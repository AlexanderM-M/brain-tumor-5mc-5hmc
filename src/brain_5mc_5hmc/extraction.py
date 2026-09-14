"""Extraction analysis stages for the brain tumour 5mC/5hmC study.

Each stage has explicit initialization and execution; importing is read-only.
See docs/workflow.md for inputs, stage order and reproduction limits.
"""

from . import runtime as _runtime

from .runtime import config as _config

from . import molecules as _molecules
from . import regional as _regional


# COUNTS


import argparse as counts_argparse, fcntl as counts_fcntl, ctypes as counts_ctypes, hashlib as counts_hashlib, json as counts_json, time as counts_time
from pathlib import Path as counts_Path
import numpy as counts_np, pysam as counts_pysam


def counts_save(path, obj):
    p = counts_Path(path).resolve()
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(counts_json.dumps(obj, indent=2))
    tmp.replace(p)


def counts_log(s):
    print(counts_time.strftime("%Y-%m-%d %H:%M:%S"), s, flush=True)


def counts_reference():
    p = counts_D / "reference/cpg_positions.npy"
    if not p.exists():
        f = counts_pysam.FastaFile(str(counts_REF))
        sites = []
        offsets = [0]
        lengths = []
        for i in range(1, 23):
            chrom = "chr" + str(i)
            a = counts_np.frombuffer(f.fetch(chrom).upper().encode(), dtype=counts_np.uint8)
            pos = counts_np.flatnonzero((a[:-1] == 67) & (a[1:] == 71)).astype(counts_np.uint32)
            sites.append(pos)
            offsets.append(offsets[-1] + len(pos))
            lengths.append(len(a))
            counts_log(chrom + ": " + str(len(pos)) + " reference CpGs")
        counts_np.save(p, counts_np.concatenate(sites))
        counts_np.save(
            counts_D / "reference/cpg_offsets.npy", counts_np.array(offsets, counts_np.uint32)
        )
        counts_np.save(
            counts_D / "reference/chromosome_lengths.npy",
            counts_np.array(lengths, counts_np.uint32),
        )
        counts_save(
            counts_D / "reference/reference.json",
            {
                "fasta": str(counts_REF),
                "fasta_bytes": counts_REF.stat().st_size,
                "fasta_mtime_ns": counts_REF.stat().st_mtime_ns,
                "fasta_fai_sha256": counts_hashlib.sha256(
                    counts_Path(str(counts_REF) + ".fai").read_bytes()
                ).hexdigest(),
                "autosomes": [f"chr{i}" for i in range(1, 23)],
                "cpg_coordinate": "0-based C position; both strands summed",
                "CpGs": offsets[-1],
            },
        )
    return tuple(
        (
            counts_np.load(counts_D / "reference" / name, mmap_mode="r")
            for name in ["cpg_positions.npy", "cpg_offsets.npy", "chromosome_lengths.npy"]
        )
    )


class counts_Engine:

    def __init__(self):
        (self.sites, self.offsets, self.lengths) = counts_reference()
        self.counts = counts_np.zeros((len(self.sites), 8), counts_np.uint32)
        self.regions = counts_np.zeros(
            (int(((self.lengths.astype(counts_np.uint64) + 99999) // 100000).sum()), 12),
            counts_np.uint64,
        )
        self.qc = counts_np.zeros(14, counts_np.uint64)
        counts_ctypes.CDLL(counts_pysam.libchtslib.__file__, mode=counts_ctypes.RTLD_GLOBAL)
        self.lib = counts_ctypes.CDLL(str(_runtime.native_path("regional_counts.so")))
        ptr = counts_ctypes.c_void_p
        self.lib.regional_create.argtypes = [ptr] * 6 + [counts_ctypes.c_uint32] * 2
        self.lib.regional_create.restype = ptr
        self.lib.regional_process.argtypes = [
            ptr,
            counts_ctypes.c_char_p,
            counts_ctypes.c_char_p,
            counts_ctypes.c_size_t,
        ]
        self.lib.regional_process.restype = counts_ctypes.c_int
        self.lib.regional_destroy.argtypes = [ptr]
        self.ctx = self.lib.regional_create(
            *[
                a.ctypes.data
                for a in [
                    self.sites,
                    self.offsets,
                    self.lengths,
                    self.counts,
                    self.regions,
                    self.qc,
                ]
            ],
            len(self.sites),
            len(self.regions),
        )
        if not self.ctx:
            raise MemoryError("Native context")

    def process(self, p):
        err = counts_ctypes.create_string_buffer(1024)
        rc = self.lib.regional_process(self.ctx, str(p).encode(), err, len(err))
        if rc:
            raise RuntimeError(err.value.decode())

    def close(self):
        if self.ctx:
            self.lib.regional_destroy(self.ctx)
            self.ctx = None


def counts_validate(bam_path=None, validation_name="native_validation"):
    e = counts_Engine()
    bam = (
        counts_Path(bam_path)
        if bam_path
        else counts_R / "provenance/validation/core/real_reads.bam"
    )
    t = counts_time.time()
    e.process(bam)
    native_sum = e.counts.sum(0, dtype=counts_np.uint64)
    expected = {}
    oracle_regions = counts_np.zeros_like(e.regions)
    fa = counts_pysam.FastaFile(str(counts_REF))
    n = 0
    ro = counts_np.concatenate(
        ([0], counts_np.cumsum((e.lengths.astype(counts_np.uint64) + 99999) // 100000))
    )

    def rejection(read):
        if read.is_secondary or read.is_supplementary:
            return True
        if read.is_unmapped or read.is_duplicate or read.is_qcfail:
            return True
        if read.reference_name not in {f"chr{i}" for i in range(1, 23)}:
            return True
        if read.mapping_quality < 20 or read.mapping_quality == 255:
            return True
        if not read.has_tag("qs") or read.get_tag("qs") < 10:
            return True
        return not read.has_tag("MM") or not read.has_tag("ML")

    with counts_pysam.AlignmentFile(str(bam), "rb") as b:
        for read in b:
            if rejection(read):
                continue
            n += 1
            rev = int(read.is_reverse)
            chrom = int(read.reference_name[3:]) - 1
            seq = read.query_sequence
            start = read.reference_start
            ref = fa.fetch(read.reference_name, max(0, start - 1), read.reference_end + 1).upper()
            origin = max(0, start - 1)
            mods = read.modified_bases
            m = dict(mods.get(("C", rev, "m"), []))
            h = dict(mods.get(("C", rev, "h"), []))
            cp = e.sites[e.offsets[chrom] : e.offsets[chrom + 1]]
            for q, r in read.get_aligned_pairs(matches_only=True):
                if q < 100 or q >= len(seq) - 100 or r <= 0 or (r >= e.lengths[chrom] - 1):
                    continue
                if (
                    not rev
                    and (seq[q] != "C" or ref[r - origin : r - origin + 2] != "CG")
                    or (rev and (seq[q] != "G" or ref[r - origin - 1 : r - origin + 1] != "CG"))
                ):
                    continue
                site = r - rev
                j = int(e.offsets[chrom]) + int(counts_np.searchsorted(cp, site))
                v = expected.setdefault(j, counts_np.zeros(8, counts_np.uint64))
                v[5] += 1
                intact = seq[q - 1] == "C" if rev else seq[q + 1] == "G"
                v[6] += int(intact)
                if q not in m or q not in h or m[q] < 0 or (h[q] < 0):
                    v[7] += int(intact)
                    continue
                pm = (m[q] + 0.5) / 256
                ph = (h[q] + 0.5) / 256
                pr = [max(0, 1 - pm - ph), pm, ph]
                win = int(counts_np.argmax(pr))
                conf = pr[win]
                v[4] += 1
                v[win if conf >= 0.8 else 3] += 1
                region = int(ro[chrom]) + site // 100000
                for ti, th in enumerate([0.7, 0.8, 0.9, 0.95]):
                    if conf >= th:
                        oracle_regions[region, ti * 3 + win] += 1
    keys = counts_np.array(sorted(expected), counts_np.int64)
    values = counts_np.array([expected[int(j)] for j in keys])
    assert counts_np.array_equal(e.counts[keys], values), "Per-site native/oracle mismatch"
    assert counts_np.array_equal(native_sum, values.sum(0)), "Extra native calls"
    assert counts_np.array_equal(e.regions, oracle_regions), "Regional thresholds mismatch"
    e.process(bam)
    assert (
        counts_np.array_equal(native_sum, e.counts.sum(0)) and int(e.qc[8]) == n
    ), "Duplicate read handling mismatch"
    result = {
        "passed": True,
        "reads": n,
        "site_fields": counts_SITE_FIELDS,
        "totals": native_sum.tolist(),
        "reference_sites_checked": len(keys),
        "native_scalar_oracle_exact_match": True,
        "four_threshold_region_match": True,
        "duplicate_reads_removed_on_repeat": int(e.qc[8]),
        "seconds": counts_time.time() - t,
        "note": "Scalar reference-sequence/aligned-pairs oracle includes missing and intact-read-CpG opportunity counts; original probability implementation separately matched modkit.",
    }
    counts_save(counts_D / "validation" / (validation_name + ".json"), result)
    counts_log(counts_json.dumps(result))
    e.close()


def counts_run_sample(entry):
    label = entry["study_label"]
    out = counts_D / "results" / label
    out.mkdir(exist_ok=True)
    done = out / "completion.json"
    if done.exists():
        counts_log(label + ": already complete")
        return
    files = counts_json.loads(counts_Path(entry["manifest"]).read_text())
    e = counts_Engine()
    start = counts_time.time()
    progress = []
    for i, f in enumerate(files):
        path = counts_Path(f["path"])
        st = path.stat()
        if st.st_size != f["bytes"] or st.st_mtime_ns != f["mtime_ns"]:
            raise ValueError("Source BAM changed before read")
        t = counts_time.time()
        e.process(path)
        after = path.stat()
        if after.st_size != st.st_size or after.st_mtime_ns != st.st_mtime_ns:
            raise ValueError("Source BAM changed during read")
        progress.append({"batch": i, "seconds": counts_time.time() - t, "bytes": f["bytes"]})
        q = dict(zip(counts_QC_FIELDS, map(int, e.qc)))
        counts_save(
            out / "progress.json",
            {
                "study_label": label,
                "completed_bams": i + 1,
                "total_bams": len(files),
                "elapsed_seconds": counts_time.time() - start,
                "qc": q,
            },
        )
        counts_log(
            f"{label}: {i + 1}/{len(files)} BAMs, {counts_time.time() - t:.1f}s, {q['qc_passing_reads']:,} retained reads"
        )
    tmp = out / "site_counts.tmp.npz"
    counts_np.savez_compressed(tmp, counts=e.counts, fields=counts_np.array(counts_SITE_FIELDS))
    tmp.replace(out / "site_counts.npz")
    counts_np.savez_compressed(
        out / "window_counts.npz",
        counts=e.regions,
        thresholds=counts_np.array([0.7, 0.8, 0.9, 0.95]),
    )
    counts_save(
        out / "qc.json",
        {
            "study_label": label,
            "diagnosis": entry["diagnosis"],
            "aligned_genome_equivalent": int(e.qc[10]) / int(e.lengths.sum()),
            "qc": dict(zip(counts_QC_FIELDS, map(int, e.qc))),
            "site_totals": dict(
                zip(counts_SITE_FIELDS, map(int, e.counts.sum(0, dtype=counts_np.uint64)))
            ),
            "batch_timings": progress,
        },
    )
    counts_save(
        done,
        {
            "study_label": label,
            "complete": True,
            "full_run": True,
            "bams": len(files),
            "elapsed_seconds": counts_time.time() - start,
            "native_source_sha256": counts_hashlib.sha256(
                _runtime.native_source("regional_counts.c").read_bytes()
            ).hexdigest(),
            "native_library_sha256": counts_hashlib.sha256(
                _runtime.native_path("regional_counts.so").read_bytes()
            ).hexdigest(),
        },
    )
    e.close()
    counts_log(label + ": FULL RUN COMPLETE")


def counts_benchmark():
    entry = next(
        (
            s
            for s in counts_json.loads((counts_D / "cohort.json").read_text())["samples"]
            if s["study_label"] == "MEN-01"
        )
    )
    record = counts_json.loads(counts_Path(entry["manifest"]).read_text())[0]
    e = counts_Engine()
    t = counts_time.time()
    e.process(record["path"])
    elapsed = counts_time.time() - t
    result = {
        "study_label": entry["study_label"],
        "bam_GB": record["bytes"] / 1000000000.0,
        "seconds": elapsed,
        "MB_per_second": record["bytes"] / 1000000.0 / elapsed,
        "qc": dict(zip(counts_QC_FIELDS, map(int, e.qc))),
        "site_totals": dict(
            zip(counts_SITE_FIELDS, map(int, e.counts.sum(0, dtype=counts_np.uint64)))
        ),
    }
    counts_save(counts_D / "validation/benchmark.json", result)
    counts_log(counts_json.dumps(result))
    e.close()


def counts_main():
    a = counts_argparse.ArgumentParser()
    a.add_argument("--validation-bam")
    a.add_argument("--validation-name", default="native_validation")
    a.add_argument("--reference", action="store_true")
    a.add_argument("--validate", action="store_true")
    a.add_argument("--benchmark", action="store_true")
    a.add_argument("--run", action="store_true")
    a.add_argument("--workers", type=int, default=4)
    args = a.parse_args()
    if args.reference:
        counts_reference()
    if args.validate:
        counts_validate(args.validation_bam, args.validation_name)
    if args.benchmark:
        counts_benchmark()
    if args.run:
        lock = (counts_D / "run.lock").open("w")
        counts_fcntl.flock(lock, counts_fcntl.LOCK_EX | counts_fcntl.LOCK_NB)
        entries = counts_json.loads((counts_D / "cohort.json").read_text())["samples"]
        entries.sort(key=lambda s: s["study_label"])
        with _runtime.process_pool(max_workers=args.workers) as pool:
            list(pool.map(counts_run_sample, entries))
        counts_save(
            counts_D / "extraction_complete.json", {"samples": len(entries), "complete": True}
        )
        counts_log("ALL FULL-RUN EXTRACTIONS COMPLETE")


def initialize_counts():
    """Initialize the counts stage once; load its declared inputs."""
    global counts_D, counts_QC_FIELDS, counts_R, counts_REF, counts_SITE_FIELDS
    if _runtime.initialized("deep_regional"):
        return
    _runtime.begin("deep_regional")
    _runtime.initialize("pooled_analysis")
    counts_R = counts_Path(_config.workspace)
    counts_D = counts_R / ".analysis"
    counts_REF = counts_Path(_config.reference)
    counts_SITE_FIELDS = [
        "C80",
        "m80",
        "h80",
        "failed80",
        "explicit_paired",
        "reference_CpG_opportunities",
        "intact_read_CpG_opportunities",
        "missing_at_intact_read_CpG",
    ]
    counts_QC_FIELDS = [
        "records",
        "secondary_supplementary",
        "unmapped",
        "duplicate_qcfail_flags",
        "low_mapq",
        "non_autosomal",
        "low_missing_qscore",
        "missing_MM_ML",
        "duplicate_read_ids",
        "qc_passing_reads",
        "aligned_bp",
        "quantization_residual_clamped",
        "reads_with_explicit_CpG",
        "reads_without_explicit_CpG",
    ]
    _runtime.finish("deep_regional")


def run_counts():
    """Execute the counts workflow stage."""
    initialize_counts()
    counts_main()


# TARGETS


from pathlib import Path as targets_Path
import fcntl as targets_fcntl, hashlib as targets_hashlib, json as targets_json, subprocess as targets_subprocess, time as targets_time, argparse as targets_argparse, concurrent.futures as targets_futures


def targets_save(p, x):
    p = targets_Path(p).resolve()
    t = p.with_suffix(p.suffix + ".tmp")
    t.write_text(targets_json.dumps(x, indent=2) + "\n")
    t.replace(p)


def targets_extract(s):
    label = s["study_label"]
    out = targets_S / "targeted_bams" / label
    out.mkdir(exist_ok=True)
    sha = targets_hashlib.sha256((targets_S / "plans/panel.bed").read_bytes()).hexdigest()
    done = out / "completion.json"
    if done.exists():
        assert targets_json.loads(done.read_text())["panel_sha256"] == sha
        return label
    plan = out / "panel_sha256.txt"
    if plan.exists():
        assert plan.read_text().strip() == sha, "Panel changed during extraction"
    else:
        plan.write_text(sha + "\n")
    files = targets_json.loads(targets_Path(s["manifest"]).read_text())
    paths = []
    t = targets_time.time()
    for i, f in enumerate(files):
        source = targets_Path(f["path"])
        st = source.stat()
        assert st.st_size == f["bytes"] and st.st_mtime_ns == f["mtime_ns"], "Source changed"
        p = out / f"batch_{i:04d}.bam"
        paths.append(p)
        if not p.exists():
            tmp = p.with_suffix(".tmp.bam")
            cmd = [
                targets_SAM,
                "view",
                "-b",
                "--no-PG",
                "-M",
                "-L",
                str(targets_S / "plans/panel.bed"),
                "-q",
                "20",
                "-F",
                "3844",
                "-o",
                str(tmp),
                str(source),
            ]
            r = targets_subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode:
                raise RuntimeError(r.stderr)
            after = source.stat()
            assert after.st_size == st.st_size and after.st_mtime_ns == st.st_mtime_ns
            tmp.replace(p)
        targets_save(
            out / "progress.json",
            {
                "study_label": label,
                "completed_bams": i + 1,
                "total_bams": len(files),
                "seconds": targets_time.time() - t,
            },
        )
    filelist = out / "merge_inputs.txt"
    filelist.write_text("".join((str(p) + "\n" for p in paths)))
    tmp = out / "reads.tmp.bam"
    final = out / "reads.bam"
    for cmd in [
        [targets_SAM, "merge", "-f", "-c", "-p", "-@", "1", "-b", str(filelist), "-o", str(tmp)],
        [targets_SAM, "quickcheck", str(tmp)],
    ]:
        r = targets_subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(r.stderr)
    tmp.replace(final)
    r = targets_subprocess.run([targets_SAM, "index", str(final)], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr)
    targets_save(
        done,
        {
            "study_label": label,
            "panel_sha256": sha,
            "bams": len(files),
            "targeted_bam_bytes": final.stat().st_size,
            "seconds": targets_time.time() - t,
        },
    )
    print(label, "TARGETED EXTRACTION COMPLETE", flush=True)
    return label


def initialize_targets():
    """Initialize the targets stage once; load its declared inputs."""
    global targets_S
    global targets_D, targets_R, targets_S, targets_SAM
    if _runtime.initialized("story_extract"):
        return
    _runtime.begin("story_extract")
    targets_R = targets_Path(_config.workspace)
    targets_D = targets_R / ".analysis"
    targets_S = targets_D / "story_analysis"
    targets_SAM = _config.samtools
    _runtime.finish("story_extract")


def run_targets():
    """Execute the targets workflow stage."""
    initialize_targets()
    global targets_args, targets_entries, targets_lock, targets_p, targets_pool
    targets_p = targets_argparse.ArgumentParser()
    targets_p.add_argument("--workers", type=int, default=20)
    targets_args = targets_p.parse_args()
    targets_lock = (targets_S / "extract.lock").open("w")
    targets_fcntl.flock(targets_lock, targets_fcntl.LOCK_EX | targets_fcntl.LOCK_NB)
    targets_entries = targets_json.loads((targets_D / "cohort.json").read_text())["samples"]
    with targets_futures.ThreadPoolExecutor(max_workers=targets_args.workers) as targets_pool:
        list(targets_pool.map(targets_extract, targets_entries))
    targets_save(
        targets_S / "results/extraction_complete.json",
        {
            "specimens": len(targets_entries),
            "completed_utc": targets_time.strftime("%Y-%m-%dT%H:%M:%SZ", targets_time.gmtime()),
        },
    )
    print("ALL TARGETED EXTRACTIONS COMPLETE", flush=True)


# FOLLOWUP TARGETS


from pathlib import Path as followup_targets_Path
import json as followup_targets_json, os as followup_targets_os, time as followup_targets_time, fcntl as followup_targets_fcntl, concurrent.futures as followup_targets_futures


def initialize_followup_targets():
    """Initialize the followup_targets stage once; load its declared inputs."""
    global targets_S
    global followup_targets_R
    if _runtime.initialized("followup_extract"):
        return
    _runtime.begin("followup_extract")
    _runtime.initialize("story_extract")
    followup_targets_R = followup_targets_Path(_config.workspace)
    targets_S = followup_targets_R / ".analysis/focused_followup"
    _runtime.finish("followup_extract")


def run_followup_targets():
    """Execute the followup_targets workflow stage."""
    initialize_followup_targets()
    global followup_targets_entries, followup_targets_lock, followup_targets_pool
    followup_targets_lock = (targets_S / "extract.lock").open("w")
    followup_targets_fcntl.flock(
        followup_targets_lock, followup_targets_fcntl.LOCK_EX | followup_targets_fcntl.LOCK_NB
    )
    targets_save(
        targets_S / "extraction_job.json",
        {
            "pid": followup_targets_os.getpid(),
            "started_utc": followup_targets_time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", followup_targets_time.gmtime()
            ),
            "workers": _config.workers,
            "log": "logs/extraction.log",
        },
    )
    followup_targets_entries = followup_targets_json.loads((targets_D / "cohort.json").read_text())[
        "samples"
    ]
    with followup_targets_futures.ThreadPoolExecutor(
        max_workers=_config.workers
    ) as followup_targets_pool:
        list(followup_targets_pool.map(targets_extract, followup_targets_entries))
    targets_save(
        targets_S / "results/extraction_complete.json",
        {
            "specimens": 20,
            "completed_utc": followup_targets_time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", followup_targets_time.gmtime()
            ),
        },
    )
    print("EXTRACTION COMPLETE", flush=True)


# UNIVERSE TARGETS


from pathlib import Path as universe_targets_Path
import json as universe_targets_json, hashlib as universe_targets_hashlib, fcntl as universe_targets_fcntl, subprocess as universe_targets_subprocess, time as universe_targets_time, os as universe_targets_os, concurrent.futures as universe_targets_futures


def universe_targets_extract(s):
    label = s["study_label"]
    out = universe_targets_S / "targeted_bams" / label
    out.mkdir(exist_ok=True)
    needed = set(
        universe_targets_json.loads(
            (universe_targets_S / "plans/needed_by_patient.json").read_text()
        )[label]
    )
    panel = [
        p
        for p in universe_targets_json.loads(
            (universe_targets_S / "plans/new_panel.json").read_text()
        )
        if p["universe_index"] in needed
    ]
    bed = out / "targets.bed"
    body = "".join(
        (
            f"{p['chromosome']}\t{p['start']}\t{p['end']}\t{p['locus_id']}\n"
            for p in sorted(panel, key=lambda p: (int(p["chromosome"][3:]), p["start"]))
        )
    )
    sha = universe_targets_hashlib.sha256(body.encode()).hexdigest()
    if (out / "completion.json").exists():
        assert (
            universe_targets_json.loads((out / "completion.json").read_text())["panel_sha256"]
            == sha
        )
        return
    if bed.exists():
        assert bed.read_text() == body
    else:
        bed.write_text(body)
    files = universe_targets_json.loads(universe_targets_Path(s["manifest"]).read_text())
    paths = []
    start = universe_targets_time.time()
    for i, f in enumerate(files):
        source = universe_targets_Path(f["path"])
        st = source.stat()
        assert st.st_size == f["bytes"] and st.st_mtime_ns == f["mtime_ns"]
        p = out / f"batch_{i:04d}.bam"
        paths.append(p)
        if not p.exists():
            tmp = p.with_suffix(".tmp.bam")
            cmd = [
                universe_targets_SAM,
                "view",
                "-b",
                "--no-PG",
                "-M",
                "-L",
                str(bed),
                "-q",
                "20",
                "-F",
                "3844",
                "-o",
                str(tmp),
                str(source),
            ]
            r = universe_targets_subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode:
                raise RuntimeError(r.stderr)
            after = source.stat()
            assert after.st_size == st.st_size and after.st_mtime_ns == st.st_mtime_ns
            tmp.replace(p)
        _regional.regional_save(
            out / "progress.json",
            {
                "study_label": label,
                "completed_bams": i + 1,
                "total_bams": len(files),
                "target_promoters": len(panel),
                "seconds": universe_targets_time.time() - start,
            },
        )
    inputs = out / "merge_inputs.txt"
    inputs.write_text("".join((str(p) + "\n" for p in paths)))
    tmp = out / "reads.tmp.bam"
    final = out / "reads.bam"
    for cmd in [
        [
            universe_targets_SAM,
            "merge",
            "-f",
            "-c",
            "-p",
            "-@",
            "1",
            "-b",
            str(inputs),
            "-o",
            str(tmp),
        ],
        [universe_targets_SAM, "quickcheck", str(tmp)],
    ]:
        rr = universe_targets_subprocess.run(cmd, capture_output=True, text=True)
        if rr.returncode:
            raise RuntimeError(rr.stderr)
    tmp.replace(final)
    universe_targets_subprocess.run([universe_targets_SAM, "index", str(final)], check=True)
    _regional.regional_save(
        out / "completion.json",
        {
            "study_label": label,
            "panel_sha256": sha,
            "target_promoters": len(panel),
            "bytes": final.stat().st_size,
            "seconds": universe_targets_time.time() - start,
        },
    )
    print(label, "EXTRACTION COMPLETE", flush=True)


def initialize_universe_targets():
    """Initialize the universe_targets stage once; load its declared inputs."""
    global universe_targets_D, universe_targets_R, universe_targets_S, universe_targets_SAM
    if _runtime.initialized("strengthen_extract"):
        return
    _runtime.begin("strengthen_extract")
    _runtime.initialize("pooled_analysis")
    universe_targets_R = universe_targets_Path(_config.workspace)
    universe_targets_D = universe_targets_R / ".analysis"
    universe_targets_S = universe_targets_D / "strengthening"
    universe_targets_SAM = _config.samtools
    _runtime.finish("strengthen_extract")


def run_universe_targets():
    """Execute the universe_targets workflow stage."""
    initialize_universe_targets()
    global universe_targets_lock, universe_targets_pool, universe_targets_samples
    universe_targets_lock = (universe_targets_S / "extract.lock").open("w")
    universe_targets_fcntl.flock(
        universe_targets_lock, universe_targets_fcntl.LOCK_EX | universe_targets_fcntl.LOCK_NB
    )
    _regional.regional_save(
        universe_targets_S / "extraction_job.json",
        {
            "pid": universe_targets_os.getpid(),
            "workers": _config.workers,
            "started": universe_targets_time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", universe_targets_time.gmtime()
            ),
        },
    )
    universe_targets_samples = universe_targets_json.loads(
        (universe_targets_D / "cohort.json").read_text()
    )["samples"]
    with universe_targets_futures.ThreadPoolExecutor(
        max_workers=_config.workers
    ) as universe_targets_pool:
        list(universe_targets_pool.map(universe_targets_extract, universe_targets_samples))
    _regional.regional_save(
        universe_targets_S / "results/extraction_complete.json", {"patients": 20}
    )
    print("EXTRACTION COMPLETE", flush=True)


# MATRICES


from pathlib import Path as matrices_Path
import json as matrices_json, time as matrices_time, os as matrices_os, fcntl as matrices_fcntl
import numpy as matrices_np, pysam as matrices_pysam


def matrices_sample(label):
    out = matrices_S / "matrices" / label
    out.mkdir(parents=True, exist_ok=True)
    if (out / "completion.json").exists():
        return matrices_json.loads((out / "completion.json").read_text())
    needed = set(
        matrices_json.loads((matrices_S / "plans/needed_by_patient.json").read_text())[label]
    )
    panel = [
        p
        for p in matrices_json.loads((matrices_S / "plans/new_panel.json").read_text())
        if p["universe_index"] in needed
    ]
    sites = matrices_np.load(matrices_D / "reference/cpg_positions.npy", mmap_mode="r")
    off = matrices_np.load(matrices_D / "reference/cpg_offsets.npy")
    common = matrices_np.load(matrices_D / "pooled_analysis/common_cpg_mask_depth5.npy")
    indices = {}
    mat = {}
    ids = {}
    quality = {}
    bychrom = {}
    for p in panel:
        ch = int(p["chromosome"][3:]) - 1
        (lo, hi) = map(int, off[ch : ch + 2])
        (a, b) = matrices_np.searchsorted(sites[lo:hi], [p["start"], p["end"]])
        k = p["locus_id"]
        indices[k] = matrices_np.arange(lo + a, lo + b)
        mat[k] = []
        ids[k] = []
        quality[k] = []
        bychrom.setdefault(p["chromosome"], []).append(p)
    intervals = {}
    for ch, pp in bychrom.items():
        pp.sort(key=lambda p: p["start"])
        intervals[ch] = (
            pp,
            matrices_np.array([p["start"] for p in pp]),
            matrices_np.maximum.accumulate([p["end"] for p in pp]),
        )
    seen = set()
    start = matrices_time.time()
    with matrices_pysam.AlignmentFile(
        str(matrices_S / "targeted_bams" / label / "reads.bam"), "rb"
    ) as bam:
        for r in bam:
            if not _molecules.blocks_eligible(r) or r.query_name in seen:
                continue
            seen.add(r.query_name)
            if r.reference_name not in intervals:
                continue
            (pp, starts, maxends) = intervals[r.reference_name]
            a = int(matrices_np.searchsorted(maxends, r.reference_start, side="right"))
            b = int(matrices_np.searchsorted(starts, r.reference_end, side="left"))
            for p in pp[a:b]:
                if p["end"] <= r.reference_start:
                    continue
                k = p["locus_id"]
                mat[k].append(_molecules.blocks_states_at(r, sites[indices[k]]))
                ids[k].append(r.query_name)
                quality[k].append(r.mapping_quality)
    full = matrices_np.load(matrices_D / "results" / label / "site_counts.npz")["counts"]
    checks = []
    for p in panel:
        k = p["locus_id"]
        ix = indices[k]
        x = matrices_np.asarray(mat[k], matrices_np.int8).reshape(-1, len(ix))
        obs = matrices_np.column_stack([(x == v).sum(0) for v in [0, 1, 2]])
        assert matrices_np.array_equal(obs, full[ix, :3]), label + " " + k + " full-count mismatch"
        tmp = out / (k + ".tmp.npz")
        matrices_np.savez_compressed(
            tmp,
            states=x,
            read_ids=matrices_np.asarray(ids[k]),
            mapq=matrices_np.asarray(quality[k], matrices_np.uint8),
            cpg_positions=sites[ix],
            common_mask=common[ix],
        )
        tmp.replace(out / (k + ".npz"))
        checks.append(
            {"locus_id": k, "CpGs": len(ix), "molecules": len(x), "exact_count_match": True}
        )
    result = {
        "study_label": label,
        "loci": len(panel),
        "unique_reads": len(seen),
        "all_site_counts_match": True,
        "seconds": matrices_time.time() - start,
        "checks": checks,
    }
    _regional.regional_save(out / "completion.json", result)
    print(label, "MATRICES COMPLETE", len(panel), round(matrices_time.time() - start), flush=True)
    return result


def initialize_matrices():
    """Initialize the matrices stage once; load its declared inputs."""
    global matrices_D, matrices_R, matrices_S
    if _runtime.initialized("strengthen_matrices"):
        return
    _runtime.begin("strengthen_matrices")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    matrices_R = matrices_Path(_config.workspace)
    matrices_D = matrices_R / ".analysis"
    matrices_S = matrices_D / "strengthening"
    _runtime.finish("strengthen_matrices")


def run_matrices():
    """Execute the matrices workflow stage."""
    initialize_matrices()
    global matrices_done, matrices_f, matrices_jobs, matrices_label, matrices_labels, matrices_lock, matrices_pending, matrices_pool, matrices_results
    matrices_lock = (matrices_S / "matrix.lock").open("w")
    matrices_fcntl.flock(matrices_lock, matrices_fcntl.LOCK_EX | matrices_fcntl.LOCK_NB)
    _regional.regional_save(
        matrices_S / "matrix_job.json", {"pid": matrices_os.getpid(), "workers": _config.workers}
    )
    matrices_labels = sorted(
        matrices_json.loads((matrices_S / "plans/needed_by_patient.json").read_text())
    )
    matrices_pending = set(matrices_labels)
    matrices_jobs = {}
    matrices_results = []
    with _runtime.process_pool(max_workers=_config.workers) as matrices_pool:
        while matrices_pending or matrices_jobs:
            for matrices_label in sorted(matrices_pending):
                if (matrices_S / "targeted_bams" / matrices_label / "completion.json").exists():
                    matrices_jobs[matrices_pool.submit(matrices_sample, matrices_label)] = (
                        matrices_label
                    )
                    matrices_pending.remove(matrices_label)
            matrices_done = [f for f in matrices_jobs if f.done()]
            for matrices_f in matrices_done:
                matrices_results.append(matrices_f.result())
                del matrices_jobs[matrices_f]
            if matrices_pending or matrices_jobs:
                matrices_time.sleep(10)
    _regional.regional_save(
        matrices_S / "results/matrices_complete.json",
        {
            "patients": len(matrices_results),
            "new_patient_loci": sum((r["loci"] for r in matrices_results)),
            "all_site_counts_match": True,
        },
    )
    print("ALL MATRICES COMPLETE", flush=True)
