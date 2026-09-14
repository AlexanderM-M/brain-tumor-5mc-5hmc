"""Molecules analysis stages for the brain tumour 5mC/5hmC study.

Each stage has explicit initialization and execution; importing is read-only.
See docs/workflow.md for inputs, stage order and reproduction limits.
"""

from . import runtime as _runtime

from .runtime import config as _config

from . import regional as _regional


# BLOCKS


from pathlib import Path as blocks_Path
import json as blocks_json, time as blocks_time, hashlib as blocks_hashlib
import numpy as blocks_np, pysam as blocks_pysam
import matplotlib as blocks_matplotlib
import matplotlib.pyplot as blocks_plt
from matplotlib.text import Text as blocks_Text
import re as blocks_re


def blocks_figure(fig, name):
    for t in fig.findobj(match=blocks_Text):
        assert not blocks_re.search("\\bN(?:20)?\\d{2}[._-]\\d+\\b", t.get_text(), blocks_re.I)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(blocks_S / "figures" / f"{name}.{ext}", dpi=180, bbox_inches="tight")
    blocks_plt.close(fig)


def blocks_eligible(r):
    return (
        not r.flag & 3844
        and 20 <= r.mapping_quality < 255
        and r.has_tag("qs")
        and blocks_np.isfinite(r.get_tag("qs"))
        and (r.get_tag("qs") >= 10)
        and r.has_tag("MM")
        and r.has_tag("ML")
    )


def blocks_states_at(r, cpg):
    """Decode explicit paired MM/ML calls at reference CpGs.

    Returns -1 for missing/failed calls, 0 for canonical C, 1 for 5mC and 2
    for 5hmC. Trims 100 query bases at each end and requires probability 0.80."""
    states = blocks_np.full(len(cpg), -1, blocks_np.int8)
    rev = int(r.is_reverse)
    seq = r.query_sequence
    mods = r.modified_bases
    m = dict(mods.get(("C", rev, "m"), []))
    h = dict(mods.get(("C", rev, "h"), []))
    q = 0
    rp = r.reference_start
    for op, length in r.cigartuples:
        if op in [0, 7, 8]:
            (lo, hi) = blocks_np.searchsorted(cpg, [rp - rev, rp + length - rev])
            for j in range(lo, hi):
                qp = q + int(cpg[j]) + rev - rp
                if qp < 100 or qp >= len(seq) - 100 or seq[qp] != ("G" if rev else "C"):
                    continue
                if qp not in m or qp not in h or m[qp] < 0 or (h[qp] < 0):
                    continue
                assert m[qp] + h[qp] <= 256
                pm = (m[qp] + 0.5) / 256
                ph = (h[qp] + 0.5) / 256
                probs = [max(0.0, 1 - pm - ph), pm, ph]
                win = int(blocks_np.argmax(probs))
                if probs[win] >= 0.8:
                    states[j] = win
            q += length
            rp += length
        elif op in [1, 4]:
            q += length
        elif op in [2, 3]:
            rp += length
    return states


def blocks_patterns_metrics(x, rng, repeats=64, n=12):
    """Compare three-CpG patterns with independent column permutations.

    Each of 64 draws uses 12 complete molecules. Returns observed/null
    entropy, entropy excess, joint-5hmC excess and multi-5hmC excess."""
    values = []
    for _ in range(repeats):
        z = x[rng.choice(len(x), n, replace=False)]
        per = blocks_np.column_stack([rng.permutation(z[:, j]) for j in range(3)])

        def metrics(a):
            code = a[:, 0] + 3 * a[:, 1] + 9 * a[:, 2]
            freq = blocks_np.bincount(code, minlength=27) / len(a)
            freq = freq[freq > 0]
            entropy = -blocks_np.sum(freq * blocks_np.log2(freq))
            h = a == 2
            joint = blocks_np.mean(
                [blocks_np.mean(h[:, i] & h[:, j]) for (i, j) in [(0, 1), (0, 2), (1, 2)]]
            )
            return (entropy, joint, blocks_np.mean(h.sum(1) >= 2))

        (eo, jo, ho) = metrics(z)
        (en, jn, hn) = metrics(per)
        values.append([eo, en, eo - en, jo - jn, ho - hn])
    return blocks_np.mean(values, axis=0)


def blocks_sample(s):
    label = s["study_label"]
    panel = blocks_json.loads((blocks_S / "plans/panel.json").read_text())
    sites = blocks_np.load(blocks_D / "reference/cpg_positions.npy", mmap_mode="r")
    off = blocks_np.load(blocks_D / "reference/cpg_offsets.npy")
    common = blocks_np.load(blocks_D / "pooled_analysis/common_cpg_mask_depth5.npy")
    matrices = {}
    ids = {}
    cp = {}
    indices = {}
    for p in panel:
        ch = int(p["chromosome"][3:]) - 1
        (lo, hi) = map(int, off[ch : ch + 2])
        (a, b) = blocks_np.searchsorted(sites[lo:hi], [p["start"], p["end"]])
        ix = blocks_np.arange(lo + a, lo + b)
        cp[p["locus_id"]] = sites[ix]
        indices[p["locus_id"]] = ix
        matrices[p["locus_id"]] = []
        ids[p["locus_id"]] = []
    seen = set()
    seen_qc = 0
    with blocks_pysam.AlignmentFile(
        str(blocks_S / "targeted_bams" / label / "reads.bam"), "rb"
    ) as b:
        for r in b:
            if not blocks_eligible(r) or r.query_name in seen:
                continue
            seen.add(r.query_name)
            seen_qc += 1
            for p in panel:
                if (
                    r.reference_name != p["chromosome"]
                    or r.reference_start >= p["end"]
                    or r.reference_end <= p["start"]
                ):
                    continue
                x = blocks_states_at(r, cp[p["locus_id"]])
                matrices[p["locus_id"]].append(x)
                ids[p["locus_id"]].append(r.query_name)
    target = blocks_S / "results" / label
    target.mkdir(exist_ok=True)
    full = blocks_np.load(blocks_D / "results" / label / "site_counts.npz")["counts"]
    blockrows = []
    locusrows = []
    checks = []
    for p in panel:
        key = p["locus_id"]
        x = blocks_np.array(matrices[key], dtype=blocks_np.int8).reshape((-1, len(cp[key])))
        ix = indices[key]
        observed = blocks_np.column_stack([(x == v).sum(0) for v in [0, 1, 2]])
        assert blocks_np.array_equal(
            observed, full[ix, :3]
        ), f"{label} {key}: read-level counts differ from full extraction"
        checks.append({"locus_id": key, "CpGs": len(ix), "exact_count_match": True})
        blocks_np.savez_compressed(
            target / (key + ".npz"),
            states=x,
            read_ids=blocks_np.array(ids[key]),
            cpg_positions=cp[key],
            common_mask=common[ix],
        )
        if p["kind"] == "amplification_screen":
            continue
        ci = blocks_np.flatnonzero(common[ix])
        starts = blocks_np.arange(0, len(ci) - 2, 3)
        starts = starts[
            blocks_np.unique(
                blocks_np.linspace(0, len(starts) - 1, min(12, len(starts))).round().astype(int)
            )
        ]
        usable = []
        for bi, start in enumerate(starts):
            col = ci[start : start + 3]
            z = x[:, col]
            z = z[(z >= 0).all(1)]
            n = len(z)
            row = {
                "study_label": label,
                "diagnosis": s["diagnosis"],
                "locus_id": key,
                "gene": p["gene"],
                "kind": p["kind"],
                "pair_id": p["pair_id"],
                "block": bi,
                "CpG_positions": ",".join(map(str, cp[key][col])),
                "fully_observed_molecules": n,
                "eligible": n >= 12,
            }
            if n >= 12:
                seed = int.from_bytes(
                    blocks_hashlib.sha256((label + key + str(bi)).encode()).digest()[:8], "little"
                )
                v = blocks_patterns_metrics(z, blocks_np.random.default_rng(seed))
                row.update(
                    entropy_bits=float(v[0]),
                    permuted_entropy_bits=float(v[1]),
                    entropy_excess_bits=float(v[2]),
                    joint_5hmC_excess_pp=float(v[3] * 100),
                    two_or_more_5hmC_excess_pp=float(v[4] * 100),
                )
                usable.append(row)
            else:
                row.update(
                    entropy_bits="",
                    permuted_entropy_bits="",
                    entropy_excess_bits="",
                    joint_5hmC_excess_pp="",
                    two_or_more_5hmC_excess_pp="",
                )
            blockrows.append(row)
        locusrows.append(
            {
                "study_label": label,
                "diagnosis": s["diagnosis"],
                "locus_id": key,
                "gene": p["gene"],
                "kind": p["kind"],
                "pair_id": p["pair_id"],
                "retained_molecules": len(x),
                "tested_blocks": len(starts),
                "eligible_blocks": len(usable),
                "eligible_locus": len(usable) >= 3,
                "mean_joint_5hmC_excess_pp": (
                    float(blocks_np.mean([r["joint_5hmC_excess_pp"] for r in usable]))
                    if len(usable) >= 3
                    else ""
                ),
                "mean_entropy_excess_bits": (
                    float(blocks_np.mean([r["entropy_excess_bits"] for r in usable]))
                    if len(usable) >= 3
                    else ""
                ),
            }
        )
    _regional.regional_table(target / "blocks.tsv", blockrows)
    _regional.regional_table(target / "loci.tsv", locusrows)
    _regional.regional_save(
        target / "validation.json",
        {
            "study_label": label,
            "exact_reference_CpG_counts_match_full_extraction": True,
            "loci": checks,
            "unique_QC_passing_reads": seen_qc,
        },
    )
    print(label, "MOLECULE ANALYSIS COMPLETE", flush=True)
    return (locusrows, blockrows)


def blocks_summarize(samples, all_loci, all_blocks):
    _regional.regional_table(blocks_S / "results/Q3_locus_metrics.tsv", all_loci)
    _regional.regional_table(blocks_S / "results/Q3_block_metrics.tsv", all_blocks)
    panel = [
        p
        for p in blocks_json.loads((blocks_S / "plans/panel.json").read_text())
        if p["kind"] != "amplification_screen"
    ]
    labels = [s["study_label"] for s in samples]
    lookup = {(r["study_label"], r["locus_id"]): r for r in all_loci}
    cov = blocks_np.array(
        [[lookup[l, p["locus_id"]]["eligible_blocks"] for l in labels] for p in panel]
    )
    (fig, ax) = blocks_plt.subplots(figsize=(12, 7))
    im = ax.imshow(cov, aspect="auto", vmin=0, vmax=12, cmap="viridis")
    ax.set_xticks(range(20), labels, rotation=60)
    ax.set_yticks(range(len(panel)), [p["locus_id"] + " " + p["gene"] for p in panel], fontsize=8)
    ax.set_title("Molecule-level coverage: fixed 3-CpG blocks with at least 12 complete reads")
    fig.colorbar(im, ax=ax, label="Eligible blocks (maximum 12)")
    fig.tight_layout()
    blocks_figure(fig, "Q3_01_molecule_coverage")
    specimen = []
    for s in samples:
        paired = []
        for n in range(1, 13):
            c = lookup[s["study_label"], f"C{n:02d}"]
            r = lookup[s["study_label"], f"R{n:02d}"]
            if c["eligible_locus"] and r["eligible_locus"]:
                paired.append((c, r))
        row = {
            "study_label": s["study_label"],
            "diagnosis": s["diagnosis"],
            "eligible_matched_pairs": len(paired),
        }
        for kind, k in [("candidate", 0), ("comparison", 1)]:
            row[kind + "_joint_excess_pp"] = (
                float(blocks_np.mean([p[k]["mean_joint_5hmC_excess_pp"] for p in paired]))
                if paired
                else None
            )
            row[kind + "_entropy_excess_bits"] = (
                float(blocks_np.mean([p[k]["mean_entropy_excess_bits"] for p in paired]))
                if paired
                else None
            )
        specimen.append(row)
    _regional.regional_table(blocks_S / "results/Q3_specimen_metrics.tsv", specimen)
    colors = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    (fig, axes) = blocks_plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, diag in zip(axes, colors):
        for r in specimen:
            if r["diagnosis"] != diag or r["eligible_matched_pairs"] < 6:
                continue
            ax.plot(
                [0, 1],
                [r["candidate_joint_excess_pp"], r["comparison_joint_excess_pp"]],
                "-o",
                color=colors[diag],
                alpha=0.6,
                ms=4,
            )
        ax.axhline(0, c="gray", ls="--", lw=0.8)
        ax.set_xticks([0, 1], ["Opposing-change\npromoters", "Matched comparison\npromoters"])
        ax.set_ylabel("Excess co-occurring 5hmC (percentage points)")
        ax.set_title(diag)
    fig.suptitle(
        "Molecular 5hmC co-occurrence beyond site-frequency-preserving permutations\nEach line is one specimen; 12 molecules per block in every resample",
        fontsize=11,
    )
    fig.tight_layout()
    blocks_figure(fig, "Q3_02_specimen_molecular_patterns")
    locus = "C02"
    p = next((p for p in panel if p["locus_id"] == locus))
    (fig, axes) = blocks_plt.subplots(4, 5, figsize=(15, 10))
    from matplotlib.colors import ListedColormap, BoundaryNorm

    cmap = ListedColormap(["#eeeeee", "#cfcfcf", "#497ba6", "#df9d29"])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], 4)
    for ax, s in zip(axes.ravel(), samples):
        a = blocks_np.load(blocks_S / "results" / s["study_label"] / (locus + ".npz"))
        x = a["states"][:, a["common_mask"]]
        x = x[(x >= 0).sum(1) >= 8]
        rng = blocks_np.random.default_rng(20260912)
        if len(x) > 60:
            x = x[blocks_np.sort(rng.choice(len(x), 60, replace=False))]
        if len(x):
            ax.imshow(x, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
        ax.set_title(s["study_label"], color=colors[s["diagnosis"]], fontsize=10)
        ax.set_xticks([])
        ax.set_ylabel("DNA molecules", fontsize=8)
    fig.suptitle(
        p["gene"]
        + " promoter: individual molecules at common CpGs\nGrey: canonical C; blue: 5mC; gold: 5hmC; pale: missing/filtered. Up to 60 randomly selected reads per specimen.",
        fontsize=11,
    )
    fig.tight_layout()
    blocks_figure(fig, "Q3_03_P2RX1_read_patterns")
    usable = [r for r in specimen if r["eligible_matched_pairs"] >= 6]
    stats = {
        "specimens": 20,
        "specimens_with_at_least_six_matched_pairs": len(usable),
        "eligible_locus_specimen_combinations": sum((r["eligible_locus"] for r in all_loci)),
        "total_locus_specimen_combinations": len(all_loci),
        "eligible_fixed_blocks": sum((r["eligible"] for r in all_blocks)),
        "total_fixed_blocks": len(all_blocks),
        "mean_candidate_minus_comparison_joint_excess_pp": {
            diag: float(
                blocks_np.mean(
                    [
                        r["candidate_joint_excess_pp"] - r["comparison_joint_excess_pp"]
                        for r in usable
                        if r["diagnosis"] == diag
                    ]
                )
            )
            for diag in colors
        },
        "all_20_read_level_count_validations_passed": True,
    }
    _regional.regional_save(blocks_S / "results/Q3_summary.json", stats)
    lines = [
        "# Q3: individual DNA molecules",
        "",
        "Analysis completed on the fixed panel of 12 opposing-change and 12 matched comparison promoters. Four extra amplification-screen intervals were extracted for later questions.",
        "",
        "## Numerical results",
        "",
        blocks_json.dumps(stats, indent=2),
        "",
        "## What this analysis answers",
        "",
        "The BAM extracts preserve molecule linkage, and every extracted locus matched the original full-run passing per-CpG counts exactly in all 20 specimens. Fixed groups of three common CpGs were examined only when at least 12 molecules had all three states observed. Up to 12 non-overlapping blocks were chosen by genomic order in each promoter. Locus summaries require at least three eligible blocks; the paired specimen comparison requires at least six eligible locus pairs.",
        "",
        "For each block, 64 resamples without replacement use exactly 12 molecules. A separate permutation of each CpG column preserves that column's modification frequencies while breaking linkage between neighbouring CpGs. Excess joint 5hmC is the observed average pairwise co-occurrence minus this permuted baseline. Positive values indicate clustering of 5hmC on the same molecules beyond local marginal frequencies. Negative entropy excess indicates patterns more ordered than the corresponding permuted states. Neither metric identifies tumour clones, cell types or dynamic conversion of 5mC to 5hmC.",
        "",
        "The specimen-level paired candidate/comparison summaries and the P2RX1 read diagrams show the molecular structure underlying the regional averages. These are descriptive in-cohort follow-up results. Selection of loci by prior group differences and limited meningioma numbers prevent treating this as independent biomarker validation. Variable DNA mixture, haplotypes and copy number remain alternative explanations to investigate in Q4.",
        "",
        "## Figures",
        "",
        "- Q3_01_molecule_coverage: which specimen/locus comparisons have sufficient molecules.",
        "- Q3_02_specimen_molecular_patterns: matched, specimen-level co-occurrence summaries.",
        "- Q3_03_P2RX1_read_patterns: molecule patterns for all specimens, using study labels.",
    ]
    (blocks_S / "results/Q3_ANSWER.md").write_text("\n".join(lines) + "\n")
    _regional.regional_save(
        blocks_S / "results/Q3_complete.json",
        {
            "completed_utc": blocks_time.strftime("%Y-%m-%dT%H:%M:%SZ", blocks_time.gmtime()),
            "analysis_and_figures_generated": True,
            "requires_visual_review_before_checklist_tick": True,
        },
    )
    print("Q3 COMPLETE", blocks_json.dumps(stats), flush=True)


def initialize_blocks():
    """Initialize the blocks stage once; load its declared inputs."""
    global blocks_S
    global blocks_D, blocks_R, blocks_S
    if _runtime.initialized("story_molecules"):
        return
    _runtime.begin("story_molecules")
    _runtime.initialize("pooled_analysis")
    blocks_matplotlib.use("Agg")
    blocks_R = blocks_Path(_config.workspace)
    blocks_D = blocks_R / ".analysis"
    blocks_S = blocks_D / "story_analysis"
    _runtime.finish("story_molecules")


def run_blocks():
    """Execute the blocks workflow stage."""
    initialize_blocks()
    global blocks_pool, blocks_results, blocks_samples
    assert (
        blocks_S / "results/extraction_complete.json"
    ).exists(), "Target extraction is incomplete"
    blocks_samples = sorted(
        blocks_json.loads((blocks_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    with _runtime.process_pool(max_workers=_config.workers) as blocks_pool:
        blocks_results = list(blocks_pool.map(blocks_sample, blocks_samples))
    blocks_summarize(
        blocks_samples,
        [r for (rs, bs) in blocks_results for r in rs],
        [r for (rs, bs) in blocks_results for r in bs],
    )


# FOLLOWUP BLOCKS


from pathlib import Path as followup_blocks_Path
import json as followup_blocks_json, hashlib as followup_blocks_hashlib, time as followup_blocks_time, os as followup_blocks_os
import numpy as followup_blocks_np, pysam as followup_blocks_pysam


def followup_blocks_joint_excess(x, rng, repeats=64):
    values = []
    for _ in range(repeats):
        z = x[rng.choice(len(x), 12, replace=False)]
        p = followup_blocks_np.column_stack([rng.permutation(z[:, j]) for j in range(3)])
        a = (z == 2).sum(1)
        b = (p == 2).sum(1)
        values.append(float(followup_blocks_np.sum(a * (a - 1) - b * (b - 1))) / 72)
    return float(followup_blocks_np.mean(values) * 100)


def followup_blocks_sample(s):
    label = s["study_label"]
    out = followup_blocks_F / "results" / label
    out.mkdir(exist_ok=True)
    if (out / "molecules_complete.json").exists():
        import csv

        with (out / "loci.tsv").open() as f:
            loci = list(csv.DictReader(f, delimiter="\t"))
        with (out / "blocks.tsv").open() as f:
            blocks = list(csv.DictReader(f, delimiter="\t"))
        return (loci, blocks)
    panel = followup_blocks_json.loads((followup_blocks_F / "plans/panel.json").read_text())
    positions = followup_blocks_np.load(
        followup_blocks_D / "reference/cpg_positions.npy", mmap_mode="r"
    )
    offset = followup_blocks_np.load(followup_blocks_D / "reference/cpg_offsets.npy")
    common = followup_blocks_np.load(
        followup_blocks_D / "pooled_analysis/common_cpg_mask_depth5.npy"
    )
    indices = {}
    mat = {}
    ids = {}
    quality = {}
    bychrom = {}
    for p in panel:
        (lo, hi) = map(int, offset[int(p["chromosome"][3:]) - 1 : int(p["chromosome"][3:]) + 1])
        (a, b) = followup_blocks_np.searchsorted(positions[lo:hi], [p["start"], p["end"]])
        ix = followup_blocks_np.arange(lo + a, lo + b)
        indices[p["locus_id"]] = ix
        mat[p["locus_id"]] = []
        ids[p["locus_id"]] = []
        quality[p["locus_id"]] = []
        bychrom.setdefault(p["chromosome"], []).append(p)
    seen = set()
    with followup_blocks_pysam.AlignmentFile(
        str(followup_blocks_F / "targeted_bams" / label / "reads.bam"), "rb"
    ) as bam:
        for r in bam:
            if not blocks_eligible(r) or r.query_name in seen:
                continue
            seen.add(r.query_name)
            for p in bychrom.get(r.reference_name, []):
                if r.reference_start >= p["end"] or r.reference_end <= p["start"]:
                    continue
                k = p["locus_id"]
                mat[k].append(blocks_states_at(r, positions[indices[k]]))
                ids[k].append(r.query_name)
                quality[k].append(r.mapping_quality)
    full = followup_blocks_np.load(followup_blocks_D / "results" / label / "site_counts.npz")[
        "counts"
    ]
    loci = []
    blocks = []
    checks = []
    for p in panel:
        k = p["locus_id"]
        ix = indices[k]
        x = followup_blocks_np.asarray(mat[k], followup_blocks_np.int8).reshape(-1, len(ix))
        obs = followup_blocks_np.column_stack([(x == state).sum(0) for state in [0, 1, 2]])
        assert followup_blocks_np.array_equal(obs, full[ix, :3]), (
            label + " " + k + " site-count mismatch"
        )
        checks.append({"locus_id": k, "exact_site_counts": True})
        followup_blocks_np.savez_compressed(
            out / (k + ".npz"),
            states=x,
            read_ids=followup_blocks_np.array(ids[k]),
            mapq=followup_blocks_np.array(quality[k], followup_blocks_np.uint8),
            cpg_positions=positions[ix],
            common_mask=common[ix],
        )
        ci = followup_blocks_np.flatnonzero(common[ix])
        starts = followup_blocks_np.arange(0, len(ci) - 2, 3)
        starts = starts[
            followup_blocks_np.unique(
                followup_blocks_np.linspace(0, len(starts) - 1, min(12, len(starts)))
                .round()
                .astype(int)
            )
        ]
        scores = []
        for bi, start in enumerate(starts):
            col = ci[start : start + 3]
            z = x[:, col]
            z = z[(z >= 0).all(1)]
            row = {
                "study_label": label,
                "diagnosis": s["diagnosis"],
                "locus_id": k,
                "block": bi,
                "CpG_positions": ",".join(map(str, positions[ix[col]])),
                "complete_molecules": len(z),
                "eligible": len(z) >= 12,
                "joint_5hmC_excess_pp": "",
            }
            if len(z) >= 12:
                seedkey = p["original_locus_id"] or k
                seed = int.from_bytes(
                    followup_blocks_hashlib.sha256((label + seedkey + str(bi)).encode()).digest()[
                        :8
                    ],
                    "little",
                )
                score = followup_blocks_joint_excess(z, followup_blocks_np.random.default_rng(seed))
                scores.append(score)
                row["joint_5hmC_excess_pp"] = score
            blocks.append(row)
        loci.append(
            {
                "study_label": label,
                "diagnosis": s["diagnosis"],
                "locus_id": k,
                "gene": p["gene"],
                "blind_panel": p["blind_panel"],
                "original_locus_id": p["original_locus_id"],
                "eligible_blocks": len(scores),
                "eligible_locus": len(scores) >= 3,
                "mean_joint_5hmC_excess_pp": (
                    float(followup_blocks_np.mean(scores)) if len(scores) >= 3 else ""
                ),
            }
        )
    _regional.regional_table(out / "loci.tsv", loci)
    _regional.regional_table(out / "blocks.tsv", blocks)
    _regional.regional_save(
        out / "molecules_complete.json",
        {"loci": len(panel), "all_site_counts_match_full_extraction": True, "checks": checks},
    )
    print(label, "MOLECULES COMPLETE", flush=True)
    return (loci, blocks)


def initialize_followup_blocks():
    """Initialize the followup_blocks stage once; load its declared inputs."""
    global followup_blocks_D, followup_blocks_F, followup_blocks_R
    if _runtime.initialized("followup_molecules"):
        return
    _runtime.begin("followup_molecules")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    followup_blocks_R = followup_blocks_Path(_config.workspace)
    followup_blocks_D = followup_blocks_R / ".analysis"
    followup_blocks_F = followup_blocks_D / "focused_followup"
    _runtime.finish("followup_molecules")


def run_followup_blocks():
    """Execute the followup_blocks workflow stage."""
    initialize_followup_blocks()
    global followup_blocks_actual, followup_blocks_expected, followup_blocks_out, followup_blocks_pool, followup_blocks_samples, followup_blocks_x
    assert (followup_blocks_F / "results/extraction_complete.json").exists()
    followup_blocks_x = followup_blocks_np.random.default_rng(7).integers(0, 3, (25, 3))
    followup_blocks_expected = (
        blocks_patterns_metrics(followup_blocks_x, followup_blocks_np.random.default_rng(42))[3]
        * 100
    )
    followup_blocks_actual = followup_blocks_joint_excess(
        followup_blocks_x, followup_blocks_np.random.default_rng(42)
    )
    assert abs(followup_blocks_expected - followup_blocks_actual) < 1e-12
    _regional.regional_save(
        followup_blocks_F / "validation/fast_metric_equivalence.json",
        {
            "original_metric": float(followup_blocks_expected),
            "optimized_metric": followup_blocks_actual,
            "exact_equivalence": True,
        },
    )
    _regional.regional_save(
        followup_blocks_F / "molecule_job.json",
        {
            "pid": followup_blocks_os.getpid(),
            "workers": _config.workers,
            "started_utc": followup_blocks_time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", followup_blocks_time.gmtime()
            ),
        },
    )
    followup_blocks_samples = sorted(
        followup_blocks_json.loads((followup_blocks_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    with _runtime.process_pool(max_workers=_config.workers) as followup_blocks_pool:
        followup_blocks_out = list(
            followup_blocks_pool.map(followup_blocks_sample, followup_blocks_samples)
        )
    _regional.regional_table(
        followup_blocks_F / "results/all_locus_metrics.tsv",
        [r for (a, b) in followup_blocks_out for r in a],
    )
    _regional.regional_table(
        followup_blocks_F / "results/all_block_metrics.tsv",
        [r for (a, b) in followup_blocks_out for r in b],
    )
    _regional.regional_save(
        followup_blocks_F / "results/molecule_analysis_complete.json",
        {
            "specimens": 20,
            "loci": len(
                followup_blocks_json.loads((followup_blocks_F / "plans/panel.json").read_text())
            ),
            "all_site_counts_validated": True,
        },
    )
    print("ALL MOLECULES COMPLETE", flush=True)


# SELECTION SUMMARY


from pathlib import Path as selection_summary_Path
import json as selection_summary_json, csv as selection_summary_csv, numpy as selection_summary_np
import matplotlib.pyplot as selection_summary_plt


def selection_summary_run():
    assert (selection_summary_F / "results/molecule_analysis_complete.json").exists()
    panel = selection_summary_json.loads((selection_summary_F / "plans/panel.json").read_text())
    folds = selection_summary_json.loads((selection_summary_F / "plans/folds.json").read_text())
    rows = list(
        selection_summary_csv.DictReader(
            (selection_summary_F / "results/all_locus_metrics.tsv").open(), delimiter="\t"
        )
    )
    lookup = {(r["study_label"], r["locus_id"]): r for r in rows}
    labels = [f["held_out"] for f in folds]
    old = {p["original_locus_id"]: p["locus_id"] for p in panel if p["original_locus_id"]}
    scores = []
    freq = {p["locus_id"]: 0 for p in panel}
    pairrows = []
    for fold in folds:
        label = fold["held_out"]
        assert label not in fold["training_labels"]
        diag = lookup[label, panel[0]["locus_id"]]["diagnosis"]
        cv = []
        original = []
        for k, pair in enumerate(fold["pairs"]):
            a = lookup[label, pair["candidate_locus_id"]]
            b = lookup[label, pair["comparison_locus_id"]]
            freq[pair["candidate_locus_id"]] += 1
            valid = a["eligible_locus"] == "True" and b["eligible_locus"] == "True"
            diff = (
                float(a["mean_joint_5hmC_excess_pp"]) - float(b["mean_joint_5hmC_excess_pp"])
                if valid
                else None
            )
            if valid:
                cv.append(diff)
            pairrows.append(
                {
                    "study_label": label,
                    "pair": k + 1,
                    "candidate": pair["candidate_locus_id"],
                    "comparison": pair["comparison_locus_id"],
                    "eligible": valid,
                    "heldout_difference_pp": diff,
                    "training_matching_distance": pair["matching_distance"],
                }
            )
        for i in range(1, 13):
            a = lookup[label, old[f"C{i:02d}"]]
            b = lookup[label, old[f"R{i:02d}"]]
            if a["eligible_locus"] == "True" and b["eligible_locus"] == "True":
                original.append(
                    float(a["mean_joint_5hmC_excess_pp"]) - float(b["mean_joint_5hmC_excess_pp"])
                )
        scores.append(
            {
                "study_label": label,
                "diagnosis": diag,
                "selected_pairs": len(fold["pairs"]),
                "eligible_heldout_pairs": len(cv),
                "eligible_heldout_score": len(cv) >= 6,
                "heldout_candidate_minus_comparison_pp": (
                    float(selection_summary_np.mean(cv)) if len(cv) >= 6 else None
                ),
                "original_candidate_minus_comparison_pp": (
                    float(selection_summary_np.mean(original)) if len(original) >= 6 else None
                ),
            }
        )
    _regional.regional_table(selection_summary_F / "results/F1_specimen_scores.tsv", scores)
    _regional.regional_table(selection_summary_F / "results/F1_pair_scores.tsv", pairrows)
    stability = [
        {
            "locus_id": p["locus_id"],
            "gene": p["gene"],
            "candidate_folds": freq[p["locus_id"]],
            "original_candidate": p["original_locus_id"].startswith("C"),
        }
        for p in panel
        if freq[p["locus_id"]] > 0
    ]
    _regional.regional_table(selection_summary_F / "results/F1_selection_stability.tsv", stability)
    blind = selection_summary_json.loads(
        (selection_summary_F / "plans/blind_panel.json").read_text()
    )
    common = [k for k in blind if all((lookup[l, k]["eligible_locus"] == "True" for l in labels))]
    background = []
    for l in labels:
        background.append(
            {
                "study_label": l,
                "diagnosis": lookup[l, panel[0]["locus_id"]]["diagnosis"],
                "common_blind_promoters": len(common),
                "mean_joint_5hmC_excess_pp": (
                    float(
                        selection_summary_np.mean(
                            [float(lookup[l, k]["mean_joint_5hmC_excess_pp"]) for k in common]
                        )
                    )
                    if common
                    else None
                ),
            }
        )
    _regional.regional_table(selection_summary_F / "results/F1_blind_panel_scores.tsv", background)
    (fig, axes) = selection_summary_plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, diag in zip(axes, selection_summary_COLORS):
        rr = [r for r in scores if r["diagnosis"] == diag and r["eligible_heldout_score"]]
        for r in rr:
            ax.plot(
                [0, 1],
                [
                    r["original_candidate_minus_comparison_pp"],
                    r["heldout_candidate_minus_comparison_pp"],
                ],
                "-o",
                c=selection_summary_COLORS[diag],
                alpha=0.7,
                ms=4,
            )
        ax.axhline(0, c="gray", lw=0.7)
        ax.set_xticks([0, 1], ["Original selected panel", "Held-out selection"])
        ax.set_ylabel("Candidate minus comparison co-occurrence (pp)")
        ax.set_title(diag + f" (n={len(rr)})")
    fig.suptitle(
        "Patient-held-out selection: each patient excluded from selecting their own promoters\nInternal evaluation; the pre-existing common-coverage universe is fixed",
        fontsize=10,
    )
    fig.tight_layout()
    blocks_figure(fig, "F1_01_heldout_molecular_scores")
    top = sorted(stability, key=lambda r: (-r["candidate_folds"], r["locus_id"]))[:30]
    (fig, ax) = selection_summary_plt.subplots(figsize=(9, 7))
    ax.barh(
        range(len(top)),
        [r["candidate_folds"] for r in top],
        color=["#df9d29" if r["original_candidate"] else "#497ba6" for r in top],
    )
    ax.set_yticks(range(len(top)), [r["gene"] for r in top], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 20)
    ax.set_xlabel("Training folds selecting this candidate (out of 20)")
    ax.set_title("Selection stability; gold marks original candidate promoters")
    fig.tight_layout()
    blocks_figure(fig, "F1_02_selection_stability")
    (fig, axes) = selection_summary_plt.subplots(1, 2, figsize=(12, 5))
    mat = selection_summary_np.array(
        [[int(lookup[l, k]["eligible_blocks"]) for l in labels] for k in blind]
    )
    im = axes[0].imshow(mat, aspect="auto", cmap="viridis", vmin=0, vmax=12)
    axes[0].set_xticks(range(20), labels, rotation=90, fontsize=7)
    axes[0].set_ylabel("Diagnosis-blind promoters")
    axes[0].set_title("Coverage across 120 frozen promoters")
    fig.colorbar(im, ax=axes[0], label="Eligible 3-CpG blocks")
    for i, diag in enumerate(selection_summary_COLORS):
        rr = [r for r in background if r["diagnosis"] == diag]
        v = [r["mean_joint_5hmC_excess_pp"] for r in rr]
        axes[1].scatter(
            i + selection_summary_np.linspace(-0.1, 0.1, len(v)),
            v,
            c=selection_summary_COLORS[diag],
        )
        axes[1].plot([i - 0.14, i + 0.14], [selection_summary_np.mean(v)] * 2, c="black")
    axes[1].set_xticks([0, 1], ["Glioblastoma", "Meningioma"])
    axes[1].set_ylabel("Mean excess molecular 5hmC co-occurrence (pp)")
    axes[1].set_title(f"Same {len(common)} eligible blind promoters in every patient")
    fig.tight_layout()
    blocks_figure(fig, "F1_03_diagnosis_blind_panel")
    summary = {
        "heldout_groups": 20,
        "unique_candidate_promoters_across_folds": len(stability),
        "blind_promoters": 120,
        "blind_promoters_eligible_in_every_specimen": len(common),
        "groups": {},
        "all_full_extraction_count_checks_passed": True,
    }
    for diag in selection_summary_COLORS:
        rr = [r for r in scores if r["diagnosis"] == diag and r["eligible_heldout_score"]]
        v = selection_summary_np.array([r["heldout_candidate_minus_comparison_pp"] for r in rr])
        summary["groups"][diag] = {
            "eligible_patients": len(rr),
            "positive_heldout_differences": int((v > 0).sum()),
            "mean_heldout_difference_pp": float(v.mean()),
            "mean_original_difference_pp": float(
                selection_summary_np.mean([r["original_candidate_minus_comparison_pp"] for r in rr])
            ),
            "mean_blind_promoter_cooccurrence_pp": float(
                selection_summary_np.mean(
                    [r["mean_joint_5hmC_excess_pp"] for r in background if r["diagnosis"] == diag]
                )
            ),
        }
    _regional.regional_save(selection_summary_F / "results/F1_summary.json", summary)
    text = [
        "# F1: selection-aware molecular evaluation",
        "",
        "## Results",
        "",
        selection_summary_json.dumps(summary, indent=2),
        "",
        "## Interpretation",
        "",
        "Each patient was excluded from computing all candidate effects, stability criteria, training-only comparison means and matching for their own evaluation. Twenty distinct patient groups are recorded in the existing internal linkage metadata. The scientific leakage test changed held-out modification values drastically and confirmed unchanged training selection. The broader 120-promoter panel was chosen by a fixed random seed using coverage eligibility and genomic separation, without diagnosis contrasts.",
        "",
        "This is a retrospective internal analysis after observing the exploratory finding. It reduces reuse of an individual patient in selection; it does not erase prior hypothesis selection or provide an independent cohort. The shared CpG/feature coverage universe was fixed from all specimens before this follow-up, so held-out coverage participates in eligibility. Overlapping training folds also make scores statistically dependent. No naive independent-fold significance test or bootstrap confidence interval is reported. Read-level counts match the full extraction at every promoter/specimen.",
        "",
        "The diagnosis-blind panel evaluates breadth of co-occurrence, not a second set of opposing-change promoters or a direct replication of the selected-panel contrast. Its abundance dependence is addressed in F2.",
        "",
        "## Figures",
        "",
        "F1_01_heldout_molecular_scores; F1_02_selection_stability; F1_03_diagnosis_blind_panel.",
    ]
    (selection_summary_F / "results/F1_ANSWER.md").write_text("\n".join(text) + "\n")
    _regional.regional_save(
        selection_summary_F / "results/F1_complete.json",
        {"analysis_and_figures_generated": True, "requires_review": True},
    )
    print("F1 COMPLETE", selection_summary_json.dumps(summary), flush=True)


def initialize_selection_summary():
    """Initialize the selection_summary stage once; load its declared inputs."""
    global blocks_S
    global selection_summary_COLORS, selection_summary_D, selection_summary_F, selection_summary_R
    if _runtime.initialized("followup_q1"):
        return
    _runtime.begin("followup_q1")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    selection_summary_R = selection_summary_Path(_config.workspace)
    selection_summary_D = selection_summary_R / ".analysis"
    selection_summary_F = selection_summary_D / "focused_followup"
    blocks_S = selection_summary_F
    selection_summary_COLORS = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    _runtime.finish("followup_q1")


def run_selection_summary():
    """Execute the selection_summary workflow stage."""
    initialize_selection_summary()
    selection_summary_run()
