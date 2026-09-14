"""Regional analysis stages for the brain tumour 5mC/5hmC study.

Each stage has explicit initialization and execution; importing is read-only.
See docs/workflow.md for inputs, stage order and reproduction limits.
"""

from . import runtime as _runtime

from .runtime import config as _config


# REGIONAL


import argparse as regional_argparse, csv as regional_csv, fcntl as regional_fcntl, gzip as regional_gzip, hashlib as regional_hashlib, json as regional_json, re as regional_re, time as regional_time
from pathlib import Path as regional_Path
import numpy as regional_np
import matplotlib as regional_matplotlib
import matplotlib.pyplot as regional_plt
from matplotlib.text import Text as regional_Text
from scipy.stats import spearmanr as regional_spearmanr


def regional_log(s):
    print(regional_time.strftime("%Y-%m-%d %H:%M:%S"), s, flush=True)


def regional_save(path, obj):
    p = regional_Path(path).resolve()
    t = p.with_suffix(p.suffix + ".tmp")
    t.write_text(regional_json.dumps(obj, indent=2) + "\n")
    t.replace(p)


def regional_table(path, rows):
    rows = list(rows)
    if not rows:
        return
    op = regional_gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as f:
        w = regional_csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def regional_figure(fig, name):
    for text in fig.findobj(match=regional_Text):
        if regional_re.search("\\bN(?:20)?\\d{2}[._-]\\d+\\b", text.get_text(), regional_re.I):
            raise ValueError("Clinical accession in figure")
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(regional_O / "figures" / (name + "." + ext), dpi=180, bbox_inches="tight")
    regional_plt.close(fig)


def regional_promoter(start, end, strand, length):
    return (
        (max(0, start - 2000), min(length, start + 500))
        if strand == "+"
        else (max(0, end - 500), min(length, end + 2000))
    )


def regional_features_from_snapshot(lengths):
    groups = {}
    for chrom in range(1, 23):
        with regional_gzip.open(
            regional_D / "reference/ucsc_refseq_snapshot" / f"chr{chrom}.json.gz", "rt"
        ) as f:
            records = regional_json.load(f)["ncbiRefSeq"]
        for r in records:
            if not r["name"].startswith(("NM_", "NR_")):
                continue
            assert (
                r["chrom"] == f"chr{chrom}" and 0 <= r["txStart"] < r["txEnd"] <= lengths[chrom - 1]
            )
            groups.setdefault((chrom, r["name2"], r["strand"]), []).append(r)
    out = []
    for (chrom, gene, strand), rs in sorted(groups.items()):
        rs.sort(key=lambda r: (r["txStart"], r["txEnd"]))
        clusters = []
        for r in rs:
            if not clusters or r["txStart"] > clusters[-1]["end"]:
                clusters.append(
                    {"start": r["txStart"], "end": r["txEnd"], "transcripts": [r["name"]]}
                )
            else:
                clusters[-1]["end"] = max(clusters[-1]["end"], r["txEnd"])
                clusters[-1]["transcripts"].append(r["name"])
        for c in clusters:
            gid = f"{gene}|chr{chrom}:{c['start']}-{c['end']}:{strand}"
            for kind, lo, hi, n in [
                ("gene_body", c["start"], c["end"], 20),
                (
                    "promoter",
                    *regional_promoter(c["start"], c["end"], strand, int(lengths[chrom - 1])),
                    10,
                ),
            ]:
                out.append(
                    {
                        "feature_id": gid + "|" + kind,
                        "gene": gene,
                        "chromosome": f"chr{chrom}",
                        "start": lo,
                        "end": hi,
                        "strand": strand,
                        "feature": kind,
                        "minimum_CpGs": n,
                        "transcript_count": len(c["transcripts"]),
                    }
                )
    return out


def regional_aggregate(v, rid, nw, left, right):
    """Average CpG percentages equally within windows and annotated features.

    nw contains the number of common CpGs in each window; left/right index
    the sorted common-CpG rows for each feature. Coverage does not weight CpGs."""
    win = regional_np.column_stack(
        [
            regional_np.divide(
                regional_np.bincount(rid, weights=v[:, k], minlength=len(nw)),
                nw,
                out=regional_np.full(len(nw), regional_np.nan),
                where=nw > 0,
            )
            for k in range(2)
        ]
    )
    feats = []
    for k in range(2):
        cs = regional_np.concatenate(
            ([0.0], regional_np.cumsum(v[:, k], dtype=regional_np.float64))
        )
        feats.append((cs[right] - cs[left]) / (right - left))
    return (win, regional_np.column_stack(feats))


def regional_prepare_context(depth, idx, sites, off, lengths, regions, all_features):
    regional_np.save(regional_O / "cache" / f"idx{depth}.npy", idx)
    rid = regional_np.empty(len(idx), regional_np.int32)
    rbase = 0
    for chrom in range(22):
        (lo, hi) = regional_np.searchsorted(idx, [off[chrom], off[chrom + 1]])
        rid[lo:hi] = sites[idx[lo:hi]] // 100000 + rbase
        rbase += int((int(lengths[chrom]) + 99999) // 100000)
    nw = regional_np.bincount(rid, minlength=len(regions))
    features = []
    left = []
    right = []
    for f in all_features:
        ch = int(f["chromosome"][3:]) - 1
        cg = sites[int(off[ch]) : int(off[ch + 1])]
        (a, b) = regional_np.searchsorted(cg, [f["start"], f["end"]]) + off[ch]
        (l, r) = regional_np.searchsorted(idx, [a, b])
        if r - l >= f["minimum_CpGs"]:
            features.append({**f, "common_CpGs": int(r - l)})
            left.append(l)
            right.append(r)
    regional_np.savez(
        regional_O / "cache" / f"context{depth}.npz",
        rid=rid,
        nw=nw,
        left=regional_np.array(left, dtype=regional_np.int64),
        right=regional_np.array(right, dtype=regional_np.int64),
    )
    regional_save(regional_O / "cache" / f"features{depth}.json", features)
    return (nw, features)


def regional_specimen(task):
    (label, depth) = task
    idx = regional_np.load(regional_O / "cache" / f"idx{depth}.npy", mmap_mode="r")
    ctx = regional_np.load(regional_O / "cache" / f"context{depth}.npz")
    a = regional_np.load(regional_D / "results" / label / "site_counts.npz")["counts"]
    if depth == 5:
        counts = a[:, :3].sum(1, dtype=regional_np.uint64)
        assert regional_np.all(counts[idx] >= 5)
        regional_np.save(regional_O / "cache" / (label + "_depth10.npy"), counts >= 10)
        del counts
    c = a[idx, :3].astype(regional_np.float64)
    del a
    den = c.sum(1)
    assert regional_np.all(den >= depth)
    v = c[:, 1:3] / den[:, None] * 100
    del c, den
    (win, feat) = regional_aggregate(v, ctx["rid"], ctx["nw"], ctx["left"], ctx["right"])
    global_mean = v.mean(0)
    regional_np.savez_compressed(
        regional_O / "cache" / f"{label}_depth{depth}.npz",
        windows=win,
        features=feat,
        global_mean=global_mean,
    )
    regional_log(f"{label}: common-CpG profiles complete at depth {depth}")
    return label


def regional_effects(beta, group):
    """Compute GBM-minus-meningioma contrasts with equal specimen weights.

    The common-CpG set is fixed when each specimen is omitted. Returns mean
    contrasts, sign stability, minimum omitted-sample effect, medians and
    cross-diagnosis pairwise agreement."""
    delta = beta[group].mean(0) - beta[~group].mean(0)
    loo = []
    for i in range(len(group)):
        keep = regional_np.arange(len(group)) != i
        loo.append(beta[keep & group].mean(0) - beta[keep & ~group].mean(0))
    loo = regional_np.asarray(loo)
    stable = regional_np.all(
        regional_np.sign(loo[:, :, 1]) == regional_np.sign(delta[:, 1]), axis=0
    )
    median = regional_np.median(beta[group], axis=0) - regional_np.median(beta[~group], axis=0)
    pairwise = regional_np.mean(
        regional_np.sign(beta[group, :, 1][:, None, :] - beta[~group, :, 1][None, :, :])
        == regional_np.sign(delta[:, 1]),
        axis=(0, 1),
    )
    return (delta, stable, regional_np.min(abs(loo[:, :, 1]), axis=0), median, pairwise)


def regional_bootstrap_intervals(beta, group, seed=20260912, replicates=2000):
    """Return pointwise 95% intervals from within-diagnosis specimen resampling.

    The seeded multinomial weights reproduce explicit bootstrap sample means.
    Intervals are descriptive and are not adjusted for multiple comparisons."""
    rng = regional_np.random.default_rng(seed)
    n = len(group)
    weights = regional_np.zeros((replicates, n))
    gi = regional_np.flatnonzero(group)
    mi = regional_np.flatnonzero(~group)
    weights[:, gi] = rng.multinomial(
        len(gi), regional_np.ones(len(gi)) / len(gi), size=replicates
    ) / len(gi)
    weights[:, mi] = -rng.multinomial(
        len(mi), regional_np.ones(len(mi)) / len(mi), size=replicates
    ) / len(mi)
    intervals = regional_np.empty((beta.shape[1], 2))
    for start in range(0, beta.shape[1], 1000):
        stop = min(start + 1000, beta.shape[1])
        dist = weights @ beta[:, start:stop, 1]
        intervals[start:stop] = regional_np.quantile(dist, [0.025, 0.975], axis=0).T
    return intervals


def regional_result_rows(info, beta, group, intervals=None):
    (delta, stable, minimum, median, pair) = regional_effects(beta, group)
    rows = []
    for i, f in enumerate(info):
        h = float(delta[i, 1])
        m = float(delta[i, 0])
        candidate = abs(h) >= 2 and bool(stable[i])
        cancel = candidate and m * h < 0 and (abs(m + h) <= 1)
        row = {
            **f,
            "GBM_minus_meningioma_5mC_pp": m,
            "GBM_minus_meningioma_5hmC_pp": h,
            "GBM_minus_meningioma_combined_pp": m + h,
            "5hmC_direction_stable_all_single_omissions": bool(stable[i]),
            "minimum_abs_leave_one_out_5hmC_pp": float(minimum[i]),
            "median_group_5hmC_difference_pp": float(median[i, 1]),
            "pairwise_comparisons_agreeing_with_mean_direction": float(pair[i]),
            "exploratory_5hmC_candidate": candidate,
            "opposing_changes_with_small_combined_difference": cancel,
        }
        if intervals is not None:
            row.update(
                bootstrap_5hmC_lower_pp=float(intervals[i, 0]),
                bootstrap_5hmC_upper_pp=float(intervals[i, 1]),
            )
        rows.append(row)
    return rows


def regional_run(workers):
    lock = (regional_O / "analysis.lock").open("w")
    regional_fcntl.flock(lock, regional_fcntl.LOCK_EX | regional_fcntl.LOCK_NB)
    samples = sorted(
        regional_json.loads((regional_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    labels = [s["study_label"] for s in samples]
    group = regional_np.array([s["diagnosis"] == "Glioblastoma" for s in samples])
    assert len(samples) == 20 and group.sum() == 15
    assert all(((regional_D / "results" / l / "completion.json").exists() for l in labels))
    sites = regional_np.load(regional_D / "reference/cpg_positions.npy", mmap_mode="r")
    off = regional_np.load(regional_D / "reference/cpg_offsets.npy")
    lengths = regional_np.load(regional_D / "reference/chromosome_lengths.npy")
    regions = []
    for chrom, length in enumerate(lengths, 1):
        for start in range(0, int(length), 100000):
            regions.append(
                {
                    "region": f"chr{chrom}:{start}-{min(start + 100000, int(length))}",
                    "chromosome": f"chr{chrom}",
                    "start": start,
                    "end": min(start + 100000, int(length)),
                }
            )
    common = regional_np.ones(len(sites), bool)
    for label in labels:
        counts = regional_np.load(regional_D / "results" / label / "site_counts.npz")["counts"]
        common &= counts[:, :3].sum(1, dtype=regional_np.uint64) >= 5
        del counts
    idx = regional_np.flatnonzero(common)
    regional_np.save(regional_O / "common_cpg_mask_depth5.npy", common)
    del common
    features = regional_features_from_snapshot(lengths)
    regional_table(regional_O / "tables" / "annotation_features.tsv.gz", features)
    regional_log(f"Pooled common CpGs: {len(idx):,}; annotation features: {len(features):,}")
    (n5, f5) = regional_prepare_context(5, idx, sites, off, lengths, regions, features)
    with _runtime.process_pool(max_workers=workers) as pool:
        list(pool.map(regional_specimen, [(l, 5) for l in labels]))
    common10 = regional_np.ones(len(sites), bool)
    for l in labels:
        common10 &= regional_np.load(regional_O / "cache" / (l + "_depth10.npy"), mmap_mode="r")
    idx10 = regional_np.flatnonzero(common10)
    assert regional_np.all(regional_np.isin(idx10, idx, assume_unique=True))
    regional_np.save(regional_O / "common_cpg_mask_depth10.npy", common10)
    del common10
    (n10, f10) = regional_prepare_context(10, idx10, sites, off, lengths, regions, features)
    with _runtime.process_pool(max_workers=workers) as pool:
        list(pool.map(regional_specimen, [(l, 10) for l in labels]))
    loaded5 = [regional_np.load(regional_O / "cache" / f"{l}_depth5.npz") for l in labels]
    loaded10 = [regional_np.load(regional_O / "cache" / f"{l}_depth10.npz") for l in labels]
    b5 = regional_np.array([a["windows"] for a in loaded5])
    b10 = regional_np.array([a["windows"] for a in loaded10])
    g5 = regional_np.array([a["features"] for a in loaded5])
    g10 = regional_np.array([a["features"] for a in loaded10])
    global5 = regional_np.array([a["global_mean"] for a in loaded5])
    valid = n5 >= 50
    valid10 = n10 >= 50
    good = regional_np.flatnonzero(valid)
    regional_np.savez_compressed(
        regional_O / "pooled_profiles.npz",
        labels=regional_np.array(labels),
        diagnoses=regional_np.array([s["diagnosis"] for s in samples]),
        windows5=b5,
        windows10=b10,
        features5=g5,
        features10=g10,
        valid5=valid,
        valid10=valid10,
        common_CpGs5=n5,
        common_CpGs10=n10,
        global5=global5,
    )
    regional_log("Calculating specimen bootstrap intervals and omission sensitivity")
    rows = regional_result_rows(
        [{**regions[j], "common_CpGs": int(n5[j])} for j in good],
        b5[:, good],
        group,
        regional_bootstrap_intervals(b5[:, good], group),
    )
    d10 = b10[group].mean(0) - b10[~group].mean(0)
    for j, row in zip(good, rows):
        row.update(
            depth10_eligible=bool(valid10[j]),
            depth10_5hmC_difference_pp=float(d10[j, 1]) if valid10[j] else "",
            depth10_same_5hmC_direction=(
                bool(
                    regional_np.sign(d10[j, 1])
                    == regional_np.sign(row["GBM_minus_meningioma_5hmC_pp"])
                )
                if valid10[j]
                else ""
            ),
        )
    regional_table(regional_O / "tables" / "pooled_100kb_regions.tsv.gz", rows)
    genes = regional_result_rows(f5, g5, group, regional_bootstrap_intervals(g5, group))
    f10idx = {x["feature_id"]: i for (i, x) in enumerate(f10)}
    gd10 = g10[group].mean(0) - g10[~group].mean(0)
    for row in genes:
        k = f10idx.get(row["feature_id"])
        row.update(
            depth10_eligible=k is not None,
            depth10_5hmC_difference_pp=float(gd10[k, 1]) if k is not None else "",
            depth10_same_5hmC_direction=(
                bool(
                    regional_np.sign(gd10[k, 1])
                    == regional_np.sign(row["GBM_minus_meningioma_5hmC_pp"])
                )
                if k is not None
                else ""
            ),
        )
    regional_table(regional_O / "tables" / "pooled_gene_promoter_results.tsv.gz", genes)
    for name, rs in [("window", rows), ("gene_promoter", genes)]:
        candidates = sorted(
            [r for r in rs if r["exploratory_5hmC_candidate"]],
            key=lambda r: -abs(r["GBM_minus_meningioma_5hmC_pp"]),
        )
        cancel = [r for r in candidates if r["opposing_changes_with_small_combined_difference"]]
        regional_table(regional_O / "tables" / f"{name}_candidates.tsv.gz", candidates)
        regional_table(regional_O / "tables" / f"{name}_opposing_changes.tsv.gz", cancel)
    qc = []
    for i, s in enumerate(samples):
        q = regional_json.loads((regional_D / "results" / s["study_label"] / "qc.json").read_text())
        t = q["site_totals"]
        qc.append(
            {
                "study_label": s["study_label"],
                "diagnosis": s["diagnosis"],
                "grade": s.get("grade", ""),
                "aligned_genome_equivalent": q["aligned_genome_equivalent"],
                "common_CpG_5mC_percent": float(global5[i, 0]),
                "common_CpG_5hmC_percent": float(global5[i, 1]),
                "common_CpG_combined_percent": float(global5[i].sum()),
                "missing_intact_read_CpG_fraction": t["missing_at_intact_read_CpG"]
                / max(1, t["intact_read_CpG_opportunities"]),
            }
        )
    regional_table(regional_O / "tables" / "specimen_summary.tsv", qc)
    regional_plot_results(samples, global5, rows, b5[:, good], genes, g5)
    ok = valid & valid10
    d5 = b5[group].mean(0) - b5[~group].mean(0)
    stats = {
        "specimens": 20,
        "glioblastomas": 15,
        "meningiomas": 5,
        "common_CpGs_depth5": len(idx),
        "common_CpGs_depth10": len(idx10),
        "eligible_windows_depth5": int(valid.sum()),
        "eligible_windows_depth10": int(valid10.sum()),
        "eligible_gene_bodies": sum((f["feature"] == "gene_body" for f in f5)),
        "eligible_promoters": sum((f["feature"] == "promoter" for f in f5)),
        "windows_meeting_exploratory_5hmC_rule": sum(
            (r["exploratory_5hmC_candidate"] for r in rows)
        ),
        "windows_with_opposing_changes_and_small_combined_difference": sum(
            (r["opposing_changes_with_small_combined_difference"] for r in rows)
        ),
        "gene_body_candidates": sum(
            (r["exploratory_5hmC_candidate"] and r["feature"] == "gene_body" for r in genes)
        ),
        "promoter_candidates": sum(
            (r["exploratory_5hmC_candidate"] and r["feature"] == "promoter" for r in genes)
        ),
        "gene_promoter_opposing_change_candidates": sum(
            (r["opposing_changes_with_small_combined_difference"] for r in genes)
        ),
        "depth5_depth10_5hmC_effect_spearman": float(
            regional_spearmanr(d5[ok, 1], d10[ok, 1]).statistic
        ),
        "global_common_CpG_means": {
            diag: {
                "5mC": float(global5[group if diag == "Glioblastoma" else ~group, 0].mean()),
                "5hmC": float(global5[group if diag == "Glioblastoma" else ~group, 1].mean()),
            }
            for diag in ["Glioblastoma", "Meningioma"]
        },
    }
    regional_save(regional_O / "summary.json", stats)
    regional_write_report(stats, rows, genes)
    regional_save(
        regional_O / "analysis_complete.json",
        {
            "completed_utc": regional_time.strftime("%Y-%m-%dT%H:%M:%SZ", regional_time.gmtime()),
            "pooled_windows": True,
            "gene_promoter_analysis": True,
            "coverage_sensitivity": True,
            "script_sha256": regional_hashlib.sha256(
                regional_Path(__file__).read_bytes()
            ).hexdigest(),
        },
    )
    regional_log("POOLED ANALYSIS COMPLETE " + regional_json.dumps(stats))


def regional_plot_results(samples, global5, rows, beta, genes, gene_beta):
    group = regional_np.array([s["diagnosis"] == "Glioblastoma" for s in samples])
    labels = [s["study_label"] for s in samples]
    rng = regional_np.random.default_rng(20260912)
    (fig, axes) = regional_plt.subplots(1, 3, figsize=(11, 4))
    for k, (ax, title) in enumerate(zip(axes, ["5mC", "5hmC", "5mC + 5hmC"])):
        vals = global5[:, k] if k < 2 else global5.sum(1)
        for gi, diag in enumerate(["Glioblastoma", "Meningioma"]):
            mask = group if gi == 0 else ~group
            v = vals[mask]
            ax.scatter(
                gi + rng.uniform(-0.13, 0.13, len(v)), v, c=regional_COLORS[diag], s=36, alpha=0.85
            )
            ax.plot([gi - 0.2, gi + 0.2], [v.mean()] * 2, c="black", lw=2)
        ax.set_xticks([0, 1], ["Glioblastoma\n(n=15)", "Meningioma\n(n=5)"])
        ax.set_title(title)
        ax.set_ylabel("Mean fraction at common CpGs (%)")
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Genome-wide modification fractions: one point per specimen")
    fig.tight_layout()
    regional_figure(fig, "01_pooled_modification_fractions")
    m = regional_np.array([r["GBM_minus_meningioma_5mC_pp"] for r in rows])
    h = regional_np.array([r["GBM_minus_meningioma_5hmC_pp"] for r in rows])
    sel = regional_np.array([r["opposing_changes_with_small_combined_difference"] for r in rows])
    (fig, axes) = regional_plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, x, xlabel in [
        (axes[0], m, "5mC difference (percentage points)"),
        (axes[1], m + h, "Combined difference (percentage points)"),
    ]:
        ax.scatter(x, h, s=3, alpha=0.13, c="#667785", rasterized=True)
        ax.scatter(
            x[sel],
            h[sel],
            s=7,
            c="#D79A19",
            alpha=0.65,
            label="Opposing-change candidates",
            rasterized=True,
        )
        ax.axhline(0, c="gray", lw=0.6)
        ax.axvline(0, c="gray", lw=0.6)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("5hmC difference (percentage points)")
    axes[1].legend(fontsize=8)
    fig.suptitle("100 kb windows: glioblastoma minus meningioma")
    fig.tight_layout()
    regional_figure(fig, "02_separate_and_combined_modifications")
    candidates = [i for (i, r) in enumerate(rows) if r["exploratory_5hmC_candidate"]]
    selected = sorted(candidates, key=lambda i: -abs(h[i]))[:40]
    if selected:
        mat = beta[:, selected, 1].T
        mat -= mat.mean(1)[:, None]
        lim = max(1, float(abs(mat).max()))
        (fig, ax) = regional_plt.subplots(figsize=(11, 8))
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim)
        ax.set_xticks(range(len(labels)), labels, rotation=60)
        ax.set_yticks(range(len(selected)), [rows[i]["region"] for i in selected], fontsize=6)
        ax.axvline(14.5, color="black", lw=1)
        ax.set_title(
            "Largest 5hmC contrasts stable to every single-specimen omission\nSelected using these same specimens; exploratory"
        )
        fig.colorbar(im, ax=ax, label="5hmC deviation from regional mean (percentage points)")
        fig.tight_layout()
        regional_figure(fig, "03_candidate_region_heatmap")
    choices = []
    seen = set()
    for i in sorted(
        range(len(genes)), key=lambda j: -abs(genes[j]["GBM_minus_meningioma_5hmC_pp"])
    ):
        r = genes[i]
        if (
            r["exploratory_5hmC_candidate"]
            and r["gene"] not in seen
            and (r["depth10_same_5hmC_direction"] is True)
        ):
            choices.append(i)
            seen.add(r["gene"])
        if len(choices) == 6:
            break
    if choices:
        (fig, axes) = regional_plt.subplots(2, 3, figsize=(12, 7))
        for ax, i in zip(axes.ravel(), choices):
            r = genes[i]
            for gi, diag in enumerate(["Glioblastoma", "Meningioma"]):
                mask = group if gi == 0 else ~group
                v = gene_beta[mask, i, 1]
                ax.scatter(
                    gi + rng.uniform(-0.12, 0.12, len(v)), v, color=regional_COLORS[diag], s=25
                )
                ax.plot([gi - 0.2, gi + 0.2], [v.mean()] * 2, color="black")
            ax.set_xticks([0, 1], ["GBM", "Meningioma"])
            ax.set_title(r["gene"] + " / " + r["feature"].replace("_", " "))
            ax.set_ylabel("Mean 5hmC (%)")
            ax.grid(axis="y", alpha=0.2)
        for ax in axes.ravel()[len(choices) :]:
            ax.set_visible(False)
        fig.suptitle("Illustrative strongest gene-associated contrasts; selected from this cohort")
        fig.tight_layout()
        regional_figure(fig, "04_gene_examples")


def regional_write_report(s, rows, genes):
    lines = [
        "# Pooled 5hmC analysis",
        "",
        "All 20 specimens were analysed together by confirmed diagnosis: 15 glioblastomas and 5 meningiomas. No software-version grouping or adjustment was performed.",
        "",
        "## Results",
        "",
    ]
    lines += [
        f"- Common CpGs with at least five passing calls in every specimen: {s['common_CpGs_depth5']:,}.",
        f"- Eligible 100 kb windows: {s['eligible_windows_depth5']:,}.",
        f"- Eligible gene bodies/promoters: {s['eligible_gene_bodies']:,} / {s['eligible_promoters']:,}.",
        f"- Windows with absolute mean 5hmC contrast >=2 percentage points, retaining direction after every single-specimen omission: {s['windows_meeting_exploratory_5hmC_rule']:,}.",
        f"- Of these, opposing 5mC/5hmC changes with absolute combined contrast <=1 percentage point: {s['windows_with_opposing_changes_and_small_combined_difference']:,}.",
        f"- Correlation of window 5hmC contrasts at the five-call versus ten-call CpG threshold: {s['depth5_depth10_5hmC_effect_spearman']:.3f}.",
        "",
        "## Methods and interpretation",
        "",
        "Identical CpGs are used in every specimen; each CpG contributes equally within a window or feature. Group differences are means of specimen values (GBM minus meningioma). Leave-one-specimen-out checks keep this common CpG set fixed. Bootstrap intervals use 2,000 resamples within diagnosis, seed 20260912. They are pointwise intervals and are not adjusted for multiple comparisons. No significance tests or classifier were fitted. The 2/1 percentage-point candidate cutoffs are exploratory choices, not validated biological thresholds. A small combined mean difference does not establish statistical equivalence.",
        "",
        "Curated NM/NR RefSeq transcripts were retrieved from the [UCSC hg38 API](https://genome.ucsc.edu/goldenPath/help/api.html); all 22 autosomal snapshots, source timestamps and checksums are stored under ../reference/ucsc_refseq_snapshot/. Overlapping transcripts with the same gene symbol, chromosome and strand are merged into a locus. Gene bodies span that locus; promoters use its upstream-most TSS, -2,000/+500 bp. Separate disjoint loci remain separate. These are annotation-defined regions, not evidence of gene expression or regulatory activity.",
        "",
        "Patients/specimens, rather than reads or CpGs, provide biological replication. Five meningiomas limit precision; the study includes 20 distinct recorded patients. Tissue composition, tumour lineage, copy number and unmeasured technical effects may contribute. Normal controls and orthogonal modification measurements are unavailable. Selected example plots and heatmaps use this same cohort and are not independent validation.",
        "",
        "## Files",
        "",
        "- tables/specimen_summary.tsv: study-labelled specimen measurements.",
        "- tables/pooled_100kb_regions.tsv.gz: all eligible windows, contrasts, intervals and robustness.",
        "- tables/pooled_gene_promoter_results.tsv.gz: gene-body/promoter results.",
        "- tables/*candidates.tsv.gz and *opposing_changes.tsv.gz: exploratory shortlists.",
        "- figures/: PNG, PDF and SVG plots, all without clinical accession numbers.",
        "- pooled_profiles.npz: specimen-level arrays for follow-up analyses.",
        "",
        "## Next interpretation step",
        "",
        "Review candidate loci and their genomic context, check whether neighbouring windows represent the same broad locus, and assess whether the opposing-change examples support the proposed paper question. Do not infer pathway activation, diagnostic accuracy or causality from these exploratory contrasts.",
    ]
    (regional_O / "POOLED_REPORT.md").write_text("\n".join(lines) + "\n")


def initialize_regional():
    """Initialize the regional stage once; load its declared inputs."""
    global regional_COLORS, regional_D, regional_O, regional_R
    if _runtime.initialized("pooled_analysis"):
        return
    _runtime.begin("pooled_analysis")
    regional_matplotlib.use("Agg")
    regional_R = regional_Path(_config.workspace)
    regional_D = regional_R / ".analysis"
    regional_O = regional_D / "pooled_analysis"
    regional_COLORS = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    _runtime.finish("pooled_analysis")


def run_regional():
    """Execute the regional workflow stage."""
    initialize_regional()
    global regional_a, regional_p
    regional_p = regional_argparse.ArgumentParser()
    regional_p.add_argument("--workers", type=int, default=6)
    regional_a = regional_p.parse_args()
    regional_run(regional_a.workers)


# INTERPRET


import csv as interpret_csv, gzip as interpret_gzip, json as interpret_json, numpy as interpret_np


def interpret_stronger(r):
    h = float(r["GBM_minus_meningioma_5hmC_pp"])
    return (
        r["depth10_same_5hmC_direction"] == "True"
        and float(r["pairwise_comparisons_agreeing_with_mean_direction"]) >= 0.8
        and (h * float(r["median_group_5hmC_difference_pp"]) > 0)
        and (float(r["bootstrap_5hmC_lower_pp"]) * float(r["bootstrap_5hmC_upper_pp"]) > 0)
    )


def initialize_interpret():
    """Initialize the interpret stage once; load its declared inputs."""
    global interpret_O, interpret_ax, interpret_axes, interpret_beta, interpret_blocks, interpret_chosen, interpret_col, interpret_delta, interpret_diag, interpret_end, interpret_f, interpret_features, interpret_fi, interpret_fig, interpret_genes, interpret_gi, interpret_group, interpret_kind, interpret_lines, interpret_mask, interpret_merged, interpret_profiles, interpret_r, interpret_rng, interpret_row, interpret_s, interpret_seen, interpret_sign, interpret_start, interpret_state, interpret_strong, interpret_v, interpret_val, interpret_values, interpret_windows, interpret_y
    if _runtime.initialized("interpret_pooled"):
        return
    _runtime.begin("interpret_pooled")
    _runtime.initialize("pooled_analysis")
    interpret_O = regional_O
    with interpret_gzip.open(
        interpret_O / "tables/pooled_100kb_regions.tsv.gz", "rt"
    ) as interpret_f:
        interpret_windows = list(interpret_csv.DictReader(interpret_f, delimiter="\t"))
    with interpret_gzip.open(
        interpret_O / "tables/gene_promoter_opposing_changes.tsv.gz", "rt"
    ) as interpret_f:
        interpret_genes = list(interpret_csv.DictReader(interpret_f, delimiter="\t"))
    interpret_strong = sorted(
        [r for r in interpret_genes if interpret_stronger(r)],
        key=lambda r: -abs(float(r["GBM_minus_meningioma_5hmC_pp"])),
    )
    regional_table(
        interpret_O / "tables/opposing_gene_features_additional_consistency.tsv.gz",
        interpret_strong,
    )
    interpret_blocks = {}
    for interpret_kind, interpret_col in [
        ("5hmC_candidates", "exploratory_5hmC_candidate"),
        ("opposing_changes", "opposing_changes_with_small_combined_difference"),
    ]:
        interpret_merged = []
        for interpret_r in interpret_windows:
            if interpret_r[interpret_col] != "True":
                continue
            (interpret_start, interpret_end) = (int(interpret_r["start"]), int(interpret_r["end"]))
            interpret_sign = int(
                interpret_np.sign(float(interpret_r["GBM_minus_meningioma_5hmC_pp"]))
            )
            if (
                interpret_merged
                and interpret_merged[-1]["chromosome"] == interpret_r["chromosome"]
                and (interpret_merged[-1]["end"] == interpret_start)
                and (interpret_merged[-1]["direction"] == interpret_sign)
            ):
                interpret_merged[-1]["end"] = interpret_end
                interpret_merged[-1]["window_count"] += 1
            else:
                interpret_merged.append(
                    {
                        "chromosome": interpret_r["chromosome"],
                        "start": interpret_start,
                        "end": interpret_end,
                        "direction": interpret_sign,
                        "window_count": 1,
                    }
                )
        regional_table(
            interpret_O / "tables" / (interpret_kind + "_adjacent_blocks.tsv"), interpret_merged
        )
        interpret_blocks[interpret_kind] = len(interpret_merged)
    interpret_features = interpret_json.loads((interpret_O / "cache/features5.json").read_text())
    interpret_fi = {x["feature_id"]: i for (i, x) in enumerate(interpret_features)}
    interpret_profiles = interpret_np.load(interpret_O / "pooled_profiles.npz")
    interpret_group = interpret_profiles["diagnoses"] == "Glioblastoma"
    interpret_beta = interpret_profiles["features5"]
    interpret_chosen = []
    interpret_seen = set()
    for interpret_r in interpret_strong:
        if int(interpret_r["common_CpGs"]) >= 50 and interpret_r["gene"] not in interpret_seen:
            interpret_chosen.append(interpret_r)
            interpret_seen.add(interpret_r["gene"])
        if len(interpret_chosen) == 3:
            break
    (interpret_fig, interpret_axes) = regional_plt.subplots(
        len(interpret_chosen), 3, figsize=(11, 3 * len(interpret_chosen)), squeeze=False
    )
    interpret_rng = interpret_np.random.default_rng(20260912)
    for interpret_row, interpret_r in enumerate(interpret_chosen):
        interpret_v = interpret_beta[:, interpret_fi[interpret_r["feature_id"]], :]
        interpret_values = [interpret_v[:, 0], interpret_v[:, 1], interpret_v.sum(1)]
        for interpret_col, (interpret_state, interpret_val) in enumerate(
            zip(["5mC", "5hmC", "5mC + 5hmC"], interpret_values)
        ):
            interpret_ax = interpret_axes[interpret_row, interpret_col]
            for interpret_gi, interpret_diag in enumerate(["Glioblastoma", "Meningioma"]):
                interpret_mask = interpret_group if interpret_gi == 0 else ~interpret_group
                interpret_y = interpret_val[interpret_mask]
                interpret_ax.scatter(
                    interpret_gi + interpret_rng.uniform(-0.12, 0.12, len(interpret_y)),
                    interpret_y,
                    c=regional_COLORS[interpret_diag],
                    s=23,
                    alpha=0.85,
                )
                interpret_ax.plot(
                    [interpret_gi - 0.2, interpret_gi + 0.2],
                    [interpret_y.mean()] * 2,
                    c="black",
                    lw=2,
                )
            interpret_delta = float(
                interpret_val[interpret_group].mean() - interpret_val[~interpret_group].mean()
            )
            interpret_ax.set_title(f"{interpret_state}: difference {interpret_delta:+.1f} pp")
            interpret_ax.set_xticks([0, 1], ["GBM", "Meningioma"])
            interpret_ax.set_ylabel("Mean fraction at common CpGs (%)")
            interpret_ax.grid(axis="y", alpha=0.2)
        interpret_axes[interpret_row, 0].text(
            -0.34,
            0.5,
            interpret_r["gene"]
            + "\n"
            + interpret_r["feature"].replace("_", " ")
            + "\n"
            + interpret_r["common_CpGs"]
            + " CpGs",
            transform=interpret_axes[interpret_row, 0].transAxes,
            ha="right",
            va="center",
            fontsize=10,
        )
    interpret_fig.suptitle(
        "Opposing 5mC and 5hmC differences with similar combined group means\nExploratory examples selected from this cohort",
        fontsize=12,
    )
    interpret_fig.tight_layout()
    regional_figure(interpret_fig, "05_opposing_changes_examples")
    regional_save(
        interpret_O / "interpretation_summary.json",
        {
            "adjacent_same_direction_blocks": interpret_blocks,
            "gene_promoter_features_with_additional_consistency": len(interpret_strong),
            "additional_consistency_rule": "Depth10 effect direction retained; >=80% GBM-MEN specimen pairwise comparisons agree with mean direction; group median agrees; pointwise bootstrap interval excludes zero. Exploratory, not independent validation or multiple-testing control.",
            "figure_example_rule": "Top three unique genes by absolute effect among additional-consistency candidates with >=50 common CpGs",
            "figure_examples": [
                {
                    k: r[k]
                    for k in [
                        "gene",
                        "feature",
                        "common_CpGs",
                        "GBM_minus_meningioma_5mC_pp",
                        "GBM_minus_meningioma_5hmC_pp",
                        "GBM_minus_meningioma_combined_pp",
                    ]
                }
                for r in interpret_chosen
            ],
        },
    )
    interpret_s = interpret_json.loads((interpret_O / "summary.json").read_text())
    interpret_lines = [
        "# Working paper outline",
        "",
        "## Working title",
        "",
        "Separate 5mC and 5hmC profiling in deeply sequenced glioblastoma and meningioma specimens",
        "",
        "## Main question",
        "",
        "Does separate measurement of 5mC and 5hmC expose regional contrasts that are small in the combined modification measure?",
        "",
        "## Current findings",
        "",
        "- 20 specimens: 15 glioblastomas and 5 meningiomas; all pooled by diagnosis.",
        f"- {interpret_s['common_CpGs_depth5']:,} common CpGs and {interpret_s['eligible_windows_depth5']:,} eligible 100 kb windows.",
        f"- {interpret_s['windows_meeting_exploratory_5hmC_rule']:,} exploratory windows satisfy the effect-size and single-omission direction rule; adjacent same-direction windows form {interpret_blocks['5hmC_candidates']:,} blocks.",
        f"- {interpret_s['windows_with_opposing_changes_and_small_combined_difference']:,} windows show opposing 5mC/5hmC differences and a combined contrast within one percentage point; these form {interpret_blocks['opposing_changes']:,} adjacent blocks.",
        f"- {len(interpret_strong):,} gene/promoter features additionally retain direction at stricter coverage, agree in >=80% of cross-diagnosis specimen pairs, have matching median direction, and have a pointwise 5hmC interval excluding zero.",
        f"- Window effect correlation between coverage thresholds: {interpret_s['depth5_depth10_5hmC_effect_spearman']:.3f}.",
        "",
        "## Figures",
        "",
        "1. Specimen-level common-CpG 5mC, 5hmC and combined fractions.",
        "2. Regional contrasts for separate and combined modifications.",
        "3. Candidate-window heatmap and coverage/omission sensitivity (supplement as appropriate).",
        "4. Gene-associated examples, emphasizing the separate-versus-combined panels in Figure 05.",
        "",
        "## Interpretation limits",
        "",
        "These are descriptive candidates selected from this cohort, not established biomarkers. Adjacent windows, overlapping genes and promoters are not independent findings. Five meningiomas and variation within each diagnosis limit precision. Similar combined group means are not an equivalence test. Bulk lineage/cell composition and copy-number differences may contribute. No gene-expression or pathway-activity conclusions follow from annotation alone.",
        "",
        "## Remaining manuscript work",
        "",
        "Inspect genomic context and annotation boundaries for the chosen examples; relate them to the relevant literature without selecting only familiar genes. Confirm reporting completeness, patient linkage where possible, and missing grade information from existing records. Draft the methods/results and assess journal fit after the biological interpretation. No further full-run BAM extraction is required.",
    ]
    (interpret_O / "MANUSCRIPT_OUTLINE.md").write_text("\n".join(interpret_lines) + "\n")
    print(
        interpret_json.dumps(
            {
                "blocks": interpret_blocks,
                "additional_consistency_features": len(interpret_strong),
                "examples": [(r["gene"], r["feature"]) for r in interpret_chosen],
            }
        )
    )
    _runtime.finish("interpret_pooled")


def run_interpret():
    """Execute the interpret workflow stage."""
    initialize_interpret()
