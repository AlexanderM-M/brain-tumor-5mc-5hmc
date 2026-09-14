"""Association analysis stages for the brain tumour 5mC/5hmC study.

Each stage has explicit initialization and execution; importing is read-only.
See docs/workflow.md for inputs, stage order and reproduction limits.
"""

from . import runtime as _runtime

from .runtime import config as _config

from . import molecules as _molecules
from . import regional as _regional


# DISTANCE


from pathlib import Path as distance_Path
import json as distance_json, gzip as distance_gzip, hashlib as distance_hashlib, time as distance_time, os as distance_os
import numpy as distance_np, pysam as distance_pysam


def distance_prepare():
    plan = {
        "distance_edges_bp": distance_BINS.tolist(),
        "pairs_per_locus_distance_bin_max": 40,
        "pair_selection": "fixed coordinate-seeded random selection before outcome analysis",
        "minimum_complete_molecules": 12,
        "fixed_depth": 12,
        "resamples": 64,
        "primary_covariance": "Exact expected joint-minus-independent excess under random subsampling without replacement to 12 reads: 11/12 times unbiased full-pair sample covariance",
        "phi": "Mean binary phi across 64 hypergeometric samples of 12 reads; >=3 positive and >=3 negative reads for each CpG in full pair; >=32 nondegenerate resamples. A conditional column-permutation null has expected phi zero.",
        "matching": "within specimen and panel, exact distance-bin x sorted two-site frequency bins x coarse immediate-flank-GC bin; overlap weights=min(candidate pair count,comparison pair count)",
        "frequency_edges": distance_FREQ.tolist(),
        "flank_GC_bin": "0-1, 2, or 3-4 G/C bases across two immediate CpG flanks at two sites",
        "overall_eligibility": "at least 100 overlap-pair weight and at least 6 unique contributing promoters in each kind",
        "distance_bin_eligibility": "at least 20 overlap-pair weight and at least 4 unique contributing promoters in each kind",
        "sensitivities": [
            "MAPQ >=60",
            "CpG separation >=25bp",
            "original panel only: both CpG dyads outside archived RepeatMasker annotations",
        ],
        "scope": "original panel, each patient-held-out panel, and 120 diagnosis-blind promoters",
        "limits": "Internal selected hypotheses; unknown repeat annotation outside original panel; phi is frequency-normalised but remains bounded by marginals; sparse strata are omitted with coverage reported; overlapping training folds are not independent replicates",
    }
    p = distance_F / "plans/F2_protocol.json"
    if p.exists():
        assert distance_json.loads(p.read_text()) == plan
    else:
        _regional.regional_save(p, plan)
    panel = distance_json.loads((distance_F / "plans/panel.json").read_text())
    positions = distance_np.load(distance_D / "reference/cpg_positions.npy", mmap_mode="r")
    offset = distance_np.load(distance_D / "reference/cpg_offsets.npy")
    common = distance_np.load(distance_D / "pooled_analysis/common_cpg_mask_depth5.npy")
    pairs = {}
    with distance_pysam.FastaFile(_config.reference) as fa:
        for p in panel:
            ch = int(p["chromosome"][3:]) - 1
            (lo, hi) = map(int, offset[ch : ch + 2])
            (a, b) = distance_np.searchsorted(positions[lo:hi], [p["start"], p["end"]])
            ix = distance_np.arange(lo + a, lo + b)
            cp = positions[ix][common[ix]]
            seq = fa.fetch(p["chromosome"], p["start"] - 2, p["end"] + 3).upper()
            flank = distance_np.array(
                [
                    int(seq[int(c) - p["start"] + 1] in "GC")
                    + int(seq[int(c) - p["start"] + 4] in "GC")
                    for c in cp
                ]
            )
            repeats = None
            if p["original_locus_id"]:
                path = (
                    distance_D
                    / "reference/story_context_snapshot"
                    / ("rmsk_" + p["original_locus_id"] + ".json.gz")
                )
                with distance_gzip.open(path, "rt") as h:
                    rr = distance_json.load(h)["rmsk"]
                repeats = distance_np.array(
                    [
                        any((r["genoStart"] < int(c) + 2 and r["genoEnd"] > int(c) for r in rr))
                        for c in cp
                    ]
                )
            (ii, jj) = distance_np.triu_indices(len(cp), 1)
            dist = cp[jj].astype(distance_np.int64) - cp[ii].astype(distance_np.int64)
            rows = []
            seed = int.from_bytes(
                distance_hashlib.sha256((p["locus_id"] + "distance").encode()).digest()[:8],
                "little",
            )
            rng = distance_np.random.default_rng(seed)
            for k, (low, high) in enumerate(zip(distance_BINS[:-1], distance_BINS[1:])):
                selected = distance_np.flatnonzero((dist >= low) & (dist < high))
                if len(selected) > 40:
                    selected = distance_np.sort(rng.choice(selected, 40, replace=False))
                for u in selected:
                    gc = int(flank[ii[u]] + flank[jj[u]])
                    rows.append(
                        {
                            "i": int(ii[u]),
                            "j": int(jj[u]),
                            "distance_bin": k,
                            "distance_bp": int(dist[u]),
                            "flank_gc_bin": 0 if gc <= 1 else 1 if gc == 2 else 2,
                            "original_nonrepeat": (
                                None
                                if repeats is None
                                else bool(not repeats[ii[u]] and (not repeats[jj[u]]))
                            ),
                        }
                    )
            pairs[p["locus_id"]] = rows
    _regional.regional_save(distance_F / "plans/distance_pairs.json", pairs)
    return pairs


def distance_metric(z, rng):
    n = len(z)
    h = z == 2
    code = h[:, 0].astype(int) + 2 * h[:, 1]
    counts = distance_np.bincount(code, minlength=4)
    p = h.mean(0)
    joint = distance_np.mean(h[:, 0] & h[:, 1])
    cov = (joint - p.prod()) * n / (n - 1) * 11 / 12 * 100
    phi = None
    valid = 0
    if (
        min(
            counts[1] + counts[3],
            counts[2] + counts[3],
            counts[0] + counts[1],
            counts[0] + counts[2],
        )
        >= 3
    ):
        draw = rng.multivariate_hypergeometric(counts, 12, size=64)
        x = (draw[:, 1] + draw[:, 3]) / 12
        y = (draw[:, 2] + draw[:, 3]) / 12
        den = distance_np.sqrt(x * (1 - x) * y * (1 - y))
        ok = den > 0
        valid = int(ok.sum())
        if valid >= 32:
            phi = float(distance_np.mean((draw[ok, 3] / 12 - x[ok] * y[ok]) / den[ok]))
    coding = distance_np.where(z == 2, 1, distance_np.where(z == 1, -1, 0))
    cm = float(distance_np.cov(coding.T, ddof=1)[0, 1] * 11 / 12)
    return {
        "complete_molecules": n,
        "h1_fraction": float(p[0]),
        "h2_fraction": float(p[1]),
        "joint_excess_pp": float(cov),
        "phi": phi,
        "phi_valid_resamples": valid,
        "h_minus_m_covariance": cm,
    }


def distance_sample(s):
    label = s["study_label"]
    out = distance_F / "results" / label
    done = out / "distance_complete.json"
    if done.exists():
        return label
    panel = distance_json.loads((distance_F / "plans/panel.json").read_text())
    pairs = distance_json.loads((distance_F / "plans/distance_pairs.json").read_text())
    rows = []
    for p in panel:
        z = distance_np.load(out / (p["locus_id"] + ".npz"))
        x = z["states"][:, z["common_mask"]]
        mapq = z["mapq"]
        for k, pair in enumerate(pairs[p["locus_id"]]):
            a = x[:, [pair["i"], pair["j"]]]
            complete = (a >= 0).all(1)
            base = None
            for scope, keep in [("standard", complete), ("MAPQ60", complete & (mapq >= 60))]:
                if keep.sum() < 12:
                    continue
                if scope == "MAPQ60" and distance_np.array_equal(keep, complete):
                    value = base.copy()
                else:
                    seed = int.from_bytes(
                        distance_hashlib.sha256(
                            (label + p["locus_id"] + str(k) + "F2").encode()
                        ).digest()[:8],
                        "little",
                    )
                    value = distance_metric(a[keep], distance_np.random.default_rng(seed))
                if scope == "standard":
                    base = value.copy()
                rows.append(
                    {
                        "study_label": label,
                        "diagnosis": s["diagnosis"],
                        "locus_id": p["locus_id"],
                        "pair_index": k,
                        "scope": scope,
                        **pair,
                        **value,
                    }
                )
    _regional.regional_table(out / "distance_pairs.tsv.gz", rows)
    _regional.regional_save(
        done,
        {
            "pair_scope_records": len(rows),
            "standard_pairs": sum((r["scope"] == "standard" for r in rows)),
            "phi_eligible_standard_pairs": sum(
                (r["scope"] == "standard" and r["phi"] is not None for r in rows)
            ),
        },
    )
    print(label, "DISTANCE COMPLETE", len(rows), flush=True)
    return label


def initialize_distance():
    """Initialize the distance stage once; load its declared inputs."""
    global distance_BINS, distance_D, distance_F, distance_FREQ, distance_R
    if _runtime.initialized("followup_distance"):
        return
    _runtime.begin("followup_distance")
    _runtime.initialize("pooled_analysis")
    distance_R = distance_Path(_config.workspace)
    distance_D = distance_R / ".analysis"
    distance_F = distance_D / "focused_followup"
    distance_BINS = distance_np.array([1, 10, 25, 50, 100, 250, 500, 1000, 2500])
    distance_FREQ = distance_np.array([0, 0.05, 0.15, 0.3, 0.5, 1.000001])
    _runtime.finish("followup_distance")


def run_distance():
    """Execute the distance workflow stage."""
    initialize_distance()
    global distance_a, distance_pool, distance_samples, distance_z
    distance_prepare()
    distance_z = distance_np.array([[2, 2]] * 12 + [[0, 0]] * 12)
    distance_a = distance_metric(distance_z, distance_np.random.default_rng(1))
    assert (
        abs(distance_a["joint_excess_pp"] - 100 * 11 / 12 * 24 / 23 * 0.25) < 1e-12
        and abs(distance_a["phi"] - 1) < 1e-12
    )
    distance_z = distance_np.array([[2, 2], [2, 0], [0, 2], [0, 0]] * 6)
    distance_a = distance_metric(distance_z, distance_np.random.default_rng(1))
    assert abs(distance_a["joint_excess_pp"]) < 1e-12
    _regional.regional_save(
        distance_F / "validation/F2_metric_checks.json",
        {
            "perfect_positive_binary_association_phi_one": True,
            "independent_balanced_pair_exact_covariance_zero": True,
            "known_fixed_depth_covariance_formula": True,
        },
    )
    _regional.regional_save(
        distance_F / "distance_job.json",
        {
            "pid": distance_os.getpid(),
            "workers": _config.workers,
            "started_utc": distance_time.strftime("%Y-%m-%dT%H:%M:%SZ", distance_time.gmtime()),
        },
    )
    distance_samples = sorted(
        distance_json.loads((distance_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    with _runtime.process_pool(max_workers=_config.workers) as distance_pool:
        list(distance_pool.map(distance_sample, distance_samples))
    _regional.regional_save(
        distance_F / "results/distance_analysis_complete.json", {"specimens": 20}
    )
    print("DISTANCE ANALYSIS COMPLETE", flush=True)


# MATCHING


from pathlib import Path as matching_Path
import json as matching_json, csv as matching_csv, gzip as matching_gzip, numpy as matching_np
import matplotlib.pyplot as matching_plt


def matching_get_arrays(label):
    with matching_gzip.open(matching_F / "results" / label / "distance_pairs.tsv.gz", "rt") as f:
        rows = list(matching_csv.DictReader(f, delimiter="\t"))
    out = {
        k: matching_np.array([r[k] for r in rows])
        for k in ["locus_id", "scope", "original_nonrepeat"]
    }
    for k in [
        "distance_bin",
        "distance_bp",
        "flank_gc_bin",
        "complete_molecules",
        "h1_fraction",
        "h2_fraction",
        "joint_excess_pp",
        "phi",
    ]:
        out[k] = matching_np.array(
            [float(r[k]) if r[k] not in ["", "None"] else matching_np.nan for r in rows]
        )
    out["low_h"] = matching_np.minimum(out["h1_fraction"], out["h2_fraction"])
    out["high_h"] = matching_np.maximum(out["h1_fraction"], out["h2_fraction"])
    lo = matching_np.searchsorted(distance_FREQ, out["low_h"], side="right") - 1
    hi = matching_np.searchsorted(distance_FREQ, out["high_h"], side="right") - 1
    out["stratum"] = (out["distance_bin"] * 75 + lo * 15 + hi * 3 + out["flank_gc_bin"]).astype(int)
    return out


def matching_match(a, candidate, comparison, keep, metric, minimum_weight, minimum_loci):
    good = keep & matching_np.isfinite(a[metric])
    cm = good & matching_np.isin(a["locus_id"], candidate)
    rm = good & matching_np.isin(a["locus_id"], comparison)
    strata = a["stratum"]
    n1 = matching_np.bincount(strata[cm], minlength=600)
    n0 = matching_np.bincount(strata[rm], minlength=600)
    weight = matching_np.minimum(n1, n0)
    total = int(weight.sum())
    shared = matching_np.flatnonzero(weight > 0)
    usec = cm & matching_np.isin(strata, shared)
    user = rm & matching_np.isin(strata, shared)
    lc = len(set(a["locus_id"][usec]))
    lr = len(set(a["locus_id"][user]))
    eligible = total >= minimum_weight and min(lc, lr) >= minimum_loci
    row = {
        "eligible": eligible,
        "overlap_pair_weight": total,
        "candidate_promoters_contributing": lc,
        "comparison_promoters_contributing": lr,
        "candidate_pairs_available": int(cm.sum()),
        "comparison_pairs_available": int(rm.sum()),
        "candidate_matched_mean": "",
        "comparison_matched_mean": "",
        "matched_difference": "",
        "unmatched_difference": (
            float(a[metric][cm].mean() - a[metric][rm].mean()) if cm.any() and rm.any() else ""
        ),
        "matched_low_h_difference_pp": "",
        "matched_high_h_difference_pp": "",
        "matched_log_distance_difference": "",
    }
    if total:

        def means(key):
            s1 = matching_np.bincount(strata[cm], weights=a[key][cm], minlength=600)
            s0 = matching_np.bincount(strata[rm], weights=a[key][rm], minlength=600)
            return (
                float(matching_np.sum(weight[shared] * s1[shared] / n1[shared]) / total),
                float(matching_np.sum(weight[shared] * s0[shared] / n0[shared]) / total),
            )

        (c, r) = means(metric)
        row.update(candidate_matched_mean=c, comparison_matched_mean=r, matched_difference=c - r)
        for key in ["low_h", "high_h"]:
            (c, r) = means(key)
            row["matched_" + key + "_difference_pp"] = (c - r) * 100
        (c, r) = means("log_distance")
        row["matched_log_distance_difference"] = c - r
    return row


def matching_sample(s):
    label = s["study_label"]
    a = matching_get_arrays(label)
    a["log_distance"] = matching_np.log(a["distance_bp"])
    fold = next(
        (
            f
            for f in matching_json.loads((matching_F / "plans/folds.json").read_text())
            if f["held_out"] == label
        )
    )
    panel = matching_json.loads((matching_F / "plans/panel.json").read_text())
    originalC = [p["locus_id"] for p in panel if p["original_locus_id"].startswith("C")]
    originalR = [p["locus_id"] for p in panel if p["original_locus_id"].startswith("R")]
    panels = {
        "heldout": (
            [p["candidate_locus_id"] for p in fold["pairs"]],
            [p["comparison_locus_id"] for p in fold["pairs"]],
        ),
        "original": (originalC, originalR),
    }
    out = []
    curves = []
    for name, (c, r) in panels.items():
        variants = ["standard", "MAPQ60", "distance_at_least_25bp"] + (
            ["original_nonrepeat"] if name == "original" else []
        )
        for variant in variants:
            keep = a["scope"] == ("MAPQ60" if variant == "MAPQ60" else "standard")
            if variant == "distance_at_least_25bp":
                keep &= a["distance_bp"] >= 25
            if variant == "original_nonrepeat":
                keep &= a["original_nonrepeat"] == "True"
            for metric in ["joint_excess_pp", "phi"]:
                row = matching_match(a, c, r, keep, metric, 100, 6)
                out.append(
                    {
                        "study_label": label,
                        "diagnosis": s["diagnosis"],
                        "panel": name,
                        "sensitivity": variant,
                        "metric": metric,
                        **row,
                    }
                )
                if variant != "standard":
                    continue
                for k in range(8):
                    row = matching_match(a, c, r, keep & (a["distance_bin"] == k), metric, 20, 4)
                    curves.append(
                        {
                            "study_label": label,
                            "diagnosis": s["diagnosis"],
                            "panel": name,
                            "metric": metric,
                            "distance_bin": k,
                            "distance_min_bp": int(distance_BINS[k]),
                            "distance_max_exclusive_bp": int(distance_BINS[k + 1]),
                            **row,
                        }
                    )
    blind = matching_json.loads((matching_F / "plans/blind_panel.json").read_text())
    background = []
    for k in range(8):
        for metric in ["joint_excess_pp", "phi"]:
            lm = []
            pairs = 0
            for locus in blind:
                good = (
                    (a["locus_id"] == locus)
                    & (a["scope"] == "standard")
                    & (a["distance_bin"] == k)
                    & matching_np.isfinite(a[metric])
                )
                n = int(good.sum())
                if n >= 5:
                    lm.append(float(a[metric][good].mean()))
                    pairs += n
            background.append(
                {
                    "study_label": label,
                    "diagnosis": s["diagnosis"],
                    "metric": metric,
                    "distance_bin": k,
                    "eligible_promoters": len(lm),
                    "eligible": len(lm) >= 30,
                    "pairs": pairs,
                    "mean_association": float(matching_np.mean(lm)) if len(lm) >= 30 else "",
                }
            )
    print(label, "DISTANCE SUMMARY COMPLETE", flush=True)
    return (out, curves, background)


def matching_run():
    assert (matching_F / "results/distance_analysis_complete.json").exists()
    samples = sorted(
        matching_json.loads((matching_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    with _runtime.process_pool(max_workers=_config.workers) as pool:
        result = list(pool.map(matching_sample, samples))
    overall = [r for (a, b, c) in result for r in a]
    curves = [r for (a, b, c) in result for r in b]
    background = [r for (a, b, c) in result for r in c]
    _regional.regional_table(matching_F / "results/F2_matched_specimen_scores.tsv", overall)
    _regional.regional_table(matching_F / "results/F2_distance_curves.tsv", curves)
    _regional.regional_table(matching_F / "results/F2_blind_distance_curves.tsv", background)
    labels = [s["study_label"] for s in samples]
    summary = []
    for panel in ["heldout", "original"]:
        for sensitivity in ["standard", "MAPQ60", "distance_at_least_25bp"] + (
            ["original_nonrepeat"] if panel == "original" else []
        ):
            for metric in ["joint_excess_pp", "phi"]:
                for diag in matching_COLORS:
                    rr = [
                        r
                        for r in overall
                        if r["panel"] == panel
                        and r["sensitivity"] == sensitivity
                        and (r["metric"] == metric)
                        and (r["diagnosis"] == diag)
                        and r["eligible"]
                    ]
                    v = [r["matched_difference"] for r in rr]
                    summary.append(
                        {
                            "panel": panel,
                            "sensitivity": sensitivity,
                            "metric": metric,
                            "diagnosis": diag,
                            "eligible_patients": len(v),
                            "positive_differences": sum((x > 0 for x in v)),
                            "mean_matched_difference": float(matching_np.mean(v)) if v else None,
                            "median_matched_difference": (
                                float(matching_np.median(v)) if v else None
                            ),
                            "max_abs_residual_low_h_difference_pp": (
                                max((abs(r["matched_low_h_difference_pp"]) for r in rr))
                                if rr
                                else None
                            ),
                            "max_abs_residual_high_h_difference_pp": (
                                max((abs(r["matched_high_h_difference_pp"]) for r in rr))
                                if rr
                                else None
                            ),
                        }
                    )
    _regional.regional_save(
        matching_F / "results/F2_summary.json",
        {
            "overall": summary,
            "distance_edges_bp": distance_BINS.tolist(),
            "no_independent_fold_p_values": True,
        },
    )
    (fig, axes) = matching_plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, metric in zip(axes, ["joint_excess_pp", "phi"]):
        for i, label in enumerate(labels):
            rr = next(
                (
                    r
                    for r in overall
                    if r["panel"] == "heldout"
                    and r["sensitivity"] == "standard"
                    and (r["metric"] == metric)
                    and (r["study_label"] == label)
                )
            )
            if rr["eligible"]:
                ax.plot(
                    [i - 0.15, i + 0.15],
                    [rr["unmatched_difference"], rr["matched_difference"]],
                    c=matching_COLORS[rr["diagnosis"]],
                    alpha=0.65,
                )
                ax.scatter(
                    [i - 0.15, i + 0.15],
                    [rr["unmatched_difference"], rr["matched_difference"]],
                    c=matching_COLORS[rr["diagnosis"]],
                    s=[12, 30],
                )
        ax.axhline(0, c="gray", lw=0.8)
        ax.set_xticks(range(20), labels, rotation=90, fontsize=7)
        ax.set_ylabel(
            "Candidate minus comparison "
            + ("excess co-occurrence (pp)" if metric == "joint_excess_pp" else "phi")
        )
        ax.set_title(
            "Fixed-depth excess co-occurrence"
            if metric == "joint_excess_pp"
            else "Frequency-normalised binary association"
        )
    fig.suptitle(
        "Held-out promoters: before and after within-patient matching\nSmall left: unmatched; large right: matched. Missing points have insufficient eligible overlap.",
        fontsize=11,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "F2_01_abundance_matched_association")
    mids = matching_np.sqrt(distance_BINS[:-1] * distance_BINS[1:])
    (fig, axes) = matching_plt.subplots(2, 2, figsize=(11, 8))
    for row, metric in enumerate(["joint_excess_pp", "phi"]):
        for col, diag in enumerate(matching_COLORS):
            ax = axes[row, col]
            group = []
            for label in [s["study_label"] for s in samples if s["diagnosis"] == diag]:
                vv = []
                for k in range(8):
                    r = next(
                        (
                            r
                            for r in curves
                            if r["panel"] == "heldout"
                            and r["metric"] == metric
                            and (r["study_label"] == label)
                            and (r["distance_bin"] == k)
                        )
                    )
                    vv.append(r["matched_difference"] if r["eligible"] else matching_np.nan)
                group.append(vv)
                ax.plot(mids, vv, color=matching_COLORS[diag], alpha=0.22, lw=0.8)
            v = matching_np.array(group)
            n = matching_np.isfinite(v).sum(0)
            mean = matching_np.divide(
                matching_np.nansum(v, axis=0),
                n,
                out=matching_np.full(8, matching_np.nan),
                where=n > 0,
            )
            ax.plot(mids, mean, "-o", c=matching_COLORS[diag], lw=2, ms=4)
            for x, y, nn in zip(mids, mean, n):
                if matching_np.isfinite(y):
                    ax.annotate(
                        "n=" + str(nn),
                        (x, y),
                        fontsize=7,
                        xytext=(0, 6),
                        textcoords="offset points",
                        ha="center",
                    )
            ax.axhline(0, c="gray", lw=0.7)
            ax.set_xscale("log")
            ax.set_xlabel("CpG separation (bp; bin midpoint)")
            ax.set_ylabel(
                "Matched candidate–comparison "
                + ("excess (pp)" if metric == "joint_excess_pp" else "phi")
            )
            ax.set_title(diag)
    fig.suptitle(
        "Physical scale of association at held-out promoters\nThin lines: patients; thick lines: means of eligible patients; missing bins remain missing",
        fontsize=11,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "F2_02_distance_dependence")
    (fig, axes) = matching_plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, metric in zip(axes, ["joint_excess_pp", "phi"]):
        for diag in matching_COLORS:
            vv = []
            for k in range(8):
                v = [
                    r["mean_association"]
                    for r in background
                    if r["diagnosis"] == diag
                    and r["metric"] == metric
                    and (r["distance_bin"] == k)
                    and r["eligible"]
                ]
                vv.append(float(matching_np.mean(v)) if v else matching_np.nan)
            ax.plot(mids, vv, "-o", c=matching_COLORS[diag], label=diag)
        ax.set_xscale("log")
        ax.set_xlim(2, 2500)
        ax.axhline(0, c="gray", lw=0.7)
        ax.set_xlabel("CpG separation (bp; bin midpoint)")
        ax.set_ylabel(
            "Mean " + ("excess co-occurrence (pp)" if metric == "joint_excess_pp" else "phi")
        )
        ax.set_title("Diagnosis-blind promoter reference")
        if any((r["eligible"] for r in background if r["metric"] == metric)):
            ax.legend(fontsize=8)
        else:
            ax.text(
                0.5,
                0.5,
                "No patient meets the requirement of\nat least 30 eligible promoters per distance bin",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=10,
                bbox=dict(facecolor="white", edgecolor="none"),
            )
    fig.suptitle(
        "Breadth of molecular association across the 120 frozen reference promoters\nEach eligible promoter contributes equally within a patient; no diagnosis-based locus selection",
        fontsize=11,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "F2_03_blind_panel_distance")
    (fig, axes) = matching_plt.subplots(1, 2, figsize=(11, 4.5))
    variants = ["standard", "MAPQ60", "distance_at_least_25bp"]
    for ax, diag in zip(axes, matching_COLORS):
        for label in [s["study_label"] for s in samples if s["diagnosis"] == diag]:
            rr = [
                next(
                    (
                        r
                        for r in overall
                        if r["panel"] == "heldout"
                        and r["sensitivity"] == v
                        and (r["metric"] == "phi")
                        and (r["study_label"] == label)
                    )
                )
                for v in variants
            ]
            y = [r["matched_difference"] if r["eligible"] else matching_np.nan for r in rr]
            ax.plot(range(3), y, "-o", c=matching_COLORS[diag], alpha=0.65, ms=4)
        ax.axhline(0, c="gray", lw=0.8)
        ax.set_xticks(range(3), ["Standard", "MAPQ ≥60", "Separation ≥25 bp"])
        ax.set_ylabel("Matched candidate minus comparison phi")
        ax.set_title(diag)
    fig.suptitle(
        "Sensitivity of frequency-normalised association at held-out promoters\nPatients with insufficient overlapping strata are not assigned a score",
        fontsize=11,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "F2_04_mapping_and_distance_sensitivity")
    _regional.regional_save(
        matching_F / "results/F2_complete.json",
        {"analysis_and_figures_generated": True, "requires_review": True},
    )
    print(
        "F2 SUMMARY COMPLETE",
        matching_json.dumps(
            [r for r in summary if r["panel"] == "heldout" and r["sensitivity"] == "standard"]
        ),
        flush=True,
    )


def initialize_matching():
    """Initialize the matching stage once; load its declared inputs."""
    global matching_COLORS, matching_D, matching_F, matching_R
    if _runtime.initialized("followup_distance_summary"):
        return
    _runtime.begin("followup_distance_summary")
    _runtime.initialize("followup_distance")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    matching_R = matching_Path(_config.workspace)
    matching_D = matching_R / ".analysis"
    matching_F = matching_D / "focused_followup"
    _molecules.blocks_S = matching_F
    matching_COLORS = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    _runtime.finish("followup_distance_summary")


def run_matching():
    """Execute the matching workflow stage."""
    initialize_matching()
    matching_run()


# ROW NULL


from pathlib import Path as row_null_Path
import json as row_null_json, csv as row_null_csv, hashlib as row_null_hashlib, os as row_null_os, numpy as row_null_np


def row_null_trade(z, rng, attempts):
    """Randomise jointly observed discordant entries of two binary rows.

    Every accepted trade preserves row totals, column totals and missingness.
    The matrix is modified in place; the return value counts nontrivial moves."""
    changed = 0
    for _ in range(attempts):
        (a, b) = rng.integers(0, len(z), 2)
        if a == b:
            continue
        idx = row_null_np.flatnonzero((z[a] >= 0) & (z[b] >= 0) & (z[a] != z[b]))
        n = int((z[a, idx] == 1).sum())
        if n == 0 or n == len(idx):
            continue
        chosen = rng.choice(len(idx), n, replace=False)
        new = row_null_np.zeros(len(idx), row_null_np.int8)
        new[chosen] = 1
        if row_null_np.array_equal(new, z[a, idx]):
            continue
        z[a, idx] = new
        z[b, idx] = 1 - new
        changed += 1
    return changed


def row_null_pair_cov(z, ii, jj, complete, n):
    a = z[:, ii] == 1
    b = z[:, jj] == 1
    p = (a & complete).sum(0) / n
    q = (b & complete).sum(0) / n
    j = (a & b & complete).sum(0) / n
    return (j - p * q) * n / (n - 1) * 11 / 12 * 100


def row_null_one_locus(label, key, pairs):
    data = row_null_np.load(row_null_F / "results" / label / (key + ".npz"))
    x = data["states"][:, data["common_mask"]]
    x = x[(x >= 0).sum(1) >= 10]
    z = row_null_np.where(x < 0, -1, (x == 2).astype(row_null_np.int8)).astype(row_null_np.int8)
    ii = row_null_np.array([p["i"] for p in pairs])
    jj = row_null_np.array([p["j"] for p in pairs])
    complete = (z[:, ii] >= 0) & (z[:, jj] >= 0)
    n = complete.sum(0)
    good = n >= 12
    if not good.any():
        return (
            [],
            {
                "study_label": label,
                "locus_id": key,
                "eligible": False,
                "reason": "no pair with 12 sufficiently observed molecules",
            },
        )
    ii = ii[good]
    jj = jj[good]
    complete = complete[:, good]
    n = n[good]
    pselected = [p for (p, g) in zip(pairs, good) if g]
    observed = row_null_pair_cov(z, ii, jj, complete, n)
    original_rows = (z == 1).sum(1)
    original_cols = (z == 1).sum(0)
    original_missing = z < 0
    null = []
    changed = []
    chain_means = []
    for chain in range(2):
        seed = int.from_bytes(
            row_null_hashlib.sha256((label + key + "rownull" + str(chain)).encode()).digest()[:8],
            "little",
        )
        rng = row_null_np.random.default_rng(seed)
        zz = z.copy()
        moves = row_null_trade(zz, rng, 50 * len(zz))
        values = []
        for _ in range(32):
            moves += row_null_trade(zz, rng, 5 * len(zz))
            values.append(row_null_pair_cov(zz, ii, jj, complete, n))
        assert (
            row_null_np.array_equal((zz == 1).sum(1), original_rows)
            and row_null_np.array_equal((zz == 1).sum(0), original_cols)
            and row_null_np.array_equal(zz < 0, original_missing)
        )
        null.extend(values)
        changed.append(moves)
        chain_means.append(row_null_np.mean(values, axis=0))
    null = row_null_np.asarray(null)
    eligible = min(changed) >= 64
    rows = []
    for k in range(8):
        ix = row_null_np.array([p["distance_bin"] == k for p in pselected])
        eligiblebin = eligible and ix.sum() >= 5
        rows.append(
            {
                "study_label": label,
                "locus_id": key,
                "distance_bin": k,
                "pairs": int(ix.sum()),
                "eligible": eligiblebin,
                "observed_excess_pp": float(observed[ix].mean()) if ix.any() else "",
                "row_column_null_excess_pp": float(null[:, ix].mean()) if ix.any() else "",
                "residual_excess_pp": (
                    float(observed[ix].mean() - null[:, ix].mean()) if eligiblebin else ""
                ),
                "chain_mean_difference_pp": (
                    float(row_null_np.mean(chain_means[0][ix] - chain_means[1][ix]))
                    if ix.any()
                    else ""
                ),
            }
        )
    return (
        rows,
        {
            "study_label": label,
            "locus_id": key,
            "eligible": eligible,
            "reason": "eligible" if eligible else "insufficient nontrivial trades",
            "retained_molecules": len(z),
            "pair_count": len(observed),
            "nontrivial_trades_per_chain": changed,
            "row_col_missingness_invariants_passed": True,
        },
    )


def row_null_sample(s):
    label = s["study_label"]
    out = row_null_F / "results" / label
    if (out / "row_null_complete.json").exists():
        with (out / "row_null_locus_bins.tsv").open() as f:
            rr = list(row_null_csv.DictReader(f, delimiter="\t"))
        return (rr, row_null_json.loads((out / "row_null_complete.json").read_text())["checks"])
    panel = row_null_json.loads((row_null_F / "plans/panel.json").read_text())
    fold = next(
        (
            f
            for f in row_null_json.loads((row_null_F / "plans/folds.json").read_text())
            if f["held_out"] == label
        )
    )
    keys = set((p["locus_id"] for p in panel if p["original_locus_id"]))
    keys.update(
        (p[k] for p in fold["pairs"] for k in ["candidate_locus_id", "comparison_locus_id"])
    )
    pairs = row_null_json.loads((row_null_F / "plans/distance_pairs.json").read_text())
    rows = []
    checks = []
    for key in sorted(keys):
        (rr, check) = row_null_one_locus(label, key, pairs[key])
        rows.extend(rr)
        checks.append(check)
    for r in rows:
        r["diagnosis"] = s["diagnosis"]
    _regional.regional_table(out / "row_null_locus_bins.tsv", rows)
    _regional.regional_save(out / "row_null_complete.json", {"checks": checks})
    print(label, "ROW NULL COMPLETE", len(keys), flush=True)
    return (rows, checks)


def initialize_row_null():
    """Initialize the row_null stage once; load its declared inputs."""
    global row_null_D, row_null_F, row_null_R
    if _runtime.initialized("followup_row_null"):
        return
    _runtime.begin("followup_row_null")
    _runtime.initialize("pooled_analysis")
    row_null_R = row_null_Path(_config.workspace)
    row_null_D = row_null_R / ".analysis"
    row_null_F = row_null_D / "focused_followup"
    _runtime.finish("followup_row_null")


def run_row_null():
    """Execute the row_null workflow stage."""
    initialize_row_null()
    global row_null_moves, row_null_original, row_null_out, row_null_p, row_null_plan, row_null_pool, row_null_samples, row_null_z
    row_null_plan = {
        "purpose": "Spatial arrangement beyond molecule-wide abundance and per-CpG abundance; sensitivity, not a calibrated significance test",
        "defined": "Before inspecting F2 matched-results outputs",
        "panels": "original and held-out promoter unions for each patient",
        "read_eligibility": "at least 10 observed common CpGs per molecule; pairs need >=12 complete molecules",
        "randomisation": "Choose two rows; among jointly observed columns with discordant h states, redistribute assignments while preserving each row count. Missing cells and each column h count remain unchanged.",
        "chains": 2,
        "burn_in_row_trades_per_read": 50,
        "samples_per_chain": 32,
        "between_sample_row_trades_per_read": 5,
        "minimum_nontrivial_trades_each_chain": 64,
        "minimum_pairs_per_locus_bin": 5,
        "patient_eligibility": "at least six matched promoter pairs for overall score; four for distance-bin score",
        "invariants": ["per-read h count", "per-CpG h count", "missingness mask"],
        "limits": "Structural missingness can restrict randomisation connectivity; mixing is not guaranteed. Chain differences and trade counts are reported. No formal permutation p values. Restricted-read subset differs from F2 primary pair analysis.",
    }
    row_null_p = row_null_F / "plans/F2_row_null_protocol.json"
    if row_null_p.exists():
        assert row_null_json.loads(row_null_p.read_text()) == row_null_plan
    else:
        _regional.regional_save(row_null_p, row_null_plan)
    row_null_z = row_null_np.array(
        [[1, 0, 0, 1], [0, 1, 1, 0], [1, -1, 0, 0], [0, -1, 1, 1]], row_null_np.int8
    )
    row_null_original = row_null_z.copy()
    row_null_moves = row_null_trade(row_null_z, row_null_np.random.default_rng(2), 1000)
    assert (
        row_null_moves > 0
        and row_null_np.array_equal((row_null_z == 1).sum(0), (row_null_original == 1).sum(0))
        and row_null_np.array_equal((row_null_z == 1).sum(1), (row_null_original == 1).sum(1))
        and row_null_np.array_equal(row_null_z < 0, row_null_original < 0)
    )
    _regional.regional_save(
        row_null_F / "validation/F2_row_trade_test.json",
        {"row_column_missingness_preserved": True, "nontrivial_trades": row_null_moves},
    )
    _regional.regional_save(
        row_null_F / "row_null_job.json", {"pid": row_null_os.getpid(), "workers": _config.workers}
    )
    row_null_samples = sorted(
        row_null_json.loads((row_null_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    with _runtime.process_pool(max_workers=_config.workers) as row_null_pool:
        row_null_out = list(row_null_pool.map(row_null_sample, row_null_samples))
    _regional.regional_table(
        row_null_F / "results/F2_row_null_locus_bins.tsv",
        [r for (rr, cc) in row_null_out for r in rr],
    )
    _regional.regional_save(
        row_null_F / "validation/F2_row_null_invariants.json",
        [r for (rr, cc) in row_null_out for r in cc],
    )
    _regional.regional_save(
        row_null_F / "results/row_null_analysis_complete.json", {"specimens": 20}
    )
    print("ROW NULL COMPLETE", flush=True)


# ROW SUMMARY


from pathlib import Path as row_summary_Path
import json as row_summary_json, csv as row_summary_csv, numpy as row_summary_np
import matplotlib.pyplot as row_summary_plt


def row_summary_run():
    assert (row_summary_F / "results/row_null_analysis_complete.json").exists()
    rows = list(
        row_summary_csv.DictReader(
            (row_summary_F / "results/F2_row_null_locus_bins.tsv").open(), delimiter="\t"
        )
    )
    lookup = {(r["study_label"], r["locus_id"], int(r["distance_bin"])): r for r in rows}
    folds = row_summary_json.loads((row_summary_F / "plans/folds.json").read_text())
    panel = row_summary_json.loads((row_summary_F / "plans/panel.json").read_text())
    old = {p["original_locus_id"]: p["locus_id"] for p in panel if p["original_locus_id"]}
    orig = [
        {"candidate_locus_id": old[f"C{i:02d}"], "comparison_locus_id": old[f"R{i:02d}"]}
        for i in range(1, 13)
    ]
    out = []
    curves = []
    for fold in folds:
        label = fold["held_out"]
        diag = next((r["diagnosis"] for r in rows if r["study_label"] == label))
        for name, pairs in [("heldout", fold["pairs"]), ("original", orig)]:
            pairmean = []
            bins = {k: [] for k in range(8)}
            for pair in pairs:
                values = []
                for k in range(8):
                    a = lookup.get((label, pair["candidate_locus_id"], k))
                    b = lookup.get((label, pair["comparison_locus_id"], k))
                    if (
                        a is None
                        or b is None
                        or a["eligible"] != "True"
                        or (b["eligible"] != "True")
                    ):
                        continue
                    v = [
                        float(a[m]) - float(b[m])
                        for m in [
                            "observed_excess_pp",
                            "row_column_null_excess_pp",
                            "residual_excess_pp",
                            "chain_mean_difference_pp",
                        ]
                    ]
                    values.append(v)
                    bins[k].append(v)
                if values:
                    pairmean.append(row_summary_np.mean(values, axis=0))
            n = len(pairmean)
            v = (
                row_summary_np.mean(pairmean, axis=0)
                if pairmean
                else row_summary_np.full(4, row_summary_np.nan)
            )
            row = {
                "study_label": label,
                "diagnosis": diag,
                "panel": name,
                "matched_promoter_pairs": n,
                "eligible": n >= 6,
                "observed_paired_difference_pp": float(v[0]),
                "null_expected_paired_difference_pp": float(v[1]),
                "residual_paired_difference_pp": float(v[2]),
                "chain1_residual_paired_difference_pp": float(v[2] - v[3] / 2),
                "chain2_residual_paired_difference_pp": float(v[2] + v[3] / 2),
            }
            assert abs(v[0] - v[1] - v[2]) < 1e-10 if n else True
            out.append(row)
            for k in range(8):
                n = len(bins[k])
                v = (
                    row_summary_np.mean(bins[k], axis=0)
                    if n
                    else row_summary_np.full(4, row_summary_np.nan)
                )
                curves.append(
                    {
                        "study_label": label,
                        "diagnosis": diag,
                        "panel": name,
                        "distance_bin": k,
                        "matched_promoter_pairs": n,
                        "eligible": n >= 4,
                        "observed_paired_difference_pp": float(v[0]),
                        "null_expected_paired_difference_pp": float(v[1]),
                        "residual_paired_difference_pp": float(v[2]),
                    }
                )
    _regional.regional_table(row_summary_F / "results/F2_row_null_specimen_scores.tsv", out)
    _regional.regional_table(row_summary_F / "results/F2_row_null_distance_curves.tsv", curves)
    summary = []
    for name in ["heldout", "original"]:
        for diag in row_summary_COLORS:
            rr = [r for r in out if r["panel"] == name and r["diagnosis"] == diag and r["eligible"]]
            summary.append(
                {
                    "panel": name,
                    "diagnosis": diag,
                    "eligible_patients": len(rr),
                    "mean_observed_paired_difference_pp": (
                        float(row_summary_np.mean([r["observed_paired_difference_pp"] for r in rr]))
                        if rr
                        else None
                    ),
                    "mean_null_expected_paired_difference_pp": (
                        float(
                            row_summary_np.mean(
                                [r["null_expected_paired_difference_pp"] for r in rr]
                            )
                        )
                        if rr
                        else None
                    ),
                    "mean_residual_paired_difference_pp": (
                        float(row_summary_np.mean([r["residual_paired_difference_pp"] for r in rr]))
                        if rr
                        else None
                    ),
                    "positive_residual_patients": sum(
                        (r["residual_paired_difference_pp"] > 0 for r in rr)
                    ),
                    "positive_in_both_chains": sum(
                        (
                            min(
                                r["chain1_residual_paired_difference_pp"],
                                r["chain2_residual_paired_difference_pp"],
                            )
                            > 0
                            for r in rr
                        )
                    ),
                    "max_abs_chain_difference_pp": (
                        max(
                            (
                                abs(
                                    r["chain1_residual_paired_difference_pp"]
                                    - r["chain2_residual_paired_difference_pp"]
                                )
                                for r in rr
                            )
                        )
                        if rr
                        else None
                    ),
                }
            )
    _regional.regional_save(
        row_summary_F / "results/F2_row_null_summary.json",
        {
            "groups": summary,
            "interpretation": "Constrained randomisation sensitivity; no formal significance test or guarantee of complete mixing",
        },
    )
    (fig, axes) = row_summary_plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, diag in zip(axes, row_summary_COLORS):
        rr = [
            r for r in out if r["panel"] == "heldout" and r["diagnosis"] == diag and r["eligible"]
        ]
        for r in rr:
            ax.plot(
                range(3),
                [
                    r["observed_paired_difference_pp"],
                    r["null_expected_paired_difference_pp"],
                    r["residual_paired_difference_pp"],
                ],
                "-o",
                c=row_summary_COLORS[diag],
                alpha=0.65,
                ms=4,
            )
        ax.axhline(0, c="gray", lw=0.8)
        ax.set_xticks(range(3), ["Observed", "Row/column null", "Residual"])
        ax.set_ylabel("Candidate minus comparison co-occurrence (pp)")
        ax.set_title(diag + f" (n={len(rr)})")
    fig.suptitle(
        "Held-out promoters: arrangement beyond molecule-wide and site-specific abundance\nRandomisation preserves each read’s and each CpG’s 5hmC counts and the missing-call mask",
        fontsize=10,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "F2_05_row_column_preserving_null")
    mids = row_summary_np.sqrt(distance_BINS[:-1] * distance_BINS[1:])
    (fig, axes) = row_summary_plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, diag in zip(axes, row_summary_COLORS):
        vectors = []
        for label in [
            f["held_out"]
            for f in folds
            if next((r["diagnosis"] for r in rows if r["study_label"] == f["held_out"])) == diag
        ]:
            y = []
            for k in range(8):
                r = next(
                    (
                        r
                        for r in curves
                        if r["study_label"] == label
                        and r["panel"] == "heldout"
                        and (r["distance_bin"] == k)
                    )
                )
                y.append(
                    r["residual_paired_difference_pp"] if r["eligible"] else row_summary_np.nan
                )
            vectors.append(y)
            ax.plot(mids, y, c=row_summary_COLORS[diag], alpha=0.25, lw=0.8)
        v = row_summary_np.array(vectors)
        n = row_summary_np.isfinite(v).sum(0)
        mean = row_summary_np.divide(
            row_summary_np.nansum(v, axis=0),
            n,
            out=row_summary_np.full(8, row_summary_np.nan),
            where=n > 0,
        )
        ax.plot(mids, mean, "-o", c=row_summary_COLORS[diag], lw=2, ms=4)
        for x, y, nn in zip(mids, mean, n):
            if row_summary_np.isfinite(y):
                ax.annotate(
                    "n=" + str(nn),
                    (x, y),
                    fontsize=7,
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center",
                )
        ax.axhline(0, c="gray", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlabel("CpG separation (bp; bin midpoint)")
        ax.set_ylabel("Residual candidate minus comparison co-occurrence (pp)")
        ax.set_title(diag)
    fig.suptitle(
        "Spatial residual after molecule- and CpG-abundance-preserving randomisation\nOnly bins with at least four matched promoter pairs per patient; changing eligibility is shown",
        fontsize=10,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "F2_06_spatial_residual_distance")
    print("ROW SUMMARY COMPLETE", row_summary_json.dumps(summary), flush=True)


def initialize_row_summary():
    """Initialize the row_summary stage once; load its declared inputs."""
    global row_summary_COLORS, row_summary_D, row_summary_F, row_summary_R
    if _runtime.initialized("followup_row_summary"):
        return
    _runtime.begin("followup_row_summary")
    _runtime.initialize("followup_distance")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    row_summary_R = row_summary_Path(_config.workspace)
    row_summary_D = row_summary_R / ".analysis"
    row_summary_F = row_summary_D / "focused_followup"
    _molecules.blocks_S = row_summary_F
    row_summary_COLORS = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    _runtime.finish("followup_row_summary")


def run_row_summary():
    """Execute the row_summary workflow stage."""
    initialize_row_summary()
    row_summary_run()


# RANGES


from pathlib import Path as ranges_Path
import csv as ranges_csv, json as ranges_json, numpy as ranges_np


def ranges_run():
    rows = list(
        ranges_csv.DictReader(
            (ranges_F / "results/F2_row_null_locus_bins.tsv").open(), delimiter="\t"
        )
    )
    lookup = {(r["study_label"], r["locus_id"], int(r["distance_bin"])): r for r in rows}
    out = []
    for fold in ranges_json.loads((ranges_F / "plans/folds.json").read_text()):
        label = fold["held_out"]
        diag = next((r["diagnosis"] for r in rows if r["study_label"] == label))
        for name, allowed in [
            ("25_to_249bp", [2, 3, 4]),
            ("at_least_25bp", list(range(2, 8))),
            ("at_least_250bp", [5, 6, 7]),
        ]:
            pairs = []
            for pair in fold["pairs"]:
                differences = []
                for k in allowed:
                    a = lookup.get((label, pair["candidate_locus_id"], k))
                    b = lookup.get((label, pair["comparison_locus_id"], k))
                    if (
                        a is not None
                        and b is not None
                        and (a["eligible"] == "True")
                        and (b["eligible"] == "True")
                    ):
                        differences.append(
                            float(a["residual_excess_pp"]) - float(b["residual_excess_pp"])
                        )
                if differences:
                    pairs.append(ranges_np.mean(differences))
            out.append(
                {
                    "study_label": label,
                    "diagnosis": diag,
                    "distance_scope": name,
                    "matched_promoter_pairs": len(pairs),
                    "eligible": len(pairs) >= 6,
                    "residual_difference_pp": (
                        float(ranges_np.mean(pairs)) if len(pairs) >= 6 else None
                    ),
                }
            )
    _regional.regional_table(ranges_F / "results/F2_row_null_distance_sensitivity.tsv", out)
    summary = []
    for name in ["25_to_249bp", "at_least_25bp", "at_least_250bp"]:
        for diag in ["Glioblastoma", "Meningioma"]:
            v = [
                r["residual_difference_pp"]
                for r in out
                if r["distance_scope"] == name and r["diagnosis"] == diag and r["eligible"]
            ]
            summary.append(
                {
                    "distance_scope": name,
                    "diagnosis": diag,
                    "eligible_patients": len(v),
                    "positive_patients": sum((x > 0 for x in v)),
                    "mean_residual_pp": float(ranges_np.mean(v)) if v else None,
                }
            )
    _regional.regional_save(
        ranges_F / "results/F2_row_null_distance_sensitivity_summary.json", summary
    )


def initialize_ranges():
    """Initialize the ranges stage once; load its declared inputs."""
    global ranges_F, ranges_R
    if _runtime.initialized("followup_distance_ranges"):
        return
    _runtime.begin("followup_distance_ranges")
    _runtime.initialize("pooled_analysis")
    ranges_R = ranges_Path(_config.workspace)
    ranges_F = ranges_R / ".analysis/focused_followup"
    _runtime.finish("followup_distance_ranges")


def run_ranges():
    """Execute the ranges workflow stage."""
    initialize_ranges()
    ranges_run()


# NATIVE


from pathlib import Path as native_Path
import sys as native_sys, json as native_json, hashlib as native_hashlib, ctypes as native_ctypes, numpy as native_np, os as native_os


def native_binary_score(x, ii, jj, encoding, seed):
    """Score complete CpG pairs under the two-chain constrained binary null.

    Encodings: 0=5hmC, 1=5mC, 2=combined. Six output columns are documented
    in native/strengthen_null.c. info reports moves per chain and invariants."""
    x = native_np.ascontiguousarray(x, native_np.int8)
    ii = native_np.ascontiguousarray(ii, native_np.int32)
    jj = native_np.ascontiguousarray(jj, native_np.int32)
    out = native_np.full((len(ii), 6), native_np.nan)
    info = native_np.zeros(3, dtype=native_np.int64)
    status = native_LIB.score(
        x.ctypes.data,
        len(x),
        x.shape[1],
        ii.ctypes.data,
        jj.ctypes.data,
        len(ii),
        encoding,
        seed,
        out.ctypes.data,
        info.ctypes.data,
    )
    assert status == 0 and info[2] == 1
    return (out, info)


def native_pairs_for(p, cp):
    if p.get("existing_locus_id"):
        return native_EXISTING_PAIRS[p["existing_locus_id"]]
    (ii, jj) = native_np.triu_indices(len(cp), 1)
    dist = cp[jj].astype(native_np.int64) - cp[ii].astype(native_np.int64)
    rng = native_np.random.default_rng(
        int.from_bytes(
            native_hashlib.sha256((p["locus_id"] + "distance").encode()).digest()[:8], "little"
        )
    )
    rows = []
    for k, (a, b) in enumerate(zip(native_EDGES[:-1], native_EDGES[1:])):
        selected = native_np.flatnonzero((dist >= a) & (dist < b))
        if len(selected) > 40:
            selected = native_np.sort(rng.choice(selected, 40, replace=False))
        rows.extend(
            (
                {"i": int(ii[j]), "j": int(jj[j]), "distance_bin": k, "distance_bp": int(dist[j])}
                for j in selected
            )
        )
    return rows


def native_score_locus(label, p, path, encodings):
    z = native_np.load(path)
    x = z["states"][:, z["common_mask"]]
    x = x[(x >= 0).sum(1) >= 10]
    cp = z["cpg_positions"][z["common_mask"]]
    pairs = native_pairs_for(p, cp)
    ii = native_np.array([a["i"] for a in pairs], native_np.int32)
    jj = native_np.array([a["j"] for a in pairs], native_np.int32)
    n = ((x[:, ii] >= 0) & (x[:, jj] >= 0)).sum(0)
    good = n >= 12
    ii = ii[good]
    jj = jj[good]
    bins = native_np.array([a["distance_bin"] for a in pairs])[good]
    rows = []
    if not len(ii):
        return ([], {"locus_id": p["locus_id"], "eligible": False, "reason": "no complete pair"})
    outs = {}
    infos = {}
    for encoding in encodings:
        seed = int.from_bytes(
            native_hashlib.sha256(
                (label + p["locus_id"] + native_ENC[encoding] + "native").encode()
            ).digest()[:8],
            "little",
        )
        (out, info) = native_binary_score(x, ii, jj, encoding, seed)
        outs[encoding] = out
        infos[encoding] = info
        for k in range(8):
            use = bins == k
            eligible = use.sum() >= 5 and min(info[:2]) >= 64
            rows.append(
                {
                    "study_label": label,
                    "universe_index": p["universe_index"],
                    "locus_id": p["locus_id"],
                    "encoding": native_ENC[encoding],
                    "support": "all_pairs",
                    "distance_bin": k,
                    "pairs": int(use.sum()),
                    "eligible": bool(eligible),
                    "observed_covariance_pp": float(out[use, 0].mean()) if use.any() else None,
                    "null_covariance_pp": float(out[use, 1].mean()) if use.any() else None,
                    "residual_pp": float((out[use, 0] - out[use, 1]).mean()) if eligible else None,
                    "chain_difference_pp": float(out[use, 2].mean()) if use.any() else None,
                    "residual_phi": None,
                }
            )
    if len(encodings) == 3:
        common = native_np.ones(len(ii), bool)
        for encoding in encodings:
            binary = x == 2 if encoding == 0 else x == 1 if encoding == 1 else x > 0
            valid = (x[:, ii] >= 0) & (x[:, jj] >= 0)
            a = (binary[:, ii] & valid).sum(0)
            b = (binary[:, jj] & valid).sum(0)
            nn = valid.sum(0)
            common &= (
                (native_np.minimum.reduce([a, b, nn - a, nn - b]) >= 3)
                & native_np.isfinite(outs[encoding][:, 3])
                & native_np.isfinite(outs[encoding][:, 4])
            )
        for encoding in encodings:
            out = outs[encoding]
            for k in range(8):
                use = (bins == k) & common
                eligible = use.sum() >= 5 and all((min(infos[e][:2]) >= 64 for e in encodings))
                rows.append(
                    {
                        "study_label": label,
                        "universe_index": p["universe_index"],
                        "locus_id": p["locus_id"],
                        "encoding": native_ENC[encoding],
                        "support": "joint_informative",
                        "distance_bin": k,
                        "pairs": int(use.sum()),
                        "eligible": bool(eligible),
                        "observed_covariance_pp": float(out[use, 0].mean()) if use.any() else None,
                        "null_covariance_pp": float(out[use, 1].mean()) if use.any() else None,
                        "residual_pp": (
                            float((out[use, 0] - out[use, 1]).mean()) if eligible else None
                        ),
                        "chain_difference_pp": float(out[use, 2].mean()) if use.any() else None,
                        "residual_phi": (
                            float((out[use, 3] - out[use, 4]).mean()) if eligible else None
                        ),
                    }
                )
    return (
        rows,
        {
            "locus_id": p["locus_id"],
            "eligible": True,
            "invariants_passed": True,
            "trades": {native_ENC[e]: infos[e][:2].tolist() for e in encodings},
            "encodings": list(map(lambda e: native_ENC[e], encodings)),
        },
    )


def native_sample(task):
    (label, panel, mode) = task
    outdir = native_S / "results" / label
    outdir.mkdir(exist_ok=True)
    rows = []
    checks = []
    for p in panel:
        out = outdir / (p["locus_id"] + "_" + mode + ".json")
        if out.exists():
            data = native_json.loads(out.read_text())
            rows.extend(data["rows"])
            checks.append(data["check"])
            continue
        path = (
            native_F / "results" / label / (p["existing_locus_id"] + ".npz")
            if p.get("existing_locus_id")
            else native_S / "matrices" / label / (p["locus_id"] + ".npz")
        )
        (rr, check) = native_score_locus(label, p, path, [0, 1, 2] if mode == "three" else [0])
        _regional.regional_save(out, {"rows": rr, "check": check})
        rows.extend(rr)
        checks.append(check)
    _regional.regional_table(outdir / (mode + "_locus_bins.tsv"), rows)
    _regional.regional_save(outdir / (mode + "_checks.json"), checks)
    print(label, mode, "SCORED", len(panel), flush=True)
    return (rows, checks)


def native_validate():
    rng = native_np.random.default_rng(7)
    x = rng.integers(0, 3, (40, 30), dtype=native_np.int8)
    x[rng.random(x.shape) < 0.15] = -1
    ii = native_np.array([0, 2, 4, 8])
    jj = native_np.array([1, 3, 9, 12])
    checks = []
    for e in range(3):
        (out, info) = native_binary_score(x, ii, jj, e, 42)
        for k, (i, j) in enumerate(zip(ii, jj)):
            z = x[:, [i, j]]
            z = z[(z >= 0).all(1)]
            a = z == 2 if e == 0 else z == 1 if e == 1 else z > 0
            expected = native_np.cov(a.T, ddof=1)[0, 1] * 11 / 12 * 100
            assert abs(out[k, 0] - expected) < 1e-10
        checks.append(
            {
                "encoding": native_ENC[e],
                "direct_covariance_matches_numpy": True,
                "invariants_passed": bool(info[2]),
                "moves": info[:2].tolist(),
            }
        )
    _regional.regional_save(native_S / "validation/native_null_validation.json", checks)


def initialize_native():
    """Initialize the native stage once; load its declared inputs."""
    global native_D, native_EDGES, native_ENC, native_EXISTING_PAIRS, native_F, native_LIB, native_R, native_S
    if _runtime.initialized("strengthen_score"):
        return
    _runtime.begin("strengthen_score")
    _runtime.initialize("pooled_analysis")
    native_R = native_Path(_config.workspace)
    native_D = native_R / ".analysis"
    native_F = native_D / "focused_followup"
    native_S = native_D / "strengthening"
    native_LIB = native_ctypes.CDLL(str(_runtime.native_path("strengthen_null.so")))
    native_LIB.score.argtypes = [
        native_ctypes.c_void_p,
        native_ctypes.c_int,
        native_ctypes.c_int,
        native_ctypes.c_void_p,
        native_ctypes.c_void_p,
        native_ctypes.c_int,
        native_ctypes.c_int,
        native_ctypes.c_uint64,
        native_ctypes.c_void_p,
        native_ctypes.c_void_p,
    ]
    native_LIB.score.restype = native_ctypes.c_int
    native_EDGES = native_np.array([1, 10, 25, 50, 100, 250, 500, 1000, 2500])
    native_ENC = ["5hmC", "5mC", "combined"]
    native_EXISTING_PAIRS = (
        native_json.loads((native_F / "plans/distance_pairs.json").read_text())
        if (native_F / "plans/distance_pairs.json").exists()
        else {}
    )
    _runtime.finish("strengthen_score")


def run_native():
    """Execute the native workflow stage."""
    initialize_native()
    global native_fold, native_folds, native_keys, native_mode, native_needed, native_old, native_panel, native_pool, native_pp, native_result, native_s, native_samples, native_tasks
    native_validate()
    native_mode = native_sys.argv[1] if len(native_sys.argv) > 1 else "three"
    native_samples = sorted(
        native_json.loads((native_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    native_tasks = []
    if native_mode == "three":
        native_old = native_json.loads((native_F / "plans/panel.json").read_text())
        native_folds = native_json.loads((native_F / "plans/folds.json").read_text())
        for native_s in native_samples:
            native_fold = next(
                (f for f in native_folds if f["held_out"] == native_s["study_label"])
            )
            native_keys = set(
                (
                    v[k]
                    for v in native_fold["pairs"]
                    for k in ["candidate_locus_id", "comparison_locus_id"]
                )
            )
            native_pp = [
                {**p, "locus_id": f"U{p['universe_index']:05d}", "existing_locus_id": p["locus_id"]}
                for p in native_old
                if p["locus_id"] in native_keys
            ]
            native_tasks.append((native_s["study_label"], native_pp, native_mode))
    else:
        native_needed = native_json.loads((native_S / "plans/needed_by_patient.json").read_text())
        native_panel = native_json.loads((native_S / "plans/panel.json").read_text())
        for native_s in native_samples:
            native_tasks.append(
                (
                    native_s["study_label"],
                    [
                        p
                        for p in native_panel
                        if p["universe_index"] in native_needed[native_s["study_label"]]
                    ],
                    native_mode,
                )
            )
    _regional.regional_save(
        native_S / (native_mode + "_job.json"),
        {"pid": native_os.getpid(), "workers": _config.workers},
    )
    with _runtime.process_pool(max_workers=_config.workers) as native_pool:
        native_result = list(native_pool.map(native_sample, native_tasks))
    _regional.regional_table(
        native_S / "results" / (native_mode + "_locus_bins.tsv"),
        [r for (rr, cc) in native_result for r in rr],
    )
    _regional.regional_save(
        native_S / "validation" / (native_mode + "_all_checks.json"),
        [c for (rr, cc) in native_result for c in cc],
    )
    _regional.regional_save(
        native_S / "results" / (native_mode + "_scoring_complete.json"), {"patients": 20}
    )
    print(native_mode, "SCORING COMPLETE", flush=True)
