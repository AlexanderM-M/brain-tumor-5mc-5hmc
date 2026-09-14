"""Robustness analysis stages for the brain tumour 5mC/5hmC study.

Each stage has explicit initialization and execution; importing is read-only.
See docs/workflow.md for inputs, stage order and reproduction limits.
"""

from . import runtime as _runtime

from .runtime import config as _config

from . import molecules as _molecules
from . import regional as _regional


# INFLUENCE


from pathlib import Path as influence_Path
import json as influence_json, csv as influence_csv, numpy as influence_np
import matplotlib.pyplot as influence_plt


def influence_run():
    rows = list(
        influence_csv.DictReader(
            (influence_F / "results/F2_row_null_locus_bins.tsv").open(), delimiter="\t"
        )
    )
    lookup = {(r["study_label"], r["locus_id"], int(r["distance_bin"])): r for r in rows}
    folds = influence_json.loads((influence_F / "plans/folds.json").read_text())
    panel = influence_json.loads((influence_F / "plans/panel.json").read_text())
    genes = {p["locus_id"]: p["gene"] for p in panel}
    labels = [f["held_out"] for f in folds]
    pairs = []
    patient = []
    deletions = []
    for scope, bins in [("all_distances", range(8)), ("25_249bp", [2, 3, 4])]:
        for f in folds:
            label = f["held_out"]
            diag = next((r["diagnosis"] for r in rows if r["study_label"] == label))
            pp = []
            for pair in f["pairs"]:
                values = []
                for k in bins:
                    a = lookup.get((label, pair["candidate_locus_id"], k))
                    b = lookup.get((label, pair["comparison_locus_id"], k))
                    if a and b and (a["eligible"] == "True") and (b["eligible"] == "True"):
                        values.append(
                            float(a["residual_excess_pp"]) - float(b["residual_excess_pp"])
                        )
                if values:
                    row = {
                        "study_label": label,
                        "diagnosis": diag,
                        "scope": scope,
                        "candidate": pair["candidate_locus_id"],
                        "candidate_gene": genes[pair["candidate_locus_id"]],
                        "comparison": pair["comparison_locus_id"],
                        "comparison_gene": genes[pair["comparison_locus_id"]],
                        "eligible_bins": len(values),
                        "pair_residual_pp": float(influence_np.mean(values)),
                    }
                    pp.append(row)
                    pairs.append(row)
            baseline = (
                influence_np.mean([p["pair_residual_pp"] for p in pp])
                if len(pp) >= 6
                else influence_np.nan
            )
            leave = []
            for j, p in enumerate(pp):
                value = (
                    float(
                        influence_np.mean(
                            [x["pair_residual_pp"] for (k, x) in enumerate(pp) if k != j]
                        )
                    )
                    if len(pp) - 1 >= 6
                    else None
                )
                deletions.append(
                    {
                        **p,
                        "remaining_pairs": len(pp) - 1,
                        "eligible_after_removal": value is not None,
                        "baseline_pp": float(baseline),
                        "after_pair_removal_pp": value,
                        "change_pp": None if value is None else value - baseline,
                    }
                )
                if value is not None:
                    leave.append(value)
            patient.append(
                {
                    "study_label": label,
                    "diagnosis": diag,
                    "scope": scope,
                    "eligible_pairs": len(pp),
                    "baseline_pp": float(baseline),
                    "eligible_single_pair_deletions": len(leave),
                    "minimum_after_removal_pp": float(min(leave)) if leave else None,
                    "maximum_after_removal_pp": float(max(leave)) if leave else None,
                    "all_eligible_deletions_positive": bool(leave and min(leave) > 0),
                    "at_least_one_sign_flip": bool(
                        leave and any((v * baseline <= 0 for v in leave))
                    ),
                }
            )
    _regional.regional_table(influence_S / "results/P1_pair_contributions.tsv", pairs)
    _regional.regional_table(influence_S / "results/P1_pair_deletions.tsv", deletions)
    _regional.regional_table(influence_S / "results/P1_patient_influence.tsv", patient)
    globalrows = []
    allloci = sorted(set((p[k] for p in pairs for k in ["candidate", "comparison"])))
    for scope in ["all_distances", "25_249bp"]:
        for locus in allloci:
            groupvalues = {"Glioblastoma": [], "Meningioma": []}
            affected = 0
            for label in labels:
                pp = [p for p in pairs if p["scope"] == scope and p["study_label"] == label]
                keep = [p for p in pp if locus not in [p["candidate"], p["comparison"]]]
                affected += len(keep) != len(pp)
                if len(keep) >= 6:
                    groupvalues[pp[0]["diagnosis"]].append(
                        float(influence_np.mean([p["pair_residual_pp"] for p in keep]))
                    )
            mg = {
                d: float(influence_np.mean(v)) if v else influence_np.nan
                for (d, v) in groupvalues.items()
            }
            globalrows.append(
                {
                    "scope": scope,
                    "removed_locus": locus,
                    "gene": genes[locus],
                    "patients_affected": affected,
                    "eligible_GBM": len(groupvalues["Glioblastoma"]),
                    "eligible_MEN": len(groupvalues["Meningioma"]),
                    "GBM_mean_after_removal_pp": mg["Glioblastoma"],
                    "MEN_mean_after_removal_pp": mg["Meningioma"],
                    "MEN_minus_GBM_after_removal_pp": mg["Meningioma"] - mg["Glioblastoma"],
                    "MEN_positive_after_removal": sum((x > 0 for x in groupvalues["Meningioma"])),
                }
            )
    _regional.regional_table(influence_S / "results/P1_global_locus_deletions.tsv", globalrows)
    colors = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    (fig, axes) = influence_plt.subplots(1, 2, figsize=(13, 4.5))
    for ax, scope in zip(axes, ["all_distances", "25_249bp"]):
        rr = [r for r in patient if r["scope"] == scope]
        for i, r in enumerate(rr):
            if r["minimum_after_removal_pp"] is not None:
                ax.plot(
                    [i, i],
                    [r["minimum_after_removal_pp"], r["maximum_after_removal_pp"]],
                    c=colors[r["diagnosis"]],
                    lw=2,
                )
            ax.scatter(i, r["baseline_pp"], c=colors[r["diagnosis"]], s=24, zorder=3)
        ax.axhline(0, c="gray", lw=0.8)
        ax.set_xticks(range(20), labels, rotation=90, fontsize=7)
        ax.set_ylabel("Candidate minus comparison residual (pp)")
        ax.set_title("All distances" if scope == "all_distances" else "Exploratory 25–249 bp range")
    fig.suptitle(
        "Promoter-pair influence: dots are full-panel scores; lines span single-pair deletions\nRanges describe sensitivity, not confidence intervals; panels are not reselected after deletion",
        fontsize=10,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P1_01_patient_influence")
    cand = sorted(set((p["candidate"] for p in pairs)))
    mat = influence_np.full((len(cand), 20), influence_np.nan)
    for p in pairs:
        if p["scope"] == "25_249bp":
            mat[cand.index(p["candidate"]), labels.index(p["study_label"])] = p["pair_residual_pp"]
    lim = influence_np.nanpercentile(abs(mat), 98)
    (fig, ax) = influence_plt.subplots(figsize=(12, 9))
    cmap = influence_plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#ededed")
    im = ax.imshow(mat, aspect="auto", vmin=-lim, vmax=lim, cmap=cmap)
    ax.set_xticks(range(20), labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(cand)), [genes[k] for k in cand], fontsize=7)
    ax.set_title(
        "Individual held-out candidate–comparison pair contributions at 25–249 bp\nComparison partners can differ between folds; grey means not selected or not eligible"
    )
    fig.colorbar(
        im, ax=ax, label="Pair residual difference (pp; colour range clipped at 98th percentile)"
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P1_02_promoter_contributions")
    rr = [r for r in globalrows if r["scope"] == "25_249bp"]
    baseline = {
        d: influence_np.mean(
            [r["baseline_pp"] for r in patient if r["scope"] == "25_249bp" and r["diagnosis"] == d]
        )
        for d in colors
    }
    base = baseline["Meningioma"] - baseline["Glioblastoma"]
    rr = sorted(rr, key=lambda r: -abs(r["MEN_minus_GBM_after_removal_pp"] - base))[:20]
    (fig, ax) = influence_plt.subplots(figsize=(9, 7))
    ax.barh(
        range(len(rr)), [r["MEN_minus_GBM_after_removal_pp"] - base for r in rr], color="#497ba6"
    )
    ax.set_yticks(range(len(rr)), [r["gene"] for r in rr], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, c="gray", lw=0.8)
    ax.set_xlabel("Change in MEN-minus-GBM mean residual after removing locus (pp)")
    ax.set_title(
        "Largest locus influences at 25–249 bp\nLocus removed wherever selected as candidate or comparison"
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P1_03_global_locus_influence")
    summary = {
        "groups": [],
        "global_deletions": len(globalrows),
        "unique_loci_deleted": len(allloci),
    }
    for scope in ["all_distances", "25_249bp"]:
        for diag in colors:
            rr = [r for r in patient if r["scope"] == scope and r["diagnosis"] == diag]
            summary["groups"].append(
                {
                    "scope": scope,
                    "diagnosis": diag,
                    "patients": len(rr),
                    "patients_positive_after_every_eligible_single_pair_deletion": sum(
                        (r["all_eligible_deletions_positive"] for r in rr)
                    ),
                    "patients_with_a_sign_flip": sum((r["at_least_one_sign_flip"] for r in rr)),
                    "mean_baseline_pp": float(influence_np.mean([r["baseline_pp"] for r in rr])),
                }
            )
    _regional.regional_save(influence_S / "results/P1_summary.json", summary)
    (influence_S / "results/P1_ANSWER.md").write_text(
        "# P1: promoter influence\n\n"
        + influence_json.dumps(summary, indent=2)
        + "\n\nSingle-pair deletions keep the selected panels fixed and require six remaining eligible pairs. Global deletions remove a locus wherever it appears as candidate or comparison. These are influence diagnostics, not repeated discovery, confidence intervals or independent validation. The 25–249bp grouping remains exploratory. Pair-level and global tables retain all results, including loss of eligibility.\n"
    )
    _regional.regional_save(
        influence_S / "results/P1_complete.json",
        {"analysis_and_figures_generated": True, "requires_review": True},
    )
    print("P1 COMPLETE", influence_json.dumps(summary), flush=True)


def initialize_influence():
    """Initialize the influence stage once; load its declared inputs."""
    global influence_D, influence_F, influence_R, influence_S
    if _runtime.initialized("strengthen_influence"):
        return
    _runtime.begin("strengthen_influence")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    influence_R = influence_Path(_config.workspace)
    influence_D = influence_R / ".analysis"
    influence_F = influence_D / "focused_followup"
    influence_S = influence_D / "strengthening"
    _molecules.blocks_S = influence_S
    _runtime.finish("strengthen_influence")


def run_influence():
    """Execute the influence workflow stage."""
    initialize_influence()
    influence_run()


# ENCODINGS


from pathlib import Path as encodings_Path
import json as encodings_json, csv as encodings_csv, numpy as encodings_np
import matplotlib.pyplot as encodings_plt


def encodings_run():
    rows = list(
        encodings_csv.DictReader(
            (encodings_S / "results/three_locus_bins.tsv").open(), delimiter="\t"
        )
    )
    lookup = {
        (
            r["study_label"],
            int(r["universe_index"]),
            r["encoding"],
            r["support"],
            int(r["distance_bin"]),
        ): r
        for r in rows
    }
    folds = encodings_json.loads((encodings_F / "plans/folds.json").read_text())
    samples = {
        s["study_label"]: s
        for s in encodings_json.loads((encodings_D / "cohort.json").read_text())["samples"]
    }
    out = []
    for fold in folds:
        label = fold["held_out"]
        diag = samples[label]["diagnosis"]
        for scope, bins in [("all_distances", range(8)), ("25_249bp", [2, 3, 4])]:
            for support in ["all_pairs", "joint_informative"]:
                values = []
                phi = []
                for pair in fold["pairs"]:
                    vi = []
                    pi = []
                    for k in bins:
                        rr = [
                            [
                                lookup.get((label, pair[kind + "_index"], enc, support, k))
                                for kind in ["candidate", "comparison"]
                            ]
                            for enc in encodings_ENC
                        ]
                        if any(
                            (
                                a is None
                                or b is None
                                or a["eligible"] != "True"
                                or (b["eligible"] != "True")
                                for (a, b) in rr
                            )
                        ):
                            continue
                        assert (
                            len({a["pairs"] for (a, b) in rr}) == 1
                            and len({b["pairs"] for (a, b) in rr}) == 1
                        )
                        vi.append(
                            [float(a["residual_pp"]) - float(b["residual_pp"]) for (a, b) in rr]
                        )
                        if support == "joint_informative":
                            pi.append(
                                [
                                    float(a["residual_phi"]) - float(b["residual_phi"])
                                    for (a, b) in rr
                                ]
                            )
                    if vi:
                        values.append(encodings_np.mean(vi, axis=0))
                    if pi:
                        phi.append(encodings_np.mean(pi, axis=0))
                for e, encoding in enumerate(encodings_ENC):
                    out.append(
                        {
                            "study_label": label,
                            "diagnosis": diag,
                            "scope": scope,
                            "support": support,
                            "encoding": encoding,
                            "matched_promoter_pairs": len(values),
                            "eligible": len(values) >= 6,
                            "residual_pp": (
                                float(encodings_np.mean(values, axis=0)[e])
                                if len(values) >= 6
                                else None
                            ),
                            "residual_phi": (
                                float(encodings_np.mean(phi, axis=0)[e]) if len(phi) >= 6 else None
                            ),
                        }
                    )
    _regional.regional_table(encodings_S / "results/P2_specimen_comparison.tsv", out)
    summary = []
    for scope in ["all_distances", "25_249bp"]:
        for support in ["all_pairs", "joint_informative"]:
            for diag in encodings_COL:
                for encoding in encodings_ENC:
                    rr = [
                        r
                        for r in out
                        if r["scope"] == scope
                        and r["support"] == support
                        and (r["diagnosis"] == diag)
                        and (r["encoding"] == encoding)
                        and r["eligible"]
                    ]
                    summary.append(
                        {
                            "scope": scope,
                            "support": support,
                            "diagnosis": diag,
                            "encoding": encoding,
                            "eligible_patients": len(rr),
                            "mean_residual_pp": (
                                float(encodings_np.mean([r["residual_pp"] for r in rr]))
                                if rr
                                else None
                            ),
                            "positive_patients": sum((r["residual_pp"] > 0 for r in rr)),
                            "mean_residual_phi": (
                                float(encodings_np.mean([r["residual_phi"] for r in rr]))
                                if rr and support == "joint_informative"
                                else None
                            ),
                        }
                    )
    _regional.regional_save(
        encodings_S / "results/P2_summary.json",
        {
            "comparisons": summary,
            "same_physical_pairs_in_primary_three_encoding_comparison": True,
            "phi_uses_full_complete_pair_counts": "Secondary normalised residual, distinct from prior F2 subsampled phi",
        },
    )
    (fig, axes) = encodings_plt.subplots(2, 2, figsize=(11, 8))
    for row, scope in enumerate(["all_distances", "25_249bp"]):
        for col, diag in enumerate(encodings_COL):
            ax = axes[row, col]
            n = 0
            for label in [
                f["held_out"] for f in folds if samples[f["held_out"]]["diagnosis"] == diag
            ]:
                rr = [
                    next(
                        (
                            r
                            for r in out
                            if r["study_label"] == label
                            and r["scope"] == scope
                            and (r["support"] == "all_pairs")
                            and (r["encoding"] == e)
                        )
                    )
                    for e in encodings_ENC
                ]
                if all((r["eligible"] for r in rr)):
                    ax.plot(
                        range(3),
                        [r["residual_pp"] for r in rr],
                        "-o",
                        c=encodings_COL[diag],
                        alpha=0.65,
                        ms=4,
                    )
                    n += 1
            ax.axhline(0, c="gray", lw=0.7)
            ax.set_xticks(range(3), encodings_ENC)
            ax.set_ylabel("Candidate minus comparison residual (pp)")
            ax.set_title(
                diag
                + f" (n={n})"
                + ("\nAll distances" if scope == "all_distances" else "\nExploratory 25–249 bp")
            )
    fig.suptitle(
        "Separate and combined modifications on shared reads and CpG pairs\nEach encoding has its own molecule/CpG-count-preserving null; residual scales still depend on prevalence",
        fontsize=10,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P2_01_encoding_comparison")
    labels = [f["held_out"] for f in folds]
    (fig, axes) = encodings_plt.subplots(1, 2, figsize=(12, 4))
    for ax, support in zip(axes, ["all_pairs", "joint_informative"]):
        m = encodings_np.array(
            [
                [
                    next(
                        (
                            r["matched_promoter_pairs"]
                            for r in out
                            if r["study_label"] == label
                            and r["scope"] == scope
                            and (r["support"] == support)
                            and (r["encoding"] == "5hmC")
                        )
                    )
                    for label in labels
                ]
                for scope in ["all_distances", "25_249bp"]
            ]
        )
        im = ax.imshow(m, aspect="auto", cmap="viridis", vmin=0, vmax=12)
        ax.set_xticks(range(20), labels, rotation=90, fontsize=7)
        ax.set_yticks([0, 1], ["All distances", "25–249 bp"])
        ax.set_title(
            "Shared raw-pair support"
            if support == "all_pairs"
            else "Shared informative support for phi"
        )
        fig.colorbar(im, ax=ax, label="Matched promoter pairs (at least 6 needed)")
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P2_02_shared_support")
    (fig, axes) = encodings_plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, diag in zip(axes, encodings_COL):
        n = 0
        for label in [f["held_out"] for f in folds if samples[f["held_out"]]["diagnosis"] == diag]:
            rr = [
                next(
                    (
                        r
                        for r in out
                        if r["study_label"] == label
                        and r["scope"] == "25_249bp"
                        and (r["support"] == "joint_informative")
                        and (r["encoding"] == e)
                    )
                )
                for e in encodings_ENC
            ]
            if all((r["eligible"] for r in rr)):
                ax.plot(
                    range(3),
                    [r["residual_phi"] for r in rr],
                    "-o",
                    c=encodings_COL[diag],
                    alpha=0.7,
                )
                n += 1
        ax.set_xticks(range(3), encodings_ENC)
        ax.axhline(0, c="gray", lw=0.7)
        ax.set_ylabel("Candidate minus comparison residual phi")
        ax.set_title(diag + f" (n={n})")
        if not n:
            ax.text(
                0.5,
                0.5,
                "Insufficient shared informative pairs",
                transform=ax.transAxes,
                ha="center",
                bbox=dict(facecolor="white", edgecolor="none"),
            )
    fig.suptitle(
        "Secondary frequency-normalised comparison at 25–249 bp\nSame informative pairs in all three encodings; full-pair phi differs from prior subsampled F2 phi",
        fontsize=10,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P2_03_normalised_comparison")
    panel = encodings_json.loads((encodings_F / "plans/panel.json").read_text())
    byid = {p["locus_id"]: p for p in panel}
    frequency = {
        p["locus_id"]: sum(
            (any((q["candidate_locus_id"] == p["locus_id"] for q in f["pairs"])) for f in folds)
        )
        for p in panel
    }
    key = min(frequency, key=lambda k: (-frequency[k], byid[k]["gene"], k))
    p = byid[key]
    from matplotlib.colors import ListedColormap, BoundaryNorm

    (fig, axes) = encodings_plt.subplots(2, 4, figsize=(13, 6))
    chosen = ["GBM-01", "GBM-02", "MEN-01", "MEN-02"]
    cmap = ListedColormap(["#eeeeee", "#cfcfcf", "#497ba6", "#df9d29"])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], 4)
    for col, label in enumerate(chosen):
        z = encodings_np.load(encodings_F / "results" / label / (key + ".npz"))
        x = z["states"][:, z["common_mask"]]
        x = x[(x >= 0).sum(1) >= 10]
        rng = encodings_np.random.default_rng(20260912)
        if len(x) > 50:
            x = x[encodings_np.sort(rng.choice(len(x), 50, replace=False))]
        axes[0, col].imshow(x, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
        collapsed = encodings_np.where(x < 0, -1, encodings_np.where(x > 0, 1, 0))
        axes[1, col].imshow(
            collapsed,
            aspect="auto",
            interpolation="nearest",
            cmap=ListedColormap(["#eeeeee", "#cfcfcf", "#147D92"]),
            vmin=-1,
            vmax=1,
        )
        for row in [0, 1]:
            axes[row, col].set_title(
                label + (" — separate" if row == 0 else " — combined"), fontsize=9
            )
            axes[row, col].set_xlabel("Same CpGs")
            axes[row, col].set_ylabel("Same molecules")
    fig.suptitle(
        p["gene"]
        + ": what combining the modifications removes from the read display\nTop: grey C, blue 5mC, gold 5hmC. Bottom: grey C, teal modified. Pale = missing. Same rows/columns.",
        fontsize=10,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P2_04_same_molecule_encoding")
    (encodings_S / "results/P2_ANSWER.md").write_text(
        "# P2: separate versus combined molecular signals\n\n"
        + encodings_json.dumps(summary, indent=2)
        + "\n\nAll three encodings use the same source molecules and physical CpG pairs. Primary comparisons require the same eligible locus/bin support in all three encodings, with separate row/column-preserving nulls. Covariance residuals retain prevalence-dependent bounds and cannot alone establish which modification carries the most information. The secondary phi comparison uses the same sufficiently polymorphic pairs with finite null phi in every encoding; limited support is explicit. It is full-pair phi, not the prior F2 12-molecule-subsampled phi. This comparison is descriptive and selected using 5mC/5hmC regional contrasts; it is not independent evidence of modification-specific tumour biology.\n"
    )
    _regional.regional_save(
        encodings_S / "results/P2_complete.json",
        {"analysis_and_figures_generated": True, "requires_review": True},
    )
    print(
        "P2 COMPLETE",
        encodings_json.dumps([r for r in summary if r["scope"] == "25_249bp"]),
        flush=True,
    )


def initialize_encodings():
    """Initialize the encodings stage once; load its declared inputs."""
    global encodings_COL, encodings_D, encodings_ENC, encodings_F, encodings_R, encodings_S
    if _runtime.initialized("strengthen_compare"):
        return
    _runtime.begin("strengthen_compare")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    encodings_R = encodings_Path(_config.workspace)
    encodings_D = encodings_R / ".analysis"
    encodings_F = encodings_D / "focused_followup"
    encodings_S = encodings_D / "strengthening"
    _molecules.blocks_S = encodings_S
    encodings_ENC = ["5hmC", "5mC", "combined"]
    encodings_COL = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    _runtime.finish("strengthen_compare")


def run_encodings():
    """Execute the encodings workflow stage."""
    initialize_encodings()
    encodings_run()


# RANDOMISATION


from pathlib import Path as randomisation_Path
import json as randomisation_json, csv as randomisation_csv, numpy as randomisation_np
import matplotlib.pyplot as randomisation_plt


def randomisation_patient(score, pairs, a, b):
    """Average matched differences over shared finite bins, then over pairs.

    At least six eligible promoter pairs are required. Each eligible pair
    receives equal weight, regardless of its number of available bins."""
    if not pairs:
        return (randomisation_np.nan, 0)
    c = score[[p["candidate_index"] for p in pairs], a : b + 1]
    r = score[[p["comparison_index"] for p in pairs], a : b + 1]
    v = c - r
    n = randomisation_np.isfinite(v).sum(1)
    ok = n > 0
    means = randomisation_np.nansum(v[ok], axis=1) / n[ok]
    return (float(means.mean()) if len(means) >= 6 else randomisation_np.nan, int(len(means)))


def randomisation_statistic(values, men):
    """Return an absolute Welch statistic for one diagnosis allocation.

    Requires at least three eligible meningiomas and ten glioblastomas.
    Unsupported allocations return zero and remain in the randomisation null."""
    mask = randomisation_np.zeros(20, bool)
    mask[men] = True
    m = values[mask & randomisation_np.isfinite(values)]
    g = values[~mask & randomisation_np.isfinite(values)]
    eligible = len(m) >= 3 and len(g) >= 10
    effect = float(m.mean() - g.mean()) if len(m) and len(g) else None
    if not eligible:
        return (0.0, effect, len(m), len(g), False)
    se = randomisation_np.sqrt(m.var(ddof=1) / len(m) + g.var(ddof=1) / len(g))
    if se == 0:
        return (float("inf") if effect else 0.0, effect, len(m), len(g), True)
    return (abs(effect) / se, effect, len(m), len(g), True)


def randomisation_run():
    assert (randomisation_S / "results/h_scoring_complete.json").exists()
    alloc = randomisation_json.loads((randomisation_S / "plans/assignments.json").read_text())
    assert (
        len(alloc) == 500
        and len({tuple(a["meningioma_indices"]) for a in alloc}) == 500
        and all((len(set(a["meningioma_indices"])) == 5 for a in alloc))
    )
    obs = randomisation_json.loads(
        (randomisation_S / "plans/permutation_folds/0000.json").read_text()
    )
    labels = [f["held_out"] for f in obs["folds"]]
    old = randomisation_json.loads((randomisation_F / "plans/folds.json").read_text())
    for f in obs["folds"]:
        prior = next((x for x in old if x["held_out"] == f["held_out"]))
        assert [(p["candidate_index"], p["comparison_index"]) for p in f["pairs"]] == [
            (p["candidate_index"], p["comparison_index"]) for p in prior["pairs"]
        ]
    nfeatures = len(
        randomisation_json.loads((randomisation_S / "plans/feature_universe.json").read_text())
    )
    scores = randomisation_np.full((20, nfeatures, 8), randomisation_np.nan)
    li = {l: i for (i, l) in enumerate(labels)}
    with (randomisation_S / "results/h_locus_bins.tsv").open() as fh:
        for r in randomisation_csv.DictReader(fh, delimiter="\t"):
            if r["eligible"] == "True":
                scores[li[r["study_label"]], int(r["universe_index"]), int(r["distance_bin"])] = (
                    float(r["residual_pp"])
                )
    allrows = []
    patientrows = []
    maxima = []
    oracle = 0
    for allocation in alloc:
        index = allocation["assignment_index"]
        aobj = randomisation_json.loads(
            (randomisation_S / "plans/permutation_folds" / f"{index:04d}.json").read_text()
        )
        assert aobj["meningioma_indices"] == allocation["meningioma_indices"]
        ar = []
        for a, b in randomisation_RANGES:
            values = []
            for fold in aobj["folds"]:
                i = fold["held_out_index"]
                (v, n) = randomisation_patient(scores[i], fold["pairs"], a, b)
                values.append(v)
                if index in [0, 1, 2] and (a, b) in [(0, 7), (2, 4), (3, 3)]:
                    brute = []
                    for pair in fold["pairs"]:
                        dd = [
                            scores[i, pair["candidate_index"], k]
                            - scores[i, pair["comparison_index"], k]
                            for k in range(a, b + 1)
                            if randomisation_np.isfinite(scores[i, pair["candidate_index"], k])
                            and randomisation_np.isfinite(scores[i, pair["comparison_index"], k])
                        ]
                        if dd:
                            brute.append(sum(dd) / len(dd))
                    expected = sum(brute) / len(brute) if len(brute) >= 6 else randomisation_np.nan
                    assert (
                        randomisation_np.isnan(v)
                        and randomisation_np.isnan(expected)
                        or randomisation_np.isclose(v, expected, atol=1e-12)
                    )
                    assert n == len(brute)
                    oracle += 1
                if index == 0:
                    patientrows.append(
                        {
                            "study_label": fold["held_out"],
                            "diagnosis": (
                                "Meningioma"
                                if i in allocation["meningioma_indices"]
                                else "Glioblastoma"
                            ),
                            "first_bin": a,
                            "last_bin": b,
                            "distance_from_bp": randomisation_EDGES[a],
                            "distance_to_bp": randomisation_EDGES[b + 1] - 1,
                            "eligible_pairs": n,
                            "eligible": bool(randomisation_np.isfinite(v)),
                            "residual_pp": float(v) if randomisation_np.isfinite(v) else None,
                        }
                    )
            (t, e, nm, ng, eligible) = randomisation_statistic(
                randomisation_np.asarray(values), allocation["meningioma_indices"]
            )
            ar.append(
                {
                    "assignment_index": index,
                    "first_bin": a,
                    "last_bin": b,
                    "distance_from_bp": randomisation_EDGES[a],
                    "distance_to_bp": randomisation_EDGES[b + 1] - 1,
                    "statistic": float(t),
                    "mean_MEN_minus_GBM_pp": e,
                    "eligible_MEN": nm,
                    "eligible_GBM": ng,
                    "eligible": eligible,
                }
            )
        allrows.extend(ar)
        maxima.append(
            {
                "assignment_index": index,
                "maximum_statistic": max((r["statistic"] for r in ar)),
                "eligible_ranges": sum((r["eligible"] for r in ar)),
                "selected_pairs": sum((len(f["pairs"]) for f in aobj["folds"])),
            }
        )
        if index % 50 == 0:
            print("STATISTICS", index, flush=True)
    null = randomisation_np.array([m["maximum_statistic"] for m in maxima[1:]])
    observed = allrows[:36]
    for row in observed:
        row["max_family_adjusted_p"] = (
            1 + int(randomisation_np.sum(null >= row["statistic"]))
        ) / 500
    omnibus = (1 + int(randomisation_np.sum(null >= maxima[0]["maximum_statistic"]))) / 500
    _regional.regional_table(randomisation_S / "results/P3_all_interval_statistics.tsv", allrows)
    _regional.regional_table(randomisation_S / "results/P3_observed_intervals.tsv", observed)
    _regional.regional_table(
        randomisation_S / "results/P3_observed_patient_scores.tsv", patientrows
    )
    _regional.regional_table(randomisation_S / "results/P3_assignment_maxima.tsv", maxima)
    summary = {
        "observed_max_statistic": maxima[0]["maximum_statistic"],
        "omnibus_p": omnibus,
        "random_assignments": 499,
        "p_resolution": 0.002,
        "random_assignments_no_eligible_range": sum(
            (m["eligible_ranges"] == 0 for m in maxima[1:])
        ),
        "random_assignments_with_eligible_range": sum(
            (m["eligible_ranges"] > 0 for m in maxima[1:])
        ),
        "observed_selected_pairs": maxima[0]["selected_pairs"],
        "best_observed_interval": max(observed, key=lambda r: r["statistic"]),
        "exploratory_25_249bp": next(
            (r for r in observed if r["first_bin"] == 2 and r["last_bin"] == 4)
        ),
        "all_distances": next((r for r in observed if r["first_bin"] == 0 and r["last_bin"] == 7)),
        "scope": "Entire selection, matching and residual-scoring procedure under diagnosis-label exchangeability; not an isolated spatial mechanism test",
        "family": "36 contiguous distance intervals; does not cover all historical study hypotheses",
    }
    _regional.regional_save(randomisation_S / "results/P3_summary.json", summary)
    _regional.regional_save(
        randomisation_S / "validation/P3_numerical_checks.json",
        {
            "unique_assignments": 500,
            "all_have_5_MEN_15_GBM": True,
            "observed_selections_match_previous_heldout_folds": True,
            "independent_patient_score_checks": oracle,
            "failed_assignments_retained": True,
            "all_500_maxima_present": len(maxima) == 500,
            "no_nan_statistics": bool(
                not randomisation_np.isnan([r["statistic"] for r in allrows]).any()
            ),
        },
    )
    (fig, ax) = randomisation_plt.subplots(figsize=(8, 5))
    finite = null[randomisation_np.isfinite(null)]
    ax.hist(finite, bins=30, color="#9aacb7", edgecolor="white")
    ax.axvline(
        maxima[0]["maximum_statistic"],
        c="#147D92",
        lw=2,
        label=f"Observed maximum; p={omnibus:.3f}",
    )
    ax.set_xlabel("Maximum absolute studentized contrast over 36 intervals")
    ax.set_ylabel("Random diagnosis assignments")
    ax.legend(loc="center right")
    ax.set_title("Full selection and matching repeated for each diagnosis assignment")
    ax.text(
        0.98,
        0.97,
        f"499 random assignments; {summary['random_assignments_no_eligible_range']} without an eligible interval\nFailures retained with statistic 0",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P3_01_selection_aware_randomisation")
    (fig, axes) = randomisation_plt.subplots(1, 2, figsize=(13, 6))
    mats = [
        randomisation_np.full((8, 8), randomisation_np.nan),
        randomisation_np.full((8, 8), randomisation_np.nan),
    ]
    for r in observed:
        mats[0][r["first_bin"], r["last_bin"]] = r["statistic"]
        mats[1][r["first_bin"], r["last_bin"]] = r["max_family_adjusted_p"]
    for ax, m, title in zip(
        axes,
        mats,
        ["Observed absolute studentized contrast", "P adjusted over all 36 distance intervals"],
    ):
        im = ax.imshow(
            randomisation_np.ma.masked_invalid(m),
            cmap="viridis" if ax is axes[0] else "viridis_r",
            vmin=0,
        )
        ax.set_xticks(range(8), [str(x - 1) for x in randomisation_EDGES[1:]], rotation=45)
        ax.set_yticks(range(8), randomisation_EDGES[:-1])
        ax.set_xlabel("Interval upper bound (bp)")
        ax.set_ylabel("Interval lower bound (bp)")
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.7)
        for a, b in randomisation_RANGES:
            ax.text(
                b,
                a,
                f"{m[a, b]:.2f}" if ax is axes[0] else f"{m[a, b]:.3f}",
                ha="center",
                va="center",
                fontsize=7,
                color=(
                    "white"
                    if ax is axes[0]
                    and m[a, b] < randomisation_np.nanmax(m) * 0.5
                    or (ax is axes[1] and m[a, b] > 0.5)
                    else "black"
                ),
            )
    fig.suptitle("Exploratory distance ranges assessed as one statistical family", fontsize=11)
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P3_02_distance_family")
    (fig, axes) = randomisation_plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, key, title in zip(
        axes,
        ["selected_pairs", "eligible_ranges"],
        [
            "Selected pairs across 20 held-out folds",
            "Distance intervals with sufficient patient support",
        ],
    ):
        ax.hist([m[key] for m in maxima[1:]], bins=20, color="#9aacb7", edgecolor="white")
        ax.axvline(maxima[0][key], c="#147D92", lw=2, label="Observed")
        ax.set_xlabel(title)
        ax.set_ylabel("Random assignments")
        ax.legend()
    fig.suptitle("Selection and support are part of the randomised procedure", fontsize=11)
    fig.tight_layout()
    _molecules.blocks_figure(fig, "P3_03_randomisation_support")
    answer = (
        "# P3: selection-aware diagnosis randomisation\n\n"
        + randomisation_json.dumps(summary, indent=2)
        + "\n\nFor each allocation, promoters and matched comparisons were selected anew in every patient-held-out training set. Molecule scores were cached by patient/locus because they do not depend on the allocation. All 36 contiguous unions of the eight distance bins contribute to the maximum-statistic family, including the previously explored 25–249 bp range. Assignments with insufficient support were retained with statistic zero. This is a test of the entire selection-and-scoring procedure under exchangeable diagnosis labels, conditional on the fixed all-specimen coverage universe and randomisation settings. It does not separate regional selection signal from additional spatial information, establish causality, validate an independent cohort, or correct every historical hypothesis explored in this project. The 499 sampled assignments limit p-value resolution to 0.002. Null chains preserve molecule counts, CpG counts and missingness; their two-chain checks cannot prove complete mixing.\n"
    )
    (randomisation_S / "results/P3_ANSWER.md").write_text(answer)
    _regional.regional_save(
        randomisation_S / "results/P3_complete.json",
        {"analysis_and_figures_generated": True, "requires_review": True},
    )
    print("P3 COMPLETE", randomisation_json.dumps(summary), flush=True)


def initialize_randomisation():
    """Initialize the randomisation stage once; load its declared inputs."""
    global randomisation_D, randomisation_EDGES, randomisation_F, randomisation_R, randomisation_RANGES, randomisation_S
    if _runtime.initialized("strengthen_permutation"):
        return
    _runtime.begin("strengthen_permutation")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    randomisation_R = randomisation_Path(_config.workspace)
    randomisation_D = randomisation_R / ".analysis"
    randomisation_F = randomisation_D / "focused_followup"
    randomisation_S = randomisation_D / "strengthening"
    _molecules.blocks_S = randomisation_S
    randomisation_EDGES = [1, 10, 25, 50, 100, 250, 500, 1000, 2500]
    randomisation_RANGES = [(a, b) for a in range(8) for b in range(a, 8)]
    _runtime.finish("strengthen_permutation")


def run_randomisation():
    """Execute the randomisation workflow stage."""
    initialize_randomisation()
    randomisation_run()
