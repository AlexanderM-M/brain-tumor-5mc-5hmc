"""Figures analysis stages for the brain tumour 5mC/5hmC study.

Each stage has explicit initialization and execution; importing is read-only.
See docs/workflow.md for inputs, stage order and reproduction limits.
"""

from . import runtime as _runtime

from .runtime import config as _config


# MAIN


from pathlib import Path as main_Path
import csv as main_csv, gzip as main_gzip, json as main_json, re as main_re
import numpy as main_np
import matplotlib as main_matplotlib
import matplotlib.pyplot as main_plt
from matplotlib.colors import (
    ListedColormap as main_ListedColormap,
    BoundaryNorm as main_BoundaryNorm,
)
from matplotlib.patches import Patch as main_Patch
from matplotlib.text import Text as main_Text


def main_table(path):
    p = main_R / path
    with main_gzip.open(p, "rt") if str(p).endswith(".gz") else p.open() as h:
        return list(main_csv.DictReader(h, delimiter="\t"))


def main_js(path):
    return main_json.loads((main_R / path).read_text())


def main_panel(ax, letter, title):
    ax.set_title(title, loc="left", pad=12)
    ax.text(
        -0.12, 1.08, letter, transform=ax.transAxes, fontweight="bold", fontsize=12, va="bottom"
    )


def main_save(fig, n):
    for t in fig.findobj(match=main_Text):
        assert not main_re.search("\\bN(?:20)?\\d{2}[._-]\\d+\\b", t.get_text(), main_re.I)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(main_OUT / f"Figure_{n}.{ext}", dpi=350)
    main_plt.close(fig)


def main_patient_points(ax, rows, key):
    for i, d in enumerate(main_COL):
        a = main_np.array([float(r[key]) for r in rows if r["diagnosis"] == d])
        j = main_np.linspace(-0.14, 0.14, len(a))
        ax.scatter(i + j, a, c=main_COL[d], s=18, alpha=0.85, zorder=3)
        ax.plot([i - 0.2, i + 0.2], [a.mean()] * 2, c="black", lw=1.5)
    ax.set_xticks([0, 1], ["GBM", "MEN"])
    ax.grid(axis="y", alpha=0.15)


def main_paired(ax, rows, keys, ticks):
    for r in rows:
        ax.plot(
            range(len(keys)),
            [float(r[k]) for k in keys],
            "-o",
            c=main_COL[r["diagnosis"]],
            alpha=0.6,
            lw=0.8,
            ms=3,
        )
    ax.axhline(0, c="#999999", lw=0.6)
    ax.set_xticks(range(len(keys)), ticks)
    ax.grid(axis="y", alpha=0.15)


def initialize_main():
    """Initialize the main stage once; load its declared inputs."""
    global main_COL, main_OUT, main_R, main_a, main_ax, main_axes, main_bottom, main_c, main_colors, main_d, main_enc, main_example, main_examples, main_f1, main_fig, main_folds, main_frequency, main_gs, main_h, main_i, main_im, main_j, main_k, main_key, main_label, main_labels, main_lookup, main_m, main_mat, main_means, main_mids, main_norm, main_null, main_observed, main_pan, main_r, main_row, main_rows, main_rr, main_sel, main_spec, main_ss, main_title, main_v, main_win, main_x, main_y, main_z
    if _runtime.initialized("manuscript_results_figures"):
        return
    _runtime.begin("manuscript_results_figures")
    main_matplotlib.use("Agg")
    main_R = main_Path(_config.workspace)
    main_OUT = main_R / "figures/main/manuscript"
    main_OUT.mkdir(exist_ok=True)
    main_plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )
    main_COL = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    main_spec = main_table("results/regional/tables/specimen_summary.tsv")
    main_win = main_table("results/regional/tables/pooled_100kb_regions.tsv.gz")
    main_fig = main_plt.figure(figsize=(7.2, 7.7))
    main_gs = main_fig.add_gridspec(3, 6, height_ratios=[1, 1.45, 1.2], hspace=0.72, wspace=1.35)
    for main_k, (main_key, main_title) in enumerate(
        [
            ("common_CpG_5mC_percent", "Genome-wide 5mC"),
            ("common_CpG_5hmC_percent", "Genome-wide 5hmC"),
            ("common_CpG_combined_percent", "Genome-wide combined"),
        ]
    ):
        main_ax = main_fig.add_subplot(main_gs[0, 2 * main_k : 2 * main_k + 2])
        main_patient_points(main_ax, main_spec, main_key)
        main_panel(main_ax, chr(65 + main_k), main_title)
        main_ax.set_ylabel("Mean fraction (%)")
    main_m = main_np.array([float(r["GBM_minus_meningioma_5mC_pp"]) for r in main_win])
    main_h = main_np.array([float(r["GBM_minus_meningioma_5hmC_pp"]) for r in main_win])
    main_sel = main_np.array(
        [r["opposing_changes_with_small_combined_difference"] == "True" for r in main_win]
    )
    assert len(main_win) == 26585 and main_sel.sum() == 593
    for main_k, (main_x, main_label, main_title) in enumerate(
        [
            (main_m, "5mC difference (pp)", "Separate modifications"),
            (main_m + main_h, "Combined difference (pp)", "Combined modification"),
        ]
    ):
        main_ax = main_fig.add_subplot(main_gs[1, 3 * main_k : 3 * main_k + 3])
        main_ax.scatter(main_x, main_h, s=2, alpha=0.12, c="#667785", rasterized=True)
        main_ax.scatter(
            main_x[main_sel],
            main_h[main_sel],
            s=4,
            c="#D79A19",
            alpha=0.8,
            rasterized=True,
            label="Opposing-change windows",
        )
        main_ax.axhline(0, c="gray", lw=0.6)
        main_ax.axvline(0, c="gray", lw=0.6)
        main_ax.set_xlabel(main_label)
        main_ax.set_ylabel("5hmC difference (pp)")
        main_panel(main_ax, chr(68 + main_k), main_title)
        if main_k:
            main_ax.legend(loc="lower left", frameon=False, fontsize=6.5)
    main_examples = main_js("results/regional/tables/interpretation_summary.json")[
        "figure_examples"
    ]
    for main_k, main_r in enumerate(main_examples):
        main_ax = main_fig.add_subplot(main_gs[2, 2 * main_k : 2 * main_k + 2])
        main_v = [
            float(main_r["GBM_minus_meningioma_" + s + "_pp"]) for s in ["5mC", "5hmC", "combined"]
        ]
        main_ax.bar(range(3), main_v, color=["#497ba6", "#df9d29", "#888888"], width=0.65)
        main_ax.axhline(0, c="gray", lw=0.6)
        main_ax.set_xticks(range(3), ["5mC", "5hmC", "Sum"])
        main_ax.set_ylim(-17, 18)
        main_ax.set_ylabel("GBM − MEN (pp)")
        for main_j, main_y in enumerate(main_v):
            main_ax.text(
                main_j,
                main_y + (1 if main_y >= 0 else -1),
                f"{main_y:+.2f}",
                ha="center",
                va="bottom" if main_y >= 0 else "top",
                fontsize=7,
            )
        main_panel(
            main_ax, chr(70 + main_k), main_r["gene"] + "\n" + main_r["feature"].replace("_", " ")
        )
    main_fig.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.07)
    main_fig.text(
        0.5,
        0.988,
        "Regional composition: separate and combined modifications",
        ha="center",
        va="top",
        fontsize=11,
        fontweight="bold",
    )
    main_save(main_fig, 1)
    main_labels = [r["study_label"] for r in main_spec]
    main_colors = main_ListedColormap(["#eeeeee", "#cfcfcf", "#497ba6", "#df9d29"])
    main_norm = main_BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], 4)
    main_fig = main_plt.figure(figsize=(7.2, 7.6))
    main_gs = main_fig.add_gridspec(
        6, 5, height_ratios=[1, 1, 1, 1, 0.22, 2.25], hspace=0.55, wspace=0.25
    )
    for main_i, main_label in enumerate(main_labels):
        main_z = main_np.load(main_R / f"data/molecules/molecular/{main_label}/C02.npz")
        main_x = main_z["states"][:, main_z["common_mask"]]
        main_x = main_x[(main_x >= 0).sum(1) >= 8]
        if len(main_x) > 60:
            main_x = main_x[
                main_np.sort(
                    main_np.random.default_rng(20260912).choice(len(main_x), 60, replace=False)
                )
            ]
        main_ax = main_fig.add_subplot(main_gs[main_i // 5, main_i % 5])
        main_ax.imshow(
            main_x, aspect="auto", interpolation="nearest", cmap=main_colors, norm=main_norm
        )
        main_ax.set_xticks([])
        main_ax.set_yticks([])
        main_ax.set_title(
            main_label, color=main_COL[main_spec[main_i]["diagnosis"]], fontsize=8, pad=3
        )
    main_fig.text(0.065, 0.95, "A", fontsize=12, fontweight="bold")
    main_fig.text(0.105, 0.951, "P2RX1 promoter: individual molecules at common CpGs", fontsize=9)
    main_fig.legend(
        handles=[
            main_Patch(color=c, label=l)
            for (c, l) in zip(main_colors.colors, ["Missing", "C", "5mC", "5hmC"])
        ],
        loc="upper center",
        bbox_to_anchor=(0.53, 0.38),
        ncol=4,
        frameon=False,
    )
    main_bottom = main_gs[5, :].subgridspec(1, 3, wspace=0.7)
    main_f1 = main_table("results/followup/tables/F1_specimen_scores.tsv")
    for main_k, main_d in enumerate(main_COL):
        main_ax = main_fig.add_subplot(main_bottom[0, main_k])
        main_ss = [r for r in main_f1 if r["diagnosis"] == main_d]
        main_paired(
            main_ax,
            main_ss,
            ["original_candidate_minus_comparison_pp", "heldout_candidate_minus_comparison_pp"],
            ["Original", "Held-out"],
        )
        main_panel(main_ax, chr(66 + main_k), main_d)
        main_ax.set_ylabel("Candidate − comparison (pp)")
    main_ax = main_fig.add_subplot(main_bottom[0, 2])
    main_patient_points(
        main_ax,
        main_table("results/followup/tables/F1_blind_panel_scores.tsv"),
        "mean_joint_5hmC_excess_pp",
    )
    main_panel(main_ax, "D", "Diagnosis-blind panel")
    main_ax.set_ylabel("Excess co-occurrence (pp)")
    main_fig.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.07)
    main_fig.text(
        0.5,
        0.992,
        "From regional selection to single-molecule patterns",
        ha="center",
        va="top",
        fontsize=11,
        fontweight="bold",
    )
    main_save(main_fig, 2)
    (main_fig, main_axes) = main_plt.subplots(3, 2, figsize=(7.2, 8.7))
    main_fig.subplots_adjust(left=0.11, right=0.98, bottom=0.08, top=0.92, hspace=0.86, wspace=0.48)
    main_ss = [
        r
        for r in main_table("results/followup/tables/F2_matched_specimen_scores.tsv")
        if r["panel"] == "heldout"
        and r["sensitivity"] == "standard"
        and (r["metric"] == "joint_excess_pp")
        and (r["eligible"] == "True")
    ]
    main_ax = main_axes[0, 0]
    main_paired(
        main_ax, main_ss, ["unmatched_difference", "matched_difference"], ["Unmatched", "Matched"]
    )
    main_panel(main_ax, "A", "Distance / abundance matching")
    main_ax.set_ylabel("Candidate − comparison (pp)")
    main_ss = [
        r
        for r in main_table("results/followup/tables/F2_row_null_specimen_scores.tsv")
        if r["panel"] == "heldout" and r["eligible"] == "True"
    ]
    main_ax = main_axes[0, 1]
    main_paired(
        main_ax,
        main_ss,
        [
            "observed_paired_difference_pp",
            "null_expected_paired_difference_pp",
            "residual_paired_difference_pp",
        ],
        ["Observed", "Null", "Residual"],
    )
    main_panel(main_ax, "B", "Read- and CpG-count control")
    main_ax.set_ylabel("Candidate − comparison (pp)")
    main_rows = [
        r
        for r in main_table("results/followup/tables/F2_row_null_distance_curves.tsv")
        if r["panel"] == "heldout"
    ]
    main_mids = main_np.sqrt(
        main_np.array([1, 10, 25, 50, 100, 250, 500, 1000])
        * main_np.array([10, 25, 50, 100, 250, 500, 1000, 2500])
    )
    main_ax = main_axes[1, 0]
    for main_d in main_COL:
        main_means = []
        for main_k in range(8):
            main_v = [
                float(r["residual_paired_difference_pp"])
                for r in main_rows
                if r["diagnosis"] == main_d
                and int(r["distance_bin"]) == main_k
                and (r["eligible"] == "True")
            ]
            main_means.append(main_np.mean(main_v) if main_v else main_np.nan)
        main_ax.plot(main_mids, main_means, "-o", c=main_COL[main_d], ms=4, label=main_d)
    main_ax.set_xscale("log")
    main_ax.axhline(0, c="gray", lw=0.6)
    main_ax.set_xlabel("CpG separation (bp; bin midpoint)")
    main_ax.set_ylabel("Mean paired residual (pp)")
    main_panel(main_ax, "C", "Residual by CpG distance")
    main_ax.legend(frameon=False, fontsize=6.5)
    main_ax = main_axes[1, 1]
    main_rr = [
        r
        for r in main_table("results/strengthening/tables/P1_patient_influence.tsv")
        if r["scope"] == "25_249bp"
    ]
    for main_i, main_r in enumerate(main_rr):
        main_c = main_COL[main_r["diagnosis"]]
        main_ax.plot(
            [main_i, main_i],
            [float(main_r["minimum_after_removal_pp"]), float(main_r["maximum_after_removal_pp"])],
            c=main_c,
            lw=1.6,
        )
        main_ax.scatter(main_i, float(main_r["baseline_pp"]), c=main_c, s=13, zorder=3)
    main_ax.axhline(0, c="gray", lw=0.6)
    main_ax.set_xticks(range(20), [r["study_label"] for r in main_rr], rotation=90, fontsize=5.5)
    main_ax.set_ylabel("Paired residual (pp)")
    main_panel(main_ax, "D", "Promoter-pair removal: 25–249 bp")
    main_rr = main_table("results/strengthening/tables/P3_assignment_maxima.tsv")
    main_null = main_np.array([float(r["maximum_statistic"]) for r in main_rr[1:]])
    main_observed = float(main_rr[0]["maximum_statistic"])
    assert len(main_null) == 499
    main_ax = main_axes[2, 0]
    main_ax.hist(main_null, bins=25, color="#9aacb7", edgecolor="white")
    main_ax.axvline(main_observed, c=main_COL["Meningioma"], lw=2)
    main_ax.set_xlabel("Maximum absolute Welch statistic")
    main_ax.set_ylabel("Random label allocations")
    main_panel(main_ax, "E", "Whole-procedure randomisation")
    main_ax.set_title("Whole-procedure randomisation", loc="left", y=1.17, pad=12)
    main_ax.texts[0].set_position((-0.12, 1.25))
    main_ax.text(
        0,
        1.015,
        "Omnibus p = 0.004\n250/499 allocations without eligible intervals",
        transform=main_ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.5,
        clip_on=False,
    )
    main_rr = main_table("results/strengthening/tables/P3_observed_intervals.tsv")
    main_mat = main_np.full((8, 8), main_np.nan)
    for main_r in main_rr:
        main_mat[int(main_r["first_bin"]), int(main_r["last_bin"])] = float(
            main_r["max_family_adjusted_p"]
        )
    main_ax = main_axes[2, 1]
    main_im = main_ax.imshow(
        main_np.ma.masked_invalid(main_mat), vmin=0, vmax=1, cmap="viridis_r", aspect="auto"
    )
    main_ax.set_xticks(range(8), [9, 24, 49, 99, 249, 499, 999, 2499], rotation=60)
    main_ax.set_yticks(range(8), [1, 10, 25, 50, 100, 250, 500, 1000])
    main_ax.set_xlabel("Interval end (bp)")
    main_ax.set_ylabel("Interval start (bp)")
    main_panel(main_ax, "F", "Distance-family adjusted p")
    for main_i in range(8):
        for main_j in range(main_i, 8):
            main_ax.text(
                main_j,
                main_i,
                f"{main_mat[main_i, main_j]:.3f}",
                ha="center",
                va="center",
                fontsize=4.5,
                color="white" if main_mat[main_i, main_j] > 0.5 else "black",
            )
    main_fig.text(
        0.5,
        0.985,
        "Abundance controls, promoter influence and statistical support",
        ha="center",
        va="top",
        fontsize=10.5,
        fontweight="bold",
    )
    main_save(main_fig, 3)
    main_fig = main_plt.figure(figsize=(7.2, 7.3))
    main_gs = main_fig.add_gridspec(3, 4, height_ratios=[2.4, 1, 1], hspace=0.8, wspace=0.35)
    main_enc = ["5hmC", "5mC", "combined"]
    main_rr = main_table("results/strengthening/tables/P2_specimen_comparison.tsv")
    for main_k, main_d in enumerate(main_COL):
        main_ax = main_fig.add_subplot(main_gs[0, 2 * main_k : 2 * main_k + 2])
        for main_label in [r["study_label"] for r in main_spec if r["diagnosis"] == main_d]:
            main_a = [
                next(
                    (
                        r
                        for r in main_rr
                        if r["study_label"] == main_label
                        and r["scope"] == "25_249bp"
                        and (r["support"] == "all_pairs")
                        and (r["encoding"] == e)
                    )
                )
                for e in main_enc
            ]
            assert all((r["eligible"] == "True" for r in main_a))
            main_ax.plot(
                range(3),
                [float(r["residual_pp"]) for r in main_a],
                "-o",
                c=main_COL[main_d],
                alpha=0.7,
                ms=4,
                lw=1,
            )
        main_ax.axhline(0, c="gray", lw=0.6)
        main_ax.set_xticks(range(3), ["5hmC", "5mC", "Combined"])
        main_ax.set_ylabel("Candidate − comparison residual (pp)")
        main_panel(main_ax, chr(65 + main_k), main_d + " (25–249 bp)")
    main_pan = main_js("data/plans/followup/panel.json")
    main_folds = main_js("data/plans/followup/folds.json")
    main_frequency = {
        p["locus_id"]: sum(
            (
                any((q["candidate_locus_id"] == p["locus_id"] for q in f["pairs"]))
                for f in main_folds
            )
        )
        for p in main_pan
    }
    main_lookup = {p["locus_id"]: p for p in main_pan}
    main_key = min(main_frequency, key=lambda k: (-main_frequency[k], main_lookup[k]["gene"], k))
    main_example = main_lookup[main_key]["gene"]
    for main_k, main_label in enumerate(["GBM-01", "GBM-02", "MEN-01", "MEN-02"]):
        main_z = main_np.load(main_R / f"data/molecules/followup/{main_label}/{main_key}.npz")
        main_x = main_z["states"][:, main_z["common_mask"]]
        main_x = main_x[(main_x >= 0).sum(1) >= 10]
        if len(main_x) > 50:
            main_x = main_x[
                main_np.sort(
                    main_np.random.default_rng(20260912).choice(len(main_x), 50, replace=False)
                )
            ]
        for main_row in [1, 2]:
            main_ax = main_fig.add_subplot(main_gs[main_row, main_k])
            if main_row == 1:
                main_ax.imshow(
                    main_x, aspect="auto", interpolation="nearest", cmap=main_colors, norm=main_norm
                )
            else:
                main_ax.imshow(
                    main_np.where(main_x < 0, -1, main_np.where(main_x > 0, 1, 0)),
                    aspect="auto",
                    interpolation="nearest",
                    cmap=main_ListedColormap(["#eeeeee", "#cfcfcf", "#147D92"]),
                    vmin=-1,
                    vmax=1,
                )
            main_ax.set_xticks([])
            main_ax.set_yticks([])
            main_ax.set_title(main_label, fontsize=8, pad=3)
            if main_k == 0:
                main_ax.set_ylabel("Separate" if main_row == 1 else "Combined")
            if main_row == 2:
                main_ax.set_xlabel("Same CpGs", fontsize=7)
    main_fig.text(0.05, 0.545, "C", fontsize=12, fontweight="bold")
    main_fig.text(
        0.09,
        0.548,
        main_example + ": the same molecules before and after combining states",
        fontsize=8.5,
    )
    main_fig.legend(
        handles=[
            main_Patch(color=c, label=l)
            for (c, l) in zip(
                ["#eeeeee", "#cfcfcf", "#497ba6", "#df9d29", "#147D92"],
                ["Missing", "C", "5mC", "5hmC", "Either modification"],
            )
        ],
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.52, 0.025),
        fontsize=7,
    )
    main_fig.subplots_adjust(left=0.11, right=0.98, top=0.9, bottom=0.13)
    main_fig.text(
        0.5,
        0.985,
        "Molecular organisation is not exclusive to 5hmC",
        ha="center",
        va="top",
        fontsize=11,
        fontweight="bold",
    )
    main_save(main_fig, 4)
    (main_OUT / "figure_metadata.json").write_text(
        main_json.dumps(
            {
                "source": "frozen numerical tables and read-linked matrices",
                "new_extraction": False,
                "new_statistical_tests": False,
                "figure_4_example_gene": main_example,
                "figure_4_example_locus": main_key,
                "labels": main_labels,
                "figure3_distance_eligible_counts": {
                    d: [
                        sum(
                            (
                                r["diagnosis"] == d
                                and int(r["distance_bin"]) == k
                                and (r["eligible"] == "True")
                                for r in main_rows
                            )
                        )
                        for k in range(8)
                    ]
                    for d in main_COL
                },
            },
            indent=2,
        )
        + "\n"
    )
    print("Rendered four figures; Figure 4 example:", main_example, flush=True)
    _runtime.finish("manuscript_results_figures")


def run_main():
    """Execute the main workflow stage."""
    initialize_main()


# SUPPLEMENTARY


from pathlib import Path as supplementary_Path
import csv as supplementary_csv, gzip as supplementary_gzip, json as supplementary_json, re as supplementary_re, io as supplementary_io, html as supplementary_html
import numpy as supplementary_np
import matplotlib as supplementary_matplotlib
import matplotlib.pyplot as supplementary_plt
from matplotlib.lines import Line2D as supplementary_Line2D
from matplotlib.text import Text as supplementary_Text
from reportlab.pdfgen import canvas as supplementary_canvas
from reportlab.lib.pagesizes import A4 as supplementary_A4
from reportlab.platypus import Paragraph as supplementary_Paragraph
from reportlab.lib.styles import ParagraphStyle as supplementary_ParagraphStyle
from reportlab.pdfbase import pdfmetrics as supplementary_pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont as supplementary_TTFont
from pypdf import (
    PdfReader as supplementary_PdfReader,
    PdfWriter as supplementary_PdfWriter,
    Transformation as supplementary_Transformation,
)


def supplementary_tab(p):
    p = supplementary_R / p
    with supplementary_gzip.open(p, "rt") if p.suffix == ".gz" else p.open() as f:
        return list(supplementary_csv.DictReader(f, delimiter="\t"))


def supplementary_num(r, k):
    return float(r[k]) if r[k] not in ("", "None") else supplementary_np.nan


def supplementary_panel(ax, l, title):
    ax.set_title(title, loc="left", pad=12)
    ax.text(
        -26 / (ax.get_position().width * ax.figure.get_figwidth() * 72),
        1.06,
        l,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
    )


def supplementary_diagkey(fig, y=0.99):
    fig.legend(
        handles=[
            supplementary_Line2D([], [], marker="o", ls="", color=c, label=d, ms=4)
            for (d, c) in supplementary_COL.items()
        ],
        loc="upper center",
        bbox_to_anchor=(0.55, y),
        ncol=2,
        frameon=False,
    )


def supplementary_patientaxis(ax):
    ax.set_xticks(range(20), supplementary_labels, rotation=90)
    ax.set_xlim(-0.7, 19.7)
    ax.axvline(14.5, color="#dddddd", lw=0.7)


def supplementary_zero(ax):
    ax.axhline(0, color="#aaaaaa", lw=0.7, zorder=0)


def supplementary_save(fig, n):
    texts = [t.get_text() for t in fig.findobj(supplementary_Text)]
    assert not any(
        (
            supplementary_re.search("\\bN(?:20)?\\d{2}[._-]\\d+\\b", s, supplementary_re.I)
            for s in texts
        )
    )
    supplementary_figtexts[str(n)] = texts
    fig.canvas.draw()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(supplementary_O / f"Figure_S{n}.{ext}", dpi=350)
    supplementary_plt.close(fig)


def initialize_supplementary():
    """Initialize the supplementary stage once; load its declared inputs."""
    global supplementary_COL, supplementary_H, supplementary_O, supplementary_R, supplementary_W, supplementary__, supplementary_assoc, supplementary_ax, supplementary_axis, supplementary_axs, supplementary_b, supplementary_buf, supplementary_c, supplementary_caps, supplementary_cax, supplementary_cb, supplementary_cv, supplementary_diag, supplementary_display, supplementary_e, supplementary_eff, supplementary_f, supplementary_fh, supplementary_fig, supplementary_figtexts, supplementary_fontdir, supplementary_fp, supplementary_fr, supplementary_fw, supplementary_gene, supplementary_gs, supplementary_hi, supplementary_i, supplementary_im, supplementary_j, supplementary_k, supplementary_key, supplementary_l, supplementary_label, supplementary_labels, supplementary_layout, supplementary_legend, supplementary_letter, supplementary_lh, supplementary_lo, supplementary_lookup, supplementary_lp, supplementary_ly, supplementary_m, supplementary_margin, supplementary_mat, supplementary_maxheight, supplementary_md, supplementary_metric, supplementary_n, supplementary_name, supplementary_names, supplementary_out, supplementary_p, supplementary_page, supplementary_ph, supplementary_pw, supplementary_px, supplementary_py, supplementary_q, supplementary_r, supplementary_rows, supplementary_rr, supplementary_scale, supplementary_sel, supplementary_state, supplementary_style, supplementary_support, supplementary_t, supplementary_th, supplementary_title, supplementary_title_style, supplementary_ty, supplementary_v, supplementary_vals, supplementary_variant, supplementary_width, supplementary_writer, supplementary_xlabel, supplementary_xp, supplementary_z
    if _runtime.initialized("manuscript_supplementary_figures"):
        return
    _runtime.begin("manuscript_supplementary_figures")
    supplementary_matplotlib.use("Agg")
    supplementary_R = supplementary_Path(_config.workspace)
    supplementary_O = supplementary_R / "supplementary/figures/manuscript"
    supplementary_O.mkdir(exist_ok=True)
    supplementary_plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )
    supplementary_COL = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    supplementary_labels = [f"GBM-{i:02}" for i in range(1, 16)] + [
        f"MEN-{i:02}" for i in range(1, 6)
    ]
    supplementary_figtexts = {}
    supplementary_rows = supplementary_tab("results/regional/tables/pooled_100kb_regions.tsv.gz")
    supplementary_z = supplementary_np.load(supplementary_R / "data/regional/pooled_profiles.npz")
    supplementary_b = supplementary_z["windows5"][:, supplementary_z["valid5"], :]
    assert (
        supplementary_b.shape == (20, len(supplementary_rows), 2)
        and supplementary_z["labels"].tolist() == supplementary_labels
    )
    supplementary_sel = sorted(
        [
            i
            for (i, r) in enumerate(supplementary_rows)
            if r["exploratory_5hmC_candidate"] == "True"
        ],
        key=lambda i: -abs(
            supplementary_num(supplementary_rows[i], "GBM_minus_meningioma_5hmC_pp")
        ),
    )[:40]
    supplementary_mat = supplementary_b[:, supplementary_sel, 1].T.copy()
    supplementary_mat -= supplementary_mat.mean(1)[:, None]
    assert supplementary_np.allclose(
        supplementary_b[:15, :, 1].mean(0) - supplementary_b[15:, :, 1].mean(0),
        [supplementary_num(r, "GBM_minus_meningioma_5hmC_pp") for r in supplementary_rows],
        atol=1e-05,
    )
    (supplementary_fig, supplementary_ax) = supplementary_plt.subplots(figsize=(7.2, 6.5))
    supplementary_fig.subplots_adjust(left=0.3, right=0.86, bottom=0.13, top=0.96)
    supplementary_im = supplementary_ax.imshow(
        supplementary_mat,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-abs(supplementary_mat).max(),
        vmax=abs(supplementary_mat).max(),
    )
    supplementary_patientaxis(supplementary_ax)
    supplementary_ax.set_yticks(
        range(40), [supplementary_rows[i]["region"] for i in supplementary_sel], fontsize=6
    )
    supplementary_ax.set_xlabel("Specimen")
    supplementary_ax.set_ylabel("100-kb window (hg38)")
    supplementary_cb = supplementary_fig.colorbar(
        supplementary_im, ax=supplementary_ax, fraction=0.035, pad=0.035
    )
    supplementary_cb.set_label("5hmC deviation from window mean (pp)")
    supplementary_save(supplementary_fig, 1)
    supplementary_r = supplementary_tab("results/followup/tables/F2_matched_specimen_scores.tsv")
    supplementary_lookup = {
        (x["study_label"], x["sensitivity"], x["metric"]): x
        for x in supplementary_r
        if x["panel"] == "heldout"
    }
    (supplementary_fig, supplementary_axs) = supplementary_plt.subplots(2, 2, figsize=(7.2, 6.2))
    supplementary_fig.subplots_adjust(
        left=0.1, right=0.98, bottom=0.1, top=0.89, wspace=0.42, hspace=0.8
    )
    supplementary_diagkey(supplementary_fig)
    for supplementary_ax, supplementary_metric, supplementary_letter, supplementary_title in zip(
        supplementary_axs[0],
        ["joint_excess_pp", "phi"],
        "AB",
        ["Excess co-occurrence", "Frequency-normalised association"],
    ):
        for supplementary_i, supplementary_label in enumerate(supplementary_labels):
            supplementary_q = supplementary_lookup[
                supplementary_label, "standard", supplementary_metric
            ]
            if supplementary_q["eligible"] == "True":
                supplementary_vals = [
                    supplementary_num(supplementary_q, "unmatched_difference"),
                    supplementary_num(supplementary_q, "matched_difference"),
                ]
                supplementary_c = supplementary_COL[supplementary_q["diagnosis"]]
                supplementary_ax.plot(
                    [supplementary_i - 0.16, supplementary_i + 0.16],
                    supplementary_vals,
                    c=supplementary_c,
                    lw=0.7,
                    alpha=0.7,
                )
                supplementary_ax.scatter(
                    [supplementary_i - 0.16, supplementary_i + 0.16],
                    supplementary_vals,
                    c=supplementary_c,
                    s=[8, 22],
                    zorder=3,
                )
        supplementary_patientaxis(supplementary_ax)
        supplementary_zero(supplementary_ax)
        supplementary_ax.set_ylabel(
            "Candidate − comparison (pp)"
            if supplementary_metric == "joint_excess_pp"
            else "Candidate − comparison φ"
        )
        supplementary_panel(supplementary_ax, supplementary_letter, supplementary_title)
    for supplementary_ax, supplementary_diag, supplementary_letter in zip(
        supplementary_axs[1], supplementary_COL, "CD"
    ):
        for supplementary_label in supplementary_labels:
            if (
                supplementary_lookup[supplementary_label, "standard", "phi"]["diagnosis"]
                != supplementary_diag
            ):
                continue
            supplementary_vals = []
            for supplementary_variant in ["standard", "MAPQ60", "distance_at_least_25bp"]:
                supplementary_q = supplementary_lookup[
                    supplementary_label, supplementary_variant, "phi"
                ]
                supplementary_vals.append(
                    supplementary_num(supplementary_q, "matched_difference")
                    if supplementary_q["eligible"] == "True"
                    else supplementary_np.nan
                )
            supplementary_ax.plot(
                range(3),
                supplementary_vals,
                "-o",
                c=supplementary_COL[supplementary_diag],
                lw=0.8,
                alpha=0.75,
                ms=3,
            )
        supplementary_zero(supplementary_ax)
        supplementary_ax.set_xticks(range(3), ["Standard", "MAPQ ≥60", "Distance\n≥25 bp"])
        supplementary_ax.set_ylabel("Matched candidate − comparison φ")
        supplementary_panel(supplementary_ax, supplementary_letter, supplementary_diag)
    supplementary_save(supplementary_fig, 2)
    supplementary_r = supplementary_tab("results/strengthening/tables/P3_assignment_maxima.tsv")
    assert (
        len(supplementary_r) == 500
        and sum((int(x["eligible_ranges"]) == 0 for x in supplementary_r[1:])) == 250
    )
    (supplementary_fig, supplementary_axs) = supplementary_plt.subplots(1, 2, figsize=(7.2, 3.7))
    supplementary_fig.subplots_adjust(left=0.1, right=0.98, bottom=0.19, top=0.86, wspace=0.35)
    for (
        supplementary_ax,
        supplementary_key,
        supplementary_letter,
        supplementary_title,
        supplementary_xlabel,
    ) in zip(
        supplementary_axs,
        ["selected_pairs", "eligible_ranges"],
        "AB",
        ["Promoter selection", "Patient support"],
        ["Selected pairs across 20 folds", "Eligible distance intervals"],
    ):
        supplementary_ax.hist(
            [int(x[supplementary_key]) for x in supplementary_r[1:]],
            bins=20,
            color="#9aacb7",
            edgecolor="white",
        )
        supplementary_ax.axvline(
            int(supplementary_r[0][supplementary_key]), c="#147D92", lw=1.8, label="Observed"
        )
        supplementary_ax.set_xlabel(supplementary_xlabel)
        supplementary_ax.set_ylabel("Random allocations")
        supplementary_ax.legend(frameon=False, loc="upper center")
        supplementary_panel(supplementary_ax, supplementary_letter, supplementary_title)
    supplementary_save(supplementary_fig, 3)
    supplementary_r = supplementary_tab("results/strengthening/tables/P2_specimen_comparison.tsv")
    supplementary_lookup = {
        (x["study_label"], x["scope"], x["support"]): x
        for x in supplementary_r
        if x["encoding"] == "5hmC"
    }
    (supplementary_fig, supplementary_axs) = supplementary_plt.subplots(2, 1, figsize=(7.2, 4.6))
    supplementary_fig.subplots_adjust(left=0.14, right=0.86, bottom=0.28, top=0.9, hspace=1.0)
    for supplementary_ax, supplementary_support, supplementary_letter, supplementary_title in zip(
        supplementary_axs,
        ["all_pairs", "joint_informative"],
        "AB",
        ["Shared CpG-pair support", "Shared informative support for φ"],
    ):
        supplementary_m = supplementary_np.array(
            [
                [
                    int(
                        supplementary_lookup[label, scope, supplementary_support][
                            "matched_promoter_pairs"
                        ]
                    )
                    for label in supplementary_labels
                ]
                for scope in ["all_distances", "25_249bp"]
            ]
        )
        supplementary_im = supplementary_ax.pcolormesh(
            supplementary_np.arange(21) - 0.5,
            supplementary_np.arange(3) - 0.5,
            supplementary_m,
            cmap="viridis",
            vmin=0,
            vmax=12,
            shading="flat",
        )
        supplementary_ax.set_ylim(1.5, -0.5)
        supplementary_patientaxis(supplementary_ax)
        supplementary_ax.set_yticks([0, 1], ["1–2,499 bp", "25–249 bp"])
        supplementary_panel(supplementary_ax, supplementary_letter, supplementary_title)
        for supplementary_j in range(2):
            for supplementary_i in range(20):
                supplementary_ax.text(
                    supplementary_i,
                    supplementary_j,
                    str(supplementary_m[supplementary_j, supplementary_i]),
                    ha="center",
                    va="center",
                    fontsize=6,
                    color=(
                        "black"
                        if supplementary_m[supplementary_j, supplementary_i] > 7
                        else "white"
                    ),
                )
    supplementary_cax = supplementary_fig.add_axes([0.9, 0.28, 0.017, 0.6])
    supplementary_fig.colorbar(
        supplementary_im,
        cax=supplementary_cax,
        ticks=[0, 3, 6, 9, 12],
        label="Matched promoter pairs",
    )
    supplementary_fig.text(0.14, 0.09, "Eligible for shared φ at 25–249 bp", fontsize=8)
    supplementary_fig.text(
        0.14, 0.045, "Glioblastoma: 0 / 15     Meningioma: 0 / 5", fontsize=8, fontweight="bold"
    )
    assert all(
        (
            supplementary_lookup[l, "25_249bp", "joint_informative"]["eligible"] == "False"
            for l in supplementary_labels
        )
    )
    supplementary_save(supplementary_fig, 4)
    supplementary_r = supplementary_tab("supplementary/results/molecular/Q5_specimen_context.tsv")
    supplementary_eff = supplementary_tab("supplementary/results/molecular/Q5_context_effects.tsv")
    for supplementary_n, supplementary_axis in [(5, "gene"), (6, "CpG")]:
        supplementary_names = list(
            dict.fromkeys(
                (x["context"] for x in supplementary_r if x["context_axis"] == supplementary_axis)
            )
        )
        supplementary_display = [
            x.replace("Gene body outside promoters", "Gene body\n(non-promoter)")
            .replace("Shore (0–2 kb)", "Shore\n(0–2 kb)")
            .replace("Shelf (2–4 kb)", "Shelf\n(2–4 kb)")
            .replace("Open sea (>4 kb)", "Open sea\n(>4 kb)")
            for x in supplementary_names
        ]
        supplementary_fig = supplementary_plt.figure(figsize=(7.2, 5.5))
        supplementary_gs = supplementary_fig.add_gridspec(
            2, 2, height_ratios=[1.25, 1], hspace=0.7, wspace=0.35
        )
        supplementary_fig.subplots_adjust(left=0.12, right=0.98, bottom=0.11, top=0.88)
        supplementary_diagkey(supplementary_fig)
        for supplementary_k, (
            supplementary_state,
            supplementary_letter,
            supplementary_title,
        ) in enumerate([("5hmC_percent", "A", "5hmC"), ("5mC_percent", "B", "5mC")]):
            supplementary_ax = supplementary_fig.add_subplot(supplementary_gs[0, supplementary_k])
            for supplementary_i, supplementary_name in enumerate(supplementary_names):
                for supplementary_j, (supplementary_diag, supplementary_c) in enumerate(
                    supplementary_COL.items()
                ):
                    supplementary_vals = [
                        supplementary_num(x, supplementary_state)
                        for x in supplementary_r
                        if x["context_axis"] == supplementary_axis
                        and x["context"] == supplementary_name
                        and (x["diagnosis"] == supplementary_diag)
                    ]
                    supplementary_xp = supplementary_i + (supplementary_j - 0.5) * 0.25
                    supplementary_ax.scatter(
                        supplementary_xp
                        + supplementary_np.linspace(-0.045, 0.045, len(supplementary_vals)),
                        supplementary_vals,
                        c=supplementary_c,
                        s=13,
                        alpha=0.8,
                    )
                    supplementary_ax.plot(
                        [supplementary_xp - 0.065, supplementary_xp + 0.065],
                        [supplementary_np.mean(supplementary_vals)] * 2,
                        color="black",
                        lw=1.5,
                    )
            supplementary_ax.set_xticks(range(len(supplementary_names)), supplementary_display)
            supplementary_ax.set_ylabel("Fraction (%)")
            supplementary_ax.grid(axis="y", alpha=0.12)
            supplementary_panel(supplementary_ax, supplementary_letter, supplementary_title)
        supplementary_ax = supplementary_fig.add_subplot(supplementary_gs[1, :])
        supplementary_e = [
            next(
                (
                    x
                    for x in supplementary_eff
                    if x["context_axis"] == supplementary_axis and x["context"] == name
                )
            )
            for name in supplementary_names
        ]
        supplementary_v = supplementary_np.array(
            [supplementary_num(x, "GBM_minus_MEN_5hmC_pp") for x in supplementary_e]
        )
        supplementary_lo = supplementary_np.array(
            [supplementary_num(x, "bootstrap_95_low") for x in supplementary_e]
        )
        supplementary_hi = supplementary_np.array(
            [supplementary_num(x, "bootstrap_95_high") for x in supplementary_e]
        )
        supplementary_ax.errorbar(
            range(len(supplementary_e)),
            supplementary_v,
            yerr=[supplementary_v - supplementary_lo, supplementary_hi - supplementary_v],
            fmt="o",
            color="#455966",
            ms=4,
            capsize=3,
            lw=1,
        )
        supplementary_zero(supplementary_ax)
        supplementary_ax.set_xticks(range(len(supplementary_e)), supplementary_display)
        supplementary_ax.set_ylabel("GBM − MEN (pp)")
        supplementary_ax.set_xlim(-0.5, len(supplementary_e) - 0.5)
        supplementary_panel(supplementary_ax, "C", "5hmC contrast")
        supplementary_save(supplementary_fig, supplementary_n)
    supplementary_r = supplementary_tab("supplementary/results/molecular/Q5_repeat_sensitivity.tsv")
    (supplementary_fig, supplementary_axs) = supplementary_plt.subplots(1, 2, figsize=(7.2, 3.8))
    supplementary_fig.subplots_adjust(left=0.11, right=0.98, bottom=0.16, top=0.84, wspace=0.4)
    for supplementary_ax, supplementary_diag, supplementary_letter in zip(
        supplementary_axs, supplementary_COL, "AB"
    ):
        supplementary_rr = [
            x
            for x in supplementary_r
            if x["diagnosis"] == supplementary_diag and x["eligible_nonrepeat_comparison"] == "True"
        ]
        assert len(supplementary_rr) == (12 if supplementary_diag == "Glioblastoma" else 5)
        for supplementary_q in supplementary_rr:
            supplementary_ax.plot(
                [0, 1],
                [
                    supplementary_num(supplementary_q, "original_candidate_minus_comparison_pp"),
                    supplementary_num(supplementary_q, "nonrepeat_candidate_minus_comparison_pp"),
                ],
                "-o",
                c=supplementary_COL[supplementary_diag],
                lw=0.8,
                ms=3,
                alpha=0.75,
            )
        supplementary_zero(supplementary_ax)
        supplementary_ax.set_xticks([0, 1], ["Original blocks", "Outside repeats"])
        supplementary_ax.set_ylabel("Candidate − comparison (pp)")
        supplementary_panel(
            supplementary_ax,
            supplementary_letter,
            f"{supplementary_diag} (n = {len(supplementary_rr)})",
        )
    supplementary_save(supplementary_fig, 7)
    supplementary_r = supplementary_tab("supplementary/results/molecular/Q4a_alignment_depth.tsv")
    supplementary_assoc = supplementary_tab(
        "supplementary/results/molecular/Q4a_within_specimen_associations.tsv"
    )
    supplementary_fig = supplementary_plt.figure(figsize=(7.2, 6.8))
    supplementary_gs = supplementary_fig.add_gridspec(
        3, 2, height_ratios=[1, 1, 1], hspace=0.72, wspace=0.42
    )
    supplementary_fig.subplots_adjust(left=0.12, right=0.98, bottom=0.14, top=0.88)
    supplementary_diagkey(supplementary_fig)
    for supplementary_i, supplementary_gene in enumerate(["EGFR", "PDGFRA", "CDK4", "MDM2"]):
        supplementary_ax = supplementary_fig.add_subplot(
            supplementary_gs[supplementary_i // 2, supplementary_i % 2]
        )
        supplementary_rr = [
            q
            for q in supplementary_r
            if q["gene"] == supplementary_gene and q["kind"] == "amplification_screen"
        ]
        assert len(supplementary_rr) == 20
        for supplementary_q in supplementary_rr:
            supplementary_ax.scatter(
                supplementary_num(supplementary_q, "ratio_to_genome_mean"),
                supplementary_num(supplementary_q, "5hmC_percent"),
                s=22,
                c=supplementary_COL[supplementary_q["diagnosis"]],
                edgecolors=(
                    "black"
                    if supplementary_q["relative_depth_enrichment_screen"] == "True"
                    else "none"
                ),
                linewidths=0.7,
                alpha=0.9,
            )
        supplementary_ax.set_xscale("log", base=2)
        supplementary_ax.axvline(1, c="gray", ls="--", lw=0.6)
        supplementary_ax.axvline(2.5, c="gray", ls=":", lw=0.6)
        supplementary_ax.set_xlabel("Depth / genome mean")
        supplementary_ax.set_ylabel("5hmC (%)")
        supplementary_panel(supplementary_ax, chr(65 + supplementary_i), supplementary_gene)
    supplementary_ax = supplementary_fig.add_subplot(supplementary_gs[2, :])
    supplementary_q = {x["study_label"]: x for x in supplementary_assoc}
    supplementary_ax.scatter(
        range(20),
        [
            supplementary_num(supplementary_q[l], "depth_molecular_5hmC_spearman")
            for l in supplementary_labels
        ],
        c=[supplementary_COL[supplementary_q[l]["diagnosis"]] for l in supplementary_labels],
        s=21,
    )
    supplementary_zero(supplementary_ax)
    supplementary_patientaxis(supplementary_ax)
    supplementary_ax.set_ylim(-1, 1)
    supplementary_ax.set_ylabel("Spearman ρ")
    supplementary_panel(supplementary_ax, "E", "Depth and promoter co-occurrence")
    supplementary_save(supplementary_fig, 8)
    supplementary_r = [
        q
        for q in supplementary_tab("supplementary/results/molecular/Q4b_allele_loci.tsv")
        if q["eligible_allele_comparison"] == "True"
    ]
    assert len(supplementary_r) == 204
    (supplementary_fig, supplementary_ax) = supplementary_plt.subplots(figsize=(7.2, 3.8))
    supplementary_fig.subplots_adjust(left=0.11, right=0.98, bottom=0.22, top=0.85)
    supplementary_diagkey(supplementary_fig)
    for supplementary_i, supplementary_l in enumerate(supplementary_labels):
        supplementary_rr = [q for q in supplementary_r if q["study_label"] == supplementary_l]
        supplementary_ax.scatter(
            supplementary_i + supplementary_np.linspace(-0.23, 0.23, len(supplementary_rr)),
            [supplementary_num(q, "alt_minus_ref_5hmC_pp") for q in supplementary_rr],
            c=[supplementary_COL[q["diagnosis"]] for q in supplementary_rr],
            s=17,
            alpha=0.8,
        )
    supplementary_zero(supplementary_ax)
    supplementary_patientaxis(supplementary_ax)
    supplementary_ax.set_ylabel("Alternate − reference 5hmC (pp)")
    supplementary_save(supplementary_fig, 9)
    supplementary_caps = [
        (
            "Regional profiles of the largest exploratory 5hmC contrasts",
            "The 40 candidate 100-kb windows with the largest absolute glioblastoma-minus-meningioma 5hmC contrasts are shown in descending effect-size order. Each column is one specimen (15 glioblastomas, GBM; five meningiomas, MEN); each row is a genomic window. Colour denotes the specimen 5hmC fraction minus the equally weighted mean across all 20 specimens for that window, in percentage points (pp). Fractions use the same reference CpGs with at least five passing calls in every specimen. Windows met the absolute 5hmC contrast threshold of 2 pp and retained their direction after each specimen omission. Coordinates refer to hg38. Selection and display use the same cohort; this heatmap is descriptive and does not estimate diagnostic performance.",
        ),
        (
            "Abundance matching and sensitivity of frequency-normalised association",
            "(A,B) Patient-held-out candidate-minus-comparison scores before and after matching CpG-distance bin, two-site 5hmC frequencies and flanking GC. Small left points are unmatched; larger right points are matched, with lines joining the same patient. A shows excess co-occurrence on the 12-molecule scale in percentage points (pp); B shows dimensionless binary phi (φ). Only patients eligible for the matched comparison are displayed. A includes 15 glioblastomas and five meningiomas; B includes four and three, respectively. (C,D) Matched phi under standard filtering, mapping quality (MAPQ) ≥60, and CpG separation ≥25 bp, for glioblastoma and meningioma. Lines join eligible measurements within patients; missing points indicate insufficient support, not zero association. Red denotes glioblastoma and teal meningioma. Eligibility requires an overlapping pair weight of at least 100 and at least six contributing promoters in each panel. Observation subsets and eligible patient counts can differ between sensitivities; vertical scales differ between panels.",
        ),
        (
            "Promoter selection and support under diagnosis randomisation",
            "(A) Total selected candidate–comparison pairs across 20 patient-held-out folds for each of 499 alternative diagnosis allocations. (B) Number of distance intervals meeting patient-support requirements among the 36 contiguous unions of eight distance bins. Histograms include all alternative allocations; vertical teal lines indicate the observed allocation. Promoter selection and matching were repeated within each allocation. Patient scores required at least six eligible promoter pairs; an interval required at least three meningiomas and ten glioblastomas. Of the 499 alternative allocations, 250 had no eligible interval and were retained with maximum statistic zero. The randomisation test therefore includes selection and evaluability as components of the complete procedure.",
        ),
        (
            "Shared support for comparisons across modification encodings",
            "(A,B) Numbers of matched candidate–comparison promoter pairs contributing to each patient at 1–2,499 bp and 25–249 bp. A uses shared physical CpG pairs and reads across the 5hmC, 5mC and combined encodings. B imposes the stricter requirement of sufficiently polymorphic pairs with finite observed and constrained-null phi (φ) for all three encodings. Cell numbers and colour both indicate promoter-pair counts; patient summaries require at least six pairs. No patient met the requirement for the shared, frequency-normalised comparison at 25–249 bp (glioblastoma, 0/15; meningioma, 0/5). This full-pair phi comparison differs from the 12-molecule-resampled phi analysis in Fig. S2. Insufficient shared support prevents ranking modification specificity and does not demonstrate absence of an effect.",
        ),
        (
            "Modification fractions across broad gene annotations",
            "(A,B) Specimen-level 5hmC and 5mC fractions for promoters, gene bodies outside promoters, and intergenic regions. Each point represents one patient (red, 15 glioblastomas; teal, five meningiomas); black lines indicate group means. CpGs are assigned hierarchically to promoters, then remaining gene bodies, then intergenic regions. The same common CpGs are used for every specimen, with equal CpG weights within each category and specimen. (C) Glioblastoma-minus-meningioma mean 5hmC contrasts in percentage points (pp). Points are contrasts and error bars are pointwise 95% percentile intervals from 5,000 within-diagnosis specimen-bootstrap resamples, without multiple-testing adjustment. All three intervals include zero. Genomic annotation membership does not measure expression or regulatory activity.",
        ),
        (
            "Modification fractions across CpG-island contexts",
            "(A,B) Specimen-level 5hmC and 5mC fractions for CpG islands, shores within 2 kb, shelves 2–4 kb away, and open sea beyond 4 kb. Each point represents one patient (red, 15 glioblastomas; teal, five meningiomas); black lines indicate group means. Categories partition the common reference CpGs, and every specimen is summarised using the same CpGs with equal weights. (C) Glioblastoma-minus-meningioma mean 5hmC contrasts in percentage points (pp). Error bars are pointwise 95% percentile intervals from 5,000 within-diagnosis specimen-bootstrap resamples, without multiple-testing adjustment. All four intervals include zero. These categories describe genomic context rather than a tumour-specific mechanism.",
        ),
        (
            "Molecular co-occurrence after exclusion of annotated repeats",
            "(A,B) Original candidate-minus-comparison excess 5hmC co-occurrence and the corresponding score after restricting the original three-CpG blocks to those whose complete CpG dyads lie outside RepeatMasker intervals. Lines join measurements from the same patient: 12 eligible glioblastomas (A) and five meningiomas (B). Eligibility requires at least three remaining blocks per promoter and six matched promoter pairs per patient. Scores use the original 12-molecule resampling and independently permuted-column comparison, expressed in percentage points (pp); they are not the later molecule/CpG-count-preserving-null residuals. All five meningiomas retain a positive paired excess. Vertical scales differ. This sensitivity was specified after inspecting repeat overlap and does not establish repeat-family enrichment.",
        ),
        (
            "Targeted relative alignment depth and 5hmC",
            "(A–D) Relative alignment depth and 5hmC at 20-kb promoter-neighbourhood targets for EGFR, PDGFRA, CDK4 and MDM2, respectively. Each point represents one patient (red, 15 glioblastomas; teal, five meningiomas). The horizontal axis is depth divided by the specimen genome mean, on a base-2 logarithmic scale. Dashed and dotted reference lines mark ratios of 1 and 2.5. Black outlines identify screen-positive target/specimen intervals: depth at least 2.5-fold both the genome mean and the GC-fitted comparison expectation, with at least 80% of 1-kb bins at least twofold the genome mean. Fifteen target/specimen intervals pass, all in glioblastomas. (E) Within-specimen Spearman correlations between GC-normalised depth and original promoter 5hmC co-occurrence, using eligible candidate and comparison promoters. Alignment depth was calculated from aligned match blocks independently of modification states, using the same quality-filtered reads. These targeted relative-depth summaries do not establish absolute or ploidy-adjusted copy number, amplicon boundaries, ecDNA, or causation.",
        ),
        (
            "Local marker-associated 5hmC differences",
            "Each point is one of 204 eligible specimen/promoter comparisons across 20 patients (red, glioblastoma; teal, meningioma). One sequence-selected marker per specimen/promoter defines local reference- and alternate-allele read groups, evaluated at shared CpGs. The vertical axis shows alternate-minus-reference mean 5hmC in percentage points (pp); horizontal jitter separates points within a patient. CD200 in GBM-05 (+10.82 pp) and CLEC18A in GBM-09 (−15.39 pp) meet the descriptive criteria of an absolute contrast of at least 5 pp and at least 80% CpG-direction agreement; both are comparison promoters. Reference/alternate direction is marker-specific and does not identify parental or somatic origin. These are exploratory local observations, not multiple-testing-adjusted discoveries, independently validated variants or chromosome-wide haplotypes.",
        ),
    ]
    assert all((len(t.split()) <= 15 and len(c.split()) <= 300 for (t, c) in supplementary_caps))
    (supplementary_O / "figure_captions.json").write_text(
        supplementary_json.dumps(
            [
                {"figure": f"S{i}", "title": t, "legend": c}
                for (i, (t, c)) in enumerate(supplementary_caps, 1)
            ],
            indent=2,
        )
        + "\n"
    )
    (supplementary_O / "figure_text_inventory.json").write_text(
        supplementary_json.dumps(supplementary_figtexts, indent=2) + "\n"
    )
    supplementary_fontdir = (
        supplementary_Path(supplementary_matplotlib.get_data_path()) / "fonts/ttf"
    )
    supplementary_pdfmetrics.registerFont(
        supplementary_TTFont("DVS", str(supplementary_fontdir / "DejaVuSans.ttf"))
    )
    supplementary_pdfmetrics.registerFont(
        supplementary_TTFont("DVS-Bold", str(supplementary_fontdir / "DejaVuSans-Bold.ttf"))
    )
    supplementary_style = supplementary_ParagraphStyle(
        "legend", fontName="DVS", fontSize=9, leading=13, textColor="#222222"
    )
    supplementary_title_style = supplementary_ParagraphStyle(
        "title", fontName="DVS-Bold", fontSize=12, leading=16
    )
    supplementary_writer = supplementary_PdfWriter()
    (supplementary_W, supplementary_H) = supplementary_A4
    supplementary_margin = 42
    supplementary_width = supplementary_W - 2 * supplementary_margin
    supplementary_layout = []
    for supplementary_i, (supplementary_title, supplementary_legend) in enumerate(
        supplementary_caps, 1
    ):
        supplementary_buf = supplementary_io.BytesIO()
        supplementary_cv = supplementary_canvas.Canvas(supplementary_buf, pagesize=supplementary_A4)
        supplementary_cv.setFont("DVS", 8)
        supplementary_cv.setFillColorRGB(0.35, 0.4, 0.44)
        supplementary_cv.drawString(
            supplementary_margin,
            supplementary_H - 28,
            "ADDITIONAL FILE 1  |  SUPPLEMENTARY FIGURES",
        )
        supplementary_cv.setStrokeColorRGB(0.8, 0.83, 0.85)
        supplementary_cv.line(
            supplementary_margin,
            supplementary_H - 35,
            supplementary_W - supplementary_margin,
            supplementary_H - 35,
        )
        supplementary_p = supplementary_Paragraph(
            f"Figure S{supplementary_i}. {supplementary_html.escape(supplementary_title)}",
            supplementary_title_style,
        )
        (supplementary__, supplementary_th) = supplementary_p.wrap(supplementary_width, 100)
        supplementary_ty = supplementary_H - 51 - supplementary_th
        supplementary_p.drawOn(supplementary_cv, supplementary_margin, supplementary_ty)
        supplementary_lp = supplementary_Paragraph(
            supplementary_html.escape(supplementary_legend), supplementary_style
        )
        (supplementary__, supplementary_lh) = supplementary_lp.wrap(
            supplementary_width, supplementary_H
        )
        supplementary_fr = supplementary_PdfReader(
            str(supplementary_O / f"Figure_S{supplementary_i}.pdf")
        )
        supplementary_fp = supplementary_fr.pages[0]
        (supplementary_fw, supplementary_fh) = (
            float(supplementary_fp.mediabox.width),
            float(supplementary_fp.mediabox.height),
        )
        supplementary_maxheight = supplementary_ty - 14 - supplementary_lh - 18 - 48
        supplementary_scale = min(
            supplementary_width / supplementary_fw, supplementary_maxheight / supplementary_fh
        )
        (supplementary_pw, supplementary_ph) = (
            supplementary_fw * supplementary_scale,
            supplementary_fh * supplementary_scale,
        )
        supplementary_px = (supplementary_W - supplementary_pw) / 2
        supplementary_py = supplementary_ty - 12 - supplementary_ph
        supplementary_ly = supplementary_py - 15 - supplementary_lh
        assert supplementary_ly >= 45 and supplementary_scale > 0
        supplementary_lp.drawOn(supplementary_cv, supplementary_margin, supplementary_ly)
        supplementary_cv.setFillColorRGB(0.45, 0.45, 0.45)
        supplementary_cv.setFont("DVS", 8)
        supplementary_cv.drawRightString(
            supplementary_W - supplementary_margin, 25, f"{supplementary_i} / 9"
        )
        supplementary_cv.showPage()
        supplementary_cv.save()
        supplementary_page = supplementary_PdfReader(
            supplementary_io.BytesIO(supplementary_buf.getvalue())
        ).pages[0]
        supplementary_page.merge_transformed_page(
            supplementary_fp,
            supplementary_Transformation()
            .scale(supplementary_scale)
            .translate(supplementary_px, supplementary_py),
        )
        supplementary_writer.add_page(supplementary_page)
        supplementary_layout.append(
            {
                "figure": supplementary_i,
                "legend_words": len(supplementary_legend.split()),
                "plot_width_pt": supplementary_pw,
                "plot_height_pt": supplementary_ph,
                "legend_bottom_pt": supplementary_ly,
            }
        )
    supplementary_writer.add_metadata(
        {
            "/Title": "Additional file 1: supplementary figures S1–S9",
            "/Author": "Alexander Miller-Michlits; Yelyzaveta Miller-Michlits; Adelheid Woehrer",
            "/Subject": "Regional 5mC and 5hmC composition and single-molecule organisation in brain tumors",
        }
    )
    supplementary_out = (
        supplementary_R / "supplementary/5hmC_Additional_File_1_Supplementary_Figures.pdf"
    )
    with supplementary_out.open("wb") as supplementary_f:
        supplementary_writer.write(supplementary_f)
    assert supplementary_out.stat().st_size < 20000000
    supplementary_md = "# Additional file 1: supplementary figures S1–S9\n\nCombined submission file: [PDF](5hmC_Additional_File_1_Supplementary_Figures.pdf).\n"
    for supplementary_i, (supplementary_t, supplementary_c) in enumerate(supplementary_caps, 1):
        supplementary_md += f"\n## Figure S{supplementary_i}. {supplementary_t}\n\n![Figure S{supplementary_i}](figures/manuscript/Figure_S{supplementary_i}.png)\n\n{supplementary_c}\n"
    (supplementary_R / "supplementary/ADDITIONAL_FILE_1.md").write_text(supplementary_md)
    (supplementary_O / "layout_validation.json").write_text(
        supplementary_json.dumps(supplementary_layout, indent=2) + "\n"
    )
    print(
        supplementary_json.dumps(
            {
                "pages": len(supplementary_writer.pages),
                "pdf_bytes": supplementary_out.stat().st_size,
                "figures": 9,
                "legends": 9,
                "layout": supplementary_layout,
            },
            indent=2,
        )
    )
    _runtime.finish("manuscript_supplementary_figures")


def run_supplementary():
    """Execute the supplementary workflow stage."""
    initialize_supplementary()
