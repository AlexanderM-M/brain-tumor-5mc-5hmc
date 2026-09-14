"""Context analysis stages for the brain tumour 5mC/5hmC study.

Each stage has explicit initialization and execution; importing is read-only.
See docs/workflow.md for inputs, stage order and reproduction limits.
"""

from . import runtime as _runtime

from .runtime import config as _config

from . import molecules as _molecules
from . import regional as _regional


# DEPTH


from pathlib import Path as depth_Path
import json as depth_json, csv as depth_csv, numpy as depth_np, pysam as depth_pysam
from scipy.stats import spearmanr as depth_spearmanr
import matplotlib.pyplot as depth_plt


def depth_overlap_bp(blocks, start, end):
    return sum((max(0, min(end, b) - max(start, a)) for (a, b) in blocks))


def depth_run():
    samples = sorted(
        depth_json.loads((depth_D / "cohort.json").read_text())["samples"],
        key=lambda x: x["study_label"],
    )
    panel = depth_json.loads((depth_S / "plans/panel.json").read_text())
    fa = depth_pysam.FastaFile(_config.reference)
    gc = []
    for p in panel:
        seq = fa.fetch(p["chromosome"], p["start"], p["end"]).upper()
        gc.append((seq.count("G") + seq.count("C")) / len(seq))
    gc = depth_np.array(gc)
    controls = depth_np.array([p["kind"] == "comparison" for p in panel])
    rows = []
    binrows = []
    for s in samples:
        label = s["study_label"]
        aligned = depth_np.zeros(len(panel), depth_np.int64)
        bins = [depth_np.arange(p["start"], p["end"] + 1, 1000) for p in panel]
        binbp = [depth_np.zeros(len(v) - 1, depth_np.int64) for v in bins]
        seen = set()
        with depth_pysam.AlignmentFile(
            str(depth_S / "targeted_bams" / label / "reads.bam"), "rb"
        ) as b:
            for r in b:
                if not _molecules.blocks_eligible(r) or r.query_name in seen:
                    continue
                seen.add(r.query_name)
                blocks = r.get_blocks()
                for k, p in enumerate(panel):
                    if (
                        r.reference_name != p["chromosome"]
                        or r.reference_start >= p["end"]
                        or r.reference_end <= p["start"]
                    ):
                        continue
                    aligned[k] += depth_overlap_bp(blocks, p["start"], p["end"])
                    for z, (a, c) in enumerate(zip(bins[k][:-1], bins[k][1:])):
                        binbp[k][z] += depth_overlap_bp(blocks, int(a), int(c))
        depths = aligned / depth_np.array([p["end"] - p["start"] for p in panel])
        global_depth = depth_json.loads((depth_D / "results" / label / "qc.json").read_text())[
            "aligned_genome_equivalent"
        ]
        control_median = depth_np.median(depths[controls])
        y = depth_np.log2(depth_np.maximum(depths[controls], 0.01))
        x = gc[controls]
        keep = depth_np.ones(len(x), bool)
        for _ in range(3):
            (slope, intercept) = depth_np.polyfit(x[keep], y[keep], 1)
            res = y - (slope * x + intercept)
            mad = depth_np.median(abs(res - depth_np.median(res)))
            new = abs(res - depth_np.median(res)) <= max(0.5, 3 * 1.4826 * mad)
            if new.sum() < 8:
                break
            keep = new
        expected = 2 ** (slope * gc + intercept)
        for k, p in enumerate(panel):
            a = depth_np.load(depth_S / "results" / label / (p["locus_id"] + ".npz"))
            mat = a["states"][:, a["common_mask"]]
            c = depth_np.column_stack([(mat == v).sum(0) for v in [0, 1, 2]])
            den = c.sum(1)
            assert depth_np.all(den >= 5)
            (m, h) = (100 * c[:, 1:3] / den[:, None]).mean(0)
            raw = depths[k] / global_depth
            corr = depths[k] / expected[k]
            br = binbp[k] / 1000 / global_depth
            enriched = raw >= 2.5 and corr >= 2.5 and (depth_np.mean(br >= 2) >= 0.8)
            rows.append(
                {
                    "study_label": label,
                    "diagnosis": s["diagnosis"],
                    "locus_id": p["locus_id"],
                    "gene": p["gene"],
                    "kind": p["kind"],
                    "GC_fraction": float(gc[k]),
                    "aligned_depth": float(depths[k]),
                    "genome_mean_aligned_depth": float(global_depth),
                    "ratio_to_genome_mean": float(raw),
                    "ratio_to_comparison_median": float(depths[k] / control_median),
                    "ratio_to_GC_fitted_comparisons": float(corr),
                    "fraction_1kb_bins_at_least_2fold_genome_mean": float(depth_np.mean(br >= 2)),
                    "relative_depth_enrichment_screen": bool(enriched),
                    "common_CpGs": len(den),
                    "5mC_percent": float(m),
                    "5hmC_percent": float(h),
                    "combined_percent": float(m + h),
                }
            )
            for z, (lo, hi) in enumerate(zip(bins[k][:-1], bins[k][1:])):
                binrows.append(
                    {
                        "study_label": label,
                        "locus_id": p["locus_id"],
                        "chromosome": p["chromosome"],
                        "start": int(lo),
                        "end": int(hi),
                        "depth": float(binbp[k][z] / 1000),
                        "ratio_to_genome_mean": float(br[z]),
                    }
                )
    _regional.regional_table(depth_S / "results/Q4a_alignment_depth.tsv", rows)
    _regional.regional_table(depth_S / "results/Q4a_depth_bins.tsv", binrows)
    genes = ["EGFR", "PDGFRA", "CDK4", "MDM2"]
    colors = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    (fig, axes) = depth_plt.subplots(1, 4, figsize=(15, 4.5))
    for ax, gene in zip(axes, genes):
        for r in rows:
            if r["kind"] != "amplification_screen" or r["gene"] != gene:
                continue
            ax.scatter(
                r["ratio_to_genome_mean"],
                r["5hmC_percent"],
                c=colors[r["diagnosis"]],
                edgecolors="black" if r["relative_depth_enrichment_screen"] else "none",
                s=40,
            )
        ax.set_xscale("log", base=2)
        ax.axvline(1, c="gray", ls="--", lw=0.7)
        ax.axvline(2.5, c="gray", ls=":", lw=0.7)
        ax.set_xlabel("Depth / specimen genome mean")
        ax.set_ylabel("5hmC at common CpGs (%)")
        ax.set_title(gene)
        ax.grid(alpha=0.15)
    fig.suptitle(
        "Targeted copy-abundance screen and 5hmC\nRed: glioblastoma; teal: meningioma. Outlined points pass the relative-depth screen.",
        fontsize=11,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "Q4a_01_depth_and_5hmC")
    candidates = [
        r
        for r in rows
        if r["kind"] == "amplification_screen" and r["relative_depth_enrichment_screen"]
    ]
    summary = []
    for gene in genes:
        rs = [
            r
            for r in rows
            if r["gene"] == gene
            and r["kind"] == "amplification_screen"
            and (r["diagnosis"] == "Glioblastoma")
        ]
        rho = depth_spearmanr(
            [r["ratio_to_genome_mean"] for r in rs], [r["5hmC_percent"] for r in rs]
        ).statistic
        summary.append(
            {
                "gene": gene,
                "GBM_cases_passing_screen": sum(
                    (r["relative_depth_enrichment_screen"] for r in rs)
                ),
                "GBM_depth_5hmC_spearman_descriptive": float(rho),
            }
        )
    q3 = list(
        depth_csv.DictReader((depth_S / "results/Q3_locus_metrics.tsv").open(), delimiter="\t")
    )
    q3lookup = {(r["study_label"], r["locus_id"]): r for r in q3}
    associations = []
    for s in samples:
        rs = [
            r
            for r in rows
            if r["study_label"] == s["study_label"]
            and r["kind"] != "amplification_screen"
            and (q3lookup[r["study_label"], r["locus_id"]]["eligible_locus"] == "True")
        ]
        xx = [depth_np.log2(r["ratio_to_GC_fitted_comparisons"]) for r in rs]
        yy = [
            float(q3lookup[r["study_label"], r["locus_id"]]["mean_joint_5hmC_excess_pp"])
            for r in rs
        ]
        rho = float(depth_spearmanr(xx, yy).statistic) if len(rs) >= 10 else None
        associations.append(
            {
                "study_label": s["study_label"],
                "diagnosis": s["diagnosis"],
                "loci": len(rs),
                "depth_molecular_5hmC_spearman": rho,
            }
        )
    _regional.regional_table(depth_S / "results/Q4a_within_specimen_associations.tsv", associations)
    (fig, ax) = depth_plt.subplots(figsize=(10, 4))
    ax.scatter(
        range(20),
        [r["depth_molecular_5hmC_spearman"] for r in associations],
        c=[colors[r["diagnosis"]] for r in associations],
    )
    ax.axhline(0, c="gray", lw=0.7)
    ax.set_ylim(-1, 1)
    ax.set_xticks(range(20), [r["study_label"] for r in associations], rotation=60)
    ax.set_ylabel("Within-specimen Spearman correlation")
    ax.set_title(
        "Alignment depth versus molecular 5hmC co-occurrence across Q3 promoters\nDescriptive association; no absolute copy-number or causal inference"
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "Q4a_02_depth_molecular_association")
    _regional.regional_save(
        depth_S / "results/Q4a_summary.json",
        {
            "amplification_target_screen": summary,
            "screen_positive_target_specimen_pairs": [
                {
                    k: r[k]
                    for k in [
                        "study_label",
                        "gene",
                        "ratio_to_genome_mean",
                        "ratio_to_GC_fitted_comparisons",
                        "5hmC_percent",
                    ]
                }
                for r in candidates
            ],
            "within_specimen_associations": associations,
        },
    )
    text = [
        "# Q4a: are molecular patterns associated with copy abundance?",
        "",
        "## Approach",
        "",
        "Alignment depth was calculated directly from aligned match blocks in the targeted BAMs, independently of modification states. The same read-level QC and exact read-ID deduplication were used as in the extraction. Mean depth in each target was normalized to the full specimen genome-wide aligned depth and separately to comparison-promoter coverage. A robust linear GC fit on the 12 comparison promoters provides a sensitivity normalization. This is a targeted relative-depth screen, not an absolute or ploidy-aware copy-number call.",
        "",
        "An oncogene target passes the exploratory screen only when its depth is >=2.5-fold both the genome mean and the GC-fitted comparison expectation, and at least 80% of its 1kb bins have >=2-fold genome-mean depth. These thresholds were set before this screening analysis. The 20kb targets sample the promoter neighbourhood; they do not define complete amplicon boundaries.",
        "",
        "## Results",
        "",
        depth_json.dumps(summary, indent=2),
        "",
        "Screen-positive target/specimen pairs: "
        + str(len(candidates))
        + ". Complete values and all within-specimen associations are in Q4a_summary.json and Q4a_alignment_depth.tsv.",
        "",
        "## Answer and limits",
        "",
        "The figures and tables test whether local DNA abundance covaries with 5hmC or with the Q3 molecular co-occurrence score. Correlations across 15 GBMs are descriptive and may be dominated by individual specimens. Twelve small comparison regions cannot establish absolute tumour ploidy or control every GC and mapping effect. Positive depth signals are candidates for formal genome-wide CNV confirmation. A negative screen does not exclude an alteration outside the sampled interval. Amplification does not establish ecDNA, and these data do not establish that copy-number changes cause the molecular patterns. Q4b tests a different explanation: local allele-associated modification.",
    ]
    (depth_S / "results/Q4a_ANSWER.md").write_text("\n".join(text) + "\n")
    _regional.regional_save(
        depth_S / "results/Q4a_complete.json",
        {
            "analysis_and_figures_generated": True,
            "scope": "Targeted copy-abundance association and feasibility; not genome-wide definitive CNV",
        },
    )
    print(
        "Q4a COMPLETE", depth_json.dumps(summary), "screen positives", len(candidates), flush=True
    )


def initialize_depth():
    """Initialize the depth stage once; load its declared inputs."""
    global depth_D, depth_R, depth_S
    if _runtime.initialized("story_copy_number"):
        return
    _runtime.begin("story_copy_number")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    depth_R = depth_Path(_config.workspace)
    depth_D = depth_R / ".analysis"
    depth_S = depth_D / "story_analysis"
    _runtime.finish("story_copy_number")


def run_depth():
    """Execute the depth workflow stage."""
    initialize_depth()
    assert depth_overlap_bp([(0, 10), (20, 30)], 5, 25) == 10
    depth_run()


# ALLELES


from pathlib import Path as alleles_Path
import json as alleles_json, csv as alleles_csv, hashlib as alleles_hashlib
import numpy as alleles_np, pysam as alleles_pysam
import matplotlib.pyplot as alleles_plt


def alleles_marker_allowed(seq, i):
    if seq[i] not in alleles_BASES:
        return False
    if seq[i] == "C" and seq[i + 1] == "G" or (seq[i] == "G" and seq[i - 1] == "C"):
        return False
    near = seq[i - 5 : i + 6]
    return not any((b * 4 in near for b in alleles_BASES))


def alleles_fraction(x, state):
    den = (x >= 0).sum(0)
    return (
        alleles_np.divide(
            (x == state).sum(0), den, out=alleles_np.full(x.shape[1], alleles_np.nan), where=den > 0
        )
        * 100
    )


def alleles_sample(s):
    label = s["study_label"]
    panel = [
        p
        for p in alleles_json.loads((alleles_S / "plans/panel.json").read_text())
        if p["kind"] != "amplification_screen"
    ]
    rows = []
    markers = []
    checks = []
    depths = {
        (r["study_label"], r["locus_id"]): r
        for r in alleles_csv.DictReader(
            (alleles_S / "results/Q4a_alignment_depth.tsv").open(), delimiter="\t"
        )
    }
    with (
        alleles_pysam.AlignmentFile(
            str(alleles_S / "targeted_bams" / label / "reads.bam"), "rb"
        ) as bam,
        alleles_pysam.FastaFile(_config.reference) as fa,
    ):
        for p in panel:
            key = p["locus_id"]
            seq = fa.fetch(p["chromosome"], p["start"] - 6, p["end"] + 6).upper()
            allowed = alleles_np.array(
                [alleles_marker_allowed(seq, i + 6) for i in range(p["end"] - p["start"])]
            )
            valid = {}
            candidates = []
            for col in bam.pileup(
                p["chromosome"],
                p["start"],
                p["end"],
                truncate=True,
                stepper="nofilter",
                min_base_quality=0,
                compute_baq=False,
                ignore_overlaps=False,
                ignore_orphans=False,
                max_depth=100000,
            ):
                pos = col.reference_pos
                if not allowed[pos - p["start"]]:
                    continue
                counts = alleles_np.zeros((4, 2), int)
                obs = {}
                for pr in col.pileups:
                    r = pr.alignment
                    rid = r.query_name
                    if rid not in valid:
                        valid[rid] = _molecules.blocks_eligible(r) and 30 <= r.mapping_quality < 255
                    if not valid[rid] or pr.is_del or pr.is_refskip or (pr.query_position is None):
                        continue
                    q = pr.query_position
                    if (
                        q < 100
                        or q >= r.query_length - 100
                        or r.query_qualities is None
                        or (r.query_qualities[q] < 25)
                    ):
                        continue
                    base = r.query_sequence[q]
                    if base not in alleles_BASES or rid in obs:
                        continue
                    obs[rid] = (base, bool(r.is_reverse))
                    counts[alleles_BASES.index(base), int(r.is_reverse)] += 1
                total = int(counts.sum())
                if total < 16:
                    continue
                ref = seq[pos - p["start"] + 6]
                ri = alleles_BASES.index(ref)
                ac = counts.sum(1).copy()
                ac[ri] = -1
                ai = int(ac.argmax())
                alt = alleles_BASES[ai]
                nr = int(counts[ri].sum())
                na = int(counts[ai].sum())
                if (
                    min(nr, na) < 8
                    or min(counts[ri].min(), counts[ai].min()) < 2
                    or (nr + na) / total < 0.95
                    or (not 0.25 <= na / (nr + na) <= 0.75)
                ):
                    continue
                rec = {
                    "study_label": label,
                    "locus_id": key,
                    "gene": p["gene"],
                    "chromosome": p["chromosome"],
                    "position_1based": pos + 1,
                    "reference": ref,
                    "alternate": alt,
                    "ref_reads": nr,
                    "alt_reads": na,
                    "ref_forward": int(counts[ri, 0]),
                    "ref_reverse": int(counts[ri, 1]),
                    "alt_forward": int(counts[ai, 0]),
                    "alt_reverse": int(counts[ai, 1]),
                    "qualified_base_depth": total,
                }
                markers.append(rec)
                candidates.append((rec, obs))
            row = {
                "study_label": label,
                "diagnosis": s["diagnosis"],
                "locus_id": key,
                "gene": p["gene"],
                "kind": p["kind"],
                "candidate_markers": len(candidates),
                "marker_found": bool(candidates),
                "eligible_allele_comparison": False,
                "reason": "no marker meeting filters",
                "chromosome": p["chromosome"],
                "position_1based": "",
                "reference": "",
                "alternate": "",
                "ref_reads": 0,
                "alt_reads": 0,
                "shared_CpGs": 0,
                "ref_5hmC_percent": "",
                "alt_5hmC_percent": "",
                "alt_minus_ref_5hmC_pp": "",
                "alt_minus_ref_5mC_pp": "",
                "CpG_direction_agreement": "",
                "bootstrap_95_low": "",
                "bootstrap_95_high": "",
                "bootstrap_valid": 0,
                "descriptive_large_concordant_difference": False,
                "depth_ratio_to_genome": float(depths[label, key]["ratio_to_genome_mean"]),
            }
            if candidates:
                (rec, obs) = sorted(
                    candidates,
                    key=lambda z: (
                        -min(z[0]["ref_reads"], z[0]["alt_reads"]),
                        -z[0]["qualified_base_depth"],
                        z[0]["position_1based"],
                    ),
                )[0]
                row.update(
                    {
                        k: rec[k]
                        for k in [
                            "position_1based",
                            "reference",
                            "alternate",
                            "ref_reads",
                            "alt_reads",
                        ]
                    }
                )
                pos = rec["position_1based"] - 1
                verified = {}
                for r in bam.fetch(p["chromosome"], pos, pos + 1):
                    if r.query_name not in obs or r.query_name in verified:
                        continue
                    qq = next(
                        (q for (q, rp) in r.get_aligned_pairs(matches_only=True) if rp == pos), None
                    )
                    assert qq is not None and r.query_sequence[qq] == obs[r.query_name][0]
                    verified[r.query_name] = True
                assert len(verified) == len(obs)
                checks.append(
                    {
                        "locus_id": key,
                        "marker_1based": pos + 1,
                        "independent_read_base_checks": len(verified),
                    }
                )
                z = alleles_np.load(alleles_S / "results" / label / (key + ".npz"))
                ids = z["read_ids"]
                states = z["states"]
                cpg = z["cpg_positions"]
                common = z["common_mask"] & (abs(cpg.astype(alleles_np.int64) - pos) > 5)
                groups = alleles_np.array(
                    [
                        (
                            0
                            if str(rid) in obs and obs[str(rid)][0] == rec["reference"]
                            else (
                                1
                                if str(rid) in obs and obs[str(rid)][0] == rec["alternate"]
                                else -1
                            )
                        )
                        for rid in ids
                    ]
                )
                x = states[groups == 0]
                y = states[groups == 1]
                shared = common & ((x >= 0).sum(0) >= 4) & ((y >= 0).sum(0) >= 4)
                row["shared_CpGs"] = int(shared.sum())
                row["reason"] = "insufficient shared CpGs"
                alleles_np.savez_compressed(
                    alleles_S / "results" / label / (key + "_allele_groups.npz"),
                    groups=groups,
                    shared_mask=shared,
                    marker_position_0based=pos,
                )
                if shared.sum() >= 5 and min(len(x), len(y)) >= 8:
                    x = x[:, shared]
                    y = y[:, shared]
                    xh = alleles_fraction(x, 2)
                    yh = alleles_fraction(y, 2)
                    effect = float(alleles_np.mean(yh - xh))
                    direction = float(alleles_np.mean((yh - xh) * alleles_np.sign(effect) > 0))
                    seed = int.from_bytes(
                        alleles_hashlib.sha256((label + key + "alleles").encode()).digest()[:8],
                        "little",
                    )
                    rng = alleles_np.random.default_rng(seed)
                    boot = []
                    for _ in range(500):
                        a = alleles_fraction(x[rng.integers(0, len(x), len(x))], 2)
                        b = alleles_fraction(y[rng.integers(0, len(y), len(y))], 2)
                        if alleles_np.isfinite(a).all() and alleles_np.isfinite(b).all():
                            boot.append(float(alleles_np.mean(b - a)))
                    ci = (
                        alleles_np.percentile(boot, [2.5, 97.5])
                        if len(boot) >= 400
                        else [alleles_np.nan, alleles_np.nan]
                    )
                    row.update(
                        eligible_allele_comparison=True,
                        reason="eligible",
                        ref_5hmC_percent=float(xh.mean()),
                        alt_5hmC_percent=float(yh.mean()),
                        alt_minus_ref_5hmC_pp=effect,
                        alt_minus_ref_5mC_pp=float(
                            alleles_np.mean(alleles_fraction(y, 1) - alleles_fraction(x, 1))
                        ),
                        CpG_direction_agreement=direction,
                        bootstrap_95_low=float(ci[0]),
                        bootstrap_95_high=float(ci[1]),
                        bootstrap_valid=len(boot),
                        descriptive_large_concordant_difference=abs(effect) >= 5
                        and direction >= 0.8,
                    )
            rows.append(row)
    (alleles_S / "results" / label / "allele_base_validation.json").write_text(
        alleles_json.dumps(checks, indent=2) + "\n"
    )
    _regional.regional_table(alleles_S / "results" / label / "allele_loci.tsv", rows)
    print(
        label,
        "ALLELES",
        sum((r["marker_found"] for r in rows)),
        sum((r["eligible_allele_comparison"] for r in rows)),
        flush=True,
    )
    return (rows, markers)


def alleles_summarize(samples, rows, markers):
    _regional.regional_table(alleles_S / "results/Q4b_allele_loci.tsv", rows)
    if markers:
        _regional.regional_table(alleles_S / "results/Q4b_candidate_markers.tsv", markers)
    panel = [
        p
        for p in alleles_json.loads((alleles_S / "plans/panel.json").read_text())
        if p["kind"] != "amplification_screen"
    ]
    labels = [s["study_label"] for s in samples]
    lookup = {(r["study_label"], r["locus_id"]): r for r in rows}
    a = alleles_np.array(
        [
            [
                (
                    2
                    if lookup[l, p["locus_id"]]["eligible_allele_comparison"]
                    else 1 if lookup[l, p["locus_id"]]["marker_found"] else 0
                )
                for l in labels
            ]
            for p in panel
        ]
    )
    from matplotlib.colors import ListedColormap

    (fig, ax) = alleles_plt.subplots(figsize=(12, 7))
    im = ax.imshow(
        a, vmin=0, vmax=2, cmap=ListedColormap(["#e7e7e7", "#dec778", "#278b94"]), aspect="auto"
    )
    ax.set_xticks(range(20), labels, rotation=60)
    ax.set_yticks(range(24), [p["locus_id"] + " " + p["gene"] for p in panel], fontsize=8)
    cb = fig.colorbar(im, ax=ax, ticks=[0, 1, 2])
    cb.ax.set_yticklabels(
        ["No eligible marker", "Marker; CpGs insufficient", "Allele comparison eligible"]
    )
    ax.set_title("Local allele-group feasibility at the fixed promoter panel")
    fig.tight_layout()
    _molecules.blocks_figure(fig, "Q4b_01_allele_feasibility")
    valid = [r for r in rows if r["eligible_allele_comparison"]]
    (fig, ax) = alleles_plt.subplots(figsize=(11, 4.5))
    colors = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    for i, l in enumerate(labels):
        rr = [r for r in valid if r["study_label"] == l]
        jitter = alleles_np.linspace(-0.23, 0.23, len(rr))
        ax.scatter(
            i + jitter,
            [r["alt_minus_ref_5hmC_pp"] for r in rr],
            c=[colors[r["diagnosis"]] for r in rr],
            s=25,
            alpha=0.8,
        )
    ax.axhline(0, c="gray", lw=0.8)
    ax.set_xticks(range(20), labels, rotation=60)
    ax.set_ylabel("Alternate minus reference 5hmC (percentage points)")
    ax.set_title(
        "Local marker-associated differences\nEach point is one eligible promoter; allele direction is arbitrary across markers"
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "Q4b_02_allele_differences")
    selected = sorted(
        valid, key=lambda r: (-abs(r["alt_minus_ref_5hmC_pp"]), r["study_label"], r["locus_id"])
    )[:4]
    if selected:
        from matplotlib.colors import BoundaryNorm

        (fig, axes) = alleles_plt.subplots(
            2, len(selected), figsize=(4 * len(selected), 6), squeeze=False
        )
        cmap = ListedColormap(["#eeeeee", "#cfcfcf", "#497ba6", "#df9d29"])
        norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], 4)
        for k, r in enumerate(selected):
            z = alleles_np.load(alleles_S / "results" / r["study_label"] / (r["locus_id"] + ".npz"))
            g = alleles_np.load(
                alleles_S / "results" / r["study_label"] / (r["locus_id"] + "_allele_groups.npz")
            )
            for group in [0, 1]:
                x = z["states"][g["groups"] == group][:, g["shared_mask"]]
                axes[group, k].imshow(
                    x, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm
                )
                axes[group, k].set_title(
                    r["study_label"]
                    + " "
                    + r["gene"]
                    + "\n"
                    + (
                        "Reference " + r["reference"]
                        if group == 0
                        else "Alternate " + r["alternate"]
                    ),
                    fontsize=9,
                )
                axes[group, k].set_xlabel("Shared CpGs")
                axes[group, k].set_ylabel("Molecules")
        fig.suptitle(
            "Largest observed local allele-group 5hmC differences (selected illustrations)\nGrey C; blue 5mC; gold 5hmC; pale missing. Selection exaggerates apparent effects.",
            fontsize=11,
        )
        fig.tight_layout()
        _molecules.blocks_figure(fig, "Q4b_03_allele_read_examples")
    summary = {
        "panel_specimen_combinations": len(rows),
        "combinations_with_candidate_marker": sum((r["marker_found"] for r in rows)),
        "eligible_allele_comparisons": len(valid),
        "specimens_with_eligible_comparison": len(set((r["study_label"] for r in valid))),
        "all_candidate_markers": len(markers),
        "descriptive_large_concordant_differences": sum(
            (r["descriptive_large_concordant_difference"] for r in valid)
        ),
        "median_absolute_5hmC_difference_pp": (
            float(alleles_np.median([abs(r["alt_minus_ref_5hmC_pp"]) for r in valid]))
            if valid
            else None
        ),
        "largest_examples": selected,
        "all_selected_marker_read_bases_independently_checked": True,
    }
    _regional.regional_save(alleles_S / "results/Q4b_summary.json", summary)
    answer = [
        "# Q4b: local marker-supported allele groups",
        "",
        "## Answer",
        "",
        f"{len(valid)} of {len(rows)} specimen/promoter combinations support an allele-group comparison, spanning {summary['specimens_with_eligible_comparison']} specimens. {summary['combinations_with_candidate_marker']} combinations have a candidate marker passing the sequence, balance and strand filters. {summary['descriptive_large_concordant_differences']} eligible comparisons show an absolute 5hmC difference of at least 5 percentage points with at least 80% of shared CpGs agreeing in direction. These are descriptive shortlist criteria, not corrected significance or validated allele-specific methylation calls.",
        "",
        "The table and read figures reveal whether molecular differences can coincide with local sequence alleles. They cannot establish that every Q3 pattern is explained by alleles. Read-depth limitations and marker absence leave individual loci unresolved. Reference/alternate direction has no common biological meaning across different variants.",
        "",
        "## Methods and validation",
        "",
        "One marker per specimen/promoter was chosen exclusively by allele read support, before examining modification differences. Candidate SNVs require MAPQ >=30, base quality >=25, 100bp read-end trimming, at least eight reads and at least two reads on each strand for each allele, allele fraction 25–75%, and >=95% reference-plus-alternate observations. Reference CpG-disrupting positions and local homopolymers were excluded. The exact marker bases on all assigned reads were independently checked with aligned pairs. Both alleles use the same reference CpGs, with at least four calls per allele at each of at least five CpGs, excluding positions within five bases of the marker.",
        "",
        "5hmC fractions are averaged equally across shared CpGs. The table includes 500 within-allele molecule bootstrap resamples; intervals are reported only if >=400 resamples retain observations at all fixed CpGs. These intervals describe read-sampling uncertainty in a specimen, not biological replication, variant-calling accuracy or multiple-testing-adjusted evidence.",
        "",
        "## Limits",
        "",
        "The SNVs are candidates inferred from these nanopore reads; no matched normal or orthogonal validation is available. Local biallelic groups are not necessarily exactly two physical chromosome copies in an aneuploid tumour. Mapping bias, sequence-context errors, tumour purity and additional alleles remain possible explanations. This analysis does not phase the genome, assign parental origin, distinguish somatic from germline variants, or identify tumour clones. Large examples were selected after observing the effects and need independent confirmation.",
        "",
        "## Outputs",
        "",
        "Q4b_allele_loci.tsv; Q4b_candidate_markers.tsv; Q4b_summary.json; figures Q4b_01–03.",
    ]
    (alleles_S / "results/Q4b_ANSWER.md").write_text("\n".join(answer) + "\n")
    _regional.regional_save(
        alleles_S / "results/Q4b_complete.json",
        {"analysis_and_figures_generated": True, "requires_output_review": True},
    )
    print(
        "Q4b COMPLETE",
        alleles_json.dumps({k: v for (k, v) in summary.items() if k != "largest_examples"}),
        flush=True,
    )


def initialize_alleles():
    """Initialize the alleles stage once; load its declared inputs."""
    global alleles_BASES, alleles_D, alleles_PLAN, alleles_R, alleles_S
    if _runtime.initialized("story_alleles"):
        return
    _runtime.begin("story_alleles")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    alleles_R = alleles_Path(_config.workspace)
    alleles_D = alleles_R / ".analysis"
    alleles_S = alleles_D / "story_analysis"
    alleles_PLAN = {
        "marker_MAPQ_min": 30,
        "marker_base_quality_min": 25,
        "read_end_trim": 100,
        "alternate_fraction": [0.25, 0.75],
        "minimum_reads_per_allele": 8,
        "minimum_reads_per_strand_per_allele": 2,
        "minimum_ref_plus_alt_fraction": 0.95,
        "exclude_reference_CpG_disrupting_marker": True,
        "exclude_homopolymer_run_at_least_4_within_5bp": True,
        "marker_selection": "maximum smaller allele count, then total count, then genomic coordinate; independent of modification effect",
        "minimum_common_CpGs": 5,
        "minimum_calls_per_allele_per_CpG": 4,
        "exclude_CpGs_within_marker_bp": 5,
        "bootstrap_read_resamples": 500,
        "descriptive_large_difference_pp": 5,
        "minimum_direction_concordance": 0.8,
        "scope": "local read groups at candidate SNVs; no validated variant, parental origin, somatic origin, tumour clone or whole-chromosome haplotype claim",
    }
    alleles_BASES = "ACGT"
    _runtime.finish("story_alleles")


def run_alleles():
    """Execute the alleles workflow stage."""
    initialize_alleles()
    global alleles_out, alleles_planpath, alleles_pool, alleles_samples
    alleles_planpath = alleles_S / "plans/Q4b_plan.json"
    if alleles_planpath.exists():
        assert alleles_json.loads(alleles_planpath.read_text()) == alleles_PLAN
    else:
        _regional.regional_save(alleles_planpath, alleles_PLAN)
    assert not alleles_marker_allowed("AAAAAACGAAAAAA", 6)
    alleles_samples = sorted(
        alleles_json.loads((alleles_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    with _runtime.process_pool(max_workers=_config.workers) as alleles_pool:
        alleles_out = list(alleles_pool.map(alleles_sample, alleles_samples))
    alleles_summarize(
        alleles_samples,
        [r for (rows, markers) in alleles_out for r in rows],
        [m for (rows, markers) in alleles_out for m in markers],
    )


# ANNOTATION


from pathlib import Path as annotation_Path
import json as annotation_json, csv as annotation_csv, gzip as annotation_gzip, time as annotation_time, urllib.request, urllib as annotation_urllib
import numpy as annotation_np
import matplotlib.pyplot as annotation_plt


def annotation_download(track, chrom, start=None, end=None, key=None):
    path = annotation_ANNOT / f"{track}_{key or chrom}.json.gz"
    if path.exists():
        with annotation_gzip.open(path, "rt") as f:
            d = annotation_json.load(f)
    else:
        if not _config.allow_annotation_download:
            raise FileNotFoundError(
                f"Archived annotation required: {path}. Configure allow_annotation_download explicitly for a new annotation snapshot."
            )
        u = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track={track};chrom={chrom};maxItemsOutput=1000000"
        if start is not None:
            u += f";start={start};end={end}"
        for attempt in range(3):
            try:
                with annotation_urllib.request.urlopen(u, timeout=45) as f:
                    d = annotation_json.load(f)
                break
            except Exception:
                if attempt == 2:
                    raise
                annotation_time.sleep(2)
        assert (
            track in d
            and (not d.get("maxItemsLimit"))
            and (len(d[track]) == d.get("itemsReturned", 0))
        )
        d["_retrieval_url"] = u
        with annotation_gzip.open(path, "wt") as f:
            annotation_json.dump(d, f)
        annotation_time.sleep(0.25)
    return d[track]


def annotation_union_bp(intervals, start, end):
    clipped = sorted(
        ((max(a, start), min(b, end)) for (a, b) in intervals if a < end and b > start)
    )
    total = 0
    stop = start
    for a, b in clipped:
        total += max(0, b - max(a, stop))
        stop = max(stop, b)
    return total


def annotation_annotations():
    pos = annotation_np.load(annotation_D / "reference/cpg_positions.npy", mmap_mode="r")
    off = annotation_np.load(annotation_D / "reference/cpg_offsets.npy")
    gc = annotation_np.full(len(pos), 2, annotation_np.int8)
    cc = annotation_np.full(len(pos), 3, annotation_np.int8)
    islands = {}
    for i in range(22):
        chrom = f"chr{i + 1}"
        rr = annotation_download("cpgIslandExt", chrom)
        islands[chrom] = [(int(r["chromStart"]), int(r["chromEnd"])) for r in rr]
        (lo, hi) = map(int, off[i : i + 2])
        p = pos[lo:hi]
        for padding, value in [(4000, 2), (2000, 1), (0, 0)]:
            for a, b in islands[chrom]:
                (x, y) = annotation_np.searchsorted(p, [max(0, a - padding), b + padding])
                cc[lo + x : lo + y] = value
        print(chrom, "CpG context", len(rr), flush=True)
    with annotation_gzip.open(
        annotation_D / "pooled_analysis/tables/annotation_features.tsv.gz", "rt"
    ) as f:
        features = list(annotation_csv.DictReader(f, delimiter="\t"))
    for ft, value in [("gene_body", 1), ("promoter", 0)]:
        for r in features:
            if r["feature"] != ft:
                continue
            ch = int(r["chromosome"][3:]) - 1
            (lo, hi) = map(int, off[ch : ch + 2])
            (x, y) = annotation_np.searchsorted(pos[lo:hi], [int(r["start"]), int(r["end"])])
            gc[lo + x : lo + y] = value
    annotation_np.savez_compressed(
        annotation_S / "results/Q5_context_masks.npz", gene_context=gc, cpg_context=cc
    )
    panel = [
        p
        for p in annotation_json.loads((annotation_S / "plans/panel.json").read_text())
        if p["kind"] != "amplification_screen"
    ]
    rows = []
    for p in panel:
        rr = annotation_download("rmsk", p["chromosome"], p["start"], p["end"], p["locus_id"])
        intervals = [(int(r["genoStart"]), int(r["genoEnd"])) for r in rr]
        width = p["end"] - p["start"]
        z = annotation_np.load(annotation_S / "results/GBM-01" / (p["locus_id"] + ".npz"))
        cp = z["cpg_positions"][z["common_mask"]]
        repeatmask = annotation_np.array([any((a <= v < b for (a, b) in intervals)) for v in cp])
        ic = annotation_np.array(
            [any((a <= v < b for (a, b) in islands[p["chromosome"]])) for v in cp]
        )
        row = {
            "locus_id": p["locus_id"],
            "pair_id": p["pair_id"],
            "gene": p["gene"],
            "kind": p["kind"],
            "interval_bp": width,
            "repeat_overlap_bp": annotation_union_bp(intervals, p["start"], p["end"]),
            "repeat_common_CpG_fraction": float(repeatmask.mean()),
            "island_overlap_bp": annotation_union_bp(
                islands[p["chromosome"]], p["start"], p["end"]
            ),
            "island_common_CpG_fraction": float(ic.mean()),
            "common_CpGs": len(cp),
            "repeat_classes": ",".join(sorted(set((r["repClass"] for r in rr)))),
        }
        row["repeat_interval_fraction"] = row["repeat_overlap_bp"] / width
        row["island_interval_fraction"] = row["island_overlap_bp"] / width
        rows.append(row)
    _regional.regional_table(annotation_S / "results/Q5_panel_context.tsv", rows)
    return rows


def annotation_sample(s):
    label = s["study_label"]
    masks = annotation_np.load(annotation_S / "results/Q5_context_masks.npz")
    g = masks["gene_context"]
    c = masks["cpg_context"]
    common = annotation_np.load(annotation_D / "pooled_analysis/common_cpg_mask_depth5.npy")
    counts = annotation_np.load(annotation_D / "results" / label / "site_counts.npz")["counts"]
    total = annotation_np.zeros((19, 3), float)
    n = annotation_np.zeros(19, annotation_np.int64)
    for lo in range(0, len(common), 1000000):
        hi = min(lo + 1000000, len(common))
        keep = common[lo:hi]
        x = counts[lo:hi, :3][keep].astype(float)
        den = x.sum(1)
        assert annotation_np.all(den >= 5)
        x = x / den[:, None] * 100
        gi = g[lo:hi][keep]
        ci = c[lo:hi][keep]
        for code, offset, size in [(gi, 0, 3), (ci, 3, 4), (gi * 4 + ci, 7, 12)]:
            n[offset : offset + size] += annotation_np.bincount(code, minlength=size)
            for state in range(3):
                total[offset : offset + size, state] += annotation_np.bincount(
                    code, weights=x[:, state], minlength=size
                )
    assert annotation_np.all(n > 0)
    means = total / n[:, None]
    assert annotation_np.allclose(means.sum(1), 100)
    assert (
        n[:3].sum() == common.sum()
        and n[3:7].sum() == common.sum()
        and (n[7:].sum() == common.sum())
    )
    rows = []
    for k in range(19):
        axis = "gene" if k < 3 else "CpG" if k < 7 else "joint"
        name = (
            annotation_GENE[k]
            if k < 3
            else (
                annotation_CPG[k - 3]
                if k < 7
                else annotation_GENE[(k - 7) // 4] + " / " + annotation_CPG[(k - 7) % 4]
            )
        )
        rows.append(
            {
                "study_label": label,
                "diagnosis": s["diagnosis"],
                "context_axis": axis,
                "context": name,
                "common_CpGs": int(n[k]),
                "canonical_C_percent": float(means[k, 0]),
                "5mC_percent": float(means[k, 1]),
                "5hmC_percent": float(means[k, 2]),
                "combined_percent": float(means[k, 1] + means[k, 2]),
            }
        )
    _regional.regional_save(
        annotation_S / "validation" / f"Q5_{label}.json",
        {
            "all_context_partitions_equal_common_mask": True,
            "all_state_fractions_sum_to_100": True,
            "common_CpGs": int(common.sum()),
        },
    )
    print(label, "CONTEXT COMPLETE", flush=True)
    return rows


def annotation_summarize(rows, panel):
    _regional.regional_table(annotation_S / "results/Q5_specimen_context.tsv", rows)
    colors = {"Glioblastoma": "#B65046", "Meningioma": "#147D92"}
    rng = annotation_np.random.default_rng(20260912)
    effects = []
    for axis, names in [("gene", annotation_GENE), ("CpG", annotation_CPG)]:
        (fig, axes) = annotation_plt.subplots(1, 2, figsize=(12, 4.5))
        for ax, state in zip(axes, ["5hmC_percent", "5mC_percent"]):
            for i, name in enumerate(names):
                for j, diag in enumerate(colors):
                    rr = [
                        r
                        for r in rows
                        if r["context_axis"] == axis
                        and r["context"] == name
                        and (r["diagnosis"] == diag)
                    ]
                    v = [r[state] for r in rr]
                    ax.scatter(
                        i + (j - 0.5) * 0.25 + annotation_np.linspace(-0.045, 0.045, len(v)),
                        v,
                        c=colors[diag],
                        s=21,
                        alpha=0.8,
                        label=diag if i == 0 else None,
                    )
                    ax.plot(
                        [i + (j - 0.5) * 0.25 - 0.065, i + (j - 0.5) * 0.25 + 0.065],
                        [annotation_np.mean(v)] * 2,
                        c="black",
                        lw=2,
                    )
            ax.set_xticks(range(len(names)), names, rotation=15, ha="right")
            ax.set_ylabel(state.replace("_percent", "") + " (%)")
            ax.grid(axis="y", alpha=0.15)
        axes[0].legend(fontsize=8)
        fig.suptitle(
            "Modification patterns by "
            + ("gene context" if axis == "gene" else "CpG-island context")
            + "\nSame CpGs in every specimen; each dot is one specimen; black lines are group means",
            fontsize=11,
        )
        fig.tight_layout()
        _molecules.blocks_figure(
            fig, "Q5_01_gene_context" if axis == "gene" else "Q5_02_CpG_context"
        )
        for name in names:
            rr = [r for r in rows if r["context_axis"] == axis and r["context"] == name]
            x = annotation_np.array(
                [r["5hmC_percent"] for r in rr if r["diagnosis"] == "Glioblastoma"]
            )
            y = annotation_np.array(
                [r["5hmC_percent"] for r in rr if r["diagnosis"] == "Meningioma"]
            )
            boot = x[rng.integers(0, len(x), (5000, len(x)))].mean(1) - y[
                rng.integers(0, len(y), (5000, len(y)))
            ].mean(1)
            ci = annotation_np.percentile(boot, [2.5, 97.5])
            effects.append(
                {
                    "context_axis": axis,
                    "context": name,
                    "common_CpGs": rr[0]["common_CpGs"],
                    "GBM_mean_5hmC_percent": float(x.mean()),
                    "MEN_mean_5hmC_percent": float(y.mean()),
                    "GBM_minus_MEN_5hmC_pp": float(x.mean() - y.mean()),
                    "bootstrap_95_low": float(ci[0]),
                    "bootstrap_95_high": float(ci[1]),
                }
            )
    _regional.regional_table(annotation_S / "results/Q5_context_effects.tsv", effects)
    (fig, axes) = annotation_plt.subplots(1, 2, figsize=(9, 4.5))
    for ax, metric, title in zip(
        axes,
        ["repeat_common_CpG_fraction", "island_common_CpG_fraction"],
        ["Common CpGs in annotated repeats", "Common CpGs in CpG islands"],
    ):
        for i in range(1, 13):
            a = next((r for r in panel if r["locus_id"] == f"C{i:02d}"))
            b = next((r for r in panel if r["locus_id"] == f"R{i:02d}"))
            ax.plot([0, 1], [a[metric] * 100, b[metric] * 100], "-o", alpha=0.6, ms=4)
        ax.set_xticks([0, 1], ["Opposing-change\npromoters", "Matched comparison\npromoters"])
        ax.set_ylabel("Common CpGs (%)")
        ax.set_title(title)
        ax.set_ylim(-3, 103)
    fig.suptitle(
        "Genomic context of the 12 fixed promoter pairs\nAnnotation overlap only; each line is one pair",
        fontsize=11,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "Q5_03_matched_panel_context")
    summary = {
        "common_CpGs": sum((e["common_CpGs"] for e in effects if e["context_axis"] == "gene")),
        "context_effects": effects,
        "panel_context_means": {
            kind: {
                metric: float(annotation_np.mean([r[metric] for r in panel if r["kind"] == kind]))
                for metric in ["repeat_common_CpG_fraction", "island_common_CpG_fraction"]
            }
            for kind in ["candidate", "comparison"]
        },
        "scope": "Genome-wide common-CpG context plus descriptive repeat annotation restricted to fixed promoter panel",
    }
    _regional.regional_save(annotation_S / "results/Q5_summary.json", summary)
    lines = [
        "# Q5: genomic context",
        "",
        "## Scope and methods",
        "",
        "All 24,097,410 common CpGs were assigned to mutually exclusive gene contexts (promoter first, then gene body, otherwise intergenic) and CpG-island contexts (island; shore within 2kb; shelf 2–4kb; open sea beyond 4kb). Gene loci use the existing curated RefSeq snapshot. CpG islands and panel RepeatMasker intervals were retrieved from the UCSC hg38 API and archived with retrieval URLs and track dates. See https://genome.ucsc.edu/goldenPath/help/api.html. Transcript/locus promoter definitions and overlapping annotations affect these descriptive labels; they are not measured regulatory activity.",
        "",
        "Each specimen has an equal-CpG-weighted modification mean for each context, using exactly the same reference CpGs across specimens. The tables include 5000 within-diagnosis specimen bootstrap resamples; intervals are descriptive, pointwise and unadjusted. Reads or CpGs are not treated as independent patients. Gene and CpG categories each partition the common mask exactly; canonical C, 5mC and 5hmC fractions sum to 100% in every specimen/context.",
        "",
        "## Results",
        "",
    ]
    for e in effects:
        lines.append(
            f"- {e['context']}: GBM mean 5hmC {e['GBM_mean_5hmC_percent']:.2f}%, MEN {e['MEN_mean_5hmC_percent']:.2f}%; GBM-minus-MEN {e['GBM_minus_MEN_5hmC_pp']:+.2f} percentage points (descriptive interval {e['bootstrap_95_low']:+.2f} to {e['bootstrap_95_high']:+.2f}); {e['common_CpGs']:,} common CpGs."
        )
    lines += [
        "",
        "## Matched promoter context",
        "",
        annotation_json.dumps(summary["panel_context_means"], indent=2),
        "",
        "Repeat analysis here measures annotation overlap at the 24 fixed promoters, not genome-wide repeat-family enrichment or repeat-specific 5hmC accuracy. Existing mapping QC does not guarantee uniquely resolved repeat copies. The small, selected, matched panel cannot support general enrichment claims.",
        "",
        "## Biological limits",
        "",
        "Genomic location contextualizes the modification signal but does not measure gene expression, pathway activity, chromatin accessibility or active demethylation. GBM-versus-meningioma differences may include lineage and tissue-composition effects; no normal-brain comparison is available. This context analysis is exploratory and does not establish a tumour-specific regulatory mechanism.",
        "",
        "## Figures",
        "",
        "Q5_01_gene_context; Q5_02_CpG_context; Q5_03_matched_panel_context.",
    ]
    (annotation_S / "results/Q5_ANSWER.md").write_text("\n".join(lines) + "\n")
    _regional.regional_save(
        annotation_S / "results/Q5_complete.json",
        {"analysis_and_figures_generated": True, "requires_output_review": True},
    )
    print("Q5 COMPLETE", annotation_json.dumps(summary), flush=True)


def initialize_annotation():
    """Initialize the annotation stage once; load its declared inputs."""
    global annotation_ANNOT, annotation_CPG, annotation_D, annotation_GENE, annotation_PLAN, annotation_R, annotation_S
    if _runtime.initialized("story_context"):
        return
    _runtime.begin("story_context")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    annotation_R = annotation_Path(_config.workspace)
    annotation_D = annotation_R / ".analysis"
    annotation_S = annotation_D / "story_analysis"
    annotation_ANNOT = annotation_D / "reference/story_context_snapshot"
    annotation_ANNOT.mkdir(exist_ok=True)
    annotation_PLAN = {
        "genome": "hg38",
        "gene_annotation": "all curated RefSeq loci in existing annotation_features.tsv.gz, not filtered by modification differences",
        "gene_context_priority": ["promoter", "gene_body_outside_promoters", "intergenic"],
        "promoter_definition": "existing strand-aware -2000/+500bp locus TSS interval",
        "CpG_context_priority": ["island", "shore_0_2kb", "shelf_2_4kb", "open_sea_beyond_4kb"],
        "main_aggregation": "equal weight per common CpG, same >=5 passing-call mask in every specimen",
        "replicate": "specimen; 15 GBM and 5 MEN",
        "uncertainty": "5000 specimen bootstrap draws within diagnosis; descriptive pointwise intervals",
        "repeat_scope": "24 fixed promoters only, paired candidate/comparison interval overlap; no genome-wide repeat enrichment claim",
        "sources": [
            "https://api.genome.ucsc.edu",
            "https://genome.ucsc.edu/goldenPath/help/api.html",
        ],
    }
    annotation_GENE = ["Promoter", "Gene body outside promoters", "Intergenic"]
    annotation_CPG = ["CpG island", "Shore (0–2 kb)", "Shelf (2–4 kb)", "Open sea (>4 kb)"]
    _runtime.finish("story_context")


def run_annotation():
    """Execute the annotation workflow stage."""
    initialize_annotation()
    global annotation_panel, annotation_pool, annotation_pp, annotation_rs, annotation_samples
    annotation_pp = annotation_S / "plans/Q5_plan.json"
    if annotation_pp.exists():
        assert annotation_json.loads(annotation_pp.read_text()) == annotation_PLAN
    else:
        _regional.regional_save(annotation_pp, annotation_PLAN)
    assert annotation_union_bp([(0, 10), (5, 20), (25, 30)], 3, 28) == 20
    annotation_panel = annotation_annotations()
    annotation_samples = sorted(
        annotation_json.loads((annotation_D / "cohort.json").read_text())["samples"],
        key=lambda s: s["study_label"],
    )
    with _runtime.process_pool(max_workers=_config.workers) as annotation_pool:
        annotation_rs = list(annotation_pool.map(annotation_sample, annotation_samples))
    annotation_summarize([r for rr in annotation_rs for r in rr], annotation_panel)


# REPEATS


from pathlib import Path as repeats_Path
import csv as repeats_csv, gzip as repeats_gzip, json as repeats_json, numpy as repeats_np
import matplotlib.pyplot as repeats_plt


def repeats_run():
    plan = {
        "analysis": "post-context sensitivity, specified after observing panel repeat overlap",
        "selection": "retain originally chosen Q3 blocks only if all three full CpG dyads avoid RepeatMasker intervals",
        "minimum_remaining_blocks_per_locus": 3,
        "minimum_matched_pairs_per_specimen": 6,
        "metric": "same precomputed 12-molecule co-occurrence score; no reselection or new permutations",
    }
    _regional.regional_save(repeats_S / "plans/Q5_repeat_sensitivity_plan.json", plan)
    panel = [
        p
        for p in repeats_json.loads((repeats_S / "plans/panel.json").read_text())
        if p["kind"] != "amplification_screen"
    ]
    intervals = {}
    for p in panel:
        with repeats_gzip.open(
            repeats_D / "reference/story_context_snapshot" / ("rmsk_" + p["locus_id"] + ".json.gz"),
            "rt",
        ) as f:
            rr = repeats_json.load(f)["rmsk"]
        intervals[p["locus_id"]] = [(r["genoStart"], r["genoEnd"]) for r in rr]
    blocks = list(
        repeats_csv.DictReader((repeats_S / "results/Q3_block_metrics.tsv").open(), delimiter="\t")
    )
    kept = []
    for b in blocks:
        positions = list(map(int, b["CpG_positions"].split(",")))
        retain = not any(
            (a < p + 2 and e > p for p in positions for (a, e) in intervals[b["locus_id"]])
        )
        b["all_CpG_dyads_outside_RepeatMasker"] = retain
        if retain and b["eligible"] == "True":
            kept.append(b)
    _regional.regional_table(repeats_S / "results/Q5_repeat_block_eligibility.tsv", blocks)
    original = list(
        repeats_csv.DictReader(
            (repeats_S / "results/Q3_specimen_metrics.tsv").open(), delimiter="\t"
        )
    )
    rows = []
    for s in original:
        label = s["study_label"]
        means = {}
        for p in panel:
            b = [r for r in kept if r["study_label"] == label and r["locus_id"] == p["locus_id"]]
            if len(b) >= 3:
                means[p["locus_id"]] = repeats_np.mean(
                    [float(r["joint_5hmC_excess_pp"]) for r in b]
                )
        pairs = [
            (means[f"C{i:02d}"], means[f"R{i:02d}"])
            for i in range(1, 13)
            if f"C{i:02d}" in means and f"R{i:02d}" in means
        ]
        difference = (
            float(repeats_np.mean([a - b for (a, b) in pairs])) if len(pairs) >= 6 else None
        )
        rows.append(
            {
                "study_label": label,
                "diagnosis": s["diagnosis"],
                "original_candidate_minus_comparison_pp": float(s["candidate_joint_excess_pp"])
                - float(s["comparison_joint_excess_pp"]),
                "nonrepeat_matched_pairs": len(pairs),
                "eligible_nonrepeat_comparison": len(pairs) >= 6,
                "nonrepeat_candidate_minus_comparison_pp": difference,
            }
        )
    _regional.regional_table(repeats_S / "results/Q5_repeat_sensitivity.tsv", rows)
    rng = repeats_np.random.default_rng(20260912)
    summary = {}
    for diag in ["Glioblastoma", "Meningioma"]:
        v = repeats_np.array(
            [
                r["nonrepeat_candidate_minus_comparison_pp"]
                for r in rows
                if r["diagnosis"] == diag and r["eligible_nonrepeat_comparison"]
            ]
        )
        ci = (
            repeats_np.percentile(v[rng.integers(0, len(v), (10000, len(v)))].mean(1), [2.5, 97.5])
            if len(v) >= 3
            else [None, None]
        )
        summary[diag] = {
            "eligible_specimens": len(v),
            "positive_paired_differences": int((v > 0).sum()),
            "mean_nonrepeat_difference_pp": float(v.mean()) if len(v) else None,
            "bootstrap_95_low": None if ci[0] is None else float(ci[0]),
            "bootstrap_95_high": None if ci[1] is None else float(ci[1]),
        }
    (fig, axes) = repeats_plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, diag, col in zip(axes, summary, ["#B65046", "#147D92"]):
        rr = [r for r in rows if r["diagnosis"] == diag and r["eligible_nonrepeat_comparison"]]
        for r in rr:
            ax.plot(
                [0, 1],
                [
                    r["original_candidate_minus_comparison_pp"],
                    r["nonrepeat_candidate_minus_comparison_pp"],
                ],
                "-o",
                c=col,
                alpha=0.7,
                ms=4,
            )
        ax.axhline(0, c="gray", lw=0.8)
        ax.set_xticks([0, 1], ["Original blocks", "Outside repeats"])
        ax.set_ylabel("Candidate minus comparison co-occurrence (pp)")
        ax.set_title(diag + f" ({len(rr)} eligible specimens)")
    fig.suptitle(
        "Sensitivity to repeat overlap: same original blocks, excluding annotated repeats\nRequires at least 3 remaining blocks per promoter and 6 matched pairs per specimen",
        fontsize=10,
    )
    fig.tight_layout()
    _molecules.blocks_figure(fig, "Q5_04_nonrepeat_molecular_sensitivity")
    _regional.regional_save(repeats_S / "results/Q5_repeat_sensitivity_summary.json", summary)
    print(repeats_json.dumps(summary), flush=True)


def initialize_repeats():
    """Initialize the repeats stage once; load its declared inputs."""
    global repeats_D, repeats_R, repeats_S
    if _runtime.initialized("story_repeat_sensitivity"):
        return
    _runtime.begin("story_repeat_sensitivity")
    _runtime.initialize("pooled_analysis")
    _runtime.initialize("story_molecules")
    repeats_R = repeats_Path(_config.workspace)
    repeats_D = repeats_R / ".analysis"
    repeats_S = repeats_D / "story_analysis"
    _runtime.finish("story_repeat_sensitivity")


def run_repeats():
    """Execute the repeats workflow stage."""
    initialize_repeats()
    repeats_run()
