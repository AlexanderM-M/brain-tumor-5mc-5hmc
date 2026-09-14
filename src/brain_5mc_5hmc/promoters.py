"""Promoters analysis stages for the brain tumour 5mC/5hmC study.

Each stage has explicit initialization and execution; importing is read-only.
See docs/workflow.md for inputs, stage order and reproduction limits.
"""

from . import runtime as _runtime

from .runtime import config as _config

from . import regional as _regional


# PANEL


from pathlib import Path as panel_Path
import csv as panel_csv, gzip as panel_gzip, json as panel_json, hashlib as panel_hashlib, numpy as panel_np, pysam as panel_pysam


def panel_strong(r):
    return (
        r["opposing_changes_with_small_combined_difference"] == "True"
        and r["depth10_same_5hmC_direction"] == "True"
        and (float(r["pairwise_comparisons_agreeing_with_mean_direction"]) >= 0.8)
        and (
            float(r["GBM_minus_meningioma_5hmC_pp"]) * float(r["median_group_5hmC_difference_pp"])
            > 0
        )
    )


def panel_separate(r, ss):
    return all(
        (
            r["chromosome"] != x["chromosome"] or abs(int(r["start"]) - int(x["start"])) >= 100000
            for x in ss
        )
    )


def initialize_panel():
    """Initialize the panel stage once; load its declared inputs."""
    global panel_D, panel_O, panel_R, panel_S, panel_annotation, panel_b, panel_c, panel_candidates, panel_controls, panel_distance, panel_f, panel_fa, panel_features, panel_fi, panel_fs, panel_i, panel_kind, panel_mid, panel_name, panel_panel, panel_plan, panel_pool, panel_prom, panel_r, panel_rows, panel_rr, panel_selected, panel_seq, panel_used
    if _runtime.initialized("story_panel"):
        return
    _runtime.begin("story_panel")
    _runtime.initialize("pooled_analysis")
    panel_R = panel_Path(_config.workspace)
    panel_D = panel_R / ".analysis"
    panel_O = panel_D / "pooled_analysis"
    panel_S = panel_D / "story_analysis"
    with panel_gzip.open(panel_O / "tables/pooled_gene_promoter_results.tsv.gz", "rt") as panel_f:
        panel_rows = list(panel_csv.DictReader(panel_f, delimiter="\t"))
    panel_features = panel_json.loads((panel_O / "cache/features5.json").read_text())
    panel_fi = {x["feature_id"]: i for (i, x) in enumerate(panel_features)}
    panel_b = panel_np.load(panel_O / "pooled_profiles.npz")["features5"]
    panel_fa = panel_pysam.FastaFile(_config.reference)
    panel_prom = [
        r
        for r in panel_rows
        if r["feature"] == "promoter"
        and int(r["common_CpGs"]) >= 40
        and (r["depth10_eligible"] == "True")
    ]
    for panel_r in panel_prom:
        panel_i = panel_fi[panel_r["feature_id"]]
        panel_seq = panel_fa.fetch(
            panel_r["chromosome"], int(panel_r["start"]), int(panel_r["end"])
        ).upper()
        panel_r["GC_fraction"] = sum((panel_seq.count(c) for c in "GC")) / len(panel_seq)
        panel_r["pooled_5hmC"] = float(panel_b[:, panel_i, 1].mean())
        panel_r["pooled_combined"] = float(panel_b[:, panel_i, :].sum(1).mean())
    panel_candidates = sorted(
        [r for r in panel_prom if panel_strong(r)],
        key=lambda r: -abs(float(r["GBM_minus_meningioma_5hmC_pp"])),
    )
    panel_selected = []
    for panel_r in panel_candidates:
        if panel_separate(panel_r, panel_selected):
            panel_selected.append(panel_r)
        if len(panel_selected) == 12:
            break
    panel_controls = [
        r
        for r in panel_prom
        if abs(float(r["GBM_minus_meningioma_5hmC_pp"])) < 1
        and abs(float(r["GBM_minus_meningioma_combined_pp"])) < 2
    ]
    panel_panel = []
    panel_used = list(panel_selected)
    for panel_i, panel_r in enumerate(panel_selected, 1):
        panel_pool = [c for c in panel_controls if panel_separate(c, panel_used)]

        def panel_distance(c):
            return (
                (panel_np.log(int(c["common_CpGs"]) / int(panel_r["common_CpGs"])) / 0.5) ** 2
                + ((c["GC_fraction"] - panel_r["GC_fraction"]) / 0.1) ** 2
                + ((c["pooled_combined"] - panel_r["pooled_combined"]) / 10) ** 2
                + ((c["pooled_5hmC"] - panel_r["pooled_5hmC"]) / 5) ** 2
            )

        panel_c = min(panel_pool, key=panel_distance)
        panel_used.append(panel_c)
        for panel_kind, panel_rr in [("candidate", panel_r), ("comparison", panel_c)]:
            panel_panel.append(
                {
                    "locus_id": f"{('C' if panel_kind == 'candidate' else 'R')}{panel_i:02d}",
                    "pair_id": panel_i,
                    "kind": panel_kind,
                    "gene": panel_rr["gene"],
                    "chromosome": panel_rr["chromosome"],
                    "start": int(panel_rr["start"]),
                    "end": int(panel_rr["end"]),
                    "strand": panel_rr["strand"],
                    "feature_id": panel_rr["feature_id"],
                    "common_CpGs": int(panel_rr["common_CpGs"]),
                    "GC_fraction": panel_rr["GC_fraction"],
                    "pooled_5hmC": panel_rr["pooled_5hmC"],
                    "pooled_combined": panel_rr["pooled_combined"],
                    "5hmC_group_difference_pp": float(panel_rr["GBM_minus_meningioma_5hmC_pp"]),
                    "matching_distance": float(panel_distance(panel_c)),
                }
            )
    panel_annotation = []
    with panel_gzip.open(panel_O / "tables/annotation_features.tsv.gz", "rt") as panel_f:
        panel_annotation = list(panel_csv.DictReader(panel_f, delimiter="\t"))
    for panel_name in ["EGFR", "PDGFRA", "CDK4", "MDM2"]:
        panel_fs = [
            f for f in panel_annotation if f["gene"] == panel_name and f["feature"] == "promoter"
        ]
        assert panel_fs, panel_name
        panel_f = max(panel_fs, key=lambda x: int(x["transcript_count"]))
        panel_mid = (int(panel_f["start"]) + int(panel_f["end"])) // 2
        panel_panel.append(
            {
                "locus_id": "A_" + panel_name,
                "pair_id": 0,
                "kind": "amplification_screen",
                "gene": panel_name,
                "chromosome": panel_f["chromosome"],
                "start": max(0, panel_mid - 10000),
                "end": panel_mid + 10000,
                "strand": panel_f["strand"],
                "feature_id": panel_f["feature_id"],
                "common_CpGs": 0,
                "GC_fraction": 0.0,
                "pooled_5hmC": 0.0,
                "pooled_combined": 0.0,
                "5hmC_group_difference_pp": 0.0,
                "matching_distance": 0.0,
            }
        )
    assert len(panel_panel) == 28
    _regional.regional_save(panel_S / "plans/panel.json", panel_panel)
    _regional.regional_table(panel_S / "plans/panel.tsv", panel_panel)
    (panel_S / "plans/panel.bed").write_text(
        "".join(
            (
                f"{r['chromosome']}\t{r['start']}\t{r['end']}\t{r['locus_id']}\n"
                for r in sorted(panel_panel, key=lambda x: (int(x["chromosome"][3:]), x["start"]))
            )
        )
    )
    panel_plan = {
        "question": "Q3 molecule-level structure behind opposing 5mC/5hmC contrasts",
        "panel": "12 candidate promoters, 12 matched comparison promoters; 4 prespecified amplification-screen intervals reserved for Q4",
        "selection": "Candidate promoters >=40 common CpGs; >=80% pairwise direction agreement; mean and median agree; depth10 direction agrees; largest absolute 5hmC difference, separated by >=100kb",
        "matching": "Comparison promoter abs 5hmC difference <1 pp, abs combined difference <2 pp; nearest weighted distance in CpG count, GC, pooled combined and pooled 5hmC; >=100kb separation",
        "read_filters": "Primary mapped autosomal reads, MAPQ 20-254, qs>=10, no duplicate/QC-fail flag, exact read-ID deduplication, MM+ML required; trim100bp; exact reference CpG; probability>=0.8",
        "molecule_analysis": "Fixed 3-CpG blocks selected by coordinates and common-site coverage, not by read patterns; >=12 completely observed molecules per block/specimen; all comparisons specimen-level; depth-matched resampling and per-site-frequency-preserving permutations",
        "interpretation": "Exploratory within the same cohort; reads are molecules, not cell identities or clones; no biological null conclusion from failed coverage",
        "downstream": "Keep targeted BAMs for allele and independent depth feasibility; do not infer copy number from modification count alone",
    }
    _regional.regional_save(panel_S / "plans/Q3_plan.json", panel_plan)
    print([(r["locus_id"], r["gene"]) for r in panel_panel])
    print(
        "Panel SHA256",
        panel_hashlib.sha256((panel_S / "plans/panel.json").read_bytes()).hexdigest(),
    )
    _runtime.finish("story_panel")


def run_panel():
    """Execute the panel workflow stage."""
    initialize_panel()


# SELECTION


from pathlib import Path as selection_Path
import json as selection_json, numpy as selection_np, pysam as selection_pysam


def selection_setup():
    a = selection_np.load(selection_O / "pooled_profiles.npz")
    f5 = selection_json.loads((selection_O / "cache/features5.json").read_text())
    f10 = selection_json.loads((selection_O / "cache/features10.json").read_text())
    m10 = {r["feature_id"]: i for (i, r) in enumerate(f10)}
    ix = [
        i
        for (i, r) in enumerate(f5)
        if r["feature"] == "promoter" and r["common_CpGs"] >= 40 and (r["feature_id"] in m10)
    ]
    features = [dict(f5[i]) for i in ix]
    b = a["features5"][:, ix, :]
    b10 = a["features10"][:, [m10[r["feature_id"]] for r in features], :]
    labels = a["labels"].tolist()
    diag = a["diagnoses"].tolist()
    with selection_pysam.FastaFile(_config.reference) as fa:
        for r in features:
            seq = fa.fetch(r["chromosome"], r["start"], r["end"]).upper()
            r["GC_fraction"] = (seq.count("G") + seq.count("C")) / len(seq)
    return (features, b, b10, labels, selection_np.array(diag) == "Glioblastoma")


def selection_separate(i, used, features):
    r = features[i]
    return all(
        (
            r["chromosome"] != features[j]["chromosome"]
            or abs(r["start"] - features[j]["start"]) >= 100000
            for j in used
        )
    )


def selection_select(features, b, b10, group, train):
    """Select and sequentially match promoters using training specimens only.

    The caller supplies a fixed coverage-qualified feature universe. Held-out
    modification values never enter candidate ranking or comparison matching."""
    gb = selection_np.array([i for i in train if group[i]])
    me = selection_np.array([i for i in train if not group[i]])
    delta = b[gb].mean(0) - b[me].mean(0)
    d10 = b10[gb].mean(0) - b10[me].mean(0)
    h = delta[:, 1]
    m = delta[:, 0]
    stable = selection_np.ones(len(h), bool)
    for drop in train:
        gg = gb[gb != drop]
        mm = me[me != drop]
        dh = b[gg, :, 1].mean(0) - b[mm, :, 1].mean(0)
        stable &= dh * h > 0
    median = selection_np.median(b[gb, :, 1], axis=0) - selection_np.median(b[me, :, 1], axis=0)
    agree = selection_np.mean(
        (b[gb, :, 1][:, None, :] - b[me, :, 1][None, :, :]) * h[None, None, :] > 0, axis=(0, 1)
    )
    strong = (
        (abs(h) >= 2)
        & stable
        & (h * m < 0)
        & (abs(h + m) <= 1)
        & (d10[:, 1] * h > 0)
        & (agree >= 0.8)
        & (median * h > 0)
    )
    rank = sorted(
        selection_np.flatnonzero(strong), key=lambda i: (-abs(h[i]), features[i]["feature_id"])
    )
    chosen = []
    for i in rank:
        if selection_separate(i, chosen, features):
            chosen.append(int(i))
        if len(chosen) == 12:
            break
    controls = selection_np.flatnonzero((abs(h) < 1) & (abs(h + m) < 2))
    used = chosen.copy()
    pooled = b[train].mean(0)
    pairs = []
    for i in chosen:
        pool = [int(j) for j in controls if selection_separate(j, used, features)]

        def distance(j):
            return (
                (selection_np.log(features[j]["common_CpGs"] / features[i]["common_CpGs"]) / 0.5)
                ** 2
                + ((features[j]["GC_fraction"] - features[i]["GC_fraction"]) / 0.1) ** 2
                + ((pooled[j].sum() - pooled[i].sum()) / 10) ** 2
                + ((pooled[j, 1] - pooled[i, 1]) / 5) ** 2
            )

        if not pool:
            continue
        j = min(pool, key=lambda j: (distance(j), features[j]["feature_id"]))
        used.append(j)
        pairs.append(
            {
                "candidate_index": i,
                "comparison_index": j,
                "matching_distance": float(distance(j)),
                "training_5hmC_difference_pp": float(h[i]),
                "training_combined_difference_pp": float(h[i] + m[i]),
            }
        )
    return pairs


def selection_run():
    (features, b, b10, labels, group) = selection_setup()
    samples = selection_json.loads((selection_D / "cohort.json").read_text())["samples"]
    lookup = {r["study_label"]: r for r in samples}
    patient_groups = [lookup[l]["patient_group"] for l in labels]
    assert len(set(patient_groups)) == 20
    plan = {
        "scope": "Only follow-up questions 1 and 2; no additional cohort/public-data validation",
        "Q1": "leave-one-recorded-patient-out; all 20 recorded groups distinct; reselect and match promoters with training modification values only",
        "fixed_QC_universe": "Existing all-specimen common CpGs and >=10-call feature availability are fixed unsupervised coverage filters; held-out coverage therefore contributes to this eligibility universe, but no held-out modification values or diagnosis contrasts enter selection. Internal validation conditional on this universe, not fully prospective validation.",
        "candidate_rule": "training abs h difference>=2pp; opposite m difference; abs combined<=1pp; direction stable after all training-only single omissions; >=80% training cross-diagnosis pairs agree; median agrees; >=10-call direction agrees; top12 separated >=100kb",
        "matching": "training-only means of h and combined, plus fixed GC and common-CpG count; same matching scale as original; control training abs h<1pp and abs combined<2pp",
        "blind_panel": "120 promoters, uniform seeded shuffle of coverage-qualified universe, >=100kb apart, no diagnosis/profile-value selection",
        "Q1_molecular_metric": "same fixed 3-CpG blocks and 12-molecule resampling as original; >=3 eligible blocks/locus and >=6 matched pairs/specimen",
        "Q1_inference": "held-out specimen scores and selection stability; no fold-independent or independent-cohort claim",
        "Q2_distance_bins_bp": [1, 10, 25, 50, 100, 250, 500, 1000, 2500],
        "Q2_sampling": "Up to 40 coordinate-selected CpG pairs per locus/distance bin; >=12 complete molecules; 64 fixed-depth resamples of 12 molecules; same pairs for diagnoses",
        "Q2_metrics": [
            "5hmC joint occurrence minus exact column-permutation expectation",
            "binary phi correlation at polymorphic pairs; degeneracy reported",
            "covariance of h minus m state coding, supplementary",
        ],
        "Q2_controls": "within-specimen distance/abundance-stratified candidate-comparison overlap weighting; minimum overlap support reported; missingness preserved by complete-pair analysis; both original and held-out panels; diagnosis-blind panel breadth",
        "Q2_limits": "phi remains bounded by marginals; distance is not causality; row-wide state mixture can generate long-range association; no clones or active demethylation claim",
        "seed": 20260912,
    }
    pp = selection_F / "plans/protocol.json"
    if pp.exists():
        assert selection_json.loads(pp.read_text()) == plan
    else:
        _regional.regional_save(pp, plan)
    rng = selection_np.random.default_rng(plan["seed"])
    blind = []
    for i in rng.permutation(len(features)):
        if selection_separate(i, blind, features):
            blind.append(int(i))
        if len(blind) == 120:
            break
    folds = []
    indices = set(blind)
    for heldout, label in enumerate(labels):
        train = [i for i in range(20) if patient_groups[i] != patient_groups[heldout]]
        pairs = selection_select(features, b, b10, group, train)
        folds.append(
            {"held_out": label, "training_labels": [labels[i] for i in train], "pairs": pairs}
        )
        indices.update((v[k] for v in pairs for k in ["candidate_index", "comparison_index"]))
        print(label, "selected", len(pairs), "pairs", flush=True)
    checks = []
    for heldout in [0, 15]:
        altered = b.copy()
        altered10 = b10.copy()
        altered[heldout] += 1234
        altered10[heldout] -= 1234
        train = [i for i in range(20) if i != heldout]
        assert (
            selection_select(features, altered, altered10, group, train) == folds[heldout]["pairs"]
        )
        checks.append(
            {
                "held_out": labels[heldout],
                "selection_unchanged_after_extreme_held_out_value_perturbation": True,
            }
        )
    original = [
        p
        for p in selection_json.loads((selection_D / "story_analysis/plans/panel.json").read_text())
        if p["kind"] != "amplification_screen"
    ]
    feature_index = {r["feature_id"]: i for (i, r) in enumerate(features)}
    for p in original:
        indices.add(feature_index[p["feature_id"]])
    panel = []
    ids = {}
    for n, i in enumerate(
        sorted(
            indices,
            key=lambda i: (
                int(features[i]["chromosome"][3:]),
                features[i]["start"],
                features[i]["feature_id"],
            ),
        ),
        1,
    ):
        p = dict(features[i])
        p.update(
            locus_id=f"P{n:04d}",
            universe_index=i,
            blind_panel=i in blind,
            original_locus_id=next(
                (r["locus_id"] for r in original if r["feature_id"] == p["feature_id"]), ""
            ),
        )
        panel.append(p)
        ids[i] = p["locus_id"]
    for fold in folds:
        for pair in fold["pairs"]:
            for kind in ["candidate", "comparison"]:
                pair[kind + "_locus_id"] = ids[pair[kind + "_index"]]
    _regional.regional_save(selection_F / "plans/panel.json", panel)
    _regional.regional_table(selection_F / "plans/panel.tsv", panel)
    _regional.regional_save(selection_F / "plans/folds.json", folds)
    _regional.regional_save(selection_F / "plans/blind_panel.json", [ids[i] for i in blind])
    _regional.regional_save(selection_F / "validation/selection_leakage_checks.json", checks)
    _regional.regional_save(
        selection_F / "validation/patient_group_audit.json",
        {
            "specimens": 20,
            "distinct_recorded_patient_groups": 20,
            "linkage_basis": sorted(set((s["patient_linkage_basis"] for s in samples))),
            "no_personal_identifiers_exported": True,
        },
    )
    (selection_F / "plans/panel.bed").write_text(
        "".join((f"{p['chromosome']}\t{p['start']}\t{p['end']}\t{p['locus_id']}\n" for p in panel))
    )
    print("TOTAL PANEL", len(panel), "BLIND", len(blind), flush=True)


def initialize_selection():
    """Initialize the selection stage once; load its declared inputs."""
    global selection_D, selection_F, selection_O, selection_R
    if _runtime.initialized("followup_select"):
        return
    _runtime.begin("followup_select")
    _runtime.initialize("pooled_analysis")
    selection_R = selection_Path(_config.workspace)
    selection_D = selection_R / ".analysis"
    selection_O = selection_D / "pooled_analysis"
    selection_F = selection_D / "focused_followup"
    _runtime.finish("followup_select")


def run_selection():
    """Execute the selection workflow stage."""
    initialize_selection()
    selection_run()


# ALLOCATION PLAN


from pathlib import Path as allocation_plan_Path
import json as allocation_plan_json, itertools as allocation_plan_itertools, time as allocation_plan_time, os as allocation_plan_os, numpy as allocation_plan_np


def allocation_plan_assignment(task):
    (index, men) = task
    path = allocation_plan_S / "plans/permutation_folds" / f"{index:04d}.json"
    if path.exists():
        return allocation_plan_json.loads(path.read_text())
    group = allocation_plan_np.ones(20, bool)
    group[list(men)] = False
    folds = []
    for held in range(20):
        pairs = selection_select(
            allocation_plan_FEATURES,
            allocation_plan_B,
            allocation_plan_B10,
            group,
            [i for i in range(20) if i != held],
        )
        folds.append(
            {"held_out_index": held, "held_out": allocation_plan_LABELS[held], "pairs": pairs}
        )
    result = {"assignment_index": index, "meningioma_indices": list(men), "folds": folds}
    _regional.regional_save(path, result)
    print("SELECTION", index, "pairs", sum((len(f["pairs"]) for f in folds)), flush=True)
    return result


def allocation_plan_run():
    global allocation_plan_FEATURES, allocation_plan_B, allocation_plan_B10, allocation_plan_LABELS
    (
        allocation_plan_FEATURES,
        allocation_plan_B,
        allocation_plan_B10,
        allocation_plan_LABELS,
        group,
    ) = selection_setup()
    observed = tuple(allocation_plan_np.flatnonzero(~group).tolist())
    allassign = [x for x in allocation_plan_itertools.combinations(range(20), 5) if x != observed]
    rng = allocation_plan_np.random.default_rng(20260913)
    chosen = [allassign[i] for i in rng.choice(len(allassign), 499, replace=False)]
    alloc = [observed] + chosen
    plan = {
        "scope": "Three authorised strengthening analyses; no final main/supplement allocation before discussion",
        "P1": "delete each candidate/comparison promoter pair; held-out panels, overall and exploratory 25-249bp residual; >=6 remaining pairs per patient; report sign changes and eligibility",
        "P2": "same physical CpG pairs, reads, missing mask and 12-molecule covariance scale for 5hmC,5mC,combined; separate binary row/column-preserving null per encoding; joint-support comparisons and limited randomisability explicit",
        "P3_assignments": 500,
        "P3_random_assignments": 499,
        "P3_seed": 20260913,
        "P3_label_space": 15504,
        "P3_reselection": "Repeat ALL 20 patient-held-out selections and comparison matching for every assignment, using original feature universe and thresholds; do not restrict to current selected promoters",
        "P3_primary_family": "5hmC row/column-null residual only; maximum absolute Welch studentized group contrast over all 36 contiguous unions of eight fixed distance bins; covers the previously explored 25-249bp interval",
        "P3_patient_score": "mean matched candidate-control residual, equal bins within pair then equal pairs; same available bins for both members; >=6 eligible pairs per patient",
        "P3_statistic_eligibility": "at least 3 pseudo-meningiomas and 10 pseudo-GBMs; otherwise statistic=0, representing failure of the selection-and-scoring procedure; no null allocations discarded",
        "P3_p_value": "(1 + number of random-assignment maxima >= observed statistic) / 500; report omnibus and max-family-adjusted values for each interval",
        "P3_null_scope": "exchangeability of diagnosis labels for the entire procedure conditional on fixed coverage universe, reads, and simulation seeds; not a causal or isolated modification-specific null",
        "null_randomisation": "two row/column-preserving binary chains; 50 trades/read burn-in; 32 retained states each separated by 5 trades/read; >=64 nontrivial trades each; >=5 CpG pairs per locus/bin; reads>=10 observed CpGs; pairs>=12 complete reads",
        "limitations": "MC resolution 0.002; retrospective hypotheses; modification encodings descriptive; strict selection may fail under null, which is explicitly part of the tested procedure; no post-hoc conditioning on successful permutations",
    }
    p = allocation_plan_S / "plans/protocol.json"
    if p.exists():
        assert allocation_plan_json.loads(p.read_text()) == plan
    else:
        _regional.regional_save(p, plan)
    (allocation_plan_S / "plans/permutation_folds").mkdir(exist_ok=True)
    _regional.regional_save(
        allocation_plan_S / "plans/feature_universe.json", allocation_plan_FEATURES
    )
    _regional.regional_save(
        allocation_plan_S / "plans/assignments.json",
        [{"assignment_index": i, "meningioma_indices": list(a)} for (i, a) in enumerate(alloc)],
    )
    _regional.regional_save(
        allocation_plan_S / "selection_job.json",
        {
            "pid": allocation_plan_os.getpid(),
            "workers": _config.workers,
            "started": allocation_plan_time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", allocation_plan_time.gmtime()
            ),
        },
    )
    with _runtime.process_pool(max_workers=_config.workers) as pool:
        results = list(pool.map(allocation_plan_assignment, list(enumerate(alloc))))
    needed = {i: set() for i in range(20)}
    for a in results:
        for fold in a["folds"]:
            for pair in fold["pairs"]:
                needed[fold["held_out_index"]].update(
                    [pair["candidate_index"], pair["comparison_index"]]
                )
    for p in allocation_plan_json.loads((allocation_plan_F / "plans/panel.json").read_text()):
        if p["original_locus_id"]:
            for i in needed:
                needed[i].add(p["universe_index"])
    allneeded = sorted(set.union(*needed.values()))
    existing = {
        p["universe_index"]: p
        for p in allocation_plan_json.loads((allocation_plan_F / "plans/panel.json").read_text())
    }
    panel = []
    for i in allneeded:
        p = dict(allocation_plan_FEATURES[i])
        p.update(
            locus_id=f"U{i:05d}",
            universe_index=i,
            existing_locus_id=existing.get(i, {}).get("locus_id", ""),
        )
        panel.append(p)
    _regional.regional_save(allocation_plan_S / "plans/panel.json", panel)
    _regional.regional_table(allocation_plan_S / "plans/panel.tsv", panel)
    _regional.regional_save(
        allocation_plan_S / "plans/needed_by_patient.json",
        {allocation_plan_LABELS[i]: sorted(v) for (i, v) in needed.items()},
    )
    new = [p for p in panel if not p["existing_locus_id"]]
    _regional.regional_save(allocation_plan_S / "plans/new_panel.json", new)
    (allocation_plan_S / "plans/panel.bed").write_text(
        "".join(
            (
                f"{p['chromosome']}\t{p['start']}\t{p['end']}\t{p['locus_id']}\n"
                for p in sorted(new, key=lambda p: (int(p["chromosome"][3:]), p["start"]))
            )
        )
    )
    _regional.regional_save(
        allocation_plan_S / "results/selection_complete.json",
        {
            "assignments": 500,
            "universe_promoters": len(allocation_plan_FEATURES),
            "selected_union_promoters": len(panel),
            "new_promoters": len(new),
            "patient_locus_scores_needed": sum(map(len, needed.values())),
        },
    )
    print(
        "SELECTION COMPLETE",
        len(panel),
        "promoters;",
        len(new),
        "new;",
        sum(map(len, needed.values())),
        "patient-locus scores",
        flush=True,
    )


def initialize_allocation_plan():
    """Initialize the allocation_plan stage once; load its declared inputs."""
    global allocation_plan_B, allocation_plan_B10, allocation_plan_D, allocation_plan_F, allocation_plan_FEATURES, allocation_plan_LABELS, allocation_plan_R, allocation_plan_S
    if _runtime.initialized("strengthen_plan"):
        return
    _runtime.begin("strengthen_plan")
    _runtime.initialize("followup_select")
    _runtime.initialize("pooled_analysis")
    allocation_plan_R = allocation_plan_Path(_config.workspace)
    allocation_plan_D = allocation_plan_R / ".analysis"
    allocation_plan_F = allocation_plan_D / "focused_followup"
    allocation_plan_S = allocation_plan_D / "strengthening"
    allocation_plan_FEATURES = allocation_plan_B = allocation_plan_B10 = allocation_plan_LABELS = (
        None
    )
    _runtime.finish("strengthen_plan")


def run_allocation_plan():
    """Execute the allocation_plan workflow stage."""
    initialize_allocation_plan()
    allocation_plan_run()
