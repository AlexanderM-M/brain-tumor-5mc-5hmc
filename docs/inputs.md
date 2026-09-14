# Inputs and output layout

All coordinates are hg38, 0-based; intervals are half-open unless an output explicitly states otherwise. Reference CpG positions identify the forward-strand C. Counts sum evidence from both strands. Fractions in aggregate profiles and tables are percentages, not values between zero and one.

## Cohort and read manifests (private)

`data/cohort/cohort.json` contains a `samples` list. Required fields depend on the entry point:

```json
{
  "samples": [
    {
      "study_label": "SYNTHETIC-01",
      "diagnosis": "Glioblastoma",
      "patient_group": "SYNTHETIC-PATIENT-01",
      "manifest": "/private/workspace/data/cohort/manifests/SYNTHETIC-01.json"
    }
  ]
}
```

This one-row illustration is a schema example, not a runnable study cohort. The original study uses stable `GBM-01`–`GBM-15` and `MEN-01`–`MEN-05` labels. Clinical accessions must not be substituted into public figures or files. Use a unique pseudonymous `patient_group` for each independent patient. BAM extraction manifests are JSON lists of records with `path`, `bytes` and `mtime_ns`. Paths and source-integrity values are local to the authorised data holder. No original cohort manifest is distributed.

Some follow-up stage provenance also records `patient_linkage_basis`; if needed, supply a non-identifying explanation of the verified linkage. Do not include accession numbers, full clinical records or patient identifiers.

## CpG counts and regional profiles

- `data/reference/cpg_positions.npy`: autosomal reference-CpG C positions, concatenated by chromosome, `uint32`.
- `data/reference/cpg_offsets.npy`: 23 chromosome offsets, `uint32`.
- `data/reference/chromosome_lengths.npy`: 22 autosomal lengths, `uint32`.
- `data/counts/{study_label}/site_counts.npz`: reference-CpG rows with eight columns, `uint32`: passing C, 5mC, 5hmC; failed 0.80 calls; explicit paired calls; reference-CpG opportunities; intact read-CpG opportunities; missing explicit calls at intact read CpGs. The array is stored under the `counts` key, with field names under `fields`.
- `data/regional/common_cpg_mask_depth5.npy` and `...depth10.npy`: identical-CpG eligibility masks across the cohort.
- `data/regional/cache/features5.json` and `features10.json`: locus definitions, feature IDs and common-CpG counts.
- `data/regional/pooled_profiles.npz`: specimen labels and diagnoses; global, window and feature modification arrays for both coverage thresholds. The last dimension is `[5mC, 5hmC]`. `valid5` and `valid10` specify eligible windows. Corresponding feature order comes from the feature JSON files.

Denominators include only passing C, 5mC and 5hmC calls. Both explicit modifications are required; absent tags are not silently interpreted as unmethylated C. ML bin midpoints are `(v + 0.5) / 256`, canonical probability is `max(0, 1 - p5mC - p5hmC)`, and the most likely state must reach 0.80. Both read ends are trimmed by 100 query bases. The native counter and Python decoder are tested together using synthetic forward and reverse reads.

## Panels and single-molecule data

`data/plans/molecular/`, `followup/` and `strengthening/` hold frozen panel definitions, matching pairs, fold assignments, distance pairs, protocol choices and the full feature universe. Feature IDs and locus IDs are distinct: retain the mappings in the panel files. Seed construction includes stable specimen and locus labels; changing these labels can change Monte Carlo draws.

Molecule `.npz` files live under `data/molecules/molecular/`, `followup/` and `strengthening/`. They contain:

- `states`: reads × reference CpGs, signed integer values `-1` (missing/failed), `0` (C), `1` (5mC), `2` (5hmC).
- `cpg_positions`: reference CpG positions in the same column order.
- `common_mask`: the fixed common-CpG eligibility mask for those columns.
- Some internal archives additionally contain `read_ids` and QC metadata. The numerical association kernels operate on the state matrices; raw read identifiers are not needed for public interpretation.

Preserve missingness and molecule membership. bedMethyl summaries lose molecule membership and therefore cannot reproduce these analyses. Any patient matrix release remains subject to the approved access arrangement.

Some original stage interfaces locate per-specimen matrices beneath their `results` alias. To use canonical matrix locations, create private per-specimen links such as `.analysis/focused_followup/results/{study_label}` → `data/molecules/followup/{study_label}`; use the same convention for `story_analysis`/`molecular`. `init-workspace` creates the parent layout; specimen-specific links must be populated with the supplied data package. Never overwrite a results directory containing unique files while creating links.

## Completed results for figures

- `results/regional/tables/`: CpG coverage, regional contrasts, feature tables and cohort summaries.
- `results/molecular/tables/`: original panel and molecule/block summaries.
- `results/followup/tables/`: held-out evaluation, distance/abundance matching, constrained null and distance curves.
- `results/strengthening/tables/`: influence, shared-encoding scores and 500-allocation randomisation summaries.
- `supplementary/results/molecular/`: contextual and allele/depth sensitivity summaries.
- Panel/feature definitions and selected matrices remain necessary for the same-molecule example panels.

The figure code lists the exact table paths it reads. It intentionally fails when required results or eligibility checks are missing; it does not substitute invented values. The combined supplementary PDF contains nine figures with legends; the manuscript text itself is not distributed.
