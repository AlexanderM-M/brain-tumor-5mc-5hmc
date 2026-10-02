# AGE manuscript reproducibility

Code for the three-figure AGE manuscript: fixed Horvath projections, molecule/patient resampling, qualified WGS-purity associations, five-CpG reference substitution, modification-resolved chemistry, local four-CpG discordance, beta-matched controls, targeted external validation and exploratory GBM survival.

This is a processed-input release. No raw patient sequencing, clinical records, accession maps or exact clinical dates are included. Restricted data require appropriate ethics approval and an executed Data Transfer Agreement (DTA) with the Medical University of Innsbruck. Figure source tables alone do not contain the molecule/count matrices needed to repeat the resampling analyses; those analyses require the appropriate de-identified processed inputs described below. Missing inputs are not reconstructed or invented.

## Layout and environment

- `workflow/AGE_ANALYSIS/`: portable copies of retained scientific scripts, in their dependency layout.
- `run.py`: installs these scripts in a separate workspace and invokes one named stage.
- `wgs_qc.py`: validates existing SAVANA processed CNA outputs and completion evidence.
- `survival_core.py`, `survival.py`: retained Cox mathematics and a duration-only command-line interface.
- `format_tables.py`: presentation-only TSV-to-XLSX export for the supplementary tables.
- `source_manifest.json`: source-relative provenance and hashes; no patient identifiers or machine paths.
- `tests/`: synthetic numerical and interface tests; no patient data.

Recorded scientific environment: Python 3.9.25, NumPy 1.26.4, pandas 2.3.3, SciPy 1.13.1, pysam 0.22.1 and Matplotlib 3.5.1. SAVANA source version was 1.3.8. Unknown historical basecaller, modification-model or upstream versions are not supplied by inference. `requirements.txt` pins the recorded Python scientific packages; openpyxl is an additional presentation dependency. Linux is required by the retained file-locking routines.

```bash
python -m pip install -r requirements.txt
python age/run.py --workspace /path/to/processed-workspace init
```

`init` only copies code and creates output directories. It refuses to overwrite modified workspace scripts. Put processed inputs in that separate workspace; never use the original/raw-data directory as the workspace. The scientific stages write derived tables and caches within it. Input checks fail if required data are absent; the launcher does not download private data.

## Controlled-workspace input contracts (not the public table schema)

All specimen keys must be the fixed pseudonyms GBM-01 through GBM-15 and MEN-01 through MEN-05. These are controlled-workspace contracts, not authorisation to distribute their inputs publicly. Public files are restricted to the privacy-minimized layer documented in release/README.md.

Paths below are relative to the workspace, not this repository.

| Component | Required processed inputs |
|---|---|
| Clock/chemistry | `data/reference/horvath2013/clock353_hg38.tsv` (353 rows; `CpGmarker`, `CoefficientTraining`, `reference_index`, hg38 coordinates); `clock353_sesame_reference.csv` (`Probe_ID`, `median`); `results/age_hox_joint/tables/{sample}_clock353_counts.npz` (`counts`, shape 353 × 3, unmodified/5mC/5hmC); `AGE_ANALYSIS/next_generation/00_audit/frozen_previous/01_tables/age_master_table.tsv` (de-identified fixed cohort, chronological age and previously verified clock/QC columns). |
| Molecule/patient bootstrap | `AGE_ANALYSIS/next_generation/intermediate/{sample}_molecules.npz`: `indices`, `row`, `col`, `state`, `read_hash`; hash labels are opaque molecule IDs, never original read names. Also the fixed triplet definitions and reference arrays required by `molecules.py`. The molecule calls must reproduce retained counts exactly. |
| WGS qualification | Pseudonymised technical-evidence table `AGE_ANALYSIS/deep_dive/04_WGS_genomic_backbone/tables/cohort_WGS_QC.tsv`; retained segment evidence in `FINAL_INTEGRATION/provenance/WGS_segment_evidence.tsv`; SAVANA fitted-purity/ploidy, ranked solutions, absolute/relative segments and preliminary allele-block stability summaries. `AGE_WGS_PROCESSED` is set through `--wgs-processed`; its per-sample folders contain `savana_cna/` processed files. No BAM/CRAM decoding occurs. Original full-decode/index/reference checks are evidence inputs, not reproducible from processed counts alone. |
| Purity associations/five CpGs | Final patient clock summary at `next_generation/final/tables/age_patient_summary.tsv`; qualified WGS QC; five-locus contribution baseline and annotation tables at `05_five_locus_signal/tables/` and `next_generation/tables/clock_locus_chemical_modes.tsv`. Frozen baseline tables are comparison inputs, not refitted clocks. |
| Local PDR/entropy | `deep_dive/.intermediate/local_epialleles/{sample}/chr{chrom}.npz` with `indices` and four-CpG 81-state `counts`; completion evidence; genomic reference CpG arrays, offsets and retained RefSeq/CpG-island annotations. The raw extractor is intentionally omitted. |
| Beta controls | Local metric caches above; existing genomic blacklist window mask; fixed five targets; same-patient coverage, quartet span, density and context arrays. Candidate/control caches can be reused if already frozen. `focused_beta.py` documents all matching gates, seeds and null construction. |
| External cohorts | Public processed metadata/target-row caches for GSE74193, GSE60274, GSE195684 and TCGA-GBM; five-locus local modification and canonical phenotype tables. `read_archive.py` accepts a verified archive with its download/checksum manifests, or `extract_targeted_rows.py` can prepare targeted caches from the public processed matrices. No new cohort or whole-methylome analysis is selected. |
| Figures | Frozen canonical five TSV tables under `FINAL_AGE_STORY/tables/`, plus `09_EXTERNAL_VALIDATION/tables/normal_brain_validation.tsv`, `GBM_validation.tsv`, and `five_locus_summary.tsv`. Figure code retains exactly three main figures. |
| Survival | Exactly five columns: `sample_id`, `diagnosis`, `DNAm_age_acceleration`, `post_sampling_survival_days`, `event`. Event is 1 for death, 0 for confirmed alive/right-censored. Supply the fixed primary cohort and, optionally, two separate duration-only sensitivity files. Clinical adjudication is performed privately upstream and is not inferred by this code. |

The existing table readers and assertions are executable schemas. Preserve canonical specimen order, feature order, count definitions, matching gates, seeds and the supplied missing-value conventions. Supplying rounded publication values instead of full-precision processed values can change estimates; use the full-precision source tables.

## Analysis order with authorised controlled inputs

These commands require the exact controlled inputs described above; they must not be run on privacy-transformed public ages. Preparation of this code release did not rerun patient analyses. The public marker makes the stage launcher refuse those inputs.

1. **Fixed clock and chemical decomposition:** run `clock-chemistry`, then `clock-bootstrap`. The clock is not trained or recalibrated. Molecule resampling and patient resampling retain their original distinction. The bootstrap writes `next_generation/tables/age_patient_summary.tsv`; copy that completed table into the corresponding `next_generation/final/tables/` location for downstream integration. Preserve the complete table, not just publication columns.
2. **WGS QC:** use `age/wgs_qc.py --sample GBM-01 --folder PATH --log PATH --exit-code 0 --output PATH` on de-identified processed SAVANA output. Supply the retained cohort technical/allele-block evidence, then run `wgs-fit-stability`, `purity-associations`, `five-cpg`, and `joint-ploidy-qc`. The last stage applies the recorded joint-ploidy qualification and repeats downstream integration with those qualified estimates. The original genotype/segmentation work is not rerun. The preliminary allele-derived range and rejected fits must remain explicitly qualified.
3. **Local four-CpG statistics:** run `local-annotation`, then `local-metrics GBM-01` for each specimen. Reuse the fixed blacklist mask and reference annotations. This computes PDR, entropy and rarefaction from processed counts.
4. **Beta-matched regions:** run `beta-setup` once in a fresh workspace (or retain its existing frozen caches), then `beta-patient --sample GBM-01` for each specimen, followed by `beta-chemistry` and `beta-summary`. These retain the original candidate sampling, beta calipers, patient eligibility and dependence sensitivities.
5. **Targeted external validation:** use `external-ingest --archive PATH` or `external-targets SOURCE --cache PATH --name STEM --clock353`, then `external-validation`. Copy its completed `normal_brain_validation.tsv`, `GBM_validation.tsv` and `five_locus_summary.tsv` from the external cache to `AGE_ANALYSIS/deep_dive/09_EXTERNAL_VALIDATION/tables/` before figure generation. Use the cache schemas read by `curate_inputs.py`; do not change sample eligibility or infer missing ages.
6. **Survival:** invoke the duration-only interface below. No age/recurrence adjustment or secondary molecular models are added to this sparse-event analysis.
7. **Figures/tables:** run `main-figures` from the frozen canonical/external tables. This generates the retained three PDFs and preview PNGs without refitting statistical models. `format_tables.py` converts publication-ready, de-identified TSVs to two professionally formatted workbook sheets without recalculating results.

Example stage invocation:

```bash
python age/run.py --workspace /path/to/processed-workspace clock-chemistry
python age/run.py --workspace /path/to/processed-workspace local-metrics GBM-01
python age/run.py --workspace /path/to/processed-workspace beta-patient --sample GBM-01
python age/survival.py --input primary_durations.tsv \
  --sensitivity-earlier sensitivity_earlier.tsv \
  --sensitivity-later sensitivity_later.tsv --output /path/to/survival-output
python age/format_tables.py --cohort cohort.tsv --survival survival.tsv \
  --output /path/to/Supplementary_Tables.xlsx
```

The primary survival predictor is continuous apparent acceleration per +10 years. Time is post-sampling survival, not survival from diagnosis. Standard Cox, leave-one-out and bootstrap rules are preserved; bootstrap CIs are withheld when more than 10% of fits fail. No optimal cutpoint or machine-learning step is present. A finite fit with three deaths does not establish prognostic independence.

## Validation and scope

Scientific kernels were compared with retained source functions. Portability changes are limited to configurable paths, CPU-affinity handling, workspace bookkeeping, and replacing private clinical integration with a duration-only interface. Release validation uses syntax/import checks and synthetic numerical examples; no new patient analysis was run. The repository also retains its earlier regional/molecular workflow under `src/brain_5mc_5hmc`; those modules support the earlier study and are not additional AGE main figures.

```bash
python -m unittest discover -s age/tests -v
```

The complete original sequencing workflow and private clinical-date adjudication cannot be reproduced from this public repository alone. Their exclusion is deliberate. Missing upstream versions and undistributed processed matrices remain explicit limitations, not implied public access commitments.

## Frozen manuscript release

See [release instructions](../release/README.md) for v1.0.1, the privacy-transformed public source tables and their restricted-input limitations. Exact original figure/table regeneration is not supported by these public tables. Reserved release DOI: 10.5281/zenodo.23104961; publication awaits institutional data-licence confirmation.
