# Frozen manuscript release v1.0.0

This public release contains the analysis code and ten de-identified publication source TSV tables. The detailed patient-level matrices remain controlled-access under the approved sharing scope. No raw sequencing, accession mapping or exact clinical date is included. The release is sufficient for regeneration of the three main figures and the two supplementary table sheets from frozen estimates. It is not a complete public rerun of the statistical analyses.

## Environment and figure/table regeneration

Use Python 3.9.25 on Linux and install `requirements.txt`. The scientific versions are recorded in `age/recorded_environment.json`. Openpyxl is an additional table-presentation dependency. From the extracted release directory:

```bash
python -m pip install -r requirements.txt
python age/run.py --workspace ../age-reproduction init
cp -R release/processed/. ../age-reproduction/
python age/run.py --workspace ../age-reproduction main-figures
python age/format_tables.py --cohort ../age-reproduction/AGE_ANALYSIS/deep_dive/FINAL_AGE_STORY/supplementary/Supplementary_Table_1_cohort.tsv --survival ../age-reproduction/AGE_ANALYSIS/deep_dive/FINAL_AGE_STORY/supplementary/Supplementary_Table_2_survival.tsv --output ../age-reproduction/Supplementary_Tables.xlsx
```

Generated PDFs are under `../age-reproduction/AGE_ANALYSIS/deep_dive/FINAL_AGE_STORY/figures/`. Input checksums are in `release/input_manifest.json`. The source TSVs preserve the reported values without refitting models; figure generation produces a derived external-cohort source table.

## Analysis coverage and restricted dependencies

| Analysis | Included implementation | Public rerun status |
|---|---|---|
| DNAm-age reconstruction and 5mC/5hmC decomposition | `age/workflow/AGE_ANALYSIS/next_generation/scripts/chemistry.py` | Requires restricted U/M/H count matrices and additional reference inputs |
| Molecule/patient resampling | `age/workflow/AGE_ANALYSIS/next_generation/scripts/molecules.py` | Requires restricted molecule matrices and linked input contracts |
| Genomic tumour-content QC and associations | `age/wgs_qc.py`, `age/workflow/AGE_ANALYSIS/deep_dive/FINAL_INTEGRATION/scripts/post_wgs_analysis.py` | Requires restricted processed SAVANA/technical evidence |
| Five-CpG reference substitution and matched sets | `post_wgs_analysis.py`, `06_independent_preparation/scripts/five_null.py` | Frozen matched-null source table included; upstream reconstruction requires restricted matrices |
| Beta-matched local analyses and PDR/entropy | `07_FIVE_LOCUS_BETA_VALIDATION/scripts/`, `06_PDR_EPY/scripts/` under the workflow | Requires restricted local epiallele matrices and genomic annotation/control inputs |
| External cohorts | `09_EXTERNAL_VALIDATION/scripts/` under the workflow | Result tables included; public target-row caches/metadata must be obtained separately using the documented cohort accessions and input contracts |
| Figure and source-data generation | `FINAL_AGE_STORY/scripts/revise_figures.py` under the workflow | Included publication source tables are sufficient; verified in an isolated workspace |
| Supplementary table generation | `age/format_tables.py` | Both source TSVs included |
| Exploratory survival | `age/survival.py`, `age/survival_core.py` | Frozen publication tables included; precise adjudicated duration inputs remain controlled-access |

See `age/README.md` for detailed analysis order, schemas, seeds and limitations. The archive contains scripts for all requested analysis categories but does not contain all inputs required to rerun the underlying analyses. Historical Dorado/model versions are not inferred. Synthetic tests check computational kernels; they are not patient-analysis validation.

## Restricted access

Raw patient-level sequencing data, clinical accession mappings and exact clinical dates are not publicly available because they may contain potentially identifying information and are subject to ethical and institutional restrictions. Access to these restricted data requires appropriate ethics approval and an executed Data Transfer Agreement (DTA) with the Medical University of Innsbruck. Detailed patient-level clock/molecule and local epiallele count matrices also remain controlled-access, with the same access requirements. Public supplementary/source tables contain fixed study pseudonyms only.

## Archival status

This is the v1.0.0 manuscript release. The archive filename and accompanying release record identify the exact Git commit. No Zenodo DOI has been obtained; deposition requires an authenticated Zenodo account. Do not cite this release as permanently archived until deposition is published and its DOI has been verified. Software retains its existing MIT licence; confirm the data licence in Zenodo metadata before publishing. The software licence does not grant rights to restricted data.
