# Frozen public release v1.0.1

Reserved Zenodo DOI for this frozen public release: **10.5281/zenodo.23104961** (https://doi.org/10.5281/zenodo.23104961). This DOI was supplied by the depositor; reservation is not evidence of publication. Final Zenodo publication is blocked until the institution confirms the patient-derived source-data licence. Software retains the MIT licence; no source-data licence is assumed or granted. See `DATA_TERMS.txt` and `zenodo_metadata.json`.

Public source tables use study pseudonyms and privacy-minimized clinical variables. Exact chronological ages derived from clinical dates, exact clinical dates, clinical accession mappings, patient-level survival times and vital status, raw patient sequencing data and detailed patient-level matrices are not included in the public release. Access to restricted patient-level data requires appropriate ethics approval and an executed Data Transfer Agreement (DTA) with the Medical University of Innsbruck.

## What is in the public data layer

Ten publication TSVs are supplied in `processed/`, including two privacy-transformed local-cohort tables and an aggregate-only survival table. Each adult age is the completed integer number of years; a minor is represented only as `<18`. The public age headers deliberately differ from the exact-input schemas. Do not recalculate any analysis using these transformed ages.

`patient_phenotype_WGS_chemistry.tsv` supplies the public Figure 1-related molecular values, but it is NOT the exact Figure 1 source table. Individual apparent acceleration, its clinical-age-offset CIs, post-reference-substitution offsets and age-offset significance flags were removed. Together with predicted DNAm ages, those fields could reveal exact chronological ages by subtraction. The public cohort table likewise omits apparent acceleration, follow-up days, vital status and survival-inclusion flags. Molecular DNAm-age predictions are model outputs, not chronological ages, and remain unchanged.

The previous patient-linked Supplementary Table 2 is absent. `Supplementary_Table_2_survival_aggregate.tsv` contains only the three already-reported primary/sensitivity Cox summaries: HR, CI, P value, sample size, event count and descriptive model labels. Named patient-omission rows were removed because differences in event counts could reveal an individual's death status. No survival model was rerun.

The other seven aggregate/external TSVs are unchanged. Public GEO/TCGA accessions are retained. Study pseudonyms remain GBM-01 through GBM-15 and MEN-01 through MEN-05. `input_manifest.json` records checksums of the public inputs. The marker `processed/PUBLIC_INPUTS.json` prevents the original stage launcher from silently accepting a privacy-transformed workspace.

## Supported public outputs

Use Python 3.9.25 and the scientific dependencies in `../requirements.txt`; `../age/recorded_environment.json` preserves the recorded versions. From the repository root:

```bash
python -m pip install -r requirements.txt
python release/public_source_data.py --output ../public-release-outputs
python -m unittest discover -s age/tests -v
```

The dedicated public command validates the public schemas and emits a privacy-labelled workbook containing the two local summary tables and aggregate survival results, plus a molecular-signal overview figure. It only displays supplied values; it does not fit or recalculate any scientific estimate and does not produce replacement manuscript results. Adult integer ages and `<18` remain text categories in the workbook.

## Exact manuscript reproduction requires controlled inputs

The original analysis, figure and source-generation scripts are retained unchanged under `age/workflow/`, along with the original supplementary-table formatter and survival implementation. They cover DNAm-age reconstruction, molecule/patient resampling, WGS QC and associations, five-CpG reference substitution, matched-set analyses, external cohorts, chemical decomposition, local PDR/entropy and figures. Their schemas require exact authorised inputs documented in `../age/README.md`; those inputs are not all public. All original scientific numerical results are unchanged.

The privacy-transformed tables cannot regenerate exact-age analyses, patient-level survival analyses, the complete original Figure 1 or the complete original Figure 2. Full exact figure/table regeneration claims from v1.0.0 do not apply to this release. The public source-data exporter is separate from the exact-input analytical workflow. No exact clinical values were inferred or recomputed for this release.

## History and deposition

v1.0.1 is a new commit/tag; v1.0.0 is not overwritten. Earlier Git commits contain now-restricted patient information. Deleting a release ZIP alone does not remove Git history. Repository-owner/institutional review and separately authorised remediation are needed; no history rewrite, remote deletion or force-push is performed by this release preparation.

Upload only the single v1.0.1 ZIP and verify its supplied SHA-256. Do not upload the previous ZIP, repository history, private audits or internal inputs. Publication must await institutional data-licence confirmation and review of the historical exposure. No licence for patient-derived data has been invented.
