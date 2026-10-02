# Brain tumour methylation analyses

Reproducibility code for regional 5mC/5hmC analyses and the AGE manuscript, **Age-concordant and discordant methylation remodelling accompany apparent epigenetic age acceleration in brain tumours**.

The AGE study includes 15 glioblastomas and five meningiomas. The code covers fixed Horvath clock reconstruction, molecule/patient resampling, WGS tumour-content qualification and associations, five-CpG reference substitution, modification-resolved chemistry, local four-CpG PDR/entropy, beta-matched controls, targeted external validation, exploratory GBM survival and the three retained manuscript figures.

Start with [AGE instructions and processed-input contracts](age/README.md). Install the recorded scientific packages with `python -m pip install -r requirements.txt`; the recorded interpreter is Python 3.9.25. Work in a separate processed-data workspace. No patient data analysis runs during installation.

- `age/`: AGE scripts, portable launcher, input documentation and synthetic tests.
- `src/brain_5mc_5hmc/`: earlier regional and single-molecule workflow.
- `docs/`: earlier workflow, input and validation documentation.
- `tests/`: synthetic tests for the earlier workflow.

Public source tables use study pseudonyms and privacy-minimized clinical variables. Exact chronological ages derived from clinical dates, exact clinical dates, clinical accession mappings, patient-level survival times and vital status, raw patient sequencing data and detailed patient-level matrices are not included in the public release. Access to restricted patient-level data requires appropriate ethics approval and an executed Data Transfer Agreement (DTA) with the Medical University of Innsbruck.

For the earlier workflow, install with `python -m pip install -e '.[test]'`, copy `config.example.json` to `config.local.json`, and see [workflow documentation](docs/workflow.md). Both test suites use synthetic inputs.

See [CITATION.cff](CITATION.cff) and the [MIT software licence](LICENSE). The reserved DOI for frozen public release v1.0.1 is [10.5281/zenodo.23104961](https://doi.org/10.5281/zenodo.23104961). Zenodo publication is pending institutional confirmation of the source-data licence; no data licence is assumed.

## Frozen manuscript release

See [release instructions](release/README.md) for v1.0.1, the privacy-transformed public tables and controlled-input limitations. These public tables do not reproduce exact clinical-age analyses or the complete original figures. Use the dedicated public-output command; do not substitute rounded ages into the scientific workflows.
