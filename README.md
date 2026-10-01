# Brain tumour methylation analyses

Reproducibility code for regional 5mC/5hmC analyses and the AGE manuscript, **Age-concordant and discordant methylation remodelling accompany apparent epigenetic age acceleration in brain tumours**.

The AGE study includes 15 glioblastomas and five meningiomas. The code covers fixed Horvath clock reconstruction, molecule/patient resampling, WGS tumour-content qualification and associations, five-CpG reference substitution, modification-resolved chemistry, local four-CpG PDR/entropy, beta-matched controls, targeted external validation, exploratory GBM survival and the three retained manuscript figures.

Start with [AGE instructions and processed-input contracts](age/README.md). Install the recorded scientific packages with `python -m pip install -r requirements.txt`; the recorded interpreter is Python 3.9.25. Work in a separate processed-data workspace. No patient data analysis runs during installation.

- `age/`: AGE scripts, portable launcher, input documentation and synthetic tests.
- `src/brain_5mc_5hmc/`: earlier regional and single-molecule workflow.
- `docs/`: earlier workflow, input and validation documentation.
- `tests/`: synthetic tests for the earlier workflow.

De-identified processed data accompany the article as Supplementary Data and Source Data. Raw patient sequencing, clinical accession maps and exact clinical dates are excluded for ethical/privacy reasons. No raw patient data or credentials are included, and no controlled-access or data-on-request arrangement is promised. Full resampling reproduction requires the appropriate processed molecule/count inputs; figure-level source tables alone are insufficient.

For the earlier workflow, install with `python -m pip install -e '.[test]'`, copy `config.example.json` to `config.local.json`, and see [workflow documentation](docs/workflow.md). Both test suites use synthetic inputs.

See [CITATION.cff](CITATION.cff) and the [MIT licence](LICENSE). No archival DOI is asserted.
