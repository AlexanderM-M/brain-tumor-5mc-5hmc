# Brain tumour 5mC and 5hmC

Analysis code for **What the sum conceals: regional 5mC and 5hmC composition and single-molecule organisation in brain tumors**.

Nanopore analysis of regional modification patterns and single-molecule organisation in 15 glioblastomas and five meningiomas. Eight modules cover CpG counting, promoter selection, molecular analysis, robustness checks and publication figures.

## Installation

Requires Linux, Python 3.9+, a C compiler and samtools for BAM extraction.

```bash
git clone https://github.com/AlexanderM-M/brain-tumor-5mc-5hmc.git
cd brain-tumor-5mc-5hmc
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
pytest -q
```

Tests use synthetic data; patient data are not included.

## Usage

Copy the example configuration and set the workspace, hg38 reference and samtools paths:

```bash
cp config.example.json config.local.json
# Edit config.local.json before continuing.
brain-5mc-5hmc --config config.local.json init-workspace
brain-5mc-5hmc --config config.local.json build-native
brain-5mc-5hmc list
brain-5mc-5hmc --config config.local.json run STAGE
```

Replace `STAGE` with a listed stage and supply its required inputs. Use a working copy of the data, as analysis stages write derived outputs.

See the [workflow](docs/workflow.md), [input formats](docs/inputs.md) and [validation](docs/validation.md) for details.

## Data access

Processed data may be requested from the corresponding authors, subject to ethical approval and institutional requirements. Patient BAM files are not shared.

## Citation and licence

See [CITATION.cff](CITATION.cff). Code is available under the [MIT licence](LICENSE).
