# Brain tumour 5mC and 5hmC

Analysis code for **What the sum conceals: regional 5mC and 5hmC composition and single-molecule organisation in brain tumors**.

The study compares 15 glioblastomas and five meningiomas using nanopore-derived cytosine modification calls. It asks whether regional 5mC and 5hmC differences can be concealed by their combined fraction, and whether selected promoters show molecular organisation beyond modification abundance.

This repository contains **code, documentation and synthetic tests only**. Patient BAMs, modification counts, molecule matrices, clinical metadata and the manuscript are not included. Code availability does not grant access to patient data.

## Code organisation

Twenty-eight analysis and figure scripts have been consolidated into eight readable modules. Stage prefixes distinguish functions with different scientific roles; `initialize_*` loads a stage's configuration and inputs, and `run_*` executes it. Importing the package does not start analyses. A source map preserves the relationship to the original analysis.

| Module | Responsibility |
| --- | --- |
| `extraction.py` | CpG counting, targeted BAM extraction and molecule matrices |
| `regional.py` | Common-CpG profiles, regional contrasts and candidate summaries |
| `promoters.py` | Candidate selection, comparison matching, held-out folds and allocation plans |
| `molecules.py` | Three-CpG molecular patterns and held-out evaluation |
| `association.py` | Distance/abundance matching, constrained nulls and three modification encodings |
| `robustness.py` | Promoter deletion, encoding comparisons and whole-procedure randomisation |
| `context.py` | Genomic context, repeat sensitivity, depth and allele-linked follow-up |
| `figures.py` | Four main figures and a nine-page supplementary figure PDF with legends |

Two small C kernels accelerate counting and constrained randomisation. They are provided as source and compiled locally. `runtime.py` and `cli.py` handle configuration and commands.

## Installation and synthetic checks

Linux, Python 3.9 or later, a C compiler and samtools are required for the complete workflow. BAM extraction uses samtools; plotting from processed results does not. Use an isolated environment:

```bash
git clone https://github.com/AlexanderM-M/brain-tumor-5mc-5hmc.git
cd brain-tumor-5mc-5hmc
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
pytest -q
brain-5mc-5hmc list
```

The tests generate artificial reads and matrices in a temporary directory. They do not need patient data or an hg38 download. They check strand-aware calls, the 0.80 threshold, duplicate handling, equal-CpG weighting, held-out selection, missingness-preserving randomisation and patient-level inference.

## Configure a private analysis workspace

```bash
cp config.example.json config.local.json
# Edit workspace, reference and samtools paths in config.local.json.
brain-5mc-5hmc --config config.local.json init-workspace
brain-5mc-5hmc --config config.local.json build-native
```

Paths are relative to the configuration file unless absolute. Keep this workspace outside the Git repository. Configuration is local and ignored by Git. The `workers` setting controls fixed-size worker pools; extraction and regional stages also accept their documented `--workers` argument. The hg38 FASTA must be indexed, with autosomes named `chr1` through `chr22`.

Populate the workspace with authorised inputs before running analysis stages. **An empty workspace is not sufficient to reproduce the patient results.** See [input formats](docs/inputs.md), [workflow and execution order](docs/workflow.md), and [validation and limitations](docs/validation.md).

To redraw manuscript figures from the retained analysis tables and matrices:

```bash
brain-5mc-5hmc --config config.local.json run manuscript_results_figures
brain-5mc-5hmc --config config.local.json run manuscript_supplementary_figures
```

Outputs are written to `figures/main/manuscript/`, `supplementary/figures/manuscript/`, and the combined supplementary PDF in `supplementary/`. Run these commands in a working copy of the authorised data layout: output stages can replace existing derived files.

## Scientific scope

Patients provide biological replication. Common CpGs receive equal weight within each specimen; specimens receive equal weight within diagnosis. Regional contrasts are glioblastoma minus meningioma in percentage points. Regional candidate thresholds are exploratory, without an FDR discovery claim. Small combined differences do not establish statistical equivalence or biochemical conversion.

Molecular scores use complete CpG pairs and preserve the stated eligibility rules. The constrained null preserves binary row totals, column totals and missingness. The 500-allocation test reselects and rematches promoters and evaluates the maximum statistic across 36 contiguous distance intervals; unsupported allocations remain at statistic zero. It tests the complete selection-and-scoring procedure. It is not independent-cohort validation or a test of a specific biological mechanism.

## Data access and citation

Processed study data may be requested from the corresponding authors, subject to ethical approval and institutional requirements. The underlying patient BAM files are not offered for sharing. Some stages require read-linked matrices; bedMethyl files alone cannot reproduce single-molecule analyses. This repository does not establish journal acceptance of the data-access arrangement.

See [CITATION.cff](CITATION.cff). No manuscript DOI or archived software DOI has been assigned in this repository. For submission, archive a tagged release and cite its persistent identifier once available.

## Licence

The code is distributed under the MIT licence. The licence does not apply to patient data or third-party reference annotations.
