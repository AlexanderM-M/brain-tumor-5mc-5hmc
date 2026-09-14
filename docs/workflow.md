# Workflow

Run each stage as a separate command:

```bash
brain-5mc-5hmc --config config.local.json run STAGE
```

The command names preserve the original stage identifiers so that archived plans and outputs can be traced. They are functions within eight modules, not separate scripts. `docs/source_map.json` records each stage, its function prefix, dependencies and the SHA-256 of its original source. Dependencies initialize helper functions and constants; they do **not** automatically run preceding analysis stages.

## Entry points

Choose the entry point that matches the authorised inputs:

1. **Raw mapped reads:** indexed BAMs with explicit MM/ML joint 5mC/5hmC calls, sample manifests and the reference are required. Counting and target extraction can be expensive. The public repository contains no patient reads.
2. **Processed analysis inputs:** common-CpG counts/profiles, the frozen annotation snapshots, selected-panel plans and read-linked state matrices can support regional and molecular reanalysis. Some historical stages also perform BAM extraction or count validation; see the matrix-only entry point below.
3. **Completed analysis outputs:** retained tables, plans and the example state matrices support figure regeneration without full extraction or permutation scoring.

The implementation retains the study's cohort assumptions: 20 independent specimens, 15 glioblastomas and five meningiomas. Adapting it to a different cohort requires reviewing cohort dimensions, selection rules and eligibility thresholds. It is research reproduction code, not a general clinical pipeline.

## Stage order

| Step | Stages, in execution order | Main inputs and products |
| --- | --- | --- |
| 1. Regional profiles | `deep_regional --reference`; `deep_regional --run --workers 4`; `pooled_analysis --workers 4`; `interpret_pooled` | Reference, BAM manifests → CpG counts, five-/ten-call profiles, regional tables |
| 2. Original promoter panel | `story_panel`; `story_extract --workers 4`; `story_molecules` | Regional results and reference → candidate/comparison panel, targeted reads, state matrices, three-CpG scores |
| 3. Held-out evaluation | `followup_select`; `followup_extract`; `followup_molecules`; `followup_q1` | Fixed coverage universe → held-out and blind panels, molecular scores |
| 4. Abundance and distance | `followup_distance`; `followup_distance_summary`; `followup_row_null`; `followup_row_summary`; `followup_distance_ranges` | Panels and molecule evidence → matching and row/column-preserving-null results |
| 5. Promoter influence | `strengthen_influence` | Held-out residuals and promoter pairs → pair/locus-deletion checks |
| 6. Three encodings | `strengthen_score three`; `strengthen_compare` | Existing held-out matrices → shared-pair 5hmC, 5mC and combined scores |
| 7. Whole-procedure randomisation | `strengthen_plan`; `strengthen_extract`; `strengthen_matrices`; `strengthen_score h`; `strengthen_permutation` | Full promoter universe and 500 allocations → reselection, matrix scoring, max-statistic results |
| 8. Contextual follow-up | `story_copy_number`; `story_alleles`; `story_context`; `story_repeat_sensitivity` | Original panel, BAMs/reference/annotations and block scores → descriptive context and sensitivity tables |
| 9. Publication figures | `manuscript_results_figures`; `manuscript_supplementary_figures` | Completed tables and example matrices → main figures and supplementary figures with legends |

Steps 5, 6 and 8 do not need to run sequentially relative to one another once their inputs exist. The three-encoding step must precede its comparison renderer. The full randomisation step requires reselection across the entire eligible promoter universe, not only previously selected candidates.

## Processed matrix entry point

For pre-existing matrices, `followup_row_null`, `strengthen_score`, influence summaries and randomisation summaries operate on processed inputs. `followup_molecules` still extracts states from follow-up BAMs and verifies counts; its kernels `blocks_patterns_metrics` and `followup_blocks_joint_excess` in `molecules.py` can be used directly with complete three-CpG matrices. `followup_distance` includes read-based MAPQ sensitivity and therefore also requires targeted BAMs. The source is included to document those analyses; public code alone cannot recreate unavailable reads.

Prepare every required matrix, plan and consolidated result before selecting a later entry point. Archived completion flags are not substitutes for actual inputs. Original execution/visual-review gates were removed from the release; numerical QC and source-integrity checks remain.

## Key method functions

- `extraction.counts_Engine`: streams BAMs through the native counting kernel.
- `molecules.blocks_states_at`: decodes explicit paired modification calls into C/5mC/5hmC states.
- `regional.regional_effects` and `regional.regional_bootstrap_intervals`: specimen-level contrasts and pointwise intervals.
- `promoters.selection_select`: training-only selection and comparison matching for a held-out fold.
- `association.row_null_trade`: Python binary row/column-preserving randomisation.
- `association.native_binary_score`: compiled constrained-null scoring for three encodings.
- `robustness.randomisation_patient` and `robustness.randomisation_statistic`: equal-pair specimen scores and the supported patient-group statistic.

Before calling stage functions directly, configure the workspace and initialize the corresponding stage with `runtime.initialize("stage_name")`. Use a fresh Python process for another workspace or workflow stage: some figure helpers deliberately share the current output directory. Worker pools explicitly use Linux fork to inherit the initialized stage state. Windows/spawn-based execution is not supported by this release.

## Reference annotations and outputs

Archived hg38 RefSeq/CpG-island/RepeatMasker snapshots are required for exact reproduction. Annotation downloads are disabled by default. Set `allow_annotation_download: true` in the local JSON configuration only when intentionally creating a new snapshot; results from a new snapshot are not automatically identical to the study.

The private workspace uses canonical `data/`, `results/`, `figures/`, `supplementary/`, `paper/` and `provenance/` directories. Relative aliases under `.analysis/` preserve the proven input/output interfaces. `init-workspace` creates empty directories and aliases without copying patient data. See `src/brain_5mc_5hmc/layout.json` for the mapping. Do not run output-producing stages directly against the only copy of frozen study results.
