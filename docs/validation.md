# Release validation

The code was reorganised into eight modules while preserving the numerical kernels. Module-level state is initialized explicitly; stage-specific names prevent collisions. VM paths, temporary dependency paths, fixed worker-pool sizes and obsolete visual-review scheduling gates were replaced with portable configuration. Extraction no longer automatically starts the regional analysis; invoke the next documented stage explicitly. Annotation downloads require an explicit configuration setting.

## Synthetic tests

Twelve tests pass on the release environment. They cover:

- strand-aware promoter coordinates and inherited stage state in a fork worker;
- equally weighted CpG fractions despite unequal call depth;
- exploratory opposing-change criteria and repeatable bootstrap intervals;
- exclusion of held-out modification values from selection;
- sequential comparison matching without reuse;
- agreement of original and optimized three-CpG metrics;
- row/column/missingness invariance during Python null trades;
- native covariance agreement with direct complete-pair calculations in all three encodings;
- equal-pair specimen scoring with shared available distance bins;
- retention of unsupported label allocations at statistic zero;
- explicit forward/reverse MM/ML decoding, the 0.80 threshold, and native/scalar count agreement including duplicate ingestion.

These tests use artificial data only. The native tests build both C kernels locally. GitHub Actions is configured to run the tests; the presence of a workflow file is not a claim that a remote CI run has already passed.

## Packaging

The source passes static undefined-name/import checks. A wheel was built and installed into a separate temporary directory; all eight modules imported, the JSON resources were present, and both packaged C sources compiled successfully.

## Checks against the retained study

The local numerical regression audit, summarized in `regression_summary.json`, found:

- all 20 held-out promoter selections exactly matched both the original implementation and saved fold selections;
- all 26,585 window contrasts and single-omission direction flags matched;
- seeded 2,000-resample bootstrap intervals matched exactly for 80 checked windows;
- the rebuilt native null matched the original library exactly on synthetic matrices in all three encodings;
- full three-encoding scores matched on two retained patient/locus matrices;
- Python null moves, their random-number streams and block scores matched;
- recalculation from the saved 500 allocation maxima reproduced omnibus p = 0.004 and short-distance family-adjusted p = 0.034.

The four main figures and nine supplementary figures were regenerated into a separate temporary workspace. Pixel comparisons against the retained PNGs are recorded in `figure_regression.json`. The supplementary renderer produced nine pages, nine legends and nine figures with the original layout dimensions.

These checks establish agreement for the calculations and outputs tested. The full multi-terabyte BAM extraction, every molecular matrix and all 500 reselection/scoring runs were not repeated during this packaging task. The original extraction and scoring procedures remain available as code, with their numerical validation checks retained.

## Reproduction limits

The public checkout is sufficient for synthetic tests, but patient results require approved data access. Exact reproduction also requires the frozen reference/annotation inputs, panel identifiers, seeds and retained eligibility masks. A new annotation snapshot or changed specimen/locus identifiers can change outputs. Public code does not remove these input requirements.

The analysis is conditional on a cohort of 20 specimens and the fixed common-CpG universe. It does not establish independent-cohort validation, causality, expression changes or a clinical classifier. The strict shared-information comparison may be infeasible; the code retains that negative result rather than manufacturing a normalized comparison.
