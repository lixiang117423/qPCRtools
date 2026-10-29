# CRAN submission comments

This submission updates qPCRtools from 1.0.1 to 1.0.2.

## Test environments

- macOS 27.0.1, R 4.6.1, aarch64-apple-darwin27.0.0

## R CMD check results

- `R CMD check --as-cran --no-manual`: 0 errors, 0 warnings, 0 notes

## Changes

- Corrects the reference-group calibrator in `CalExp2ddCt()`.
- Corrects the standard error to SD/sqrt(n).
- Back-transforms `CalExpCurve()` output to concentration and aggregates
  technical replicates within biological replicates.
- Corrects RqPCR replicate aggregation and geNorm stability input.
- Corrects Tukey compact-letter grouping and removes deprecated plot syntax.
- Removes artificial red markers from boxplots and fixes stacked bar heights.
- Adds `RawPValue`, `PValue`, and `PAdjustMethod` columns to the expression
  result tables; ANOVA/Tukey output now reports actual adjusted P values.

## Downstream checks

The release includes regression tests for the affected functions. A separate
impact assessment documents the expected changes to the 2022 paper's example
figures and will accompany the GitHub release.
