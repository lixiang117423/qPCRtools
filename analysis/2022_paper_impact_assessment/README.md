# qPCRtools 2022 paper impact assessment

Compared qPCRtools 1.0.1 with the corrected 1.0.2 implementation on the
example datasets distributed with the 2022 Frontiers in Genetics paper.

## Overall finding

The numerical results feeding Figures 2–4 are affected. The paper’s
conclusions must be re-run from the corrected package before those figures
or source-data sheets are reused.

This assessment covers the reproducibility impact on the 2022 article's
example analyses. It does not estimate how many external papers may have
used the affected 1.0.1 calculations; that would require a separate
bibliographic and author-contact audit.

## Publication action

- Figure 1 / CalCurve can be retained because all numerical outputs are
  unchanged.
- Figures 2–4 and their source-data sheets must be regenerated with
  qPCRtools 1.0.2; do not reuse the 1.0.1 figure panels.
- Re-run any conclusions or supplementary statistics that depended on the
  affected group means, standard errors, significance calls or geNorm
  input.
- Use the numeric tables in `results/` as the traceable comparison record
  for a correction or repository release note.

## Output row counts

| Case | 1.0.1 rows | 1.0.2 rows |
| --- | ---: | ---: |
| calcurve | 4 | 4 |
| calexpcurve | 80 | 40 |
| calexp2ddct | 40 | 40 |
| calexprqpcr | 32 | 32 |

## Figure 1 / CalCurve

All numeric standard-curve entries are unchanged (maximum absolute difference `0.00e+00`). The figure code changed only to remove deprecated `..rr.label..` syntax.

## Figure 2 / CalExpCurve

The corrected version back-transforms log concentration to concentration
and averages technical wells within each biological replicate. The result
table falls from 80 to 40 rows; plotted group values and significance
calls must be regenerated.

| Gene / treatment | 1.0.1 mean | 1.0.2 mean | Change | Significance |
| --- | ---: | ---: | ---: | --- |
| OSPOX8 / CK | 1.036 | 1.488 | +43.7% | — → — |
| OSPOX8 / Treatment | 0.9551 | 0.5327 | -44.2% | NS → NS |
| OsCeBip / CK | 0.8421 | 0.5465 | -35.1% | — → — |
| OsCeBip / Treatment | 0.8478 | 0.497 | -41.4% | NS → NS |
| OsPR10 / CK | 1.249 | 34.1 | +2630.0% | — → — |
| OsPR10 / Treatment | 1.056 | 2.915 | +176.2% | NS → NS |
| OsRBBI2 / CK | 0.5399 | 0.01111 | -97.9% | — → — |
| OsRBBI2 / Treatment | 0.5294 | 0.007104 | -98.7% | NS → NS |
| OsWAK91 / CK | 0.7632 | 0.2112 | -72.3% | — → — |
| OsWAK91 / Treatment | 0.9844 | 0.9893 | +0.5% | * → NS |

## Figure 3 / CalExp2ddCt

The corrected version calibrates ΔCt within the control group and reports
SE = SD/sqrt(n). Group means and error bars change; the bundled example
keeps the same significance-symbol pattern.

| Gene / treatment | 1.0.1 mean | 1.0.2 mean | Change | Significance |
| --- | ---: | ---: | ---: | --- |
| OSPOX8 / CK | 1.224 | 1.014 | -17.2% | — → — |
| OSPOX8 / Treatment | 0.8381 | 0.6941 | -17.2% | * → * |
| OsCeBip / CK | 1.189 | 1.179 | -0.8% | — → — |
| OsCeBip / Treatment | 1.113 | 1.104 | -0.8% | NS → NS |
| OsPR10 / CK | 2.033 | 1.316 | -35.2% | — → — |
| OsPR10 / Treatment | 0.6896 | 0.4466 | -35.2% | NS → NS |
| OsRBBI2 / CK | 1.162 | 1.009 | -13.1% | — → — |
| OsRBBI2 / Treatment | 0.8939 | 0.7765 | -13.1% | NS → NS |
| OsWAK91 / CK | 0.7399 | 1.212 | +63.8% | — → — |
| OsWAK91 / Treatment | 1.683 | 2.757 | +63.8% | * → * |

## Figure 4 / CalExpRqPCR

The corrected version no longer merges equal numeric expression values and
uses relative quantities for geNorm. Group means and standard errors change;
the published figure must be regenerated with 1.0.2.

| Gene / group | 1.0.1 mean | 1.0.2 mean | Change |
| --- | ---: | ---: | ---: |
| OSPOX8 / CK | 1.394 | 1.265 | -9.3% |
| OSPOX8 / Treatment | 1 | 1 | +0.0% |
| OsCeBip / CK | 1.124 | 1 | -11.0% |
| OsCeBip / Treatment | 1 | 1.034 | +3.4% |
| OsPR10 / CK | 1.947 | 1.947 | +0.0% |
| OsPR10 / Treatment | 1 | 1 | +0.0% |
| OsWAK91 / CK | 1 | 1 | +0.0% |
| OsWAK91 / Treatment | 4.22 | 3.745 | -11.2% |

## Reproduction

Run `01_run_version_case.R` separately with R libraries containing
qPCRtools 1.0.1 and 1.0.2, then run `02_compare_versions.py`.

The complete numeric tables are retained in `results/`.
