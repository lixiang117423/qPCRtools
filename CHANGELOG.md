# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [2.0.0] - 2026-06-10

### BREAKING

- Rename all dot.case parameters to snake_case per development guide:
  `cq.table` → `cq_table`, `design.table` → `design_table`, `curve.table` → `curve_table`,
  `concen.table` → `concen_table`, `ref.gene` → `ref_gene`, `ref.group` → `ref_group`,
  `stat.method` → `stat_method`, `fig.type` → `fig_type`, `fig.ncol` → `fig_ncol`,
  `highest.concen` → `highest_concen`, `lowest.concen` → `lowest_concen`,
  `by.mean` → `by_mean`, `RNA.weight` → `rna_weight`, `remove.outliers` → `remove_outliers`

### Changed

- Rewrite README.md with correct parameter names and usage examples

## [1.1.0] - 2026-06-10

### Changed

- Conform to `DEVELOP_GUIDE.md` development standards
- Replace `exportPattern` in NAMESPACE with roxygen2-generated explicit exports
- Restructure roxygen2 documentation: `@title`, `@description`, `@param`, `@return`, `@importFrom`, `@export`, `@examples`, `@author`
- Wrap all `@examples` in `\dontrun{}`
- Fix `parse = T` to `parse = TRUE`
- Replace bare `sd()` with `stats::sd()`
- Clean up `globalVariables()` lists (remove function parameter names, keep only data column names)
- Extract internal helpers with `@keywords internal`: `cal_curve_by_mean`, `cal_curve_by_raw`, `find_outlier`, `cal_stat_test`, `build_exp_plot`, etc.
- Add `Encoding: UTF-8` and `VignetteBuilder: knitr` to DESCRIPTION
- Update `.Rbuildignore` and `.gitignore` with proper entries
- Vignette: use `read.table()` instead of `data.table::fread()`

### Fixed

- Fix wrong parameter names in `CalCurve` example (`dilu`/`by` → `dilution`/`by.mean`)
- Fix wrong parameter names in vignette for `CalCurve`
- Fix missing example data files in `inst/examples/`
- Remove stale `man/hello.Rd`

### Added

- Add input validation (`stop()`) to all exported functions
- Add proper package description in DESCRIPTION
- Add CHANGELOG.md

## [1.0.1] - 2022-11-11

### Fixed

- Fix Issue #2
- Add `remove.outliers` parameter to `CalExp2ddCt`

## [1.0.0] - 2022-06-20

### Added

- Initial CRAN release
- Functions: `CalRTable`, `CalCurve`, `CalExpCurve`, `CalExp2dCt`, `CalExp2ddCt`, `CalExpRqPCR`
