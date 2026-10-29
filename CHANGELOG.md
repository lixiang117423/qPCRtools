# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [2.0.1] - 2026-10-29

### Released

- GitHub Release:
  <https://github.com/lixiang117423/qPCRtools/releases/tag/v2.0.1>
- Assets: `qPCRtools_2.0.1.tar.gz`, `qPCRtools_1.0.2.tar.gz`, and the
  2022-paper impact assessment workbook/archive.

### Fixed

- `CalExp2ddCt()` now calculates the calibration ΔCt from the reference group
  only, so the reference group is centred at 1.
- `CalExp2ddCt()` reports `SE = SD / sqrt(n)` instead of `mean / sqrt(n)`.
- `CalExpCurve()` now back-transforms the fitted value to concentration,
  aggregates technical replicates within biological replicates before testing,
  and warns when `BioRep` is absent.
- `CalExpRqPCR()` keeps independent biological replicates with identical
  numeric expression values and uses relative quantities for the geNorm
  stability calculation.
- ANOVA/Tukey letter grouping now uses `level = 0.05` and the correctly
  spelled `decreasing` argument.
- Boxplots no longer draw the artificial large red dots used to reserve space.
- Bar heights are plotted on the group-mean scale rather than divided by the
  replicate count.
- `CalCurve()` no longer uses deprecated `..rr.label..` aesthetics.
- `CalExp2ddCt()`, `CalExpCurve()`, and `CalExpRqPCR()` now return the
  unadjusted (`RawPValue`) and reported (`PValue`) P value for every
  non-reference group. ANOVA/Tukey results report both the unadjusted
  pairwise P value and the Tukey-adjusted P value, with `PAdjustMethod`.
- Significance symbols now use explicit `< 0.001`, `< 0.01`, and `< 0.05`
  thresholds.

## [1.0.2] - 2026-10-29

### Released

- Published on CRAN on 2026-10-29; the CRAN package page now serves version
  1.0.2.
- Identical source tarball attached to the GitHub `v2.0.1` release.
- CRAN pre-tests requested one documentation fix: the six package manuals now
  contain `\usage` sections, and technical abbreviations in `Description` are
  quoted for the incoming spell check. The corrected tarball passes local
  `R CMD check --as-cran --no-manual` with `Status: OK` and was resubmitted to
  CRAN on 2026-10-29.

### Fixed

- Backported the same numerical and plotting fixes to the CRAN-source tree in
  `CRAN/qPCRtools/`.
- Replaced the `PKG_DESC.` placeholder with a complete package description.
- Removed the remaining deprecated `..rr.label..` use in `CalCurve()`.
- Added raw and reported P-value columns to the three expression-analysis
  result tables; ANOVA/Tukey tables now retain the actual adjusted P values
  alongside compact-letter grouping.
- Added `\usage` sections to all exported-function manuals and quoted
  `dCt`, `ddCt` and `RqPCR` in `Description` to satisfy CRAN incoming checks.
- Added a regression test with unequal `BioRep` counts and an IQR outlier;
  this closes the final open issue (#1).

## [Unreleased] — 2026-10-15（针对 CRAN 目录 `CRAN/qPCRtools/` 的 1.0.1 版）

> 本条记录的是 2026-10-15 在 `CRAN/qPCRtools/` 下做的修改（**尚未提交**；该工作区此前已含一批
> 未提交的改动），来源是导师 2026-10-14 的审稿意见与本仓库 issue #1 / #5。

### Added（2026-10-15 第二轮）

- `tests/01_verify_calexpcurve.R`：随 `R CMD check` 运行的验证脚本，覆盖 ① 反对数（表达量还原为浓度）、
  ② 技术重复层级（技术孔先归并到生物学重复、组间检验以生物学重复为单位；判别用例：按重复 p = 0.074 → NS，
  而按孔 p = 1.2e-4 → ***）、③ 缺 `BioRep` 时必须告警。

### Fixed（2026-10-15 第二轮，`CalExpCurve()` 第 5 条）

- **反对数**：`expre` 由 `(Cq - Intercept) / Slope`（对数浓度）改为
  `dilution^((Cq - Intercept) / Slope)`（浓度）。原实现把两个对数值相除，结果没有意义
  （`CalCurve()` 内部先取 `log(Conc, base = dilution)`，故新增 `dilution = 4` 参数与之一致，并同步
  roxygen 与 `man/CalExpCurve.Rd`）。
- **技术重复层级**：先在 `(Gene, Treatment, BioRep)` 内平均技术孔，再以**生物学重复**为单位做组间检验与作图
  （原实现按孔统计，同一 BioRep 的多个技术孔被当成独立样本 → 伪重复、p 值偏小）；设计表缺 `BioRep` 列时
  给出明确 warning，提示需要补列才能按生物学重复检验。
- `fig.type = "bar"` 的柱高由 `mean.expre / n` 改为 `mean.expre`（除以重复数没有意义）。

### Verified（2026-10-15 第二轮）

- `Rscript tests/01_verify_calexpcurve.R` 通过（0.8 s）。修正后示例数据的定性结论需在 R 包发布修正版后重跑。

### Fixed

- `CalExp2ddCt()`：ΔΔCt 的校准基准改为**只取对照组**（Livak & Schmittgen 2001）。
  原实现（含此前只加了 `df.sub.ck` 却没换用的半成品）用全部样本的平均 ΔCt，导致对照组 ≠ 1、
  换 `ref_group` 结果不变（issue #5 / #1）。
- `CalExp2ddCt()`：`se.expression` 由 `mean.expression / sqrt(n.biorep)` 改为
  `sd.expression / sqrt(n.biorep)`（SE = SD/√n）。
- `CalExp2ddCt()` / `CalExpCurve()` / `CalExpRqPCR()`：`multcomp::cld()` 的 `level` 由 `0.95` 改为
  `0.05`（该参数是显著性水平，0.95 等于用 α = 0.95 分组，会把无差异的组标成不同字母），
  并修正 `ddecreasing` → `decreasing`（原拼写被忽略）。
- `CalExpRqPCR()`：组内汇总不再按数值去重 —— `mean(unique(expression))` /
  `sd(unique(expression))` 改为 `mean(expression)` / `sd(expression)`，`se = sd/√n`；
  数值碰巧相同的生物学重复仍是独立观测（Hilliou & Tran 2013 §2.2.4）。
- `CalExpRqPCR()`：geNorm 稳定性 M 改用**相对量** `Q = 2^(minCq − Cq)`（原用原始 Cq 相除
  `log2(Cq_j / Cq_k)`）；实现上用 `sweep()` 按列减最小值。
- `CalExpRqPCR()`：`geometric.mean()` 从"自动选内参"分支提到函数开头 —— 原来显式给 `ref_gene` 时
  会在后半段报 `没有 geometric.mean 这个函数`。
- `CalExpRqPCR()`：删除自动选内参分支中一段从未被使用的 M 计算（与逐步剔除循环里的重复实现）。

### Verified

- 三组无差异数据上复现并确认字母分组问题：`cld(level = 0.95)` → a,a,b（假阳性）；`cld()` → a,a,a。
- 修正后 `CalExpRqPCR()` 自动选内参由 `OSPOX8, OsRBBI2` 变为 `OsUBQ, OsRBBI2`，与 C++ 版一致；
  显式 `ref_gene = "OsUBQ"` 时的组均值与 C++ 修正口径一致到 3.9e-11。

### 待办（同一批审查意见，尚未修）

- GSL 缺失时的占位 p 值（意见第 1 条）。
- `CalExpCurve()` 技术重复层级：应在生物学重复内先合并技术孔再做组间检验（意见第 5 条）。
- 本仓库根目录的 2.0.0 开发版（`R/`，snake_case 命名）存在同样的问题，尚未同步这批修正。

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
