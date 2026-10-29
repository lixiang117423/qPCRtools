#!/usr/bin/env python3
"""Compare qPCRtools 1.0.1 and 1.0.2 outputs on the 2022 paper examples."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"


def load(version: str, case: str) -> pd.DataFrame:
    return pd.read_csv(OUT / f"{version}_{case}.csv")


def summarize(
    old: pd.DataFrame,
    new: pd.DataFrame,
    keys: list[str],
    mean_col: str,
    signif_col: str | None,
) -> pd.DataFrame:
    cols = keys + [mean_col]
    if signif_col:
        cols.append(signif_col)
    old = old[cols].drop_duplicates()
    new = new[cols].drop_duplicates()
    old = old.rename(columns={mean_col: "old_mean"})
    new = new.rename(columns={mean_col: "new_mean"})
    if signif_col:
        old = old.rename(columns={signif_col: "old_significance"})
        new = new.rename(columns={signif_col: "new_significance"})
    merged = old.merge(new, on=keys, how="outer", validate="one_to_one")
    merged["absolute_change"] = merged["new_mean"] - merged["old_mean"]
    merged["percent_change"] = (
        100 * merged["absolute_change"] / merged["old_mean"].replace(0, np.nan)
    )
    if signif_col:
        merged["old_significance"] = merged["old_significance"].fillna("—")
        merged["new_significance"] = merged["new_significance"].fillna("—")
    return merged


def main() -> None:
    curve_summary = summarize(
        load("v1.0.1", "calexpcurve"),
        load("v1.0.2", "calexpcurve"),
        ["Gene", "Treatment"],
        "mean.expre",
        "signif",
    )
    ddct_summary = summarize(
        load("v1.0.1", "calexp2ddct"),
        load("v1.0.2", "calexp2ddct"),
        ["gene", "Treatment"],
        "mean.expre",
        "signif",
    )
    rq_summary = summarize(
        load("v1.0.1", "calexprqpcr"),
        load("v1.0.2", "calexprqpcr"),
        ["gene", "group"],
        "Expression",
        None,
    )
    for name, frame in (
        ("summary_calexpcurve.csv", curve_summary),
        ("summary_calexp2ddct.csv", ddct_summary),
        ("summary_calexprqpcr.csv", rq_summary),
    ):
        frame.to_csv(OUT / name, index=False)

    old_curve = load("v1.0.1", "calcurve")
    new_curve = load("v1.0.2", "calcurve")
    numeric_cols = [
        c for c in old_curve.columns
        if c in new_curve.columns and pd.api.types.is_numeric_dtype(old_curve[c])
    ]
    curve_max_abs = float(
        np.nanmax(
            np.abs(
                old_curve[numeric_cols].to_numpy()
                - new_curve[numeric_cols].to_numpy()
            )
        )
    )

    report = [
        "# qPCRtools 2022 paper impact assessment",
        "",
        "Compared qPCRtools 1.0.1 with the corrected 1.0.2 implementation on the",
        "example datasets distributed with the 2022 Frontiers in Genetics paper.",
        "",
        "## Overall finding",
        "",
        "The numerical results feeding Figures 2–4 are affected. The paper’s",
        "conclusions must be re-run from the corrected package before those figures",
        "or source-data sheets are reused.",
        "",
        "## Output row counts",
        "",
        "| Case | 1.0.1 rows | 1.0.2 rows |",
        "| --- | ---: | ---: |",
    ]
    for case in ("calcurve", "calexpcurve", "calexp2ddct", "calexprqpcr"):
        report.append(
            f"| {case} | {len(load('v1.0.1', case))} | "
            f"{len(load('v1.0.2', case))} |"
        )
    report.extend(
        [
            "",
            "## Figure 1 / CalCurve",
            "",
            "All numeric standard-curve entries are unchanged "
            f"(maximum absolute difference `{curve_max_abs:.2e}`). The figure code "
            "changed only to remove deprecated `..rr.label..` syntax.",
            "",
            "## Figure 2 / CalExpCurve",
            "",
            "The corrected version back-transforms log concentration to concentration",
            "and averages technical wells within each biological replicate. The result",
            "table falls from 80 to 40 rows; plotted group values and significance",
            "calls must be regenerated.",
            "",
            "| Gene / treatment | 1.0.1 mean | 1.0.2 mean | Change | Significance |",
            "| --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in curve_summary.itertuples():
        report.append(
            f"| {row.Gene} / {row.Treatment} | {row.old_mean:.4g} | "
            f"{row.new_mean:.4g} | {row.percent_change:+.1f}% | "
            f"{row.old_significance} → {row.new_significance} |"
        )
    report.extend(
        [
            "",
            "## Figure 3 / CalExp2ddCt",
            "",
            "The corrected version calibrates ΔCt within the control group and reports",
            "SE = SD/sqrt(n). Group means and error bars change; the bundled example",
            "keeps the same significance-symbol pattern.",
            "",
            "| Gene / treatment | 1.0.1 mean | 1.0.2 mean | Change | Significance |",
            "| --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in ddct_summary.itertuples():
        report.append(
            f"| {row.gene} / {row.Treatment} | {row.old_mean:.4g} | "
            f"{row.new_mean:.4g} | {row.percent_change:+.1f}% | "
            f"{row.old_significance} → {row.new_significance} |"
        )
    report.extend(
        [
            "",
            "## Figure 4 / CalExpRqPCR",
            "",
            "The corrected version no longer merges equal numeric expression values and",
            "uses relative quantities for geNorm. Group means and standard errors change;",
            "the published figure must be regenerated with 1.0.2.",
            "",
            "| Gene / group | 1.0.1 mean | 1.0.2 mean | Change |",
            "| --- | ---: | ---: | ---: |",
        ]
    )
    for row in rq_summary.itertuples():
        report.append(
            f"| {row.gene} / {row.group} | {row.old_mean:.4g} | "
            f"{row.new_mean:.4g} | {row.percent_change:+.1f}% |"
        )
    report.extend(
        [
            "",
            "## Reproduction",
            "",
            "Run `01_run_version_case.R` separately with R libraries containing",
            "qPCRtools 1.0.1 and 1.0.2, then run `02_compare_versions.py`.",
            "",
            "The complete numeric tables are retained in `results/`.",
            "",
        ]
    )
    (ROOT / "README.md").write_text("\n".join(report), encoding="utf-8")
    print("wrote", ROOT / "README.md")


if __name__ == "__main__":
    main()
