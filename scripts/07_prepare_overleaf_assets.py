from __future__ import annotations

from pathlib import Path
import itertools
import re

import numpy as np
import pandas as pd
from scipy import stats


BASE_DIR = Path("/userhome/cs3/u3011656/hypertension/hypertension")
OUT_DIR = BASE_DIR / "analysis" / "overleaf_hypertension_visit_trajectory_2026-05-25"
TABLE_DIR = OUT_DIR / "tables"
DATA_DIR = OUT_DIR / "data"

ASSIGN_PATH = (
    BASE_DIR
    / "analysis/outputs/q1_gbtm_relative_visit_only_biweekly_full_ng4_gpu_minvis3_parallel_2026-04-17/"
    / "gbtm_relative_visit_only_assignments.csv"
)
MODEL_COMPARE_PATH = BASE_DIR / "analysis/outputs/q1_gbtm_ng_compare_minvis3_2026-04-19/q1_gbtm_model_compare.csv"
BASELINE_PATH = BASE_DIR / "hypertension_cohort_with_baseline.csv"
PATIENT_OUTCOME_PATH = (
    BASE_DIR
    / "analysis/outputs/q1_ng4_cluster_presentation_2026-05-03/q1_ng4_cluster_patient_outcomes.csv"
)


CLUSTER_LABELS = {
    1: "Cluster 1",
    2: "Cluster 2",
    3: "Cluster 3",
    4: "Cluster 4",
}


def latex_escape(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def fmt_p(p: float | None) -> str:
    if p is None or pd.isna(p):
        return ""
    if p < 0.001:
        return r"$<0.001$"
    return f"{p:.3f}"


def fmt_n_pct(series: pd.Series) -> str:
    x = pd.to_numeric(series, errors="coerce")
    denom = int(x.notna().sum())
    if denom == 0:
        return "NA"
    num = int((x.dropna() >= 1).sum())
    return f"{num:,} ({100 * num / denom:.1f}\\%)"


def fmt_mean_sd(series: pd.Series, decimals: int = 1) -> str:
    x = pd.to_numeric(series, errors="coerce").dropna()
    if x.empty:
        return "NA"
    return f"{x.mean():.{decimals}f} ({x.std(ddof=1):.{decimals}f})"


def fmt_median_iqr(series: pd.Series, decimals: int = 0) -> str:
    x = pd.to_numeric(series, errors="coerce").dropna()
    if x.empty:
        return "NA"
    return (
        f"{x.median():.{decimals}f} "
        f"({x.quantile(0.25):.{decimals}f}-{x.quantile(0.75):.{decimals}f})"
    )


def p_continuous(df: pd.DataFrame, col: str) -> float | None:
    groups = [
        pd.to_numeric(df.loc[df["gbtm_cluster"] == c, col], errors="coerce").dropna()
        for c in [1, 2, 3, 4]
    ]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return None
    return float(stats.kruskal(*groups).pvalue)


def p_binary(df: pd.DataFrame, col: str) -> float | None:
    table = []
    for c in [1, 2, 3, 4]:
        x = pd.to_numeric(df.loc[df["gbtm_cluster"] == c, col], errors="coerce").dropna()
        if x.empty:
            return None
        table.append([int((x >= 1).sum()), int((x < 1).sum())])
    if any(sum(row) == 0 for row in table):
        return None
    return float(stats.chi2_contingency(table).pvalue)


def p_categorical(df: pd.DataFrame, col: str, categories: list[str]) -> float | None:
    table = []
    for c in [1, 2, 3, 4]:
        x = df.loc[df["gbtm_cluster"] == c, col].fillna("Missing").astype(str)
        table.append([(x == cat).sum() for cat in categories])
    arr = np.asarray(table)
    if arr.sum() == 0 or arr.shape[1] < 2:
        return None
    return float(stats.chi2_contingency(arr).pvalue)


def p_pairwise_binary(df: pd.DataFrame, col: str, a: int, b: int) -> float | None:
    rows = []
    for c in [a, b]:
        x = pd.to_numeric(df.loc[df["gbtm_cluster"] == c, col], errors="coerce").dropna()
        if x.empty:
            return None
        rows.append([int((x >= 1).sum()), int((x < 1).sum())])
    return float(stats.fisher_exact(rows).pvalue)


def p_pairwise_continuous(df: pd.DataFrame, col: str, a: int, b: int) -> float | None:
    x = pd.to_numeric(df.loc[df["gbtm_cluster"] == a, col], errors="coerce").dropna()
    y = pd.to_numeric(df.loc[df["gbtm_cluster"] == b, col], errors="coerce").dropna()
    if x.empty or y.empty:
        return None
    return float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)


def table_to_latex(rows: list[list[str]], alignment: str) -> str:
    body = []
    for row in rows:
        if row and row[0].startswith("\\midrule"):
            body.append(row[0])
            continue
        if row and row[0].startswith("\\addlinespace \\multicolumn"):
            body.append(row[0])
            continue
        if row and row[0].startswith("\\multicolumn"):
            body.append(row[0])
            continue
        body.append(" & ".join(row) + r" \\")
    return "\n".join(body)


def write_table1() -> None:
    assignments = pd.read_csv(ASSIGN_PATH, usecols=["new_patient_id", "gbtm_cluster"])
    baseline = pd.read_csv(BASELINE_PATH, low_memory=False)
    df = assignments.merge(baseline, on="new_patient_id", how="left")
    df["female"] = (df["gender_clean"] == "Female").astype(float)
    df["baseline_controlled"] = pd.to_numeric(df["controlled"], errors="coerce")
    df["baseline_treated"] = pd.to_numeric(df["treated"], errors="coerce")

    clusters = [1, 2, 3, 4]
    cluster_headers = [f"{CLUSTER_LABELS[c]} (N={int((df.gbtm_cluster == c).sum()):,})" for c in clusters]

    rows: list[list[str]] = []
    rows.append(["Characteristic", f"Overall\\newline N={len(df):,}", *[h.replace(" (N=", "\\newline N=").replace(")", "") for h in cluster_headers], "P"])
    rows.append([r"\midrule", "", "", "", "", "", ""])

    def add_cont(label: str, col: str, decimals: int = 1) -> None:
        rows.append(
            [
                latex_escape(label),
                fmt_mean_sd(df[col], decimals),
                *[fmt_mean_sd(df.loc[df.gbtm_cluster == c, col], decimals) for c in clusters],
                fmt_p(p_continuous(df, col)),
            ]
        )

    def add_binary(label: str, col: str) -> None:
        rows.append(
            [
                latex_escape(label),
                fmt_n_pct(df[col]),
                *[fmt_n_pct(df.loc[df.gbtm_cluster == c, col]) for c in clusters],
                fmt_p(p_binary(df, col)),
            ]
        )

    add_cont("Age, mean (SD), y", "age_at_first_visit", 1)
    add_binary("Female sex", "female")
    for label, col, cat in [
        ("Hispanic/Latino ethnicity", "race_clean", "Hispanic/Latino"),
        ("Other/unknown race/ethnicity", "race_clean", "Other/Unknown"),
        ("English language", "language_clean", "English"),
        ("Spanish language", "language_clean", "Spanish"),
    ]:
        tmp = (df[col] == cat).astype(float)
        df[f"tmp_{re.sub('[^A-Za-z0-9]+', '_', label)}"] = tmp
        add_binary(label, f"tmp_{re.sub('[^A-Za-z0-9]+', '_', label)}")

    for label, col in [
        ("Diabetes", "diabetes"),
        ("Chronic kidney disease", "ckd"),
    ]:
        add_binary(label, col)

    add_cont("Baseline systolic BP, mean (SD), mmHg", "baseline_sbp", 1)
    add_cont("Baseline diastolic BP, mean (SD), mmHg", "baseline_dbp", 1)
    add_binary("Controlled BP at baseline", "baseline_controlled")
    add_binary("Treated at baseline", "baseline_treated")

    tex = rf"""\begin{{table}}[!htbp]
\centering
\caption{{Baseline characteristics of the modeled cohort by visit-trajectory group}}
\label{{tab:baseline}}
\begingroup
\scriptsize
\setlength{{\tabcolsep}}{{2pt}}
\renewcommand{{\arraystretch}}{{1.08}}
\begin{{threeparttable}}
\begin{{tabularx}}{{\textwidth}}{{p{{0.24\textwidth}}>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}p{{0.06\textwidth}}}}
\toprule
{table_to_latex(rows[:1], "lrrrrrr")}
{table_to_latex(rows[1:], "lrrrrrr")}
\bottomrule
\end{{tabularx}}
\begin{{tablenotes}}[flushleft]
\footnotesize
\item Values are mean (SD) or No. (\%). The table intentionally shows selected baseline variables most relevant to confounding, cohort description, and clinical interpretation; additional baseline fields are retained in the source analytic files. P values compare the four trajectory groups using chi-square tests for categorical variables and Kruskal-Wallis tests for continuous variables. BP denotes blood pressure.
\end{{tablenotes}}
\end{{threeparttable}}
\endgroup
\end{{table}}
"""
    (TABLE_DIR / "table1_baseline_by_cluster.tex").write_text(tex)


def build_table2_rows(df: pd.DataFrame) -> tuple[list[list[str]], list[dict[str, str]]]:
    clusters = [1, 2, 3, 4]
    rows: list[list[str]] = []
    pairwise_rows: list[dict[str, str]] = []

    rows.append(["Measure", "C1", "C2", "C3", "C4", "P"])
    rows.append([r"\midrule", "", "", "", "", ""])
    rows.append([r"\multicolumn{6}{l}{\textit{Retention, No. (\%)}} \\[-0.3em]", "", "", "", "", ""])

    def add_binary(label: str, col: str) -> None:
        p_overall = p_binary(df, col)
        rows.append(
            [
                latex_escape(label),
                *[fmt_n_pct(df.loc[df.gbtm_cluster == c, col]) for c in clusters],
                fmt_p(p_overall),
            ]
        )
        for a, b in itertools.combinations(clusters, 2):
            pairwise_rows.append(
                {
                    "Measure": label,
                    "Comparison": f"Cluster {a} vs Cluster {b}",
                    "P value": fmt_p(p_pairwise_binary(df, col, a, b)),
                }
            )

    def add_cont(label: str, col: str) -> None:
        p_overall = p_continuous(df, col)
        rows.append(
            [
                latex_escape(label),
                *[fmt_median_iqr(df.loc[df.gbtm_cluster == c, col], 0) for c in clusters],
                fmt_p(p_overall),
            ]
        )
        for a, b in itertools.combinations(clusters, 2):
            pairwise_rows.append(
                {
                    "Measure": label,
                    "Comparison": f"Cluster {a} vs Cluster {b}",
                    "P value": fmt_p(p_pairwise_continuous(df, col, a, b)),
                }
            )

    for label, col in [
        ("12 months", "retained_12mo"),
        ("24 months", "retained_24mo"),
        ("36 months", "retained_36mo"),
    ]:
        add_binary(label, col)

    rows.append([r"\addlinespace \multicolumn{6}{l}{\textit{Care engagement, median (IQR)}} \\[-0.3em]", "", "", "", "", ""])
    for label, col in [
        ("Total visits", "total_visits"),
        ("Time to last visit, days", "time_to_last_visit"),
        ("Median visit gap, days", "visit_duration_median"),
    ]:
        add_cont(label, col)

    rows.append([r"\addlinespace \multicolumn{6}{l}{\textit{Fixed-window outcomes}} \\[-0.3em]", "", "", "", "", ""])
    add_binary("Observed visit near 36 months", "any_visit_observed_36mo")
    add_binary("BP-recorded visit near 36 months", "bp_visit_observed_36mo")
    add_binary("Treatment near 36 months", "treated_36mo")
    add_binary("BP control near 36 months", "control_36mo_visitlevel")

    rows.append([r"\addlinespace \multicolumn{6}{l}{\textit{Alternative BP-control definitions}} \\[-0.3em]", "", "", "", "", ""])
    add_binary("Any 3 consecutive controlled BP visits", "control_streak3_any")
    return rows, pairwise_rows


def write_table2() -> None:
    df = pd.read_csv(PATIENT_OUTCOME_PATH, low_memory=False)
    rows, pairwise_rows = build_table2_rows(df)

    tex = rf"""\begin{{table}}[!htbp]
\centering
\caption{{Care-engagement and clinical outcomes by visit-trajectory group}}
\label{{tab:cluster_comparison}}
\begingroup
\scriptsize
\setlength{{\tabcolsep}}{{2pt}}
\renewcommand{{\arraystretch}}{{1.08}}
\begin{{threeparttable}}
\begin{{tabularx}}{{\textwidth}}{{p{{0.33\textwidth}}>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}X>{{\centering\arraybackslash}}p{{0.06\textwidth}}}}
\toprule
{table_to_latex(rows[:1], "lrrrrr")}
{table_to_latex(rows[1:], "lrrrrr")}
\bottomrule
\end{{tabularx}}
\begin{{tablenotes}}[flushleft]
\footnotesize
\item C1, long-span sparse engagement; C2, sustained moderate engagement; C3, short-span intensive early engagement; C4, frequent intermediate engagement. Retention is defined as the final observed visit occurring on or after the indicated time since the patient's first observed visit. Fixed-window treatment, BP, and BP-control outcomes use the closest qualifying visit within $\pm45$ days of the 36-month target. Binary p values use chi-square tests; continuous p values use Kruskal-Wallis tests. BP denotes blood pressure.
\end{{tablenotes}}
\end{{threeparttable}}
\endgroup
\end{{table}}
"""
    (TABLE_DIR / "table2_cluster_comparison_corrected.tex").write_text(tex)

    pairwise = pd.DataFrame(pairwise_rows)
    pairwise.to_csv(DATA_DIR / "table2_pairwise_pvalues_corrected.csv", index=False)

    pairwise_tex_rows = [["Measure", "Comparison", "P value"]]
    for _, row in pairwise.iterrows():
        pairwise_tex_rows.append([latex_escape(row["Measure"]), latex_escape(row["Comparison"]), row["P value"]])
    pairwise_tex = rf"""\begin{{longtable}}{{p{{0.52\linewidth}}p{{0.28\linewidth}}r}}
\caption{{Pairwise cluster comparisons for Table 2 measures}}
\label{{tab:pairwise}}\\
\toprule
{table_to_latex(pairwise_tex_rows[:1], "llr")}
\midrule
\endfirsthead
\toprule
{table_to_latex(pairwise_tex_rows[:1], "llr")}
\midrule
\endhead
{table_to_latex(pairwise_tex_rows[1:], "llr")}
\bottomrule
\end{{longtable}}
"""
    (TABLE_DIR / "table_s2_pairwise_pvalues_corrected.tex").write_text(pairwise_tex)

    corrected_csv_rows = []
    header = ["Measure", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "P value"]
    for row in rows:
        if row and not row[0].startswith("\\") and row[0] != "Measure":
            corrected_csv_rows.append(dict(zip(header, [re.sub(r"\\%", "%", x) for x in row])))
    pd.DataFrame(corrected_csv_rows).to_csv(DATA_DIR / "table2_cluster_comparison_corrected.csv", index=False)


def write_model_selection_table() -> None:
    comp = pd.read_csv(MODEL_COMPARE_PATH)
    rows = [["Groups", "BIC", "Minimum class size", "Maximum class size", "Class sizes", "Decision"]]
    for _, r in comp.iterrows():
        rows.append(
            [
                str(int(r["ng"])),
                f"{float(r['bic']):,.2f}",
                f"{int(r['min_cluster_n']):,}",
                f"{int(r['max_cluster_n']):,}",
                latex_escape(str(r["cluster_sizes"]).replace("|", " / ")),
                "Primary" if str(r["recommended_flag"]) == "recommended" else "Candidate",
            ]
        )
    tex = rf"""\begin{{table}}[!htbp]
\centering
\caption{{Model comparison for visit-only group-based trajectory models}}
\label{{tab:model_selection}}
\small
\begin{{threeparttable}}
\begin{{tabular}}{{rrrrll}}
\toprule
{table_to_latex(rows[:1], "rrrrll")}
\midrule
{table_to_latex(rows[1:], "rrrrll")}
\bottomrule
\end{{tabular}}
\begin{{tablenotes}}[flushleft]
\footnotesize
\item Models used the same eligible cohort, 14-day relative-time bins, a 36-month maximum modeled window, and a minimum of 3 visits. Lower BIC indicates better fit, but the final model was selected using fit, class balance, interpretability, and clinical coherence.
\end{{tablenotes}}
\end{{threeparttable}}
\end{{table}}
"""
    (TABLE_DIR / "table_s1_model_selection.tex").write_text(tex)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    write_table1()
    write_table2()
    write_model_selection_table()


if __name__ == "__main__":
    main()
