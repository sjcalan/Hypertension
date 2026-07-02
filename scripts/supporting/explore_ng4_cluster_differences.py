#!/usr/bin/env python3
"""Explore ng=4 cluster differences beyond retention trajectories."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


BASE_DIR = Path("/userhome/cs3/u3011656/hypertension/hypertension")
OUT_DIR = BASE_DIR / "analysis" / "outputs" / "q1_ng4_cluster_presentation_2026-05-03"
FIG_DIR = OUT_DIR / "ng4_cluster_difference_exploration"
PATIENT_OUTCOMES = OUT_DIR / "q1_ng4_cluster_patient_outcomes.csv"
VISIT_LEVEL = BASE_DIR / "visit_level_cohort.csv"
QUARTER_DAYS = 90
QUARTERLY_CI_Z = 1.6448536269514722  # 90% confidence interval; narrower than 95%.

CLUSTER_COLORS = {
    1: "#d94b4b",
    2: "#2c92c9",
    3: "#9a6b00",
    4: "#18b84f",
}


@dataclass(frozen=True)
class BinaryMetric:
    label: str
    column: str
    denominator: str = "observed"


@dataclass(frozen=True)
class ContinuousMetric:
    label: str
    column: str


BINARY_TRAJECTORY_METRICS = {
    "Control rate": [
        BinaryMetric("12 mo", "control_12mo_visitlevel"),
        BinaryMetric("24 mo", "control_24mo_visitlevel"),
        BinaryMetric("36 mo", "control_36mo_visitlevel"),
    ],
    "Treatment rate": [
        BinaryMetric("12 mo", "treated_12mo"),
        BinaryMetric("24 mo", "treated_24mo"),
        BinaryMetric("36 mo", "treated_36mo"),
    ],
}

CONTINUOUS_TRAJECTORY_METRICS = {
    "SBP, mmHg": [
        ContinuousMetric("12 mo", "sbp_12mo"),
        ContinuousMetric("24 mo", "sbp_24mo"),
        ContinuousMetric("36 mo", "sbp_36mo"),
    ],
    "DBP, mmHg": [
        ContinuousMetric("12 mo", "dbp_12mo"),
        ContinuousMetric("24 mo", "dbp_24mo"),
        ContinuousMetric("36 mo", "dbp_36mo"),
    ],
}

CONTROL_DEFINITION_METRICS = [
    BinaryMetric("12 mo control", "control_12mo_visitlevel"),
    BinaryMetric("24 mo control", "control_24mo_visitlevel"),
    BinaryMetric("36 mo control", "control_36mo_visitlevel"),
    BinaryMetric("Last BP visit controlled", "last_bp_controlled"),
    BinaryMetric("Any 3 controlled visits in a row", "control_streak3_any"),
]

SCREEN_BINARY_METRICS = [
    BinaryMetric("Retained at 6 mo", "retained_6mo"),
    BinaryMetric("Retained at 12 mo", "retained_12mo"),
    BinaryMetric("Retained at 18 mo", "retained_18mo"),
    BinaryMetric("Retained at 24 mo", "retained_24mo"),
    BinaryMetric("Retained at 30 mo", "retained_30mo"),
    BinaryMetric("Retained at 36 mo", "retained_36mo"),
    BinaryMetric("Any visit around 12 mo", "any_visit_observed_12mo"),
    BinaryMetric("Any visit around 24 mo", "any_visit_observed_24mo"),
    BinaryMetric("Any visit around 36 mo", "any_visit_observed_36mo"),
    BinaryMetric("BP visit around 12 mo", "bp_visit_observed_12mo"),
    BinaryMetric("BP visit around 24 mo", "bp_visit_observed_24mo"),
    BinaryMetric("BP visit around 36 mo", "bp_visit_observed_36mo"),
    BinaryMetric("Treatment at 12 mo", "treated_12mo"),
    BinaryMetric("Treatment at 24 mo", "treated_24mo"),
    BinaryMetric("Treatment at 36 mo", "treated_36mo"),
    BinaryMetric("Control at 12 mo", "control_12mo_visitlevel"),
    BinaryMetric("Control at 24 mo", "control_24mo_visitlevel"),
    BinaryMetric("Control at 36 mo", "control_36mo_visitlevel"),
    BinaryMetric("Last BP visit controlled", "last_bp_controlled"),
    BinaryMetric("Any 3 controlled visits in a row", "control_streak3_any"),
]

SCREEN_CONTINUOUS_METRICS = [
    ContinuousMetric("Total visits", "total_visits"),
    ContinuousMetric("Time to last visit, days", "time_to_last_visit"),
    ContinuousMetric("Median visit gap, days", "visit_duration_median"),
    ContinuousMetric("Time to first gap >6 mo, days", "time_to_gap_6mo"),
    ContinuousMetric("Time to first gap >12 mo, days", "time_to_gap_12mo"),
    ContinuousMetric("Time to first gap >18 mo, days", "time_to_gap_18mo"),
    ContinuousMetric("SBP at 12 mo", "sbp_12mo"),
    ContinuousMetric("SBP at 24 mo", "sbp_24mo"),
    ContinuousMetric("SBP at 36 mo", "sbp_36mo"),
    ContinuousMetric("DBP at 12 mo", "dbp_12mo"),
    ContinuousMetric("DBP at 24 mo", "dbp_24mo"),
    ContinuousMetric("DBP at 36 mo", "dbp_36mo"),
]


def wilson_interval(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    p = successes / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt((p * (1 - p) / n) + (z**2 / (4 * n**2))) / denom
    return (max(0.0, centre - margin), min(1.0, centre + margin))


def clean_cluster(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["gbtm_cluster"] = out["gbtm_cluster"].astype(int)
    return out


def binary_summary(df: pd.DataFrame, metric: BinaryMetric) -> pd.DataFrame:
    rows = []
    for cluster, sub in df.groupby("gbtm_cluster", sort=True):
        vals = sub[metric.column].dropna().astype(int)
        n = int(vals.shape[0])
        successes = int(vals.sum())
        prop = successes / n if n else np.nan
        lo, hi = wilson_interval(successes, n)
        rows.append(
            {
                "cluster": int(cluster),
                "metric": metric.label,
                "column": metric.column,
                "n": n,
                "successes": successes,
                "proportion": prop,
                "percent": prop * 100,
                "ci_low": lo,
                "ci_high": hi,
                "ci_low_percent": lo * 100,
                "ci_high_percent": hi * 100,
            }
        )
    return pd.DataFrame(rows)


def continuous_summary(df: pd.DataFrame, metric: ContinuousMetric) -> pd.DataFrame:
    rows = []
    for cluster, sub in df.groupby("gbtm_cluster", sort=True):
        vals = sub[metric.column].dropna()
        n = int(vals.shape[0])
        mean = vals.mean()
        sd = vals.std()
        se = sd / np.sqrt(n) if n > 0 else np.nan
        ci_margin = 1.96 * se if n > 1 else np.nan
        rows.append(
            {
                "cluster": int(cluster),
                "metric": metric.label,
                "column": metric.column,
                "n": n,
                "mean": mean,
                "sd": sd,
                "ci_low": mean - ci_margin if pd.notna(ci_margin) else np.nan,
                "ci_high": mean + ci_margin if pd.notna(ci_margin) else np.nan,
                "median": vals.median(),
                "q1": vals.quantile(0.25),
                "q3": vals.quantile(0.75),
            }
        )
    return pd.DataFrame(rows)


def cramer_v_for_binary(df: pd.DataFrame, column: str) -> tuple[float, float, int]:
    sub = df[["gbtm_cluster", column]].dropna().copy()
    sub[column] = sub[column].astype(int)
    table = pd.crosstab(sub["gbtm_cluster"], sub[column])
    if table.shape[0] < 2 or table.shape[1] < 2:
        return (np.nan, np.nan, int(sub.shape[0]))
    chi2, p, _, _ = stats.chi2_contingency(table)
    n = table.to_numpy().sum()
    v = np.sqrt(chi2 / (n * min(table.shape[0] - 1, table.shape[1] - 1)))
    return (float(p), float(v), int(n))


def kruskal_effect(df: pd.DataFrame, column: str) -> tuple[float, float, int]:
    groups = [sub[column].dropna().to_numpy() for _, sub in df.groupby("gbtm_cluster", sort=True)]
    groups = [g for g in groups if len(g) > 0]
    n = int(sum(len(g) for g in groups))
    if len(groups) < 2 or n <= len(groups):
        return (np.nan, np.nan, n)
    h, p = stats.kruskal(*groups)
    eps_sq = max(0.0, (h - len(groups) + 1) / (n - len(groups)))
    return (float(p), float(eps_sq), n)


def p_text(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def derive_visit_level_metrics(patient_df: pd.DataFrame) -> pd.DataFrame:
    usecols = [
        "new_patient_id",
        "visit_sequence",
        "encounter_date",
        "systolic_bp",
        "diastolic_bp",
        "controlled",
        "treated",
        "active_med_count",
    ]
    visits = pd.read_csv(VISIT_LEVEL, usecols=usecols)
    clusters = patient_df[["new_patient_id", "gbtm_cluster"]].drop_duplicates()
    visits = visits.merge(clusters, on="new_patient_id", how="inner")
    visits = visits.sort_values(["new_patient_id", "visit_sequence"])

    def first_nonmissing(s: pd.Series) -> float:
        s = s.dropna()
        return s.iloc[0] if len(s) else np.nan

    def last_nonmissing(s: pd.Series) -> float:
        s = s.dropna()
        return s.iloc[-1] if len(s) else np.nan

    patient_visit = (
        visits.groupby(["new_patient_id", "gbtm_cluster"], as_index=False)
        .agg(
            visit_count=("visit_sequence", "count"),
            bp_visit_count=("systolic_bp", lambda s: int(s.notna().sum())),
            controlled_visit_pct=("controlled", "mean"),
            treated_visit_pct=("treated", "mean"),
            mean_active_med_count=("active_med_count", "mean"),
            first_sbp=("systolic_bp", first_nonmissing),
            last_sbp=("systolic_bp", last_nonmissing),
            first_dbp=("diastolic_bp", first_nonmissing),
            last_dbp=("diastolic_bp", last_nonmissing),
        )
    )
    patient_visit["sbp_change_last_minus_first"] = patient_visit["last_sbp"] - patient_visit["first_sbp"]
    patient_visit["dbp_change_last_minus_first"] = patient_visit["last_dbp"] - patient_visit["first_dbp"]
    return patient_visit


def normal_interval(mean: float, sd: float, n: int, z: float = 1.96) -> tuple[float, float]:
    if n <= 1 or pd.isna(mean) or pd.isna(sd):
        return (np.nan, np.nan)
    margin = z * sd / np.sqrt(n)
    return (mean - margin, mean + margin)


def derive_quarterly_clinical_summary(patient_df: pd.DataFrame) -> pd.DataFrame:
    usecols = [
        "new_patient_id",
        "visit_sequence",
        "encounter_dateD",
        "systolic_bp",
        "diastolic_bp",
        "controlled",
        "treated",
    ]
    visits = pd.read_csv(VISIT_LEVEL, usecols=usecols)
    visits = visits.merge(patient_df[["new_patient_id", "gbtm_cluster"]], on="new_patient_id", how="inner")
    visits = visits.dropna(subset=["encounter_dateD", "visit_sequence"]).copy()
    visits["encounter_dateD"] = pd.to_numeric(visits["encounter_dateD"], errors="coerce")
    visits["visit_sequence"] = pd.to_numeric(visits["visit_sequence"], errors="coerce")
    visits = visits.dropna(subset=["encounter_dateD", "visit_sequence"]).copy()
    first_day = visits.groupby("new_patient_id", as_index=False)["encounter_dateD"].min().rename(columns={"encounter_dateD": "first_day"})
    visits = visits.merge(first_day, on="new_patient_id", how="left")
    visits["days_since_first"] = visits["encounter_dateD"] - visits["first_day"]
    visits = visits.loc[(visits["days_since_first"] >= 0) & (visits["days_since_first"] < 12 * QUARTER_DAYS)].copy()
    visits["quarter"] = (visits["days_since_first"] // QUARTER_DAYS).astype(int) + 1

    any_last = (
        visits.sort_values(["new_patient_id", "quarter", "days_since_first", "visit_sequence"])
        .groupby(["new_patient_id", "gbtm_cluster", "quarter"], as_index=False)
        .last()[["new_patient_id", "gbtm_cluster", "quarter", "treated"]]
    )
    bp_last = (
        visits.loc[
            visits["systolic_bp"].notna()
            & visits["diastolic_bp"].notna()
            & visits["controlled"].notna()
        ]
        .sort_values(["new_patient_id", "quarter", "days_since_first", "visit_sequence"])
        .groupby(["new_patient_id", "gbtm_cluster", "quarter"], as_index=False)
        .last()[["new_patient_id", "gbtm_cluster", "quarter", "systolic_bp", "diastolic_bp", "controlled"]]
    )

    rows = []
    clusters = sorted(patient_df["gbtm_cluster"].unique())
    for cluster in clusters:
        for quarter in range(1, 13):
            any_sub = any_last.loc[(any_last["gbtm_cluster"] == cluster) & (any_last["quarter"] == quarter)]
            bp_sub = bp_last.loc[(bp_last["gbtm_cluster"] == cluster) & (bp_last["quarter"] == quarter)]

            treated_vals = pd.to_numeric(any_sub["treated"], errors="coerce").dropna().astype(int)
            treated_n = int(treated_vals.shape[0])
            treated_success = int(treated_vals.sum())
            treated_prop = treated_success / treated_n if treated_n else np.nan
            treated_lo, treated_hi = wilson_interval(treated_success, treated_n, z=QUARTERLY_CI_Z)

            control_vals = pd.to_numeric(bp_sub["controlled"], errors="coerce").dropna().astype(int)
            control_n = int(control_vals.shape[0])
            control_success = int(control_vals.sum())
            control_prop = control_success / control_n if control_n else np.nan
            control_lo, control_hi = wilson_interval(control_success, control_n, z=QUARTERLY_CI_Z)

            sbp_vals = pd.to_numeric(bp_sub["systolic_bp"], errors="coerce").dropna()
            dbp_vals = pd.to_numeric(bp_sub["diastolic_bp"], errors="coerce").dropna()
            sbp_mean = sbp_vals.mean()
            dbp_mean = dbp_vals.mean()
            sbp_lo, sbp_hi = normal_interval(sbp_mean, sbp_vals.std(), int(sbp_vals.shape[0]), z=QUARTERLY_CI_Z)
            dbp_lo, dbp_hi = normal_interval(dbp_mean, dbp_vals.std(), int(dbp_vals.shape[0]), z=QUARTERLY_CI_Z)

            rows.append(
                {
                    "cluster": int(cluster),
                    "quarter": quarter,
                    "month": quarter * 3,
                    "control_n": control_n,
                    "control_success": control_success,
                    "control_rate": control_prop,
                    "control_ci_low": control_lo,
                    "control_ci_high": control_hi,
                    "treatment_n": treated_n,
                    "treatment_success": treated_success,
                    "treatment_rate": treated_prop,
                    "treatment_ci_low": treated_lo,
                    "treatment_ci_high": treated_hi,
                    "bp_n": int(sbp_vals.shape[0]),
                    "sbp_mean": sbp_mean,
                    "sbp_ci_low": sbp_lo,
                    "sbp_ci_high": sbp_hi,
                    "dbp_mean": dbp_mean,
                    "dbp_ci_low": dbp_lo,
                    "dbp_ci_high": dbp_hi,
                }
            )
    return pd.DataFrame(rows)


def derive_quarterly_retention_summary(patient_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    clusters = sorted(patient_df["gbtm_cluster"].unique())
    for cluster in clusters:
        sub = patient_df.loc[patient_df["gbtm_cluster"] == cluster]
        total_n = int(sub["new_patient_id"].nunique())
        for quarter in range(1, 13):
            threshold_days = quarter * QUARTER_DAYS
            retained = (pd.to_numeric(sub["time_to_last_visit"], errors="coerce") >= threshold_days).astype(int)
            retained_n = int(retained.sum())
            prop = retained_n / total_n if total_n else np.nan
            lo, hi = wilson_interval(retained_n, total_n, z=QUARTERLY_CI_Z)
            rows.append(
                {
                    "cluster": int(cluster),
                    "quarter": quarter,
                    "month": quarter * 3,
                    "retained_n": retained_n,
                    "denominator": total_n,
                    "retention_rate": prop,
                    "retention_ci_low": lo,
                    "retention_ci_high": hi,
                }
            )
    return pd.DataFrame(rows)

def plot_quarterly_metric(
    summary: pd.DataFrame,
    value: str,
    lo: str,
    hi: str,
    ylabel: str,
    title: str,
    filename: str,
    *,
    proportion: bool = False,
    zoom: bool = False,
    uncertainty: str = "ribbon",
    ribbon_alpha: float = 0.10,
) -> None:
    del title
    fig, ax = plt.subplots(figsize=(9.2, 6.6))
    for cluster, sub in summary.groupby("cluster", sort=True):
        sub = sub.sort_values("quarter")
        x = sub["quarter"].to_numpy(dtype=float)
        y = sub[value].to_numpy(dtype=float)
        low = sub[lo].to_numpy(dtype=float)
        high = sub[hi].to_numpy(dtype=float)
        color = CLUSTER_COLORS[int(cluster)]
        if uncertainty == "ribbon":
            ax.fill_between(
                x,
                low,
                high,
                color=color,
                alpha=ribbon_alpha,
                linewidth=0,
                zorder=1,
            )
        ax.plot(
            x,
            y,
            color=color,
            linewidth=2.7,
            solid_capstyle="round",
            solid_joinstyle="round",
            label=f"Cluster {cluster}",
        )
        ax.scatter(x, y, color=color, edgecolor="white", linewidth=0.8, s=45, zorder=4)
        if uncertainty == "errorbar":
            yerr = np.vstack([y - low, high - y])
            yerr = np.where(np.isfinite(yerr), yerr, np.nan)
            ax.errorbar(
                x,
                y,
                yerr=yerr,
                fmt="none",
                ecolor=color,
                elinewidth=1.1,
                capsize=3.2,
                capthick=1.0,
                alpha=0.65,
                zorder=3,
            )

    ax.set_xlabel("Quarters since first visit")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.8, 12.2)
    ax.set_xticks(range(1, 13))
    ax.grid(True, color="#e8e8e8")
    ax.set_axisbelow(True)
    if proportion:
        vals = summary[value].dropna()
        if zoom and len(vals):
            bottom = max(0.0, vals.min() - 0.10)
            top = min(1.0, vals.max() + 0.10)
            ax.set_ylim(bottom, top)
        else:
            ax.set_ylim(0, 1.02)
    else:
        vals = pd.concat([summary[value], summary[lo], summary[hi]], axis=0).dropna()
        if len(vals):
            pad = (vals.max() - vals.min()) * 0.12
            ax.set_ylim(vals.min() - pad, vals.max() + pad)
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.09), frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIG_DIR / f"{filename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_quarterly_clinical_trajectories(summary: pd.DataFrame) -> None:
    plot_quarterly_metric(
        summary,
        value="control_rate",
        lo="control_ci_low",
        hi="control_ci_high",
        ylabel="Control rate",
        title="Quarterly Blood Pressure Control Rate by ng=4 GBTM Cluster",
        filename="quarterly_control_rate_ng4",
        proportion=True,
        zoom=True,
        uncertainty="ribbon",
        ribbon_alpha=0.09,
    )
    plot_quarterly_metric(
        summary,
        value="treatment_rate",
        lo="treatment_ci_low",
        hi="treatment_ci_high",
        ylabel="Treatment rate",
        title="Quarterly Treatment Rate by ng=4 GBTM Cluster",
        filename="quarterly_treatment_rate_ng4",
        proportion=True,
        uncertainty="ribbon",
        ribbon_alpha=0.09,
    )
    plot_quarterly_metric(
        summary,
        value="sbp_mean",
        lo="sbp_ci_low",
        hi="sbp_ci_high",
        ylabel="Mean SBP, mmHg",
        title="Quarterly Mean SBP by ng=4 GBTM Cluster",
        filename="quarterly_sbp_ng4",
        uncertainty="ribbon",
        ribbon_alpha=0.09,
    )
    plot_quarterly_metric(
        summary,
        value="dbp_mean",
        lo="dbp_ci_low",
        hi="dbp_ci_high",
        ylabel="Mean DBP, mmHg",
        title="Quarterly Mean DBP by ng=4 GBTM Cluster",
        filename="quarterly_dbp_ng4",
        uncertainty="ribbon",
        ribbon_alpha=0.09,
    )


def plot_quarterly_retention_trajectory(summary: pd.DataFrame) -> None:
    plot_quarterly_metric(
        summary,
        value="retention_rate",
        lo="retention_ci_low",
        hi="retention_ci_high",
        ylabel="Retention rate",
        title="Quarterly Retention Rate by ng=4 GBTM Cluster",
        filename="quarterly_retention_rate_ng4",
        proportion=True,
    )


def make_screen_table(patient_df: pd.DataFrame, visit_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in SCREEN_BINARY_METRICS:
        p, effect, n = cramer_v_for_binary(patient_df, metric.column)
        rows.append(
            {
                "metric_type": "binary",
                "metric": metric.label,
                "column": metric.column,
                "n_observed": n,
                "overall_p": p,
                "effect_size": effect,
                "effect_size_name": "Cramer's V",
            }
        )

    continuous_source = pd.concat(
        [
            patient_df,
            visit_metrics.drop(columns=["gbtm_cluster"], errors="ignore"),
        ],
        axis=1,
    )
    extra_continuous = [
        ContinuousMetric("Across all visits: percent BP visits controlled", "controlled_visit_pct"),
        ContinuousMetric("Across all visits: percent visits treated", "treated_visit_pct"),
        ContinuousMetric("Across all visits: mean active medication count", "mean_active_med_count"),
        ContinuousMetric("Across all visits: SBP change, last-first", "sbp_change_last_minus_first"),
        ContinuousMetric("Across all visits: DBP change, last-first", "dbp_change_last_minus_first"),
        ContinuousMetric("Across all visits: BP-recorded visit count", "bp_visit_count"),
    ]
    for metric in SCREEN_CONTINUOUS_METRICS + extra_continuous:
        p, effect, n = kruskal_effect(continuous_source, metric.column)
        rows.append(
            {
                "metric_type": "continuous",
                "metric": metric.label,
                "column": metric.column,
                "n_observed": n,
                "overall_p": p,
                "effect_size": effect,
                "effect_size_name": "Kruskal epsilon-squared",
            }
        )

    screen = pd.DataFrame(rows)
    screen = screen.sort_values(["effect_size", "n_observed"], ascending=[False, False])
    screen["overall_p_text"] = screen["overall_p"].map(p_text)
    return screen


def plot_clinical_trajectories(patient_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=False)
    axes = axes.ravel()

    for ax, (title, metrics) in zip(axes[:2], BINARY_TRAJECTORY_METRICS.items()):
        plot_rows = []
        for metric in metrics:
            tmp = binary_summary(patient_df, metric)
            tmp["time"] = metric.label
            plot_rows.append(tmp)
        data = pd.concat(plot_rows, ignore_index=True)
        for cluster, sub in data.groupby("cluster", sort=True):
            x = np.arange(len(metrics))
            y = sub["proportion"].to_numpy()
            lo = sub["ci_low"].to_numpy()
            hi = sub["ci_high"].to_numpy()
            color = CLUSTER_COLORS[int(cluster)]
            ax.fill_between(x, lo, hi, color=color, alpha=0.09, linewidth=0)
            ax.plot(x, y, color=color, linewidth=2.6, marker="o", markersize=6, label=f"Cluster {cluster}")
        ax.set_title(title)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("Proportion")
        ax.set_xticks(np.arange(len(metrics)))
        ax.set_xticklabels([m.label for m in metrics])
        ax.grid(True, color="#e8e8e8")

    for ax, (title, metrics) in zip(axes[2:], CONTINUOUS_TRAJECTORY_METRICS.items()):
        plot_rows = []
        for metric in metrics:
            tmp = continuous_summary(patient_df, metric)
            tmp["time"] = metric.label
            plot_rows.append(tmp)
        data = pd.concat(plot_rows, ignore_index=True)
        for cluster, sub in data.groupby("cluster", sort=True):
            x = np.arange(len(metrics))
            y = sub["mean"].to_numpy()
            lo = sub["ci_low"].to_numpy()
            hi = sub["ci_high"].to_numpy()
            color = CLUSTER_COLORS[int(cluster)]
            ax.fill_between(x, lo, hi, color=color, alpha=0.09, linewidth=0)
            ax.plot(x, y, color=color, linewidth=2.6, marker="o", markersize=6, label=f"Cluster {cluster}")
        ax.set_title(title)
        ax.set_ylabel("Mean")
        ax.set_xticks(np.arange(len(metrics)))
        ax.set_xticklabels([m.label for m in metrics])
        ax.grid(True, color="#e8e8e8")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Clinical Outcome Trajectories by ng=4 GBTM Cluster", y=1.08, fontsize=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "clinical_outcome_trajectories_ng4.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "clinical_outcome_trajectories_ng4.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_control_definition_bars(patient_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in CONTROL_DEFINITION_METRICS:
        rows.append(binary_summary(patient_df, metric))
    data = pd.concat(rows, ignore_index=True)

    fig, ax = plt.subplots(figsize=(12.5, 6.2))
    metrics = [m.label for m in CONTROL_DEFINITION_METRICS]
    clusters = sorted(data["cluster"].unique())
    x = np.arange(len(metrics), dtype=float)
    width = 0.18
    offsets = np.linspace(-1.5 * width, 1.5 * width, len(clusters))

    for offset, cluster in zip(offsets, clusters):
        sub = data.loc[data["cluster"] == cluster].set_index("metric").loc[metrics]
        y = sub["proportion"].to_numpy()
        yerr = np.vstack(
            [
                y - sub["ci_low"].to_numpy(),
                sub["ci_high"].to_numpy() - y,
            ]
        )
        ax.bar(
            x + offset,
            y,
            width=width,
            color=CLUSTER_COLORS[int(cluster)],
            edgecolor="white",
            linewidth=0.8,
            label=f"Cluster {cluster}",
        )
        ax.errorbar(
            x + offset,
            y,
            yerr=yerr,
            fmt="none",
            ecolor="#303030",
            elinewidth=0.9,
            capsize=2.5,
            capthick=0.9,
        )

    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Proportion")
    ax.set_xlabel("")
    ax.set_title("Blood Pressure Control Definitions by ng=4 GBTM Cluster", fontweight="bold", pad=34)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=18, ha="right")
    ax.grid(True, axis="y", color="#e8e8e8")
    ax.set_axisbelow(True)
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.12), frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / "control_definition_comparison_ng4.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "control_definition_comparison_ng4.pdf", bbox_inches="tight")
    plt.close(fig)
    return data


def plot_engagement_distributions(patient_df: pd.DataFrame, visit_metrics: pd.DataFrame) -> None:
    merged = patient_df.merge(
        visit_metrics[
            [
                "new_patient_id",
                "controlled_visit_pct",
                "treated_visit_pct",
                "mean_active_med_count",
                "sbp_change_last_minus_first",
                "dbp_change_last_minus_first",
            ]
        ],
        on="new_patient_id",
        how="left",
    )
    long_specs = [
        ("Total visits", "total_visits"),
        ("Time to last visit, days", "time_to_last_visit"),
        ("Median visit gap, days", "visit_duration_median"),
        ("Percent BP visits controlled", "controlled_visit_pct"),
        ("Percent visits treated", "treated_visit_pct"),
        ("Mean active medication count", "mean_active_med_count"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8.4))
    axes = axes.ravel()
    for ax, (title, col) in zip(axes, long_specs):
        sub = merged[["gbtm_cluster", col]].dropna().copy()
        if "Percent" in title:
            sub[col] = sub[col] * 100
        sns.boxplot(
            data=sub,
            x="gbtm_cluster",
            y=col,
            hue="gbtm_cluster",
            palette=CLUSTER_COLORS,
            showfliers=False,
            ax=ax,
            legend=False,
        )
        ax.set_title(title)
        ax.set_xlabel("Cluster")
        ax.set_ylabel("")
        ax.grid(True, axis="y", color="#e8e8e8")
        ax.set_axisbelow(True)
    fig.suptitle("Care Engagement and Visit-Level Clinical Profiles", y=1.02, fontsize=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "care_engagement_and_visit_profile_ng4.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "care_engagement_and_visit_profile_ng4.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_cluster_fingerprint(patient_df: pd.DataFrame, visit_metrics: pd.DataFrame) -> None:
    merged = patient_df.merge(
        visit_metrics[
            [
                "new_patient_id",
                "controlled_visit_pct",
                "treated_visit_pct",
                "mean_active_med_count",
                "sbp_change_last_minus_first",
                "dbp_change_last_minus_first",
            ]
        ],
        on="new_patient_id",
        how="left",
    )
    specs = [
        ("Retained 36 mo", "retained_36mo", "mean"),
        ("Any 3 controlled visits", "control_streak3_any", "mean"),
        ("Last BP controlled", "last_bp_controlled", "mean"),
        ("Treatment 36 mo", "treated_36mo", "mean"),
        ("All BP visits controlled", "controlled_visit_pct", "median"),
        ("All visits treated", "treated_visit_pct", "median"),
        ("Total visits", "total_visits", "median"),
        ("Time to last visit", "time_to_last_visit", "median"),
        ("Median visit gap", "visit_duration_median", "median"),
    ]

    rows = []
    for label, col, agg in specs:
        if agg == "mean":
            vals = merged.groupby("gbtm_cluster")[col].mean()
        else:
            vals = merged.groupby("gbtm_cluster")[col].median()
        for cluster, value in vals.items():
            rows.append({"metric": label, "cluster": int(cluster), "value": value})
    profile = pd.DataFrame(rows)

    wide = profile.pivot(index="metric", columns="cluster", values="value").loc[[s[0] for s in specs]]
    color_values = wide.copy()
    for idx in color_values.index:
        row = color_values.loc[idx]
        row_min = row.min(skipna=True)
        row_max = row.max(skipna=True)
        if pd.isna(row_min) or pd.isna(row_max) or row_max == row_min:
            color_values.loc[idx] = 0.5
        else:
            color_values.loc[idx] = (row - row_min) / (row_max - row_min)

    def format_heatmap_value(metric: str, value: float) -> str:
        if pd.isna(value):
            return ""
        if metric in {
            "Retained 36 mo",
            "Any 3 controlled visits",
            "Last BP controlled",
            "Treatment 36 mo",
            "All BP visits controlled",
            "All visits treated",
        }:
            return f"{value * 100:.1f}%"
        if metric in {"Total visits"}:
            return f"{value:.0f}"
        return f"{value:.0f}"

    annot = wide.astype(object).copy()
    for metric in annot.index:
        annot.loc[metric] = [format_heatmap_value(metric, v) for v in wide.loc[metric]]

    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    sns.heatmap(
        color_values,
        cmap="YlGnBu",
        vmin=0,
        vmax=1,
        linewidths=0.7,
        linecolor="white",
        cbar_kws={"label": "Relative value within each row"},
        annot=annot,
        fmt="",
        ax=ax,
    )
    ax.set_xlabel("Cluster")
    ax.set_ylabel("")
    ax.set_title("ng=4 Cluster Fingerprint Across Retention, Control, and Engagement", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cluster_fingerprint_heatmap_ng4.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "cluster_fingerprint_heatmap_ng4.pdf", bbox_inches="tight")
    plt.close(fig)

    wide.to_csv(FIG_DIR / "cluster_fingerprint_raw_values_ng4.csv")


def write_recommendations(screen: pd.DataFrame) -> None:
    top = screen.head(12).copy()
    lines = [
        "# ng=4 Cluster Difference Plot Recommendations",
        "",
        "The simple retention trajectory is statistically different, but it is monotonic rather than crossing like the reference paper figure.",
        "For the paper, the strongest visual candidates are:",
        "",
        "1. `care_engagement_and_visit_profile_ng4.png`: best for showing that the clusters represent different care-use patterns.",
        "2. `control_definition_comparison_ng4.png`: best for showing clinically interpretable BP-control differences, especially durable control.",
        "3. `clinical_outcome_trajectories_ng4.png`: best if you want a trajectory-style clinical figure with control rate, treatment rate, SBP, and DBP.",
        "4. `cluster_fingerprint_heatmap_ng4.png`: best compact summary figure if space is limited.",
        "",
        "Interval note:",
        "",
        "- Quarterly binary trajectory panels use Wilson 90% pointwise intervals around observed proportions.",
        "- Quarterly SBP and DBP trajectory panels use normal-approximation 90% pointwise intervals around observed means.",
        "- Clinical trajectory panels draw light shaded ribbons for pointwise intervals at observed quarters. The ribbons are intentionally light because many BP intervals overlap.",
        "",
        "Top variables by screening effect size:",
        "",
    ]
    for _, row in top.iterrows():
        lines.append(
            f"- {row['metric']}: p={row['overall_p_text']}, "
            f"{row['effect_size_name']}={row['effect_size']:.3f}, n={int(row['n_observed']):,}"
        )
    lines.append("")
    lines.append("My recommendation: use the care-engagement profile as the main replacement/companion for retention, and use the control-definition comparison as the clinical outcome figure.")
    (FIG_DIR / "plot_recommendations_ng4.md").write_text("\n".join(lines), encoding="ascii")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.25)
    patient_df = clean_cluster(pd.read_csv(PATIENT_OUTCOMES))
    visit_metrics = derive_visit_level_metrics(patient_df)
    quarterly_summary = derive_quarterly_clinical_summary(patient_df)
    quarterly_retention_summary = derive_quarterly_retention_summary(patient_df)

    screen = make_screen_table(patient_df, visit_metrics)
    screen.to_csv(FIG_DIR / "ng4_cluster_difference_metric_screen.csv", index=False)

    visit_metrics.to_csv(FIG_DIR / "ng4_patient_visit_level_derived_metrics.csv", index=False)
    quarterly_summary.to_csv(FIG_DIR / "quarterly_clinical_trajectory_summary_ng4.csv", index=False)
    quarterly_retention_summary.to_csv(FIG_DIR / "quarterly_retention_trajectory_summary_ng4.csv", index=False)
    plot_clinical_trajectories(patient_df)
    plot_quarterly_clinical_trajectories(quarterly_summary)
    plot_quarterly_retention_trajectory(quarterly_retention_summary)
    control_data = plot_control_definition_bars(patient_df)
    control_data.to_csv(FIG_DIR / "control_definition_comparison_ng4_data.csv", index=False)
    plot_engagement_distributions(patient_df, visit_metrics)
    plot_cluster_fingerprint(patient_df, visit_metrics)
    write_recommendations(screen)

    print(f"Wrote exploration outputs to {FIG_DIR}")
    print(screen.head(12)[["metric", "overall_p_text", "effect_size_name", "effect_size", "n_observed"]].to_string(index=False))


if __name__ == "__main__":
    main()
