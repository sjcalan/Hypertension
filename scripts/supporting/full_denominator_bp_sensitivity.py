#!/usr/bin/env python3
"""Create retention-denominator sensitivity plots for BP outcomes.

These plots are exploratory. They apply the cluster-size denominator logic used
for retention to quarterly BP outcomes so we can inspect how the trajectories
look when denominator-driven uncertainty is reduced.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import PchipInterpolator


BASE_DIR = Path("/userhome/cs3/u3011656/hypertension/hypertension")
OUT_DIR = BASE_DIR / "analysis" / "outputs" / "q1_ng4_cluster_presentation_2026-05-03"
FIG_DIR = OUT_DIR / "ng4_cluster_difference_exploration"
PATIENT_OUTCOMES = OUT_DIR / "q1_ng4_cluster_patient_outcomes.csv"
QUARTERLY_SUMMARY = FIG_DIR / "quarterly_clinical_trajectory_summary_ng4.csv"
Z_90 = 1.6448536269514722

CLUSTER_COLORS = {
    1: "#d94b4b",
    2: "#2c92c9",
    3: "#9a6b00",
    4: "#18b84f",
}


def wilson_interval(successes: int, n: int, z: float = Z_90) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    p = successes / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * np.sqrt((p * (1 - p) / n) + (z**2 / (4 * n**2))) / denom
    return (max(0.0, centre - margin), min(1.0, centre + margin))


def smooth_series(
    x: np.ndarray,
    y: np.ndarray,
    grid: np.ndarray,
    *,
    ymin: float | None = None,
    ymax: float | None = None,
) -> np.ndarray:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) == 0:
        out = np.full_like(grid, np.nan, dtype=float)
    elif len(x) == 1:
        out = np.full_like(grid, y[0], dtype=float)
    else:
        out = PchipInterpolator(x, y)(grid)
    if ymin is not None or ymax is not None:
        out = np.clip(out, ymin if ymin is not None else -np.inf, ymax if ymax is not None else np.inf)
    return out


def build_sensitivity_summary() -> pd.DataFrame:
    patient_df = pd.read_csv(PATIENT_OUTCOMES, usecols=["new_patient_id", "gbtm_cluster"])
    cluster_sizes = patient_df.groupby("gbtm_cluster")["new_patient_id"].nunique().astype(int)
    summary = pd.read_csv(QUARTERLY_SUMMARY)

    rows: list[dict[str, float | int]] = []
    for _, row in summary.iterrows():
        cluster = int(row["cluster"])
        quarter = int(row["quarter"])
        cluster_n = int(cluster_sizes.loc[cluster])
        control_success = int(row["control_success"])
        control_rate_full = control_success / cluster_n
        control_lo, control_hi = wilson_interval(control_success, cluster_n)

        bp_n = int(row["bp_n"])

        def full_denominator_mean_interval(prefix: str) -> tuple[float, float]:
            mean = float(row[f"{prefix}_mean"])
            observed_low = float(row[f"{prefix}_ci_low"])
            observed_high = float(row[f"{prefix}_ci_high"])
            observed_margin = (observed_high - observed_low) / 2
            if bp_n <= 0 or cluster_n <= 0 or not np.isfinite(observed_margin):
                return (np.nan, np.nan)
            full_margin = observed_margin * np.sqrt(bp_n / cluster_n)
            return (mean - full_margin, mean + full_margin)

        sbp_lo, sbp_hi = full_denominator_mean_interval("sbp")
        dbp_lo, dbp_hi = full_denominator_mean_interval("dbp")

        rows.append(
            {
                "cluster": cluster,
                "quarter": quarter,
                "month": int(row["month"]),
                "cluster_n": cluster_n,
                "bp_observed_n": bp_n,
                "control_success": control_success,
                "control_rate_full_denominator": control_rate_full,
                "control_ci_low_full_denominator": control_lo,
                "control_ci_high_full_denominator": control_hi,
                "control_rate_full_denominator_percent": control_rate_full * 100,
                "control_ci_low_full_denominator_percent": control_lo * 100,
                "control_ci_high_full_denominator_percent": control_hi * 100,
                "sbp_mean_observed": float(row["sbp_mean"]),
                "sbp_ci_low_full_denominator": sbp_lo,
                "sbp_ci_high_full_denominator": sbp_hi,
                "dbp_mean_observed": float(row["dbp_mean"]),
                "dbp_ci_low_full_denominator": dbp_lo,
                "dbp_ci_high_full_denominator": dbp_hi,
            }
        )
    return pd.DataFrame(rows)


def plot_metric(
    df: pd.DataFrame,
    value_col: str,
    low_col: str,
    high_col: str,
    ylabel: str,
    filename: str,
    *,
    proportion: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 6.6))
    for cluster, sub in df.groupby("cluster", sort=True):
        sub = sub.sort_values("quarter")
        x = sub["quarter"].to_numpy(dtype=float)
        y = sub[value_col].to_numpy(dtype=float)
        lo = sub[low_col].to_numpy(dtype=float)
        hi = sub[high_col].to_numpy(dtype=float)
        grid = np.linspace(1, 12, 240)
        color = CLUSTER_COLORS[int(cluster)]
        ymin = 0 if proportion else None
        ymax = 1 if proportion else None
        ax.fill_between(x, lo, hi, color=color, alpha=0.10, linewidth=0, zorder=1)
        ax.plot(
            grid,
            smooth_series(x, y, grid, ymin=ymin, ymax=ymax),
            color=color,
            linewidth=2.7,
            solid_capstyle="round",
            label=f"Cluster {cluster}",
        )
        ax.scatter(x, y, color=color, edgecolor="white", linewidth=0.8, s=45, zorder=4)

    ax.set_xlabel("Quarters since first visit")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.8, 12.2)
    ax.set_xticks(range(1, 13))
    ax.grid(True, color="#e8e8e8")
    ax.set_axisbelow(True)
    if proportion:
        vals = df[value_col].dropna()
        bottom = max(0.0, vals.min() - 0.04)
        top = min(1.0, vals.max() + 0.04)
        ax.set_ylim(bottom, top)
    else:
        vals = pd.concat([df[value_col], df[low_col], df[high_col]], axis=0).dropna()
        pad = (vals.max() - vals.min()) * 0.12 if len(vals) else 1
        ax.set_ylim(vals.min() - pad, vals.max() + pad)
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.09), frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(FIG_DIR / f"{filename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.25)
    sensitivity = build_sensitivity_summary()
    sensitivity.to_csv(FIG_DIR / "quarterly_full_denominator_bp_sensitivity_ng4.csv", index=False)

    plot_metric(
        sensitivity,
        "control_rate_full_denominator_percent",
        "control_ci_low_full_denominator_percent",
        "control_ci_high_full_denominator_percent",
        "Control rate (%)",
        "quarterly_control_rate_full_denominator_ng4",
    )
    plot_metric(
        sensitivity,
        "sbp_mean_observed",
        "sbp_ci_low_full_denominator",
        "sbp_ci_high_full_denominator",
        "Mean SBP, mmHg",
        "quarterly_sbp_full_denominator_ci_ng4",
    )
    plot_metric(
        sensitivity,
        "dbp_mean_observed",
        "dbp_ci_low_full_denominator",
        "dbp_ci_high_full_denominator",
        "Mean DBP, mmHg",
        "quarterly_dbp_full_denominator_ci_ng4",
    )
    print(f"Wrote full-denominator sensitivity plots to {FIG_DIR}")


if __name__ == "__main__":
    main()
