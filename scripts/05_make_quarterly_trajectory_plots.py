#!/usr/bin/env python3
"""Regenerate quarterly clinical trajectory plots for the ng=4 manuscript.

These are the manuscript-facing clinical panels. They intentionally use the
same complete-case mean/rate denominator logic as Table 2:

- treatment rate: patients with an observed visit and nonmissing treatment
  status in the quarter
- BP control rate: patients with a BP-recorded visit in the quarter
- SBP/DBP means: patients with a BP-recorded visit in the quarter
- SBP/DBP ribbons: cluster-size scaled standard-error bands, used only to
  avoid over-emphasizing sparse observed-BP denominators in the visual display

The full-cluster-denominator sensitivity plots are useful for checking
measurement coverage, but they should not be labeled as BP control rates.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_DIR = Path("/userhome/cs3/u3011656/hypertension/hypertension")
PAPER_DIR = BASE_DIR / "analysis/outputs/q1_ng4_cluster_paper_assets_2026-05-15"
PRESENTATION_DIR = BASE_DIR / "analysis/outputs/q1_ng4_cluster_presentation_2026-05-03/ng4_cluster_difference_exploration"
OVERLEAF_FIG_DIR = BASE_DIR / "analysis/overleaf_hypertension_visit_trajectory_2026-05-25/figures"
SUMMARY_CSV = PAPER_DIR / "quarterly_clinical_trajectory_summary_ng4.csv"
RETENTION_SUMMARY_CSV = PRESENTATION_DIR / "quarterly_retention_trajectory_summary_ng4.csv"
PATIENT_OUTCOMES = BASE_DIR / "analysis/outputs/q1_ng4_cluster_presentation_2026-05-03/q1_ng4_cluster_patient_outcomes.csv"
BACKUP_DIR = PAPER_DIR / "superseded_quarterly_full_denominator_active_2026-07-02"
CI_AUDIT_CSV = PAPER_DIR / "quarterly_sbp_dbp_cluster_n_scaled_ci_audit_2026-07-02.csv"

CLUSTER_COLORS = {
    1: "#d94b4b",
    2: "#2c92c9",
    3: "#9a6b00",
    4: "#18b84f",
}

RIBBON_ALPHA = 0.22
LINE_WIDTH = 2.7
MEAN_CI_Z = 1.6448536269514722

def backup_existing_outputs(stems: list[str]) -> None:
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    for directory, prefix in [
        (PAPER_DIR, "paper_assets"),
        (PRESENTATION_DIR, "presentation"),
        (OVERLEAF_FIG_DIR, "overleaf_figures"),
    ]:
        for stem in stems:
            for ext in ("png", "pdf"):
                src = directory / f"{stem}.{ext}"
                if src.exists():
                    dest = BACKUP_DIR / f"{prefix}_{stem}.{ext}"
                    if not dest.exists():
                        shutil.copy2(src, dest)


def add_cluster_size_scaled_bp_ci(summary: pd.DataFrame) -> pd.DataFrame:
    """Use full cluster size for SBP/DBP ribbon SEs, while keeping observed means.

    The source summary stores conventional complete-case CIs:
    mean +/- z * SD / sqrt(bp_n). We infer SD from those CIs and then rescale
    the standard error using the trajectory-group size. This is a visualization
    choice requested for narrower ribbons; the means remain complete-case means
    among BP-recorded patient-quarters.
    """
    out = summary.copy()
    patient = pd.read_csv(PATIENT_OUTCOMES, usecols=["new_patient_id", "gbtm_cluster"])
    cluster_n = (
        patient.groupby("gbtm_cluster")["new_patient_id"]
        .nunique()
        .rename("cluster_n")
        .reset_index()
        .rename(columns={"gbtm_cluster": "cluster"})
    )
    out = out.merge(cluster_n, on="cluster", how="left")

    audit_rows = []
    for prefix in ["sbp", "dbp"]:
        mean_col = f"{prefix}_mean"
        lo_col = f"{prefix}_ci_low"
        hi_col = f"{prefix}_ci_high"
        scaled_lo_col = f"{prefix}_ci_low_cluster_n_scaled"
        scaled_hi_col = f"{prefix}_ci_high_cluster_n_scaled"

        old_half_width = (out[hi_col] - out[lo_col]) / 2.0
        implied_sd = old_half_width * np.sqrt(out["bp_n"]) / MEAN_CI_Z
        new_half_width = MEAN_CI_Z * implied_sd / np.sqrt(out["cluster_n"])

        out[scaled_lo_col] = out[mean_col] - new_half_width
        out[scaled_hi_col] = out[mean_col] + new_half_width

        tmp = out[
            [
                "cluster",
                "quarter",
                "bp_n",
                "cluster_n",
                mean_col,
                lo_col,
                hi_col,
                scaled_lo_col,
                scaled_hi_col,
            ]
        ].copy()
        tmp["metric"] = prefix.upper()
        tmp["old_half_width_observed_bp_n"] = old_half_width
        tmp["new_half_width_cluster_n_scaled"] = new_half_width
        tmp["half_width_ratio_new_over_old"] = new_half_width / old_half_width
        audit_rows.append(tmp)

    audit = pd.concat(audit_rows, ignore_index=True)
    audit.to_csv(CI_AUDIT_CSV, index=False)
    return out


def plot_metric(
    summary: pd.DataFrame,
    *,
    value: str,
    lo: str,
    hi: str,
    ylabel: str,
    filename: str,
    percent: bool = False,
    reference_y: float | None = None,
    reference_label: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 6.6))
    factor = 100.0 if percent else 1.0

    all_values: list[pd.Series] = []
    for cluster, sub in summary.groupby("cluster", sort=True):
        sub = sub.sort_values("quarter")
        x = sub["quarter"].to_numpy(dtype=float)
        y = sub[value].to_numpy(dtype=float) * factor
        low = sub[lo].to_numpy(dtype=float) * factor
        high = sub[hi].to_numpy(dtype=float) * factor
        color = CLUSTER_COLORS[int(cluster)]

        ax.fill_between(x, low, high, color=color, alpha=RIBBON_ALPHA, linewidth=0, zorder=1)
        ax.plot(
            x,
            y,
            color=color,
            linewidth=LINE_WIDTH,
            solid_capstyle="round",
            solid_joinstyle="round",
            label=f"Cluster {cluster}",
            zorder=3,
        )
        ax.scatter(x, y, color=color, edgecolor="white", linewidth=0.8, s=45, zorder=4)
        all_values.extend([pd.Series(y), pd.Series(low), pd.Series(high)])

    if reference_y is not None:
        ax.axhline(reference_y, color="#555555", linestyle="--", linewidth=1.2, alpha=0.75, zorder=2)
        if reference_label:
            ax.text(
                12.15,
                reference_y,
                reference_label,
                ha="right",
                va="bottom",
                fontsize=8,
                color="#444444",
            )
        all_values.append(pd.Series([reference_y]))

    ax.set_xlabel("Quarters since first visit")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.8, 12.2)
    ax.set_xticks(range(1, 13))
    ax.grid(True, color="#e8e8e8")
    ax.set_axisbelow(True)

    vals = pd.concat(all_values, ignore_index=True).dropna()
    if len(vals):
        pad = (vals.max() - vals.min()) * 0.12
        if percent:
            ax.set_ylim(max(0.0, vals.min() - pad), min(100.0, vals.max() + pad))
        else:
            ax.set_ylim(vals.min() - pad, vals.max() + pad)

    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.09), frameon=False)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    for out_dir in [PAPER_DIR, PRESENTATION_DIR, OVERLEAF_FIG_DIR]:
        fig.savefig(out_dir / f"{filename}.png", dpi=300, bbox_inches="tight")
        fig.savefig(out_dir / f"{filename}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    stems = [
        "quarterly_treatment_rate_ng4",
        "quarterly_control_rate_ng4",
        "quarterly_sbp_ng4",
        "quarterly_dbp_ng4",
        "quarterly_retention_rate_ng4",
    ]
    backup_existing_outputs(stems)

    summary = add_cluster_size_scaled_bp_ci(pd.read_csv(SUMMARY_CSV))
    retention_summary = pd.read_csv(RETENTION_SUMMARY_CSV)

    plot_metric(
        retention_summary,
        value="retention_rate",
        lo="retention_ci_low",
        hi="retention_ci_high",
        ylabel="Retention rate (%)",
        filename="quarterly_retention_rate_ng4",
        percent=True,
    )

    plot_metric(
        summary,
        value="treatment_rate",
        lo="treatment_ci_low",
        hi="treatment_ci_high",
        ylabel="Treatment rate among observed visits (%)",
        filename="quarterly_treatment_rate_ng4",
        percent=True,
    )
    plot_metric(
        summary,
        value="control_rate",
        lo="control_ci_low",
        hi="control_ci_high",
        ylabel="BP control rate among BP-recorded visits (%)",
        filename="quarterly_control_rate_ng4",
        percent=True,
    )
    plot_metric(
        summary,
        value="sbp_mean",
        lo="sbp_ci_low_cluster_n_scaled",
        hi="sbp_ci_high_cluster_n_scaled",
        ylabel="Mean SBP among BP-recorded visits, mmHg",
        filename="quarterly_sbp_ng4",
    )
    plot_metric(
        summary,
        value="dbp_mean",
        lo="dbp_ci_low_cluster_n_scaled",
        hi="dbp_ci_high_cluster_n_scaled",
        ylabel="Mean DBP among BP-recorded visits, mmHg",
        filename="quarterly_dbp_ng4",
    )

    print(f"Regenerated quarterly clinical plots in {PAPER_DIR}")
    print(f"Backed up replaced active files in {BACKUP_DIR}")
    print(f"Saved SBP/DBP CI scaling audit to {CI_AUDIT_CSV}")


if __name__ == "__main__":
    main()
