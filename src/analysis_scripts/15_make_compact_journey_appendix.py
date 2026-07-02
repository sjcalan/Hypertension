from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache-codex")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import pandas as pd


BASE_DIR = Path("/userhome/cs3/u3011656/hypertension/hypertension")
ASSIGN_PATH = (
    BASE_DIR
    / "analysis/outputs/q1_gbtm_relative_visit_only_biweekly_full_ng4_gpu_minvis3_parallel_2026-04-17/"
    / "gbtm_relative_visit_only_assignments.csv"
)
VISIT_PATH_CANDIDATES = [
    BASE_DIR / "visit_level_cohort.csv",
    BASE_DIR.parent / "_cleanup_archive_2026-06-30/hypertension/visit_level_cohort.csv",
    BASE_DIR.parent / "_cleanup_archive_2026-06-30/visit_level_cohort.csv",
]
OUT_DIR = BASE_DIR / "analysis/overleaf_hypertension_visit_trajectory_2026-05-25/figures"
PAPER_ASSET_DIR = BASE_DIR / "analysis/outputs/q1_ng4_cluster_paper_assets_2026-05-15"
OUT_FILE = OUT_DIR / "appendix_compact_journey_2x2_lines_only.png"
MAX_DAYS = 1095
CONTROLLED_COLOR = "#0077BB"
UNCONTROLLED_COLOR = "#EE7733"
CHANGED_COLOR = "#777777"


def find_visit_path() -> Path:
    for path in VISIT_PATH_CANDIDATES:
        if path.exists():
            return path
    candidates = "\n".join(str(path) for path in VISIT_PATH_CANDIDATES)
    raise FileNotFoundError(f"Could not find visit_level_cohort.csv. Checked:\n{candidates}")


def add_cluster_journey(
    ax: plt.Axes,
    bp_visits: pd.DataFrame,
    assignments: pd.DataFrame,
    cluster: int,
    line_width: float = 0.45,
    alpha: float = 0.72,
) -> None:
    cluster_all_n = int((assignments["gbtm_cluster"] == cluster).sum())
    sub = bp_visits.loc[bp_visits["gbtm_cluster"] == cluster].copy()
    sub = sub.loc[sub["days_since_first_visit"].between(0, MAX_DAYS, inclusive="both")].copy()
    sub = sub.sort_values(["new_patient_id", "days_since_first_visit", "visit_sequence"])

    patient_counts = (
        sub.groupby("new_patient_id", as_index=False)
        .agg(
            total_bp_visits=("visit_sequence", "size"),
            max_days=("days_since_first_visit", "max"),
        )
        # Order patients by observed BP follow-up length, longest to shortest.
        # Sorting by visit count first creates repeated "peaks" in the right edge.
        .sort_values(["max_days", "total_bp_visits", "new_patient_id"], ascending=[False, False, True])
        .reset_index(drop=True)
    )
    bp_patient_count = len(patient_counts)
    patient_counts = patient_counts.loc[patient_counts["total_bp_visits"] >= 2].reset_index(drop=True)
    patient_order = patient_counts["new_patient_id"].tolist()
    y_lookup = {pid: i for i, pid in enumerate(patient_order)}
    sub = sub.loc[sub["new_patient_id"].isin(y_lookup)].copy()

    cc_segments = []
    uu_segments = []
    changed_segments = []
    initial_carried_segments = []

    for pid, person in sub.groupby("new_patient_id", sort=False):
        person = person.sort_values(["days_since_first_visit", "visit_sequence"])
        y = y_lookup[pid]
        times = person["days_since_first_visit"].to_numpy()
        ctrl = person["controlled"].astype(int).to_numpy()
        if len(times) > 1 and times[0] > 0:
            # The journey plot is anchored at each patient's first overall visit.
            # If the first BP-recorded encounter happened later, carry its first
            # observed BP status back to day 0 for display so visible journeys do
            # not begin with a blank lead-in.
            initial_carried_segments.append(([(0, y), (times[0], y)], ctrl[0]))
        for j in range(len(times) - 1):
            segment = [(times[j], y), (times[j + 1], y)]
            if ctrl[j] == 1 and ctrl[j + 1] == 1:
                cc_segments.append(segment)
            elif ctrl[j] == 0 and ctrl[j + 1] == 0:
                uu_segments.append(segment)
            else:
                changed_segments.append(segment)

    for segment, first_status in initial_carried_segments:
        color = CONTROLLED_COLOR if first_status == 1 else UNCONTROLLED_COLOR
        ax.add_collection(LineCollection([segment], colors=color, linewidths=line_width, alpha=alpha, zorder=1))
    for segments, color in [
        (cc_segments, CONTROLLED_COLOR),
        (uu_segments, UNCONTROLLED_COLOR),
        (changed_segments, CHANGED_COLOR),
    ]:
        if segments:
            ax.add_collection(LineCollection(segments, colors=color, linewidths=line_width, alpha=alpha, zorder=1))

    n_patients = len(patient_order)
    n_visits = int(len(sub))
    control_rate = sub["controlled"].astype(int).mean() * 100 if n_visits else 0
    ax.set_title(
        f"Cluster {cluster}: N={cluster_all_n:,}; line patients={n_patients:,}\n"
        f"BP patients={bp_patient_count:,}; controlled={control_rate:.1f}%",
        fontsize=9,
        fontweight="bold",
    )
    ax.set_xlim(0, MAX_DAYS)
    ax.set_ylim(max(n_patients - 0.5, 0.5), -0.5)
    ax.set_box_aspect(1)
    if n_patients <= 1:
        ax.set_yticks([0])
        ax.set_yticklabels(["1"])
    else:
        mid = (n_patients - 1) / 2
        ax.set_yticks([0, mid, n_patients - 1])
        ax.set_yticklabels(["1", f"{int(round(mid + 1)):,}", f"{n_patients:,}"], fontsize=8)
    ax.set_xticks([0, 200, 400, 600, 800, 1000])
    ax.set_xticklabels(["0", "200", "400", "600", "800", "1000"], fontsize=8)
    ax.grid(axis="x", color="#dddddd", linewidth=0.5)
    ax.set_xlabel("Time (days)", fontsize=8)
    ax.set_ylabel("Patient rank", fontsize=8)


def save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PAPER_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    assignments = pd.read_csv(ASSIGN_PATH, usecols=["new_patient_id", "gbtm_cluster"], low_memory=False)
    assignments["new_patient_id"] = pd.to_numeric(assignments["new_patient_id"], errors="coerce")
    assignments["gbtm_cluster"] = pd.to_numeric(assignments["gbtm_cluster"], errors="coerce")
    assignments = assignments.dropna().astype({"new_patient_id": int, "gbtm_cluster": int})

    visit_path = find_visit_path()
    visit = pd.read_csv(
        visit_path,
        usecols=["new_patient_id", "visit_sequence", "encounter_dateD", "systolic_bp", "diastolic_bp", "controlled"],
        low_memory=False,
    )
    for col in ["new_patient_id", "visit_sequence", "encounter_dateD", "systolic_bp", "diastolic_bp", "controlled"]:
        visit[col] = pd.to_numeric(visit[col], errors="coerce")
    visit = visit.dropna(subset=["new_patient_id", "visit_sequence", "encounter_dateD"]).copy()
    visit["new_patient_id"] = visit["new_patient_id"].astype(int)
    visit["visit_sequence"] = visit["visit_sequence"].astype(int)
    visit = visit.merge(assignments, on="new_patient_id", how="inner")
    first_day = visit.groupby("new_patient_id", as_index=False)["encounter_dateD"].min().rename(columns={"encounter_dateD": "first_day"})
    visit = visit.merge(first_day, on="new_patient_id", how="left")
    visit["days_since_first_visit"] = visit["encounter_dateD"] - visit["first_day"]

    bp_visits = visit.loc[
        visit["systolic_bp"].notna() & visit["diastolic_bp"].notna() & visit["controlled"].notna()
    ].copy()

    for cluster in [1, 2, 3, 4]:
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        add_cluster_journey(ax, bp_visits, assignments, cluster, line_width=0.55, alpha=0.80)
        fig.subplots_adjust(left=0.11, right=0.98, top=0.82, bottom=0.20)
        legend_elements = [
            Line2D([0], [0], color=CONTROLLED_COLOR, linewidth=2, label="Controlled to controlled"),
            Line2D([0], [0], color=UNCONTROLLED_COLOR, linewidth=2, label="Uncontrolled to uncontrolled"),
            Line2D([0], [0], color=CHANGED_COLOR, linewidth=2, label="Status changed"),
        ]
        fig.legend(handles=legend_elements, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.02))
        individual_path = OUT_DIR / f"figure3_cluster{cluster}_journey_time_lines_only.png"
        save_figure(fig, individual_path)
        save_figure(fig, PAPER_ASSET_DIR / individual_path.name)
        plt.close(fig)
        print(f"Saved {individual_path}")

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 10.5))
    for ax, cluster in zip(axes.flat, [1, 2, 3, 4]):
        add_cluster_journey(ax, bp_visits, assignments, cluster)

    legend_elements = [
        Line2D([0], [0], color=CONTROLLED_COLOR, linewidth=2, label="Controlled to controlled"),
        Line2D([0], [0], color=UNCONTROLLED_COLOR, linewidth=2, label="Uncontrolled to uncontrolled"),
        Line2D([0], [0], color=CHANGED_COLOR, linewidth=2, label="Status changed"),
    ]
    fig.subplots_adjust(left=0.06, right=0.98, top=0.97, bottom=0.09, hspace=0.26, wspace=0.12)
    fig.legend(handles=legend_elements, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.01))
    save_figure(fig, OUT_FILE)
    save_figure(fig, PAPER_ASSET_DIR / OUT_FILE.name)
    plt.close(fig)
    print(f"Saved {OUT_FILE}")
    print(f"Read visit-level data from {visit_path}")


if __name__ == "__main__":
    main()
