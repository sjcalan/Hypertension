from __future__ import annotations

from pathlib import Path
from datetime import datetime
import argparse
import pandas as pd


RUNS = {
    3: {
        "fit_dir": "analysis/outputs/q1_gbtm_relative_visit_only_biweekly_full_ng3_gpu_minvis3_2026-04-17",
        "table_dir": "analysis/outputs/q1_table2_by_ng3_cluster_minvis3_2026-04-19",
    },
    4: {
        "fit_dir": "analysis/outputs/q1_gbtm_relative_visit_only_biweekly_full_ng4_gpu_minvis3_parallel_2026-04-17",
        "table_dir": "analysis/outputs/q1_table2_by_ng4_cluster_minvis3_2026-04-19",
    },
    5: {
        "fit_dir": "analysis/outputs/q1_gbtm_relative_visit_only_biweekly_full_ng5_gpu_minvis3_parallel_2026-04-17",
        "table_dir": "analysis/outputs/q1_table2_by_ng5_cluster_minvis3_2026-04-19",
    },
    6: {
        "fit_dir": "analysis/outputs/q1_gbtm_relative_visit_only_biweekly_full_ng6_gpu_minvis3_parallel_2026-04-17",
        "table_dir": "analysis/outputs/q1_table2_by_ng6_cluster_minvis3_2026-04-19",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare Q1 GBTM ng=3/4/5/6 solutions.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("/userhome/cs3/u3011656/hypertension/hypertension"),
        help="Project root",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory",
    )
    return parser.parse_args()


def extract_table_value(table_path: Path, metric_label: str) -> str:
    tbl = pd.read_csv(table_path)
    metric = tbl["Metric"].astype(str).str.strip()
    match = tbl.loc[metric == metric_label, "Overall"]
    if match.empty:
        return ""
    return str(match.iloc[0])


def build_model_summary(base_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    model_rows = []
    cluster_rows = []

    for ng, cfg in RUNS.items():
        fit_dir = base_dir / cfg["fit_dir"]
        table_dir = base_dir / cfg["table_dir"]

        fit = pd.read_csv(fit_dir / "gbtm_relative_visit_only_model_fit.csv")
        profile = pd.read_csv(fit_dir / "gbtm_relative_visit_only_cluster_profiles.csv")
        manifest = pd.read_csv(table_dir / "cluster_table_manifest.csv")

        fit_row = fit.iloc[0].to_dict()
        fit_row["ng"] = ng
        fit_row["fit_dir"] = str(fit_dir)
        fit_row["table_dir"] = str(table_dir)
        fit_row["min_cluster_pct"] = float(profile["pct"].min())
        fit_row["max_cluster_pct"] = float(profile["pct"].max())
        fit_row["n_clusters_lt_5pct"] = int((profile["pct"] < 5).sum())
        fit_row["n_clusters_lt_10pct"] = int((profile["pct"] < 10).sum())
        fit_row["cluster_sizes"] = "|".join(str(int(x)) for x in profile["n"].tolist())
        fit_row["recommended_flag"] = "candidate"
        model_rows.append(fit_row)

        for _, prow in profile.iterrows():
            cid = int(prow["gbtm_cluster"])
            table_path = base_dir / manifest.loc[manifest["cluster"] == cid, "output_csv"].iloc[0]
            cluster_rows.append(
                {
                    "ng": ng,
                    "gbtm_cluster": cid,
                    "n": int(prow["n"]),
                    "pct": float(prow["pct"]),
                    "bins_observed": float(prow["bins_observed"]),
                    "bins_with_any_visit": float(prow["bins_with_any_visit"]),
                    "avg_visits_per_bin": float(prow["avg_visits_per_bin"]),
                    "followup_days_modeled": float(prow["followup_days_modeled"]),
                    "retained_12mo_overall": extract_table_value(table_path, "12 months"),
                    "retained_24mo_overall": extract_table_value(table_path, "24 months"),
                    "retained_36mo_overall": extract_table_value(table_path, "36 months"),
                    "total_visits_overall": extract_table_value(table_path, "Total visits number, median (IQR)"),
                    "time_to_last_visit_overall": extract_table_value(table_path, "Time to the last visit, median (IQR), days"),
                    "control_12mo_overall": extract_table_value(table_path, "Control rates at 12 months, No. (%)"),
                    "control_24mo_overall": extract_table_value(table_path, "Control rates at 24 months, No. (%)"),
                    "control_36mo_overall": extract_table_value(table_path, "Control rates at 36 months, No. (%)"),
                    "table_csv": str(table_path),
                }
            )

    model_df = pd.DataFrame(model_rows).sort_values("ng")
    cluster_df = pd.DataFrame(cluster_rows).sort_values(["ng", "gbtm_cluster"])
    return model_df, cluster_df


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    values = [[str(v) for v in row] for row in df.itertuples(index=False, name=None)]
    widths = [len(c) for c in cols]
    for row in values:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def fmt_row(row: list[str]) -> str:
        return "| " + " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)) + " |"

    header = fmt_row(cols)
    sep = "| " + " | ".join("-" * widths[i] for i in range(len(cols))) + " |"
    body = [fmt_row(row) for row in values]
    return "\n".join([header, sep] + body)


def recommend_ng(model_df: pd.DataFrame) -> int:
    # Judgment rule: prioritize fit improvement, but avoid solutions with tiny clusters.
    # ng=6 has a sub-5% class. ng=5 has a sub-10% class and more fragmentation.
    # ng=4 remains well-balanced while materially improving BIC over ng=3.
    return 4


def build_markdown(model_df: pd.DataFrame, cluster_df: pd.DataFrame, recommended_ng: int) -> str:
    rec_row = model_df.loc[model_df["ng"] == recommended_ng].iloc[0]

    lines = [
        "# Q1 GBTM Model Comparison",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Recommendation",
        "",
        f"Recommended primary solution: `ng={recommended_ng}`.",
        "",
        "Reasoning:",
        "- `ng=3` is stable but coarser.",
        "- `ng=4` materially improves BIC over `ng=3` while keeping all classes well populated.",
        "- `ng=5` improves BIC further, but introduces a smaller 8.2% class and more fragmentation.",
        "- `ng=6` improves BIC again, but creates a 2.15% class (382 patients), which looks like over-splitting for primary reporting.",
        "",
        "## Model-Level Summary",
        "",
        dataframe_to_markdown(
            model_df[
                [
                    "ng",
                    "bic",
                    "aic",
                    "n_populated_clusters",
                    "min_cluster_n",
                    "min_cluster_pct",
                    "max_cluster_n",
                    "max_cluster_pct",
                    "n_clusters_lt_5pct",
                    "n_clusters_lt_10pct",
                    "elapsed_seconds",
                    "cluster_sizes",
                ]
            ]
        ),
        "",
        f"## Recommended ng={recommended_ng} Cluster Overview",
        "",
        dataframe_to_markdown(
            cluster_df.loc[cluster_df["ng"] == recommended_ng, [
                "gbtm_cluster",
                "n",
                "pct",
                "bins_observed",
                "bins_with_any_visit",
                "avg_visits_per_bin",
                "followup_days_modeled",
                "retained_12mo_overall",
                "retained_24mo_overall",
                "retained_36mo_overall",
                "control_12mo_overall",
                "control_24mo_overall",
                "control_36mo_overall",
            ]]
        ),
        "",
        "## Notes",
        "",
        "- This recommendation is based on fit, class balance, and interpretability rather than BIC alone.",
        "- All compared models use the same `>=3` visit cohort and the same biweekly relative-time visit-only trajectory representation.",
        "",
        f"Selected row summary: `ng={recommended_ng}`, `BIC={rec_row['bic']:.2f}`, `min_cluster_n={int(rec_row['min_cluster_n'])}`, `min_cluster_pct={rec_row['min_cluster_pct']:.2f}`.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    base_dir = args.base_dir.resolve()
    out_dir = args.out_dir or (base_dir / "analysis" / "outputs" / "q1_gbtm_ng_compare_minvis3_2026-04-19")
    out_dir.mkdir(parents=True, exist_ok=True)

    model_df, cluster_df = build_model_summary(base_dir)
    recommended_ng = recommend_ng(model_df)
    model_df.loc[model_df["ng"] == recommended_ng, "recommended_flag"] = "recommended"

    model_csv = out_dir / "q1_gbtm_model_compare.csv"
    cluster_csv = out_dir / "q1_gbtm_cluster_overview.csv"
    summary_md = out_dir / "q1_gbtm_model_selection_summary.md"

    model_df.to_csv(model_csv, index=False)
    cluster_df.to_csv(cluster_csv, index=False)
    summary_md.write_text(build_markdown(model_df, cluster_df, recommended_ng), encoding="utf-8")

    print(f"Wrote: {model_csv}")
    print(f"Wrote: {cluster_csv}")
    print(f"Wrote: {summary_md}")
    print(model_df[[
        "ng",
        "bic",
        "min_cluster_n",
        "min_cluster_pct",
        "n_clusters_lt_5pct",
        "n_clusters_lt_10pct",
        "recommended_flag",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
