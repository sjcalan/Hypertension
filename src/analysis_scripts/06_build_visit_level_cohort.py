from __future__ import annotations

from collections import defaultdict
from heapq import heappop, heappush
from pathlib import Path
import json
import numpy as np
import pandas as pd


ORIGIN_DATE = "1960-01-01"


def pick_cohort_file(base_dir: Path) -> Path:
    candidates = [
        base_dir / "hypertension_cohort.csv",
        base_dir / "hypertension_cohort_with_baseline.csv",
        base_dir / "hypertension_cohort_with_metrics.csv",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("No hypertension cohort file found in project root.")


def load_patient_profile(base_dir: Path, cohort_path: Path) -> tuple[set[int], pd.DataFrame]:
    wanted = [
        "new_patient_id",
        "pt_dob",
        "pt_ethnicity",
        "pt_gender",
        "pt_language",
        "pt_city",
        "pt_state",
        "pt_race",
        "height_cm",
        "weight_kg",
        "bmi",
    ]
    cols = pd.read_csv(cohort_path, nrows=0).columns.tolist()
    usecols = [c for c in wanted if c in cols]
    cohort = pd.read_csv(cohort_path, usecols=usecols, low_memory=False)
    cohort["new_patient_id"] = normalize_patient_id(cohort["new_patient_id"])
    cohort = cohort[cohort["new_patient_id"].notna()].copy()
    cohort["new_patient_id"] = cohort["new_patient_id"].astype(int)
    cohort_ids = set(cohort["new_patient_id"].tolist())

    # Ensure required demographic columns exist even if source file is missing any.
    for c in ["pt_dob", "pt_ethnicity", "pt_gender", "pt_language", "pt_city", "pt_state", "pt_race"]:
        if c not in cohort.columns:
            cohort[c] = np.nan

    rename_map = {}
    if "height_cm" in cohort.columns:
        rename_map["height_cm"] = "baseline_height_cm"
    if "weight_kg" in cohort.columns:
        rename_map["weight_kg"] = "baseline_weight_kg"
    if "bmi" in cohort.columns:
        rename_map["bmi"] = "baseline_bmi"
    cohort = cohort.rename(columns=rename_map)

    # One patient-level row from cohort.
    cohort = cohort.sort_values("new_patient_id").drop_duplicates("new_patient_id", keep="first")

    # Pull demographics directly from lds_demographics (primary source), then fallback to cohort.
    demo_path = base_dir / "lds_demographics.csv"
    if demo_path.exists():
        demo_cols = ["new_patient_id", "pt_dob", "pt_ethnicity", "pt_gender", "pt_language", "pt_city", "pt_state", "pt_race"]
        d = pd.read_csv(demo_path, usecols=demo_cols, low_memory=False)
        d["new_patient_id"] = normalize_patient_id(d["new_patient_id"])
        d = d[d["new_patient_id"].isin(cohort_ids)].copy()
        d["new_patient_id"] = d["new_patient_id"].astype(int)
        d = d.sort_values("new_patient_id").drop_duplicates("new_patient_id", keep="first")

        merged = cohort.merge(d, on="new_patient_id", how="left", suffixes=("_cohort", "_lds"))
        for c in ["pt_dob", "pt_ethnicity", "pt_gender", "pt_language", "pt_city", "pt_state", "pt_race"]:
            merged[c] = merged[f"{c}_lds"].combine_first(merged[f"{c}_cohort"])
        drop_cols = [f"{c}_cohort" for c in ["pt_dob", "pt_ethnicity", "pt_gender", "pt_language", "pt_city", "pt_state", "pt_race"]]
        drop_cols += [f"{c}_lds" for c in ["pt_dob", "pt_ethnicity", "pt_gender", "pt_language", "pt_city", "pt_state", "pt_race"]]
        cohort = merged.drop(columns=drop_cols)

    return cohort_ids, cohort


def normalize_patient_id(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def parse_daynum(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def daynum_to_date(daynum: pd.Series) -> pd.Series:
    return pd.to_datetime(daynum, unit="D", origin=ORIGIN_DATE, errors="coerce")


def build_visit_base(base_dir: Path, cohort_ids: set[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    enc_path = base_dir / "lds_encounter_diagnoses.csv"
    chunks = []
    usecols = ["new_patient_id", "encounter_dateD", "dx_name"]

    for chunk in pd.read_csv(enc_path, usecols=usecols, chunksize=400_000, low_memory=False):
        chunk["new_patient_id"] = normalize_patient_id(chunk["new_patient_id"])
        chunk = chunk[chunk["new_patient_id"].isin(cohort_ids)]
        if chunk.empty:
            continue
        chunk["encounter_dateD"] = parse_daynum(chunk["encounter_dateD"])
        chunk = chunk[chunk["encounter_dateD"].between(0, 30000, inclusive="both")]
        if chunk.empty:
            continue
        chunk["dx_name"] = chunk["dx_name"].fillna("").astype(str)
        chunks.append(chunk)

    if not chunks:
        raise RuntimeError("No encounter rows found for selected hypertension cohort.")

    enc = pd.concat(chunks, ignore_index=True)
    enc["new_patient_id"] = enc["new_patient_id"].astype(int)
    enc["encounter_dateD"] = enc["encounter_dateD"].astype(int)
    enc["dx_upper"] = enc["dx_name"].str.upper()

    # Visit-level diagnosis flags.
    enc["dx_hypertension"] = enc["dx_upper"].str.contains("HYPERTENSI", na=False)
    enc["dx_diabetes"] = enc["dx_upper"].str.contains("DIABET", na=False)
    enc["dx_heart_failure"] = enc["dx_upper"].str.contains(r"HEART FAILURE|CHF", regex=True, na=False)
    enc["dx_ckd"] = enc["dx_upper"].str.contains(r"CHRONIC KIDNEY DISEASE| CKD|^CKD|CKD ", regex=True, na=False)
    enc["dx_stroke"] = enc["dx_upper"].str.contains(r"STROKE|CVA|CEREBROVASCULAR ACCIDENT", regex=True, na=False)
    enc["dx_heart_attack"] = enc["dx_upper"].str.contains(r"HEART ATTACK|MYOCARDIAL INFARCTION", regex=True, na=False)

    gcols = ["new_patient_id", "encounter_dateD"]
    visit_flags = (
        enc.groupby(gcols)[
            ["dx_hypertension", "dx_diabetes", "dx_heart_failure", "dx_ckd", "dx_stroke", "dx_heart_attack"]
        ]
        .max()
        .reset_index()
    )
    visit_flags = visit_flags.rename(
        columns={
            "dx_hypertension": "hypertension_dx_visit",
            "dx_diabetes": "diabetes_dx_visit",
            "dx_heart_failure": "heart_failure_dx_visit",
            "dx_ckd": "ckd_dx_visit",
            "dx_stroke": "stroke_dx_visit",
            "dx_heart_attack": "heart_attack_dx_visit",
        }
    )

    visit_dx = (
        enc.loc[enc["dx_name"] != "", gcols + ["dx_name"]]
        .drop_duplicates()
        .groupby(gcols)["dx_name"]
        .agg(lambda s: " | ".join(s.astype(str).head(40)))
        .reset_index()
        .rename(columns={"dx_name": "dx_names_visit"})
    )

    visits = visit_flags.merge(visit_dx, on=gcols, how="left")

    # Patient-level comorbidity flags derived from any visit.
    any_cols = [
        "diabetes_dx_visit",
        "heart_failure_dx_visit",
        "ckd_dx_visit",
        "stroke_dx_visit",
        "heart_attack_dx_visit",
    ]
    patient_any = visits.groupby("new_patient_id")[any_cols].max().reset_index()
    patient_any = patient_any.rename(columns={c: c.replace("_visit", "_any") for c in any_cols})

    visits["encounter_date"] = daynum_to_date(visits["encounter_dateD"])
    visits = visits.sort_values(["new_patient_id", "encounter_dateD"]).reset_index(drop=True)
    visits["visit_sequence"] = visits.groupby("new_patient_id").cumcount() + 1
    return visits, patient_any


def build_vitals_wide(base_dir: Path, cohort_ids: set[int]) -> pd.DataFrame:
    vit_path = base_dir / "lds_vitals.csv"
    keep_map = {
        "SYSTOLIC BP": "systolic_bp",
        "DIASTOLIC BP": "diastolic_bp",
        "RESPIRATION RATE": "respiration_rate",
        "PULSE": "pulse",
        "WEIGHT": "weight",
        "HEIGHT": "height",
        "TEMPERATURE": "temperature",
    }

    rows = []
    usecols = ["aa_patient_id", "record_dated", "measure_name", "measure_value"]
    for chunk in pd.read_csv(vit_path, usecols=usecols, chunksize=500_000, low_memory=False):
        chunk["new_patient_id"] = normalize_patient_id(chunk["aa_patient_id"])
        chunk = chunk[chunk["new_patient_id"].isin(cohort_ids)]
        if chunk.empty:
            continue
        chunk["encounter_dateD"] = parse_daynum(chunk["record_dated"])
        chunk = chunk[chunk["encounter_dateD"].between(0, 30000, inclusive="both")]
        if chunk.empty:
            continue
        upper_name = chunk["measure_name"].astype(str).str.upper().str.strip()
        chunk["measure_key"] = upper_name.map(keep_map)
        chunk = chunk[chunk["measure_key"].notna()]
        if chunk.empty:
            continue
        chunk["measure_value"] = pd.to_numeric(chunk["measure_value"], errors="coerce")
        chunk = chunk[chunk["measure_value"].notna()]
        if chunk.empty:
            continue
        rows.append(chunk[["new_patient_id", "encounter_dateD", "measure_key", "measure_value"]])

    if not rows:
        return pd.DataFrame(columns=["new_patient_id", "encounter_dateD"])

    v = pd.concat(rows, ignore_index=True)
    v["new_patient_id"] = v["new_patient_id"].astype(int)
    v["encounter_dateD"] = v["encounter_dateD"].astype(int)

    wide = (
        v.pivot_table(
            index=["new_patient_id", "encounter_dateD"],
            columns="measure_key",
            values="measure_value",
            aggfunc="median",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    return wide


def build_med_intervals(base_dir: Path, cohort_ids: set[int]) -> dict[int, list[tuple[int, int]]]:
    med_path = base_dir / "lds_current_medications.csv"
    usecols = ["new_patient_id", "encounter_dateD", "med_startdate", "med_stopdate", "med_name"]
    intervals: dict[int, list[tuple[int, int]]] = defaultdict(list)

    for chunk in pd.read_csv(med_path, usecols=usecols, chunksize=350_000, low_memory=False):
        chunk["new_patient_id"] = normalize_patient_id(chunk["new_patient_id"])
        chunk = chunk[chunk["new_patient_id"].isin(cohort_ids)]
        if chunk.empty:
            continue

        med_name = chunk["med_name"].fillna("").astype(str).str.upper()
        chunk = chunk[~med_name.str.contains("NO CURRENT MEDICATIONS", na=False)]
        if chunk.empty:
            continue

        enc_d = parse_daynum(chunk["encounter_dateD"])
        start_date = pd.to_datetime(chunk["med_startdate"], errors="coerce", utc=True).dt.tz_localize(None)
        stop_date = pd.to_datetime(chunk["med_stopdate"], errors="coerce", utc=True).dt.tz_localize(None)

        start_d = ((start_date - pd.Timestamp(ORIGIN_DATE)).dt.days).astype("float64")
        stop_d = ((stop_date - pd.Timestamp(ORIGIN_DATE)).dt.days).astype("float64")

        # Fallbacks if explicit interval missing.
        start_d = start_d.fillna(enc_d)
        stop_d = stop_d.fillna(start_d)

        valid = start_d.notna() & stop_d.notna()
        c = chunk.loc[valid, ["new_patient_id"]].copy()
        c["start_d"] = start_d[valid].astype(int)
        c["stop_d"] = stop_d[valid].astype(int)

        swap_mask = c["stop_d"] < c["start_d"]
        if swap_mask.any():
            svals = c.loc[swap_mask, "start_d"].copy()
            c.loc[swap_mask, "start_d"] = c.loc[swap_mask, "stop_d"].values
            c.loc[swap_mask, "stop_d"] = svals.values

        for pid, s, e in c.itertuples(index=False):
            intervals[int(pid)].append((int(s), int(e)))

    for pid in list(intervals.keys()):
        intervals[pid].sort(key=lambda x: (x[0], x[1]))
    return intervals


def assign_treatment(visits: pd.DataFrame, med_intervals: dict[int, list[tuple[int, int]]]) -> pd.DataFrame:
    visits = visits.copy()
    visits["treated"] = 0
    visits["active_med_count"] = 0

    for pid, idx in visits.groupby("new_patient_id").groups.items():
        intervals = med_intervals.get(int(pid), [])
        if not intervals:
            continue

        days = visits.loc[idx, "encounter_dateD"].to_numpy(dtype=int)
        order = np.argsort(days)
        sorted_days = days[order]
        treated = np.zeros(len(sorted_days), dtype=np.int8)
        counts = np.zeros(len(sorted_days), dtype=np.int16)

        ptr = 0
        active: list[int] = []
        for i, d in enumerate(sorted_days):
            while ptr < len(intervals) and intervals[ptr][0] <= d:
                heappush(active, intervals[ptr][1])
                ptr += 1
            while active and active[0] < d:
                heappop(active)
            if active:
                treated[i] = 1
                counts[i] = len(active)

        back = np.empty_like(order)
        back[order] = np.arange(len(order))
        visits.loc[idx, "treated"] = treated[back]
        visits.loc[idx, "active_med_count"] = counts[back]

    return visits


def add_bmi_and_controlled(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["weight"] = pd.to_numeric(out.get("weight"), errors="coerce")
    out["height"] = pd.to_numeric(out.get("height"), errors="coerce")

    # LDS vitals are mostly lbs/in; keep robust conversion for mixed unit encodings.
    out["weight_kg_visit"] = np.where(out["weight"] > 0, out["weight"] * 0.45359237, np.nan)
    out["height_m_visit"] = np.select(
        [
            out["height"].between(36, 96, inclusive="both"),   # likely inches
            out["height"].between(100, 250, inclusive="both"), # likely centimeters
        ],
        [
            out["height"] * 0.0254,
            out["height"] / 100.0,
        ],
        default=np.nan,
    )

    # Carry forward/backward within patient for visit-derived anthropometrics.
    out = out.sort_values(["new_patient_id", "encounter_dateD"]).reset_index(drop=True)
    out["weight_kg"] = out.groupby("new_patient_id")["weight_kg_visit"].transform(lambda s: s.ffill().bfill())
    out["height_m"] = out.groupby("new_patient_id")["height_m_visit"].transform(lambda s: s.ffill().bfill())

    # Fallback to patient-level baseline values if still missing.
    if "baseline_weight_kg" in out.columns:
        out["baseline_weight_kg"] = pd.to_numeric(out["baseline_weight_kg"], errors="coerce")
        out["weight_kg"] = out["weight_kg"].fillna(out["baseline_weight_kg"])
    if "baseline_height_cm" in out.columns:
        out["baseline_height_cm"] = pd.to_numeric(out["baseline_height_cm"], errors="coerce")
        out["height_m"] = out["height_m"].fillna(out["baseline_height_cm"] / 100.0)

    out["bmi"] = np.where(
        (out["height_m"] > 0) & out["weight_kg"].notna(),
        out["weight_kg"] / (out["height_m"] ** 2),
        np.nan,
    )
    if "baseline_bmi" in out.columns:
        out["baseline_bmi"] = pd.to_numeric(out["baseline_bmi"], errors="coerce")
        out["bmi"] = out["bmi"].fillna(out["baseline_bmi"])

    # Keep only one final anthropometric pair for output:
    # weight in kg and height in meters.
    out["weight"] = out["weight_kg"]
    out["height"] = out["height_m"]

    sbp = pd.to_numeric(out["systolic_bp"], errors="coerce")
    dbp = pd.to_numeric(out["diastolic_bp"], errors="coerce")
    out["controlled"] = np.where(
        sbp.notna() & dbp.notna(),
        np.where((sbp >= 140) | (dbp >= 90), 0, 1),
        np.nan,
    )
    return out


def main() -> None:
    base_dir = Path(r"G:\hypertension")
    analysis_dir = base_dir / "analysis"
    out_dir = analysis_dir / "outputs" / "visit_level"
    out_dir.mkdir(parents=True, exist_ok=True)

    cohort_path = pick_cohort_file(base_dir)
    cohort_ids, patient_profile = load_patient_profile(base_dir, cohort_path)

    visits, patient_any = build_visit_base(base_dir, cohort_ids)
    vitals_wide = build_vitals_wide(base_dir, cohort_ids)
    merged = visits.merge(vitals_wide, on=["new_patient_id", "encounter_dateD"], how="left")
    merged = merged.merge(patient_any, on="new_patient_id", how="left")
    merged = merged.merge(patient_profile, on="new_patient_id", how="left")

    for c in ["pt_ethnicity", "pt_gender", "pt_language", "pt_city", "pt_state", "pt_race"]:
        if c in merged.columns:
            merged[c] = merged[c].astype(str).replace({".": np.nan, "": np.nan, "nan": np.nan}).fillna("Unknown")

    med_intervals = build_med_intervals(base_dir, cohort_ids)
    merged = assign_treatment(merged, med_intervals)
    merged = add_bmi_and_controlled(merged)

    # Column order for easy downstream use.
    col_order = [
        "new_patient_id",
        "visit_sequence",
        "encounter_dateD",
        "encounter_date",
        "pt_dob",
        "pt_ethnicity",
        "pt_gender",
        "pt_language",
        "pt_city",
        "pt_state",
        "pt_race",
        "dx_names_visit",
        "hypertension_dx_visit",
        "diabetes_dx_visit",
        "heart_failure_dx_visit",
        "ckd_dx_visit",
        "stroke_dx_visit",
        "heart_attack_dx_visit",
        "diabetes_dx_any",
        "heart_failure_dx_any",
        "ckd_dx_any",
        "stroke_dx_any",
        "heart_attack_dx_any",
        "systolic_bp",
        "diastolic_bp",
        "pulse",
        "respiration_rate",
        "temperature",
        "weight",
        "height",
        "bmi",
        "controlled",
        "treated",
        "active_med_count",
    ]
    existing_cols = [c for c in col_order if c in merged.columns]
    merged = merged[existing_cols].sort_values(["new_patient_id", "encounter_dateD"]).reset_index(drop=True)

    out_main = base_dir / "visit_level_cohort.csv"
    out_copy = out_dir / "visit_level_cohort.csv"
    merged.to_csv(out_main, index=False)
    merged.to_csv(out_copy, index=False)

    summary = {
        "cohort_file_used": str(cohort_path),
        "n_unique_patients_input": int(len(cohort_ids)),
        "n_visit_rows_output": int(len(merged)),
        "n_unique_patients_output": int(merged["new_patient_id"].nunique()),
        "n_with_vitals_sbp_dbp": int((merged["systolic_bp"].notna() & merged["diastolic_bp"].notna()).sum()),
        "controlled_non_missing": int(merged["controlled"].notna().sum()),
        "treated_positive_rows": int((merged["treated"] == 1).sum()),
        "output_paths": [str(out_main), str(out_copy)],
    }
    with open(out_dir / "visit_level_summary.json", "w", encoding="ascii") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
