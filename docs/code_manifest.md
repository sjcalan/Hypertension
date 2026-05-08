# Code Manifest

This manifest records the main scripts from the original working directory that belong to the **non follow-up vs follow-up** part of the project.

## Main entry-point scripts

- `/userhome/cs3/u3011656/hypertension/hypertension/step1_patient_level_clean.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/step2_types_clean.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/step2_type4_type5_clean.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/table1_baseline_characteristics.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/table2_retention_patterns.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/figure1_visit_distribution.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/figure2_blood_pressure_analysis.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/figure3_patient_journey.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/figure4_kaplan_meier.py`
- `/userhome/cs3/u3011656/hypertension/hypertension/figure5_forest_plot.py`

## Expected mapping into this repo

- `step1_patient_level_clean.py`
  - `src/data_prep/01_patient_selection/`
- diagnosis-derived comorbidity logic from the early cohort build
  - `src/data_prep/02_comorbidity_flags/`
- visit-level merging and cohort assembly logic
  - `src/cohort_build/01_visit_level_merge/`
- BP, BMI, treatment, and control enrichment logic
  - `src/cohort_build/02_bp_bmi_enrichment/`
- follow-up group creation logic
  - `src/cohort_build/03_followup_groups/`
- `table1_baseline_characteristics.py`
  - `src/tables/table1_baseline/`
- `table2_retention_patterns.py`
  - `src/tables/table2_retention/`
- `figure1_visit_distribution.py`
  - `src/figures/figure1_visit_distribution/`
- `figure2_blood_pressure_analysis.py`
  - `src/figures/figure2_blood_pressure/`
- `figure3_patient_journey.py`
  - `src/figures/figure3_patient_journey/`
- `figure4_kaplan_meier.py`
  - `src/figures/figure4_kaplan_meier/`
- `figure5_forest_plot.py`
  - `src/figures/figure5_forest_plot/`
- survival / hazard-style follow-up analyses
  - `src/survival/`

## Important note

The code itself has **not** yet been copied into this repo structure.
This file only records what belongs here when migration starts.
