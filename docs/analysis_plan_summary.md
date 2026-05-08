# Analysis Plan Summary

This summary is based on the original analysis plan in:

- `/userhome/cs3/u3011656/hypertension/hypertension/Analysis Plan.docx`

and on the completed project work tracked in:

- `/userhome/cs3/u3011656/hypertension/hypertension/analysis/docs/project_status_report_2026-04-24.md`

## First-Part Workflow Covered by This Repo

### 1. Hypertension cohort selection

- start from `lds_encounter_diagnoses`
- identify patients with diagnosis strings containing `hypertensi*`
- exclude records with missing `new_patient_id`

### 2. Patient-level enrichment

- merge demographics
- add comorbidity flags derived from diagnosis text:
  - diabetes
  - heart failure
  - chronic kidney disease
  - stroke
  - heart attack

### 3. Visit-level enrichment

- merge with `lds_vitals`
- merge visit-level diagnosis and medication context
- include core vitals such as:
  - systolic BP
  - diastolic BP
  - pulse
  - respiration rate
  - temperature
  - height
  - weight
- calculate BMI where possible

### 4. Visit-level derived indicators

- create `controlled` at each visit from systolic and diastolic BP
- create `treated` at each visit from medication exposure timing

### 5. Follow-up grouping

The first analysis block focused on comparing patients by follow-up pattern, especially:

- no follow-up
- at least one follow-up

The related scripts also generated more detailed visit-count groupings in the working directory.

### 6. First-phase outputs

The completed first-phase scripts produced:

- visit-number summary outputs
- baseline comparison tables
- retention-pattern tables
- blood-pressure-related figures
- patient journey plots
- Kaplan-Meier style follow-up analyses
- forest plot output

## Why This Repo Exists

This repo is a clean structure for that first analysis block only. It is not meant to package the later trajectory-modeling phase yet.
