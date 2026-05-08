# Hypertension

This repository is a **structure-first export** for the **first part** of the hypertension project:

- comparison of patients with **no follow-up**
- comparison of patients with **at least one follow-up visit**

It is intentionally limited to the earlier cohort-building, descriptive, retention, and outcome-comparison workflow. The later latent-trajectory work such as GBTM, LCMM, and BCClong is **not** the focus of this repo.

## What This Repo Is For

This GitHub repo is meant to give a clean, viewable project skeleton for:

- brief reporting
- collaboration
- gradual code migration from the working directory

The current source-of-truth scripts still live in the main project workspace, and this repo currently provides the folder structure that matches the first-phase analysis plan.

## Current Structure

- `docs/`
  - project notes, scope, analysis-plan summary, migration notes, code manifest
- `src/`
  - empty source-code structure for cohort construction, tables, figures, and survival analyses
- `data_reference/`
  - placeholders for raw-source and derived-cohort references
- `outputs_placeholder/`
  - empty output structure for tables, figures, and survival/model summaries

## Workflow Reflected Here

Based on the original analysis plan and the completed first-phase work, this repo structure is organized around:

1. selecting the hypertension cohort
2. merging demographics, diagnoses, vitals, and medication data
3. defining per-visit `controlled` and `treated` status
4. splitting patients by follow-up pattern
5. generating summary tables and figures
6. summarizing retention and outcome differences

## Important Notes

- This repository currently contains **structure only**
- code cleanup and copying can be done later
- empty folders are made visible with placeholder `README.md` files so the structure can be reviewed on GitHub

See:

- `docs/analysis_scope.md`
- `docs/analysis_plan_summary.md`
- `docs/code_manifest.md`
- `docs/migration_notes.md`
