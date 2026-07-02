# Hypertension Visit-Trajectory Analysis

This repository contains the analysis code for a hypertension visit-trajectory study using longitudinal electronic health record data. The project identifies visit-based care-retention trajectories with group-based trajectory modeling and compares the resulting groups on retention, treatment, blood pressure control, SBP, and DBP outcomes.

## Repository Structure

- `scripts/`
  - main analysis scripts
- `scripts/supporting/`
  - exploratory and sensitivity-analysis scripts
- `jobs/`
  - SLURM job scripts used for longer model-fitting runs
- `docs/`
  - workflow notes and code manifest
- `data/`
  - data availability and local-data setup note
- `results/`
  - results

## Main Workflow

1. Build the visit-level hypertension cohort.
2. Fit visit-only group-based trajectory models.
3. Compare candidate trajectory models.
4. Generate manuscript-facing cluster comparison tables.
5. Generate quarterly clinical and retention trajectory plots.
6. Generate individual patient journey plots.
7. Prepare Overleaf manuscript assets.
8. Export a collaborator-review manuscript draft.

See `scripts/README.md` for the script-by-script map.

## Data And Privacy

The study uses sensitive EHR-derived data. These data are not included in this repository.

## Current Primary Model

The primary manuscript analysis uses a 4-group visit-only GBTM solution among patients with at least 3 visits. Visits are aligned to each patient's first observed visit and represented using 14-day relative-time bins over a 36-month window.
