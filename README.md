# Hypertension

This repository contains the source-code structure for the hypertension visit-trajectory analysis. It is intended for collaborators to review the analysis workflow and reproduce the code organization without exposing raw EHR data or patient-level derived files.

## Current Structure

- `src/analysis_scripts/`
  - current scripts for cohort construction, GBTM fitting, model comparison, tables, figures, manuscript assets, and sensitivity checks
- `docs/`
  - workflow summary and code manifest
- `data_reference/`
  - placeholders describing expected raw and derived data locations
- `outputs_placeholder/`
  - placeholder folders for generated tables, figures, and model outputs

## Data Policy

Raw EHR files, patient-level analytic datasets, generated CSV outputs, and manuscript figures are not committed to this repository.
