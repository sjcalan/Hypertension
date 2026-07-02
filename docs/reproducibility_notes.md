# Reproducibility Notes

This repository contains the analysis code but not the confidential data required to execute the full workflow.

## Local Paths

Some scripts contain absolute paths from the secure HKU analysis environment. Before rerunning the code elsewhere, update the base path constants or pass command-line arguments where supported.

## Data Requirements

The workflow expects local source or derived files containing:

- encounter diagnoses
- vitals, including systolic and diastolic blood pressure
- medication exposure information
- demographics and social history
- derived visit-level and patient-level cohort files
- GBTM model outputs and posterior cluster assignments for downstream tables and plots

## Output Policy

Generated tables, figures, model outputs, and manuscript exports are excluded from GitHub because they may contain derived patient-level information or large generated artifacts.

For manuscript submission, only de-identified aggregate tables and figures should be exported into the journal-specific supplement.
