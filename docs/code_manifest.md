# Code Manifest

This manifest lists the source files currently included for collaborator review. Data files and generated outputs are intentionally excluded.

## Main Workflow Scripts

- `scripts/01_build_visit_level_cohort.py`
  - constructs the visit-level hypertension cohort
- `scripts/02_fit_visit_only_gbtm.R`
  - fits visit-only group-based trajectory models
- `scripts/03_compare_gbtm_models.py`
  - compares candidate trajectory model solutions
- `scripts/04_make_manuscript_tables.R`
  - generates manuscript-facing cluster comparison tables
- `scripts/05_make_quarterly_trajectory_plots.py`
  - generates quarterly retention, treatment, control, SBP, and DBP plots
- `scripts/06_make_individual_journey_plots.py`
  - generates individual patient journey plots by trajectory group
- `scripts/07_prepare_overleaf_assets.py`
  - prepares selected tables and figures for Overleaf
- `scripts/08_create_manuscript_docx.py`
  - exports a collaborator-review manuscript draft

## Cluster Job Scripts

- `jobs/run_gbtm_ng4.slurm`
  - SLURM batch job for the primary ng=4 model run

## Supporting Scripts

- `scripts/supporting/explore_ng4_cluster_differences.py`
  - exploratory cluster-difference summaries and earlier plot variants
- `scripts/supporting/full_denominator_bp_sensitivity.py`
  - sensitivity plots using full cluster denominators for BP-related outcomes

## Excluded Files

- raw EHR tables
- patient-level analytic datasets
- intermediate model objects
- generated CSV tables
- generated PDF/PNG figures
- manuscript DOCX/ZIP exports
