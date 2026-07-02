# Code Manifest

This repository contains the current source scripts for the hypertension visit-trajectory analysis. Data files and generated outputs are intentionally excluded.

## Active Scripts

- `src/analysis_scripts/06_build_visit_level_cohort.py`
- `src/analysis_scripts/07i_q1_gbtm_relative_visit_only_shorttest.R`
- `src/analysis_scripts/run_q1_gbtm_relative_visit_only_full_ng4_gpu_minvis3_parallel.slurm`
- `src/analysis_scripts/10_compare_gbtm_ng_models.py`
- `src/analysis_scripts/13_ng4_explore_cluster_differences.py`
- `src/analysis_scripts/14_prepare_overleaf_manuscript_assets.py`
- `src/analysis_scripts/15_make_compact_journey_appendix.py`
- `src/analysis_scripts/16_prepare_overleaf_revision_tables.R`
- `src/analysis_scripts/17_create_google_doc_docx.py`
- `src/analysis_scripts/17_ng4_full_denominator_bp_sensitivity.py`
- `src/analysis_scripts/18_regenerate_quarterly_clinical_plots.py`

## Coverage

- Cohort construction
- Visit-only GBTM model fitting
- SLURM job configuration for the ng=4 model run
- Model comparison across candidate group counts
- Cluster-level table generation
- Quarterly clinical and retention trajectory plotting
- Individual journey plotting
- Manuscript table and figure export
- BP denominator sensitivity checks
