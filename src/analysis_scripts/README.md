# Current Analysis Scripts

This folder contains the active source scripts used for the hypertension visit-trajectory project.

## Workflow

1. Build the visit-level hypertension cohort.
2. Fit visit-only group-based trajectory models.
3. Compare candidate `ng` solutions.
4. Generate ng=4 cluster summaries, tables, and figures.
5. Prepare manuscript-facing Overleaf and Google Doc assets.
6. Run sensitivity checks for BP denominator choices and clinical trajectories.

## Script Map

- `06_build_visit_level_cohort.py`: Builds the visit-level cohort from diagnosis, vitals, medication, and demographic sources.
- `07i_q1_gbtm_relative_visit_only_shorttest.R`: Fits visit-only GBTM models using relative-time binned visit activity.
- `run_q1_gbtm_relative_visit_only_full_ng4_gpu_minvis3_parallel.slurm`: SLURM job script for the ng=4 GBTM run.
- `10_compare_gbtm_ng_models.py`: Compares GBTM model outputs across candidate numbers of groups.
- `13_ng4_explore_cluster_differences.py`: Generates ng=4 cluster summaries, exploratory clinical comparisons, and quarterly trajectories.
- `14_prepare_overleaf_manuscript_assets.py`: Prepares selected manuscript tables and assets for Overleaf.
- `15_make_compact_journey_appendix.py`: Generates compact individual journey plots by cluster.
- `16_prepare_overleaf_revision_tables.R`: Builds revised manuscript tables for the Overleaf version.
- `17_create_google_doc_docx.py`: Creates a Google-Doc-style manuscript export.
- `17_ng4_full_denominator_bp_sensitivity.py`: Generates BP outcome sensitivity plots using full cluster denominators.
- `18_regenerate_quarterly_clinical_plots.py`: Regenerates the final quarterly clinical plots using straight dot-to-dot lines.

## Not Included

Raw EHR data, patient-level analytic datasets, generated tables, and generated figures are intentionally excluded.
