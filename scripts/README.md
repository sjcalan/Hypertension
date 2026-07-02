# Scripts

Run the scripts in numeric order for the main analysis workflow.

1. `01_build_visit_level_cohort.py`
   - builds the visit-level hypertension cohort from diagnosis, vitals, demographics, and medication sources
2. `02_fit_visit_only_gbtm.R`
   - fits visit-only group-based trajectory models using relative-time 14-day visit bins
3. `03_compare_gbtm_models.py`
   - compares candidate trajectory solutions across group counts
4. `04_make_manuscript_tables.R`
   - creates manuscript-facing cluster comparison tables
5. `05_make_quarterly_trajectory_plots.py`
   - generates quarterly retention, treatment, BP control, SBP, and DBP trajectory plots
6. `06_make_individual_journey_plots.py`
   - generates individual patient journey plots by trajectory group
7. `07_prepare_overleaf_assets.py`
   - prepares figures and tables for the Overleaf manuscript package
8. `08_create_manuscript_docx.py`
   - creates a manuscript draft export for collaborator review

The `supporting/` folder contains exploratory and sensitivity-analysis scripts that informed the final manuscript figures but are not part of the shortest main workflow.
