### 1. Hypertension cohort construction

- start from encounter diagnosis records
- identify patients with diagnosis strings containing hypertension-related text
- exclude records with missing patient identifiers
- merge demographics, social history, vitals, diagnosis context, and medication information
- construct a visit-level analytic file with one row per patient visit

### 2. Visit-level derived indicators

- derive systolic and diastolic blood pressure variables
- define visit-level blood pressure control using the planned threshold:
  - uncontrolled if systolic BP is at least 140 mmHg or diastolic BP is at least 90 mmHg
  - controlled if both systolic and diastolic BP are observed and below those thresholds
- derive treatment indicators from medication exposure information
- derive visit timing, visit number, and patient-level follow-up summaries

### 3. Visit-only trajectory construction

- align every patient to their own first observed visit as day 0
- use relative follow-up time rather than absolute calendar date
- divide follow-up into 14-day bins over the modeled 36-month window
- code visit activity within the observed span between the first and last visit
- avoid padding time after a patient's last observed follow-up with artificial no-visit records

### 4. Group-based trajectory modeling

- fit visit-only group-based trajectory models using binned visit activity
- exclude patients with only 2 visits from the primary trajectory analysis
- compare candidate models with 3, 4, 5, and 6 groups
- select the 4-group model as the primary solution based on fit, class balance, interpretability, and clinical coherence

### 5. Cluster comparison

- compare the 4 visit-trajectory groups on retention, visit spacing, treatment, blood pressure, and blood pressure control outcomes
- report overall cluster-comparison p values and pairwise cluster comparisons where needed
- evaluate alternative blood pressure control definitions, including last observed BP control and durable control across consecutive controlled visits

### 6. Manuscript-facing outputs

- generate the main cluster-comparison table
- generate quarterly retention, treatment, BP control, SBP, and DBP trajectory plots
- generate compact individual patient journey plots by trajectory group
- prepare selected tables and figures for the Overleaf manuscript package
