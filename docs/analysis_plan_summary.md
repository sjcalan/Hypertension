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

- no follow-up
- at least one follow-up

### 6. First-phase outputs

- visit-number summary outputs
- baseline comparison tables
- retention-pattern tables
- blood-pressure-related figures
- patient journey plots
- Kaplan-Meier style follow-up analyses
- forest plot output
