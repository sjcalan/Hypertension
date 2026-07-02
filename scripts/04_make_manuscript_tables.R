parse_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

latex_escape <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("\\$", "\\\\$", x)
  x <- gsub("#", "\\\\#", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x
}

fmt_p <- function(p) {
  if (length(p) == 0L || is.na(p)) return("")
  if (p < 0.001) return("$<0.001$")
  sprintf("%.3f", p)
}

fmt_n <- function(n) {
  format(as.integer(n), big.mark = ",", scientific = FALSE)
}

fmt_n_pct <- function(x) {
  y <- parse_num(x)
  denom <- sum(!is.na(y))
  if (denom == 0L) return("NA")
  num <- sum(y >= 1, na.rm = TRUE)
  sprintf("%s (%.1f\\%%)", fmt_n(num), 100 * num / denom)
}

fmt_n_denom_pct <- function(x) {
  y <- parse_num(x)
  denom <- sum(!is.na(y))
  if (denom == 0L) return("NA")
  num <- sum(y >= 1, na.rm = TRUE)
  sprintf("%s/%s (%.1f\\%%)", fmt_n(num), fmt_n(denom), 100 * num / denom)
}

fmt_count_pct <- function(num, denom) {
  if (is.na(denom) || denom == 0L) return("NA")
  sprintf("%s (%.1f\\%%)", fmt_n(num), 100 * num / denom)
}

fmt_count_denom_pct <- function(num, denom) {
  if (is.na(denom) || denom == 0L) return("NA")
  sprintf("%s/%s (%.1f\\%%)", fmt_n(num), fmt_n(denom), 100 * num / denom)
}

fmt_mean_sd <- function(x, digits = 1L) {
  y <- parse_num(x)
  y <- y[!is.na(y)]
  if (!length(y)) return("NA")
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(y), sd(y))
}

fmt_median_iqr <- function(x, digits = 0L) {
  y <- parse_num(x)
  y <- y[!is.na(y)]
  if (!length(y)) return("NA")
  qs <- quantile(y, probs = c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    qs[[2]], qs[[1]], qs[[3]]
  )
}

fmt_date_iqr <- function(x) {
  d <- as.Date(x)
  d <- d[!is.na(d)]
  if (!length(d)) return("NA")
  qs <- quantile(as.numeric(d), probs = c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
  sprintf(
    "%s (%s-%s)",
    as.character(as.Date(round(qs[[2]]), origin = "1970-01-01")),
    as.character(as.Date(round(qs[[1]]), origin = "1970-01-01")),
    as.character(as.Date(round(qs[[3]]), origin = "1970-01-01"))
  )
}

p_binary <- function(df, col, group_col = "gbtm_cluster") {
  x <- parse_num(df[[col]])
  g <- df[[group_col]]
  keep <- !is.na(x) & !is.na(g)
  if (sum(keep) == 0L || length(unique(g[keep])) < 2L) return(NA_real_)
  tab <- table(factor(g[keep]), factor(ifelse(x[keep] >= 1, 1, 0), levels = c(1, 0)))
  if (any(rowSums(tab) == 0L) || ncol(tab) < 2L) return(NA_real_)
  suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
}

p_cont <- function(df, col, group_col = "gbtm_cluster") {
  y <- parse_num(df[[col]])
  g <- df[[group_col]]
  keep <- !is.na(y) & !is.na(g)
  if (sum(keep) == 0L || length(unique(g[keep])) < 2L) return(NA_real_)
  suppressWarnings(kruskal.test(y[keep] ~ as.factor(g[keep]))$p.value)
}

p_date <- function(df, col, group_col = "gbtm_cluster") {
  d <- as.numeric(as.Date(df[[col]]))
  g <- df[[group_col]]
  keep <- !is.na(d) & !is.na(g)
  if (sum(keep) == 0L || length(unique(g[keep])) < 2L) return(NA_real_)
  suppressWarnings(kruskal.test(d[keep] ~ as.factor(g[keep]))$p.value)
}

table_body <- function(rows) {
  out <- character()
  for (row in rows) {
    if (length(row) == 1L) {
      out <- c(out, row)
    } else {
      out <- c(out, paste0(paste(row, collapse = " & "), " \\\\"))
    }
  }
  out
}

write_table1 <- function(df, table_dir) {
  clusters <- 1:4
  df$female <- as.numeric(df$gender_clean == "Female")
  df$hispanic <- as.numeric(df$race_clean == "Hispanic/Latino")
  df$other_unknown <- as.numeric(df$race_clean == "Other/Unknown")
  df$english <- as.numeric(df$language_clean == "English")
  df$spanish <- as.numeric(df$language_clean == "Spanish")
  df$baseline_controlled <- parse_num(df$controlled)
  df$baseline_treated <- parse_num(df$treated)

  add_cont <- function(label, col, digits = 1L) {
    c(
      latex_escape(label),
      fmt_mean_sd(df[[col]], digits),
      sapply(clusters, function(k) fmt_mean_sd(df[df$gbtm_cluster == k, col], digits)),
      fmt_p(p_cont(df, col))
    )
  }
  add_bin <- function(label, col, denom = FALSE) {
    formatter <- if (denom) fmt_n_denom_pct else fmt_n_pct
    c(
      latex_escape(label),
      formatter(df[[col]]),
      sapply(clusters, function(k) formatter(df[df$gbtm_cluster == k, col])),
      fmt_p(p_binary(df, col))
    )
  }

  rows <- list(
    c(
      "Characteristic",
      sprintf("Overall\\newline N=%s", fmt_n(nrow(df))),
      sapply(clusters, function(k) sprintf("Cluster %d\\newline N=%s", k, fmt_n(sum(df$gbtm_cluster == k)))),
      "P"
    ),
    "\\midrule",
    add_cont("Age, mean (SD), y", "age_at_first_visit", 1L),
    add_bin("Female sex", "female"),
    add_bin("Hispanic/Latino ethnicity", "hispanic"),
    add_bin("Other/unknown race/ethnicity", "other_unknown"),
    add_bin("English language", "english"),
    add_bin("Spanish language", "spanish"),
    add_bin("Diabetes", "diabetes"),
    add_bin("Chronic kidney disease", "ckd"),
    add_cont("Baseline systolic BP, mean (SD), mmHg", "baseline_sbp", 1L),
    add_cont("Baseline diastolic BP, mean (SD), mmHg", "baseline_dbp", 1L),
    add_bin("Controlled BP at baseline, No./available (%)", "baseline_controlled", TRUE),
    add_bin("Treated at baseline", "baseline_treated")
  )
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Baseline characteristics of the modeled cohort by visit-trajectory group}",
    "\\label{tab:baseline}",
    "\\begingroup",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2pt}",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\begin{threeparttable}",
    "\\begin{tabularx}{\\textwidth}{p{0.24\\textwidth}>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}p{0.06\\textwidth}}",
    "\\toprule",
    table_body(rows),
    "\\bottomrule",
    "\\end{tabularx}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item Values are mean (SD), No. (\\%), or No./available (\\%). The denominator for baseline BP control is restricted to patients with nonmissing baseline SBP and DBP, which explains why the control denominator is smaller than the modeled cohort size. Treated at baseline indicates antihypertensive medication captured in the structured current-medication field at the first observed visit; it should not be interpreted as ever treated, prescribed outside the available EHR, or clinically eligible for treatment. P values compare the four trajectory groups using chi-square tests for categorical variables and Kruskal-Wallis tests for continuous variables. BP denotes blood pressure.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\endgroup",
    "\\end{table}"
  )
  writeLines(lines, file.path(table_dir, "table1_baseline_by_cluster.tex"))
}

compute_bp_metrics <- function(visits, assignments) {
  v <- merge(visits, assignments, by = "new_patient_id")
  v$systolic_bp <- parse_num(v$systolic_bp)
  v$diastolic_bp <- parse_num(v$diastolic_bp)
  v$controlled <- parse_num(v$controlled)
  v$visit_sequence <- parse_num(v$visit_sequence)
  v$encounter_dateD <- parse_num(v$encounter_dateD)
  bp <- v[!is.na(v$systolic_bp) & !is.na(v$diastolic_bp) & !is.na(v$controlled), ]
  bp <- bp[order(bp$new_patient_id, bp$encounter_dateD, bp$visit_sequence), ]
  ids <- unique(assignments$new_patient_id)
  out <- data.frame(
    new_patient_id = ids,
    bp_visit_count = 0L,
    control_streak3_any_calc = NA_real_
  )
  chunks <- split(bp$controlled, bp$new_patient_id)
  counts <- lengths(chunks)
  out$bp_visit_count <- as.integer(counts[as.character(out$new_patient_id)])
  out$bp_visit_count[is.na(out$bp_visit_count)] <- 0L
  out$control_streak3_any_calc <- sapply(out$new_patient_id, function(id) {
    vals <- chunks[[as.character(id)]]
    if (is.null(vals) || length(vals) < 3L) return(NA_real_)
    as.numeric(any(stats::filter(vals == 1, rep(1, 3), sides = 1) == 3, na.rm = TRUE))
  })
  out
}

write_table2 <- function(outcomes, bp_metrics, table_dir) {
  clusters <- 1:4
  df <- merge(outcomes, bp_metrics, by = "new_patient_id", all.x = TRUE)
  df$eligible_bp3 <- as.numeric(df$bp_visit_count >= 3L)
  df$control_streak3_any_calc <- parse_num(df$control_streak3_any_calc)

  add_bin <- function(label, col, denom = FALSE) {
    formatter <- if (denom) fmt_n_denom_pct else fmt_n_pct
    c(
      latex_escape(label),
      sapply(clusters, function(k) formatter(df[df$gbtm_cluster == k, col])),
      fmt_p(p_binary(df, col))
    )
  }
  add_cont <- function(label, col, digits = 0L, data = df) {
    c(
      latex_escape(label),
      sapply(clusters, function(k) fmt_median_iqr(data[data$gbtm_cluster == k, col], digits)),
      fmt_p(p_cont(data, col))
    )
  }
  add_section <- function(label) {
    sprintf("\\addlinespace \\multicolumn{6}{l}{\\textit{%s}} \\\\[-0.3em]", label)
  }

  eligible <- df[df$eligible_bp3 == 1, ]
  rows <- list(
    c("Measure", "C1", "C2", "C3", "C4", "P"),
    "\\midrule",
    "\\multicolumn{6}{l}{\\textit{Retention, No. (\\%)}} \\\\[-0.3em]",
    add_bin("12 months", "retained_12mo"),
    add_bin("24 months", "retained_24mo"),
    add_bin("36 months", "retained_36mo"),
    add_section("Care engagement, median (IQR)"),
    add_cont("Total visits", "total_visits", 0L),
    add_cont("Time to last visit, days", "time_to_last_visit", 0L),
    add_cont("Median visit gap, days", "visit_duration_median", 0L),
    add_section("Fixed-window observation, No. (\\%)"),
    add_bin("Observed visit near 12 months", "any_visit_observed_12mo"),
    add_bin("Observed visit near 24 months", "any_visit_observed_24mo"),
    add_bin("Observed visit near 36 months", "any_visit_observed_36mo"),
    add_bin("BP-recorded visit near 12 months", "bp_visit_observed_12mo"),
    add_bin("BP-recorded visit near 24 months", "bp_visit_observed_24mo"),
    add_bin("BP-recorded visit near 36 months", "bp_visit_observed_36mo"),
    add_section("Fixed-window treatment, No./available (\\%)"),
    add_bin("Treatment near 12 months", "treated_12mo", TRUE),
    add_bin("Treatment near 24 months", "treated_24mo", TRUE),
    add_bin("Treatment near 36 months", "treated_36mo", TRUE),
    add_section("Fixed-window BP control, No./available (\\%)"),
    add_bin("BP control near 12 months", "control_12mo_visitlevel", TRUE),
    add_bin("BP control near 24 months", "control_24mo_visitlevel", TRUE),
    add_bin("BP control near 36 months", "control_36mo_visitlevel", TRUE),
    add_section("Longitudinal BP-control outcome and measurement opportunity"),
    add_bin("Eligible for durable BP control (>=3 BP-recorded visits)", "eligible_bp3"),
    add_cont("BP-recorded visits among eligible patients", "bp_visit_count", 0L, eligible),
    add_bin("Any 3 consecutive controlled BP visits among eligible", "control_streak3_any_calc", TRUE)
  )
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Care-engagement and clinical outcomes by visit-trajectory group}",
    "\\label{tab:cluster_comparison}",
    "\\begingroup",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2pt}",
    "\\renewcommand{\\arraystretch}{1.03}",
    "\\begin{threeparttable}",
    "\\begin{tabularx}{\\textwidth}{p{0.34\\textwidth}>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}p{0.06\\textwidth}}",
    "\\toprule",
    table_body(rows),
    "\\bottomrule",
    "\\end{tabularx}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item C1, long-span sparse engagement; C2, sustained moderate engagement; C3, short-span intensive early engagement; C4, frequent intermediate engagement. Retention is defined as the final observed visit occurring on or after the indicated time since the patient's first observed visit. Fixed-window treatment uses the closest observed visit within $\\pm45$ days of the target date; fixed-window BP and BP-control outcomes use the closest BP-recorded visit within $\\pm45$ days. Treatment and BP-control cells are shown as No./available (\\%) to make the changing complete-case denominators explicit. Durable BP control was evaluated only among patients with at least 3 BP-recorded visits, and BP measurement counts are shown because the outcome is opportunity-dependent. Binary p values use chi-square tests; continuous p values use Kruskal-Wallis tests. BP denotes blood pressure.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\endgroup",
    "\\end{table}"
  )
  writeLines(lines, file.path(table_dir, "table2_cluster_comparison_corrected.tex"))

  display_rows <- rows[sapply(rows, length) == 6L]
  display <- as.data.frame(do.call(rbind, display_rows), stringsAsFactors = FALSE)
  names(display) <- display[1, ]
  display <- display[-1, , drop = FALSE]
  display[] <- lapply(display, function(x) {
    x <- gsub("\\\\%", "%", x)
    x <- gsub("\\$<0.001\\$", "<0.001", x)
    x
  })
  write.csv(display, file.path(data_dir, "table2_cluster_comparison_corrected.csv"), row.names = FALSE)

  invisible(df)
}

write_excluded_table <- function(baseline, assignments, table_dir) {
  baseline$modeled <- as.numeric(baseline$new_patient_id %in% assignments$new_patient_id)
  baseline$female <- as.numeric(baseline$gender_clean == "Female")
  baseline$hispanic <- as.numeric(baseline$race_clean == "Hispanic/Latino")
  baseline$other_unknown <- as.numeric(baseline$race_clean == "Other/Unknown")
  baseline$english <- as.numeric(baseline$language_clean == "English")
  baseline$spanish <- as.numeric(baseline$language_clean == "Spanish")
  baseline$baseline_controlled <- parse_num(baseline$controlled)
  baseline$baseline_treated <- parse_num(baseline$treated)
  baseline$visit_count_num <- parse_num(baseline$visit_count)
  baseline$one_two_visits <- as.numeric(baseline$visit_count_num <= 2L)

  modeled <- baseline[baseline$modeled == 1, ]
  excluded <- baseline[baseline$modeled == 0, ]
  p2bin <- function(col) fmt_p(p_binary(transform(baseline, group = modeled), col, "group"))
  p2cont <- function(col) fmt_p(p_cont(transform(baseline, group = modeled), col, "group"))
  add_bin <- function(label, col, denom = FALSE) {
    formatter <- if (denom) fmt_n_denom_pct else fmt_n_pct
    c(latex_escape(label), formatter(excluded[[col]]), formatter(modeled[[col]]), p2bin(col))
  }
  add_cont <- function(label, col, digits = 1L, median = FALSE) {
    formatter <- if (median) fmt_median_iqr else fmt_mean_sd
    c(latex_escape(label), formatter(excluded[[col]], digits), formatter(modeled[[col]], digits), p2cont(col))
  }

  rows <- list(
    c("Characteristic", sprintf("Not modeled\\newline N=%s", fmt_n(nrow(excluded))), sprintf("Modeled\\newline N=%s", fmt_n(nrow(modeled))), "P"),
    "\\midrule",
    add_cont("Total visits, median (IQR)", "visit_count_num", 0L, TRUE),
    add_bin("One or two total visits", "one_two_visits"),
    add_cont("Age, mean (SD), y", "age_at_first_visit", 1L),
    add_bin("Female sex", "female"),
    add_bin("Hispanic/Latino ethnicity", "hispanic"),
    add_bin("Other/unknown race/ethnicity", "other_unknown"),
    add_bin("English language", "english"),
    add_bin("Spanish language", "spanish"),
    add_bin("Diabetes", "diabetes"),
    add_bin("Chronic kidney disease", "ckd"),
    add_cont("Baseline systolic BP, mean (SD), mmHg", "baseline_sbp", 1L),
    add_cont("Baseline diastolic BP, mean (SD), mmHg", "baseline_dbp", 1L),
    add_bin("Controlled BP at baseline, No./available (%)", "baseline_controlled", TRUE),
    add_bin("Treated at baseline", "baseline_treated")
  )
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Patients not included in the final trajectory model compared with the modeled cohort}",
    "\\label{tab:excluded_modeled}",
    "\\begingroup",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\begin{threeparttable}",
    "\\begin{tabularx}{\\textwidth}{p{0.42\\textwidth}>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}p{0.08\\textwidth}}",
    "\\toprule",
    table_body(rows),
    "\\bottomrule",
    "\\end{tabularx}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    sprintf("\\item The not-modeled group includes all %s patients from the full visit-level cohort who were not assigned to the final 4-group model. Of these, %s had one or two total visits in the baseline cohort file; the remainder had at least 3 visits overall but did not contribute sufficient relative-time panel information for the final modeled dataset. Values are mean (SD), median (IQR), No. (\\%%), or No./available (\\%%). BP denotes blood pressure.", fmt_n(nrow(excluded)), fmt_n(sum(excluded$visit_count_num <= 2L, na.rm = TRUE))),
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\endgroup",
    "\\end{table}"
  )
  writeLines(lines, file.path(table_dir, "table_s2_excluded_vs_modeled.tex"))
}

write_admin_table <- function(baseline, assignments, visits, table_dir) {
  df <- merge(assignments, baseline[, c("new_patient_id", "first_visit_date")], by = "new_patient_id", all.x = TRUE)
  admin_end <- max(as.Date(visits$encounter_date), na.rm = TRUE)
  df$first_visit_date <- as.Date(df$first_visit_date)
  df$possible_followup_days <- as.numeric(admin_end - df$first_visit_date)
  df$observable_24mo <- as.numeric(df$possible_followup_days >= 720)
  df$observable_36mo <- as.numeric(df$possible_followup_days >= 1080)
  clusters <- 1:4
  rows <- list(
    c("Measure", "C1", "C2", "C3", "C4", "P"),
    "\\midrule",
    c("First observed visit date, median (IQR)", sapply(clusters, function(k) fmt_date_iqr(df[df$gbtm_cluster == k, "first_visit_date"])), fmt_p(p_date(df, "first_visit_date"))),
    c("Maximum administratively observable follow-up, days, median (IQR)", sapply(clusters, function(k) fmt_median_iqr(df[df$gbtm_cluster == k, "possible_followup_days"], 0L)), fmt_p(p_cont(df, "possible_followup_days"))),
    c("Administratively observable for >=24 months", sapply(clusters, function(k) fmt_n_pct(df[df$gbtm_cluster == k, "observable_24mo"])), fmt_p(p_binary(df, "observable_24mo"))),
    c("Administratively observable for >=36 months", sapply(clusters, function(k) fmt_n_pct(df[df$gbtm_cluster == k, "observable_36mo"])), fmt_p(p_binary(df, "observable_36mo")))
  )
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    sprintf("\\caption{Calendar-time entry and administrative follow-up opportunity by visit-trajectory group through the data cutoff (%s)}", as.character(admin_end)),
    "\\label{tab:administrative_followup}",
    "\\begingroup",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\begin{threeparttable}",
    "\\begin{tabularx}{\\textwidth}{p{0.36\\textwidth}>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}p{0.07\\textwidth}}",
    "\\toprule",
    table_body(rows),
    "\\bottomrule",
    "\\end{tabularx}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item Maximum administratively observable follow-up is the interval from each patient's first observed visit to the last encounter date available in the analytic visit-level file. This table assesses calendar-time truncation that can remain even after relative-time alignment. P values are descriptive and compare trajectory groups.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\endgroup",
    "\\end{table}"
  )
  writeLines(lines, file.path(table_dir, "table_s3_administrative_followup_by_cluster.tex"))
  invisible(df)
}

write_bp_opportunity_table <- function(bp_df, table_dir) {
  clusters <- 1:4
  bp_df$has_bp <- as.numeric(bp_df$bp_visit_count >= 1L)
  bp_df$eligible_bp3 <- as.numeric(bp_df$bp_visit_count >= 3L)
  eligible <- bp_df[bp_df$eligible_bp3 == 1L, ]
  add_bin <- function(label, col, data = bp_df, denom = FALSE) {
    formatter <- if (denom) fmt_n_denom_pct else fmt_n_pct
    c(latex_escape(label), sapply(clusters, function(k) formatter(data[data$gbtm_cluster == k, col])), fmt_p(p_binary(data, col)))
  }
  add_cont <- function(label, col, data = bp_df) {
    c(latex_escape(label), sapply(clusters, function(k) fmt_median_iqr(data[data$gbtm_cluster == k, col], 0L)), fmt_p(p_cont(data, col)))
  }
  rows <- list(
    c("Measure", "C1", "C2", "C3", "C4", "P"),
    "\\midrule",
    add_bin("At least 1 BP-recorded visit", "has_bp"),
    add_bin("At least 3 BP-recorded visits", "eligible_bp3"),
    add_cont("BP-recorded visits among all patients", "bp_visit_count"),
    add_cont("BP-recorded visits among patients with >=3 BP visits", "bp_visit_count", eligible),
    add_bin("Any 3 consecutive controlled BP visits among patients with >=3 BP visits", "control_streak3_any_calc", eligible, TRUE)
  )
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{BP measurement opportunity for longitudinal BP-control analyses}",
    "\\label{tab:bp_opportunity}",
    "\\begingroup",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\begin{threeparttable}",
    "\\begin{tabularx}{\\textwidth}{p{0.38\\textwidth}>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}X>{\\centering\\arraybackslash}p{0.07\\textwidth}}",
    "\\toprule",
    table_body(rows),
    "\\bottomrule",
    "\\end{tabularx}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item BP-recorded visits require nonmissing systolic and diastolic BP values. The durable-control outcome is conditional on having at least 3 BP-recorded visits and is therefore sensitive to differential measurement opportunity across trajectory groups.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\endgroup",
    "\\end{table}"
  )
  writeLines(lines, file.path(table_dir, "table_s4_bp_observation_opportunity.tex"))
}

write_diagnostic_status_table <- function(table_dir) {
  rows <- list(
    c("Diagnostic", "Current status", "Needed for final reporting"),
    "\\midrule",
    c("Average posterior probability (APP)", "Not recoverable from the current hard-assignment export", "Rerun the primary GBTM after exporting posterior membership probabilities"),
    c("Odds of correct classification (OCC)", "Not recoverable from the current hard-assignment export", "Posterior probabilities plus assigned class proportions"),
    c("Relative entropy", "The exported model-fit file contains a missing entropy field", "Saved model object or package output that retains entropy"),
    c("Prior/posterior probability ratio", "Not recoverable from the current hard-assignment export", "Posterior class probabilities and assigned class proportions")
  )
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Status of GBTM classification diagnostics in the current exported analysis files}",
    "\\label{tab:gbtm_diagnostics_status}",
    "\\begingroup",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{1.12}",
    "\\begin{threeparttable}",
    "\\begin{tabularx}{\\textwidth}{p{0.26\\textwidth}p{0.34\\textwidth}X}",
    "\\toprule",
    table_body(rows),
    "\\bottomrule",
    "\\end{tabularx}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    "\\item The trajectory script has been revised to save the fitted model object and posterior-probability-derived diagnostics in future runs. These diagnostics cannot be reconstructed from the currently saved assignment file because it contains only hard class labels.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\endgroup",
    "\\end{table}"
  )
  writeLines(lines, file.path(table_dir, "table_s5_gbtm_diagnostic_status.tex"))
}

base_dir <- "/userhome/cs3/u3011656/hypertension/hypertension"
out_dir <- file.path(base_dir, "analysis", "overleaf_hypertension_visit_trajectory_2026-05-25")
table_dir <- file.path(out_dir, "tables")
data_dir <- file.path(out_dir, "data")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

assignments <- read.csv(file.path(base_dir, "analysis/outputs/q1_gbtm_relative_visit_only_biweekly_full_ng4_gpu_minvis3_parallel_2026-04-17/gbtm_relative_visit_only_assignments.csv"), stringsAsFactors = FALSE)
baseline <- read.csv(file.path(base_dir, "hypertension_cohort_with_baseline.csv"), stringsAsFactors = FALSE)
outcomes <- read.csv(file.path(base_dir, "analysis/outputs/q1_ng4_cluster_presentation_2026-05-03/q1_ng4_cluster_patient_outcomes.csv"), stringsAsFactors = FALSE)
visits <- read.csv(file.path(base_dir, "analysis/outputs/visit_level/visit_level_cohort.csv"), stringsAsFactors = FALSE)

baseline_model <- merge(assignments, baseline, by = "new_patient_id", all.x = TRUE)
bp_metrics <- compute_bp_metrics(visits, assignments)
bp_with_cluster <- merge(assignments, bp_metrics, by = "new_patient_id", all.x = TRUE)

write_table1(baseline_model, table_dir)
table2_df <- write_table2(outcomes, bp_metrics, table_dir)
write_excluded_table(baseline, assignments, table_dir)
admin_df <- write_admin_table(baseline, assignments, visits, table_dir)
write_bp_opportunity_table(bp_with_cluster, table_dir)
write_diagnostic_status_table(table_dir)

write.csv(table2_df, file.path(data_dir, "table2_cluster_comparison_revision_source.csv"), row.names = FALSE)
write.csv(bp_with_cluster, file.path(data_dir, "bp_observation_opportunity_by_patient.csv"), row.names = FALSE)
write.csv(admin_df, file.path(data_dir, "administrative_followup_by_patient.csv"), row.names = FALSE)
