parse_args <- function(args) {
  out <- list(
    base_dir = NULL,
    out_dir = NULL,
    ng = 3L,
    nstart = 5L,
    ncores = 1L,
    maxit = 150L,
    max_months = 36L,
    bin_days = 14L,
    min_visits = 4L,
    sample_n = 1000L,
    seed = 2080L
  )
  if (length(args) == 0L) return(out)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--base-dir" && i < length(args)) {
      out$base_dir <- args[[i + 1L]]
      i <- i + 2L
      next
    }
    if (a == "--out-dir" && i < length(args)) {
      out$out_dir <- args[[i + 1L]]
      i <- i + 2L
      next
    }
    if (a == "--ng" && i < length(args)) {
      out$ng <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (a == "--nstart" && i < length(args)) {
      out$nstart <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (a == "--ncores" && i < length(args)) {
      out$ncores <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (a == "--maxit" && i < length(args)) {
      out$maxit <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (a == "--max-months" && i < length(args)) {
      out$max_months <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (a == "--bin-days" && i < length(args)) {
      out$bin_days <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (a == "--min-visits" && i < length(args)) {
      out$min_visits <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (a == "--sample-n" && i < length(args)) {
      out$sample_n <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (a == "--seed" && i < length(args)) {
      out$seed <- as.integer(args[[i + 1L]])
      i <- i + 2L
      next
    }
    i <- i + 1L
  }
  out
}

detect_base_dir <- function(cli_base_dir) {
  if (!is.null(cli_base_dir) && dir.exists(cli_base_dir)) return(cli_base_dir)
  env_base <- Sys.getenv("HYPERTENSION_BASE_DIR", unset = "")
  if (nzchar(env_base) && dir.exists(env_base)) return(env_base)
  for (p in c("E:/hypertension", "G:/hypertension")) {
    if (dir.exists(p)) return(p)
  }
  stop("Could not detect project base directory. Pass --base-dir.")
}

stamp <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
base_dir <- detect_base_dir(args$base_dir)
lib_dir <- file.path(base_dir, "analysis", "Rlib")
if (dir.exists(lib_dir)) {
  # Prefer working environment packages for core dependencies while still
  # allowing gbmt to load from the project library.
  .libPaths(c(.libPaths(), lib_dir))
}

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
  library(jsonlite)
  library(parallel)
})

out_dir <- args$out_dir
if (is.null(out_dir)) {
  out_dir <- file.path(
    base_dir,
    "analysis",
    "outputs",
    paste0("q1_gbtm_relative_visit_only_shorttest_", format(Sys.Date(), "%Y-%m-%d"), "_ng", args$ng)
  )
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

visit_file <- file.path(base_dir, "visit_level_cohort.csv")

stamp("Reading visit-level cohort: ", visit_file)
dt <- fread(
  visit_file,
  select = c("new_patient_id", "encounter_dateD", "visit_sequence")
)
dt[, new_patient_id := as.integer(new_patient_id)]
dt[, encounter_dateD := as.numeric(encounter_dateD)]
dt[, visit_sequence := as.integer(visit_sequence)]
dt <- dt[!is.na(new_patient_id) & !is.na(encounter_dateD)]
setorder(dt, new_patient_id, encounter_dateD, visit_sequence)

max_days <- args$max_months * 30.4375
dt[, first_day := min(encounter_dateD), by = new_patient_id]
dt[, day_since_index := encounter_dateD - first_day]
dt <- dt[day_since_index >= 0 & day_since_index <= max_days]

visit_counts <- dt[, .(n_visits = .N), by = new_patient_id]
eligible_ids <- visit_counts[n_visits >= args$min_visits, new_patient_id]

set.seed(args$seed)
if (!is.na(args$sample_n) && args$sample_n > 0L && length(eligible_ids) > args$sample_n) {
  use_ids <- sort(sample(eligible_ids, size = args$sample_n, replace = FALSE))
} else {
  use_ids <- sort(eligible_ids)
}

work <- dt[new_patient_id %in% use_ids]
work[, time_bin := as.integer(floor(day_since_index / args$bin_days))]

bin_agg <- work[, .(
  visit_count_bin = as.integer(.N)
), by = .(new_patient_id, time_bin)]

last_bin <- bin_agg[, .(last_time_bin = max(time_bin)), by = new_patient_id]
last_bin <- last_bin[last_time_bin >= 1L]
panel <- last_bin[, .(time_bin = 0:last_time_bin), by = new_patient_id]
panel <- merge(panel, bin_agg, by = c("new_patient_id", "time_bin"), all.x = TRUE)
panel[is.na(visit_count_bin), visit_count_bin := 0L]
panel[, time_year := (time_bin * args$bin_days) / 365.25]

gbtm_long <- panel[, .(new_patient_id, time_year, visit_count_bin)]
fwrite(gbtm_long, file.path(out_dir, "q1_relative_biweekly_visit_only_panel.csv"))

stamp(
  "Panel rows: ", nrow(gbtm_long),
  "; eligible patients: ", length(eligible_ids),
  "; patients with >=2 bins: ", uniqueN(last_bin$new_patient_id),
  "; modeled patients: ", uniqueN(gbtm_long$new_patient_id),
  "; bin_days=", args$bin_days,
  "; ng=", args$ng,
  "; nstart=", args$nstart,
  "; ncores=", args$ncores,
  "; maxit=", args$maxit
)

t0 <- proc.time()[["elapsed"]]
model_df <- as.data.frame(gbtm_long)

run_one_start <- function(start_idx, data_df, ng, maxit, lib_dir) {
  if (!is.null(lib_dir) && dir.exists(lib_dir)) {
    .libPaths(c(.libPaths(), lib_dir))
  }
  suppressPackageStartupMessages(library(gbmt))
  set.seed(2080L + as.integer(start_idx))
  gbmt(
    x.names = "visit_count_bin",
    unit = "new_patient_id",
    time = "time_year",
    ng = as.integer(ng),
    d = 1,
    data = data_df,
    scaling = 0,
    pruning = FALSE,
    delete.empty = TRUE,
    nstart = 1L,
    maxit = as.integer(maxit),
    quiet = TRUE
  )
}

ncores <- max(1L, min(as.integer(args$ncores), as.integer(args$nstart), parallel::detectCores()))

if (ncores > 1L && .Platform$OS.type != "windows") {
  fits <- mclapply(
    seq_len(as.integer(args$nstart)),
    function(i) run_one_start(i, model_df, args$ng, args$maxit, lib_dir),
    mc.cores = ncores
  )
} else {
  fits <- lapply(seq_len(as.integer(args$nstart)), function(i) {
    run_one_start(i, model_df, args$ng, args$maxit, lib_dir)
  })
}

bics <- sapply(fits, function(f) as.numeric(f$ic["bic"]))
best_idx <- which.min(bics)
fit <- fits[[best_idx]]
elapsed <- proc.time()[["elapsed"]] - t0
fit_rds_path <- file.path(out_dir, "gbtm_relative_visit_only_fit.rds")
saveRDS(fit, fit_rds_path)

assign <- data.table(
  new_patient_id = as.integer(names(fit$assign)),
  gbtm_cluster = as.integer(fit$assign)
)
assign <- assign[!is.na(gbtm_cluster)]
sizes <- assign[, .N, by = gbtm_cluster][order(gbtm_cluster)]

fwrite(assign, file.path(out_dir, "gbtm_relative_visit_only_assignments.csv"))

postprob_path <- NA_character_
diagnostics_path <- NA_character_
if (!is.null(fit$postprob)) {
  post <- as.data.table(fit$postprob)
  prob_cols <- names(post)
  if (!("new_patient_id" %in% prob_cols)) {
    id_source <- rownames(fit$postprob)
    if (is.null(id_source) || anyNA(suppressWarnings(as.integer(id_source)))) {
      id_source <- names(fit$assign)
    }
    post[, new_patient_id := as.integer(id_source)]
  }
  prob_cols <- setdiff(names(post), "new_patient_id")
  setnames(post, prob_cols, paste0("posterior_prob_g", seq_along(prob_cols)))
  post <- merge(assign, post, by = "new_patient_id", all.x = TRUE)
  postprob_path <- file.path(out_dir, "gbtm_relative_visit_only_posterior_probabilities.csv")
  fwrite(post, postprob_path)

  n_assigned <- nrow(assign)
  diagnostics <- rbindlist(lapply(sort(unique(assign$gbtm_cluster)), function(g) {
    prob_col <- paste0("posterior_prob_g", g)
    assigned_prob <- post[gbtm_cluster == g, get(prob_col)]
    app <- mean(assigned_prob, na.rm = TRUE)
    prior <- nrow(post[gbtm_cluster == g]) / n_assigned
    posterior <- mean(post[[prob_col]], na.rm = TRUE)
    occ <- (app / (1 - app)) / (prior / (1 - prior))
    data.table(
      gbtm_cluster = as.integer(g),
      n_assigned = as.integer(sum(assign$gbtm_cluster == g)),
      prior_probability = prior,
      posterior_probability = posterior,
      prior_posterior_ratio = posterior / prior,
      average_posterior_probability = app,
      odds_correct_classification = occ
    )
  }))
  diagnostics_path <- file.path(out_dir, "gbtm_relative_visit_only_classification_diagnostics.csv")
  fwrite(diagnostics, diagnostics_path)
}

fit_row <- data.table(
  ng = as.integer(args$ng),
  degree = 1L,
  scaling = 0L,
  pruning = FALSE,
  nstart = as.integer(args$nstart),
  ncores = as.integer(ncores),
  maxit = as.integer(args$maxit),
  max_months = as.integer(args$max_months),
  bin_days = as.integer(args$bin_days),
  min_visits = as.integer(args$min_visits),
  sample_n = as.integer(args$sample_n),
  aic = as.numeric(fit$ic["aic"]),
  bic = as.numeric(fit$ic["bic"]),
  ssbic = as.numeric(fit$ic["ssbic"]),
  hqic = as.numeric(fit$ic["hqic"]),
  entropy = as.numeric(fit$ic["entropy"]),
  n_eligible_patients = as.integer(length(eligible_ids)),
  n_patients_ge_2_bins = as.integer(uniqueN(last_bin$new_patient_id)),
  n_modeled_patients = as.integer(uniqueN(assign$new_patient_id)),
  n_populated_clusters = as.integer(nrow(sizes)),
  min_cluster_n = as.integer(min(sizes$N)),
  max_cluster_n = as.integer(max(sizes$N)),
  elapsed_seconds = round(elapsed, 2)
)
fwrite(fit_row, file.path(out_dir, "gbtm_relative_visit_only_model_fit.csv"))

patient_pattern <- panel[, .(
  bins_observed = as.integer(.N),
  bins_with_any_visit = as.integer(sum(visit_count_bin > 0L)),
  avg_visits_per_bin = as.numeric(mean(visit_count_bin, na.rm = TRUE)),
  p90_visits_per_bin = as.numeric(quantile(visit_count_bin, probs = 0.9, na.rm = TRUE)),
  followup_days_modeled = as.numeric(max(time_bin, na.rm = TRUE) * args$bin_days)
), by = new_patient_id]

cluster_profile <- merge(assign, patient_pattern, by = "new_patient_id", all.x = TRUE)[, .(
  n = .N,
  pct = round(100 * .N / nrow(assign), 2),
  bins_observed = round(mean(bins_observed, na.rm = TRUE), 2),
  bins_with_any_visit = round(mean(bins_with_any_visit, na.rm = TRUE), 2),
  avg_visits_per_bin = round(mean(avg_visits_per_bin, na.rm = TRUE), 4),
  p90_visits_per_bin = round(mean(p90_visits_per_bin, na.rm = TRUE), 4),
  followup_days_modeled = round(mean(followup_days_modeled, na.rm = TRUE), 2)
), by = gbtm_cluster][order(gbtm_cluster)]
fwrite(cluster_profile, file.path(out_dir, "gbtm_relative_visit_only_cluster_profiles.csv"))

summary_list <- list(
  method = "q1_gbtm_relative_visit_only_shorttest",
  input_file = visit_file,
  n_patients_total = as.integer(uniqueN(dt$new_patient_id)),
  n_patients_eligible = as.integer(length(eligible_ids)),
  n_patients_ge_2_bins = as.integer(uniqueN(last_bin$new_patient_id)),
  n_patients_modeled = as.integer(uniqueN(assign$new_patient_id)),
  modeling_design = list(
    alignment = "relative_to_each_patients_first_visit",
    outcome = "visit_count_bin",
    time_variable = "time_year",
    bin_days = as.integer(args$bin_days),
    censoring_rule = "patient_specific_panel_ends_at_last_observed_bin"
  ),
  gbtm_ng = as.integer(args$ng),
  gbtm_degree = 1L,
  gbtm_scaling = 0L,
  gbtm_pruning = FALSE,
  gbtm_nstart = as.integer(args$nstart),
  gbtm_ncores = as.integer(ncores),
  gbtm_maxit = as.integer(args$maxit),
  min_visits = as.integer(args$min_visits),
  sampling = list(
    requested_sample_n = as.integer(args$sample_n),
    modeled_n = as.integer(uniqueN(assign$new_patient_id)),
    seed = as.integer(args$seed)
  ),
  outputs = list(
    panel = file.path(out_dir, "q1_relative_biweekly_visit_only_panel.csv"),
    model_fit = file.path(out_dir, "gbtm_relative_visit_only_model_fit.csv"),
    model_rds = fit_rds_path,
    assignments = file.path(out_dir, "gbtm_relative_visit_only_assignments.csv"),
    posterior_probabilities = postprob_path,
    classification_diagnostics = diagnostics_path,
    cluster_profiles = file.path(out_dir, "gbtm_relative_visit_only_cluster_profiles.csv")
  )
)

write_json(summary_list, file.path(out_dir, "q1_gbtm_relative_visit_only_summary.json"), pretty = TRUE, auto_unbox = TRUE)
stamp("Finished relative-time visit-only short test.")
cat(toJSON(summary_list, pretty = TRUE, auto_unbox = TRUE), "\n")
