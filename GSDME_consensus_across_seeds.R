#!/usr/bin/env Rscript

# =============================================================================
# GSDME consensus analysis across random seeds
# Version 1.4.0
#
# This script combines complete (non-quick) executions produced by
# GSDME_systems_oncology_pipeline.R. It quantifies Monte Carlo variability;
# it does not estimate biological variability and does not provide causal
# experimental validation.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

parse_arguments <- function(args) {
  cfg <- list(
    parent = getwd(),
    pattern = "^resultados_finais_GSDME_seed_[0-9]+$",
    folders = NA_character_,
    out = "resultados_consenso_GSDME",
    auto_install = TRUE
  )
  if (!length(args)) return(cfg)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--help") {
      cat(paste0(
        "Usage:\n",
        "  Rscript GSDME_consensus_across_seeds.R [options]\n\n",
        "Options:\n",
        "  --parent PATH      Parent directory containing seed folders\n",
        "  --pattern REGEX    Folder-name pattern used for discovery\n",
        "  --folders LIST     Comma-separated folder names or paths\n",
        "  --out PATH         Consensus output directory\n",
        "  --no-install       Do not install missing CRAN packages\n"
      ))
      quit(save = "no", status = 0)
    }
    if (key == "--no-install") {
      cfg$auto_install <- FALSE
      i <- i + 1L
      next
    }
    if (i == length(args)) stop("Missing value after ", key)
    value <- args[[i + 1L]]
    if (key == "--parent") cfg$parent <- value
    else if (key == "--pattern") cfg$pattern <- value
    else if (key == "--folders") cfg$folders <- value
    else if (key == "--out") cfg$out <- value
    else stop("Unknown argument: ", key)
    i <- i + 2L
  }
  cfg
}

CFG <- parse_arguments(commandArgs(trailingOnly = TRUE))
CFG$parent <- normalizePath(CFG$parent, mustWork = TRUE)
if (!grepl("^(/|[A-Za-z]:[/\\\\])", CFG$out)) {
  CFG$out <- file.path(CFG$parent, CFG$out)
}
CFG$out <- normalizePath(CFG$out, mustWork = FALSE)

if (!is.na(CFG$folders)) {
  run_dirs <- trimws(strsplit(CFG$folders, ",", fixed = TRUE)[[1]])
  run_dirs <- vapply(run_dirs, function(x) {
    if (grepl("^(/|[A-Za-z]:[/\\\\])", x)) x else file.path(CFG$parent, x)
  }, character(1))
} else {
  candidates <- list.dirs(CFG$parent, recursive = FALSE, full.names = TRUE)
  run_dirs <- candidates[grepl(CFG$pattern, basename(candidates), perl = TRUE)]
}

run_dirs <- sort(unique(normalizePath(run_dirs, mustWork = FALSE)))
if (!length(run_dirs)) {
  stop(
    "No seed folders found. Expected names such as ",
    "resultados_finais_GSDME_seed_101 inside: ", CFG$parent
  )
}
if (any(!dir.exists(run_dirs))) {
  stop("Folder(s) not found: ", paste(run_dirs[!dir.exists(run_dirs)], collapse = ", "))
}

extract_seed <- function(path) {
  hit <- regexec("(?:^|_)seed_([0-9]+)(?:$|_)", basename(path), perl = TRUE)
  value <- regmatches(basename(path), hit)[[1]]
  if (length(value) < 2L) return(NA_integer_)
  as.integer(value[[2]])
}

folder_seeds <- vapply(run_dirs, extract_seed, integer(1))
if (anyNA(folder_seeds)) {
  stop("Every folder name must contain '_seed_NUMBER': ",
       paste(basename(run_dirs)[is.na(folder_seeds)], collapse = ", "))
}
if (anyDuplicated(folder_seeds)) stop("Duplicated seed numbers were found.")

dir.create(CFG$out, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(CFG$out, "tables")
figure_dir <- file.path(CFG$out, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(CFG$out, "consensus.log")
log_message <- function(...) {
  msg <- paste0(...)
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

write_csv <- function(x, filename) {
  utils::write.csv(x, file.path(table_dir, filename), row.names = FALSE, na = "")
}

required_files <- c(
  "01_model_nodes_and_rules.csv",
  "02_model_edges.csv",
  "05_perturbation_activation_frequencies.csv",
  "06_structural_control_ranking.csv",
  "07_minimum_driver_set_search.csv",
  "11_reinforcement_learning_sequence.csv",
  "12_output_manifest.csv",
  "15_geo_node_coverage.csv",
  "16_geo_perturbation_concordance.csv",
  "17_geo_group_fate_probabilities.csv",
  "19_geo_cross_dataset_summary.csv",
  "20_RL_pyroptosis_target_summary.csv",
  "22_exact_stable_state_phenotypes_DDR_ON.csv",
  "23_phenotype_reachability_DDR_ON.csv",
  "26_GEO_fixed_perturbation_independent_phenotypes.csv",
  "Figure_01_GSDME_logical_network_600dpi.png",
  "Figure_01_GSDME_logical_network.pdf",
  "Figure_06_GEO_projected_fate_composition_600dpi.png",
  "Figure_06b_seven_priority_nodes_RNA_detection_600dpi.png",
  "pipeline.log",
  "sessionInfo.txt"
)

locate_file <- function(run_dir, filename) {
  direct <- file.path(run_dir, filename)
  if (file.exists(direct)) return(direct)
  hits <- list.files(run_dir, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == filename]
  if (!length(hits)) return(NA_character_)
  hits[[1]]
}

read_run_csv <- function(run_dir, filename, seed) {
  path <- locate_file(run_dir, filename)
  if (is.na(path)) stop("Missing ", filename, " in ", run_dir)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  x$seed <- seed
  x$run_folder <- basename(run_dir)
  x
}

safe_unique <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (!length(x)) NA_character_ else paste(x, collapse = ";")
}

# --------------------------- Quality control --------------------------------

qc_rows <- lapply(seq_along(run_dirs), function(i) {
  run_dir <- run_dirs[[i]]
  seed <- folder_seeds[[i]]
  found <- vapply(required_files, function(f) !is.na(locate_file(run_dir, f)), logical(1))

  manifest_path <- locate_file(run_dir, "12_output_manifest.csv")
  manifest_seed <- NA_integer_
  pipeline_version <- NA_character_
  possible_quick <- TRUE
  if (!is.na(manifest_path)) {
    manifest <- utils::read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
    if ("seed" %in% names(manifest)) {
      manifest_seed <- suppressWarnings(as.integer(unique(manifest$seed)[1]))
    }
    if ("pipeline_version" %in% names(manifest)) {
      pipeline_version <- safe_unique(as.character(manifest$pipeline_version))
    }
    if ("quick" %in% names(manifest) && nrow(manifest)) {
      possible_quick <- any(tolower(as.character(manifest$quick)) != "false")
    }
  }

  log_path <- locate_file(run_dir, "pipeline.log")
  completed <- FALSE
  if (!is.na(log_path)) {
    completed <- any(grepl("Pipeline completed", readLines(log_path, warn = FALSE), fixed = TRUE))
  }

  coverage_path <- locate_file(run_dir, "15_geo_node_coverage.csv")
  max_geo_cells <- NA_real_
  min_geo_cells <- NA_real_
  if (!is.na(coverage_path)) {
    coverage <- utils::read.csv(coverage_path, check.names = FALSE, stringsAsFactors = FALSE)
    if ("n_cells" %in% names(coverage) && any(is.finite(coverage$n_cells))) {
      max_geo_cells <- max(coverage$n_cells, na.rm = TRUE)
      min_geo_cells <- min(coverage$n_cells, na.rm = TRUE)
    }
  }

  nodes_path <- locate_file(run_dir, "01_model_nodes_and_rules.csv")
  edges_path <- locate_file(run_dir, "02_model_edges.csv")
  nodes_md5 <- if (!is.na(nodes_path)) unname(tools::md5sum(nodes_path)) else NA_character_
  edges_md5 <- if (!is.na(edges_path)) unname(tools::md5sum(edges_path)) else NA_character_

  # Real HCC subsets can contain fewer than 500 cells even in a full run.
  # The seed manifest records whether quick mode was actually used.
  all_required <- all(found)
  seed_matches <- !is.na(manifest_seed) && identical(seed, manifest_seed)

  data.frame(
    seed = seed,
    run_folder = basename(run_dir),
    pipeline_completed = completed,
    all_required_files = all_required,
    missing_files = if (all_required) "" else paste(required_files[!found], collapse = ";"),
    manifest_seed = manifest_seed,
    seed_matches_folder = seed_matches,
    pipeline_version = pipeline_version,
    min_geo_cells = min_geo_cells,
    max_geo_cells = max_geo_cells,
    possible_quick_run = possible_quick,
    model_nodes_md5 = nodes_md5,
    model_edges_md5 = edges_md5,
    stringsAsFactors = FALSE
  )
})

qc <- do.call(rbind, qc_rows)
same_version <- length(unique(qc$pipeline_version[!is.na(qc$pipeline_version)])) == 1L
same_nodes <- length(unique(qc$model_nodes_md5[!is.na(qc$model_nodes_md5)])) == 1L
same_edges <- length(unique(qc$model_edges_md5[!is.na(qc$model_edges_md5)])) == 1L
qc$same_pipeline_version <- same_version
qc$same_model_nodes <- same_nodes
qc$same_model_edges <- same_edges
qc$quality_control_pass <- with(
  qc,
  pipeline_completed & all_required_files & seed_matches_folder &
    !possible_quick_run & same_pipeline_version & same_model_nodes & same_model_edges
)
write_csv(qc, "01_run_quality_control.csv")

if (any(!qc$quality_control_pass)) {
  stop(
    "Quality control failed. Open tables/01_run_quality_control.csv. ",
    "Do not combine incomplete, possible --quick, mismatched-version, or mismatched-model runs."
  )
}
if (nrow(qc) < 5L) {
  warning("Only ", nrow(qc), " seeds were found. Five or more are recommended for the final analysis.")
}
log_message("Quality control passed for seeds: ", paste(folder_seeds, collapse = ", "))

# Figure 1 is structural and therefore seed-independent. Because the model
# tables were verified as byte-identical above, copy it from the first run.
network_png <- locate_file(run_dirs[[1]], "Figure_01_GSDME_logical_network_600dpi.png")
network_pdf <- locate_file(run_dirs[[1]], "Figure_01_GSDME_logical_network.pdf")
file.copy(network_png, file.path(figure_dir, basename(network_png)), overwrite = TRUE)
file.copy(network_pdf, file.path(figure_dir, basename(network_pdf)), overwrite = TRUE)

# -------------------------- Generic summaries -------------------------------

bind_seed_tables <- function(filename) {
  do.call(rbind, lapply(seq_along(run_dirs), function(i) {
    read_run_csv(run_dirs[[i]], filename, folder_seeds[[i]])
  }))
}

mean_ci <- function(x, conf = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(c(n = 0, mean = NA, sd = NA, se = NA, ci_low = NA,
                   ci_high = NA, min = NA, max = NA))
  avg <- mean(x)
  sdev <- if (n > 1L) stats::sd(x) else NA_real_
  se <- if (n > 1L) sdev / sqrt(n) else NA_real_
  margin <- if (n > 1L) stats::qt(1 - (1 - conf) / 2, df = n - 1L) * se else NA_real_
  c(n = n, mean = avg, sd = sdev, se = se,
    ci_low = if (n > 1L) avg - margin else NA_real_,
    ci_high = if (n > 1L) avg + margin else NA_real_,
    min = min(x), max = max(x))
}

group_summary <- function(data, keys, metrics) {
  missing_columns <- setdiff(c(keys, metrics, "seed"), names(data))
  if (length(missing_columns)) {
    stop("Missing column(s): ", paste(missing_columns, collapse = ", "))
  }
  key_text <- do.call(paste, c(lapply(data[keys], as.character), sep = "\036"))
  groups <- split(seq_len(nrow(data)), key_text, drop = TRUE)
  rows <- lapply(groups, function(idx) {
    first <- data[idx[1], keys, drop = FALSE]
    first$n_seeds_present <- length(unique(data$seed[idx]))
    for (metric in metrics) {
      stats <- mean_ci(data[[metric]][idx])
      for (stat_name in names(stats)) {
        first[[paste(metric, stat_name, sep = "_")]] <- unname(stats[[stat_name]])
      }
    }
    first
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

clamp01 <- function(x) pmin(1, pmax(0, x))

# ----------------------- Perturbation consensus ------------------------------

perturb_all <- bind_seed_tables("05_perturbation_activation_frequencies.csv")
perturb_consensus <- group_summary(
  perturb_all,
  c("perturbation", "node"),
  c("activation_frequency", "stable_fraction")
)
perturb_consensus$activation_frequency_ci_low <- clamp01(perturb_consensus$activation_frequency_ci_low)
perturb_consensus$activation_frequency_ci_high <- clamp01(perturb_consensus$activation_frequency_ci_high)
perturb_consensus$stable_fraction_ci_low <- clamp01(perturb_consensus$stable_fraction_ci_low)
perturb_consensus$stable_fraction_ci_high <- clamp01(perturb_consensus$stable_fraction_ci_high)
write_csv(perturb_consensus, "02_perturbation_consensus.csv")

# Direct GSDME dependency contrasts. These within-seed comparisons use DDR ON
# and independent simulated trajectories per condition; matched initial-state
# transitions are evaluated separately in the epistasis script.
dependency_pairs <- data.frame(
  comparison = c("CASP3 E1", "miR-204-5p E1",
                 "SIRT1 KO", "MALAT1 KO", "BAX E1"),
  gsdme_endogenous = c("CASP3 OE", "miR-204-5p OE", "SIRT1 KO",
                       "MALAT1 KO", "BAX E1"),
  gsdme_ko = c("CASP3 OE + GSDME KO", "miR-204-5p OE + GSDME KO",
               "SIRT1 KO + GSDME KO", "MALAT1 KO + GSDME KO",
               "BAX E1 + GSDME KO"),
  stringsAsFactors = FALSE
)
dependency_phenotypes <- c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE",
                           "PROLIFERATION", "CELL_CYCLE_ARREST")
expected_conditions <- unique(c(dependency_pairs$gsdme_endogenous,
                                dependency_pairs$gsdme_ko))
if (!all(expected_conditions %in% perturb_all$perturbation) ||
    !all(dependency_phenotypes %in% perturb_all$node)) {
  stop("Missing GSDME-dependency conditions or phenotypes. Re-run every seed ",
       "using the updated seed script before computing consensus.")
}
dependency_rows <- lapply(seq_len(nrow(dependency_pairs)), function(index) {
  pair <- dependency_pairs[index, ]
  endogenous <- perturb_all[
    perturb_all$perturbation == pair$gsdme_endogenous &
      perturb_all$node %in% dependency_phenotypes,
    c("seed", "node", "activation_frequency")
  ]
  knockout <- perturb_all[
    perturb_all$perturbation == pair$gsdme_ko &
      perturb_all$node %in% dependency_phenotypes,
    c("seed", "node", "activation_frequency")
  ]
  names(endogenous)[[3]] <- "endogenous_frequency"
  names(knockout)[[3]] <- "gsdme_ko_frequency"
  result <- merge(endogenous, knockout, by = c("seed", "node"))
  if (nrow(result) != length(folder_seeds) * length(dependency_phenotypes)) {
    stop("Incomplete GSDME dependency pair: ", pair$comparison)
  }
  result$comparison <- pair$comparison
  result$endogenous_condition <- pair$gsdme_endogenous
  result$knockout_condition <- pair$gsdme_ko
  result$delta_endogenous_minus_ko <- result$endogenous_frequency -
    result$gsdme_ko_frequency
  result
})
dependency_seed_contrasts <- do.call(rbind, dependency_rows)
write_csv(dependency_seed_contrasts, "15_GSDME_dependency_contrasts_by_seed.csv")
dependency_contrasts <- group_summary(
  dependency_seed_contrasts, c("comparison", "node"),
  c("endogenous_frequency", "gsdme_ko_frequency", "delta_endogenous_minus_ko")
)
write_csv(dependency_contrasts, "16_GSDME_dependency_contrasts_consensus.csv")

# --------------------- Structural-control consensus --------------------------

structural_all <- bind_seed_tables("06_structural_control_ranking.csv")
structural_metrics <- intersect(
  c("indegree", "outdegree", "betweenness", "pagerank", "in_feedback_scc", "score"),
  names(structural_all)
)
structural_consensus <- group_summary(structural_all, "node", structural_metrics)
write_csv(structural_consensus, "03_structural_control_consensus.csv")

# -------------------------- Driver consensus ---------------------------------

driver_all <- bind_seed_tables("07_minimum_driver_set_search.csv")
driver_all$robust_numeric <- as.numeric(as.logical(driver_all$robust))
driver_consensus <- group_summary(
  driver_all,
  c("size", "intervention"),
  c("success_rate", "robust_numeric")
)
driver_consensus$selection_frequency <- driver_consensus$n_seeds_present / length(folder_seeds)
driver_consensus$success_rate_ci_low <- clamp01(driver_consensus$success_rate_ci_low)
driver_consensus$success_rate_ci_high <- clamp01(driver_consensus$success_rate_ci_high)
names(driver_consensus)[names(driver_consensus) == "robust_numeric_mean"] <- "robust_seed_fraction"
write_csv(driver_consensus, "04_driver_set_consensus.csv")

parse_driver_nodes <- function(intervention) {
  pieces <- unlist(strsplit(intervention, "\\s*(?:\\+|;|,)\\s*", perl = TRUE))
  pieces <- trimws(pieces)
  pieces[nzchar(pieces)]
}

driver_presence <- do.call(rbind, lapply(split(driver_all, driver_all$seed), function(d) {
  tokens <- unique(unlist(lapply(d$intervention[d$success_rate >= 0.95], parse_driver_nodes)))
  if (!length(tokens)) return(NULL)
  data.frame(seed = unique(d$seed)[1], driver_action = tokens, stringsAsFactors = FALSE)
}))
if (is.null(driver_presence) || !nrow(driver_presence)) {
  driver_node_frequency <- data.frame(
    driver_action = character(), seeds_selected = integer(),
    selection_frequency = numeric(), stringsAsFactors = FALSE
  )
} else {
  driver_node_frequency <- stats::aggregate(
    seed ~ driver_action, driver_presence,
    function(x) length(unique(x))
  )
  names(driver_node_frequency)[2] <- "seeds_selected"
  driver_node_frequency$selection_frequency <- driver_node_frequency$seeds_selected / length(folder_seeds)
  driver_node_frequency <- driver_node_frequency[order(-driver_node_frequency$selection_frequency,
                                                       driver_node_frequency$driver_action), ]
}
write_csv(driver_node_frequency, "05_driver_action_frequency.csv")

# --------------------------- GEO consensus -----------------------------------

coverage_all <- bind_seed_tables("15_geo_node_coverage.csv")
coverage_consensus <- group_summary(
  coverage_all,
  c("accession", "experiment", "platform", "validation_scope", "node", "measured", "matched_aliases"),
  c("detection_rate", "mean_activation_probability", "n_cells")
)
coverage_consensus$detection_rate_ci_low <- clamp01(coverage_consensus$detection_rate_ci_low)
coverage_consensus$detection_rate_ci_high <- clamp01(coverage_consensus$detection_rate_ci_high)
write_csv(coverage_consensus, "06_GEO_node_coverage_consensus.csv")

concordance_all <- bind_seed_tables("16_geo_perturbation_concordance.csv")
concordance_consensus <- group_summary(
  concordance_all,
  c("accession", "experiment", "validation_group", "validation_scope", "perturbation"),
  c("concordance_score", "weighted_concordance_score", "spearman_rank_correlation",
    "nodes_compared", "node_coverage", "n_cells")
)
concordance_consensus$concordance_score_ci_low <- clamp01(concordance_consensus$concordance_score_ci_low)
concordance_consensus$concordance_score_ci_high <- clamp01(concordance_consensus$concordance_score_ci_high)
concordance_consensus$weighted_concordance_score_ci_low <- clamp01(concordance_consensus$weighted_concordance_score_ci_low)
concordance_consensus$weighted_concordance_score_ci_high <- clamp01(concordance_consensus$weighted_concordance_score_ci_high)
write_csv(concordance_consensus, "07_GEO_perturbation_concordance_consensus.csv")

fate_all <- bind_seed_tables("17_geo_group_fate_probabilities.csv")
fate_consensus <- group_summary(
  fate_all,
  c("accession", "experiment", "validation_group", "validation_scope", "fate"),
  c("mean_probability", "n_cells")
)
fate_consensus$mean_probability_ci_low <- clamp01(fate_consensus$mean_probability_ci_low)
fate_consensus$mean_probability_ci_high <- clamp01(fate_consensus$mean_probability_ci_high)
write_csv(fate_consensus, "08_GEO_fate_consensus.csv")

fixed_phenotypes_all <- bind_seed_tables(
  "26_GEO_fixed_perturbation_independent_phenotypes.csv"
)
fixed_phenotypes_consensus <- group_summary(
  fixed_phenotypes_all,
  c("accession", "experiment", "validation_group", "validation_scope",
    "condition", "phenotype"),
  c("mean_probability", "n_stable_trajectories", "n_trajectories", "n_cells")
)
fixed_phenotypes_consensus$mean_probability_ci_low <- clamp01(
  fixed_phenotypes_consensus$mean_probability_ci_low
)
fixed_phenotypes_consensus$mean_probability_ci_high <- clamp01(
  fixed_phenotypes_consensus$mean_probability_ci_high
)
write_csv(fixed_phenotypes_consensus,
          "14_GEO_fixed_perturbation_independent_phenotypes_consensus.csv")

cross_all <- bind_seed_tables("19_geo_cross_dataset_summary.csv")
cross_consensus <- group_summary(
  cross_all,
  c("accession", "perturbation", "validation_role", "causal_confirmation"),
  c("mean_concordance", "mean_weighted_concordance", "mean_node_coverage")
)
write_csv(cross_consensus, "09_GEO_cross_dataset_consensus.csv")

# ----------------------- Reinforcement-learning consensus --------------------

rl_summary_all <- bind_seed_tables("20_RL_pyroptosis_target_summary.csv")
rl_summary_all$reached_numeric <- as.numeric(as.logical(rl_summary_all$reached))
rl_consensus <- group_summary(
  rl_summary_all,
  "target",
  c("reached_numeric", "steps_used", "horizon")
)
names(rl_consensus)[names(rl_consensus) == "reached_numeric_mean"] <- "success_fraction"
rl_consensus$successful_seeds <- sum(rl_summary_all$reached_numeric == 1, na.rm = TRUE)
rl_consensus$total_seeds <- length(folder_seeds)
rl_consensus$final_fates <- safe_unique(rl_summary_all$final_fate)
write_csv(rl_consensus, "10_RL_target_consensus.csv")

rl_sequence_all <- bind_seed_tables("11_reinforcement_learning_sequence.csv")
max_step <- max(rl_sequence_all$step, na.rm = TRUE)
action_grid <- expand.grid(
  seed = folder_seeds,
  step = seq_len(max_step),
  stringsAsFactors = FALSE
)
action_grid <- merge(
  action_grid,
  rl_sequence_all[c("seed", "step", "action")],
  by = c("seed", "step"), all.x = TRUE, sort = FALSE
)
action_grid$action[is.na(action_grid$action)] <- "No recorded step"
action_counts <- stats::aggregate(
  seed ~ step + action, action_grid,
  function(x) length(unique(x))
)
names(action_counts)[3] <- "seeds"
action_counts$frequency_all_seeds <- action_counts$seeds / length(folder_seeds)
action_counts <- action_counts[order(action_counts$step, -action_counts$frequency_all_seeds,
                                     action_counts$action), ]
write_csv(action_counts, "11_RL_action_stability.csv")

modal_action <- do.call(rbind, lapply(split(action_counts, action_counts$step), function(d) {
  d <- d[order(-d$frequency_all_seeds, d$action), ]
  d[1, c("step", "action", "seeds", "frequency_all_seeds"), drop = FALSE]
}))
rownames(modal_action) <- NULL
write_csv(modal_action, "12_RL_modal_action_by_step.csv")

# Seed-level endpoints for transparent supplementary reporting.
phenotype_nodes <- intersect(
  c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE", "PROLIFERATION",
    "CELL_CYCLE_ARREST"),
  unique(perturb_all$node)
)
seed_endpoints <- perturb_all[perturb_all$node %in% phenotype_nodes,
                              c("seed", "perturbation", "node", "activation_frequency",
                                "stable_fraction", "run_folder")]
write_csv(seed_endpoints, "13_seed_level_phenotype_endpoints.csv")

# Exact fixed points and seed-level empirical convergence are distinct metrics.
# Fixed points are enumerated from the Boolean rules; reachability reflects
# random initial states and asynchronous updates under DDR ON.
exact_seed_states <- bind_seed_tables("22_exact_stable_state_phenotypes_DDR_ON.csv")
reachability_seed <- bind_seed_tables("23_phenotype_reachability_DDR_ON.csv")
if (!all(exact_seed_states$DDR == "DDR ON") ||
    !all(reachability_seed$DDR == "DDR ON")) {
  stop("DDR OFF results cannot be included: regenerate every seed with DDR ON")
}
exact_keys <- c("condition", "phenotype")
exact_groups <- split(exact_seed_states,
                      interaction(exact_seed_states[exact_keys], drop = TRUE))
if (any(vapply(exact_groups, function(group) {
  length(unique(group$n_stable_states)) != 1L ||
    length(unique(group$states_with_phenotype)) != 1L ||
    length(unique(group$exact_enumeration_complete)) != 1L
}, logical(1)))) {
  stop("Exact stable states differ across seeds; inspect model and rules")
}
exact_summary <- exact_seed_states[
  exact_seed_states$seed == folder_seeds[[1]],
  c("DDR", "condition", "phenotype", "exact_enumeration_complete",
    "n_stable_states", "states_with_phenotype")
]
write_csv(exact_summary, "17_exact_stable_states_consensus_DDR_ON.csv")
convergence_summary <- group_summary(
  reachability_seed, c("DDR", "condition", "phenotype"),
  c("reachability", "n_stable_trajectories")
)
convergence_summary$reachability_ci_low <- clamp01(
  convergence_summary$reachability_ci_low
)
convergence_summary$reachability_ci_high <- clamp01(
  convergence_summary$reachability_ci_high
)
write_csv(convergence_summary, "18_convergence_consensus_DDR_ON.csv")

# ------------------------------- Figures -------------------------------------

ensure_package <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) return(invisible(TRUE))
  if (!CFG$auto_install) stop("Missing package: ", package)
  install.packages(package, repos = "https://cloud.r-project.org")
  if (!requireNamespace(package, quietly = TRUE)) stop("Could not install package: ", package)
}
ensure_package("ggplot2")
ensure_package("scales")

save_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    filename = file.path(figure_dir, paste0(stem, "_600dpi.png")),
    plot = plot, width = width, height = height, units = "in", dpi = 600,
    bg = "white", limitsize = FALSE
  )
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(
    filename = file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot, width = width, height = height, units = "in",
    device = pdf_device, bg = "white", limitsize = FALSE
  )
}

pretty_node <- function(x) {
  labels <- c(
    PYROPTOSIS = "GSDME-mediated pyroptosis",
    CELL_CYCLE_ARREST = "Cell-cycle arrest",
    PROLIFERATION = "Proliferation",
    RESISTANCE = "Resistance",
    APOPTOSIS = "Apoptosis",
    OTHER = "Other"
  )
  out <- unname(labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

pretty_condition <- function(x) {
  x <- gsub(" OE", " E1", x, fixed = TRUE)
  x[x == "CASP3 E1 + GSDME KO"] <- "GSDME KO + CASP3 E1"
  x
}

pretty_regulator <- function(x) {
  labels <- c(
    DFNA5 = "GSDME (DFNA5)",
    TP53 = "p53 (TP53)",
    MALAT1 = "MALAT1",
    MIR204 = "miR-204-5p (MIR204)",
    `(CCND1 ou CCND2 ou CCND3) E (CDK4 ou CDK6)` = "Cyclin D-CDK4/6"
  )
  out <- unname(labels[x])
  out[is.na(out)] <- gsub("_", " ", x[is.na(out)], fixed = TRUE)
  out
}

plot_theme <- ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(colour = "#444444"),
    plot.margin = ggplot2::margin(12, 18, 12, 18)
  )

# C1: Mean phenotype activation across seeds.
heat <- perturb_consensus[perturb_consensus$node %in% phenotype_nodes, ]
heat$node_label <- factor(pretty_node(heat$node), levels = pretty_node(phenotype_nodes))
heat$perturbation <- factor(pretty_condition(heat$perturbation),
                           levels = rev(pretty_condition(unique(heat$perturbation))))
p1 <- ggplot2::ggplot(heat, ggplot2::aes(node_label, perturbation,
                                         fill = activation_frequency_mean)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", activation_frequency_mean)),
                     size = 3.2) +
  ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#B2182B", limits = c(0, 1),
                               name = "Mean\nfrequency") +
  ggplot2::labs(
    title = "Consensus in silico perturbation screen",
    subtitle = paste0("Mean across ", length(folder_seeds),
                      " computational seeds; intervals are provided in Table 02"),
    x = NULL, y = NULL
  ) +
  plot_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 42, hjust = 1))
save_plot(p1, "Figure_02_in_silico_perturbation_consensus", 12.5, 8)

# Exact stable states and estimated stochastic convergence, strictly DDR ON.
convergence_focus <- c(
  "Reference: DDR ON, GSDME endogenous", "MALAT1 OE", "MALAT1 KO",
  "miR-204-5p OE", "miR-204-5p KO", "SIRT1 OE", "SIRT1 KO",
  "CASP3 OE", "CASP3 OE + GSDME KO", "GSDME E1 + CASP3 E1",
  "miR-204-5p OE + GSDME KO",
  "p53 E1 + miR-204-5p E1 + CASP3 KO + p21 E1",
  "MALAT1 E1 + miR-204-5p KO + SIRT1 E1",
  "SIRT1 KO + GSDME KO + CASP3 E1"
)
stable_convergence <- merge(
  exact_summary, convergence_summary,
  by = c("DDR", "condition", "phenotype"), sort = FALSE
)
stable_convergence <- stable_convergence[
  stable_convergence$condition %in% convergence_focus &
    stable_convergence$exact_enumeration_complete &
    !is.na(stable_convergence$n_stable_states) &
    stable_convergence$n_stable_states > 0L, , drop = FALSE
]
write_csv(stable_convergence, "20_stable_states_and_convergence_DDR_ON.csv")
if (nrow(stable_convergence)) {
  stable_convergence$condition_label <- factor(
    pretty_condition(stable_convergence$condition),
    levels = rev(pretty_condition(convergence_focus))
  )
  stable_convergence$phenotype_label <- factor(
    pretty_node(stable_convergence$phenotype),
    levels = pretty_node(phenotype_nodes)
  )
  stable_convergence$cell_label <- with(
    stable_convergence,
    sprintf("%d/%d fixed\n%.0f%% paths", states_with_phenotype,
            n_stable_states, 100 * reachability_mean)
  )
  p_stable_convergence <- ggplot2::ggplot(
    stable_convergence,
    ggplot2::aes(x = phenotype_label, y = condition_label,
                 fill = reachability_mean)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 2.55,
                       lineheight = 0.85) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#0072B2",
                                 limits = c(0, 1),
                                 labels = scales::percent_format()) +
    ggplot2::labs(
      title = "Stable states and estimated convergence in the miR-204-5p/MALAT1/SIRT1 network",
      subtitle = "DDR ON | exact fixed-point counts and mean across five stochastic seeds",
      caption = paste("Top: phenotype-ON exact fixed points / all fixed points.",
                      "Bottom: fraction of random starts converging to that output; not an exact state probability."),
      x = NULL, y = NULL, fill = "Paths") + plot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 32, hjust = 1))
  save_plot(p_stable_convergence,
            "Figure_03a_stable_states_and_convergence_DDR_ON", 18, 10)
}

loop_conditions_consensus <- c(
  "Reference: DDR ON, GSDME endogenous", "miR-204-5p OE",
  "miR-204-5p KO", "MALAT1 OE", "MALAT1 KO",
  "SIRT1 OE", "SIRT1 KO", "p53 E1", "p53 KO", "p21 E1",
  "p53 E1 + miR-204-5p E1 + CASP3 KO + p21 E1",
  "MALAT1 E1 + miR-204-5p KO + SIRT1 E1"
)
loop_nodes_consensus <- c("MALAT1", "MIR204", "SIRT1", "TP53",
                          "CDKN1A", "BAX", "CASP3", "DFNA5",
                          "PYROPTOSIS", "APOPTOSIS")
loop_consensus <- perturb_consensus[
  perturb_consensus$perturbation %in% loop_conditions_consensus &
    perturb_consensus$node %in% loop_nodes_consensus, , drop = FALSE
]
write_csv(loop_consensus, "19_feedback_loop_consensus_DDR_ON.csv")
loop_consensus$condition_label <- factor(
  pretty_condition(loop_consensus$perturbation),
  levels = rev(pretty_condition(loop_conditions_consensus))
)
loop_labels <- c(MALAT1 = "MALAT1", MIR204 = "miR-204-5p", SIRT1 = "SIRT1",
                 TP53 = "p53", CDKN1A = "p21", BAX = "BAX", CASP3 = "CASP3",
                 DFNA5 = "GSDME", PYROPTOSIS = "Pyroptosis",
                 APOPTOSIS = "Apoptosis")
loop_consensus$node_label <- factor(
  unname(loop_labels[loop_consensus$node]),
  levels = unname(loop_labels[loop_nodes_consensus])
)
p_loop_consensus <- ggplot2::ggplot(
  loop_consensus, ggplot2::aes(node_label, condition_label,
                             fill = activation_frequency_mean)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(
    label = sprintf("%.0f%%", 100 * activation_frequency_mean)), size = 2.8) +
  ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#0072B2",
                               limits = c(0, 1),
                               labels = scales::percent_format()) +
  ggplot2::labs(
    title = "Perturbation outcomes in miR-204-5p/SIRT1/p53 and p53/MALAT1 loops",
    subtitle = "DDR ON | independent node activities across five simulation seeds",
    caption = "E1 = Boolean ON; KO = Boolean OFF. Outputs are not mutually exclusive.",
    x = NULL, y = NULL, fill = "Node ON") + plot_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 32, hjust = 1))
save_plot(p_loop_consensus, "Figure_03b_feedback_loop_outcomes_DDR_ON", 17, 8)

# Detailed GSDME panel, with independent scores for all five phenotype nodes.
priority_conditions <- c(
  "Reference: DDR ON, GSDME endogenous", "GSDME OE", "GSDME KO",
  "CASP3 OE", "CASP3 OE + GSDME KO", "CASP3 KO",
  "miR-204-5p OE", "miR-204-5p OE + GSDME KO",
  "miR-204-5p OE + CASP3 KO", "SIRT1 KO", "SIRT1 KO + GSDME KO",
  "MALAT1 KO", "MALAT1 KO + GSDME KO",
  "GSDME E1 + CASP3 E1", "SIRT1 KO + GSDME KO + CASP3 E1",
  "p53 E1 + miR-204-5p E1 + CASP3 KO + p21 E1",
  "MALAT1 E1 + miR-204-5p KO + SIRT1 E1",
  "BAX E1", "BAX E1 + GSDME KO"
)
dep_heat <- perturb_consensus[
  perturb_consensus$perturbation %in% priority_conditions &
    perturb_consensus$node %in% dependency_phenotypes, , drop = FALSE
]
dep_heat$perturbation <- factor(pretty_condition(dep_heat$perturbation),
                                levels = rev(pretty_condition(priority_conditions)))
dep_heat$phenotype_label <- factor(pretty_node(dep_heat$node),
                                   levels = pretty_node(dependency_phenotypes))
p_gsdme_heat <- ggplot2::ggplot(
  dep_heat, ggplot2::aes(x = phenotype_label, y = perturbation,
                        fill = activation_frequency_mean)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.0f%%", 100 * activation_frequency_mean)),
    size = 3.2
  ) +
  ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#0072B2",
                               limits = c(0, 1), name = "Mean output ON") +
  ggplot2::labs(
    title = "Perturbation based validation of GSDME dependence in silico",
    subtitle = paste0("Five independent phenotype outputs across ",
                      length(folder_seeds), " computational seeds"),
    caption = paste("Outcomes can coexist. GSDME E1 is a Boolean DFNA5 clamp,",
                    "not a cleavage measurement."), x = NULL, y = NULL
  ) + plot_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 32, hjust = 1))
save_plot(p_gsdme_heat, "Figure_02c_GSDME_dependency_five_phenotypes_consensus",
          15, 9)

dep_delta <- dependency_contrasts[
  dependency_contrasts$node == "PYROPTOSIS", , drop = FALSE
]
dep_delta$comparison <- factor(dep_delta$comparison,
                               levels = rev(dependency_pairs$comparison))
p_gsdme_delta <- ggplot2::ggplot(
  dep_delta, ggplot2::aes(x = delta_endogenous_minus_ko_mean, y = comparison)
) +
  ggplot2::geom_vline(xintercept = 0, colour = "#777777") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = delta_endogenous_minus_ko_ci_low,
                 xmax = delta_endogenous_minus_ko_ci_high),
    orientation = "y", width = 0.15, na.rm = TRUE
  ) +
  ggplot2::geom_point(size = 3.5, colour = "#0072B2") +
  ggplot2::scale_x_continuous(limits = c(-1, 1),
                              labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "Loss of projected pyroptosis after GSDME knockout",
    subtitle = "DDR ON | endogenous GSDME minus matched intervention with GSDME KO",
    caption = "Intervals measure variability across simulation seeds; conditions used independent trajectory draws.",
    x = "Difference in pyroptosis frequency", y = NULL
  ) + plot_theme
save_plot(p_gsdme_delta, "Figure_02d_GSDME_dependency_contrast_consensus", 11, 5)

# C2: Pyroptosis endpoint with Monte Carlo 95% intervals.
pyro <- perturb_consensus[perturb_consensus$node == "PYROPTOSIS", ]
pyro <- pyro[order(pyro$activation_frequency_mean), ]
pyro$perturbation <- factor(pretty_condition(pyro$perturbation),
                            levels = pretty_condition(pyro$perturbation))
p2 <- ggplot2::ggplot(pyro, ggplot2::aes(activation_frequency_mean, perturbation)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = activation_frequency_ci_low, xmax = activation_frequency_ci_high),
    width = 0.18, linewidth = 0.6, colour = "#555555", na.rm = TRUE,
    orientation = "y"
  ) +
  ggplot2::geom_point(size = 3.2, colour = "#0072B2") +
  ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "Robustness of the GSDME-mediated pyroptosis endpoint",
    subtitle = "Points are seed means; bars are 95% t intervals across computational seeds",
    x = "Endpoint activation frequency", y = NULL
  ) + plot_theme
save_plot(p2, "Supplementary_Figure_S1_pyroptosis_seed_robustness", 11.5, 7.5)

# C3: Driver-set success and reproducibility.
driver_plot <- driver_consensus[order(driver_consensus$success_rate_mean), ]
driver_plot$intervention <- factor(driver_plot$intervention, levels = driver_plot$intervention)
p3 <- ggplot2::ggplot(driver_plot, ggplot2::aes(success_rate_mean, intervention,
                                                colour = selection_frequency)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = success_rate_ci_low, xmax = success_rate_ci_high),
    width = 0.18, linewidth = 0.6, colour = "#666666", na.rm = TRUE,
    orientation = "y"
  ) +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::geom_vline(xintercept = 0.95, linetype = "dashed", colour = "#B2182B") +
  ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_colour_gradient(low = "#D9D9D9", high = "#009E73", limits = c(0, 1),
                                 name = "Seed selection\nfrequency") +
  ggplot2::labs(
    title = "Consensus attractor-control interventions",
    subtitle = "Success intervals represent computational-seed variability; dashed line = 95%",
    x = "Mean control success rate", y = NULL
  ) + plot_theme
save_plot(p3, "Figure_03_minimum_driver_nodes_consensus", 12, 7.5)

# Figure 4: Mean node-detection rates across seeds.
coverage_plot <- coverage_consensus[as.logical(coverage_consensus$measured), ]
coverage_plot$node_label <- pretty_regulator(coverage_plot$node)
coverage_plot$node_label <- factor(coverage_plot$node_label,
                                   levels = unique(coverage_plot$node_label))
coverage_plot$experiment <- factor(coverage_plot$experiment,
                                   levels = rev(unique(coverage_plot$experiment)))
p4 <- ggplot2::ggplot(coverage_plot, ggplot2::aes(node_label, experiment,
                                                  fill = detection_rate_mean)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(detection_rate_mean, accuracy = 1)), size = 3.0
  ) +
  ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", limits = c(0, 1),
                               name = "Mean\ndetected") +
  ggplot2::labs(
    title = "Consensus single-cell detection of GSDME-network components",
    subtitle = "Detection rates are averaged across computational seed-specific cell samples",
    x = NULL, y = NULL
  ) + plot_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_plot(p4, "Figure_04_GEO_node_detection_consensus", 13, 5.8)

# Figure 5: Mean weighted observational concordance.
geo_heat <- concordance_consensus
geo_heat$experiment <- factor(geo_heat$experiment, levels = rev(unique(geo_heat$experiment)))
geo_heat$perturbation <- factor(pretty_condition(geo_heat$perturbation),
                                levels = pretty_condition(unique(geo_heat$perturbation)))
p5 <- ggplot2::ggplot(geo_heat, ggplot2::aes(perturbation, experiment,
                                             fill = weighted_concordance_score_mean)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", weighted_concordance_score_mean)), size = 3.1
  ) +
  ggplot2::scale_fill_gradient2(low = "#E66101", mid = "#F7F7F7", high = "#0571B0",
                                midpoint = 0.5, limits = c(0, 1),
                                name = "Weighted\nconcordance") +
  ggplot2::labs(
    title = "Consensus single-cell concordance with simulated endpoints",
    subtitle = "Observational state support only; not causal perturbation validation",
    x = NULL, y = NULL
  ) + plot_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 42, hjust = 1))
save_plot(p5, "Figure_05_GEO_weighted_perturbation_concordance_consensus", 14, 6.2)

# Assess observational agreement with the UNPERTURBED network separately.
# RNA-based concordance cannot be used to validate unperformed interventions.
reference_geo <- cross_consensus[
  grepl("^Reference:", cross_consensus$perturbation), , drop = FALSE
]
if (nrow(reference_geo)) {
  write_csv(reference_geo, "21_HCC_reference_model_observational_concordance.csv")
  reference_geo$dataset <- factor(reference_geo$accession,
                                  levels = c("GSE125449", "GSE189903"))
  p_reference_geo <- ggplot2::ggplot(
    reference_geo,
    ggplot2::aes(x = dataset, y = mean_weighted_concordance_mean)
  ) +
    ggplot2::geom_col(width = 0.57, fill = "#0072B2") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_weighted_concordance_ci_low,
                   ymax = mean_weighted_concordance_ci_high),
      width = 0.17, na.rm = TRUE
    ) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("%.2f", mean_weighted_concordance_mean)),
      vjust = -0.8, size = 4
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1.08)) +
    ggplot2::labs(
      title = "Observational HCC single-cell concordance with the DDR-ON model",
      subtitle = "Unperturbed network only | HCC-selected samples in both GEO series",
      caption = paste("Transcript-based agreement is observational and cannot",
                      "confirm CASP3 activity, GSDME cleavage or a cell-death phenotype."),
      x = NULL, y = "Reference network concordance") + plot_theme
  save_plot(p_reference_geo, "Figure_01b_HCC_observational_model_concordance", 10, 5)
}

# Figure 6: Independent Boolean outputs for the professor's fixed condition.
fate_plot <- fixed_phenotypes_consensus
fate_plot$fate_label <- factor(pretty_node(fate_plot$phenotype),
                               levels = pretty_node(c("PYROPTOSIS", "APOPTOSIS",
                                 "RESISTANCE", "PROLIFERATION", "CELL_CYCLE_ARREST")))
fate_plot$group_label <- paste(gsub("_", " ", fate_plot$experiment),
                               fate_plot$validation_group, sep = " / ")
fate_plot$condition <- factor(pretty_condition(fate_plot$condition),
                              levels = pretty_condition(c("Unperturbed",
                                  "miR-204-5p OE + GSDME KO")))
p6 <- ggplot2::ggplot(fate_plot,
                       ggplot2::aes(fate_label, group_label,
                                    fill = mean_probability_mean)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * mean_probability_mean)),
                     size = 2.9) +
  ggplot2::facet_wrap(~ condition, nrow = 1) +
  ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                               labels = scales::percent_format(accuracy = 1),
                               name = "Reachability") +
  ggplot2::labs(
    title = "Five independent model phenotypes across five computational seeds",
    subtitle = "DDR ON; fixed miR-204-5p E1 + GSDME KO versus unperturbed model",
    caption = paste0("Outputs can coexist and do not sum to 100%. Mean of seed-level fixed-point ",
                     "reachabilities; seed confidence intervals are in Table 14. ",
                     "The GEO samples were not experimentally perturbed."),
    x = NULL, y = NULL
  ) + plot_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_plot(p6, "Figure_06_GEO_projected_fate_consensus", 18, 10)

priority_nodes <- c("MIR204", "SIRT1", "MALAT1", "TP53", "BAX", "DFNA5", "CASP3")
priority_plot <- coverage_consensus[coverage_consensus$node %in% priority_nodes, ]
priority_labels <- c(MIR204 = "miR-204-5p (MIR204)", SIRT1 = "SIRT1",
                     MALAT1 = "MALAT1", TP53 = "p53 (TP53)", BAX = "BAX",
                     DFNA5 = "GSDME (DFNA5)", CASP3 = "CASP3")
priority_plot$node_label <- factor(unname(priority_labels[priority_plot$node]),
                                   levels = unname(priority_labels[priority_nodes]))
priority_plot$display_rate <- ifelse(as.logical(priority_plot$measured),
                                     priority_plot$detection_rate_mean, NA_real_)
priority_plot$label <- ifelse(as.logical(priority_plot$measured),
                              sprintf("%.0f%%", 100 * priority_plot$detection_rate_mean),
                              "Not measured")
p6b <- ggplot2::ggplot(priority_plot,
                        ggplot2::aes(node_label, experiment, fill = display_rate)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 3.0) +
  ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#3A9D84", limits = c(0, 1),
                               na.value = "grey88", labels = scales::percent_format(),
                               name = "RNA detected") +
  ggplot2::labs(
    title = "Detection of the seven priority network components",
    subtitle = "Mean RNA detection across computational seeds in selected HCC specimens",
    caption = "RNA detection cannot confirm p53 or CASP3 activation, miRNA activity, or GSDME cleavage.",
    x = NULL, y = NULL
  ) + plot_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_plot(p6b, "Figure_06b_seven_priority_nodes_RNA_detection_consensus", 15, 5.5)

# Figure 7: Stability of learned actions by step.
action_plot <- action_counts[action_counts$action != "No recorded step", ]
p7 <- ggplot2::ggplot(action_plot, ggplot2::aes(step, action,
                                                fill = frequency_all_seeds)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.55) +
  ggplot2::geom_text(ggplot2::aes(label = scales::percent(frequency_all_seeds, accuracy = 1)),
                     size = 3.1) +
  ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                               name = "Seed\nfrequency") +
  ggplot2::scale_x_continuous(breaks = seq_len(max_step)) +
  ggplot2::labs(
    title = "Stability of Q-learning-derived sequential interventions",
    subtitle = paste0("Pyroptosis target reached in ",
                      sum(rl_summary_all$reached_numeric, na.rm = TRUE),
                      " of ", length(folder_seeds), " computational seeds"),
    x = "Simulation step", y = NULL
  ) + plot_theme
save_plot(p7, "Figure_07_RL_action_stability_consensus", 12, 6.8)

# -------------------------- Reproducibility report ---------------------------

readme_lines <- c(
  "GSDME CONSENSUS ACROSS COMPUTATIONAL SEEDS",
  "==========================================",
  "",
  paste0("Seeds included: ", paste(folder_seeds, collapse = ", ")),
  paste0("Number of seeds: ", length(folder_seeds)),
  paste0("Pipeline version: ", unique(qc$pipeline_version)),
  "",
  "Interpretation:",
  "- Means, SDs and 95% t intervals quantify Monte Carlo/seed variability.",
  "- Seeds are computational replicates, not biological replicates.",
  "- GEO concordance is observational support and is not causal validation.",
  "- A driver/intervention should not be called robust from a single seed.",
  "- Biological causality requires experimental perturbation and rescue.",
  "",
  "Primary tables:",
  "01_run_quality_control.csv",
  "02_perturbation_consensus.csv",
  "04_driver_set_consensus.csv",
  "05_driver_action_frequency.csv",
  "07_GEO_perturbation_concordance_consensus.csv",
  "08_GEO_fate_consensus.csv",
  "10_RL_target_consensus.csv",
  "11_RL_action_stability.csv",
  "13_seed_level_phenotype_endpoints.csv",
  "15_GSDME_dependency_contrasts_by_seed.csv",
  "16_GSDME_dependency_contrasts_consensus.csv",
  "17_exact_stable_states_consensus_DDR_ON.csv",
  "18_convergence_consensus_DDR_ON.csv",
  "19_feedback_loop_consensus_DDR_ON.csv",
  "20_stable_states_and_convergence_DDR_ON.csv",
  "",
  "Figures:",
  "- Figure 01 is copied after verifying an identical network model.",
  "- Figures 02-07 summarize all included computational seeds.",
  "- Supplementary Figure S1 shows pyroptosis endpoint uncertainty.",
  "",
  "Publication note:",
  "Report the exact seeds, full-run parameters, pipeline version, and the",
  "distinction between computational robustness and biological evidence."
)
writeLines(readme_lines, file.path(CFG$out, "README_consensus_results.txt"))

capture.output(sessionInfo(), file = file.path(CFG$out, "sessionInfo_consensus.txt"))
manifest <- data.frame(
  file = list.files(CFG$out, recursive = TRUE, full.names = FALSE),
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
  consensus_script_version = "1.2.0",
  seeds = paste(folder_seeds, collapse = ";"),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(CFG$out, "output_manifest_consensus.csv"), row.names = FALSE)

log_message("Consensus analysis completed: ", CFG$out)
cat("\nCompleted successfully.\nConsensus results: ", CFG$out, "\n", sep = "")
