#!/usr/bin/env Rscript

# =============================================================================
# GSDME Systems-Oncology Boolean Network Pipeline
# =============================================================================
# Modules:
#   1. GINsim import, network reconstruction and attractor exploration
#   2. Sustained attractor control / minimal driver-node search
#   3. Probabilistic scRNA-seq binarisation and cell-specific simulations
#   4. GSE125449/GSE140228 observational validation with scope safeguards
#   5. Multi-omic soft-clipping for patient-specific network instances
#   6. Tabular Q-learning for sequential interventions
#   7. Publication figures (PDF + 600-dpi PNG) and CSV tables
#
# The script has a reproducible demo mode. Demo scRNA-seq and multi-omic data
# are simulated and MUST NOT be presented as experimental validation.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

PIPELINE_VERSION <- "1.1.0"

# ------------------------------- CLI -----------------------------------------

parse_cli <- function(args) {
  out <- list(
    model = "modelo_GINsim_GSDME_available.zginml",
    outdir = "GSDME_pipeline_results",
    scrna = NA_character_,
    geo = character(0),
    geo_dir = "GEO_scRNA_data",
    geo_platform = "droplet",
    geo_max_cells = 2000L,
    multiomics = NA_character_,
    seed = 204L,
    demo = TRUE,
    auto_install = TRUE,
    quick = FALSE
  )
  if (!length(args)) return(out)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      cat(paste0(
        "GSDME Systems-Oncology Boolean Network Pipeline\n\n",
        "Usage:\n",
        "  Rscript GSDME_systems_oncology_pipeline.R \\\n",
        "    --model modelo.zginml --out results [options]\n\n",
        "Options:\n",
        "  --scrna PATH        CSV/TSV/RDS matrix (genes x cells)\n",
        "  --geo ACCESSIONS    Comma-separated GEO accessions (supported: GSE125449,GSE140228)\n",
        "  --geo-dir PATH      Download/cache directory for GEO files (default GEO_scRNA_data)\n",
        "  --geo-platform X    GSE140228 platform: droplet, smartseq2 or all (default droplet)\n",
        "  --geo-max-cells N   Maximum analysed cells per GEO experiment (default 2000)\n",
        "  --multiomics PATH   Long CSV/TSV patient multi-omics table\n",
        "  --seed INTEGER      Random seed (default 204)\n",
        "  --no-demo           Require real input data\n",
        "  --no-install        Do not install missing CRAN packages\n",
        "  --quick             Smaller run for testing\n"
      ))
      quit(save = "no", status = 0)
    }
    if (key %in% c("--no-demo", "--no-install", "--quick")) {
      if (key == "--no-demo") out$demo <- FALSE
      if (key == "--no-install") out$auto_install <- FALSE
      if (key == "--quick") out$quick <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) stop("Missing value after ", key)
    value <- args[[i + 1L]]
    if (key == "--model") out$model <- value
    else if (key == "--out") out$outdir <- value
    else if (key == "--scrna") out$scrna <- value
    else if (key == "--geo") {
      out$geo <- unique(toupper(trimws(strsplit(value, ",", fixed = TRUE)[[1]])))
      out$geo <- out$geo[nzchar(out$geo)]
    }
    else if (key == "--geo-dir") out$geo_dir <- value
    else if (key == "--geo-platform") out$geo_platform <- tolower(value)
    else if (key == "--geo-max-cells") out$geo_max_cells <- as.integer(value)
    else if (key == "--multiomics") out$multiomics <- value
    else if (key == "--seed") out$seed <- as.integer(value)
    else stop("Unknown argument: ", key)
    i <- i + 2L
  }
  out
}

CFG <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!CFG$geo_platform %in% c("droplet", "smartseq2", "all")) {
  stop("--geo-platform must be droplet, smartseq2 or all")
}
if (is.na(CFG$geo_max_cells) || CFG$geo_max_cells < 20L) {
  stop("--geo-max-cells must be an integer >= 20")
}
set.seed(CFG$seed)

# ---------------------------- Dependencies -----------------------------------

required_packages <- c("xml2", "igraph", "ggplot2")
if (length(CFG$geo)) required_packages <- unique(c(required_packages, "Matrix", "data.table"))
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) && CFG$auto_install) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    ". Install them with install.packages()."
  )
}

library(xml2)
library(igraph)
library(ggplot2)

# ------------------------------ Output ---------------------------------------

dir.create(CFG$outdir, recursive = TRUE, showWarnings = FALSE)
FIG_DIR <- file.path(CFG$outdir, "figures")
TAB_DIR <- file.path(CFG$outdir, "tables")
LOG_DIR <- file.path(CFG$outdir, "logs")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(LOG_DIR, "pipeline.log"), append = TRUE)
}

save_publication_plot <- function(plot, stem, width = 11, height = 8) {
  pdf_path <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIG_DIR, paste0(stem, "_600dpi.png"))
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height,
                  units = "in", device = pdf_device, bg = "white")
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  units = "in", dpi = 600, bg = "white", limitsize = FALSE)
  invisible(c(pdf = pdf_path, png = png_path))
}

write_table <- function(x, filename) {
  utils::write.csv(x, file.path(TAB_DIR, filename), row.names = FALSE, na = "")
}

# --------------------------- GINsim parsing ----------------------------------

read_ginsim_model <- function(path) {
  if (!file.exists(path)) stop("GINsim file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  xml_path <- path
  cleanup <- NULL
  if (ext == "zginml") {
    entries <- utils::unzip(path, list = TRUE)$Name
    candidate <- entries[grepl("regulatoryGraph\\.ginml$", entries)]
    if (!length(candidate)) stop("No regulatoryGraph.ginml found inside ", path)
    td <- tempfile("ginsim_")
    dir.create(td)
    utils::unzip(path, files = candidate[[1]], exdir = td)
    xml_path <- file.path(td, candidate[[1]])
    cleanup <- td
  }
  on.exit(if (!is.null(cleanup)) unlink(cleanup, recursive = TRUE), add = TRUE)
  doc <- xml2::read_xml(xml_path, options = c("RECOVER", "NONET", "NOBLANKS"))
  graph_node <- xml2::xml_find_first(doc, ".//graph")
  node_xml <- xml2::xml_find_all(graph_node, "./node")
  edge_xml <- xml2::xml_find_all(graph_node, "./edge")

  nodes <- data.frame(
    id = xml2::xml_attr(node_xml, "id"),
    label = ifelse(is.na(xml2::xml_attr(node_xml, "name")),
                   xml2::xml_attr(node_xml, "id"), xml2::xml_attr(node_xml, "name")),
    maxvalue = as.integer(xml2::xml_attr(node_xml, "maxvalue")),
    input = xml2::xml_attr(node_xml, "input") == "true",
    rule = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(node_xml)) {
    exp_node <- xml2::xml_find_first(node_xml[[i]], "./value/exp")
    if (!inherits(exp_node, "xml_missing")) nodes$rule[[i]] <- xml2::xml_attr(exp_node, "str")
  }
  nodes$maxvalue[is.na(nodes$maxvalue)] <- 1L

  edges <- data.frame(
    id = xml2::xml_attr(edge_xml, "id"),
    from = xml2::xml_attr(edge_xml, "from"),
    to = xml2::xml_attr(edge_xml, "to"),
    sign = xml2::xml_attr(edge_xml, "sign"),
    stringsAsFactors = FALSE
  )
  model <- list(
    id = xml2::xml_attr(graph_node, "id"),
    nodes = nodes,
    edges = edges,
    node_ids = nodes$id,
    input_ids = nodes$id[nodes$input]
  )
  validate_model(model)
  compiled <- lapply(nodes$rule, function(rule) {
    if (is.na(rule) || !nzchar(rule)) NULL else parse(text = rule)[[1]]
  })
  names(compiled) <- nodes$id
  model$compiled_rules <- compiled
  model
}

validate_model <- function(model) {
  stopifnot(!anyDuplicated(model$nodes$id))
  unknown_edges <- setdiff(unique(c(model$edges$from, model$edges$to)), model$nodes$id)
  if (length(unknown_edges)) stop("Edges refer to unknown nodes: ", paste(unknown_edges, collapse = ", "))
  for (i in which(!is.na(model$nodes$rule))) {
    vars <- all.vars(parse(text = model$nodes$rule[[i]]))
    unknown <- setdiff(vars, model$nodes$id)
    if (length(unknown)) {
      stop("Rule for ", model$nodes$id[[i]], " refers to unknown node(s): ",
           paste(unknown, collapse = ", "))
    }
  }
  invisible(TRUE)
}

MODEL <- read_ginsim_model(CFG$model)
log_message("Loaded model '", MODEL$id, "': ", nrow(MODEL$nodes), " nodes, ",
            nrow(MODEL$edges), " edges, input(s): ", paste(MODEL$input_ids, collapse = ", "))
write_table(MODEL$nodes, "01_model_nodes_and_rules.csv")
write_table(MODEL$edges, "02_model_edges.csv")

# ------------------------ Boolean simulation core ----------------------------

empty_state <- function(model, value = 0L) {
  stats::setNames(rep(as.integer(value), length(model$node_ids)), model$node_ids)
}

normalise_state <- function(state, model) {
  out <- empty_state(model)
  common <- intersect(names(state), model$node_ids)
  out[common] <- as.integer(state[common] > 0)
  out
}

apply_clamp <- function(state, clamp = integer(0)) {
  if (length(clamp)) state[names(clamp)] <- as.integer(clamp)
  state
}

evaluate_rule <- function(rule_expression, state_environment) {
  if (is.null(rule_expression)) return(0L)
  as.integer(isTRUE(eval(rule_expression, envir = state_environment)))
}

logical_targets <- function(state, model, clamp = integer(0)) {
  state <- apply_clamp(normalise_state(state, model), clamp)
  target <- state
  state_environment <- list2env(as.list(as.logical(state)), parent = baseenv())
  for (i in seq_len(nrow(model$nodes))) {
    node <- model$nodes$id[[i]]
    if (node %in% names(clamp)) {
      target[[node]] <- as.integer(clamp[[node]])
    } else if (model$nodes$input[[i]]) {
      target[[node]] <- state[[node]]
    } else {
      target[[node]] <- evaluate_rule(model$compiled_rules[[node]], state_environment)
    }
  }
  apply_clamp(target, clamp)
}

step_boolean <- function(state, model, clamp = integer(0),
                         mode = c("synchronous", "asynchronous")) {
  mode <- match.arg(mode)
  state <- apply_clamp(normalise_state(state, model), clamp)
  target <- logical_targets(state, model, clamp)
  if (mode == "synchronous") return(target)
  unstable <- setdiff(names(state)[state != target], names(clamp))
  if (!length(unstable)) return(state)
  selected <- sample(unstable, 1L)
  state[[selected]] <- target[[selected]]
  apply_clamp(state, clamp)
}

encode_state <- function(state, model) paste0(state[model$node_ids], collapse = "")

simulate_to_attractor <- function(initial_state, model, clamp = integer(0),
                                  mode = "asynchronous", max_steps = 300L,
                                  return_trajectory = FALSE) {
  state <- apply_clamp(normalise_state(initial_state, model), clamp)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  trajectory <- if (return_trajectory) list(state) else NULL
  for (step in seq_len(max_steps)) {
    key <- encode_state(state, model)
    if (exists(key, seen, inherits = FALSE)) {
      first_seen <- get(key, seen, inherits = FALSE)
      return(list(state = state, steps = step - 1L, stable = FALSE,
                  cycle_length = step - first_seen,
                  trajectory = trajectory))
    }
    assign(key, step, seen)
    next_state <- step_boolean(state, model, clamp, mode)
    if (identical(unname(next_state), unname(state))) {
      return(list(state = state, steps = step - 1L, stable = TRUE,
                  cycle_length = 1L, trajectory = trajectory))
    }
    state <- next_state
    if (return_trajectory) trajectory[[length(trajectory) + 1L]] <- state
  }
  list(state = state, steps = max_steps, stable = FALSE,
       cycle_length = NA_integer_, trajectory = trajectory)
}

phenotype_nodes <- intersect(
  c("PYROPTOSIS_GSDME", "APOPTOSIS", "SURVIVAL", "RESISTANCE",
    "PROLIFERATION", "CELL_CYCLE_ARREST"),
  MODEL$node_ids
)

classify_fate <- function(state) {
  is_on <- function(x) x %in% names(state) && state[[x]] == 1L
  if (is_on("PYROPTOSIS_GSDME")) return("Pyroptosis")
  if (is_on("APOPTOSIS")) return("Apoptosis")
  if (is_on("RESISTANCE")) return("Resistance")
  if (is_on("PROLIFERATION")) return("Proliferation")
  if (is_on("CELL_CYCLE_ARREST")) return("Cell-cycle arrest")
  if (is_on("SURVIVAL")) return("Survival")
  "Other"
}

make_reference_state <- function(model, type = c("malignant", "death_primed")) {
  type <- match.arg(type)
  s <- empty_state(model)
  s[intersect(c("DDR_fixed_ON", "GSDME_availability", "ROS"), names(s))] <- 1L
  if (type == "malignant") {
    s[intersect(c("lncRNA_MALAT1", "SIRT1", "PGC1A", "SURVIVAL",
                  "RESISTANCE", "CyclinD_CDK46", "E2F1"), names(s))] <- 1L
  } else {
    s[intersect(c("miR_204_5p", "p53_ACTIVE", "PUMA", "MITO_DAMAGE", "BAX",
                  "MOMP", "CYTOCHROME_C", "CASP9", "CASP3", "GSDME_N",
                  "PYROPTOSIS_GSDME", "p21", "RB1"), names(s))] <- 1L
  }
  s
}

BASE_CLAMP <- c(DDR_fixed_ON = 1L)
BASE_CLAMP <- BASE_CLAMP[names(BASE_CLAMP) %in% MODEL$node_ids]
MALIGNANT_STATE <- simulate_to_attractor(
  make_reference_state(MODEL, "malignant"), MODEL,
  clamp = c(BASE_CLAMP, GSDME_availability = 1L), mode = "synchronous"
)$state
DEATH_STATE <- simulate_to_attractor(
  make_reference_state(MODEL, "death_primed"), MODEL,
  clamp = c(BASE_CLAMP, GSDME_availability = 1L), mode = "synchronous"
)$state

reference_table <- rbind(
  data.frame(reference = "Malignant/resistant", t(MALIGNANT_STATE), check.names = FALSE),
  data.frame(reference = "Death/pyroptosis", t(DEATH_STATE), check.names = FALSE)
)
write_table(reference_table, "03_reference_attractors.csv")

# ----------------------- Publication network plot ----------------------------

pretty_labels <- c(
  GSDME_availability = "GSDME availability",
  GSDME_N = "GSDME-N",
  CASP3 = "CASP3", CASP9 = "CASP9", CYTOCHROME_C = "Cytochrome c",
  MOMP = "MOMP", BAX = "BAX", BCL2 = "BCL2", PUMA = "PUMA",
  p53_ACTIVE = "Active p53", miR_204_5p = "miR-204-5p",
  lncRNA_MALAT1 = "MALAT1", ROS = "ROS",
  MITO_DAMAGE = "Mitochondrial damage", PGC1A = "PGC-1alpha",
  p21 = "p21", CyclinD_CDK46 = "Cyclin D-CDK4/6", RB1 = "RB1",
  E2F1 = "E2F1", DDR_fixed_ON = "DDR (fixed ON)",
  PYROPTOSIS_GSDME = "GSDME-mediated\npyroptosis", SURVIVAL = "Survival",
  RESISTANCE = "Resistance", PROLIFERATION = "Proliferation",
  APOPTOSIS = "Apoptosis", CELL_CYCLE_ARREST = "Cell-cycle arrest",
  SIRT1 = "SIRT1"
)

manual_layout <- data.frame(
  id = c("GSDME_availability", "GSDME_N", "PYROPTOSIS_GSDME", "APOPTOSIS",
         "CASP3", "CASP9", "CYTOCHROME_C", "MOMP", "BAX", "PUMA", "BCL2",
         "lncRNA_MALAT1", "miR_204_5p", "SIRT1", "PGC1A", "MITO_DAMAGE",
         "ROS", "DDR_fixed_ON", "p53_ACTIVE", "p21", "CyclinD_CDK46", "RB1",
         "E2F1", "SURVIVAL", "RESISTANCE", "PROLIFERATION", "CELL_CYCLE_ARREST"),
  x = c(5,5,2.7,7.3,5,5,5,5,5,3.8,2.6,7.5,7.5,7.5,7.5,6.5,6.5,9.3,9.3,
        11,11,11,11,9.5,9.5,7.2,4.1),
  y = c(12,11,10.4,10.4,9.7,8.6,7.5,6.4,5.3,4.4,5.3,10.1,9.0,7.9,6.8,5.3,
        4.0,5.3,7.8,7.0,5.9,4.8,3.7,9.7,2.2,2.2,2.2),
  stringsAsFactors = FALSE
)

plot_network <- function(model) {
  layout <- merge(data.frame(id = model$node_ids), manual_layout, by = "id", all.x = TRUE)
  missing <- which(is.na(layout$x) | is.na(layout$y))
  if (length(missing)) {
    fallback <- igraph::layout_with_fr(
      igraph::graph_from_data_frame(model$edges[, c("from", "to")],
                                    vertices = data.frame(name = model$nodes$id))
    )
    rownames(fallback) <- model$nodes$id
    layout$x[missing] <- scales::rescale(fallback[layout$id[missing], 1], c(1, 12))
    layout$y[missing] <- scales::rescale(fallback[layout$id[missing], 2], c(1, 12))
  }
  layout$label <- ifelse(layout$id %in% names(pretty_labels),
                         unname(pretty_labels[layout$id]), layout$id)
  outputs <- c("PYROPTOSIS_GSDME", "APOPTOSIS", "SURVIVAL", "RESISTANCE",
               "PROLIFERATION", "CELL_CYCLE_ARREST")
  layout$class <- ifelse(layout$id %in% model$input_ids, "Input",
                         ifelse(layout$id == "DDR_fixed_ON", "Fixed condition",
                                ifelse(layout$id %in% outputs, "Phenotype", "Regulator")))
  edges <- merge(model$edges, layout[, c("id", "x", "y")], by.x = "from", by.y = "id")
  names(edges)[names(edges) %in% c("x", "y")] <- c("x", "y")
  edges <- merge(edges, layout[, c("id", "x", "y")], by.x = "to", by.y = "id",
                 suffixes = c("", "end"))
  edges$edge_type <- ifelse(edges$sign == "negative", "Inhibition", "Activation")

  pos <- edges[edges$edge_type == "Activation", , drop = FALSE]
  neg <- edges[edges$edge_type == "Inhibition", , drop = FALSE]
  p <- ggplot() +
    geom_curve(data = pos,
               aes(x = x, y = y, xend = xend, yend = yend, colour = edge_type),
               curvature = 0.04, linewidth = 0.55,
               arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed")) +
    geom_curve(data = neg,
               aes(x = x, y = y, xend = xend, yend = yend, colour = edge_type),
               curvature = -0.04, linewidth = 0.65,
               arrow = grid::arrow(length = grid::unit(2.7, "mm"), type = "open", angle = 90)) +
    geom_label(data = layout,
               aes(x = x, y = y, label = label, fill = class),
               size = 3.0, label.size = 0.35, label.padding = grid::unit(0.17, "lines"),
               colour = "#17212B", fontface = "bold") +
    scale_fill_manual(values = c(Input = "#7FE7E7", `Fixed condition` = "#F1C75B",
                                 Phenotype = "#D7D7D7", Regulator = "white")) +
    scale_colour_manual(values = c(Activation = "#17823B", Inhibition = "#C43131")) +
    coord_equal(clip = "off") +
    labs(title = "Logical network of GSDME-mediated cell-fate control",
         subtitle = "Green arrows: activation; red T-ended edges: inhibition",
         fill = NULL, colour = NULL) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "bottom",
          plot.margin = margin(12, 18, 12, 18))
  p
}

network_plot <- plot_network(MODEL)
save_publication_plot(network_plot, "Figure_01_GSDME_logical_network", 13, 10)

# -------------------- Literature-grounded perturbations ----------------------

scientific_perturbations <- list(
  `Reference: DDR ON, GSDME ON` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L),
  `GSDME KO` = c(DDR_fixed_ON = 1L, GSDME_availability = 0L),
  `CASP3 KO` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L, CASP3 = 0L),
  `miR-204-5p OE` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L, miR_204_5p = 1L),
  `miR-204-5p KO` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L, miR_204_5p = 0L),
  `MALAT1 OE` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L, lncRNA_MALAT1 = 1L),
  `MALAT1 KO` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L, lncRNA_MALAT1 = 0L),
  `SIRT1 OE` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L, SIRT1 = 1L),
  `SIRT1 KO` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L, SIRT1 = 0L),
  `BCL2 KO` = c(DDR_fixed_ON = 1L, GSDME_availability = 1L, BCL2 = 0L)
)

perturbation_evidence <- data.frame(
  perturbation = names(scientific_perturbations),
  expected_direction = c(
    "Reference condition", "Pyroptosis-to-apoptosis switch when CASP3 is active",
    "Reduced GSDME-N and pyroptosis", "SIRT1 suppression; anti-tumour direction",
    "SIRT1 release; pro-survival direction", "miR-204 sequestration and SIRT1 release",
    "miR-204 release and SIRT1 suppression", "Pro-survival/resistance direction",
    "Reduced survival; mitochondrial death sensitisation", "BAX disinhibition"
  ),
  evidence_scope = c(
    "Model reference", "Direct experimental evidence; cancer, non-HCC",
    "Direct experimental evidence in HepG2", "Direct experimental evidence in HCC",
    "Mechanistic inverse inferred from HCC evidence", "Direct experimental evidence in HCC",
    "Mechanistic inverse inferred from HCC evidence", "Axis-supported model perturbation",
    "Axis-supported model perturbation", "Canonical mitochondrial apoptosis mechanism"
  ),
  reference = c(
    "Model assumption", "PMID:28459430; DOI:10.1038/nature22393",
    "PMID:35747157; DOI:10.3892/etm.2022.11383",
    "PMID:27748572; DOI:10.1002/cbf.3223",
    "PMID:27748572; DOI:10.1002/cbf.3223",
    "PMID:28720061; DOI:10.1177/1010428317718135",
    "PMID:28720061; DOI:10.1177/1010428317718135",
    "PMID:27748572; PMID:28720061", "PMID:27748572; PMID:28720061",
    "Network mechanism; experimental validation required in selected HCC system"
  ),
  stringsAsFactors = FALSE
)
write_table(perturbation_evidence, "04_scientific_perturbation_evidence.csv")

random_initial_state <- function(model, clamp = integer(0), p = 0.5) {
  s <- stats::setNames(stats::rbinom(length(model$node_ids), 1L, p), model$node_ids)
  apply_clamp(s, clamp)
}

screen_perturbations <- function(model, perturbations, n_trajectories = 160L,
                                 max_steps = 350L) {
  rows <- vector("list", length(perturbations))
  for (i in seq_along(perturbations)) {
    clamp <- perturbations[[i]]
    clamp <- clamp[names(clamp) %in% model$node_ids]
    states <- matrix(0L, nrow = n_trajectories, ncol = length(model$node_ids),
                     dimnames = list(NULL, model$node_ids))
    stable <- logical(n_trajectories)
    for (j in seq_len(n_trajectories)) {
      initial <- random_initial_state(model, clamp)
      result <- simulate_to_attractor(initial, model, clamp, "asynchronous", max_steps)
      states[j, ] <- result$state
      stable[j] <- result$stable
    }
    means <- colMeans(states)
    rows[[i]] <- data.frame(
      perturbation = names(perturbations)[[i]],
      node = names(means), activation_frequency = as.numeric(means),
      stable_fraction = mean(stable), stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

n_perturb <- if (CFG$quick) 40L else 160L
perturbation_results <- screen_perturbations(MODEL, scientific_perturbations, n_perturb)
write_table(perturbation_results, "05_perturbation_activation_frequencies.csv")

plot_heatmap_long <- function(data, row_var, col_var, value_var, title, subtitle = NULL,
                              low = "#F7FBFF", high = "#B2182B") {
  data[[row_var]] <- factor(data[[row_var]], levels = unique(data[[row_var]]))
  data[[col_var]] <- factor(data[[col_var]], levels = unique(data[[col_var]]))
  ggplot(data, aes(x = .data[[col_var]], y = .data[[row_var]], fill = .data[[value_var]])) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", .data[[value_var]])), size = 2.6) +
    scale_fill_gradient(low = low, high = high, limits = c(0, 1), name = "Frequency") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"), legend.position = "right")
}

perturbation_pheno <- perturbation_results[perturbation_results$node %in% phenotype_nodes, ]
perturbation_pheno$node <- factor(perturbation_pheno$node, levels = phenotype_nodes)
p_perturb <- plot_heatmap_long(
  perturbation_pheno, "perturbation", "node", "activation_frequency",
  "In silico perturbation screen",
  "Asynchronous endpoint activation frequencies; literature evidence is reported separately"
)
save_publication_plot(p_perturb, "Figure_02_validated_perturbation_heatmap", 11.5, 7.5)

# ---------------------- Structural control / drivers -------------------------

structural_ranking <- function(model) {
  g <- igraph::graph_from_data_frame(model$edges[, c("from", "to")], directed = TRUE,
                                     vertices = data.frame(name = model$nodes$id))
  scc <- igraph::components(g, mode = "strong")
  scc_size <- scc$csize[scc$membership]
  cyc <- as.integer(scc_size > 1L)
  data.frame(
    node = igraph::V(g)$name,
    indegree = igraph::degree(g, mode = "in"),
    outdegree = igraph::degree(g, mode = "out"),
    betweenness = igraph::betweenness(g, directed = TRUE, normalized = TRUE),
    pagerank = igraph::page_rank(g, directed = TRUE)$vector,
    in_feedback_scc = cyc,
    stringsAsFactors = FALSE
  )
}

structural_scores <- structural_ranking(MODEL)
structural_scores$score <- with(structural_scores,
  2 * in_feedback_scc + scale(outdegree)[, 1] + scale(pagerank)[, 1]
)
structural_scores <- structural_scores[order(-structural_scores$score), ]
write_table(structural_scores, "06_structural_control_ranking.csv")

candidate_actions <- c(
  lncRNA_MALAT1 = 0L, miR_204_5p = 1L, SIRT1 = 0L, PGC1A = 0L,
  BCL2 = 0L, BAX = 1L, p53_ACTIVE = 1L, CASP3 = 1L
)
candidate_actions <- candidate_actions[names(candidate_actions) %in% MODEL$node_ids]

control_success <- function(model, clamp, start_states, target_node = "PYROPTOSIS_GSDME",
                            forbidden = c("RESISTANCE", "PROLIFERATION"),
                            n_rep = 30L, max_steps = 350L) {
  success <- logical(length(start_states) * n_rep)
  k <- 1L
  for (s in start_states) {
    for (r in seq_len(n_rep)) {
      result <- simulate_to_attractor(s, model, clamp, "asynchronous", max_steps)
      target_ok <- target_node %in% names(result$state) && result$state[[target_node]] == 1L
      forbid_ok <- !any(result$state[intersect(forbidden, names(result$state))] == 1L)
      success[[k]] <- target_ok && forbid_ok
      k <- k + 1L
    }
  }
  mean(success)
}

jitter_state <- function(state, probability = 0.08, protected = names(BASE_CLAMP)) {
  s <- state
  eligible <- setdiff(names(s), protected)
  flip <- eligible[stats::runif(length(eligible)) < probability]
  s[flip] <- 1L - s[flip]
  s
}

find_minimum_driver_sets <- function(model, candidate_actions, max_size = 3L,
                                     robustness = 0.95, n_start = 12L, n_rep = 12L) {
  starts <- c(list(MALIGNANT_STATE),
              replicate(n_start - 1L, jitter_state(MALIGNANT_STATE), simplify = FALSE))
  all_results <- list()
  idx <- 1L
  minimal_size <- NA_integer_
  for (size in seq_len(min(max_size, length(candidate_actions)))) {
    combos <- utils::combn(names(candidate_actions), size, simplify = FALSE)
    for (nodes in combos) {
      clamp <- c(BASE_CLAMP, GSDME_availability = 1L, candidate_actions[nodes])
      clamp <- clamp[!duplicated(names(clamp), fromLast = TRUE)]
      rate <- control_success(model, clamp, starts, n_rep = n_rep)
      all_results[[idx]] <- data.frame(
        size = size,
        intervention = paste0(nodes, "=", candidate_actions[nodes], collapse = " + "),
        success_rate = rate,
        robust = rate >= robustness,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
    current <- do.call(rbind, all_results)
    if (any(current$size == size & current$robust)) {
      minimal_size <- size
      break
    }
  }
  result <- do.call(rbind, all_results)
  attr(result, "minimum_size") <- minimal_size
  result[order(result$size, -result$success_rate, result$intervention), ]
}

driver_results <- find_minimum_driver_sets(
  MODEL, candidate_actions,
  max_size = if (CFG$quick) 2L else 3L,
  n_start = if (CFG$quick) 5L else 12L,
  n_rep = if (CFG$quick) 4L else 12L
)
write_table(driver_results, "07_minimum_driver_set_search.csv")

top_driver <- head(driver_results[order(-driver_results$success_rate, driver_results$size), ], 20)
top_driver$intervention <- factor(top_driver$intervention,
                                  levels = rev(top_driver$intervention))
p_driver <- ggplot(top_driver, aes(x = success_rate, y = intervention, fill = factor(size))) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0.95, linetype = 2, colour = "#B2182B") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2", name = "Driver-set size") +
  labs(title = "Candidate attractor-control interventions",
       subtitle = "Success from resistant and locally perturbed initial states; dashed line = 95%",
       x = "Control success rate", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))
save_publication_plot(p_driver, "Figure_03_minimum_driver_nodes", 10.5, 7.5)

# ---------------- Probabilistic scRNA-seq integration ------------------------

gene_map <- list(
  GSDME_availability = c("GSDME", "DFNA5"), CASP3 = "CASP3", CASP9 = "CASP9",
  CYTOCHROME_C = c("CYCS", "CYTOCHROME_C"), BAX = "BAX", BCL2 = "BCL2",
  PUMA = c("BBC3", "PUMA"), p53_ACTIVE = c("TP53", "P53"),
  miR_204_5p = c("MIR204", "MIR204-5P", "MIR_204_5P"),
  lncRNA_MALAT1 = "MALAT1", PGC1A = c("PPARGC1A", "PGC1A"),
  p21 = c("CDKN1A", "P21"), CyclinD_CDK46 = c("CCND1", "CDK4", "CDK6"),
  RB1 = "RB1", E2F1 = "E2F1", SIRT1 = "SIRT1"
)

load_expression_matrix <- function(path) {
  if (!file.exists(path)) stop("Expression file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    x <- readRDS(path)
    x <- as.matrix(x)
  } else {
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    x <- utils::read.table(path, header = TRUE, row.names = 1, sep = sep,
                           check.names = FALSE, comment.char = "", quote = "\"")
    x <- as.matrix(x)
  }
  storage.mode(x) <- "numeric"
  if (any(x < 0, na.rm = TRUE)) stop("scRNA-seq counts must be non-negative")
  x
}

make_demo_scrna <- function(n_cells = 240L, seed = 204L) {
  set.seed(seed)
  genes <- unique(unlist(gene_map))
  cell_group <- rep(c("Resistant-like", "Pyroptosis-primed", "Apoptosis-primed"),
                    length.out = n_cells)
  mat <- matrix(0, nrow = length(genes), ncol = n_cells,
                dimnames = list(genes, paste0("Cell_", seq_len(n_cells))))
  baseline <- stats::runif(length(genes), 0.2, 1.1)
  names(baseline) <- genes
  for (j in seq_len(n_cells)) {
    mu <- exp(baseline)
    names(mu) <- genes
    if (cell_group[[j]] == "Resistant-like") {
      mu[intersect(c("MALAT1", "SIRT1", "PPARGC1A", "BCL2", "E2F1"), genes)] <- 9
      mu[intersect(c("MIR204", "BBC3", "BAX", "CASP3", "GSDME"), genes)] <- 0.5
    } else if (cell_group[[j]] == "Pyroptosis-primed") {
      mu[intersect(c("MIR204", "BBC3", "BAX", "CASP9", "CASP3", "GSDME"), genes)] <- 10
      mu[intersect(c("MALAT1", "SIRT1", "BCL2"), genes)] <- 0.4
    } else {
      mu[intersect(c("MIR204", "BBC3", "BAX", "CASP9", "CASP3"), genes)] <- 8
      mu[intersect(c("GSDME", "DFNA5"), genes)] <- 0.2
      mu[intersect(c("MALAT1", "SIRT1"), genes)] <- 0.5
    }
    mat[, j] <- stats::rnbinom(length(genes), mu = mu, size = 1.4)
  }
  attr(mat, "demo_group") <- cell_group
  mat
}

fit_two_gaussian_posterior <- function(x, max_iter = 100L, tol = 1e-6) {
  x <- as.numeric(x)
  if (length(unique(x)) < 3L || stats::sd(x) < 1e-8) {
    return(ifelse(x > stats::median(x), 0.9, ifelse(x == 0, 0.05, 0.5)))
  }
  mu <- as.numeric(stats::quantile(x, c(0.25, 0.75), names = FALSE))
  sigma <- rep(max(stats::sd(x), 0.2), 2L)
  pi_k <- c(0.5, 0.5)
  old_ll <- -Inf
  for (iter in seq_len(max_iter)) {
    dens <- cbind(
      pi_k[[1]] * stats::dnorm(x, mu[[1]], sigma[[1]]),
      pi_k[[2]] * stats::dnorm(x, mu[[2]], sigma[[2]])
    )
    denom <- rowSums(dens) + 1e-300
    z <- dens / denom
    nk <- colSums(z) + 1e-8
    pi_k <- nk / length(x)
    mu <- colSums(z * x) / nk
    sigma <- sqrt(pmax(colSums(z * (x - rep(mu, each = length(x)))^2) / nk, 0.05^2))
    ll <- sum(log(denom))
    if (abs(ll - old_ll) < tol) break
    old_ll <- ll
  }
  high <- which.max(mu)
  pmin(pmax(z[, high], 0.001), 0.999)
}

probabilistic_binarise_scrna <- function(counts, gene_map) {
  lib <- attr(counts, "library_size", exact = TRUE)
  if (is.null(lib) || length(lib) != ncol(counts)) lib <- colSums(counts)
  lib[lib <= 0] <- 1
  log_norm <- log1p(t(t(counts) / lib * 1e4))
  gene_post <- matrix(NA_real_, nrow = nrow(log_norm), ncol = ncol(log_norm),
                      dimnames = dimnames(log_norm))
  for (i in seq_len(nrow(log_norm))) {
    gene_post[i, ] <- fit_two_gaussian_posterior(log_norm[i, ])
  }
  node_post <- matrix(NA_real_, nrow = length(gene_map), ncol = ncol(counts),
                      dimnames = list(names(gene_map), colnames(counts)))
  upper_rows <- toupper(rownames(gene_post))
  for (node in names(gene_map)) {
    idx <- which(upper_rows %in% toupper(gene_map[[node]]))
    if (length(idx) == 1L) node_post[node, ] <- gene_post[idx, ]
    if (length(idx) > 1L) node_post[node, ] <- apply(gene_post[idx, , drop = FALSE], 2, max)
  }
  list(node_probability = node_post, gene_probability = gene_post,
       log_normalised = log_norm)
}

simulate_single_cells <- function(node_probability, model, n_draws = 5L,
                                  max_steps = 300L) {
  fate_levels <- c("Pyroptosis", "Apoptosis", "Resistance", "Proliferation",
                   "Cell-cycle arrest", "Survival", "Other")
  out <- matrix(0, nrow = ncol(node_probability), ncol = length(fate_levels),
                dimnames = list(colnames(node_probability), fate_levels))
  for (cell in seq_len(ncol(node_probability))) {
    fate <- character(n_draws)
    for (draw in seq_len(n_draws)) {
      s <- random_initial_state(model, BASE_CLAMP)
      for (node in intersect(rownames(node_probability), model$node_ids)) {
        p <- node_probability[node, cell]
        if (!is.na(p)) s[[node]] <- stats::rbinom(1L, 1L, p)
      }
      result <- simulate_to_attractor(s, model, BASE_CLAMP, "asynchronous", max_steps)
      fate[[draw]] <- classify_fate(result$state)
    }
    out[cell, ] <- tabulate(match(fate, fate_levels), nbins = length(fate_levels)) / n_draws
  }
  out
}

run_generic_scrna <- !is.na(CFG$scrna) || CFG$demo
if (run_generic_scrna) {
  if (!is.na(CFG$scrna)) {
    scrna_counts <- load_expression_matrix(CFG$scrna)
    scrna_source <- "User-provided"
  } else {
    scrna_counts <- make_demo_scrna(if (CFG$quick) 60L else 240L, CFG$seed)
    scrna_source <- "SIMULATED DEMO DATA"
  }
  scrna_binary <- probabilistic_binarise_scrna(scrna_counts, gene_map)
  scrna_fates <- simulate_single_cells(
    scrna_binary$node_probability, MODEL,
    n_draws = if (CFG$quick) 2L else 5L
  )

  scrna_fate_table <- data.frame(cell = rownames(scrna_fates), scrna_fates,
                                 check.names = FALSE)
  scrna_fate_table$data_source <- scrna_source
  write_table(scrna_fate_table, "08_scrna_cell_fate_probabilities.csv")

  dominant_fate <- colnames(scrna_fates)[max.col(scrna_fates, ties.method = "first")]
  ord <- order(dominant_fate, -apply(scrna_fates, 1, max))
  show_cells <- ord[seq_len(min(length(ord), 120L))]
  scrna_long <- do.call(rbind, lapply(seq_along(show_cells), function(i) {
    cell <- rownames(scrna_fates)[show_cells[[i]]]
    data.frame(cell_order = i, cell = cell, fate = colnames(scrna_fates),
               probability = as.numeric(scrna_fates[cell, ]), stringsAsFactors = FALSE)
  }))
  p_scrna <- ggplot(scrna_long, aes(x = cell_order, y = fate, fill = probability)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("#F7FBFF", "#6BAED6", "#08306B"), limits = c(0, 1)) +
    labs(title = "Single-cell-informed Boolean outcomes",
         subtitle = paste0(scrna_source, "; posterior sampling followed by asynchronous simulation"),
         x = "Cells ordered by dominant simulated fate", y = NULL, fill = "Probability") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(face = "bold"))
  save_publication_plot(p_scrna, "Figure_04_scRNA_cell_fate_heatmap", 12, 5.8)

  fate_counts <- as.data.frame(table(fate = dominant_fate), stringsAsFactors = FALSE)
  fate_counts$fraction <- fate_counts$Freq / sum(fate_counts$Freq)
  p_fate <- ggplot(fate_counts, aes(x = reorder(fate, fraction), y = fraction, fill = fate)) +
    geom_col(show.legend = FALSE) + coord_flip() +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
    labs(title = "Dominant simulated cell fate", subtitle = scrna_source,
         x = NULL, y = "Fraction of cells") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))
  save_publication_plot(p_fate, "Figure_05_scRNA_fate_composition", 8, 5)
} else {
  log_message("Generic scRNA-seq module skipped: --no-demo used without --scrna")
}

# ---------------------- GEO validation: liver cancer ------------------------

# These two studies answer different questions. GSE125449 contains malignant
# and stromal/immune cells from primary liver cancers and is used here only for
# tumour-cell-intrinsic state concordance. GSE140228 contains sorted CD45+
# immune cells and is therefore restricted to microenvironmental support.
# Neither dataset contains a controlled perturbation arm; the output must not
# be described as causal confirmation of a knockout or overexpression.

geo_registry <- list(
  GSE125449 = list(
    title = "Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer",
    citation = "Ma et al., Cancer Cell 2019; PMID:31588021; DOI:10.1016/j.ccell.2019.08.007",
    scope = "tumour-cell-intrinsic observational concordance",
    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/",
    files = c(
      "GSE125449_Set1_barcodes.tsv.gz", "GSE125449_Set1_genes.tsv.gz",
      "GSE125449_Set1_matrix.mtx.gz", "GSE125449_Set1_samples.txt.gz",
      "GSE125449_Set2_barcodes.tsv.gz", "GSE125449_Set2_genes.tsv.gz",
      "GSE125449_Set2_matrix.mtx.gz", "GSE125449_Set2_samples.txt.gz"
    )
  ),
  GSE140228 = list(
    title = "Landscape and Dynamics of Single Immune Cells in Hepatocellular Carcinoma",
    citation = "Zhang et al., Cell 2019; PMID:31675496; DOI:10.1016/j.cell.2019.10.003",
    scope = "CD45+ immune-microenvironment observational support only",
    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140228/suppl/",
    droplet_files = c(
      "GSE140228_UMI_counts_Droplet.mtx.gz",
      "GSE140228_UMI_counts_Droplet_barcodes.tsv.gz",
      "GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz",
      "GSE140228_UMI_counts_Droplet_genes.tsv.gz"
    ),
    smartseq2_files = c(
      "GSE140228_cell_info_Smartseq2.tsv.gz",
      "GSE140228_gene_info_Smartseq2.tsv.gz",
      "GSE140228_read_counts_Smartseq2.csv.gz"
    )
  )
)

geo_requested_files <- function(accession, platform = "droplet") {
  entry <- geo_registry[[accession]]
  if (is.null(entry)) stop("Unsupported GEO accession: ", accession)
  if (accession == "GSE125449") return(entry$files)
  if (platform == "droplet") return(entry$droplet_files)
  if (platform == "smartseq2") return(entry$smartseq2_files)
  c(entry$droplet_files, entry$smartseq2_files)
}

download_geo_bundle <- function(accession, cache_dir, platform = "droplet") {
  entry <- geo_registry[[accession]]
  files <- geo_requested_files(accession, platform)
  target_dir <- file.path(cache_dir, accession)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    filename <- files[[i]]
    destination <- file.path(target_dir, filename)
    was_cached <- file.exists(destination) && file.info(destination)$size > 0
    if (!was_cached) {
      log_message("Downloading ", accession, ": ", filename)
      ok <- tryCatch({
        utils::download.file(paste0(entry$base_url, filename), destination,
                             mode = "wb", method = "libcurl", quiet = FALSE)
        TRUE
      }, error = function(e) {
        if (file.exists(destination)) unlink(destination)
        stop("GEO download failed for ", filename, ": ", conditionMessage(e))
      })
      if (!ok || !file.exists(destination) || file.info(destination)$size <= 0) {
        stop("Downloaded GEO file is missing or empty: ", destination)
      }
    }
    rows[[i]] <- data.frame(
      accession = accession, file = filename,
      local_path = normalizePath(destination, mustWork = FALSE),
      size_bytes = file.info(destination)$size,
      source = if (was_cached) "cache" else "downloaded",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

read_gz_table <- function(path, header = TRUE, sep = "\t") {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.table(con, header = header, sep = sep, quote = "", comment.char = "",
                    check.names = FALSE, fill = TRUE, stringsAsFactors = FALSE)
}

read_gz_vector <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  trimws(scan(con, what = character(), sep = "\n", quiet = TRUE))
}

read_sparse_mtx_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  methods::as(Matrix::readMM(con), "dgCMatrix")
}

extract_network_gene_counts <- function(mat, symbols, gene_map) {
  if (nrow(mat) != length(symbols)) stop("Gene annotation and matrix row counts differ")
  symbols <- trimws(as.character(symbols))
  upper_symbols <- toupper(symbols)
  wanted <- unique(toupper(unlist(gene_map, use.names = FALSE)))
  aliases <- intersect(wanted, unique(upper_symbols))
  if (!length(aliases)) stop("None of the network genes was found in the expression matrix")
  out <- matrix(0, nrow = length(aliases), ncol = ncol(mat),
                dimnames = list(aliases, colnames(mat)))
  for (i in seq_along(aliases)) {
    idx <- which(upper_symbols == aliases[[i]])
    out[i, ] <- as.numeric(Matrix::colSums(mat[idx, , drop = FALSE]))
  }
  attr(out, "library_size") <- as.numeric(Matrix::colSums(mat))
  out
}

geo_file <- function(manifest, pattern) {
  hit <- manifest$local_path[grepl(pattern, manifest$file, ignore.case = TRUE)]
  if (length(hit) != 1L) stop("Expected exactly one GEO file matching: ", pattern)
  hit[[1]]
}

read_gse125449 <- function(manifest) {
  experiments <- list()
  for (set_name in c("Set1", "Set2")) {
    genes <- read_gz_table(geo_file(manifest, paste0(set_name, "_genes\\.tsv")),
                           header = FALSE)
    barcodes <- read_gz_vector(geo_file(manifest, paste0(set_name, "_barcodes\\.tsv")))
    metadata <- read_gz_table(geo_file(manifest, paste0(set_name, "_samples\\.txt")))
    mat <- read_sparse_mtx_gz(geo_file(manifest, paste0(set_name, "_matrix\\.mtx")))
    if (ncol(mat) != length(barcodes) || nrow(metadata) != ncol(mat)) {
      stop("GSE125449 ", set_name, ": matrix/barcode/metadata dimensions differ")
    }
    cell_ids <- make.unique(paste(metadata$Sample, barcodes, sep = "::"))
    colnames(mat) <- cell_ids
    metadata$cell_id <- cell_ids
    metadata$sample_id <- metadata$Sample
    metadata$cell_type <- metadata$Type
    metadata$tissue <- "Tumour"
    metadata$cohort <- set_name
    rownames(metadata) <- cell_ids
    symbols <- if (ncol(genes) >= 2L) genes[[2]] else genes[[1]]
    target <- extract_network_gene_counts(mat, symbols, gene_map)
    experiments[[paste0("GSE125449_", set_name)]] <- list(
      accession = "GSE125449", experiment = paste0("GSE125449_", set_name),
      platform = "10x", scope = geo_registry$GSE125449$scope,
      counts = target, metadata = metadata
    )
    rm(mat); invisible(gc())
  }
  experiments
}

align_geo_metadata <- function(metadata, barcodes, accession, platform) {
  barcode_column <- names(metadata)[tolower(names(metadata)) == "barcode"]
  if (!length(barcode_column)) stop(accession, " ", platform, ": Barcode column not found")
  barcodes <- trimws(barcodes)
  metadata[[barcode_column[[1]]]] <- trimws(metadata[[barcode_column[[1]]]])
  idx <- match(barcodes, metadata[[barcode_column[[1]]]])
  if (anyNA(idx)) stop(accession, " ", platform, ": some matrix barcodes lack metadata")
  metadata <- metadata[idx, , drop = FALSE]
  metadata$cell_id <- make.unique(barcodes)
  metadata$sample_id <- if ("Sample" %in% names(metadata)) metadata$Sample else metadata$cell_id
  metadata$cell_type <- if ("celltype_global" %in% names(metadata)) metadata$celltype_global else "CD45+"
  metadata$cohort <- platform
  rownames(metadata) <- metadata$cell_id
  metadata
}

read_gse140228_droplet <- function(manifest) {
  genes <- read_gz_table(geo_file(manifest, "Droplet_genes\\.tsv"))
  barcodes <- read_gz_vector(geo_file(manifest, "Droplet_barcodes\\.tsv"))
  metadata <- read_gz_table(geo_file(manifest, "Droplet_cellinfo\\.tsv"))
  mat <- read_sparse_mtx_gz(geo_file(manifest, "Droplet\\.mtx"))
  if (ncol(mat) != length(barcodes)) stop("GSE140228 Droplet: matrix/barcode dimensions differ")
  metadata <- align_geo_metadata(metadata, barcodes, "GSE140228", "Droplet")
  colnames(mat) <- metadata$cell_id
  symbol_column <- names(genes)[toupper(names(genes)) == "SYMBOL"]
  symbols <- if (length(symbol_column)) genes[[symbol_column[[1]]]] else genes[[2]]
  target <- extract_network_gene_counts(mat, symbols, gene_map)
  rm(mat); invisible(gc())
  list(GSE140228_Droplet = list(
    accession = "GSE140228", experiment = "GSE140228_Droplet", platform = "Droplet",
    scope = geo_registry$GSE140228$scope, counts = target, metadata = metadata
  ))
}

read_gse140228_smartseq2 <- function(manifest) {
  count_path <- geo_file(manifest, "read_counts_Smartseq2\\.csv")
  info <- read_gz_table(geo_file(manifest, "gene_info_Smartseq2\\.tsv"))
  metadata <- read_gz_table(geo_file(manifest, "cell_info_Smartseq2\\.tsv"))
  raw <- data.table::fread(count_path, data.table = FALSE, check.names = FALSE)
  gene_ids <- as.character(raw[[1]])
  mat <- as.matrix(raw[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  barcodes <- colnames(mat)
  metadata <- align_geo_metadata(metadata, barcodes, "GSE140228", "Smartseq2")
  colnames(mat) <- metadata$cell_id
  ens_column <- names(info)[toupper(names(info)) == "ENSEMBL"]
  symbol_column <- names(info)[toupper(names(info)) == "SYMBOL"]
  match_index <- if (length(ens_column)) match(sub("\\..*$", "", gene_ids), info[[ens_column[[1]]]]) else NA
  if (length(symbol_column) && sum(!is.na(match_index)) > nrow(mat) / 2) {
    symbols <- info[[symbol_column[[1]]]][match_index]
  } else if (nrow(info) == nrow(mat) && length(symbol_column)) {
    symbols <- info[[symbol_column[[1]]]]
  } else {
    symbols <- gene_ids
  }
  sparse_mat <- methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")
  target <- extract_network_gene_counts(sparse_mat, symbols, gene_map)
  rm(mat, sparse_mat, raw); invisible(gc())
  list(GSE140228_Smartseq2 = list(
    accession = "GSE140228", experiment = "GSE140228_Smartseq2", platform = "Smartseq2",
    scope = geo_registry$GSE140228$scope, counts = target, metadata = metadata
  ))
}

read_geo_experiments <- function(accession, manifest, platform) {
  if (accession == "GSE125449") return(read_gse125449(manifest))
  out <- list()
  if (platform %in% c("droplet", "all")) out <- c(out, read_gse140228_droplet(manifest))
  if (platform %in% c("smartseq2", "all")) out <- c(out, read_gse140228_smartseq2(manifest))
  out
}

stratified_cell_subset <- function(metadata, max_cells, seed, group_columns) {
  if (nrow(metadata) <= max_cells) return(seq_len(nrow(metadata)))
  group_columns <- intersect(group_columns, names(metadata))
  if (!length(group_columns)) {
    set.seed(seed)
    return(sort(sample(seq_len(nrow(metadata)), max_cells)))
  }
  group <- interaction(metadata[, group_columns, drop = FALSE], drop = TRUE, lex.order = TRUE)
  members <- split(seq_len(nrow(metadata)), group)
  set.seed(seed)
  quota <- max(1L, floor(max_cells / length(members)))
  chosen <- unlist(lapply(members, function(x) sample(x, min(length(x), quota))), use.names = FALSE)
  if (length(chosen) < max_cells) {
    remainder <- setdiff(seq_len(nrow(metadata)), chosen)
    chosen <- c(chosen, sample(remainder, min(length(remainder), max_cells - length(chosen))))
  }
  sort(chosen[seq_len(min(length(chosen), max_cells))])
}

prepare_geo_experiment <- function(experiment, max_cells, seed) {
  metadata <- experiment$metadata
  counts <- experiment$counts
  if (experiment$accession == "GSE125449") {
    malignant <- grepl("^Malignant", metadata$cell_type, ignore.case = TRUE)
    if (any(malignant)) {
      metadata <- metadata[malignant, , drop = FALSE]
      lib <- attr(counts, "library_size", exact = TRUE)[malignant]
      counts <- counts[, malignant, drop = FALSE]
      attr(counts, "library_size") <- lib
    }
    group_columns <- c("cohort")
  } else {
    group_columns <- c("Tissue", "celltype_global")
  }
  max_cells <- if (CFG$quick) min(max_cells, 500L) else max_cells
  keep <- stratified_cell_subset(metadata, max_cells, seed, group_columns)
  lib <- attr(counts, "library_size", exact = TRUE)[keep]
  counts <- counts[, keep, drop = FALSE]
  attr(counts, "library_size") <- lib
  metadata <- metadata[keep, , drop = FALSE]
  if (experiment$accession == "GSE125449") {
    metadata$validation_group <- metadata$cohort
  } else {
    metadata$validation_group <- metadata$Tissue
  }
  binary <- probabilistic_binarise_scrna(counts, gene_map)
  fates <- simulate_single_cells(binary$node_probability, MODEL,
                                 n_draws = if (CFG$quick) 1L else 3L)
  experiment$counts <- counts
  experiment$metadata <- metadata
  experiment$binary <- binary
  experiment$fates <- fates
  experiment
}

summarise_geo_coverage <- function(experiment) {
  counts <- experiment$counts
  node_probability <- experiment$binary$node_probability
  do.call(rbind, lapply(names(gene_map), function(node) {
    aliases <- intersect(toupper(gene_map[[node]]), toupper(rownames(counts)))
    detected <- if (length(aliases)) {
      colSums(counts[aliases, , drop = FALSE]) > 0
    } else rep(FALSE, ncol(counts))
    data.frame(
      accession = experiment$accession, experiment = experiment$experiment,
      platform = experiment$platform, validation_scope = experiment$scope,
      node = node, measured = length(aliases) > 0,
      matched_aliases = paste(aliases, collapse = ";"),
      detection_rate = mean(detected),
      mean_activation_probability = if (all(is.na(node_probability[node, ]))) NA_real_ else
        mean(node_probability[node, ], na.rm = TRUE),
      n_cells = ncol(counts), stringsAsFactors = FALSE
    )
  }))
}

summarise_geo_fates <- function(experiment) {
  groups <- unique(experiment$metadata$validation_group)
  do.call(rbind, lapply(groups, function(group) {
    idx <- which(experiment$metadata$validation_group == group)
    data.frame(
      accession = experiment$accession, experiment = experiment$experiment,
      validation_group = group, validation_scope = experiment$scope,
      fate = colnames(experiment$fates),
      mean_probability = colMeans(experiment$fates[idx, , drop = FALSE]),
      n_cells = length(idx), stringsAsFactors = FALSE
    )
  }))
}

summarise_geo_perturbation_concordance <- function(experiment, perturbation_results) {
  groups <- unique(experiment$metadata$validation_group)
  expected <- perturbation_results[perturbation_results$node %in% names(gene_map), ]
  do.call(rbind, lapply(groups, function(group) {
    idx <- which(experiment$metadata$validation_group == group)
    observed <- rowMeans(experiment$binary$node_probability[, idx, drop = FALSE], na.rm = TRUE)
    observed[!is.finite(observed)] <- NA_real_
    do.call(rbind, lapply(unique(expected$perturbation), function(perturbation) {
      d <- expected[expected$perturbation == perturbation, ]
      exp_state <- stats::setNames(d$activation_frequency, d$node)
      nodes <- intersect(names(exp_state), names(observed)[is.finite(observed)])
      concordance <- if (length(nodes)) 1 - mean(abs(observed[nodes] - exp_state[nodes])) else NA_real_
      rank_correlation <- if (length(nodes) >= 3L && stats::sd(observed[nodes]) > 0 &&
                              stats::sd(exp_state[nodes]) > 0) {
        suppressWarnings(stats::cor(observed[nodes], exp_state[nodes], method = "spearman"))
      } else NA_real_
      data.frame(
        accession = experiment$accession, experiment = experiment$experiment,
        validation_group = group, validation_scope = experiment$scope,
        perturbation = perturbation, concordance_score = concordance,
        spearman_rank_correlation = rank_correlation,
        nodes_compared = length(nodes), node_coverage = length(nodes) / length(gene_map),
        compared_nodes = paste(nodes, collapse = ";"), n_cells = length(idx),
        claim_level = if (experiment$accession == "GSE125449")
          "observational tumour-state concordance; not causal perturbation validation" else
          "immune-context support only; not tumour-cell perturbation validation",
        stringsAsFactors = FALSE
      )
    }))
  }))
}

summarise_geo_edge_support <- function(experiment, model) {
  node_probability <- experiment$binary$node_probability
  mapped_edges <- model$edges[model$edges$from %in% rownames(node_probability) &
                                model$edges$to %in% rownames(node_probability), ]
  rows <- lapply(seq_len(nrow(mapped_edges)), function(i) {
    edge <- mapped_edges[i, ]
    x <- node_probability[edge$from, ]
    y <- node_probability[edge$to, ]
    ok <- is.finite(x) & is.finite(y)
    rho <- p_value <- NA_real_
    if (sum(ok) >= 10L && stats::sd(x[ok]) > 0 && stats::sd(y[ok]) > 0) {
      test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      rho <- unname(test$estimate)
      p_value <- test$p.value
    }
    expected_direction <- if (edge$sign == "negative") -1 else 1
    data.frame(
      accession = experiment$accession, experiment = experiment$experiment,
      validation_scope = experiment$scope,
      edge = paste0(edge$from, ifelse(expected_direction < 0, " -| ", " -> "), edge$to),
      from = edge$from, to = edge$to, expected_sign = expected_direction,
      spearman_rho = rho, p_value = p_value,
      sign_concordant = if (is.na(rho) || rho == 0) NA else sign(rho) == expected_direction,
      n_cells = sum(ok), stringsAsFactors = FALSE
    )
  })
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) out$fdr_bh <- stats::p.adjust(out$p_value, method = "BH")
  out
}

if (length(CFG$geo)) {
  unsupported <- setdiff(CFG$geo, names(geo_registry))
  if (length(unsupported)) stop("Unsupported GEO accession(s): ", paste(unsupported, collapse = ", "))
  dir.create(CFG$geo_dir, recursive = TRUE, showWarnings = FALSE)
  geo_manifests <- lapply(CFG$geo, download_geo_bundle,
                          cache_dir = CFG$geo_dir, platform = CFG$geo_platform)
  geo_download_manifest <- do.call(rbind, geo_manifests)
  write_table(geo_download_manifest, "13_geo_download_manifest.csv")

  geo_experiments <- list()
  for (accession in CFG$geo) {
    geo_experiments <- c(
      geo_experiments,
      read_geo_experiments(accession, geo_manifests[[match(accession, CFG$geo)]], CFG$geo_platform)
    )
  }
  for (i in seq_along(geo_experiments)) {
    geo_experiments[[i]] <- prepare_geo_experiment(
      geo_experiments[[i]], CFG$geo_max_cells, CFG$seed + i
    )
  }

  geo_scope <- do.call(rbind, lapply(geo_experiments, function(x) data.frame(
    accession = x$accession, experiment = x$experiment, platform = x$platform,
    title = geo_registry[[x$accession]]$title,
    citation = geo_registry[[x$accession]]$citation,
    validation_scope = x$scope, analysed_cells = ncol(x$counts),
    causal_perturbation_labels = FALSE,
    permitted_claim = if (x$accession == "GSE125449")
      "tumour-cell state concordance with simulated endpoints" else
      "CD45+ immune-context support only",
    stringsAsFactors = FALSE
  )))
  geo_coverage <- do.call(rbind, lapply(geo_experiments, summarise_geo_coverage))
  geo_fates <- do.call(rbind, lapply(geo_experiments, summarise_geo_fates))
  geo_concordance <- do.call(rbind, lapply(
    geo_experiments, summarise_geo_perturbation_concordance,
    perturbation_results = perturbation_results
  ))
  edge_tables <- lapply(geo_experiments, summarise_geo_edge_support, model = MODEL)
  edge_tables <- edge_tables[vapply(edge_tables, nrow, integer(1)) > 0]
  geo_edges <- if (length(edge_tables)) do.call(rbind, edge_tables) else data.frame()

  write_table(geo_scope, "14_geo_dataset_scope_and_claims.csv")
  write_table(geo_coverage, "15_geo_node_coverage.csv")
  write_table(geo_concordance, "16_geo_perturbation_concordance.csv")
  write_table(geo_fates, "17_geo_group_fate_probabilities.csv")
  if (nrow(geo_edges)) write_table(geo_edges, "18_geo_regulatory_edge_support.csv")

  consensus <- stats::aggregate(
    cbind(concordance_score, node_coverage) ~ accession + perturbation,
    data = geo_concordance, FUN = mean, na.rm = TRUE
  )
  names(consensus)[names(consensus) == "concordance_score"] <- "mean_concordance"
  names(consensus)[names(consensus) == "node_coverage"] <- "mean_node_coverage"
  consensus$validation_role <- ifelse(
    consensus$accession == "GSE125449", "tumour-state concordance",
    "immune-context support only"
  )
  consensus$causal_confirmation <- FALSE
  write_table(consensus, "19_geo_cross_dataset_summary.csv")

  coverage_plot_data <- geo_coverage[geo_coverage$measured, ]
  p_geo_coverage <- ggplot(coverage_plot_data,
                           aes(x = node, y = experiment, fill = detection_rate)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * detection_rate)), size = 2.4) +
    scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1),
                        labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "GEO single-cell coverage of GSDME-network nodes",
         subtitle = "Detection rates are calculated separately for each experiment/platform",
         x = NULL, y = NULL, fill = "Detected") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
  save_publication_plot(p_geo_coverage, "Figure_08_GEO_node_detection_heatmap", 12, 5.5)

  concordance_plot_data <- stats::aggregate(
    concordance_score ~ experiment + perturbation, data = geo_concordance,
    FUN = mean, na.rm = TRUE
  )
  p_geo_concordance <- ggplot(concordance_plot_data,
                              aes(x = perturbation, y = experiment, fill = concordance_score)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", concordance_score)), size = 2.4) +
    scale_fill_gradientn(colours = c("#B2182B", "#F7F7F7", "#2166AC"), limits = c(0, 1)) +
    labs(title = "Observed single-cell concordance with simulated perturbation endpoints",
         subtitle = "Observational state matching only; the GEO studies contain no perturbation arms",
         x = NULL, y = NULL, fill = "Concordance") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
  save_publication_plot(p_geo_concordance, "Figure_09_GEO_perturbation_concordance", 13, 5.8)

  p_geo_fates <- ggplot(geo_fates,
                        aes(x = validation_group, y = mean_probability, fill = fate)) +
    geom_col(width = 0.8) + facet_wrap(~ experiment, scales = "free_x") +
    scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "Boolean fate probabilities informed by GEO single-cell profiles",
         subtitle = "GSE125449: malignant cells; GSE140228: CD45+ immune cells interpreted as context only",
         x = NULL, y = "Mean probability", fill = "Fate") +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"), legend.position = "bottom")
  save_publication_plot(p_geo_fates, "Figure_10_GEO_predicted_fate_composition", 14, 7.5)

  if (nrow(geo_edges)) {
    p_geo_edges <- ggplot(geo_edges, aes(x = edge, y = experiment, fill = spearman_rho)) +
      geom_tile(colour = "white", linewidth = 0.35) +
      geom_point(aes(shape = sign_concordant), size = 2.2, na.rm = TRUE) +
      scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                           midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
      labs(title = "Exploratory transcriptional support for mapped regulatory edges",
           subtitle = "Spearman correlation; dots indicate agreement with the model edge sign",
           x = NULL, y = NULL, fill = "Spearman rho", shape = "Sign concordant") +
      theme_minimal(base_size = 9) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"), legend.position = "bottom")
    save_publication_plot(p_geo_edges, "Figure_11_GEO_regulatory_edge_support", 13, 6.5)
  }
}

# ------------------ Multi-omic patient soft-clipping -------------------------

load_multiomics <- function(path) {
  if (!file.exists(path)) stop("Multi-omics file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  x <- utils::read.table(path, header = TRUE, sep = sep, check.names = FALSE,
                         comment.char = "", quote = "\"")
  required <- c("patient_id", "node", "expression_z", "methylation_beta",
                "cnv_log2", "mutation_effect")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing multi-omics column(s): ", paste(missing, collapse = ", "))
  if (anyDuplicated(paste(x$patient_id, x$node, sep = "::"))) {
    stop("The multi-omics table must contain one row per patient_id/node pair")
  }
  x
}

make_demo_multiomics <- function(model, n_patients = 24L, seed = 204L) {
  set.seed(seed + 1L)
  nodes <- intersect(c("GSDME_availability", "lncRNA_MALAT1", "miR_204_5p",
                       "SIRT1", "PGC1A", "BCL2", "BAX", "p53_ACTIVE",
                       "CASP3", "RB1", "E2F1"), model$node_ids)
  grid <- expand.grid(patient_id = paste0("Patient_", seq_len(n_patients)),
                      node = nodes, stringsAsFactors = FALSE)
  n <- nrow(grid)
  grid$expression_z <- stats::rnorm(n)
  grid$methylation_beta <- stats::rbeta(n, 2, 3)
  grid$cnv_log2 <- stats::rnorm(n, 0, 0.45)
  grid$mutation_effect <- sample(c(-1L, 0L, 1L), n, replace = TRUE,
                                 prob = c(0.03, 0.94, 0.03))
  # Structured demo contrast: first third is more resistant-like.
  resistant <- grid$patient_id %in% paste0("Patient_", seq_len(ceiling(n_patients / 3)))
  grid$expression_z[resistant & grid$node %in% c("lncRNA_MALAT1", "SIRT1", "BCL2")] <-
    grid$expression_z[resistant & grid$node %in% c("lncRNA_MALAT1", "SIRT1", "BCL2")] + 1.5
  grid$expression_z[resistant & grid$node %in% c("GSDME_availability", "miR_204_5p")] <-
    grid$expression_z[resistant & grid$node %in% c("GSDME_availability", "miR_204_5p")] - 1.2
  grid
}

omics_activation_prior <- function(x, weights = c(expression = 0.9, methylation = 2.0,
                                                   cnv = 1.1, mutation = 4.0)) {
  eta <- weights[["expression"]] * x$expression_z -
    weights[["methylation"]] * (x$methylation_beta - 0.5) +
    weights[["cnv"]] * x$cnv_log2 +
    weights[["mutation"]] * x$mutation_effect
  stats::plogis(eta)
}

step_soft_clipped <- function(state, model, priors, base_clamp = integer(0),
                              soft_strength = 0.35) {
  state <- apply_clamp(normalise_state(state, model), base_clamp)
  logical <- logical_targets(state, model, base_clamp)
  p <- as.numeric(logical)
  names(p) <- names(logical)
  common <- intersect(names(priors), names(p))
  p[common] <- (1 - soft_strength) * p[common] + soft_strength * priors[common]
  next_state <- stats::setNames(stats::rbinom(length(p), 1L, pmin(pmax(p, 0), 1)), names(p))
  hard_low <- common[priors[common] <= 0.02]
  hard_high <- common[priors[common] >= 0.98]
  next_state[hard_low] <- 0L
  next_state[hard_high] <- 1L
  apply_clamp(next_state, base_clamp)
}

simulate_digital_twins <- function(multiomics, model, n_rep = 80L, n_steps = 35L,
                                   soft_strength = 0.35) {
  multiomics$activation_prior <- omics_activation_prior(multiomics)
  patients <- unique(multiomics$patient_id)
  fate_levels <- c("Pyroptosis", "Apoptosis", "Resistance", "Proliferation",
                   "Cell-cycle arrest", "Survival", "Other")
  out <- matrix(0, nrow = length(patients), ncol = length(fate_levels),
                dimnames = list(patients, fate_levels))
  for (patient in patients) {
    d <- multiomics[multiomics$patient_id == patient, ]
    priors <- stats::setNames(d$activation_prior, d$node)
    priors <- priors[names(priors) %in% model$node_ids]
    fates <- character(n_rep)
    for (r in seq_len(n_rep)) {
      replicate_clamp <- BASE_CLAMP
      if ("GSDME_availability" %in% names(priors)) {
        replicate_clamp <- c(replicate_clamp,
                             GSDME_availability = stats::rbinom(1L, 1L, priors[["GSDME_availability"]]))
      }
      s <- random_initial_state(model, replicate_clamp)
      common <- intersect(names(priors), names(s))
      s[common] <- stats::rbinom(length(common), 1L, priors[common])
      s <- apply_clamp(s, replicate_clamp)
      for (step in seq_len(n_steps)) {
        s <- step_soft_clipped(s, model, priors, replicate_clamp, soft_strength)
      }
      fates[[r]] <- classify_fate(s)
    }
    out[patient, ] <- tabulate(match(fates, fate_levels), nbins = length(fate_levels)) / n_rep
  }
  list(phenotypes = out, annotated_multiomics = multiomics)
}

if (!is.na(CFG$multiomics)) {
  multiomics <- load_multiomics(CFG$multiomics)
  multiomics_source <- "User-provided"
} else if (CFG$demo) {
  multiomics <- make_demo_multiomics(MODEL, if (CFG$quick) 8L else 24L, CFG$seed)
  multiomics_source <- "SIMULATED DEMO DATA"
} else {
  multiomics <- NULL
  log_message("Multi-omics module skipped: --no-demo used without --multiomics")
}

if (!is.null(multiomics)) {
  digital_twins <- simulate_digital_twins(
    multiomics, MODEL,
    n_rep = if (CFG$quick) 20L else 80L,
    n_steps = if (CFG$quick) 15L else 35L
  )
  write_table(digital_twins$annotated_multiomics, "09_patient_omics_activation_priors.csv")
  digital_twin_table <- data.frame(patient_id = rownames(digital_twins$phenotypes),
                                   digital_twins$phenotypes, check.names = FALSE)
  digital_twin_table$data_source <- multiomics_source
  write_table(digital_twin_table, "10_digital_twin_phenotype_probabilities.csv")

  twin_long <- do.call(rbind, lapply(seq_len(nrow(digital_twins$phenotypes)), function(i) {
    data.frame(patient_id = rownames(digital_twins$phenotypes)[[i]],
               phenotype = colnames(digital_twins$phenotypes),
               probability = as.numeric(digital_twins$phenotypes[i, ]),
               stringsAsFactors = FALSE)
  }))
  p_twins <- plot_heatmap_long(
    twin_long, "patient_id", "phenotype", "probability",
    "Patient-specific multi-omic network instances",
    paste0(multiomics_source, "; probabilistic soft-clipping of Boolean update rules"),
    low = "#FFF7EC", high = "#7F0000"
  )
  save_publication_plot(p_twins, "Figure_06_multiomic_digital_twins", 10.5, 8)
}

# ---------------- Q-learning for sequential therapy --------------------------

rl_actions <- list(
  NONE = integer(0),
  MALAT1_KO = c(lncRNA_MALAT1 = 0L),
  miR204_OE = c(miR_204_5p = 1L),
  SIRT1_KO = c(SIRT1 = 0L),
  PGC1A_KO = c(PGC1A = 0L),
  BCL2_KO = c(BCL2 = 0L)
)
rl_actions <- rl_actions[vapply(rl_actions, function(x) all(names(x) %in% MODEL$node_ids), logical(1))]

rl_reward <- function(state, action_name) {
  value <- 0
  if (state[["PYROPTOSIS_GSDME"]] == 1L) value <- value + 10
  if (state[["APOPTOSIS"]] == 1L) value <- value + 7
  if (state[["RESISTANCE"]] == 1L) value <- value - 6
  if (state[["PROLIFERATION"]] == 1L) value <- value - 4
  if (state[["SURVIVAL"]] == 1L) value <- value - 2
  value <- value - ifelse(action_name == "NONE", 0.02, 0.25)
  value
}

q_key <- function(state, model, step) paste0(step, ":", encode_state(state, model))

train_q_learning <- function(model, initial_state, actions, episodes = 5000L,
                             horizon = 12L, alpha = 0.15, gamma = 0.92,
                             epsilon_start = 0.9, epsilon_end = 0.05) {
  q <- new.env(hash = TRUE, parent = emptyenv())
  action_names <- names(actions)
  get_q <- function(key) {
    if (!exists(key, q, inherits = FALSE)) assign(key, rep(0, length(action_names)), q)
    get(key, q, inherits = FALSE)
  }
  for (episode in seq_len(episodes)) {
    state <- jitter_state(initial_state, probability = 0.04, protected = names(BASE_CLAMP))
    epsilon <- epsilon_end + (epsilon_start - epsilon_end) * (1 - episode / episodes)
    for (step in seq_len(horizon)) {
      key <- q_key(state, model, step)
      q_values <- get_q(key)
      if (stats::runif(1) < epsilon) action_index <- sample(seq_along(actions), 1L)
      else action_index <- sample(which(q_values == max(q_values)), 1L)
      action_name <- action_names[[action_index]]
      transient_clamp <- c(BASE_CLAMP, GSDME_availability = 1L, actions[[action_name]])
      transient_clamp <- transient_clamp[!duplicated(names(transient_clamp), fromLast = TRUE)]
      next_state <- step_boolean(state, model, transient_clamp, "synchronous")
      reward <- rl_reward(next_state, action_name)
      next_key <- q_key(next_state, model, min(step + 1L, horizon))
      next_q <- if (step == horizon) rep(0, length(actions)) else get_q(next_key)
      q_values[[action_index]] <- q_values[[action_index]] +
        alpha * (reward + gamma * max(next_q) - q_values[[action_index]])
      assign(key, q_values, q)
      state <- next_state
      if (state[["PYROPTOSIS_GSDME"]] == 1L || state[["APOPTOSIS"]] == 1L) break
    }
  }
  list(q = q, action_names = action_names, actions = actions, horizon = horizon,
       get_q = get_q)
}

rollout_policy <- function(agent, model, initial_state) {
  state <- initial_state
  rows <- list()
  for (step in seq_len(agent$horizon)) {
    key <- q_key(state, model, step)
    q_values <- agent$get_q(key)
    action_index <- which.max(q_values)
    action_name <- agent$action_names[[action_index]]
    clamp <- c(BASE_CLAMP, GSDME_availability = 1L, agent$actions[[action_name]])
    clamp <- clamp[!duplicated(names(clamp), fromLast = TRUE)]
    next_state <- step_boolean(state, model, clamp, "synchronous")
    rows[[step]] <- data.frame(
      step = step, action = action_name, reward = rl_reward(next_state, action_name),
      fate = classify_fate(next_state), t(next_state), check.names = FALSE
    )
    state <- next_state
    if (state[["PYROPTOSIS_GSDME"]] == 1L || state[["APOPTOSIS"]] == 1L) break
  }
  do.call(rbind, rows)
}

agent <- train_q_learning(
  MODEL, MALIGNANT_STATE, rl_actions,
  episodes = if (CFG$quick) 600L else 5000L,
  horizon = if (CFG$quick) 8L else 12L
)
rl_sequence <- rollout_policy(agent, MODEL, MALIGNANT_STATE)
write_table(rl_sequence, "11_reinforcement_learning_sequence.csv")

rl_nodes <- intersect(c("lncRNA_MALAT1", "miR_204_5p", "SIRT1", "PGC1A",
                        "p53_ACTIVE", "PUMA", "BAX", "CASP3", "GSDME_N",
                        phenotype_nodes), names(rl_sequence))
rl_long <- do.call(rbind, lapply(seq_len(nrow(rl_sequence)), function(i) {
  data.frame(step = rl_sequence$step[[i]], action = rl_sequence$action[[i]],
             node = rl_nodes, state = as.numeric(rl_sequence[i, rl_nodes]),
             stringsAsFactors = FALSE)
}))
p_rl <- ggplot(rl_long, aes(x = step, y = node, fill = state)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), breaks = c(0, 1)) +
  scale_x_continuous(
    breaks = unique(rl_long$step),
    labels = stats::setNames(paste0(rl_sequence$step, "\n", rl_sequence$action), rl_sequence$step)
  ) +
  labs(title = "Q-learning-derived sequential intervention",
       subtitle = "Exploratory policy learned on the logical model; experimental validation required",
       x = "Simulation step", y = NULL, fill = "State") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(face = "bold"))
save_publication_plot(p_rl, "Figure_07_RL_sequential_therapy", 11, 7)

# ----------------------------- Manifest --------------------------------------

capture.output(sessionInfo(), file = file.path(LOG_DIR, "sessionInfo.txt"))
manifest <- data.frame(
  file = list.files(CFG$outdir, recursive = TRUE, full.names = FALSE),
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
  pipeline_version = PIPELINE_VERSION,
  seed = CFG$seed,
  stringsAsFactors = FALSE
)
write_table(manifest, "12_output_manifest.csv")

log_message("Pipeline completed. Results: ", normalizePath(CFG$outdir, mustWork = FALSE))
cat("\nCompleted successfully.\nResults: ", normalizePath(CFG$outdir, mustWork = FALSE), "\n", sep = "")
