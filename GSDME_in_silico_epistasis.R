#!/usr/bin/env Rscript

# =============================================================================
# In silico epistasis and rescue analysis for the GSDME logical network
# Version 1.3.0
# =============================================================================
#
# Purpose
#   Test model-predicted ordering of the MALAT1/miR-204-5p/SIRT1/p53/CASP3/
#   GSDME axis with matched single and double Boolean perturbations.
#
# Scope
#   Results are causal predictions within the encoded logical model. They are
#   not biological proof and must be described as "in silico epistasis" or
#   "model-based rescue". GSE125449 is used only to provide observational,
#   tumour-cell-derived initial states; it contains no controlled perturbation.
#
# Principal outputs
#   Figure_08_in_silico_epistasis_and_rescue.pdf + 600-dpi PNG
#   Supplementary_Figure_S2_epistasis_seed_robustness.pdf + 600-dpi PNG
#   CSV tables containing seed-level and consensus results
#
# Example full run
#   Rscript GSDME_in_silico_epistasis.R \
#     --model GINsim-miR_204_GSDME_Pyroptosis.zginml \
#     --seeds 101,204,307,509,811 \
#     --trajectories 500 \
#     --geo-dir GEO_scRNA_data \
#     --geo-max-cells 2000 \
#     --out resultados_epistasia_GSDME
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

SCRIPT_VERSION <- "1.3.0"

# ------------------------------- CLI -----------------------------------------

parse_cli <- function(args) {
  cfg <- list(
    model = "GINsim-miR_204_GSDME_Pyroptosis.zginml",
    outdir = "resultados_epistasia_GSDME",
    seeds = c(101L, 204L, 307L, 509L, 811L),
    trajectories = 500L,
    max_steps = 350L,
    jitter_probability = 0.08,
    geo = TRUE,
    geo_dir = "GEO_scRNA_data",
    geo_max_cells = 2000L,
    geo_cell_draws = 1L,
    geo_subset_seed = 204L,
    auto_install = TRUE,
    quick = FALSE
  )

  if (!length(args)) return(cfg)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      cat(paste0(
        "GSDME in silico epistasis and rescue analysis\n\n",
        "Usage:\n",
        "  Rscript GSDME_in_silico_epistasis.R --model MODEL --out OUT [options]\n\n",
        "Options:\n",
        "  --model PATH            GINsim .zginml/.ginml model\n",
        "  --out PATH              Output directory\n",
        "  --seeds LIST            Comma-separated seeds (default 101,204,307,509,811)\n",
        "  --trajectories N        Trajectories per condition and initial-state stratum\n",
        "  --max-steps N           Maximum asynchronous updates per trajectory\n",
        "  --jitter-probability X  Flip probability around resistant attractor (default 0.08)\n",
        "  --geo-dir PATH          Existing/download cache for GSE125449\n",
        "  --geo-max-cells N       Maximum malignant cells per GSE125449 set\n",
        "  --geo-cell-draws N      State draws per cell and seed (default 1)\n",
        "  --geo-subset-seed N     Fixed seed for identical cell subset across runs\n",
        "  --no-geo                Skip single-cell-informed epistasis\n",
        "  --no-install            Do not install missing CRAN packages\n",
        "  --quick                 Small diagnostic run; never use its numbers in a paper\n"
      ))
      quit(save = "no", status = 0)
    }

    if (key %in% c("--no-geo", "--no-install", "--quick")) {
      if (key == "--no-geo") cfg$geo <- FALSE
      if (key == "--no-install") cfg$auto_install <- FALSE
      if (key == "--quick") cfg$quick <- TRUE
      i <- i + 1L
      next
    }

    if (i == length(args)) stop("Missing value after ", key)
    value <- args[[i + 1L]]
    if (key == "--model") cfg$model <- value
    else if (key == "--out") cfg$outdir <- value
    else if (key == "--seeds") {
      cfg$seeds <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1]])))
    }
    else if (key == "--trajectories") cfg$trajectories <- as.integer(value)
    else if (key == "--max-steps") cfg$max_steps <- as.integer(value)
    else if (key == "--jitter-probability") cfg$jitter_probability <- as.numeric(value)
    else if (key == "--geo-dir") cfg$geo_dir <- value
    else if (key == "--geo-max-cells") cfg$geo_max_cells <- as.integer(value)
    else if (key == "--geo-cell-draws") cfg$geo_cell_draws <- as.integer(value)
    else if (key == "--geo-subset-seed") cfg$geo_subset_seed <- as.integer(value)
    else stop("Unknown argument: ", key)
    i <- i + 2L
  }
  cfg
}

CFG <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!length(CFG$seeds) || anyNA(CFG$seeds) || any(CFG$seeds < 1L)) {
  stop("--seeds must contain positive integers")
}
CFG$seeds <- unique(CFG$seeds)
if (is.na(CFG$trajectories) || CFG$trajectories < 20L) {
  stop("--trajectories must be an integer >= 20")
}
if (is.na(CFG$max_steps) || CFG$max_steps < 20L) {
  stop("--max-steps must be an integer >= 20")
}
if (!is.finite(CFG$jitter_probability) || CFG$jitter_probability < 0 ||
    CFG$jitter_probability > 0.5) {
  stop("--jitter-probability must be between 0 and 0.5")
}
if (is.na(CFG$geo_max_cells) || CFG$geo_max_cells < 20L) {
  stop("--geo-max-cells must be an integer >= 20")
}
if (is.na(CFG$geo_cell_draws) || CFG$geo_cell_draws < 1L) {
  stop("--geo-cell-draws must be an integer >= 1")
}
if (CFG$quick) {
  CFG$trajectories <- min(CFG$trajectories, 60L)
  CFG$geo_max_cells <- min(CFG$geo_max_cells, 200L)
  CFG$geo_cell_draws <- 1L
}

# ---------------------------- Dependencies -----------------------------------

required_packages <- c("xml2", "ggplot2", "patchwork")
if (CFG$geo) required_packages <- c(required_packages, "Matrix")
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

suppressPackageStartupMessages({
  library(xml2)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------ Output ---------------------------------------

dir.create(CFG$outdir, recursive = TRUE, showWarnings = FALSE)
FIG_DIR <- file.path(CFG$outdir, "figures")
TAB_DIR <- file.path(CFG$outdir, "tables")
LOG_DIR <- file.path(CFG$outdir, "logs")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "epistasis.log")
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

write_table <- function(x, filename) {
  utils::write.csv(x, file.path(TAB_DIR, filename), row.names = FALSE, na = "")
}

save_publication_plot <- function(plot, stem, width = 14, height = 10) {
  pdf_path <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIG_DIR, paste0(stem, "_600dpi.png"))
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(
    pdf_path, plot = plot, width = width, height = height,
    units = "in", device = pdf_device, bg = "white"
  )
  ggplot2::ggsave(
    png_path, plot = plot, width = width, height = height,
    units = "in", dpi = 600, bg = "white", limitsize = FALSE
  )
  invisible(c(pdf = pdf_path, png = png_path))
}

log_message(
  "Starting version ", SCRIPT_VERSION,
  if (CFG$quick) " [QUICK DIAGNOSTIC MODE]" else " [FULL MODE]"
)

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

  input_attribute <- xml2::xml_attr(node_xml, "input")
  nodes <- data.frame(
    id = xml2::xml_attr(node_xml, "id"),
    name = xml2::xml_attr(node_xml, "name"),
    maxvalue = as.integer(xml2::xml_attr(node_xml, "maxvalue")),
    input = !is.na(input_attribute) & tolower(input_attribute) == "true",
    rule = NA_character_, rule_name = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(node_xml)) {
    exp_node <- xml2::xml_find_first(node_xml[[i]], "./value/exp")
    if (!inherits(exp_node, "xml_missing")) nodes$rule[[i]] <- xml2::xml_attr(exp_node, "str")
  }
  nodes$name[is.na(nodes$name) | !nzchar(nodes$name)] <-
    nodes$id[is.na(nodes$name) | !nzchar(nodes$name)]
  nodes$label <- nodes$name
  nodes$maxvalue[is.na(nodes$maxvalue)] <- 1L

  # O arquivo GINsim escreve cada regra com ids. A partir daqui os estados,
  # clamps, arestas, transcritos e nomes das colunas usam somente name.
  id_to_name <- stats::setNames(nodes$name, nodes$id)
  rename_rule <- function(expr) {
    if (is.symbol(expr)) {
      id <- as.character(expr)
      if (!id %in% names(id_to_name)) stop("Unknown rule id: ", id)
      return(as.name(id_to_name[[id]]))
    }
    if (is.call(expr)) {
      operator <- as.character(expr[[1L]])
      # O parser de R representa (A | B) como uma chamada a `(`.
      # Os parenteses agrupam a regra; nao sao um operador Booleano.
      if (operator == "(" && length(expr) == 2L)
        return(rename_rule(expr[[2L]]))
      if (!operator %in% c("&", "|", "!"))
        stop("Unsupported Boolean operator in GINsim rule: ", deparse(expr))
      return(as.call(c(list(expr[[1L]]), lapply(as.list(expr)[-1L], rename_rule))))
    }
    expr
  }
  compiled <- lapply(nodes$rule, function(rule) {
    if (is.na(rule) || !nzchar(rule)) return(NULL)
    expr <- parse(text = rule)
    if (length(expr) != 1L) stop("Multiple expressions in a GINsim rule")
    rename_rule(expr[[1L]])
  })
  nodes$rule_name <- vapply(compiled, function(rule) if (is.null(rule)) ""
    else paste(deparse(rule), collapse = ""), "")

  edges <- data.frame(
    id = xml2::xml_attr(edge_xml, "id"),
    from_id = xml2::xml_attr(edge_xml, "from"),
    to_id = xml2::xml_attr(edge_xml, "to"),
    sign = xml2::xml_attr(edge_xml, "sign"),
    stringsAsFactors = FALSE
  )
  edges$from <- unname(id_to_name[edges$from_id])
  edges$to <- unname(id_to_name[edges$to_id])
  model <- list(
    id = xml2::xml_attr(graph_node, "id"),
    nodes = nodes,
    edges = edges,
    node_ids = nodes$name,
    input_ids = nodes$name[which(nodes$input)],
    id_to_name = id_to_name
  )
  validate_model(model)
  names(compiled) <- nodes$name
  model$compiled_rules <- compiled
  model
}

validate_model <- function(model) {
  if (anyNA(model$nodes$id) || any(!nzchar(model$nodes$id))) {
    stop("Every GINsim node must have a non-empty ID")
  }
  stopifnot(!anyDuplicated(model$nodes$id))
  if (anyNA(model$nodes$name) || any(!nzchar(model$nodes$name)) ||
      anyDuplicated(model$nodes$name)) stop("Every node needs a unique non-empty name")
  if (anyNA(model$nodes$maxvalue) || any(model$nodes$maxvalue != 1L))
    stop("Only Boolean GINsim nodes (maxvalue=1) are supported")
  unknown_edges <- setdiff(unique(c(model$edges$from, model$edges$to)), model$nodes$name)
  if (length(unknown_edges)) stop("Edges refer to unknown nodes: ", paste(unknown_edges, collapse = ", "))
  for (i in which(!is.na(model$nodes$rule))) {
    vars <- all.vars(parse(text = model$nodes$rule[[i]]))
    unknown <- setdiff(vars, model$nodes$id)
    if (length(unknown)) {
      stop("Rule for ", model$nodes$name[[i]], " refers to unknown id(s): ",
           paste(unknown, collapse = ", "))
    }
  }
  invisible(TRUE)
}

MODEL <- read_ginsim_model(CFG$model)
required_nodes <- c(
  "DDR", "DFNA5", "MALAT1", "MIR204",
  "SIRT1", "TP53", "CASP3", "PYROPTOSIS", "APOPTOSIS",
  "RESISTANCE", "PROLIFERATION", "CELL_CYCLE_ARREST"
)
missing_model_nodes <- setdiff(required_nodes, MODEL$node_ids)
if (length(missing_model_nodes)) {
  stop("Model lacks required node(s): ", paste(missing_model_nodes, collapse = ", "))
}
log_message(
  "Loaded model '", MODEL$id, "': ", length(MODEL$node_ids),
  " nodes and ", nrow(MODEL$edges), " edges"
)

# ------------------------ Boolean simulation core ----------------------------

empty_state <- function(model, value = 0L) {
  stats::setNames(rep(as.integer(value), length(model$node_ids)), model$node_ids)
}

normalise_state <- function(state, model) {
  output <- empty_state(model)
  common <- intersect(names(state), model$node_ids)
  output[common] <- as.integer(state[common] > 0)
  output
}

merge_clamps <- function(...) {
  clamp <- unlist(list(...), use.names = TRUE)
  if (!length(clamp)) return(integer(0))
  if (is.null(names(clamp)) || anyNA(names(clamp)) || any(!nzchar(names(clamp)))) {
    stop("All clamp values must have node names")
  }
  clamp <- clamp[!duplicated(names(clamp), fromLast = TRUE)]
  stats::setNames(as.integer(clamp), names(clamp))
}

apply_clamp <- function(state, clamp = integer(0)) {
  if (length(clamp)) {
    unknown <- setdiff(names(clamp), names(state))
    if (length(unknown)) {
      stop("Clamp refers to unknown node(s): ", paste(unknown, collapse = ", "))
    }
    state[names(clamp)] <- as.integer(clamp)
  }
  state
}

evaluate_rule <- function(rule_expression, state_environment) {
  if (is.null(rule_expression)) return(0L)
  as.integer(isTRUE(eval(rule_expression, envir = state_environment)))
}

logical_targets <- function(state, model, clamp = integer(0)) {
  state <- apply_clamp(normalise_state(state, model), clamp)
  target <- state
  state_values <- as.list(as.logical(unname(state)))
  names(state_values) <- names(state)
  state_environment <- list2env(state_values, parent = baseenv())
  for (i in seq_len(nrow(model$nodes))) {
    node <- model$nodes$name[[i]]
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
                                  max_steps = 350L) {
  state <- apply_clamp(normalise_state(initial_state, model), clamp)
  for (step in 0:max_steps) {
    target <- logical_targets(state, model, clamp)
    unstable <- setdiff(names(state)[state != target], names(clamp))
    if (!length(unstable)) {
      return(list(state = state, steps = step, stable = TRUE, cycle_length = 1L))
    }
    if (step == max_steps) break
    node <- sample(unstable, 1L)
    state[[node]] <- target[[node]]
  }
  list(state = state, steps = max_steps, stable = FALSE, cycle_length = NA_integer_)
}

simulate_synchronous_reference <- function(initial_state, model, clamp = integer(0),
                                           max_steps = 350L) {
  state <- apply_clamp(normalise_state(initial_state, model), clamp)
  seen <- character(0)
  for (step in seq_len(max_steps)) {
    key <- encode_state(state, model)
    if (key %in% seen) return(state)
    seen <- c(seen, key)
    next_state <- step_boolean(state, model, clamp, "synchronous")
    if (identical(unname(next_state), unname(state))) return(state)
    state <- next_state
  }
  state
}

phenotype_nodes <- intersect(
  c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE",
    "PROLIFERATION", "CELL_CYCLE_ARREST"),
  MODEL$node_ids
)

fate_levels <- c(
  "Pyroptosis", "Apoptosis", "Resistance", "Proliferation",
  "Cell-cycle arrest", "Other"
)

classify_fate <- function(state) {
  is_on <- function(node) node %in% names(state) && state[[node]] == 1L
  if (is_on("PYROPTOSIS")) return("Pyroptosis")
  if (is_on("APOPTOSIS")) return("Apoptosis")
  if (is_on("RESISTANCE")) return("Resistance")
  if (is_on("PROLIFERATION")) return("Proliferation")
  if (is_on("CELL_CYCLE_ARREST")) return("Cell-cycle arrest")
  "Other"
}

BASE_CONTEXT <- c(DDR = 1L)

make_resistant_reference <- function(model) {
  state <- empty_state(model)
  state[intersect(names(BASE_CONTEXT), names(state))] <- BASE_CONTEXT[
    intersect(names(BASE_CONTEXT), names(state))
  ]
  state[intersect(
    c("MALAT1", "SIRT1", "RESISTANCE",
      "(CCND1 ou CCND2 ou CCND3) E (CDK4 ou CDK6)", "E2F1"),
    names(state)
  )] <- 1L
  simulate_synchronous_reference(
    state, model, BASE_CONTEXT, max_steps = CFG$max_steps
  )
}

RESISTANT_REFERENCE <- make_resistant_reference(MODEL)

make_initial_states <- function(model, n, seed, stratum,
                                jitter_probability = 0.08) {
  set.seed(seed)
  if (stratum == "Global random") {
    matrix_values <- stats::rbinom(n * length(model$node_ids), 1L, 0.5)
    states <- matrix(
      matrix_values, nrow = n, ncol = length(model$node_ids),
      dimnames = list(NULL, model$node_ids)
    )
  } else if (stratum == "Resistant-local") {
    states <- matrix(
      rep(RESISTANT_REFERENCE, each = n), nrow = n,
      dimnames = list(NULL, model$node_ids)
    )
    protected <- names(BASE_CONTEXT)
    eligible <- setdiff(model$node_ids, protected)
    flips <- matrix(
      stats::runif(n * length(eligible)) < jitter_probability,
      nrow = n, ncol = length(eligible)
    )
    states[, eligible] <- abs(states[, eligible, drop = FALSE] - flips)
  } else {
    stop("Unknown initial-state stratum: ", stratum)
  }
  states[, "DDR"] <- 1L
  storage.mode(states) <- "integer"
  states
}

# ----------------------- Perturbation conditions -----------------------------

condition_definitions <- list(
  Reference = integer(0),
  MALAT1_OFF = c(MALAT1 = 0L),
  miR204_OFF = c(MIR204 = 0L),
  MALAT1_OFF__miR204_OFF = c(MALAT1 = 0L, MIR204 = 0L),
  SIRT1_ON = c(SIRT1 = 1L),
  MALAT1_OFF__SIRT1_ON = c(MALAT1 = 0L, SIRT1 = 1L),
  miR204_ON = c(MIR204 = 1L),
  miR204_ON__GSDME_OFF = c(MIR204 = 1L, DFNA5 = 0L),
  miR204_ON__CASP3_OFF = c(MIR204 = 1L, CASP3 = 0L),
  miR204_ON__SIRT1_ON = c(MIR204 = 1L, SIRT1 = 1L),
  SIRT1_OFF = c(SIRT1 = 0L),
  p53_OFF = c(TP53 = 0L),
  SIRT1_OFF__p53_OFF = c(SIRT1 = 0L, TP53 = 0L),
  p53_ON = c(TP53 = 1L),
  CASP3_OFF = c(CASP3 = 0L),
  p53_ON__CASP3_OFF = c(TP53 = 1L, CASP3 = 0L),
  CASP3_ON = c(CASP3 = 1L),
  GSDME_ON = c(DFNA5 = 1L),
  GSDME_OFF = c(DFNA5 = 0L),
  CASP3_ON__GSDME_OFF = c(CASP3 = 1L, DFNA5 = 0L),
  MALAT1_OFF__GSDME_OFF = c(MALAT1 = 0L, DFNA5 = 0L),
  SIRT1_OFF__GSDME_OFF = c(SIRT1 = 0L, DFNA5 = 0L),
  GSDME_ON__CASP3_ON = c(DFNA5 = 1L, CASP3 = 1L),
  p53_ON__miR204_ON__CASP3_OFF__p21_ON =
    c(TP53 = 1L, MIR204 = 1L, CASP3 = 0L, CDKN1A = 1L),
  MALAT1_ON__miR204_OFF__SIRT1_ON =
    c(MALAT1 = 1L, MIR204 = 0L, SIRT1 = 1L),
  SIRT1_OFF__GSDME_OFF__CASP3_ON =
    c(SIRT1 = 0L, DFNA5 = 0L, CASP3 = 1L),
  BAX_ON = c(BAX = 1L),
  BAX_ON__GSDME_OFF = c(BAX = 1L, DFNA5 = 0L)
)

condition_labels <- c(
  Reference = "Reference: DDR ON, endogenous GSDME",
  MALAT1_OFF = "MALAT1 inhibition",
  miR204_OFF = "miR-204-5p inhibition",
  MALAT1_OFF__miR204_OFF = "MALAT1 inhibition + miR-204-5p inhibition",
  SIRT1_ON = "SIRT1 activation",
  MALAT1_OFF__SIRT1_ON = "MALAT1 inhibition + SIRT1 rescue",
  miR204_ON = "miR-204-5p activation",
  miR204_ON__GSDME_OFF = "miR-204-5p OE + GSDME KO",
  miR204_ON__CASP3_OFF = "miR-204-5p OE + CASP3 KO",
  miR204_ON__SIRT1_ON = "miR-204-5p activation + SIRT1 rescue",
  SIRT1_OFF = "SIRT1 inhibition",
  p53_OFF = "Active p53 inhibition",
  SIRT1_OFF__p53_OFF = "SIRT1 inhibition + active p53 inhibition",
  p53_ON = "Active p53 activation",
  CASP3_OFF = "CASP3 inhibition",
  p53_ON__CASP3_OFF = "Active p53 activation + CASP3 inhibition",
  CASP3_ON = "CASP3 activation",
  GSDME_ON = "GSDME OE (Boolean DFNA5 ON)",
  GSDME_OFF = "GSDME KO",
  CASP3_ON__GSDME_OFF = "CASP3 activation + GSDME KO",
  MALAT1_OFF__GSDME_OFF = "MALAT1 inhibition + GSDME KO",
  SIRT1_OFF__GSDME_OFF = "SIRT1 KO + GSDME KO",
  GSDME_ON__CASP3_ON = "GSDME E1 + CASP3 E1",
  p53_ON__miR204_ON__CASP3_OFF__p21_ON =
    "p53 E1 + miR-204-5p E1 + CASP3 KO + p21 E1",
  MALAT1_ON__miR204_OFF__SIRT1_ON =
    "MALAT1 E1 + miR-204-5p KO + SIRT1 E1",
  SIRT1_OFF__GSDME_OFF__CASP3_ON =
    "SIRT1 KO + GSDME KO + CASP3 E1",
  BAX_ON = "BAX E1",
  BAX_ON__GSDME_OFF = "BAX E1 + GSDME KO"
)
condition_labels <- gsub(" activation", " E1", condition_labels, fixed = TRUE)
condition_labels <- gsub(" inhibition", " KO", condition_labels, fixed = TRUE)
condition_labels <- gsub(" OE", " E1", condition_labels, fixed = TRUE)
condition_labels["CASP3_ON__GSDME_OFF"] <- "GSDME KO + CASP3 E1"
if (!identical(MODEL$input_ids, "DDR") || !identical(BASE_CONTEXT[["DDR"]], 1L))
  stop("This project requires DDR as its sole input, fixed ON (DDR=1)")

condition_table <- data.frame(
  condition = names(condition_definitions),
  label = unname(condition_labels[names(condition_definitions)]),
  clamp = vapply(condition_definitions, function(x) {
    if (!length(x)) return("DDR=1")
    effective <- merge_clamps(BASE_CONTEXT, x)
    paste0(names(effective), "=", effective, collapse = "; ")
  }, character(1)),
  interpretation_scope = "Model-based intervention; not experimental validation",
  stringsAsFactors = FALSE
)
write_table(condition_table, "01_epistasis_conditions.csv")

epistasis_pairs <- data.frame(
  pair_id = c(
    "MALAT1_to_miR204", "miR204_to_SIRT1", "SIRT1_to_p53",
    "p53_to_CASP3", "CASP3_to_GSDME", "MALAT1_to_GSDME_gate",
    "miR204_to_GSDME_gate"
  ),
  upstream = c(
    "MALAT1_OFF", "miR204_ON", "SIRT1_OFF",
    "p53_ON", "CASP3_ON", "MALAT1_OFF", "miR204_ON"
  ),
  downstream = c(
    "miR204_OFF", "SIRT1_ON", "p53_OFF",
    "CASP3_OFF", "GSDME_OFF", "GSDME_OFF", "GSDME_OFF"
  ),
  double = c(
    "MALAT1_OFF__miR204_OFF", "miR204_ON__SIRT1_ON",
    "SIRT1_OFF__p53_OFF", "p53_ON__CASP3_OFF",
    "CASP3_ON__GSDME_OFF", "MALAT1_OFF__GSDME_OFF",
    "miR204_ON__GSDME_OFF"
  ),
  expected_interpretation = c(
    "miR-204-5p inhibition rescues MALAT1 inhibition",
    "SIRT1 activation rescues miR-204-5p activation",
    "Active p53 inhibition blocks SIRT1 inhibition",
    "CASP3 inhibition blocks active p53 activation",
    "GSDME loss blocks CASP3-driven pyroptosis",
    "GSDME loss changes the terminal fate after MALAT1 inhibition",
    "GSDME loss changes the terminal fate after miR-204-5p activation"
  ),
  stringsAsFactors = FALSE
)
write_table(epistasis_pairs, "02_epistasis_pair_definitions.csv")

# ----------------------- Network-wide simulations ----------------------------

simulation_strata <- c("Global random", "Resistant-local")

matched_update_seed <- function(master_seed, stratum_index, replicate_index) {
  value <- (as.double(master_seed) * 100003 + stratum_index * 1009 + replicate_index) %%
    (.Machine$integer.max - 1)
  as.integer(value + 1)
}

simulate_condition_set <- function(model, seeds, conditions, trajectories,
                                   max_steps, jitter_probability) {
  total_blocks <- length(seeds) * length(simulation_strata) * length(conditions)
  output <- vector("list", total_blocks)
  block <- 1L

  for (seed in seeds) {
    log_message("Network epistasis: seed ", seed)
    for (stratum_index in seq_along(simulation_strata)) {
      stratum <- simulation_strata[[stratum_index]]
      initial_states <- make_initial_states(
        model, trajectories, seed + 10000L * stratum_index,
        stratum, jitter_probability
      )

      for (condition in names(conditions)) {
        clamp <- merge_clamps(BASE_CONTEXT, conditions[[condition]])
        states <- matrix(
          0L, nrow = trajectories, ncol = length(model$node_ids),
          dimnames = list(NULL, model$node_ids)
        )
        stable <- logical(trajectories)
        steps <- integer(trajectories)
        cycle_length <- rep(NA_integer_, trajectories)
        fate <- character(trajectories)

        for (replicate_index in seq_len(trajectories)) {
          set.seed(matched_update_seed(seed, stratum_index, replicate_index))
          result <- simulate_to_attractor(
            initial_states[replicate_index, ], model, clamp, max_steps
          )
          states[replicate_index, ] <- if (result$stable) result$state else 0L
          stable[[replicate_index]] <- result$stable
          steps[[replicate_index]] <- result$steps
          cycle_length[[replicate_index]] <- result$cycle_length
          fate[[replicate_index]] <- if (result$stable)
            classify_fate(result$state) else "Other"
        }

        output[[block]] <- data.frame(
          seed = seed,
          initialisation = stratum,
          condition = condition,
          replicate = seq_len(trajectories),
          stable = stable,
          steps = steps,
          cycle_length = cycle_length,
          fate = fate,
          states[, phenotype_nodes, drop = FALSE],
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        block <- block + 1L
      }
    }
  }
  do.call(rbind, output)
}

trajectory_results <- simulate_condition_set(
  MODEL, CFG$seeds, condition_definitions, CFG$trajectories,
  CFG$max_steps, CFG$jitter_probability
)
write_table(trajectory_results, "03_epistasis_trajectory_endpoints.csv")

# ------------------------------ Summaries ------------------------------------

mean_ci <- function(x, conf = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) {
    return(c(n = 0, mean = NA, sd = NA, se = NA, ci_low = NA, ci_high = NA,
             min = NA, max = NA))
  }
  average <- mean(x)
  standard_deviation <- if (n > 1L) stats::sd(x) else NA_real_
  standard_error <- if (n > 1L) standard_deviation / sqrt(n) else NA_real_
  margin <- if (n > 1L) {
    stats::qt(1 - (1 - conf) / 2, df = n - 1L) * standard_error
  } else {
    NA_real_
  }
  c(
    n = n, mean = average, sd = standard_deviation, se = standard_error,
    ci_low = if (n > 1L) average - margin else NA_real_,
    ci_high = if (n > 1L) average + margin else NA_real_,
    min = min(x), max = max(x)
  )
}

summarise_metric <- function(data, group_columns, metric) {
  key <- interaction(data[, group_columns, drop = FALSE], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(data)), key)
  rows <- lapply(groups, function(indices) {
    output <- data[indices[[1]], group_columns, drop = FALSE]
    statistics <- mean_ci(data[[metric]][indices])
    for (name in names(statistics)) output[[paste(metric, name, sep = "_")]] <- statistics[[name]]
    output
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

fate_seed_grid <- expand.grid(
  seed = CFG$seeds,
  initialisation = simulation_strata,
  condition = names(condition_definitions),
  fate = fate_levels,
  stringsAsFactors = FALSE
)
fate_counts <- stats::aggregate(
  replicate ~ seed + initialisation + condition + fate,
  data = trajectory_results,
  FUN = length
)
names(fate_counts)[names(fate_counts) == "replicate"] <- "count"
fate_by_seed <- merge(
  fate_seed_grid, fate_counts,
  by = c("seed", "initialisation", "condition", "fate"), all.x = TRUE
)
fate_by_seed$count[is.na(fate_by_seed$count)] <- 0L
fate_by_seed$frequency <- fate_by_seed$count / CFG$trajectories
fate_by_seed$condition_label <- unname(condition_labels[fate_by_seed$condition])
write_table(fate_by_seed, "04_fate_frequencies_by_seed.csv")

fate_consensus <- summarise_metric(
  fate_by_seed,
  c("initialisation", "condition", "condition_label", "fate"),
  "frequency"
)
fate_consensus$frequency_ci_low <- pmax(0, fate_consensus$frequency_ci_low)
fate_consensus$frequency_ci_high <- pmin(1, fate_consensus$frequency_ci_high)
write_table(fate_consensus, "05_fate_consensus_across_seeds.csv")

phenotype_by_seed <- do.call(rbind, lapply(phenotype_nodes, function(node) {
  aggregate_formula <- stats::as.formula(paste(node, "~ seed + initialisation + condition"))
  output <- stats::aggregate(aggregate_formula, data = trajectory_results, FUN = mean)
  names(output)[names(output) == node] <- "activation_frequency"
  output$node <- node
  output
}))
phenotype_by_seed$condition_label <- unname(condition_labels[phenotype_by_seed$condition])
write_table(phenotype_by_seed, "06_phenotype_activation_by_seed.csv")

phenotype_consensus <- summarise_metric(
  phenotype_by_seed,
  c("initialisation", "condition", "condition_label", "node"),
  "activation_frequency"
)
phenotype_consensus$activation_frequency_ci_low <- pmax(
  0, phenotype_consensus$activation_frequency_ci_low
)
phenotype_consensus$activation_frequency_ci_high <- pmin(
  1, phenotype_consensus$activation_frequency_ci_high
)
write_table(phenotype_consensus, "07_phenotype_activation_consensus.csv")

get_fate_frequency <- function(data, seed, initialisation, condition, fate) {
  value <- data$frequency[
    data$seed == seed & data$initialisation == initialisation &
      data$condition == condition & data$fate == fate
  ]
  if (length(value) != 1L) return(NA_real_)
  value[[1]]
}

epistasis_seed_rows <- list()
row_index <- 1L
for (seed in CFG$seeds) {
  for (initialisation in simulation_strata) {
    for (pair_index in seq_len(nrow(epistasis_pairs))) {
      pair <- epistasis_pairs[pair_index, ]
      p_reference <- get_fate_frequency(fate_by_seed, seed, initialisation, "Reference", "Pyroptosis")
      p_upstream <- get_fate_frequency(
        fate_by_seed, seed, initialisation, pair$upstream, "Pyroptosis"
      )
      p_downstream <- get_fate_frequency(
        fate_by_seed, seed, initialisation, pair$downstream, "Pyroptosis"
      )
      p_double <- get_fate_frequency(
        fate_by_seed, seed, initialisation, pair$double, "Pyroptosis"
      )
      upstream_effect <- p_upstream - p_reference
      rescue_fraction <- if (abs(upstream_effect) < 1e-12) {
        NA_real_
      } else {
        (p_upstream - p_double) / upstream_effect
      }
      distance_to_upstream <- abs(p_double - p_upstream)
      distance_to_downstream <- abs(p_double - p_downstream)
      dominance_denominator <- distance_to_upstream + distance_to_downstream
      downstream_dominance <- if (dominance_denominator < 1e-12) {
        NA_real_
      } else {
        distance_to_upstream / dominance_denominator
      }

      epistasis_seed_rows[[row_index]] <- data.frame(
        seed = seed,
        initialisation = initialisation,
        pair_id = pair$pair_id,
        upstream = pair$upstream,
        downstream = pair$downstream,
        double = pair$double,
        expected_interpretation = pair$expected_interpretation,
        reference_pyroptosis = p_reference,
        upstream_pyroptosis = p_upstream,
        downstream_pyroptosis = p_downstream,
        double_pyroptosis = p_double,
        upstream_effect = upstream_effect,
        rescue_fraction = rescue_fraction,
        downstream_dominance = downstream_dominance,
        downstream_dominant = is.finite(downstream_dominance) && downstream_dominance >= 0.75,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
}
epistasis_by_seed <- do.call(rbind, epistasis_seed_rows)
write_table(epistasis_by_seed, "08_epistasis_and_rescue_by_seed.csv")

epistasis_consensus <- summarise_metric(
  epistasis_by_seed,
  c(
    "initialisation", "pair_id", "upstream", "downstream", "double",
    "expected_interpretation"
  ),
  "rescue_fraction"
)
dominance_consensus <- summarise_metric(
  epistasis_by_seed,
  c("initialisation", "pair_id"),
  "downstream_dominance"
)
epistasis_consensus <- merge(
  epistasis_consensus, dominance_consensus,
  by = c("initialisation", "pair_id"), all.x = TRUE
)
write_table(epistasis_consensus, "09_epistasis_and_rescue_consensus.csv")

gate_pairs <- data.frame(
  gate_test = c("MALAT1 inhibition", "SIRT1 inhibition", "CASP3 activation",
                "miR-204-5p activation", "BAX activation"),
  gsdme_on = c("MALAT1_OFF", "SIRT1_OFF", "CASP3_ON", "miR204_ON",
               "BAX_ON"),
  gsdme_off = c(
    "MALAT1_OFF__GSDME_OFF", "SIRT1_OFF__GSDME_OFF",
    "CASP3_ON__GSDME_OFF", "miR204_ON__GSDME_OFF",
    "BAX_ON__GSDME_OFF"
  ),
  stringsAsFactors = FALSE
)

gate_rows <- list()
row_index <- 1L
for (seed in CFG$seeds) {
  for (initialisation in simulation_strata) {
    for (gate_index in seq_len(nrow(gate_pairs))) {
      gate <- gate_pairs[gate_index, ]
      on_pyro <- get_fate_frequency(
        fate_by_seed, seed, initialisation, gate$gsdme_on, "Pyroptosis"
      )
      off_pyro <- get_fate_frequency(
        fate_by_seed, seed, initialisation, gate$gsdme_off, "Pyroptosis"
      )
      on_apoptosis <- get_fate_frequency(
        fate_by_seed, seed, initialisation, gate$gsdme_on, "Apoptosis"
      )
      off_apoptosis <- get_fate_frequency(
        fate_by_seed, seed, initialisation, gate$gsdme_off, "Apoptosis"
      )
      on_resistance <- get_fate_frequency(
        fate_by_seed, seed, initialisation, gate$gsdme_on, "Resistance"
      )
      off_resistance <- get_fate_frequency(
        fate_by_seed, seed, initialisation, gate$gsdme_off, "Resistance"
      )
      gate_rows[[row_index]] <- data.frame(
        seed = seed,
        initialisation = initialisation,
        gate_test = gate$gate_test,
        gsdme_on_condition = gate$gsdme_on,
        gsdme_off_condition = gate$gsdme_off,
        pyroptosis_with_GSDME = on_pyro,
        pyroptosis_without_GSDME = off_pyro,
        pyroptosis_drop = on_pyro - off_pyro,
        apoptosis_with_GSDME = on_apoptosis,
        apoptosis_without_GSDME = off_apoptosis,
        apoptosis_gain = off_apoptosis - on_apoptosis,
        resistance_with_GSDME = on_resistance,
        resistance_without_GSDME = off_resistance,
        resistance_gain = off_resistance - on_resistance,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
}
gate_by_seed <- do.call(rbind, gate_rows)
write_table(gate_by_seed, "10_GSDME_gate_by_seed.csv")

gate_consensus <- Reduce(
  function(x, y) merge(x, y, by = c("initialisation", "gate_test"), all = TRUE),
  lapply(c("pyroptosis_drop", "apoptosis_gain", "resistance_gain"), function(metric) {
    summarise_metric(gate_by_seed, c("initialisation", "gate_test"), metric)
  })
)
write_table(gate_consensus, "11_GSDME_gate_consensus.csv")

# Matched initial states and random-update seeds across the GSDME pairs.
# Each post-KO output is scored independently among trajectories that reached
# pyroptosis with endogenous GSDME and a fixed point in both conditions.
paired_rows <- list()
paired_index <- 1L
for (seed in CFG$seeds) {
  for (stratum in simulation_strata) {
    for (pair_index in seq_len(nrow(gate_pairs))) {
      gate <- gate_pairs[pair_index, ]
      before <- trajectory_results[
        trajectory_results$seed == seed &
          trajectory_results$initialisation == stratum &
          trajectory_results$condition == gate$gsdme_on, , drop = FALSE
      ]
      after <- trajectory_results[
        trajectory_results$seed == seed &
          trajectory_results$initialisation == stratum &
          trajectory_results$condition == gate$gsdme_off, , drop = FALSE
      ]
      order_after <- match(before$replicate, after$replicate)
      if (anyNA(order_after) || nrow(before) != nrow(after)) {
        stop("Matched GSDME pair lacks replicates: ", gate$gate_test)
      }
      after <- after[order_after, , drop = FALSE]
      eligible <- before$stable & after$stable & before$PYROPTOSIS == 1L
      n_eligible <- sum(eligible)
      for (node in phenotype_nodes) {
        paired_rows[[paired_index]] <- data.frame(
          seed = seed, initialisation = stratum, gate_test = gate$gate_test,
          gsdme_endogenous_condition = gate$gsdme_on,
          gsdme_ko_condition = gate$gsdme_off,
          phenotype = node,
          matched_trajectories = nrow(before),
          pyroptotic_with_GSDME_and_both_stable = n_eligible,
          post_ko_activation_fraction = if (n_eligible > 0L)
            mean(after[[node]][eligible] == 1L) else NA_real_,
          stringsAsFactors = FALSE
        )
        paired_index <- paired_index + 1L
      }
    }
  }
}
paired_gsdme <- do.call(rbind, paired_rows)
write_table(paired_gsdme, "18_GSDME_dependent_trajectories_by_seed.csv")
paired_gsdme_consensus <- summarise_metric(
  paired_gsdme, c("initialisation", "gate_test", "phenotype"),
  "post_ko_activation_fraction"
)
eligible_consensus <- summarise_metric(
  paired_gsdme, c("initialisation", "gate_test", "phenotype"),
  "pyroptotic_with_GSDME_and_both_stable"
)
paired_gsdme_consensus <- merge(
  paired_gsdme_consensus, eligible_consensus,
  by = c("initialisation", "gate_test", "phenotype"), all.x = TRUE
)
write_table(paired_gsdme_consensus,
            "19_GSDME_dependent_trajectories_consensus.csv")

# -------------------- GSE125449 single-cell module ---------------------------

gene_map <- list(
  DFNA5 = c("GSDME", "DFNA5"),
  CASP3 = "CASP3", CASP9 = "CASP9",
  BAX = "BAX", BCL2 = "BCL2", BBC3 = c("BBC3", "PUMA"),
  TP53 = c("TP53", "P53"),
  MIR204 = c("MIR204", "MIR204-5P", "MIR_204_5P"),
  MALAT1 = "MALAT1", CDKN1A = c("CDKN1A", "P21"),
  RB1 = "RB1", E2F1 = "E2F1", SIRT1 = "SIRT1"
)
observational_only_nodes <- c("DFNA5", "CASP3", "CASP9", "TP53", "MIR204")

gse125449_files <- c(
  "GSE125449_Set1_barcodes.tsv.gz", "GSE125449_Set1_genes.tsv.gz",
  "GSE125449_Set1_matrix.mtx.gz", "GSE125449_Set1_samples.txt.gz",
  "GSE125449_Set2_barcodes.tsv.gz", "GSE125449_Set2_genes.tsv.gz",
  "GSE125449_Set2_matrix.mtx.gz", "GSE125449_Set2_samples.txt.gz"
)
gse125449_base_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/"
)

download_gse125449 <- function(cache_dir) {
  target_dir <- file.path(cache_dir, "GSE125449")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", length(gse125449_files))
  for (i in seq_along(gse125449_files)) {
    filename <- gse125449_files[[i]]
    path <- file.path(target_dir, filename)
    cached <- file.exists(path) && file.info(path)$size > 0
    if (!cached) {
      log_message("Downloading GSE125449: ", filename)
      tryCatch(
        utils::download.file(
          paste0(gse125449_base_url, filename), path,
          mode = "wb", method = "libcurl", quiet = FALSE
        ),
        error = function(error) {
          if (file.exists(path)) unlink(path)
          stop("GSE125449 download failed for ", filename, ": ", conditionMessage(error))
        }
      )
    }
    if (!file.exists(path) || file.info(path)$size <= 0) {
      stop("Missing or empty GSE125449 file: ", path)
    }
    rows[[i]] <- data.frame(
      accession = "GSE125449", file = filename,
      local_path = normalizePath(path, mustWork = TRUE),
      size_bytes = file.info(path)$size,
      source = if (cached) "cache" else "downloaded",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

read_gz_table <- function(path, header = TRUE, sep = "\t") {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  utils::read.table(
    connection, header = header, sep = sep, quote = "", comment.char = "",
    check.names = FALSE, fill = TRUE, stringsAsFactors = FALSE
  )
}

read_gz_vector <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  trimws(scan(connection, what = character(), sep = "\n", quiet = TRUE))
}

read_sparse_mtx_gz <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  methods::as(Matrix::readMM(connection), "dgCMatrix")
}

manifest_file <- function(manifest, pattern) {
  matches <- manifest$local_path[grepl(pattern, manifest$file, ignore.case = TRUE)]
  if (length(matches) != 1L) stop("Expected one GEO file matching: ", pattern)
  matches[[1]]
}

extract_network_gene_counts <- function(mat, symbols, mapping) {
  if (nrow(mat) != length(symbols)) stop("Gene annotation and matrix rows differ")
  symbols <- trimws(as.character(symbols))
  upper_symbols <- toupper(symbols)
  wanted <- unique(toupper(unlist(mapping, use.names = FALSE)))
  aliases <- intersect(wanted, unique(upper_symbols))
  if (!length(aliases)) stop("No GSDME-network genes found in GSE125449")
  output <- matrix(
    0, nrow = length(aliases), ncol = ncol(mat),
    dimnames = list(aliases, colnames(mat))
  )
  for (i in seq_along(aliases)) {
    indices <- which(upper_symbols == aliases[[i]])
    output[i, ] <- as.numeric(Matrix::colSums(mat[indices, , drop = FALSE]))
  }
  attr(output, "library_size") <- as.numeric(Matrix::colSums(mat))
  output
}

read_gse125449 <- function(manifest) {
  experiments <- list()
  for (set_name in c("Set1", "Set2")) {
    genes <- read_gz_table(
      manifest_file(manifest, paste0(set_name, "_genes\\.tsv")), header = FALSE
    )
    barcodes <- read_gz_vector(
      manifest_file(manifest, paste0(set_name, "_barcodes\\.tsv"))
    )
    metadata <- read_gz_table(
      manifest_file(manifest, paste0(set_name, "_samples\\.txt"))
    )
    sparse_matrix <- read_sparse_mtx_gz(
      manifest_file(manifest, paste0(set_name, "_matrix\\.mtx"))
    )
    if (ncol(sparse_matrix) != length(barcodes) || nrow(metadata) != ncol(sparse_matrix)) {
      stop("GSE125449 ", set_name, ": matrix/barcode/metadata dimensions differ")
    }
    cell_ids <- make.unique(paste(metadata$Sample, barcodes, sep = "::"))
    colnames(sparse_matrix) <- cell_ids
    metadata$cell_id <- cell_ids
    metadata$sample_id <- metadata$Sample
    metadata$cell_type <- metadata$Type
    metadata$cohort <- set_name
    rownames(metadata) <- cell_ids
    symbols <- if (ncol(genes) >= 2L) genes[[2]] else genes[[1]]
    counts <- extract_network_gene_counts(sparse_matrix, symbols, gene_map)
    experiments[[paste0("GSE125449_", set_name)]] <- list(
      accession = "GSE125449",
      experiment = paste0("GSE125449_", set_name),
      counts = counts,
      metadata = metadata
    )
    rm(sparse_matrix)
    invisible(gc())
  }
  experiments
}

stratified_subset <- function(metadata, max_cells, seed) {
  if (nrow(metadata) <= max_cells) return(seq_len(nrow(metadata)))
  set.seed(seed)
  grouping_columns <- intersect(c("cohort", "sample_id"), names(metadata))
  if (!length(grouping_columns)) return(sort(sample(seq_len(nrow(metadata)), max_cells)))
  groups <- interaction(
    metadata[, grouping_columns, drop = FALSE], drop = TRUE, lex.order = TRUE
  )
  members <- split(seq_len(nrow(metadata)), groups)
  quota <- max(1L, floor(max_cells / length(members)))
  chosen <- unlist(
    lapply(members, function(indices) sample(indices, min(length(indices), quota))),
    use.names = FALSE
  )
  if (length(chosen) < max_cells) {
    remaining <- setdiff(seq_len(nrow(metadata)), chosen)
    chosen <- c(chosen, sample(remaining, min(length(remaining), max_cells - length(chosen))))
  }
  sort(chosen[seq_len(min(length(chosen), max_cells))])
}

fit_two_gaussian_posterior <- function(x, max_iter = 100L, tolerance = 1e-6) {
  x <- as.numeric(x)
  if (length(unique(x)) < 3L || stats::sd(x) < 1e-8) {
    return(ifelse(x > stats::median(x), 0.9, ifelse(x == 0, 0.05, 0.5)))
  }
  means <- as.numeric(stats::quantile(x, c(0.25, 0.75), names = FALSE))
  standard_deviations <- rep(max(stats::sd(x), 0.2), 2L)
  proportions <- c(0.5, 0.5)
  old_likelihood <- -Inf
  posterior <- matrix(0, nrow = length(x), ncol = 2L)
  for (iteration in seq_len(max_iter)) {
    densities <- cbind(
      proportions[[1]] * stats::dnorm(x, means[[1]], standard_deviations[[1]]),
      proportions[[2]] * stats::dnorm(x, means[[2]], standard_deviations[[2]])
    )
    denominator <- rowSums(densities) + 1e-300
    posterior <- densities / denominator
    effective_n <- colSums(posterior) + 1e-8
    proportions <- effective_n / length(x)
    means <- colSums(posterior * x) / effective_n
    standard_deviations <- sqrt(pmax(
      colSums(posterior * (x - rep(means, each = length(x)))^2) / effective_n,
      0.05^2
    ))
    likelihood <- sum(log(denominator))
    if (abs(likelihood - old_likelihood) < tolerance) break
    old_likelihood <- likelihood
  }
  high_component <- which.max(means)
  pmin(pmax(posterior[, high_component], 0.001), 0.999)
}

probabilistic_binarise <- function(counts, mapping) {
  library_size <- attr(counts, "library_size", exact = TRUE)
  if (is.null(library_size) || length(library_size) != ncol(counts)) {
    library_size <- colSums(counts)
  }
  library_size[library_size <= 0] <- 1
  log_normalised <- log1p(t(t(counts) / library_size * 1e4))
  gene_posterior <- matrix(
    NA_real_, nrow = nrow(log_normalised), ncol = ncol(log_normalised),
    dimnames = dimnames(log_normalised)
  )
  for (i in seq_len(nrow(log_normalised))) {
    gene_posterior[i, ] <- fit_two_gaussian_posterior(log_normalised[i, ])
  }
  node_posterior <- matrix(
    NA_real_, nrow = length(mapping), ncol = ncol(counts),
    dimnames = list(names(mapping), colnames(counts))
  )
  upper_rows <- toupper(rownames(gene_posterior))
  for (node in names(mapping)) {
    indices <- which(upper_rows %in% toupper(mapping[[node]]))
    if (length(indices) == 1L) node_posterior[node, ] <- gene_posterior[indices, ]
    if (length(indices) > 1L) {
      node_posterior[node, ] <- apply(gene_posterior[indices, , drop = FALSE], 2, max)
    }
  }
  node_posterior
}

prepare_gse125449_experiment <- function(experiment, max_cells, subset_seed) {
  metadata <- experiment$metadata
  counts <- experiment$counts
  malignant <- grepl("^Malignant", metadata$cell_type, ignore.case = TRUE)
  if (!any(malignant)) stop(experiment$experiment, ": no malignant cells found")
  metadata <- metadata[malignant, , drop = FALSE]
  library_size <- attr(counts, "library_size", exact = TRUE)[malignant]
  counts <- counts[, malignant, drop = FALSE]
  attr(counts, "library_size") <- library_size

  selected <- stratified_subset(metadata, max_cells, subset_seed)
  metadata <- metadata[selected, , drop = FALSE]
  library_size <- attr(counts, "library_size", exact = TRUE)[selected]
  counts <- counts[, selected, drop = FALSE]
  attr(counts, "library_size") <- library_size
  posterior <- probabilistic_binarise(counts, gene_map)

  experiment$metadata <- metadata
  experiment$counts <- counts
  experiment$node_probability <- posterior
  experiment
}

geo_condition_ids <- c(
  "Reference", "MALAT1_OFF", "MALAT1_OFF__miR204_OFF",
  "MALAT1_OFF__SIRT1_ON", "MALAT1_OFF__GSDME_OFF",
  "miR204_ON", "miR204_ON__SIRT1_ON", "miR204_ON__GSDME_OFF"
)

geo_clamp <- function(condition) {
  # The current model has no GSDME-availability input. DFNA5 is only clamped
  # in the explicit GSDME-KO conditions; RNA does not establish cleavage.
  merge_clamps(c(DDR = 1L), condition_definitions[[condition]])
}

simulate_geo_epistasis <- function(experiment, model, seeds, cell_draws, max_steps) {
  cells <- colnames(experiment$node_probability)
  unit_cell_index <- rep(seq_along(cells), each = cell_draws)
  unit_draw <- rep(seq_len(cell_draws), times = length(cells))
  n_units <- length(unit_cell_index)
  output <- vector("list", length(seeds) * length(geo_condition_ids))
  block_index <- 1L

  for (seed in seeds) {
    log_message(experiment$experiment, ": single-cell epistasis seed ", seed)
    initial_states <- matrix(
      0L, nrow = n_units, ncol = length(model$node_ids),
      dimnames = list(NULL, model$node_ids)
    )
    state_seeds <- integer(n_units)
    for (unit_index in seq_len(n_units)) {
      cell_index <- unit_cell_index[[unit_index]]
      draw <- unit_draw[[unit_index]]
      cell <- cells[[cell_index]]
      state_seed <- as.integer(
        (as.double(seed) * 100003 + cell_index * 101 + draw * 17) %%
          (.Machine$integer.max - 1) + 1
      )
      state_seeds[[unit_index]] <- state_seed
      set.seed(state_seed)
      initial_state <- stats::setNames(
        stats::rbinom(length(model$node_ids), 1L, 0.5), model$node_ids
      )
      for (node in setdiff(intersect(rownames(experiment$node_probability),
                                     model$node_ids), observational_only_nodes)) {
        probability <- experiment$node_probability[node, cell]
        if (is.finite(probability)) {
          initial_state[[node]] <- stats::rbinom(1L, 1L, probability)
        }
      }
      initial_state[["DDR"]] <- 1L
      initial_states[unit_index, ] <- initial_state
    }

    for (condition in geo_condition_ids) {
      clamp <- geo_clamp(condition)
      fate <- character(n_units)
      stable <- logical(n_units)
      steps <- integer(n_units)
      for (unit_index in seq_len(n_units)) {
        set.seed(state_seeds[[unit_index]])
        result <- simulate_to_attractor(
          initial_states[unit_index, ], model, clamp, max_steps
        )
        fate[[unit_index]] <- if (result$stable) classify_fate(result$state) else "Other"
        stable[[unit_index]] <- result$stable
        steps[[unit_index]] <- result$steps
      }
      selected_cells <- cells[unit_cell_index]
      metadata_rows <- experiment$metadata[selected_cells, , drop = FALSE]
      output[[block_index]] <- data.frame(
        accession = experiment$accession,
        experiment = experiment$experiment,
        seed = seed,
        draw = unit_draw,
        cell = selected_cells,
        sample_id = as.character(metadata_rows$sample_id),
        cell_type = as.character(metadata_rows$cell_type),
        condition = condition,
        condition_label = unname(condition_labels[[condition]]),
        fate = fate,
        stable = stable,
        steps = steps,
        GSDME_transcript_probability = as.numeric(
          experiment$node_probability["DFNA5", selected_cells]
        ),
        stringsAsFactors = FALSE
      )
      block_index <- block_index + 1L
    }
  }
  do.call(rbind, output)
}

geo_results <- NULL
geo_response_summary <- NULL

if (CFG$geo) {
  manifest <- download_gse125449(CFG$geo_dir)
  write_table(manifest, "12_GSE125449_download_manifest.csv")
  geo_experiments <- read_gse125449(manifest)
  for (i in seq_along(geo_experiments)) {
    geo_experiments[[i]] <- prepare_gse125449_experiment(
      geo_experiments[[i]], CFG$geo_max_cells, CFG$geo_subset_seed + i
    )
  }
  geo_results <- do.call(rbind, lapply(
    geo_experiments,
    simulate_geo_epistasis,
    model = MODEL,
    seeds = CFG$seeds,
    cell_draws = CFG$geo_cell_draws,
    max_steps = CFG$max_steps
  ))
  write_table(geo_results, "13_GSE125449_cell_level_epistasis.csv")

  geo_key <- unique(geo_results[, c("experiment", "seed", "draw", "cell")])
  geo_wide <- geo_key
  for (condition in geo_condition_ids) {
    subset <- geo_results[geo_results$condition == condition, c("experiment", "seed", "draw", "cell", "fate")]
    names(subset)[names(subset) == "fate"] <- condition
    geo_wide <- merge(
      geo_wide, subset,
      by = c("experiment", "seed", "draw", "cell"), all.x = TRUE
    )
  }
  geo_wide$MALAT1_induced_pyroptosis <- with(
    geo_wide, MALAT1_OFF == "Pyroptosis" & Reference != "Pyroptosis"
  )
  geo_wide$anti_miR204_rescue <- with(
    geo_wide,
    MALAT1_OFF == "Pyroptosis" & MALAT1_OFF__miR204_OFF != "Pyroptosis"
  )
  geo_wide$SIRT1_rescue <- with(
    geo_wide,
    MALAT1_OFF == "Pyroptosis" & MALAT1_OFF__SIRT1_ON != "Pyroptosis"
  )
  geo_wide$GSDME_dependency <- with(
    geo_wide,
    MALAT1_OFF == "Pyroptosis" & MALAT1_OFF__GSDME_OFF != "Pyroptosis"
  )
  geo_wide$apoptosis_switch_without_GSDME <- with(
    geo_wide,
    MALAT1_OFF == "Pyroptosis" & MALAT1_OFF__GSDME_OFF == "Apoptosis"
  )
  geo_wide$miR204_GSDME_dependency <- with(
    geo_wide,
    miR204_ON == "Pyroptosis" & miR204_ON__GSDME_OFF != "Pyroptosis"
  )
  write_table(geo_wide, "14_GSE125449_paired_cell_responses.csv")

  response_metrics <- c(
    "MALAT1_induced_pyroptosis", "anti_miR204_rescue", "SIRT1_rescue",
    "GSDME_dependency", "apoptosis_switch_without_GSDME",
    "miR204_GSDME_dependency"
  )
  geo_response_by_seed <- do.call(rbind, lapply(response_metrics, function(metric) {
    output <- stats::aggregate(
      stats::as.formula(paste(metric, "~ experiment + seed")),
      data = geo_wide,
      FUN = mean
    )
    names(output)[names(output) == metric] <- "fraction"
    output$metric <- metric
    output
  }))
  write_table(geo_response_by_seed, "15_GSE125449_response_fractions_by_seed.csv")
  geo_response_summary <- summarise_metric(
    geo_response_by_seed, c("experiment", "metric"), "fraction"
  )
  geo_response_summary$fraction_ci_low <- pmax(0, geo_response_summary$fraction_ci_low)
  geo_response_summary$fraction_ci_high <- pmin(1, geo_response_summary$fraction_ci_high)
  geo_response_summary$claim_scope <- paste(
    "Model-projected response from observational malignant-cell states;",
    "not causal single-cell validation"
  )
  write_table(geo_response_summary, "16_GSE125449_response_consensus.csv")
}

# --------------------------- Publication figures -----------------------------

fate_colours <- c(
  Pyroptosis = "#0072B2", Apoptosis = "#D55E00", Resistance = "#6F4AB8",
  Proliferation = "#009E73", `Cell-cycle arrest` = "#E69F00",
  Other = "#8A8A8A"
)

main_conditions <- c(
  "Reference", "MALAT1_OFF", "MALAT1_OFF__miR204_OFF",
  "MALAT1_OFF__SIRT1_ON", "MALAT1_OFF__GSDME_OFF",
  "miR204_ON", "miR204_ON__SIRT1_ON", "miR204_ON__GSDME_OFF",
  "miR204_ON__CASP3_OFF", "SIRT1_OFF", "SIRT1_OFF__GSDME_OFF",
  "SIRT1_OFF__p53_OFF", "p53_ON", "p53_ON__CASP3_OFF",
  "CASP3_ON", "CASP3_ON__GSDME_OFF", "GSDME_ON", "GSDME_OFF",
  "GSDME_ON__CASP3_ON", "SIRT1_OFF__GSDME_OFF__CASP3_ON",
  "p53_ON__miR204_ON__CASP3_OFF__p21_ON",
  "MALAT1_ON__miR204_OFF__SIRT1_ON", "BAX_ON", "BAX_ON__GSDME_OFF"
)

heatmap_data <- phenotype_consensus[
  phenotype_consensus$condition %in% main_conditions &
    phenotype_consensus$initialisation == "Resistant-local",
]
heatmap_data$condition_label <- factor(
  heatmap_data$condition_label,
  levels = rev(unname(condition_labels[main_conditions]))
)
heatmap_data$phenotype_label <- factor(
  unname(c(PYROPTOSIS = "Pyroptosis", APOPTOSIS = "Apoptosis",
           RESISTANCE = "Resistance", PROLIFERATION = "Proliferation",
           CELL_CYCLE_ARREST = "Cell-cycle arrest")[heatmap_data$node]),
  levels = c("Pyroptosis", "Apoptosis", "Resistance",
             "Proliferation", "Cell-cycle arrest")
)

panel_a <- ggplot(
  heatmap_data,
  aes(x = phenotype_label, y = condition_label, fill = activation_frequency_mean)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", activation_frequency_mean)), size = 2.6) +
  scale_fill_gradient(
    low = "#F7FBFF", high = "#B2182B", limits = c(0, 1), name = "Frequency"
  ) +
  labs(
    title = "A  Five independent phenotype outputs",
    subtitle = "Matched asynchronous perturbations; outputs may coexist",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold")
  )

rescue_plot_data <- epistasis_consensus[
  epistasis_consensus$initialisation == "Resistant-local",
]
rescue_plot_data$pair_label <- factor(
  rescue_plot_data$pair_id,
  levels = rev(epistasis_pairs$pair_id),
  labels = rev(c(
    "MALAT1 -> miR-204-5p", "miR-204-5p -> SIRT1", "SIRT1 -> active p53",
    "Active p53 -> CASP3", "CASP3 -> GSDME", "MALAT1 -> GSDME gate",
    "miR-204-5p -> GSDME gate"
  ))
)

panel_b <- ggplot(
  rescue_plot_data,
  aes(x = rescue_fraction_mean, y = pair_label)
) +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.4) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#0072B2", linewidth = 0.5) +
  geom_errorbarh(
    aes(xmin = rescue_fraction_ci_low, xmax = rescue_fraction_ci_high),
    height = 0.18, colour = "#333333", na.rm = TRUE
  ) +
  geom_point(size = 3.0, colour = "#009E73") +
  coord_cartesian(xlim = c(-0.25, 1.25)) +
  labs(
    title = "B  Model-based rescue",
    subtitle = "Mean and 95% CI across random seeds; 1 indicates complete rescue",
    x = "Rescue fraction", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

gate_plot_conditions <- c("MALAT1_OFF", "MALAT1_OFF__GSDME_OFF")
gate_plot_data <- fate_consensus[
  fate_consensus$initialisation == "Resistant-local" &
    fate_consensus$condition %in% gate_plot_conditions,
]
gate_plot_data$condition_label <- factor(
  gate_plot_data$condition_label,
  levels = unname(condition_labels[gate_plot_conditions])
)

panel_c <- ggplot(
  gate_plot_data,
  aes(x = condition_label, y = frequency_mean, fill = fate)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = fate_colours, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    title = "C  GSDME terminal-gate test",
    subtitle = "Loss of GSDME tests a pyroptosis-to-apoptosis/resistance switch",
    x = NULL, y = "Mean fate frequency", fill = "Fate"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom", plot.title = element_text(face = "bold")
  )

if (!is.null(geo_response_summary) && nrow(geo_response_summary)) {
  response_labels <- c(
    MALAT1_induced_pyroptosis = "MALAT1 inhibition\ninduces pyroptosis",
    anti_miR204_rescue = "Anti-miR-204-5p\nrescue",
    SIRT1_rescue = "SIRT1\nrescue",
    GSDME_dependency = "GSDME\ndependency",
    apoptosis_switch_without_GSDME = "Apoptosis switch\nwithout GSDME",
    miR204_GSDME_dependency = "miR-204-5p / GSDME\ndependency"
  )
  geo_response_summary$metric_label <- unname(response_labels[geo_response_summary$metric])
  geo_response_summary$experiment_label <- gsub("_", " ", geo_response_summary$experiment)
  panel_d <- ggplot(
    geo_response_summary,
    aes(x = metric_label, y = fraction_mean, fill = experiment_label)
  ) +
    geom_col(position = position_dodge(width = 0.78), width = 0.70) +
    geom_errorbar(
      aes(ymin = fraction_ci_low, ymax = fraction_ci_high),
      position = position_dodge(width = 0.78), width = 0.18, na.rm = TRUE
    ) +
    scale_y_continuous(
      limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")
    ) +
    scale_fill_manual(values = c("#0072B2", "#E69F00")) +
    labs(
      title = "D  Single-cell-informed response potential",
      subtitle = "GSE125449 malignant-cell states; observational support only",
      x = NULL, y = "Projected cell fraction", fill = "Dataset"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "bottom", plot.title = element_text(face = "bold")
    )
  main_figure <- (panel_a / panel_b) | (panel_c / panel_d)
} else {
  main_figure <- panel_a / (panel_b | panel_c)
}

main_figure <- main_figure + patchwork::plot_annotation(
  title = "In silico epistasis and rescue analysis of the GSDME cell-fate network",
  subtitle = paste(
    "Model-based causal predictions across", length(CFG$seeds),
    "random seeds; experimental validation remains required"
  ),
  theme = theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 11)
  )
)
save_publication_plot(
  main_figure, "Figure_08_in_silico_epistasis_and_rescue", 16, 12
)

fixed_conditions <- c("Reference", "miR204_ON", "miR204_ON__GSDME_OFF")
fixed_plot <- phenotype_consensus[
  phenotype_consensus$condition %in% fixed_conditions, , drop = FALSE
]
fixed_plot$condition_label <- factor(
  fixed_plot$condition_label,
  levels = unname(condition_labels[fixed_conditions])
)
phenotype_labels <- c(
  PYROPTOSIS = "GSDME-mediated pyroptosis", APOPTOSIS = "Apoptosis",
  RESISTANCE = "Resistance", PROLIFERATION = "Proliferation",
  CELL_CYCLE_ARREST = "Cell-cycle arrest"
)
fixed_plot$phenotype_label <- factor(
  unname(phenotype_labels[fixed_plot$node]),
  levels = unname(phenotype_labels[phenotype_nodes])
)
p_fixed <- ggplot(fixed_plot,
                  aes(x = phenotype_label, y = condition_label,
                      fill = activation_frequency_mean)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * activation_frequency_mean)), size = 3.6) +
  facet_wrap(~ initialisation, nrow = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                      labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "miR-204-5p E1 + GSDME KO across five Boolean outputs",
       subtitle = "DDR ON; average of five computational seeds under matched initial states",
       caption = paste0("Outputs are independently scored and may coexist; seed intervals are in Table 07. ",
                        "The compound intervention is simulated, not experimentally applied to GEO cells."),
       x = NULL, y = NULL, fill = "Output ON") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))
save_publication_plot(
  p_fixed, "Figure_09_miR204_OE_GSDME_KO_five_phenotypes", 15, 6
)

# Direct dependency panel: the same five outputs under matched perturbations.
gsdme_condition_ids <- c(
  "Reference", "GSDME_ON", "GSDME_OFF", "CASP3_ON",
  "CASP3_ON__GSDME_OFF", "CASP3_OFF", "miR204_ON",
  "miR204_ON__GSDME_OFF", "miR204_ON__CASP3_OFF",
  "SIRT1_OFF", "SIRT1_OFF__GSDME_OFF",
  "MALAT1_OFF", "MALAT1_OFF__GSDME_OFF",
  "GSDME_ON__CASP3_ON", "SIRT1_OFF__GSDME_OFF__CASP3_ON",
  "p53_ON__miR204_ON__CASP3_OFF__p21_ON",
  "MALAT1_ON__miR204_OFF__SIRT1_ON", "BAX_ON", "BAX_ON__GSDME_OFF"
)
dependency_data <- phenotype_consensus[
  phenotype_consensus$condition %in% gsdme_condition_ids, , drop = FALSE
]
dependency_data$condition_label <- factor(
  unname(condition_labels[dependency_data$condition]),
  levels = rev(unname(condition_labels[gsdme_condition_ids]))
)
dependency_data$phenotype_label <- factor(
  unname(phenotype_labels[dependency_data$node]),
  levels = unname(phenotype_labels[phenotype_nodes])
)
p_gsdme_dependency <- ggplot(
  dependency_data,
  aes(x = phenotype_label, y = condition_label,
      fill = activation_frequency_mean)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * activation_frequency_mean)),
            size = 2.9) +
  facet_wrap(~ initialisation, nrow = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                      labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "GSDME requirement for model-predicted pyroptosis",
       subtitle = "DDR ON | five independent Boolean outputs | five computational seeds",
       caption = paste("Paired initial states; seed uncertainty is in Table 07.",
                       "GSDME E1 is a direct DFNA5 clamp, not a cleavage assay."),
       x = NULL, y = NULL, fill = "Output ON") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))
save_publication_plot(p_gsdme_dependency,
                      "Figure_10_GSDME_dependency_DDR_ON_five_phenotypes", 18, 9)

paired_plot <- paired_gsdme_consensus[
  paired_gsdme_consensus$post_ko_activation_fraction_n > 0L, , drop = FALSE
]
if (nrow(paired_plot)) {
  paired_plot$phenotype_label <- factor(
    unname(phenotype_labels[paired_plot$phenotype]),
    levels = unname(phenotype_labels[phenotype_nodes])
  )
  paired_plot$gate_test <- factor(
    paired_plot$gate_test, levels = rev(gate_pairs$gate_test)
  )
  p_paired_gsdme <- ggplot(
    paired_plot,
    aes(x = phenotype_label, y = gate_test,
        fill = post_ko_activation_fraction_mean)
  ) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * post_ko_activation_fraction_mean)),
              size = 3.3) +
    facet_wrap(~ initialisation, nrow = 1) +
    scale_fill_gradient(low = "#F7FBFF", high = "#D55E00", limits = c(0, 1),
                        labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "After GSDME KO: fates of matched pyroptotic trajectories",
         subtitle = "DDR ON | conditional on pyroptosis before KO and stable endpoints in both runs",
         caption = paste("Outputs may coexist. Denominators and number of eligible seeds are in Table 19.",
                         "No eligible trajectories are omitted, not scored as zero."),
         x = NULL, y = NULL, fill = "Output ON") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold"))
  save_publication_plot(p_paired_gsdme,
                        "Figure_11_GSDME_KO_matched_trajectory_fates", 15, 6)
}

seed_robustness <- fate_by_seed[
  fate_by_seed$initialisation == "Resistant-local" &
    fate_by_seed$fate == "Pyroptosis" &
    fate_by_seed$condition %in% main_conditions,
]
seed_robustness$condition_label <- factor(
  seed_robustness$condition_label,
  levels = rev(unname(condition_labels[main_conditions]))
)
seed_robustness$seed_label <- factor(
  paste0("Seed ", seed_robustness$seed),
  levels = paste0("Seed ", CFG$seeds)
)

supplementary_figure <- ggplot(
  seed_robustness,
  aes(x = seed_label, y = condition_label, fill = frequency)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", frequency)), size = 2.8) +
  scale_fill_gradient(
    low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
    name = "Pyroptosis\nfrequency"
  ) +
  labs(
    title = "Seed robustness of model-predicted GSDME-mediated pyroptosis",
    subtitle = "Resistant-local initial states; computational seeds are not biological replicates",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(12, 18, 12, 18)
  )
save_publication_plot(
  supplementary_figure,
  "Supplementary_Figure_S2_epistasis_seed_robustness",
  11, 8
)

# ------------------------- Manifest and provenance ---------------------------

manifest <- data.frame(
  script_version = SCRIPT_VERSION,
  model = normalizePath(CFG$model, mustWork = TRUE),
  model_md5 = unname(tools::md5sum(CFG$model)),
  seeds = paste(CFG$seeds, collapse = ","),
  trajectories_per_condition_and_stratum = CFG$trajectories,
  max_steps = CFG$max_steps,
  jitter_probability = CFG$jitter_probability,
  quick_mode = CFG$quick,
  geo_enabled = CFG$geo,
  geo_max_cells_per_set = if (CFG$geo) CFG$geo_max_cells else NA_integer_,
  geo_cell_draws = if (CFG$geo) CFG$geo_cell_draws else NA_integer_,
  geo_claim_scope = if (CFG$geo) {
    "Observational malignant-cell initial states; not perturbational validation"
  } else {
    "GEO module disabled"
  },
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
write_table(manifest, "17_analysis_manifest.csv")
capture.output(sessionInfo(), file = file.path(LOG_DIR, "sessionInfo.txt"))

log_message(
  "Completed. Main figure: ",
  file.path(FIG_DIR, "Figure_08_in_silico_epistasis_and_rescue_600dpi.png")
)
if (CFG$quick) {
  log_message("WARNING: quick-mode numerical results must not be reported in the manuscript")
}
