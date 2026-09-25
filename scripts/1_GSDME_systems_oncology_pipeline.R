#!/usr/bin/env Rscript

# =============================================================================
# GSDME Systems-Oncology Boolean Network Pipeline
# =============================================================================
# Modules:
#   1. GINsim import, network reconstruction and attractor exploration
#   2. Sustained attractor control / minimal driver-node search
#   3. Probabilistic scRNA-seq binarisation and cell-specific simulations
#   4. GSE125449/GSE189903 observational validation with HCC-only safeguards
#   5. Multi-omic soft-clipping for patient-specific network instances
#   6. Tabular Q-learning for sequential interventions
#   7. Publication figures (PDF + 600-dpi PNG) and CSV tables
#
# The script has a reproducible demo mode. Demo scRNA-seq and multi-omic data
# are simulated and MUST NOT be presented as experimental validation.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

PIPELINE_VERSION <- "1.6.1-HCC-GSDME-switch-DDR-ON"

# ------------------------------- CLI -----------------------------------------

parse_cli <- function(args) {
  out <- list(
    model = "GINsim-miR_204_GSDME_Pyroptosis.zginml",
    outdir = NULL,
    scrna = NA_character_,
    geo = c("GSE125449", "GSE189903"),
    geo_dir = "GEO_scRNA_data",
    geo_platform = "droplet",
    geo_max_cells = 2000L,
    multiomics = NA_character_,
    seed = getOption("gsdme.seed", 101L),
    demo = FALSE,
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
        "  --geo ACCESSIONS    Comma-separated GEO accessions (supported: GSE125449,GSE189903)\n",
        "  --geo-dir PATH      Download/cache directory for GEO files (default GEO_scRNA_data)\n",
        "  --geo-platform X    Compatibility option; ignored for the two HCC series\n",
        "  --geo-max-cells N   Maximum analysed cells per GEO experiment (default 2000)\n",
        "  --multiomics PATH   Long CSV/TSV patient multi-omics table\n",
        "  --seed INTEGER      Random seed (default 101; also options(gsdme.seed = 101L) in R)\n",
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
# Ao usar source() no R, procurar o modelo ao lado do script se necessario.
if (!file.exists(CFG$model)) {
  source_files <- lapply(sys.frames(), function(frame) frame$ofile)
  source_files <- Filter(function(path) !is.null(path) && length(path) == 1L,
                         source_files)
  if (length(source_files)) {
    local_model <- file.path(dirname(normalizePath(tail(source_files, 1L)[[1]])), CFG$model)
    if (file.exists(local_model)) CFG$model <- local_model
  }
}
if (is.null(CFG$outdir)) CFG$outdir <- paste0("resultados_finais_GSDME_seed_", CFG$seed)
if (!setequal(CFG$geo, c("GSE125449", "GSE189903")))
  stop("A analise HCC exige apenas GSE125449 e GSE189903 juntos")
if (!CFG$geo_platform %in% c("droplet", "smartseq2", "all")) {
  stop("--geo-platform must be droplet, smartseq2 or all")
}
if (is.na(CFG$geo_max_cells) || CFG$geo_max_cells < 20L) {
  stop("--geo-max-cells must be an integer >= 20")
}
set.seed(CFG$seed)

# ---------------------------- Dependencies -----------------------------------

required_packages <- c("xml2", "igraph", "ggplot2")
if (length(CFG$geo)) required_packages <- unique(c(required_packages, "Matrix"))
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
log_message("Loaded model '", MODEL$id, "': ", nrow(MODEL$nodes), " nodes, ",
            nrow(MODEL$edges), " edges, input(s): ", paste(MODEL$input_ids, collapse = ", "))
write_table(MODEL$nodes, "01_model_nodes_and_rules.csv")
write_table(MODEL$edges, "02_model_edges.csv")
if (!identical(MODEL$input_ids, "DDR"))
  stop("Expected DDR as the sole input; review the current GINsim model")
g_rule <- MODEL$nodes$rule[match("DFNA5", MODEL$nodes$name)]
a_rule <- MODEL$nodes$rule[match("APOPTOSIS", MODEL$nodes$name)]
if (identical(g_rule, "CASP3") && identical(a_rule, "CASP3 & !GSDME"))
  log_message("MODEL ALERT: DFNA5 = CASP3 and APOPTOSIS = CASP3 & !GSDME; ",
              "without DFNA5 intervention, apoptosis is zero at a Boolean fixed point.")

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
  if (length(clamp)) {
    if (is.null(names(clamp)) || anyNA(names(clamp)) || any(!nzchar(names(clamp)))) {
      stop("Every clamped value must have a valid node name")
    }
    unknown <- setdiff(names(clamp), names(state))
    if (length(unknown)) stop("Clamp refers to unknown node(s): ", paste(unknown, collapse = ", "))
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
                                  mode = "asynchronous", max_steps = 300L,
                                  return_trajectory = FALSE) {
  state <- apply_clamp(normalise_state(initial_state, model), clamp)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  trajectory <- if (return_trajectory) list(state) else NULL
  for (step in seq_len(max_steps)) {
    key <- encode_state(state, model)
    if (mode == "synchronous" && exists(key, seen, inherits = FALSE)) {
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
  c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE", "PROLIFERATION", "CELL_CYCLE_ARREST"),
  MODEL$node_ids
)

classify_fate <- function(state) {
  is_on <- function(x) x %in% names(state) && state[[x]] == 1L
  if (is_on("PYROPTOSIS")) return("Pyroptosis")
  if (is_on("APOPTOSIS")) return("Apoptosis")
  if (is_on("RESISTANCE")) return("Resistance")
  if (is_on("PROLIFERATION")) return("Proliferation")
  if (is_on("CELL_CYCLE_ARREST")) return("Cell-cycle arrest")
  "Other"
}

make_reference_state <- function(model, type = c("malignant", "death_primed")) {
  type <- match.arg(type)
  s <- empty_state(model)
  s[intersect(c("DDR"), names(s))] <- 1L
  if (type == "malignant") {
    s[intersect(c("MALAT1", "SIRT1",
                  "RESISTANCE", "(CCND1 ou CCND2 ou CCND3) E (CDK4 ou CDK6)", "E2F1"), names(s))] <- 1L
  } else {
    s[intersect(c("MIR204", "p53-Arrest", "BBC3", "BAX",
                  "CASP9", "CASP3", "DFNA5",
                  "PYROPTOSIS", "CDKN1A", "RB1"), names(s))] <- 1L
  }
  s
}

BASE_CLAMP <- c(DDR = 1L)
BASE_CLAMP <- BASE_CLAMP[names(BASE_CLAMP) %in% MODEL$node_ids]
MALIGNANT_STATE <- simulate_to_attractor(
  make_reference_state(MODEL, "malignant"), MODEL,
  clamp = BASE_CLAMP, mode = "synchronous"
)$state
DEATH_STATE <- simulate_to_attractor(
  make_reference_state(MODEL, "death_primed"), MODEL,
  clamp = BASE_CLAMP, mode = "synchronous"
)$state

reference_table <- rbind(
  data.frame(reference = "Malignant/resistant", t(MALIGNANT_STATE), check.names = FALSE),
  data.frame(reference = "Death/pyroptosis", t(DEATH_STATE), check.names = FALSE)
)
write_table(reference_table, "03_reference_attractors.csv")

# ----------------------- Publication network plot ----------------------------

pretty_labels <- c(
  DFNA5 = "GSDME (DFNA5)",
  CASP3 = "CASP3", CASP9 = "CASP9",
  MOMP = "MOMP", BAX = "BAX", BCL2 = "BCL2", BBC3 = "BBC3",
  TP53 = "p53 (TP53)", MIR204 = "miR-204-5p (MIR204)",
  MALAT1 = "MALAT1", PRKAA1 = "PRKAA1 (AMPK)", PPM1D = "PPM1D (Wip1)",
  CDKN1A = "CDKN1A", `(CCND1 ou CCND2 ou CCND3) E (CDK4 ou CDK6)` = "Cyclin D-CDK4/6", RB1 = "RB1",
  E2F1 = "E2F1", DDR = "DDR",
  PYROPTOSIS = "GSDME-mediated\npyroptosis",
  RESISTANCE = "Resistance", PROLIFERATION = "Proliferation",
  APOPTOSIS = "Apoptosis", CELL_CYCLE_ARREST = "Cell-cycle arrest",
  SIRT1 = "SIRT1", `p53-Arrest` = "p53-Arrest", `p53-Killer` = "p53-Killer"
)

display_node <- function(x) {
  out <- ifelse(x %in% names(pretty_labels), unname(pretty_labels[x]), gsub("_", " ", x))
  gsub("\n", " ", out, fixed = TRUE)
}

display_action <- function(x) {
  labels <- c(
    NONE = "No intervention", MALAT1_KO = "MALAT1 inhibition",
    miR204_OE = "miR-204-5p activation", SIRT1_KO = "SIRT1 inhibition",
    BCL2_KO = "BCL2 inhibition"
  )
  ifelse(x %in% names(labels), unname(labels[x]), gsub("_", " ", x))
}

display_perturbation <- function(x) {
  labels <- c(
    `Reference: DDR ON, GSDME endogenous` = "Reference\nDDR ON (GSDME endogenous)",
    `GSDME KO` = "GSDME loss\n(KO)", `CASP3 KO` = "CASP3 inhibition\n(KO)",
    `GSDME OE` = "GSDME activation\n(E1)", `CASP3 OE` = "CASP3 activation\n(E1)",
    `miR-204-5p OE + GSDME KO` = "miR-204-5p E1\n+ GSDME KO",
    `miR-204-5p OE` = "miR-204-5p activation\n(E1)",
    `miR-204-5p KO` = "miR-204-5p inhibition\n(KO)",
    `MALAT1 OE` = "MALAT1 activation\n(OE)",
    `MALAT1 KO` = "MALAT1 inhibition\n(KO)",
    `SIRT1 OE` = "SIRT1 activation\n(OE)",
    `SIRT1 KO` = "SIRT1 inhibition\n(KO)",
    `BCL2 KO` = "BCL2 inhibition\n(KO)"
  )
  gsub("OE", "E1", ifelse(x %in% names(labels), unname(labels[x]), x), fixed = TRUE)
}

display_condition_e1 <- function(x) {
  x <- gsub(" OE", " E1", x, fixed = TRUE)
  x[x == "CASP3 E1 + GSDME KO"] <- "GSDME KO + CASP3 E1"
  x
}

display_driver_intervention <- function(x) {
  parts <- strsplit(x, " \\+ ")[[1]]
  formatted <- vapply(parts, function(part) {
    fields <- strsplit(part, "=", fixed = TRUE)[[1]]
    node <- display_node(fields[[1]])
    state <- if (length(fields) > 1L) fields[[2]] else ""
    paste0(node, if (identical(state, "1")) " activation (ON)" else " inhibition (OFF)")
  }, character(1))
  paste(formatted, collapse = " + ")
}

manual_layout <- data.frame(
  id = c("DFNA5", "PYROPTOSIS", "APOPTOSIS", "CASP3", "CASP9", "MOMP",
         "BAX", "BBC3", "BCL2", "MALAT1", "MIR204", "SIRT1", "DDR",
         "TP53", "CDKN1A", "(CCND1 ou CCND2 ou CCND3) E (CDK4 ou CDK6)",
         "RB1", "E2F1", "RESISTANCE", "PROLIFERATION", "CELL_CYCLE_ARREST"),
  x = c(5, 2.7, 7.3, 5, 5, 5, 5, 3.8, 2.6, 7.5, 7.5, 7.5, 9.3,
        9.3, 11, 11, 11, 11, 9.5, 7.2, 4.1),
  y = c(11, 10.4, 10.4, 9.7, 8.6, 6.4, 5.3, 4.4, 5.3, 10.1, 9.0, 7.9,
        5.3, 7.8, 7.0, 5.9, 4.8, 3.7, 2.2, 2.2, 2.2),
  stringsAsFactors = FALSE
)

plot_network <- function(model) {
  layout <- merge(data.frame(id = model$node_ids), manual_layout, by = "id", all.x = TRUE)
  missing <- which(is.na(layout$x) | is.na(layout$y))
  if (length(missing)) {
    fallback <- igraph::layout_with_fr(
      igraph::graph_from_data_frame(model$edges[, c("from", "to")],
                                    vertices = data.frame(name = model$nodes$name))
    )
    rownames(fallback) <- model$nodes$name
    range_to_layout <- function(x) if (diff(range(x)) == 0) rep(6, length(x)) else
      1 + 11 * (x - min(x)) / diff(range(x))
    layout$x[missing] <- range_to_layout(fallback[layout$id[missing], 1])
    layout$y[missing] <- range_to_layout(fallback[layout$id[missing], 2])
  }
  layout$label <- ifelse(layout$id %in% names(pretty_labels),
                         unname(pretty_labels[layout$id]), layout$id)
  outputs <- c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE",
               "PROLIFERATION", "CELL_CYCLE_ARREST")
  priority_axis <- c("MIR204", "SIRT1", "MALAT1", "TP53", "BAX", "DFNA5", "CASP3")
  layout$class <- ifelse(layout$id %in% priority_axis, "Priority node",
                         ifelse(layout$id %in% model$input_ids, "Input",
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
               size = 3.0, linewidth = 0.35, label.padding = grid::unit(0.17, "lines"),
               colour = "#17212B", fontface = "bold") +
    scale_fill_manual(values = c(Input = "#73DDE1", `Fixed condition` = "#F2C14E",
                                 `Priority node` = "#F4B46B",
                                 Phenotype = "#D7D7D7", Regulator = "white")) +
    scale_colour_manual(values = c(Activation = "#007A3D", Inhibition = "#D55E00")) +
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
save_publication_plot(network_plot, "Figure_01_GSDME_logical_network", 14, 10)

# -------------------- Literature-grounded perturbations ----------------------

scientific_perturbations <- list(
  `Reference: DDR ON, GSDME endogenous` = c(DDR = 1L),
  `GSDME KO` = c(DDR = 1L, DFNA5 = 0L),
  `CASP3 KO` = c(DDR = 1L, CASP3 = 0L),
  `miR-204-5p OE` = c(DDR = 1L, MIR204 = 1L),
  `miR-204-5p KO` = c(DDR = 1L, MIR204 = 0L),
  `MALAT1 OE` = c(DDR = 1L, MALAT1 = 1L),
  `MALAT1 KO` = c(DDR = 1L, MALAT1 = 0L),
  `SIRT1 OE` = c(DDR = 1L, SIRT1 = 1L),
  `SIRT1 KO` = c(DDR = 1L, SIRT1 = 0L),
  `BCL2 KO` = c(DDR = 1L, BCL2 = 0L),
  `GSDME OE` = c(DDR = 1L, DFNA5 = 1L),
  `CASP3 OE` = c(DDR = 1L, CASP3 = 1L),
  `miR-204-5p OE + GSDME KO` = c(DDR = 1L, MIR204 = 1L, DFNA5 = 0L),
  `CASP3 OE + GSDME KO` = c(DDR = 1L, CASP3 = 1L, DFNA5 = 0L),
  `SIRT1 KO + GSDME KO` = c(DDR = 1L, SIRT1 = 0L, DFNA5 = 0L),
  `MALAT1 KO + GSDME KO` = c(DDR = 1L, MALAT1 = 0L, DFNA5 = 0L),
  `miR-204-5p OE + CASP3 KO` = c(DDR = 1L, MIR204 = 1L, CASP3 = 0L),
  `GSDME E1 + CASP3 E1` = c(DDR = 1L, DFNA5 = 1L, CASP3 = 1L),
  `p53 E1 + miR-204-5p E1 + CASP3 KO + p21 E1` =
    c(DDR = 1L, TP53 = 1L, MIR204 = 1L, CASP3 = 0L, CDKN1A = 1L),
  `MALAT1 E1 + miR-204-5p KO + SIRT1 E1` =
    c(DDR = 1L, MALAT1 = 1L, MIR204 = 0L, SIRT1 = 1L),
  `SIRT1 KO + GSDME KO + CASP3 E1` =
    c(DDR = 1L, SIRT1 = 0L, DFNA5 = 0L, CASP3 = 1L),
  `BAX E1` = c(DDR = 1L, BAX = 1L),
  `BAX E1 + GSDME KO` = c(DDR = 1L, BAX = 1L, DFNA5 = 0L),
  `p53 E1` = c(DDR = 1L, TP53 = 1L),
  `p53 KO` = c(DDR = 1L, TP53 = 0L),
  `p21 E1` = c(DDR = 1L, CDKN1A = 1L)
)
if (!all(vapply(scientific_perturbations, function(clamp)
  "DDR" %in% names(clamp) && identical(unname(clamp[["DDR"]]), 1L),
  logical(1)))) stop("Every perturbation must hold DDR ON (DDR=1)")

perturbation_evidence <- data.frame(
  perturbation = names(scientific_perturbations),
  expected_direction = c(
    "Reference condition (endogenous DFNA5)", "Direct DFNA5 loss suppresses model pyroptosis",
    "Reduced model DFNA5 and pyroptosis", "SIRT1 suppression; anti-tumour direction",
    "SIRT1 release; pro-survival direction", "miR-204 sequestration and SIRT1 release",
    "miR-204 release and SIRT1 suppression", "Pro-survival/resistance direction",
    "Reduced survival; mitochondrial death sensitisation", "BAX disinhibition",
    "DFNA5 clamped ON; model projects GSDME-dependent pyroptosis",
    "CASP3 clamped ON; model projects DFNA5 activation and pyroptosis",
    "miR-204-5p forced ON with GSDME Boolean node forced OFF",
    "CASP3 forced ON with GSDME Boolean node forced OFF",
    "SIRT1 forced OFF with GSDME Boolean node forced OFF",
    "MALAT1 forced OFF with GSDME Boolean node forced OFF",
    "miR-204-5p forced ON with CASP3 Boolean node forced OFF",
    "DFNA5 and CASP3 Boolean nodes forced ON",
    "TP53 and MIR204 ON with CASP3 OFF and CDKN1A ON",
    "MALAT1 and SIRT1 ON with MIR204 OFF",
    "SIRT1 and DFNA5 OFF with CASP3 ON",
    "BAX Boolean node forced ON",
    "BAX forced ON with DFNA5 OFF",
    "TP53 Boolean node forced ON",
    "TP53 Boolean node forced OFF",
    "CDKN1A Boolean node forced ON"
  ),
  evidence_scope = c(
    "Model reference", "Direct experimental evidence; cancer, non-HCC",
    "Direct experimental evidence in HepG2", "Direct experimental evidence in HCC",
    "Mechanistic inverse inferred from HCC evidence", "Direct experimental evidence in HCC",
    "Mechanistic inverse inferred from HCC evidence", "Axis-supported model perturbation",
    "Axis-supported model perturbation", "Canonical mitochondrial apoptosis mechanism",
    "In silico forced Boolean state; not direct evidence of GSDME cleavage",
    "In silico forced Boolean state; not direct evidence of CASP3 activity",
    "In silico combined intervention; not experimentally observed in GEO",
    rep("In silico dependency control; not experimentally applied in GEO", 13L)
  ),
  reference = c(
    "Model assumption", "PMID:28459430; DOI:10.1038/nature22393",
    "PMID:35747157; DOI:10.3892/etm.2022.11383",
    "PMID:27748572; DOI:10.1002/cbf.3223",
    "PMID:27748572; DOI:10.1002/cbf.3223",
    "PMID:28720061; DOI:10.1177/1010428317718135",
    "PMID:28720061; DOI:10.1177/1010428317718135",
    "PMID:27748572; PMID:28720061", "PMID:27748572; PMID:28720061",
    "Network mechanism; experimental validation required in selected HCC system",
    "Model rules: PYROPTOSIS = GSDME; OE is a node clamp",
    "Model rules: GSDME = CASP3; OE is a node clamp",
    "Model-defined miR-204-5p OE plus GSDME KO; requires a new simulation",
    rep("Model-defined combined Boolean perturbation; no direct GEO perturbation", 13L)
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
    if (any(!names(clamp) %in% model$node_ids))
      stop("Unknown perturbation name: ", paste(setdiff(names(clamp), model$node_ids), collapse = ", "))
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
perturbation_pheno$perturbation_label <- display_perturbation(perturbation_pheno$perturbation)
perturbation_pheno$node_label <- display_node(perturbation_pheno$node)
p_perturb <- plot_heatmap_long(
  perturbation_pheno, "perturbation_label", "node_label", "activation_frequency",
  "In silico perturbation screen",
  "Asynchronous endpoint activation frequencies; literature evidence is reported separately"
)
save_publication_plot(p_perturb, "Figure_02_in_silico_perturbation_heatmap", 12.5, 8)

# Comparacao dedicada das quatro intervencoes pedidas; usa as mesmas
# trajetorias e frequencias calculadas para a Figura 02.
gsdme_casp3_conditions <- c("GSDME OE", "GSDME KO", "CASP3 OE", "CASP3 KO")
gsdme_casp3_nodes <- c("CASP3", "DFNA5", "PYROPTOSIS", "APOPTOSIS", "RESISTANCE")
gsdme_casp3_plot <- perturbation_results[
  perturbation_results$perturbation %in% gsdme_casp3_conditions &
    perturbation_results$node %in% gsdme_casp3_nodes, , drop = FALSE
]
gsdme_casp3_plot$condition_label <- factor(
  display_condition_e1(gsdme_casp3_plot$perturbation),
  levels = rev(display_condition_e1(gsdme_casp3_conditions))
)
gsdme_casp3_plot$node_label <- factor(
  display_node(gsdme_casp3_plot$node),
  levels = display_node(gsdme_casp3_nodes)
)
p_gsdme_casp3 <- ggplot(gsdme_casp3_plot,
                         aes(x = node_label, y = condition_label,
                             fill = activation_frequency)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * activation_frequency)), size = 4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                      labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "GSDME and CASP3: activation versus knockout",
       subtitle = "Boolean node clamps (GSDME = DFNA5); asynchronous simulation endpoints",
       x = NULL, y = NULL, fill = "Endpoints\nON") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 22, hjust = 1),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(12, 18, 12, 18))
save_publication_plot(p_gsdme_casp3, "Figure_02b_GSDME_CASP3_OE_KO_comparison", 11, 5.5)

# DDR ON in every row. These five phenotype outputs are scored independently.
dependency_conditions <- c(
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
dependency_plot <- perturbation_results[
  perturbation_results$perturbation %in% dependency_conditions &
    perturbation_results$node %in% c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE",
                                      "PROLIFERATION", "CELL_CYCLE_ARREST"), , drop = FALSE
]
dependency_plot$condition_label <- factor(
  display_condition_e1(dependency_plot$perturbation),
  levels = rev(display_condition_e1(dependency_conditions))
)
dependency_plot$phenotype_label <- factor(
  display_node(dependency_plot$node),
  levels = display_node(c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE",
                          "PROLIFERATION", "CELL_CYCLE_ARREST"))
)
p_dependency <- ggplot(dependency_plot,
                       aes(x = phenotype_label, y = condition_label,
                           fill = activation_frequency)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * activation_frequency)), size = 3) +
  scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                      labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Perturbation based validation of GSDME dependence in silico",
       subtitle = "DDR ON in all conditions | asynchronous fixed-point endpoints",
       caption = "Outputs may coexist. GSDME E1 clamps DFNA5 ON; it does not measure GSDME cleavage.",
       x = NULL, y = NULL, fill = "Output ON") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))
save_publication_plot(p_dependency, "Figure_02c_GSDME_dependency_DDR_ON", 15, 9)

# Feedback-loop panel: activity of upstream and terminal nodes, keeping
# apoptosis and pyroptosis as separate (possibly coexisting) model outputs.
loop_conditions <- c(
  "Reference: DDR ON, GSDME endogenous", "miR-204-5p OE",
  "miR-204-5p KO", "MALAT1 OE", "MALAT1 KO", "SIRT1 OE",
  "SIRT1 KO", "p53 E1 + miR-204-5p E1 + CASP3 KO + p21 E1",
  "p53 E1", "p53 KO", "p21 E1",
  "MALAT1 E1 + miR-204-5p KO + SIRT1 E1",
  "SIRT1 KO + GSDME KO + CASP3 E1"
)
loop_nodes <- c("MALAT1", "MIR204", "SIRT1", "TP53", "CDKN1A",
                "BAX", "CASP3", "DFNA5", "PYROPTOSIS", "APOPTOSIS")
loop_data <- perturbation_results[
  perturbation_results$perturbation %in% loop_conditions &
    perturbation_results$node %in% loop_nodes, , drop = FALSE
]
loop_data$DDR <- "ON"
write_table(loop_data, "28_feedback_loop_perturbation_outcomes_DDR_ON.csv")
loop_data$condition_label <- factor(
  display_condition_e1(loop_data$perturbation),
  levels = rev(display_condition_e1(loop_conditions))
)
loop_data$node_label <- factor(
  unname(c(MALAT1 = "MALAT1", MIR204 = "miR-204-5p", SIRT1 = "SIRT1",
           TP53 = "p53", CDKN1A = "p21", BAX = "BAX", CASP3 = "CASP3",
           DFNA5 = "GSDME", PYROPTOSIS = "Pyroptosis",
           APOPTOSIS = "Apoptosis")[loop_data$node]),
  levels = c("MALAT1", "miR-204-5p", "SIRT1", "p53", "p21", "BAX",
             "CASP3", "GSDME", "Pyroptosis", "Apoptosis")
)
p_loops <- ggplot(loop_data,
                  aes(x = node_label, y = condition_label,
                      fill = activation_frequency)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * activation_frequency)), size = 2.8) +
  scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1)) +
  labs(title = "Perturbation outcomes of the miR-204-5p/SIRT1/p53 and p53/MALAT1 loops",
       subtitle = "DDR ON | upstream nodes, BAX/CASP3/GSDME axis and death outputs",
       caption = "E1 forces a Boolean node ON; KO forces it OFF. Outputs are model predictions.",
       x = NULL, y = NULL, fill = "Node ON") +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"))
save_publication_plot(p_loops, "Figure_02e_feedback_loop_outcomes_DDR_ON", 17, 8)

# ----------- Exact stable states and phenotype reachability ------------------

# The input DDR is fixed at 1 in every analysis. All conditions below are
# applied to the GINsim 'name' attributes, not to the original node IDs.
reachability_interventions <- list(
  `Unperturbed` = integer(0),
  `CASP3 OE + GSDME KO` = c(CASP3 = 1L, DFNA5 = 0L),
  `miR-204-5p OE + GSDME KO` = c(MIR204 = 1L, DFNA5 = 0L),
  `miR-204-5p OE + CASP3 KO` = c(MIR204 = 1L, CASP3 = 0L),
  `SIRT1 KO + GSDME KO` = c(SIRT1 = 0L, DFNA5 = 0L),
  `MALAT1 KO + GSDME KO` = c(MALAT1 = 0L, DFNA5 = 0L),
  `GSDME OE` = c(DFNA5 = 1L),
  `GSDME KO` = c(DFNA5 = 0L),
  `CASP3 OE` = c(CASP3 = 1L),
  `CASP3 KO` = c(CASP3 = 0L),
  `SIRT1 OE` = c(SIRT1 = 1L),
  `SIRT1 KO` = c(SIRT1 = 0L),
  `MALAT1 OE` = c(MALAT1 = 1L),
  `MALAT1 KO` = c(MALAT1 = 0L),
  `miR-204-5p OE` = c(MIR204 = 1L),
  `miR-204-5p KO` = c(MIR204 = 0L),
  `GSDME E1 + CASP3 E1` = c(DFNA5 = 1L, CASP3 = 1L),
  `p53 E1 + miR-204-5p E1 + CASP3 KO + p21 E1` =
    c(TP53 = 1L, MIR204 = 1L, CASP3 = 0L, CDKN1A = 1L),
  `MALAT1 E1 + miR-204-5p KO + SIRT1 E1` =
    c(MALAT1 = 1L, MIR204 = 0L, SIRT1 = 1L),
  `SIRT1 KO + GSDME KO + CASP3 E1` = c(SIRT1 = 0L, DFNA5 = 0L, CASP3 = 1L),
  `BAX E1` = c(BAX = 1L),
  `BAX E1 + GSDME KO` = c(BAX = 1L, DFNA5 = 0L),
  `p53 E1` = c(TP53 = 1L),
  `p53 KO` = c(TP53 = 0L),
  `p21 E1` = c(CDKN1A = 1L)
)
five_phenotypes <- c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE",
                     "PROLIFERATION", "CELL_CYCLE_ARREST")
if (!all(five_phenotypes %in% MODEL$node_ids))
  stop("The five phenotype outputs are required for stable-state comparisons")

# Removing a feedback vertex set leaves a directed acyclic dependency graph.
# Every assignment to those feedback nodes uniquely determines the remaining
# nodes by topological evaluation. Testing all assignments finds ALL Boolean
# fixed points without enumerating 2^31 full states or sampling attractors.
enumerate_stable_states <- function(model, clamp, max_feedback = 18L) {
  if (!all(model$input_ids %in% names(clamp)))
    stop("Fix every GINsim input before exact stable-state enumeration")
  unknown <- setdiff(names(clamp), model$node_ids)
  if (length(unknown)) stop("Unknown intervention name(s): ", paste(unknown, collapse = ", "))
  remaining <- setdiff(model$node_ids, names(clamp))
  dependencies <- lapply(remaining, function(target) {
    expr <- model$compiled_rules[[target]]
    sources <- if (is.null(expr)) character(0) else
      intersect(all.vars(expr), remaining)
    data.frame(from = sources, to = rep(target, length(sources)),
               stringsAsFactors = FALSE)
  })
  dependency_edges <- if (length(dependencies)) do.call(rbind, dependencies) else
    data.frame(from = character(0), to = character(0))
  graph <- igraph::graph_from_data_frame(
    dependency_edges, directed = TRUE, vertices = data.frame(name = remaining)
  )
  feedback <- character(0)
  acyclic <- graph
  while (!igraph::is_dag(acyclic)) {
    components <- igraph::components(acyclic, mode = "strong")
    cycle_nodes <- igraph::V(acyclic)$name[
      components$csize[components$membership] > 1L
    ]
    edges <- igraph::as_data_frame(acyclic, what = "edges")
    cycle_nodes <- unique(c(cycle_nodes, edges$from[edges$from == edges$to]))
    if (!length(cycle_nodes)) stop("Could not identify a feedback node")
    internal_degree <- vapply(cycle_nodes, function(node) {
      sum(edges$from == node & edges$to %in% cycle_nodes) +
        sum(edges$to == node & edges$from %in% cycle_nodes)
    }, integer(1))
    selected <- cycle_nodes[[which.max(internal_degree)]]
    feedback <- c(feedback, selected)
    acyclic <- igraph::delete_vertices(
      acyclic, which(igraph::V(acyclic)$name == selected)
    )
  }
  if (length(feedback) > max_feedback) {
    return(list(states = list(), complete = FALSE, feedback_nodes = feedback,
                assignments_checked = 0L))
  }
  topological_order <- igraph::as_ids(igraph::topo_sort(acyclic, mode = "out"))
  assignments <- if (length(feedback)) {
    expand.grid(rep(list(0:1), length(feedback)), KEEP.OUT.ATTRS = FALSE)
  } else data.frame(.single = 1L)
  stable <- list()
  for (i in seq_len(nrow(assignments))) {
    state <- apply_clamp(empty_state(model), clamp)
    if (length(feedback))
      state[feedback] <- as.integer(unlist(assignments[i, , drop = FALSE], use.names = FALSE))
    for (node in topological_order) {
      values <- as.list(as.logical(unname(state)))
      names(values) <- names(state)
      state[[node]] <- evaluate_rule(
        model$compiled_rules[[node]], list2env(values, parent = baseenv())
      )
    }
    # Enforce the original equations, including the enumerated feedback nodes.
    target <- logical_targets(state, model, clamp)
    if (identical(unname(state), unname(target))) {
      stable[[length(stable) + 1L]] <- state
    }
  }
  list(states = stable, complete = TRUE, feedback_nodes = feedback,
       assignments_checked = nrow(assignments))
}

# A stochastic asynchronous path may revisit a transient state. Such a
# revisit is NOT proof of a cyclic attractor. Only a state whose every rule
# agrees with its current value is counted as a reached stable state here.
reach_stable_state <- function(model, clamp, max_steps = 350L) {
  state <- random_initial_state(model, clamp)
  for (step in 0:max_steps) {
    target <- logical_targets(state, model, clamp)
    unstable <- setdiff(names(state)[state != target], names(clamp))
    if (!length(unstable))
      return(list(state = state, stable = TRUE, steps = step))
    if (step == max_steps) break
    node <- sample(unstable, 1L)
    state[[node]] <- target[[node]]
  }
  list(state = state, stable = FALSE, steps = max_steps)
}

wilson_interval <- function(success, total, z = stats::qnorm(0.975)) {
  p <- success / total
  denominator <- 1 + z^2 / total
  centre <- (p + z^2 / (2 * total)) / denominator
  margin <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  c(max(0, centre - margin), min(1, centre + margin))
}

stable_state_rows <- list()
stable_summary_rows <- list()
reachability_rows <- list()
row_index <- 1L
summary_index <- 1L
for (ddr_value in 1L) {
  for (condition_name in names(reachability_interventions)) {
    clamp <- c(DDR = ddr_value, reachability_interventions[[condition_name]])
    exact <- enumerate_stable_states(MODEL, clamp)
    ddr_label <- "DDR ON"
    if (!exact$complete)
      log_message("Exact enumeration skipped for ", ddr_label, ", ",
                  condition_name, ": too many feedback nodes")
    for (i in seq_along(exact$states)) {
      stable_state_rows[[length(stable_state_rows) + 1L]] <- data.frame(
        DDR_setting = ddr_label, condition = condition_name,
        stable_state = paste0("SS", i), t(exact$states[[i]]),
        check.names = FALSE
      )
    }
    for (phenotype in five_phenotypes) {
      stable_summary_rows[[summary_index]] <- data.frame(
        DDR = ddr_label, condition = condition_name, phenotype = phenotype,
        exact_enumeration_complete = exact$complete,
        n_stable_states = if (exact$complete) length(exact$states) else NA_integer_,
        states_with_phenotype = if (exact$complete) sum(vapply(
          exact$states, function(state) state[[phenotype]] == 1L, logical(1)
        )) else NA_integer_,
        feedback_nodes = paste(exact$feedback_nodes, collapse = ";"),
        assignments_checked = exact$assignments_checked,
        stringsAsFactors = FALSE
      )
      summary_index <- summary_index + 1L
    }
    reached <- matrix(0L, nrow = n_perturb, ncol = length(five_phenotypes),
                      dimnames = list(NULL, five_phenotypes))
    converged <- logical(n_perturb)
    for (trial in seq_len(n_perturb)) {
      # Pair random initial states within a DDR panel for all conditions.
      set.seed(as.integer(CFG$seed + 100000L * ddr_value + trial))
      result <- reach_stable_state(MODEL, clamp)
      converged[[trial]] <- result$stable
      if (result$stable)
        reached[trial, ] <- result$state[five_phenotypes]
    }
    for (phenotype in five_phenotypes) {
      count <- sum(reached[, phenotype])
      interval <- wilson_interval(count, n_perturb)
      reachability_rows[[row_index]] <- data.frame(
        DDR = ddr_label, condition = condition_name, phenotype = phenotype,
        reached_stable_state_with_phenotype = count,
        n_stable_trajectories = sum(converged), n_trajectories = n_perturb,
        reachability = count / n_perturb,
        ci95_low = interval[[1]], ci95_high = interval[[2]],
        max_steps = 350L, seed = CFG$seed,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
    log_message("DDR=", ddr_value, " | ", condition_name,
                " | exact fixed points: ", if (exact$complete) length(exact$states) else "skipped",
                " | stochastic trajectories reaching a fixed point: ",
                sum(converged), "/", n_perturb)
  }
}
stable_state_table <- if (length(stable_state_rows)) do.call(rbind, stable_state_rows) else
  data.frame(DDR_setting = character(0), condition = character(0), stable_state = character(0))
stable_summary <- do.call(rbind, stable_summary_rows)
phenotype_reachability <- do.call(rbind, reachability_rows)
write_table(stable_state_table, "21_exact_stable_states_DDR_ON.csv")
write_table(stable_summary, "22_exact_stable_state_phenotypes_DDR_ON.csv")
write_table(phenotype_reachability, "23_phenotype_reachability_DDR_ON.csv")

literature_context <- data.frame(
  module = c("CASP3/GSDME", "MALAT1/miR-204/SIRT1", "miR-204-5p/SIRT1"),
  primary_reference = c(
    "Wang et al., Nature (2017)", "Hou et al., Tumour Biology (2017)",
    "Jiang et al., Cell Biochemistry and Function (2016)"
  ),
  doi = c("10.1038/nature22393", "10.1177/1010428317718135",
          "10.1002/cbf.3223"),
  relevance = c(
    "Caspase-3 cleavage of GSDME links apoptosis to pyroptosis in cancer models",
    "MALAT1 sponges miR-204 and releases SIRT1 in HCC models",
    "miR-204-5p targets SIRT1 and affects HCC cell phenotypes"
  ),
  limitation = "Published interactions support components, not the full simulated stable states or DDR ON/OFF predictions",
  stringsAsFactors = FALSE
)
write_table(literature_context, "25_primary_literature_context.csv")

phenotype_labels <- c(
  PYROPTOSIS = "GSDME-dependent pyroptosis", APOPTOSIS = "Apoptosis",
  RESISTANCE = "Resistance", PROLIFERATION = "Proliferation",
  CELL_CYCLE_ARREST = "Cell-cycle arrest"
)
unperturbed_states <- stable_state_table[
  stable_state_table$condition == "Unperturbed", , drop = FALSE
]
if (nrow(unperturbed_states)) {
  state_heatmap <- do.call(rbind, lapply(seq_len(nrow(unperturbed_states)), function(i) {
    data.frame(DDR = unperturbed_states$DDR_setting[[i]],
               stable_state = unperturbed_states$stable_state[[i]],
               phenotype = five_phenotypes,
               value = as.integer(unlist(
                 unperturbed_states[i, five_phenotypes, drop = FALSE],
                 use.names = FALSE
               )),
               stringsAsFactors = FALSE)
  }))
  state_heatmap$phenotype_label <- factor(
    unname(phenotype_labels[state_heatmap$phenotype]),
    levels = rev(unname(phenotype_labels[five_phenotypes]))
  )
  p_exact <- ggplot(state_heatmap,
                    aes(x = stable_state, y = phenotype_label, fill = factor(value))) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = ifelse(value == 1L, "ON", "OFF")), size = 4) +
    facet_wrap(~ DDR, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c(`0` = "#E8EEF2", `1` = "#3A9D84"),
                      labels = c(`0` = "OFF", `1` = "ON")) +
    labs(title = "Exact stable states of the unperturbed Boolean network",
         subtitle = "Complete fixed-point enumeration; DDR is the only clamped input",
         caption = "Columns are distinct fixed points. Phenotype outputs may coexist within one state.",
         x = "Stable state", y = NULL, fill = "Output") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
  save_publication_plot(p_exact, "Figure_08_exact_stable_states_DDR_ON", 11, 6)
}

exact_condition_plot <- stable_summary[stable_summary$exact_enumeration_complete, ]
if (nrow(exact_condition_plot)) {
  exact_condition_plot$stable_fraction <- with(exact_condition_plot,
    ifelse(n_stable_states > 0L, states_with_phenotype / n_stable_states, NA_real_)
  )
  exact_condition_plot$cell_label <- with(exact_condition_plot,
    ifelse(n_stable_states > 0L,
           paste0(states_with_phenotype, "/", n_stable_states), "No fixed points")
  )
  exact_condition_plot$condition_label <- factor(
    display_condition_e1(exact_condition_plot$condition),
    levels = rev(display_condition_e1(names(reachability_interventions)))
  )
  exact_condition_plot$phenotype_label <- factor(
    unname(phenotype_labels[exact_condition_plot$phenotype]),
    levels = unname(phenotype_labels[five_phenotypes])
  )
  p_exact_conditions <- ggplot(exact_condition_plot,
                                aes(x = phenotype_label, y = condition_label,
                                    fill = stable_fraction)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = cell_label), size = 3) +
    facet_wrap(~ DDR, nrow = 1) +
    scale_fill_gradient(low = "#F7FBFF", high = "#3A9D84", limits = c(0, 1),
                        na.value = "grey90",
                        labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "Exact stable-state phenotypes under each intervention",
         subtitle = "Each cell shows fixed points with the output ON / all fixed points under the same DDR input",
         caption = paste0("Every fixed point is counted once. These ratios are NOT probabilities of reaching the states.\n",
                          "GSDME E1 forces the DFNA5 Boolean node ON and is not an assay of protein cleavage."),
         x = NULL, y = NULL, fill = "Fixed points") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 32, hjust = 1),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
  save_publication_plot(p_exact_conditions,
                        "Figure_08b_exact_stable_states_by_condition_DDR_ON", 16, 9)
}

# Direct DDR-ON fixed-point check of the GSDME dependency hypothesis.
# Count fixed points, not stochastic trajectory frequencies: do not interpret
# an ON/total ratio as the probability of entering an attractor.
dependency_exact <- stable_summary[
  stable_summary$DDR == "DDR ON" &
    stable_summary$condition %in% c(
      "Unperturbed", "GSDME OE", "GSDME KO", "CASP3 OE",
      "CASP3 OE + GSDME KO", "CASP3 KO", "miR-204-5p OE",
      "miR-204-5p OE + GSDME KO", "miR-204-5p OE + CASP3 KO",
      "SIRT1 KO", "SIRT1 KO + GSDME KO", "GSDME E1 + CASP3 E1",
      "p53 E1 + miR-204-5p E1 + CASP3 KO + p21 E1",
      "MALAT1 E1 + miR-204-5p KO + SIRT1 E1",
      "SIRT1 KO + GSDME KO + CASP3 E1", "BAX E1",
      "BAX E1 + GSDME KO"
    ), , drop = FALSE
]
write_table(dependency_exact, "27_GSDME_dependency_exact_stable_states_DDR_ON.csv")
dependency_exact <- dependency_exact[
  dependency_exact$exact_enumeration_complete &
    !is.na(dependency_exact$n_stable_states) &
    dependency_exact$n_stable_states > 0L, , drop = FALSE
]
if (nrow(dependency_exact)) {
  dependency_exact$on_fraction <- with(
    dependency_exact, states_with_phenotype / n_stable_states
  )
  dependency_exact$label <- with(
    dependency_exact, paste0(states_with_phenotype, "/", n_stable_states)
  )
  dependency_exact$condition <- factor(
    display_condition_e1(dependency_exact$condition),
    levels = rev(display_condition_e1(unique(as.character(dependency_exact$condition))))
  )
  dependency_exact$phenotype <- factor(
    unname(phenotype_labels[dependency_exact$phenotype]),
    levels = unname(phenotype_labels[five_phenotypes])
  )
  p_dependency_exact <- ggplot(dependency_exact,
                               aes(x = phenotype, y = condition,
                                   fill = on_fraction)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = label), size = 3.3) +
    scale_fill_gradient(low = "#F7FBFF", high = "#3A9D84", limits = c(0, 1)) +
    labs(title = "Exact GSDME-dependent phenotypes at DDR ON",
         subtitle = "Numerators count fixed points with each Boolean output ON",
         caption = paste("Each cell shows ON / all exact fixed points for that condition.",
                         "These fractions are not cell or trajectory probabilities."),
         x = NULL, y = NULL, fill = "Fixed points ON") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold"))
  save_publication_plot(p_dependency_exact,
                        "Figure_08c_GSDME_dependency_exact_DDR_ON", 17, 9)
}

unperturbed_reach <- phenotype_reachability[
  phenotype_reachability$condition == "Unperturbed", , drop = FALSE
]
unperturbed_reach$phenotype_label <- factor(
  unname(phenotype_labels[unperturbed_reach$phenotype]),
  levels = unname(phenotype_labels[five_phenotypes])
)
p_unperturbed <- ggplot(unperturbed_reach,
                        aes(x = phenotype_label, y = reachability, fill = DDR)) +
  geom_col(width = 0.68) +
  geom_errorbar(aes(ymin = ci95_low, ymax = ci95_high), width = 0.20,
                linewidth = 0.5) +
  geom_text(aes(y = pmin(reachability + 0.065, 1.1),
                label = sprintf("%.0f%%", 100 * reachability)), size = 3.4) +
  facet_wrap(~ DDR, nrow = 1) +
  scale_fill_manual(values = c(`DDR ON` = "#0072B2")) +
  scale_y_continuous(limits = c(0, 1.17), breaks = seq(0, 1, 0.25),
                     labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Unperturbed network: all five phenotype outputs",
       subtitle = paste0("Asynchronous trajectories reaching fixed points; ",
                         n_perturb, " random initial states with DDR ON"),
       caption = paste0("Bars show the fraction of all trials reaching a fixed point with each output ON; 95% Wilson intervals.\n",
                        "Outputs may coexist; zero values are shown for all five outputs."),
       x = NULL, y = "Empirical fixed-point reachability", fill = "Input") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 32, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "none", plot.title = element_text(face = "bold"))
save_publication_plot(p_unperturbed, "Figure_09_unperturbed_five_phenotypes", 15, 7)

phenotype_reachability$condition_label <- factor(
  display_condition_e1(phenotype_reachability$condition),
  levels = rev(display_condition_e1(names(reachability_interventions)))
)
phenotype_reachability$phenotype_label <- factor(
  unname(phenotype_labels[phenotype_reachability$phenotype]),
  levels = unname(phenotype_labels[five_phenotypes])
)
p_reach <- ggplot(phenotype_reachability,
                  aes(x = phenotype_label, y = condition_label, fill = reachability)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * reachability)), size = 2.9) +
  facet_wrap(~ DDR, nrow = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                      labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Phenotype reachability with DDR ON",
       subtitle = "Unperturbed and perturbed conditions; asynchronous trajectories reaching verified fixed points",
       caption = paste0("Percentages use all random starts as denominator; non-converged runs are reported in the companion CSV.\n",
                        "Outputs may coexist. GSDME E1 forces the DFNA5 Boolean node ON, not protein expression alone."),
       x = NULL, y = NULL, fill = "Reachability") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 32, hjust = 1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))
save_publication_plot(p_reach, "Figure_10_phenotype_reachability_DDR_ON", 16, 9)

baselines <- phenotype_reachability[
  phenotype_reachability$condition == "Unperturbed",
  c("DDR", "phenotype", "reachability")
]
names(baselines)[names(baselines) == "reachability"] <- "unperturbed_reachability"
perturbed_comparison <- merge(
  phenotype_reachability[phenotype_reachability$condition != "Unperturbed", ],
  baselines, by = c("DDR", "phenotype"), sort = FALSE
)
perturbed_comparison$delta <- perturbed_comparison$reachability -
  perturbed_comparison$unperturbed_reachability
write_table(perturbed_comparison, "24_perturbed_minus_unperturbed_DDR_ON.csv")
perturbed_comparison$condition_label <- factor(
  display_condition_e1(perturbed_comparison$condition),
  levels = rev(display_condition_e1(setdiff(names(reachability_interventions), "Unperturbed")))
)
perturbed_comparison$phenotype_label <- factor(
  unname(phenotype_labels[perturbed_comparison$phenotype]),
  levels = unname(phenotype_labels[five_phenotypes])
)
p_delta <- ggplot(perturbed_comparison,
                  aes(x = phenotype_label, y = condition_label, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%+.0f pp", 100 * delta)), size = 2.8) +
  facet_wrap(~ DDR, nrow = 1) +
  scale_fill_gradient2(low = "#D55E00", mid = "#F7F7F7", high = "#0072B2",
                       midpoint = 0, limits = c(-1, 1),
                       labels = function(x) paste0(round(100 * x), " pp")) +
  labs(title = "Perturbed versus unperturbed phenotype reachability",
       subtitle = "Differences in percentage points with DDR ON; model predictions only",
       caption = paste0(
         "Mechanistic context: Wang et al., Nature 2017 (doi:10.1038/nature22393); ",
         "Hou et al., Tumour Biol 2017 (doi:10.1177/1010428317718135).\n",
         "Jiang et al., Cell Biochem Funct 2016 (doi:10.1002/cbf.3223). ",
         "These papers do not validate the simulated DDR comparisons."
       ),
       x = NULL, y = NULL, fill = "Change") +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 32, hjust = 1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(size = 8.2, hjust = 0))
save_publication_plot(p_delta, "Figure_11_perturbed_vs_unperturbed_DDR_ON", 16, 9)

# ---------------------- Structural control / drivers -------------------------

structural_ranking <- function(model) {
  g <- igraph::graph_from_data_frame(model$edges[, c("from", "to")], directed = TRUE,
                                     vertices = data.frame(name = model$nodes$name))
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
  MALAT1 = 0L, MIR204 = 1L, SIRT1 = 0L,
  BCL2 = 0L, BAX = 1L, TP53 = 1L, CASP3 = 1L
)
candidate_actions <- candidate_actions[names(candidate_actions) %in% MODEL$node_ids]

control_success <- function(model, clamp, start_states, target_node = "PYROPTOSIS",
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
      clamp <- c(BASE_CLAMP, candidate_actions[nodes])
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
top_driver$display_intervention <- vapply(
  top_driver$intervention, display_driver_intervention, character(1)
)
top_driver$display_intervention <- factor(
  top_driver$display_intervention, levels = rev(top_driver$display_intervention)
)
minimum_driver_size <- attr(driver_results, "minimum_size")
driver_subtitle <- if (is.na(minimum_driver_size)) {
  "No intervention reached the 95% threshold in the tested search space"
} else {
  paste0("Minimum driver-set size = ", minimum_driver_size,
         "; dashed line indicates the 95% robustness threshold")
}
p_driver <- ggplot(top_driver, aes(x = success_rate, y = display_intervention)) +
  geom_col(width = 0.70, fill = "#009E73") +
  geom_vline(xintercept = 0.95, linetype = 2, colour = "#B2182B") +
  geom_text(aes(x = success_rate + 0.015,
                label = paste0(round(100 * success_rate), "%")),
            hjust = 0, colour = "#17212B", fontface = "bold", size = 3.2) +
  scale_x_continuous(breaks = seq(0, 1, 0.25),
                     expand = ggplot2::expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Minimum attractor-control interventions",
       subtitle = driver_subtitle,
       x = "Control success rate", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), legend.position = "none",
        plot.margin = margin(12, 20, 12, 12))
save_publication_plot(p_driver, "Figure_03_minimum_driver_nodes", 12, 7.5)

# ---------------- Probabilistic scRNA-seq integration ------------------------

gene_map <- list(
  MALAT1 = "MALAT1", SIRT1 = "SIRT1", MYC = "MYC", E2F1 = "E2F1",
  AKT1 = "AKT1", RB1 = "RB1", BCL2 = "BCL2", BAX = "BAX",
  Mdm2 = "MDM2", CDC25A = "CDC25A", CDKN1A = "CDKN1A",
  BBC3 = "BBC3", PRKAA1 = "PRKAA1", PPM1D = "PPM1D", TP53INP1 = "TP53INP1",
  DFNA5 = c("GSDME", "DFNA5"), CASP3 = "CASP3", CASP9 = "CASP9",
  TP53 = "TP53", MIR204 = "MIR204"
)
# mRNA de marcadores de clivagem/atividade nao equivale ao estado Booleano.
observational_only_nodes <- c("DFNA5", "CASP3", "CASP9", "TP53", "MIR204")

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
      mu[intersect(c("MALAT1", "SIRT1", "BCL2", "E2F1"), genes)] <- 9
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
                   "Cell-cycle arrest", "Other")
  out <- matrix(0, nrow = ncol(node_probability), ncol = length(fate_levels),
                dimnames = list(colnames(node_probability), fate_levels))
  for (cell in seq_len(ncol(node_probability))) {
    fate <- character(n_draws)
    for (draw in seq_len(n_draws)) {
      s <- random_initial_state(model, BASE_CLAMP)
      for (node in setdiff(intersect(rownames(node_probability), model$node_ids),
                           observational_only_nodes)) {
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
  save_publication_plot(p_scrna, "Supplementary_Figure_S2_generic_scRNA_cell_fates", 12, 5.8)

  fate_counts <- as.data.frame(table(fate = dominant_fate), stringsAsFactors = FALSE)
  fate_counts$fraction <- fate_counts$Freq / sum(fate_counts$Freq)
  p_fate <- ggplot(fate_counts, aes(x = reorder(fate, fraction), y = fraction, fill = fate)) +
    geom_col(show.legend = FALSE) + coord_flip() +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
    labs(title = "Dominant simulated cell fate", subtitle = scrna_source,
         x = NULL, y = "Fraction of cells") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))
  save_publication_plot(p_fate, "Supplementary_Figure_S3_generic_scRNA_composition", 8, 5)
} else {
  log_message("Generic scRNA-seq module skipped: --no-demo used without --scrna")
}

# ---------------------- GEO validation: liver cancer ------------------------

# GSE125449: somente amostras de HCC com anotacao maligna.
# GSE189903: apenas HCC de regiao tumoral; se nao houver anotacao maligna,
# registrar explicitamente que a amostra pode conter outros tipos celulares.
# Nao ha perturbacoes controladas em nenhum dos dois datasets.

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
  GSE189903 = list(
    title = "Multiregional single-cell dissection of tumor and immune cells in liver cancer (HCC subset)",
    citation = "Ma et al., Nature Communications 2022; DOI:10.1038/s41467-022-35291-5; GEO GSE189903",
    scope = "HCC tumour-region observational concordance; malignant status annotated separately",
    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189903/suppl/",
    files = c("GSE189903_barcodes.tsv.gz", "GSE189903_genes.tsv.gz",
              "GSE189903_matrix.mtx.gz", "GSE189903_Info.txt.gz")
  )
)

geo_requested_files <- function(accession, platform = "droplet") {
  entry <- geo_registry[[accession]]
  if (is.null(entry)) stop("Unsupported GEO accession: ", accession)
  entry$files
}

gzip_geo_valid <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) return(FALSE)
  gzip <- Sys.which("gzip")
  if (nzchar(gzip)) return(identical(as.integer(suppressWarnings(system2(
    gzip, c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE))), 0L))
  tryCatch({
    con <- gzfile(path, "rb")
    on.exit(close(con))
    repeat if (!length(readBin(con, "raw", n = 1048576L))) break
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
}

download_geo_bundle <- function(accession, cache_dir, platform = "droplet") {
  entry <- geo_registry[[accession]]
  files <- geo_requested_files(accession, platform)
  target_dir <- file.path(cache_dir, accession)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout", 60)
  options(timeout = max(old_timeout, 1800))
  on.exit(options(timeout = old_timeout), add = TRUE)
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    filename <- files[[i]]
    destination <- file.path(target_dir, filename)
    was_cached <- gzip_geo_valid(destination)
    if (!was_cached) {
      partial <- paste0(destination, ".part")
      if (file.exists(destination) && !file.exists(partial))
        if (!file.rename(destination, partial)) stop("Cannot preserve incomplete GEO download: ", filename)
      url <- paste0(entry$base_url, filename)
      curl <- Sys.which("curl")
      for (attempt in seq_len(5L)) {
        log_message("Downloading ", accession, ": ", filename, " attempt ", attempt, "/5")
        if (nzchar(curl)) {
          status <- suppressWarnings(tryCatch(system2(curl, c(
            "--location", "--fail", "--show-error", "--connect-timeout", "30",
            "--speed-time", "120", "--speed-limit", "1024", "--continue-at", "-",
            "--output", shQuote(partial), shQuote(url)), stdout = FALSE, stderr = FALSE),
            error = function(e) 99L))
          if (identical(as.integer(status), 33L)) unlink(partial)
        } else {
          if (file.exists(partial)) unlink(partial)
          status <- suppressWarnings(tryCatch(utils::download.file(
            url, partial, mode = "wb", method = "libcurl", quiet = FALSE),
            error = function(e) 99L))
        }
        if (identical(as.integer(status), 0L) && gzip_geo_valid(partial)) {
          if (file.exists(destination)) unlink(destination)
          if (!file.rename(partial, destination)) stop("Cannot finalize GEO file: ", filename)
          break
        }
        if (identical(as.integer(status), 0L)) unlink(partial)
        if (attempt < 5L) Sys.sleep(min(2^attempt, 15))
      }
      if (!gzip_geo_valid(destination)) stop("Incomplete GEO download: ", filename,
        ". Re-run the script to resume its .part file.")
    }
    rows[[i]] <- data.frame(accession = accession, file = filename,
      local_path = normalizePath(destination, mustWork = FALSE),
      size_bytes = file.info(destination)$size,
      source = if (was_cached) "cache" else "downloaded", stringsAsFactors = FALSE)
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

geo_metadata_column <- function(metadata, options) {
  idx <- match(tolower(options), tolower(names(metadata)), nomatch = 0L)
  idx <- idx[idx > 0L]
  if (!length(idx)) NULL else names(metadata)[idx[[1L]]]
}

align_geo_metadata <- function(metadata, barcodes, accession) {
  candidates <- grep("barcode|cell.?id|^cell$", names(metadata),
                     ignore.case = TRUE, value = TRUE)
  candidates <- unique(c(candidates, head(names(metadata), 2L)))
  matches <- vapply(candidates, function(col) sum(barcodes %in% as.character(metadata[[col]])), integer(1))
  if (!length(matches) || max(matches) != length(barcodes))
    stop(accession, ": metadata/barcode mismatch. Metadata columns: ",
         paste(names(metadata), collapse = ", "))
  col <- candidates[which.max(matches)]
  if (anyDuplicated(as.character(metadata[[col]])))
    stop(accession, ": duplicate metadata barcode in ", col)
  metadata[match(barcodes, as.character(metadata[[col]])), , drop = FALSE]
}

hcc_gse125449 <- c(
  "S16_P10_LCP18", "S02_P01_LCP21", "S10_P05_LCP23",
  "S07_P02_LCP28", "S12_P07_LCP30", "S21_P13_LCP37",
  "S15_P09_LCP38", "S351_P10_LCP34", "S364_P21_LCP65"
)

read_gse125449 <- function(manifest) {
  experiments <- list()
  for (set_name in c("Set1", "Set2")) {
    genes <- read_gz_table(geo_file(manifest, paste0(set_name, "_genes\\.tsv")),
                           header = FALSE)
    barcodes <- read_gz_vector(geo_file(manifest, paste0(set_name, "_barcodes\\.tsv")))
    metadata <- read_gz_table(geo_file(manifest, paste0(set_name, "_samples\\.txt")))
    metadata <- align_geo_metadata(metadata, barcodes, "GSE125449 ")
    mat <- read_sparse_mtx_gz(geo_file(manifest, paste0(set_name, "_matrix\\.mtx")))
    if (ncol(mat) != length(barcodes)) stop("GSE125449 ", set_name, ": matrix/barcode mismatch")
    sample_col <- geo_metadata_column(metadata, c("Sample"))
    type_col <- geo_metadata_column(metadata, c("Type"))
    if (is.null(sample_col) || is.null(type_col))
      stop("GSE125449 ", set_name, ": Sample/Type annotation not available")
    sample_id <- as.character(metadata[[sample_col]])
    # Admitir ambas as grafias encontradas para o identificador LCP65.
    hcc <- sample_id %in% c(hcc_gse125449, "364_P21_LCP65")
    malignant <- grepl("malignan|cancer cell|tumou?r cell", metadata[[type_col]], ignore.case = TRUE)
    keep <- which(!is.na(hcc) & hcc & !is.na(malignant) & malignant)
    if (!length(keep)) {
      log_message("GSE125449 ", set_name, ": no HCC malignant cells; this set is skipped")
      next
    }
    cell_ids <- make.unique(paste(sample_id[keep], barcodes[keep], sep = "::"))
    metadata <- metadata[keep, , drop = FALSE]
    metadata$cell_id <- cell_ids
    metadata$sample_id <- sample_id[keep]
    metadata$cell_type <- as.character(metadata[[type_col]])
    metadata$cohort <- set_name
    metadata$annotation_scope <- "HCC; malignant cell annotated"
    rownames(metadata) <- cell_ids
    colnames(mat) <- barcodes
    symbols <- if (ncol(genes) >= 2L) genes[[2]] else genes[[1]]
    target <- extract_network_gene_counts(mat[, keep, drop = FALSE], symbols, gene_map)
    colnames(target) <- cell_ids
    experiments[[paste0("GSE125449_", set_name)]] <- list(
      accession = "GSE125449", experiment = paste0("GSE125449_", set_name),
      platform = "10x", scope = geo_registry$GSE125449$scope,
      counts = target, metadata = metadata)
    rm(mat); invisible(gc())
  }
  if (!length(experiments)) stop("No HCC malignant cells found in GSE125449")
  experiments
}

read_gse189903 <- function(manifest) {
  genes <- read_gz_table(geo_file(manifest, "_genes\\.tsv"), header = FALSE)
  barcodes <- read_gz_vector(geo_file(manifest, "_barcodes\\.tsv"))
  metadata <- read_gz_table(geo_file(manifest, "_Info\\.txt"))
  metadata <- align_geo_metadata(metadata, barcodes, "GSE189903")
  sample_col <- geo_metadata_column(metadata, c("Sample", "sample_id", "sample", "Orig.ident", "patient_region"))
  sample_id <- if (is.null(sample_col)) rep("", length(barcodes)) else as.character(metadata[[sample_col]])
  hist_col <- geo_metadata_column(metadata, c("Cancer type", "cancer_type", "Histology", "histology"))
  region_col <- geo_metadata_column(metadata, c("Tissue type", "tissue_type", "Region", "region"))
  type_col <- geo_metadata_column(metadata, c("Cell type", "cell_type", "CellType", "celltype", "Type", "annotation"))
  types <- if (is.null(type_col)) rep("unannotated", length(barcodes)) else as.character(metadata[[type_col]])
  tag <- grepl("[1-4]H(T[0-9]+|B[0-9]*)", paste(sample_id, barcodes), ignore.case = TRUE)
  if (is.null(hist_col) && !any(tag))
    stop("GSE189903: cannot confirm HCC by sample code or histology")
  hcc <- if (is.null(hist_col)) tag else
    grepl("hepatocellular|^HCC$", as.character(metadata[[hist_col]]), ignore.case = TRUE)
  if (is.null(region_col) && !any(tag)) stop("GSE189903: tumour region not identifiable")
  tumour <- if (is.null(region_col)) tag else {
    region <- as.character(metadata[[region_col]])
    grepl("tumou?r|border|^T[123]$|^B$", region, ignore.case = TRUE) &
      !grepl("non.?tumou?r|adjacent|normal|healthy", region, ignore.case = TRUE)
  }
  selected <- !is.na(hcc) & hcc & !is.na(tumour) & tumour
  malignant <- grepl("malignan|cancer cell|tumou?r cell", types, ignore.case = TRUE)
  malignant[is.na(malignant)] <- FALSE
  if (any(selected & malignant)) selected <- selected & malignant
  keep <- which(selected)
  if (!length(keep)) stop("GSE189903: no HCC tumour-region cells after filtering")
  if (is.null(sample_col)) {
    sample_id <- sub(".*([1-4]H(T[0-9]+|B[0-9]*)).*", "\\1", barcodes)
    if (any(!grepl("^[1-4]H(T[0-9]+|B[0-9]*)$", sample_id[keep], ignore.case = TRUE)))
      stop("GSE189903: selected cells lack verifiable sample origin")
  }
  if (anyNA(sample_id[keep]) || any(!nzchar(trimws(sample_id[keep]))))
    stop("GSE189903: selected HCC cells have no verifiable sample identifier")
  mat <- read_sparse_mtx_gz(geo_file(manifest, "_matrix\\.mtx"))
  if (ncol(mat) != length(barcodes)) stop("GSE189903: matrix/barcode mismatch")
  cell_ids <- make.unique(barcodes[keep])
  colnames(mat) <- barcodes
  metadata <- metadata[keep, , drop = FALSE]
  metadata$cell_id <- cell_ids
  metadata$sample_id <- sample_id[keep]
  metadata$cell_type <- types[keep]
  metadata$cohort <- "HCC tumour region"
  metadata$annotation_scope <- if (any(selected & malignant))
    "HCC; malignant cell annotated" else "HCC; mixed tumour tissue; malignancy not confirmed"
  rownames(metadata) <- cell_ids
  symbols <- if (ncol(genes) >= 2L) genes[[2]] else genes[[1]]
  target <- extract_network_gene_counts(mat[, keep, drop = FALSE], symbols, gene_map)
  colnames(target) <- cell_ids
  rm(mat); invisible(gc())
  list(GSE189903_HCC = list(
    accession = "GSE189903", experiment = "GSE189903_HCC", platform = "10x",
    scope = geo_registry$GSE189903$scope, counts = target, metadata = metadata))
}

read_geo_experiments <- function(accession, manifest, platform) {
  if (accession == "GSE125449") return(read_gse125449(manifest))
  if (accession == "GSE189903") return(read_gse189903(manifest))
  stop("Unsupported GEO series: ", accession)
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
  # Ambos os leitores ja restringiram o material a HCC antes desta etapa.
  group_columns <- c("sample_id")
  max_cells <- if (CFG$quick) min(max_cells, 500L) else max_cells
  keep <- stratified_cell_subset(metadata, max_cells, seed, group_columns)
  lib <- attr(counts, "library_size", exact = TRUE)[keep]
  counts <- counts[, keep, drop = FALSE]
  attr(counts, "library_size") <- lib
  metadata <- metadata[keep, , drop = FALSE]
  metadata$validation_group <- metadata$sample_id
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

# The five Boolean outputs are independent: one fixed point can activate more
# than one phenotype. Pair each RNA-informed starting state across conditions.
summarise_geo_fixed_phenotypes <- function(experiment) {
  outputs <- c("PYROPTOSIS", "APOPTOSIS", "RESISTANCE",
               "PROLIFERATION", "CELL_CYCLE_ARREST")
  conditions <- list(
    `Unperturbed` = c(DDR = 1L),
    `miR-204-5p OE + GSDME KO` = c(DDR = 1L, MIR204 = 1L, DFNA5 = 0L)
  )
  n_draws <- if (CFG$quick) 1L else 3L
  n_cells <- ncol(experiment$binary$node_probability)
  counts <- array(0L, dim = c(length(conditions), n_cells, length(outputs)),
                  dimnames = list(names(conditions), NULL, outputs))
  converged <- matrix(0L, nrow = length(conditions), ncol = n_cells,
                      dimnames = list(names(conditions), NULL))
  walk <- function(initial, clamp) {
    state <- apply_clamp(initial, clamp)
    for (step in 0:350) {
      target <- logical_targets(state, MODEL, clamp)
      unstable <- setdiff(names(state)[state != target], names(clamp))
      if (!length(unstable)) return(state)
      if (step == 350L) break
      node <- sample(unstable, 1L)
      state[[node]] <- target[[node]]
    }
    NULL
  }
  for (cell in seq_len(n_cells)) {
    for (draw in seq_len(n_draws)) {
      matched_seed <- as.integer((as.double(CFG$seed) * 100003 +
        cell * 101 + draw * 17) %% (.Machine$integer.max - 1)) + 1L
      set.seed(matched_seed)
      initial <- random_initial_state(MODEL, c(DDR = 1L))
      for (node in setdiff(intersect(rownames(experiment$binary$node_probability),
                                     MODEL$node_ids), observational_only_nodes)) {
        p <- experiment$binary$node_probability[node, cell]
        if (is.finite(p)) initial[[node]] <- stats::rbinom(1L, 1L, p)
      }
      for (condition in names(conditions)) {
        set.seed(matched_seed + 500000L)
        endpoint <- walk(initial, conditions[[condition]])
        if (!is.null(endpoint)) {
          converged[condition, cell] <- converged[condition, cell] + 1L
          counts[condition, cell, ] <- counts[condition, cell, ] + endpoint[outputs]
        }
      }
    }
  }
  do.call(rbind, lapply(unique(experiment$metadata$validation_group), function(group) {
    idx <- which(experiment$metadata$validation_group == group)
    do.call(rbind, lapply(names(conditions), function(condition) {
      data.frame(
        accession = experiment$accession, experiment = experiment$experiment,
        validation_group = group, validation_scope = experiment$scope,
        condition = condition, phenotype = outputs,
        mean_probability = vapply(outputs, function(node)
          sum(counts[condition, idx, node]), integer(1)) /
          (length(idx) * n_draws),
        n_stable_trajectories = sum(converged[condition, idx]),
        n_trajectories = length(idx) * n_draws,
        n_cells = length(idx), stringsAsFactors = FALSE
      )
    }))
  }))
}

perturbation_node_weights <- function(perturbation, compared_nodes, perturbations, model) {
  clamp <- perturbations[[perturbation]]
  # DDR=1 e contexto basal compartilhado; nao e intervencao.
  direct_nodes <- setdiff(names(clamp), "DDR")
  direct_nodes <- unique(direct_nodes)
  downstream_nodes <- unique(model$edges$to[model$edges$from %in% direct_nodes])
  weights <- stats::setNames(rep(1, length(compared_nodes)), compared_nodes)
  weights[compared_nodes %in% downstream_nodes] <- 2
  weights[compared_nodes %in% direct_nodes] <- 4
  list(weights = weights, direct_nodes = direct_nodes,
       downstream_nodes = downstream_nodes)
}

summarise_geo_perturbation_concordance <- function(experiment, perturbation_results,
                                                   perturbations, model) {
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
      weight_info <- perturbation_node_weights(perturbation, nodes, perturbations, model)
      weighted_concordance <- if (length(nodes)) {
        1 - stats::weighted.mean(abs(observed[nodes] - exp_state[nodes]),
                                 weight_info$weights)
      } else NA_real_
      direct_observed <- intersect(weight_info$direct_nodes, nodes)
      downstream_observed <- intersect(weight_info$downstream_nodes, nodes)
      rank_correlation <- if (length(nodes) >= 3L && stats::sd(observed[nodes]) > 0 &&
                              stats::sd(exp_state[nodes]) > 0) {
        suppressWarnings(stats::cor(observed[nodes], exp_state[nodes], method = "spearman"))
      } else NA_real_
      data.frame(
        accession = experiment$accession, experiment = experiment$experiment,
        validation_group = group, validation_scope = experiment$scope,
        perturbation = perturbation, concordance_score = concordance,
        weighted_concordance_score = weighted_concordance,
        spearman_rank_correlation = rank_correlation,
        nodes_compared = length(nodes), node_coverage = length(nodes) / length(gene_map),
        direct_target_observed = length(direct_observed) > 0,
        direct_nodes = paste(weight_info$direct_nodes, collapse = ";"),
        observed_direct_nodes = paste(direct_observed, collapse = ";"),
        observed_first_order_nodes = paste(downstream_observed, collapse = ";"),
        compared_nodes = paste(nodes, collapse = ";"), n_cells = length(idx),
        evidence_grade = if (!length(direct_observed)) {
          "Indirect only: perturbed node not measured"
        } else if (length(nodes) / length(gene_map) < 0.5) {
          "Low coverage observational concordance"
        } else if (is.finite(weighted_concordance) && weighted_concordance >= 0.75) {
          "Moderate observational concordance"
        } else {
          "Limited observational concordance"
        },
        claim_level = paste(experiment$metadata$annotation_scope[[1]],
          "observational concordance; not causal perturbation validation"),
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
    permitted_claim = paste(unique(x$metadata$annotation_scope),
                            "observational concordance with simulated endpoints"),
    stringsAsFactors = FALSE
  )))
  geo_coverage <- do.call(rbind, lapply(geo_experiments, summarise_geo_coverage))
  geo_fates <- do.call(rbind, lapply(geo_experiments, summarise_geo_fates))
  geo_fixed_phenotypes <- do.call(rbind, lapply(
    geo_experiments, summarise_geo_fixed_phenotypes
  ))
  geo_concordance <- do.call(rbind, lapply(
    geo_experiments, summarise_geo_perturbation_concordance,
    perturbation_results = perturbation_results,
    perturbations = scientific_perturbations, model = MODEL
  ))
  edge_tables <- lapply(geo_experiments, summarise_geo_edge_support, model = MODEL)
  edge_tables <- edge_tables[vapply(edge_tables, nrow, integer(1)) > 0]
  geo_edges <- if (length(edge_tables)) do.call(rbind, edge_tables) else data.frame()

  write_table(geo_scope, "14_geo_dataset_scope_and_claims.csv")
  write_table(geo_coverage, "15_geo_node_coverage.csv")
  write_table(geo_concordance, "16_geo_perturbation_concordance.csv")
  write_table(geo_fates, "17_geo_group_fate_probabilities.csv")
  write_table(geo_fixed_phenotypes,
              "26_GEO_fixed_perturbation_independent_phenotypes.csv")
  if (nrow(geo_edges)) write_table(geo_edges, "18_geo_regulatory_edge_support.csv")

  consensus <- stats::aggregate(
    cbind(concordance_score, weighted_concordance_score, node_coverage) ~ accession + perturbation,
    data = geo_concordance, FUN = mean, na.rm = TRUE
  )
  names(consensus)[names(consensus) == "concordance_score"] <- "mean_concordance"
  names(consensus)[names(consensus) == "weighted_concordance_score"] <-
    "mean_weighted_concordance"
  names(consensus)[names(consensus) == "node_coverage"] <- "mean_node_coverage"
  consensus$validation_role <- ifelse(
    consensus$accession == "GSE125449", "annotated HCC malignant cells",
    "HCC tumour-region cells; inspect malignancy annotation"
  )
  consensus$causal_confirmation <- FALSE
  write_table(consensus, "19_geo_cross_dataset_summary.csv")

  coverage_plot_data <- geo_coverage[geo_coverage$measured, ]
  coverage_plot_data$node_label <- display_node(coverage_plot_data$node)
  coverage_plot_data$experiment_label <- gsub("_", " ", coverage_plot_data$experiment)
  p_geo_coverage <- ggplot(coverage_plot_data,
                           aes(x = node_label, y = experiment_label, fill = detection_rate)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * detection_rate)), size = 3.0) +
    scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                        labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "Single-cell detection of GSDME-network components",
         subtitle = "GSE125449: HCC malignant cells; GSE189903: HCC tumour-region cells",
         x = NULL, y = NULL, fill = "Detected") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"),
          plot.margin = margin(12, 18, 12, 18))
  save_publication_plot(p_geo_coverage, "Figure_04_GEO_node_detection_heatmap", 13, 5.8)

  concordance_plot_data <- stats::aggregate(
    weighted_concordance_score ~ experiment + perturbation, data = geo_concordance,
    FUN = mean, na.rm = TRUE
  )
  concordance_plot_data$experiment_label <- gsub("_", " ", concordance_plot_data$experiment)
  concordance_plot_data$perturbation_label <- display_perturbation(
    concordance_plot_data$perturbation
  )
  p_geo_concordance <- ggplot(concordance_plot_data,
                              aes(x = perturbation_label, y = experiment_label,
                                  fill = weighted_concordance_score)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", weighted_concordance_score)), size = 3.0) +
    scale_fill_gradient2(low = "#D55E00", mid = "#F7F7F7", high = "#0072B2",
                         midpoint = 0.5, limits = c(0, 1)) +
    labs(title = "Weighted single-cell concordance with simulated endpoints",
         subtitle = "Directly perturbed nodes receive 4x weight; observational support, not causal validation",
         x = NULL, y = NULL, fill = "Weighted\nconcordance") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"),
          plot.margin = margin(12, 18, 12, 18))
  save_publication_plot(p_geo_concordance, "Figure_05_GEO_weighted_perturbation_concordance", 14, 6.2)

  geo_fixed_phenotypes$phenotype_label <- factor(
    unname(c(PYROPTOSIS = "Pyroptosis", APOPTOSIS = "Apoptosis",
             RESISTANCE = "Resistance", PROLIFERATION = "Proliferation",
             CELL_CYCLE_ARREST = "Cell-cycle arrest")[geo_fixed_phenotypes$phenotype]),
    levels = c("Pyroptosis", "Apoptosis", "Resistance",
               "Proliferation", "Cell-cycle arrest")
  )
  geo_fixed_phenotypes$sample_label <- paste(
    gsub("_", " ", geo_fixed_phenotypes$experiment),
    geo_fixed_phenotypes$validation_group, sep = " / "
  )
  geo_fixed_phenotypes$condition <- factor(
    geo_fixed_phenotypes$condition,
    levels = c("Unperturbed", "miR-204-5p OE + GSDME KO")
  )
  p_geo_fates <- ggplot(geo_fixed_phenotypes,
                        aes(x = phenotype_label, y = sample_label,
                            fill = mean_probability)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * mean_probability)), size = 2.9) +
    facet_wrap(~ condition, nrow = 1) +
    scale_fill_gradient(low = "#F7FBFF", high = "#0072B2", limits = c(0, 1),
                        labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "Independent Boolean phenotypes under a fixed double perturbation",
         subtitle = "DDR ON; GSE125449 malignant HCC cells and GSE189903 HCC tumour-region cells",
         caption = paste0(
           "Each value is the fraction of matched asynchronous trials reaching a fixed point with the output ON. ",
           "Outputs may coexist and need not sum to 100%.\n",
           "miR-204-5p OE + GSDME KO is simulated, not an experimental GEO condition; ",
           "non-converged trials remain in the denominator (Table 26)."
         ),
         x = NULL, y = NULL, fill = "Reachability") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 32, hjust = 1),
          plot.title = element_text(face = "bold"),
          plot.margin = margin(12, 18, 12, 18))
  save_publication_plot(p_geo_fates, "Figure_06_GEO_projected_fate_composition", 18, 10)

  priority_names <- c("MIR204", "SIRT1", "MALAT1", "TP53", "BAX", "DFNA5", "CASP3")
  priority_detection <- geo_coverage[geo_coverage$node %in% priority_names, ]
  priority_detection$node_label <- factor(
    unname(c(MIR204 = "miR-204-5p (MIR204)", SIRT1 = "SIRT1",
             MALAT1 = "MALAT1", TP53 = "p53 (TP53)", BAX = "BAX",
             DFNA5 = "GSDME (DFNA5)", CASP3 = "CASP3")[priority_detection$node]),
    levels = c("miR-204-5p (MIR204)", "SIRT1", "MALAT1", "BAX",
               "p53 (TP53)", "GSDME (DFNA5)", "CASP3")
  )
  priority_detection$display_rate <- ifelse(
    priority_detection$measured, priority_detection$detection_rate, NA_real_
  )
  priority_detection$label <- ifelse(
    priority_detection$measured,
    sprintf("%.0f%%", 100 * priority_detection$detection_rate), "Not measured"
  )
  p_priority <- ggplot(priority_detection,
                       aes(x = node_label, y = experiment, fill = display_rate)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = label), size = 3.0) +
    scale_fill_gradient(low = "#F7FBFF", high = "#3A9D84", limits = c(0, 1),
                        na.value = "grey88",
                        labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "Single-cell detection of the seven priority network components",
         subtitle = "Detection of mapped RNA features in selected HCC specimens",
         caption = "RNA detection does not establish p53 or CASP3 activity, miRNA function, or GSDME cleavage.",
         x = NULL, y = NULL, fill = "RNA detected") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(face = "bold"))
  save_publication_plot(p_priority, "Figure_06b_seven_priority_nodes_RNA_detection", 15, 5.5)

  if (nrow(geo_edges)) {
    evaluable_edges <- geo_edges[is.finite(geo_edges$spearman_rho), , drop = FALSE]
    if (nrow(evaluable_edges) >= 3L) {
      evaluable_edges$edge_label <- paste0(display_node(evaluable_edges$from),
                                            ifelse(evaluable_edges$expected_sign < 0, " -| ", " -> "),
                                            display_node(evaluable_edges$to))
      evaluable_edges$experiment_label <- gsub("_", " ", evaluable_edges$experiment)
      p_geo_edges <- ggplot(evaluable_edges,
                            aes(x = spearman_rho, y = edge_label,
                                colour = sign_concordant)) +
        geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
        geom_point(size = 3.2) + facet_wrap(~ experiment_label) +
        scale_x_continuous(limits = c(-1, 1)) +
        scale_colour_manual(values = c(`TRUE` = "#007A3D", `FALSE` = "#D55E00"),
                            na.value = "grey60") +
        labs(title = "Evaluable transcriptional correlations for regulatory edges",
             subtitle = "Exploratory support only; protein activity and causal direction are not measured",
             x = "Spearman correlation", y = NULL, colour = "Expected sign") +
        theme_minimal(base_size = 10) +
        theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"),
              legend.position = "bottom")
      save_publication_plot(
        p_geo_edges, "Supplementary_Figure_S1_GEO_regulatory_edge_support", 12, 7
      )
    } else {
      log_message("Supplementary edge-correlation figure skipped: fewer than three evaluable edges; table 18 retained")
    }
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
  nodes <- intersect(c("DFNA5", "MALAT1", "MIR204",
                       "SIRT1", "PRKAA1", "BCL2", "BAX", "TP53",
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
  grid$expression_z[resistant & grid$node %in% c("MALAT1", "SIRT1", "BCL2")] <-
    grid$expression_z[resistant & grid$node %in% c("MALAT1", "SIRT1", "BCL2")] + 1.5
  grid$expression_z[resistant & grid$node %in% c("DFNA5", "MIR204")] <-
    grid$expression_z[resistant & grid$node %in% c("DFNA5", "MIR204")] - 1.2
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
                   "Cell-cycle arrest", "Other")
  out <- matrix(0, nrow = length(patients), ncol = length(fate_levels),
                dimnames = list(patients, fate_levels))
  for (patient in patients) {
    d <- multiomics[multiomics$patient_id == patient, ]
    priors <- stats::setNames(d$activation_prior, d$node)
    priors <- priors[names(priors) %in% model$node_ids]
    fates <- character(n_rep)
    for (r in seq_len(n_rep)) {
      replicate_clamp <- BASE_CLAMP
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
  save_publication_plot(p_twins, "Supplementary_Figure_S4_multiomic_digital_twins", 10.5, 8)
}

# ---------------- Q-learning for sequential therapy --------------------------

rl_actions <- list(
  NONE = integer(0),
  MALAT1_KO = c(MALAT1 = 0L),
  miR204_OE = c(MIR204 = 1L),
  SIRT1_KO = c(SIRT1 = 0L),
  BCL2_KO = c(BCL2 = 0L)
)
rl_actions <- rl_actions[vapply(rl_actions, function(x) all(names(x) %in% MODEL$node_ids), logical(1))]

rl_reward <- function(state, action_name) {
  value <- 0
  # Reward shaping is specific to GSDME-mediated pyroptosis. Intermediate
  # pathway activation receives smaller rewards, while apoptosis without
  # pyroptosis is not considered a terminal therapeutic success.
  if (state[["PYROPTOSIS"]] == 1L) value <- value + 25
  if (state[["DFNA5"]] == 1L) value <- value + 6
  if (state[["CASP3"]] == 1L) value <- value + 3
  if (state[["BAX"]] == 1L) value <- value + 1
  if (state[["APOPTOSIS"]] == 1L && state[["PYROPTOSIS"]] == 0L) value <- value - 1
  if (state[["RESISTANCE"]] == 1L) value <- value - 8
  if (state[["PROLIFERATION"]] == 1L) value <- value - 5
  value <- value - ifelse(action_name == "NONE", 0.01, 0.20)
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
      transient_clamp <- c(BASE_CLAMP, actions[[action_name]])
      transient_clamp <- transient_clamp[!duplicated(names(transient_clamp), fromLast = TRUE)]
      next_state <- step_boolean(state, model, transient_clamp, "synchronous")
      reward <- rl_reward(next_state, action_name)
      terminal <- next_state[["PYROPTOSIS"]] == 1L
      next_key <- q_key(next_state, model, min(step + 1L, horizon))
      next_q <- if (step == horizon || terminal) rep(0, length(actions)) else get_q(next_key)
      q_values[[action_index]] <- q_values[[action_index]] +
        alpha * (reward + gamma * max(next_q) - q_values[[action_index]])
      assign(key, q_values, q)
      state <- next_state
      if (terminal) break
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
    clamp <- c(BASE_CLAMP, agent$actions[[action_name]])
    clamp <- clamp[!duplicated(names(clamp), fromLast = TRUE)]
    next_state <- step_boolean(state, model, clamp, "synchronous")
    rows[[step]] <- data.frame(
      step = step, action = action_name, reward = rl_reward(next_state, action_name),
      fate = classify_fate(next_state),
      pyroptosis_target_reached = next_state[["PYROPTOSIS"]] == 1L,
      t(next_state), check.names = FALSE
    )
    state <- next_state
    if (state[["PYROPTOSIS"]] == 1L) break
  }
  do.call(rbind, rows)
}

agent <- train_q_learning(
  MODEL, MALIGNANT_STATE, rl_actions,
  episodes = if (CFG$quick) 1500L else 12000L,
  horizon = if (CFG$quick) 12L else 18L
)
rl_sequence <- rollout_policy(agent, MODEL, MALIGNANT_STATE)
write_table(rl_sequence, "11_reinforcement_learning_sequence.csv")
rl_target_reached <- any(rl_sequence$pyroptosis_target_reached)
write_table(data.frame(
  target = "GSDME-mediated pyroptosis", reached = rl_target_reached,
  steps_used = nrow(rl_sequence), horizon = agent$horizon,
  final_fate = tail(rl_sequence$fate, 1), stringsAsFactors = FALSE
), "20_RL_pyroptosis_target_summary.csv")

rl_pathway_order <- c(
  "MALAT1", "MIR204", "SIRT1", "PRKAA1", "TP53",
  "BBC3", "BAX", "CASP3", "DFNA5", "APOPTOSIS",
  "PYROPTOSIS", "CELL_CYCLE_ARREST", "RESISTANCE",
  "PROLIFERATION"
)
rl_nodes <- intersect(rl_pathway_order, names(rl_sequence))
active_rl_nodes <- rl_nodes[vapply(rl_nodes, function(node) {
  any(as.numeric(rl_sequence[[node]]) == 1L)
}, logical(1))]
rl_nodes <- unique(c(active_rl_nodes, intersect("PYROPTOSIS", rl_nodes)))
rl_long <- do.call(rbind, lapply(seq_len(nrow(rl_sequence)), function(i) {
  data.frame(step = rl_sequence$step[[i]], action = rl_sequence$action[[i]],
             node = rl_nodes, state = as.numeric(rl_sequence[i, rl_nodes]),
             stringsAsFactors = FALSE)
}))
rl_long$node_label <- display_node(rl_long$node)
rl_long$node_label <- factor(
  rl_long$node_label, levels = rev(display_node(rl_nodes))
)
rl_long$state_label <- factor(rl_long$state, levels = c(0, 1), labels = c("OFF", "ON"))
rl_step_labels <- stats::setNames(
  paste0("Step ", rl_sequence$step, "\n", display_action(rl_sequence$action)),
  rl_sequence$step
)
p_rl <- ggplot(rl_long, aes(x = step, y = node_label, fill = state_label)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  scale_fill_manual(values = c(OFF = "#F7FBFF", ON = "#0072B2"), drop = FALSE) +
  scale_x_continuous(
    breaks = unique(rl_long$step),
    labels = rl_step_labels
  ) +
  labs(title = "Pyroptosis-oriented sequential intervention learned by Q-learning",
       subtitle = if (rl_target_reached) {
         "The trajectory terminates only after GSDME-mediated pyroptosis is reached"
       } else {
         "Pyroptosis target was not reached within the tested horizon; no success claim is made"
       },
       x = "Simulation step", y = NULL, fill = "State") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        plot.margin = margin(12, 18, 12, 18), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
save_publication_plot(p_rl, "Figure_07_pyroptosis_oriented_RL_therapy", 13, 7.5)

# ----------------------------- Manifest --------------------------------------

capture.output(sessionInfo(), file = file.path(LOG_DIR, "sessionInfo.txt"))
manifest <- data.frame(
  file = list.files(CFG$outdir, recursive = TRUE, full.names = FALSE),
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
  pipeline_version = PIPELINE_VERSION,
  seed = CFG$seed,
  quick = CFG$quick,
  geo_max_cells_requested = CFG$geo_max_cells,
  stringsAsFactors = FALSE
)
write_table(manifest, "12_output_manifest.csv")

log_message("Pipeline completed. Results: ", normalizePath(CFG$outdir, mustWork = FALSE))
cat("\nCompleted successfully.\nResults: ", normalizePath(CFG$outdir, mustWork = FALSE), "\n", sep = "")
