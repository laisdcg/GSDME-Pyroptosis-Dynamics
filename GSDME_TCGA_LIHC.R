# GSDME project: independent TCGA-LIHC expression module, version 1.0
# Does not source, edit, or rerun the Boolean simulation pipeline.
# Public GDC harmonized STAR counts + BCGSC miRBase21 isomiRs.

.tcga_packages <- function() {
  p <- c('httr', 'jsonlite', 'xml2', 'ggplot2')
  missing <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Instale os pacotes: install.packages(c(',
    paste(sprintf('"%s"', missing), collapse = ', '), '))', call. = FALSE)
}
.tcga_one <- function(x, default = '') if (length(x)) as.character(x[[1]]) else default
.tcga_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = '')
.tcga_bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
.tcga_api <- function(endpoint, query = list()) {
  last <- NULL
  for (attempt in 1:4) {
    ans <- tryCatch({
      r <- httr::GET(paste0('https://api.gdc.cancer.gov/', endpoint),
                     query = query, httr::accept_json(), httr::timeout(180))
      httr::stop_for_status(r)
      jsonlite::fromJSON(httr::content(r, 'text', encoding = 'UTF-8'), simplifyVector = FALSE)
    }, error = function(e) { last <<- conditionMessage(e); NULL })
    if (!is.null(ans)) return(ans)
    if (attempt < 4) Sys.sleep(2^attempt)
  }
  stop('Falha na API GDC: ', last)
}

tcga_config <- function(model, out = 'resultados_TCGA_LIHC', aliases = NULL) {
  .tcga_packages()
  if (!file.exists(model)) stop('Modelo nao encontrado: ', model)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, mustWork = TRUE)
  for (d in c('cache', 'tables', 'figures')) dir.create(file.path(out, d), showWarnings = FALSE)
  cfg <- list(model = normalizePath(model), out = out, aliases = aliases,
              version = '1.0', project = 'TCGA-LIHC', timeout = 1200,
              attempts = 4L, minimum_pairs = 5L, minimum_correlation_n = 10L)
  saveRDS(cfg, file.path(out, 'config.rds'))
  cfg
}

.tcga_complex_components <- function(display_name) {
  # Extracts gene-symbol-like tokens from a composite Boolean node label,
  # e.g. "(CCND1 ou CCND2 ou CCND3) E (CDK4 ou CDK6)" -> CCND1, CCND2, CCND3, CDK4, CDK6.
  # Generic and language-agnostic: drops common PT/EN Boolean connectors and
  # keeps tokens that look like gene symbols (start with a letter, length >= 2).
  connectors <- c('E','OU','OR','AND','NOT','NAO','COM','SEM','DE','DO','DA','A','O')
  raw <- unlist(strsplit(display_name, '[^A-Za-z0-9]+'))
  raw <- toupper(trimws(raw)); raw <- raw[nzchar(raw)]
  raw <- raw[!raw %in% connectors]
  raw <- raw[grepl('^[A-Z][A-Z0-9]+$', raw)]
  unique(raw)
}

tcga_read_model <- function(cfg) {
  path <- cfg$model
  if (grepl('\\.zginml$', path, ignore.case = TRUE)) {
    z <- utils::unzip(path, list = TRUE)$Name
    z <- z[grepl('(^|/)ginml$', z, ignore.case = TRUE) | grepl('\\.ginml$', z, ignore.case = TRUE)]
    if (length(z) != 1L) stop('Nao foi possivel identificar um unico XML GINML no modelo.')
    con <- unz(path, z); on.exit(close(con), add = TRUE)
    doc <- xml2::read_xml(paste(readLines(con, warn = FALSE), collapse = '\n'))
  } else doc <- xml2::read_xml(path)
  nodes <- xml2::xml_find_all(doc, './/node')
  ids <- xml2::xml_attr(nodes, 'id'); names <- xml2::xml_attr(nodes, 'name')
  if (!length(ids) || anyNA(ids) || anyDuplicated(ids)) stop('IDs de nos invalidos/duplicados no modelo.')
  explicit <- !is.na(names) & nzchar(trimws(names))
  names[!explicit] <- ids[!explicit]
  symbol <- toupper(trimws(names))
  # Names take precedence; the few obsolete symbols below are explicitly audited.
  default_aliases <- c(DFNA5 = 'GSDME', P53 = 'TP53', 'MIR-204-5P' = 'MIR204',
                       'MIR_204_5P' = 'MIR204', 'HSA-MIR-204-5P' = 'MIR204')
  if (!is.null(cfg$aliases)) {
    a <- utils::read.csv(cfg$aliases, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c('model_name', 'symbol') %in% names(a))) stop('Aliases: colunas model_name e symbol obrigatorias.')
    k <- toupper(trimws(a$model_name))
    if (anyDuplicated(k)) stop('Aliases duplicados.')
    default_aliases[k] <- toupper(trimws(a$symbol))
  }
  hit <- symbol %in% names(default_aliases)
  symbol[hit] <- unname(default_aliases[symbol[hit]])
  status <- ifelse(symbol == 'MIR204', 'mature_miRNA_assay', 'RNAseq_candidate')
  conceptual <- c('DDR', 'MOMP', 'P53-ARREST', 'P53-KILLER', 'PROLIFERATION',
                  'RESISTANCE', 'PYROPTOSIS', 'APOPTOSIS', 'CELL_CYCLE_ARREST')
  is_conceptual <- symbol %in% conceptual
  is_complex <- !is_conceptual & grepl('[ ()&|]', symbol)
  status[is_conceptual] <- 'not_directly_measurable_conceptual_state'
  status[is_complex] <- 'not_directly_measurable_boolean_complex'
  base <- data.frame(model_id = ids, model_name = names,
    name_source = ifelse(explicit, 'name', 'id_fallback_missing_name'),
    assay_symbol = symbol, alias_applied = hit, status = status,
    component_of_model_id = NA_character_, component_of_label = NA_character_,
    stringsAsFactors = FALSE)
  # A composite node (logical AND/OR of several genes, e.g. a protein complex)
  # has no single transcript. Instead of leaving an empty panel, expand it into
  # its individual component genes so each can be plotted on its own if measured.
  extra <- lapply(which(is_complex), function(i) {
    comps <- .tcga_complex_components(names[i])
    comps <- setdiff(comps, base$assay_symbol)
    if (!length(comps)) return(NULL)
    data.frame(model_id = paste0(ids[i], '::', comps), model_name = comps,
      name_source = 'complex_component', assay_symbol = comps, alias_applied = FALSE,
      status = 'RNAseq_candidate', component_of_model_id = ids[i],
      component_of_label = names[i], stringsAsFactors = FALSE)
  })
  extra <- .tcga_bind(extra[!vapply(extra, is.null, logical(1))])
  if (nrow(extra)) base <- rbind(base, extra)
  base
}

tcga_prepare <- function(cfg, refresh = FALSE) {
  .tcga_packages()
  map <- tcga_read_model(cfg)
  .tcga_csv(map, file.path(cfg$out, 'tables', '01_model_name_mapping.csv'))
  rawpath <- file.path(cfg$out, 'GDC_query_snapshot.rds')
  if (file.exists(rawpath) && !refresh) {
    snapshot <- readRDS(rawpath)
  } else {
    filters <- list(op = 'and', content = list(
      list(op = 'in', content = list(field = 'cases.project.project_id', value = list('TCGA-LIHC'))),
      list(op = 'in', content = list(field = 'access', value = list('open'))),
      list(op = 'in', content = list(field = 'data_type', value = list(
        'Gene Expression Quantification', 'Isoform Expression Quantification')))))
    fields <- paste(c('file_id','file_name','file_size','md5sum','data_type','analysis.workflow_type',
      'cases.case_id','cases.submitter_id','cases.diagnoses.primary_diagnosis',
      'cases.samples.sample_id','cases.samples.submitter_id','cases.samples.sample_type'), collapse = ',')
    answer <- .tcga_api('files', list(filters = jsonlite::toJSON(filters, auto_unbox = TRUE),
                                    fields = fields, size = 10000))
    if (answer$data$pagination$total > length(answer$data$hits)) stop('Consulta GDC truncada; nao prossiga.')
    snapshot <- list(retrieved_utc = format(Sys.time(), tz = 'UTC'), response = answer,
                     gdc_status = .tcga_api('status'))
    saveRDS(snapshot, rawpath)
    jsonlite::write_json(snapshot, file.path(cfg$out, 'GDC_query_snapshot.json'),
                         auto_unbox = TRUE, pretty = TRUE)
  }
  rows <- lapply(snapshot$response$data$hits, function(h) {
    cs <- h$cases
    c1 <- if (length(cs)) cs[[1]] else list()
    ss <- c1$samples
    s1 <- if (length(ss)) ss[[1]] else list()
    dx <- unique(vapply(c1$diagnoses, function(x) .tcga_one(x$primary_diagnosis), character(1)))
    workflow <- .tcga_one(h$analysis$workflow_type)
    assay <- if (.tcga_one(h$data_type) == 'Gene Expression Quantification') 'RNAseq' else 'miRNAseq'
    reason <- 'eligible'
    if (assay == 'RNAseq' && workflow != 'STAR - Counts') reason <- 'excluded_workflow'
    if (assay == 'miRNAseq' && (workflow != 'BCGSC miRNA Profiling' ||
        !grepl('mirbase21.isoforms.quantification', h$file_name, fixed = TRUE))) reason <- 'excluded_workflow_or_annotation'
    # Conventional HCC: explicitly exclude combined and fibrolamellar diagnoses.
    if (!any(grepl('^Hepatocellular carcinoma', dx))) reason <- 'excluded_no_explicit_HCC_diagnosis'
    if (any(grepl('combined|cholangiocarcinoma|fibrolamellar', dx, ignore.case = TRUE))) reason <- 'excluded_combined_or_fibrolamellar'
    if (!.tcga_one(s1$sample_type) %in% c('Primary Tumor', 'Solid Tissue Normal')) reason <- 'excluded_sample_type'
    if (length(cs) != 1 || length(ss) != 1) reason <- 'excluded_ambiguous_sample_link'
    data.frame(file_id = h$file_id, file_name = h$file_name, file_size = h$file_size,
      md5sum = h$md5sum, assay = assay, workflow = workflow,
      case_id = .tcga_one(c1$case_id), case_barcode = .tcga_one(c1$submitter_id),
      sample_id = .tcga_one(s1$sample_id), sample_barcode = .tcga_one(s1$submitter_id),
      sample_type = .tcga_one(s1$sample_type), diagnoses = paste(dx, collapse = '; '),
      selection = reason, stringsAsFactors = FALSE)
  })
  audit <- .tcga_bind(rows)
  if (!nrow(audit)) stop('Nenhum arquivo retornado pelo GDC.')
  e <- audit[audit$selection == 'eligible', , drop = FALSE]
  # One sample per case/type/assay, prefer the same sample across the two assays.
  selected <- character()
  groups <- split(e, paste(e$case_id, e$sample_type, sep = '|'))
  for (g in groups) {
    shared <- intersect(g$sample_barcode[g$assay == 'RNAseq'], g$sample_barcode[g$assay == 'miRNAseq'])
    for (assay in unique(g$assay)) {
      a <- g[g$assay == assay, , drop = FALSE]
      if (length(shared)) a <- a[a$sample_barcode == sort(shared)[1], , drop = FALSE]
      a <- a[order(a$sample_barcode, a$file_name, a$file_id), , drop = FALSE]
      selected <- c(selected, a$file_id[1])
    }
  }
  audit$selection[audit$selection == 'eligible'] <- 'excluded_duplicate_case_type_assay'
  audit$selection[audit$file_id %in% selected] <- 'selected'
  manifest <- audit[audit$selection == 'selected', , drop = FALSE]
  manifest <- manifest[order(manifest$assay, manifest$case_barcode, manifest$sample_type), ]
  if (!all(c('RNAseq','miRNAseq') %in% manifest$assay)) stop('Falta um dos ensaios na coorte selecionada.')
  .tcga_csv(audit, file.path(cfg$out, 'tables', '02_all_files_selection_audit.csv'))
  .tcga_csv(manifest, file.path(cfg$out, 'tables', '03_download_manifest.csv'))
  job <- list(cfg = cfg, mapping = map, manifest = manifest, retrieved_utc = snapshot$retrieved_utc,
              model_md5 = unname(tools::md5sum(cfg$model)))
  saveRDS(job, file.path(cfg$out, 'TCGA_job.rds'))
  print(with(manifest, table(assay, sample_type)))
  message('Download previsto: ', round(sum(manifest$file_size)/1024^3, 2), ' GiB. Manifesto salvo.')
  job
}

.tcga_file <- function(job, row) file.path(job$cfg$out, 'cache', paste0(row$file_id, '.tsv'))
.tcga_valid <- function(path, row) {
  file.exists(path) && isTRUE(file.info(path)$size == as.numeric(row$file_size)) &&
    identical(tolower(unname(tools::md5sum(path))), tolower(as.character(row$md5sum)))
}
tcga_download <- function(job) {
  .tcga_packages()
  m <- job$manifest
  for (i in seq_len(nrow(m))) {
    row <- m[i, ]; path <- .tcga_file(job, row)
    if (.tcga_valid(path, row)) next
    message('[', i, '/', nrow(m), '] ', row$assay, ' ', row$sample_barcode)
    last <- ''
    ok <- FALSE
    for (attempt in seq_len(job$cfg$attempts)) {
      partial <- paste0(path, '.part')
      ok <- tryCatch({
        r <- httr::GET(paste0('https://api.gdc.cancer.gov/data/', row$file_id),
                       httr::timeout(job$cfg$timeout), httr::write_disk(partial, overwrite = TRUE))
        httr::stop_for_status(r)
        if (!.tcga_valid(partial, row)) stop('Tamanho/MD5 do arquivo nao confere.')
        if (file.exists(path)) unlink(path)
        if (!file.rename(partial, path)) stop('Nao foi possivel finalizar o arquivo local.')
        TRUE
      }, error = function(e) { last <<- conditionMessage(e); FALSE })
      if (ok) break
      unlink(partial)
      if (attempt < job$cfg$attempts) Sys.sleep(2^attempt)
    }
    if (!ok) stop('Falha em ', row$file_id, ': ', last,
      '\nExecute tcga_download(job) novamente; arquivos completos serao reutilizados.')
  }
  message('Download completo e MD5 conferido.')
  invisible(job)
}

.tcga_read_rna <- function(path, symbols) {
  d <- utils::read.delim(path, comment.char = '#', quote = '', check.names = FALSE,
                        stringsAsFactors = FALSE)
  required <- c('gene_id','gene_name','unstranded','tpm_unstranded')
  if (!all(required %in% names(d))) stop('STAR Counts: colunas inesperadas em ', path)
  labels <- toupper(d$gene_name)
  labels[labels == 'DFNA5'] <- 'GSDME'
  .tcga_bind(lapply(symbols, function(s) {
    ix <- which(labels == s)
    # Ambiguous gene symbols are not silently summed across different genes.
    status <- if (!length(ix)) 'not_in_annotation' else if (length(ix) > 1) 'ambiguous_gene_symbol' else 'measured'
    data.frame(feature = s, normalized = if (length(ix) == 1) as.numeric(d$tpm_unstranded[ix]) else NA_real_,
      read_count = if (length(ix) == 1) as.numeric(d$unstranded[ix]) else NA_real_,
      unit = 'TPM', measurement_status = status,
      annotation_id = paste(d$gene_id[ix], collapse = ';'), crossmapped_excluded = NA_real_,
      stringsAsFactors = FALSE)
  }))
}
.tcga_read_mirna <- function(path) {
  d <- utils::read.delim(path, quote = '', check.names = FALSE, stringsAsFactors = FALSE)
  required <- c('miRNA_ID','isoform_coords','read_count','reads_per_million_miRNA_mapped','cross-mapped','miRNA_region')
  if (!all(required %in% names(d))) stop('Isoformas: colunas inesperadas em ', path)
  if (!any(grepl('mature,MIMAT', d$miRNA_region, fixed = TRUE))) stop('Anotacao de miRNA maduro ausente: ', path)
  ix <- grepl('(^|,)MIMAT0000265(,|$)', d$miRNA_region) & d$miRNA_ID == 'hsa-mir-204'
  if (any(ix & !d[['cross-mapped']] %in% c('N','Y'))) stop('Flag cross-mapped desconhecida.')
  keep <- ix & d[['cross-mapped']] == 'N'
  data.frame(feature = 'miR-204-5p', normalized = sum(d$reads_per_million_miRNA_mapped[keep]),
    read_count = sum(d$read_count[keep]), unit = 'RPM',
    measurement_status = if (any(keep)) 'measured_mature_isomiRs' else 'not_detected_unique_mature_isomiRs',
    annotation_id = 'MIMAT0000265', crossmapped_excluded = sum(d$read_count[ix & !keep]),
    stringsAsFactors = FALSE)
}

tcga_extract <- function(job) {
  .tcga_packages()
  if (!identical(unname(tools::md5sum(job$cfg$model)), job$model_md5))
    stop('Modelo mudou desde tcga_prepare(). Prepare novamente; o cache sera preservado.')
  symbols <- unique(job$mapping$assay_symbol[job$mapping$status == 'RNAseq_candidate'])
  mir <- any(job$mapping$status == 'mature_miRNA_assay')
  if (!length(symbols)) stop('Nao ha nomes de genes candidatos no modelo.')
  rows <- vector('list', nrow(job$manifest))
  for (i in seq_len(nrow(job$manifest))) {
    m <- job$manifest[i, ]; path <- .tcga_file(job, m)
    if (!.tcga_valid(path, m)) stop('Arquivo ausente/invalido. Execute tcga_download(job): ', m$file_id)
    if (m$assay == 'miRNAseq' && !mir) next
    d <- if (m$assay == 'RNAseq') .tcga_read_rna(path, symbols) else .tcga_read_mirna(path)
    rows[[i]] <- cbind(m[rep(1, nrow(d)), c('case_id','case_barcode','sample_id','sample_barcode','sample_type','assay','file_id')], d)
  }
  x <- .tcga_bind(rows[!vapply(rows, is.null, logical(1))])
  if (any(!is.finite(x$normalized[!is.na(x$normalized)])) || any(x$normalized < 0, na.rm = TRUE)) stop('Valores de expressao invalidos.')
  x$log_expression <- log2(x$normalized + 1)
  .tcga_csv(x, file.path(job$cfg$out, 'tables', '04_expression_long.csv'))
  map <- job$mapping
  measured <- unique(x$feature[!is.na(x$normalized)])
  candidate <- map$status == 'RNAseq_candidate'
  map$status[candidate] <- ifelse(map$assay_symbol[candidate] %in% measured, 'RNAseq_measured', 'RNAseq_unmapped_or_ambiguous')
  .tcga_csv(map, file.path(job$cfg$out, 'tables', '01_model_name_mapping.csv'))
  saveRDS(x, file.path(job$cfg$out, 'TCGA_expression.rds'))
  x
}

.tcga_pairs <- function(x) {
  a <- x[x$sample_type == 'Primary Tumor', c('case_id','feature','log_expression')]
  b <- x[x$sample_type == 'Solid Tissue Normal', c('case_id','feature','log_expression')]
  names(a)[3] <- 'tumor'; names(b)[3] <- 'normal'
  merge(a, b, by = c('case_id', 'feature'))
}
.tcga_saveplot <- function(p, cfg, name, w = 11, h = 7) {
  for (ext in c('pdf','png')) ggplot2::ggsave(file.path(cfg$out, 'figures', paste0(name, '.', ext)),
    p, width = w, height = h, units = 'in', dpi = 600, bg = 'white', limitsize = FALSE)
}

# Regenerate only the expression figures from the saved expression table.
# Every model node is represented in the full view; missing measurements are labelled.
tcga_plot_expression <- function(job, expression = NULL) {
  .tcga_packages()
  if (!identical(unname(tools::md5sum(job$cfg$model)), job$model_md5))
    stop('Modelo mudou. Execute tcga_prepare() e tcga_extract() novamente.')
  if (is.null(expression)) {
    saved <- file.path(job$cfg$out, 'TCGA_expression.rds')
    if (!file.exists(saved)) stop('Expressao salva nao encontrada. Execute tcga_extract(job).')
    expression <- readRDS(saved)
  }
  x <- expression[is.finite(expression$log_expression), , drop = FALSE]
  mapping <- job$mapping
  mapping$feature <- ifelse(mapping$assay_symbol == 'MIR204', 'miR-204-5p', mapping$assay_symbol)
  mapping$measurable <- mapping$status %in% c('RNAseq_candidate','RNAseq_measured','mature_miRNA_assay')
  mapping$display <- mapping$model_name
  for (i in seq_len(nrow(mapping))) {
    if (mapping$measurable[i]) {
      s <- mapping$feature[i]
      mapping$display[i] <- if (s == 'TP53') 'p53 (TP53)' else s
    }
  }
  # Retain separate model nodes even when they map to the same measured transcript.
  mapping$panel <- make.unique(mapping$display, sep = ' / node ')
  mapping$available <- mapping$measurable & mapping$feature %in% x$feature
  .tcga_csv(mapping, file.path(job$cfg$out, 'tables', '09_expression_figure_coverage.csv'))
  not_shown <- mapping[!mapping$available, , drop = FALSE]
  if (nrow(not_shown)) {
    message('Paineis sem expressao quantificavel (nao plotados, ver tabela 09): ',
            paste(not_shown$model_id, collapse = ', '))
  }
  six <- c('miR-204-5p','MALAT1','SIRT1','TP53','CASP3','GSDME')
  selected <- mapping[match(six, mapping$feature), , drop = FALSE]
  selected <- selected[!is.na(selected$model_id), , drop = FALSE]
  # Preferred pathway order for the full-network figure; anything not listed here
  # (e.g. newly expanded complex components) keeps the model's original order,
  # appended right after the complex node group it was expanded from.
  pathway_order <- c('DDR','ATM','PRKAA1','PPM1D','MDM2','TP53','TP53INP1','CDKN1A','MYC',
    'MALAT1','miR-204-5p','SIRT1','GSDME','E2F1','AKT1',
    'CDK4','CDK6','CCND1','CCND2','CCND3','RB1','CDC25A',
    'BBC3','BCL2','BAX','CASP9','CASP3')
  mapping$order_key <- match(mapping$feature, pathway_order)
  mapping$order_key[is.na(mapping$order_key)] <- match(mapping$component_of_model_id,
    job$mapping$model_id)[is.na(mapping$order_key)]
  mapping_ordered <- mapping[order(is.na(mapping$order_key), mapping$order_key,
                                   seq_len(nrow(mapping))), ]
  plot_set <- function(spec, title, filename, ncol, width, note = NULL) {
    spec <- spec[spec$available, , drop = FALSE]
    if (!nrow(spec)) { message('Nada quantificavel para ', filename, '; figura nao gerada.'); return(invisible(NULL)) }
    spec$panel <- factor(spec$panel, levels = unique(spec$panel))
    rows <- lapply(seq_len(nrow(spec)), function(i) {
      d <- x[x$feature == spec$feature[i], , drop = FALSE]
      data.frame(panel = spec$panel[i], group = ifelse(d$sample_type == 'Primary Tumor', 2, 1),
                 log_expression = d$log_expression)
    })
    d <- .tcga_bind(rows)
    units <- ifelse(spec$feature == 'miR-204-5p', 'RPM', 'TPM')
    labels <- paste0(as.character(spec$panel), ' [', units, ']')
    labels <- vapply(labels, function(z) paste(strwrap(z, width = 32), collapse = '\n'), character(1))
    label_map <- setNames(labels, as.character(spec$panel))
    d$group_label <- factor(d$group, levels = c(1,2), labels = c('Adjacent non-tumor','Primary HCC'))
    n_panel <- nlevels(spec$panel)
    p <- ggplot2::ggplot(d, ggplot2::aes(group, log_expression, fill = group_label)) +
      ggplot2::geom_boxplot(ggplot2::aes(group = group_label), width = .55, outlier.shape = NA, linewidth = .35) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = .13, height = 0, seed = 101),
                          size = .6, alpha = .25) +
      ggplot2::facet_wrap(~panel, scales = 'free_y', ncol = ncol,
                          labeller = ggplot2::as_labeller(label_map)) +
      ggplot2::scale_x_continuous(breaks = c(1,2), labels = c('Adjacent non-tumor','Primary HCC'), limits = c(.5,2.5)) +
      ggplot2::scale_fill_manual(values = c('#7F8C8D','#0072B2')) +
      ggplot2::labs(title = title,
        subtitle = 'Conventional HCC cases; one sample per participant and tissue type; only network components with a direct transcript measurement are shown', x = NULL,
        y = 'log2(normalized expression + 1)',
        caption = paste(c('TPM: RNA-seq. RPM: mature miR-204-5p isomiRs. Expression does not establish protein activation or cell death.',
          note), collapse = '\n')) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = 'none',
        strip.background = ggplot2::element_rect(fill = '#EDF2F5'),
        axis.text.x = ggplot2::element_text(angle = 20, hjust = 1), panel.grid.minor = ggplot2::element_blank())
    .tcga_saveplot(p, job$cfg, filename, width, 2.6*ceiling(n_panel/ncol)+1.7)
    invisible(p)
  }
  a <- plot_set(selected, 'TCGA-LIHC expression of network components',
                'Figure_TCGA_01_six_components_expression', 3, 12)
  b <- plot_set(mapping_ordered, 'TCGA-LIHC expression across the model network',
                'Figure_TCGA_01b_all_model_components_expression', 4, 15,
                note = 'CDK4, CDK6, CCND1, CCND2 and CCND3 are shown individually: the model node CDK4_6_CyclinD is a logical AND/OR complex with no single transcript.')
  message('Figuras de expressao atualizadas em: ', file.path(job$cfg$out, 'figures'))
  invisible(list(six_components = a, all_components = b))
}

tcga_analyze <- function(job, expression = NULL) {
  .tcga_packages()
  if (!identical(unname(tools::md5sum(job$cfg$model)), job$model_md5))
    stop('Modelo mudou. Execute tcga_prepare() e tcga_extract() novamente.')
  if (is.null(expression)) expression <- tcga_extract(job)
  if (anyDuplicated(expression[c('case_id','sample_type','feature')]))
    stop('Expressao contem mais de uma observacao por participante/tipo/feature.')
  x <- expression[is.finite(expression$log_expression), , drop = FALSE]
  if (!nrow(x)) stop('Nenhuma expressao quantificavel.')
  pairs <- .tcga_pairs(x)
  .tcga_csv(pairs, file.path(job$cfg$out, 'tables', '05_matched_tumor_normal.csv'))
  tests <- lapply(sort(unique(x$feature)), function(f) {
    p <- pairs[pairs$feature == f & is.finite(pairs$tumor) & is.finite(pairs$normal), ]
    delta <- p$tumor - p$normal
    pv <- NA_real_
    if (nrow(p) >= job$cfg$minimum_pairs) pv <- if (all(delta == 0)) 1 else
      stats::wilcox.test(p$tumor, p$normal, paired = TRUE, exact = FALSE)$p.value
    data.frame(feature = f, paired_cases = nrow(p),
      median_paired_log2_difference = if (length(delta)) median(delta) else NA_real_,
      p_value = pv, status = if (nrow(p) < job$cfg$minimum_pairs) 'insufficient_pairs' else 'tested')
  })
  tests <- .tcga_bind(tests); tests$FDR_BH <- p.adjust(tests$p_value, 'BH')
  .tcga_csv(tests, file.path(job$cfg$out, 'tables', '06_paired_expression_tests.csv'))
  descriptive <- .tcga_bind(lapply(split(x, paste(x$feature, x$sample_type)), function(d)
    data.frame(feature = d$feature[1], sample_type = d$sample_type[1], n_cases = nrow(d),
      detected_n = sum(d$read_count > 0), median = median(d$normalized),
      q25 = unname(quantile(d$normalized, .25)), q75 = unname(quantile(d$normalized, .75)), unit = d$unit[1])))
  .tcga_csv(descriptive, file.path(job$cfg$out, 'tables', '07_expression_summary.csv'))
  priority <- c('miR-204-5p','MALAT1','SIRT1','TP53','BAX','CASP3','GSDME')
  tcga_plot_expression(job, expression)
  eligible <- tests[is.finite(tests$FDR_BH), ]
  if (nrow(eligible)) {
    eligible$label <- sprintf('%s (n=%d)', eligible$feature, eligible$paired_cases)
    eligible$significant <- factor(ifelse(eligible$FDR_BH < .05,'FDR < 0.05','FDR >= 0.05'))
    p <- ggplot2::ggplot(eligible, ggplot2::aes(median_paired_log2_difference,
         reorder(label, median_paired_log2_difference), color = significant)) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, color = 'grey60') + ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = c('FDR < 0.05' = '#0072B2','FDR >= 0.05' = '#888888')) +
      ggplot2::labs(title = 'Paired tumor-adjacent tissue expression differences',
        subtitle = 'Paired Wilcoxon tests; Benjamini-Hochberg correction across measured network features',
        x = 'Median paired difference in log2(normalized expression + 1)', y = NULL, color = NULL,
        caption = 'Positive: higher in tumor. This transformed difference is not a count-based differential-expression log2 fold change.') +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = 'bottom')
    .tcga_saveplot(p, job$cfg, 'Figure_TCGA_02_paired_differences', 12, max(5, .27*nrow(eligible)+2))
  }
  # Supplementary: per-patient tumor-vs-normal trajectories for the priority axis.
  # Complements Figure 02 (which only shows the summary median) by making
  # inter-patient heterogeneity around that median visible.
  traj <- pairs[pairs$feature %in% priority & is.finite(pairs$tumor) & is.finite(pairs$normal), ]
  if (nrow(traj)) {
    traj$feature <- factor(traj$feature, levels = priority[priority %in% traj$feature])
    long <- rbind(
      data.frame(case_id = traj$case_id, feature = traj$feature, group = 1, log_expression = traj$normal),
      data.frame(case_id = traj$case_id, feature = traj$feature, group = 2, log_expression = traj$tumor))
    long$group_label <- factor(long$group, levels = c(1,2), labels = c('Adjacent non-tumor','Primary HCC'))
    p <- ggplot2::ggplot(long, ggplot2::aes(group, log_expression, group = case_id)) +
      ggplot2::geom_line(alpha = .18, linewidth = .3, color = '#333333') +
      ggplot2::geom_point(ggplot2::aes(color = group_label), size = .9, alpha = .5) +
      ggplot2::facet_wrap(~feature, scales = 'free_y', ncol = 4) +
      ggplot2::scale_x_continuous(breaks = c(1,2), labels = levels(long$group_label), limits = c(.7,2.3)) +
      ggplot2::scale_color_manual(values = c('#7F8C8D','#0072B2')) +
      ggplot2::labs(title = 'Per-patient paired trajectories across the priority axis',
        subtitle = 'Each line links the same participant\'s adjacent non-tumor and primary HCC samples',
        x = NULL, y = 'log2(normalized expression + 1)', color = NULL,
        caption = 'Same paired cases as Figure 02; lines show individual heterogeneity around the paired median shift.') +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = 'bottom',
        strip.background = ggplot2::element_rect(fill = '#EDF2F5'),
        axis.text.x = ggplot2::element_text(angle = 20, hjust = 1), panel.grid.minor = ggplot2::element_blank())
    .tcga_saveplot(p, job$cfg, 'Figure_TCGA_02b_paired_trajectories', 12,
                   2.6*ceiling(nlevels(traj$feature)/4)+1.7)
  }
  # Correlations must refer to the SAME biological sample, never just the same patient.
  tumors <- x[x$sample_type == 'Primary Tumor' & x$feature %in% priority, ]
  features <- priority[priority %in% tumors$feature]
  correlations <- data.frame()
  if (length(features) > 1) correlations <- .tcga_bind(lapply(combn(features, 2, simplify = FALSE), function(fs) {
    a <- tumors[tumors$feature == fs[1], c('case_id','sample_id','log_expression')]
    b <- tumors[tumors$feature == fs[2], c('case_id','sample_id','log_expression')]
    z <- merge(a, b, by = c('case_id','sample_id'))
    ok <- nrow(z) >= job$cfg$minimum_correlation_n &&
      length(unique(z$log_expression.x)) > 1 && length(unique(z$log_expression.y)) > 1
    test <- if (ok) stats::cor.test(z$log_expression.x, z$log_expression.y, method = 'spearman', exact = FALSE) else NULL
    data.frame(feature_1 = fs[1], feature_2 = fs[2], matched_tumors = nrow(z),
      rho = if (ok) unname(test$estimate) else NA_real_, p_value = if (ok) test$p.value else NA_real_,
      status = if (ok) 'tested' else 'insufficient_samples_or_constant')
  }))
  if (nrow(correlations)) {
    correlations$FDR_BH <- p.adjust(correlations$p_value, 'BH')
    .tcga_csv(correlations, file.path(job$cfg$out, 'tables', '08_tumor_correlations.csv'))
    cplot <- correlations[is.finite(correlations$rho), ]
    if (nrow(cplot)) {
      cplot$label <- sprintf('%.2f%s\nn=%d', cplot$rho, ifelse(cplot$FDR_BH < .05, '*', ''), cplot$matched_tumors)
      p <- ggplot2::ggplot(cplot, ggplot2::aes(factor(feature_1, levels = priority),
         factor(feature_2, levels = rev(priority)), fill = rho)) +
        ggplot2::geom_tile(color = 'white') + ggplot2::geom_text(ggplot2::aes(label = label), size = 3.3) +
        ggplot2::scale_fill_gradient2(low = '#B35806', mid = 'white', high = '#2166AC', limits = c(-1,1)) +
        ggplot2::coord_fixed() + ggplot2::labs(title = 'Expression associations in TCGA-LIHC primary tumors',
          subtitle = 'Spearman correlations; RNA and miRNA matched by biological sample', x = NULL, y = NULL,
          fill = 'Spearman rho', caption = '* BH FDR < 0.05 across tested pairs. Unadjusted observational associations; no causal inference.') +
        ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(panel.grid = ggplot2::element_blank())
      .tcga_saveplot(p, job$cfg, 'Figure_TCGA_03_tumor_associations', 10, 8)
    }
  }
  capture.output(sessionInfo(), file = file.path(job$cfg$out, 'sessionInfo.txt'))
  writeLines(c('TCGA-LIHC expression support; not causal validation of Boolean rules.',
    paste('Module version:', job$cfg$version), paste('Model MD5:', job$model_md5),
    paste('GDC retrieval UTC:', job$retrieved_utc),
    'Simulations, DDR input, seeds, consensus and epistasis have not been modified.',
    'No expression threshold is interpreted as Boolean activation.',
    'No survival, protein cleavage, cell death or tumor purity measurements are inferred.',
    'Mature miRNA: MIMAT0000265, sum of non-cross-mapped isomiR RPM from miRBase21 files.',
    'Absent qualifying isomiRs in a valid file are recorded as not detected (zero reads), not biological absence.',
    'Normal samples are adjacent non-tumor tissue from HCC participants, not healthy controls.',
    'See selection audit for histology filters, additional diagnoses and duplicate exclusions.'),
    file.path(job$cfg$out, 'analysis_notes.txt'))
  message('Concluido. Tabelas e figuras em: ', job$cfg$out)
  invisible(list(expression = expression, paired_tests = tests, correlations = correlations))
}
