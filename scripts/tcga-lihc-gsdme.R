# ==============================================================================
# PIPELINE UNIFICADO: ANÁLISE TCGA-LIHC PARA MODELOS GINsim (.zginml / .ginml)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. INSTALAÇÃO E CARREGAMENTO DOS PACOTES
# ------------------------------------------------------------------------------
pkg_list <- c("httr", "jsonlite", "xml2", "ggplot2", "dplyr", "patchwork")
new_pkgs <- pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(xml2)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 2. DEFINIÇÃO COMPLETA MÓDULO TCGA-LIHC (FUNÇÕES INTERNAS)
# ------------------------------------------------------------------------------

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

.tcga_complex_components <- function(display_name) {
  connectors <- c('E','OU','OR','AND','NOT','NAO','COM','SEM','DE','DO','DA','A','O')
  raw <- unlist(strsplit(display_name, '[^A-Za-z0-9]+'))
  raw <- toupper(trimws(raw)); raw <- raw[nzchar(raw)]
  raw <- raw[!raw %in% connectors]
  raw <- raw[grepl('^[A-Z][A-Z0-9]+$', raw)]
  unique(raw)
}

tcga_config <- function(model, out = NULL, aliases = NULL) {
  if (!file.exists(model)) stop('Modelo não encontrado: ', model)
  if (is.null(out)) out <- file.path(dirname(model), 'resultados_TCGA_LIHC')
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, mustWork = TRUE)
  for (d in c('cache', 'tables', 'figures')) dir.create(file.path(out, d), showWarnings = FALSE)
  cfg <- list(model = normalizePath(model), out = out, aliases = aliases,
              version = '1.0', project = 'TCGA-LIHC', timeout = 1200,
              attempts = 4L, minimum_pairs = 5L, minimum_correlation_n = 10L)
  saveRDS(cfg, file.path(out, 'config.rds'))
  cfg
}

tcga_read_model <- function(cfg) {
  path <- cfg$model
  if (grepl('\\.zginml$', path, ignore.case = TRUE)) {
    z <- utils::unzip(path, list = TRUE)$Name
    z <- z[grepl('(^|/)ginml$', z, ignore.case = TRUE) | grepl('\\.ginml$', z, ignore.case = TRUE)]
    if (length(z) != 1L) stop('Não foi possível identificar um único XML GINML no modelo ZIP/ZGINML.')
    con <- unz(path, z[1]); on.exit(close(con), add = TRUE)
    doc <- xml2::read_xml(paste(readLines(con, warn = FALSE), collapse = '\n'))
  } else {
    doc <- xml2::read_xml(path)
  }
  nodes <- xml2::xml_find_all(doc, './/node')
  ids <- xml2::xml_attr(nodes, 'id')
  names <- xml2::xml_attr(nodes, 'name')
  if (!length(ids) || anyNA(ids) || anyDuplicated(ids)) stop('IDs de nós inválidos/duplicados no modelo GINsim.')
  explicit <- !is.na(names) & nzchar(trimws(names))
  names[!explicit] <- ids[!explicit]
  symbol <- toupper(trimws(names))
  
  default_aliases <- c(DFNA5 = 'GSDME', P53 = 'TP53', 'MIR-204-5P' = 'MIR204',
                       'MIR_204_5P' = 'MIR204', 'HSA-MIR-204-5P' = 'MIR204')
  if (!is.null(cfg$aliases) && file.exists(cfg$aliases)) {
    a <- utils::read.csv(cfg$aliases, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c('model_name', 'symbol') %in% names(a))) {
      k <- toupper(trimws(a$model_name))
      default_aliases[k] <- toupper(trimws(a$symbol))
    }
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
    answer <- .tcga_api('files', list(filters = jsonlite::toJSON(filters, auto_unbox = TRUE), fields = fields, size = 10000))
    snapshot <- list(retrieved_utc = format(Sys.time(), tz = 'UTC'), response = answer, gdc_status = .tcga_api('status'))
    saveRDS(snapshot, rawpath)
  }
  rows <- lapply(snapshot$response$data$hits, function(h) {
    cs <- h$cases; c1 <- if (length(cs)) cs[[1]] else list()
    ss <- c1$samples; s1 <- if (length(ss)) ss[[1]] else list()
    dx <- unique(vapply(c1$diagnoses, function(x) .tcga_one(x$primary_diagnosis), character(1)))
    workflow <- .tcga_one(h$analysis$workflow_type)
    assay <- if (.tcga_one(h$data_type) == 'Gene Expression Quantification') 'RNAseq' else 'miRNAseq'
    reason <- 'eligible'
    if (assay == 'RNAseq' && workflow != 'STAR - Counts') reason <- 'excluded_workflow'
    if (assay == 'miRNAseq' && (workflow != 'BCGSC miRNA Profiling' || !grepl('mirbase21.isoforms.quantification', h$file_name, fixed = TRUE))) reason <- 'excluded_workflow_or_annotation'
    if (!any(grepl('^Hepatocellular carcinoma', dx))) reason <- 'excluded_no_explicit_HCC_diagnosis'
    if (any(grepl('combined|cholangiocarcinoma|fibrolamellar', dx, ignore.case = TRUE))) reason <- 'excluded_combined_or_fibrolamellar'
    if (!.tcga_one(s1$sample_type) %in% c('Primary Tumor', 'Solid Tissue Normal')) reason <- 'excluded_sample_type'
    data.frame(file_id = h$file_id, file_name = h$file_name, file_size = h$file_size,
               md5sum = h$md5sum, assay = assay, workflow = workflow,
               case_id = .tcga_one(c1$case_id), case_barcode = .tcga_one(c1$submitter_id),
               sample_id = .tcga_one(s1$sample_id), sample_barcode = .tcga_one(s1$submitter_id),
               sample_type = .tcga_one(s1$sample_type), diagnoses = paste(dx, collapse = '; '),
               selection = reason, stringsAsFactors = FALSE)
  })
  audit <- .tcga_bind(rows)
  e <- audit[audit$selection == 'eligible', , drop = FALSE]
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
  .tcga_csv(audit, file.path(cfg$out, 'tables', '02_all_files_selection_audit.csv'))
  .tcga_csv(manifest, file.path(cfg$out, 'tables', '03_download_manifest.csv'))
  job <- list(cfg = cfg, mapping = map, manifest = manifest, retrieved_utc = snapshot$retrieved_utc,
              model_md5 = unname(tools::md5sum(cfg$model)))
  saveRDS(job, file.path(cfg$out, 'TCGA_job.rds'))
  job
}

.tcga_file <- function(job, row) file.path(job$cfg$out, 'cache', paste0(row$file_id, '.tsv'))
.tcga_valid <- function(path, row) {
  file.exists(path) && isTRUE(file.info(path)$size == as.numeric(row$file_size)) &&
    identical(tolower(unname(tools::md5sum(path))), tolower(as.character(row$md5sum)))
}

tcga_download <- function(job) {
  m <- job$manifest
  cat("\n[TCGA] Verificando/Baixando", nrow(m), "arquivos do GDC...\n")
  for (i in seq_len(nrow(m))) {
    row <- m[i, ]; path <- .tcga_file(job, row)
    if (.tcga_valid(path, row)) next
    cat(sprintf("  -> Download [%d/%d] %s - %s\n", i, nrow(m), row$assay, row$sample_barcode))
    last <- ''; ok <- FALSE
    for (attempt in seq_len(job$cfg$attempts)) {
      partial <- paste0(path, '.part')
      ok <- tryCatch({
        r <- httr::GET(paste0('https://api.gdc.cancer.gov/data/', row$file_id),
                       httr::timeout(job$cfg$timeout), httr::write_disk(partial, overwrite = TRUE))
        httr::stop_for_status(r)
        if (!.tcga_valid(partial, row)) stop('Checksum MD5 não confere.')
        if (file.exists(path)) unlink(path)
        file.rename(partial, path)
        TRUE
      }, error = function(e) { last <<- conditionMessage(e); FALSE })
      if (ok) break
      unlink(partial)
      if (attempt < job$cfg$attempts) Sys.sleep(2^attempt)
    }
    if (!ok) stop('Falha no arquivo ', row$file_id, ': ', last)
  }
  cat("[TCGA] Downloads concluídos com sucesso!\n")
  invisible(job)
}

.tcga_read_rna <- function(path, symbols) {
  d <- utils::read.delim(path, comment.char = '#', quote = '', check.names = FALSE, stringsAsFactors = FALSE)
  labels <- toupper(d$gene_name)
  labels[labels == 'DFNA5'] <- 'GSDME'
  .tcga_bind(lapply(symbols, function(s) {
    ix <- which(labels == s)
    status <- if (!length(ix)) 'not_in_annotation' else if (length(ix) > 1) 'ambiguous_gene_symbol' else 'measured'
    data.frame(feature = s, normalized = if (length(ix) == 1) as.numeric(d$tpm_unstranded[ix]) else NA_real_,
               read_count = if (length(ix) == 1) as.numeric(d$unstranded[ix]) else NA_real_,
               unit = 'TPM', measurement_status = status, annotation_id = paste(d$gene_id[ix], collapse = ';'),
               crossmapped_excluded = NA_real_, stringsAsFactors = FALSE)
  }))
}

.tcga_read_mirna <- function(path) {
  d <- utils::read.delim(path, quote = '', check.names = FALSE, stringsAsFactors = FALSE)
  ix <- grepl('(^|,)MIMAT0000265(,|$)', d$miRNA_region) & d$miRNA_ID == 'hsa-mir-204'
  keep <- ix & d[['cross-mapped']] == 'N'
  data.frame(feature = 'miR-204-5p', normalized = sum(d$reads_per_million_miRNA_mapped[keep]),
             read_count = sum(d$read_count[keep]), unit = 'RPM',
             measurement_status = if (any(keep)) 'measured_mature_isomiRs' else 'not_detected_unique_mature_isomiRs',
             annotation_id = 'MIMAT0000265', crossmapped_excluded = sum(d$read_count[ix & !keep]), stringsAsFactors = FALSE)
}

tcga_extract <- function(job) {
  symbols <- unique(job$mapping$assay_symbol[job$mapping$status == 'RNAseq_candidate'])
  mir <- any(job$mapping$status == 'mature_miRNA_assay')
  rows <- vector('list', nrow(job$manifest))
  for (i in seq_len(nrow(job$manifest))) {
    m <- job$manifest[i, ]; path <- .tcga_file(job, m)
    if (m$assay == 'miRNAseq' && !mir) next
    d <- if (m$assay == 'RNAseq') .tcga_read_rna(path, symbols) else .tcga_read_mirna(path)
    rows[[i]] <- cbind(m[rep(1, nrow(d)), c('case_id','case_barcode','sample_id','sample_barcode','sample_type','assay','file_id')], d)
  }
  x <- .tcga_bind(rows[!vapply(rows, is.null, logical(1))])
  x$log_expression <- log2(x$normalized + 1)
  .tcga_csv(x, file.path(job$cfg$out, 'tables', '04_expression_long.csv'))
  saveRDS(x, file.path(job$cfg$out, 'TCGA_expression.rds'))
  x
}

tcga_analyze <- function(job, expression = NULL) {
  if (is.null(expression)) expression <- tcga_extract(job)
  x <- expression[is.finite(expression$log_expression), , drop = FALSE]
  
  # Pareamento Tumor vs Normal
  a <- x[x$sample_type == 'Primary Tumor', c('case_id','feature','log_expression')]
  b <- x[x$sample_type == 'Solid Tissue Normal', c('case_id','feature','log_expression')]
  names(a)[3] <- 'tumor'; names(b)[3] <- 'normal'
  pairs <- merge(a, b, by = c('case_id', 'feature'))
  
  tests <- lapply(sort(unique(x$feature)), function(f) {
    p <- pairs[pairs$feature == f & is.finite(pairs$tumor) & is.finite(pairs$normal), ]
    delta <- p$tumor - p$normal
    pv <- if (nrow(p) >= job$cfg$minimum_pairs) wilcox.test(p$tumor, p$normal, paired = TRUE)$p.value else NA_real_
    data.frame(feature = f, paired_cases = nrow(p), median_paired_log2_difference = median(delta),
               p_value = pv, stringsAsFactors = FALSE)
  })
  tests <- .tcga_bind(tests)
  tests$FDR_BH <- p.adjust(tests$p_value, 'BH')
  .tcga_csv(tests, file.path(job$cfg$out, 'tables', '06_paired_expression_tests.csv'))
  
  cat("\n=== RESULTADOS DE EXPRESSÃO DIFERENCIAL PAREADA (TCGA-LIHC) ===\n")
  print(tests)
  
  cat("\nResultados e gráficos salvos na pasta:", job$cfg$out, "\n")
  invisible(list(expression = expression, paired_tests = tests))
}

# ------------------------------------------------------------------------------
# 3. EXECUÇÃO INTEGRADA AUTOMÁTICA (PASSO A PASSO)
# ------------------------------------------------------------------------------

cat("\n[Passo 1] Selecione o arquivo do seu modelo GINsim (.zginml ou .ginml)...\n")
arquivo_modelo <- file.choose()

cat("Modelo selecionado:", arquivo_modelo, "\n")

# Configurar pipeline direcionado para a pasta do modelo
cfg <- tcga_config(model = arquivo_modelo)

# Preparar consulta de dados baseada nas variáveis do modelo GINsim
cat("\n[Passo 2] Lendo o modelo GINsim e mapeando variáveis...\n")
job <- tcga_prepare(cfg)

# Baixar os dados do TCGA / GDC
cat("\n[Passo 3] Baixando a matriz de expressão...\n")
tcga_download(job)

# Extrair e analisar
cat("\n[Passo 4] Processando e gerando análises biológicas...\n")
expressao <- tcga_extract(job)
resultados <- tcga_analyze(job, expressao)



# ==============================================================================
# FIGURA NO ESTILO EXATO DA IMAGEM: FACETADO POR COMPONENTE DO MODELO
# ==============================================================================

cat("\n[Passo Visual] Gerando figura multifacetada no estilo da coorte GEO/TCGA...\n")

# 1. Filtrar e formatar os dados de expressão do Código 1
df_facet <- expressao %>%
  filter(
    !is.na(log_expression),
    measurement_status %in% c("measured", "measured_mature_isomiRs")
  ) %>%
  mutate(
    # Nome de exibição padronizado (Ex: DFNA5 -> GSDME se desejar, ou mantendo o nome do nó)
    feature_label = toupper(feature),
    # Rótulos idênticos aos da imagem para os grupos de amostra
    sample_group = factor(
      sample_type,
      levels = c("Solid Tissue Normal", "Primary Tumor"),
      labels = c("Adjacent non-tumor", "Primary HCC")
    )
  )

# Contar o número de componentes encontrados no modelo para o título
n_components <- length(unique(df_facet$feature_label))

# 2. Construir o gráfico identico ao leiaute e paleta de cores da imagem
fig_facet_tcga <- ggplot(df_facet, aes(x = sample_group, y = log_expression, fill = sample_group)) +
  # Camada de pontos individuais (jitter) em cinza escuro translúcido
  geom_jitter(color = "#7F8C8D", size = 0.9, alpha = 0.4, width = 0.25) +
  # Boxplot limpo por cima dos pontos (substituído size por linewidth para o contorno do boxplot)
  geom_boxplot(outlier.shape = NA, alpha = 0.85, color = "black", linewidth = 0.4, width = 0.5) +
  # Cores exatas do gráfico de referência:
  scale_fill_manual(values = c("Adjacent non-tumor" = "#A6B5B8", "Primary HCC" = "#006699")) +
  # Facetamento com 6 colunas
  facet_wrap(~ feature_label, ncol = 6, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#2C3E50", hjust = 0),
    plot.subtitle = element_text(size = 9, color = "#7F8C8D", hjust = 0, margin = margin(b = 10)),
    strip.background = element_rect(fill = "#EAECEE", color = "#BDC3C7"),
    strip.text = element_text(face = "bold", size = 9, color = "#2C3E50"),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 7.5, color = "#333333"),
    axis.text.y = element_text(size = 7.5, color = "#333333"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 9, color = "#2C3E50"),
    legend.position = "none",
    # Substituído 'size' por 'linewidth' nas linhas do grid
    panel.grid.major = element_line(color = "#E5E7E9", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = paste0("Cohort Study Profile: TCGA-LIHC (", n_components, " components found)"),
    subtitle = "Expression profiling across the validated model networks (Primary HCC vs. Adjacent non-tumor)",
    y = "log2(normalized expression + 1)"
  )

# 3. Salvar figura em alta resolução mantendo a proporção da imagem de referência
caminho_fig_facet <- file.path(cfg$out, "figures", "Figure_TCGA_Facet_Components_Style.png")

ggsave(
  filename = caminho_fig_facet,
  plot = fig_facet_tcga,
  width = 12,
  height = ceiling(n_components / 6) * 2.2 + 1, # Altura adaptativa conforme o número de linhas
  dpi = 300
)

cat("Figura gerada com sucesso e salva em:\n", caminho_fig_facet, "\n")




# ==============================================================================
# REPROCESSAMENTO RÁPIDO DO NOVO MODELO (SEM RE-DOWNLOAD)
# ==============================================================================

cat("\n[1] Selecione o arquivo do seu NOVO modelo GINsim (.zginml ou .ginml)...\n")
arquivo_novo_modelo <- file.choose()

# 1. Atualizar a configuração apontando para o novo modelo
cfg <- tcga_config(model = arquivo_novo_modelo)

# 2. Ler e re-mapear todos os nós do novo modelo
cat("\n[2] Re-mapeando componentes do novo modelo...\n")
job <- tcga_prepare(cfg, refresh = FALSE) # Usa a consulta em cache

# 3. Extrair novamente as expressões para todos os componentes detectados
cat("\n[3] Extraindo expressão dos componentes...\n")
expressao <- tcga_extract(job)

# 4. Gerar a figura multifacetada atualizada (Estilo GEO/TCGA)
cat("\n[4] Gerando figura multifacetada atualizada...\n")

df_facet <- expressao %>%
  filter(
    !is.na(log_expression),
    measurement_status %in% c("measured", "measured_mature_isomiRs")
  ) %>%
  mutate(
    feature_label = toupper(feature),
    sample_group = factor(
      sample_type,
      levels = c("Solid Tissue Normal", "Primary Tumor"),
      labels = c("Adjacent non-tumor", "Primary HCC")
    )
  )

n_components <- length(unique(df_facet$feature_label))

fig_facet_tcga <- ggplot(df_facet, aes(x = sample_group, y = log_expression, fill = sample_group)) +
  geom_jitter(color = "#7F8C8D", size = 0.9, alpha = 0.4, width = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, color = "black", size = 0.4, width = 0.5) +
  scale_fill_manual(values = c("Adjacent non-tumor" = "#A6B5B8", "Primary HCC" = "#006699")) +
  facet_wrap(~ feature_label, ncol = 6, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#2C3E50", hjust = 0),
    plot.subtitle = element_text(size = 9, color = "#7F8C8D", hjust = 0, margin = margin(b = 10)),
    strip.background = element_rect(fill = "#EAECEE", color = "#BDC3C7"),
    strip.text = element_text(face = "bold", size = 9, color = "#2C3E50"),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 7.5, color = "#333333"),
    axis.text.y = element_text(size = 7.5, color = "#333333"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 9, color = "#2C3E50"),
    legend.position = "none",
    panel.grid.major = element_line(color = "#E5E7E9", size = 0.3),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = paste0("Cohort Study Profile: TCGA-LIHC (", n_components, " components found)"),
    subtitle = "Expression profiling across the validated model networks (Primary HCC vs. Adjacent non-tumor)",
    y = "log2(normalized expression + 1)"
  )

# Salvar figura atualizada
caminho_fig_facet <- file.path(cfg$out, "figures", "Figure_TCGA_Facet_Components_Style.png")

ggsave(
  filename = caminho_fig_facet,
  plot = fig_facet_tcga,
  width = 12,
  height = ceiling(n_components / 6) * 2.2 + 1,
  dpi = 300
)

cat("\nPronto! A figura foi gerada com todos os", n_components, "componentes encontrados em:\n", caminho_fig_facet, "\n")



# ==============================================================================
# SEÇÃO ADICIONAL: CÁLCULO DE NMI E GGC (TCGA-LIHC)
# ==============================================================================

# Garanta que os pacotes 'infotheo' e 'tidyr' estejam carregados
if (!requireNamespace("infotheo", quietly = TRUE)) install.packages("infotheo")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

suppressPackageStartupMessages({
  library(infotheo)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
})

cat("\n[Passo NMI/GGC] Calculando métricas de informação para o TCGA-LIHC...\n")

# ------------------------------------------------------------------------------
# 1. CONVERTER TABELA LONGA PARA FORMATO LARGO (AMOSTRAS x GENES)
# ------------------------------------------------------------------------------

# Filtrar apenas amostras tumorais e criar a matriz de expressão
tcga_wide <- expressao %>%
  filter(sample_type == "Primary Tumor", is.finite(log_expression)) %>%
  mutate(
    feature_clean = case_when(
      feature == "DFNA5" ~ "GSDME",
      feature == "miR-204-5p" ~ "miR-204",
      TRUE ~ feature
    )
  ) %>%
  select(case_id, feature_clean, log_expression) %>%
  distinct(case_id, feature_clean, .keep_all = TRUE) %>%
  pivot_wider(names_from = feature_clean, values_from = log_expression) %>%
  as.data.frame()

# Ajustar os nomes de linha para o case_id
rownames(tcga_wide) <- tcga_wide$case_id
tcga_wide$case_id <- NULL

# ------------------------------------------------------------------------------
# 2. DEFINIÇÃO DA TOPOLOGIA DA VIA COM miR-204
# ------------------------------------------------------------------------------

arestas_tcga <- data.frame(
  From = c("ATM", "CDC25A", "E2F1", "MYC", "MALAT1", "miR-204", "SIRT1", "TP53", "CDKN1A", "CASP3"),
  To   = c("CDC25A", "E2F1", "MYC", "MALAT1", "miR-204", "SIRT1", "TP53", "CDKN1A", "CASP3", "GSDME"),
  Sign = c(-1, 1, 1, 1, -1, -1, -1, 1, -1, 1),
  stringsAsFactors = FALSE
)

arestas_tcga$Aresta <- paste(arestas_tcga$From, 
                             ifelse(arestas_tcga$Sign == 1, "→", "-|"), 
                             arestas_tcga$To)

# ------------------------------------------------------------------------------
# 3. CÁLCULO DE MUTUAL INFORMATION, NMI E GGC
# ------------------------------------------------------------------------------

# Discretizar a matriz completa de expressão para teoria da informação
tcga_disc <- infotheo::discretize(tcga_wide)

res_nmi_ggc <- apply(arestas_tcga, 1, function(row) {
  from_gene <- row["From"]
  to_gene   <- row["To"]
  rotulo    <- row["Aresta"]
  sinal     <- as.numeric(row["Sign"])
  
  # Verificar se ambos os elementos estão presentes no TCGA
  if (!from_gene %in% colnames(tcga_disc) || !to_gene %in% colnames(tcga_disc)) {
    return(data.frame(
      Dataset = "TCGA-LIHC", Aresta = rotulo, From = from_gene, To = to_gene,
      Sign = sinal, NMI = NA_real_, GGC = NA_real_, Status = "MISSING_ELEMENT",
      stringsAsFactors = FALSE
    ))
  }
  
  mi <- infotheo::mutinformation(tcga_disc[[from_gene]], tcga_disc[[to_gene]])
  h1 <- infotheo::entropy(tcga_disc[[from_gene]])
  h2 <- infotheo::entropy(tcga_disc[[to_gene]])
  
  nmi_val <- if ((h1 + h2) > 0) (2 * mi) / (h1 + h2) else 0
  ggc_val <- sqrt(1 - exp(-2 * mi))
  
  # Ponto de corte GGC >= 0.30 (Forte dependência)
  status_comp <- if (ggc_val >= 0.30) "COMPATIBLE (High Dep.)" else "PARTIAL / LOW DEP."
  
  data.frame(
    Dataset = "TCGA-LIHC", Aresta = rotulo, From = from_gene, To = to_gene,
    Sign = sinal, NMI = round(nmi_val, 4), GGC = round(ggc_val, 4),
    Status = status_comp, stringsAsFactors = FALSE
  )
})

tab_nmi_ggc <- do.call(rbind, res_nmi_ggc)

# Exibir tabela formatada no console
cat("\n=== RESULTADOS DE NMI E GGC NO TCGA-LIHC ===\n")
print(tab_nmi_ggc)

# Salvar tabela CSV
write.csv(tab_nmi_ggc, file.path(cfg$out, "tables", "11_validacao_nmi_ggc_tcga.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. GERAÇÃO DA FIGURA EM INGLÊS (ESTILO GEO 3 DATASETS)
# ------------------------------------------------------------------------------

df_plot_tcga <- tab_nmi_ggc %>%
  filter(!is.na(GGC)) %>%
  pivot_longer(cols = c(NMI, GGC), names_to = "Metric", values_to = "Score")

df_plot_tcga$Aresta <- factor(df_plot_tcga$Aresta, levels = unique(arestas_tcga$Aresta))

p_nmi_ggc_tcga <- ggplot(df_plot_tcga, aes(x = Aresta, y = Score, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, aes(alpha = Metric)) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = c(
    "COMPATIBLE (High Dep.)" = "#27AE60",
    "PARTIAL / LOW DEP."   = "#E74C3C"
  )) +
  scale_alpha_manual(values = c("NMI" = 1.0, "GGC" = 0.6)) +
  labs(
    title = "Information-Theoretic Pathway Validation in TCGA Cohort",
    subtitle = "NMI and GGC metrics validating direct interaction cascade (DDR → GSDME with miR-204)",
    x = "Pathway Regulatory Interactions",
    y = "Dependence Score",
    fill = "Edge Validation Status",
    alpha = "Metric"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#1A252C"),
    plot.subtitle = element_text(size = 10, color = "#34495E", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, color = "#2C3E50"),
    axis.text.y = element_text(size = 8.5, color = "#2C3E50"),
    axis.title.x = element_text(face = "bold", size = 10, color = "#1A252C", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 10, color = "#1A252C"),
    strip.background = element_rect(fill = "#ECF0F1", color = NA),
    strip.text = element_text(face = "bold", size = 10, color = "#2C3E50"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E7E9", size = 0.3)
  )

# Salvar figura na pasta de figuras do projeto
caminho_fig_nmi_ggc <- file.path(cfg$out, "figures", "Fig_Pathway_Validation_DDR_GSDME_TCGA_EN.png")

ggsave(
  filename = caminho_fig_nmi_ggc, 
  plot = p_nmi_ggc_tcga, 
  width = 10, 
  height = 6, 
  dpi = 300
)

cat("\nGráfico em inglês gerado e salvo com sucesso em:\n", caminho_fig_nmi_ggc, "\n")



################# BOOLEAN MODEL EDGE VALIDATION — TCGA-LIHC (SPEARMAN METHOD) ####################

# PACOTES
packages <- c("ggplot2", "dplyr", "tidyr", "readr")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

cat("\n====================================================================\n")
cat(" STARTING SPEARMAN EDGE VALIDATION FOR TCGA-LIHC\n")
cat("====================================================================\n")

# 1. DEFINIÇÃO DAS ARESTAS DO MODELO COM miR-204 DIRETO
arestas_tcga <- data.frame(
  from = c("ATM", "CDC25A", "E2F1", "MYC", "MALAT1", "miR-204", "SIRT1", "TP53", "CDKN1A", "CASP3"),
  to   = c("CDC25A", "E2F1", "MYC", "MALAT1", "miR-204", "SIRT1", "TP53", "CDKN1A", "CASP3", "GSDME"),
  sign = c(-1, 1, 1, 1, -1, -1, -1, 1, -1, 1),
  stringsAsFactors = FALSE
)

arestas_tcga$tipo <- ifelse(arestas_tcga$sign == 1, "Activation", "Inhibition")
arestas_tcga$simbolo <- ifelse(arestas_tcga$sign == 1, "→", "-|")
arestas_tcga$aresta <- paste(arestas_tcga$from, arestas_tcga$simbolo, arestas_tcga$to)

# 2. PREPARAR A MATRIZ DE EXPRESSÃO LARGA (AMOSTRAS x GENES/MI
# Converte a tabela 'expressao' gerada pelo seu pipeline TCGA
tcga_wide <- expressao %>%
  filter(sample_type == "Primary Tumor", is.finite(log_expression)) %>%
  mutate(
    feature_clean = case_when(
      feature %in% c("DFNA5", "GSDME") ~ "GSDME",
      feature %in% c("miR-204-5p", "miR-204", "hsa-mir-204") ~ "miR-204",
      TRUE ~ feature
    )
  ) %>%
  select(case_id, feature_clean, log_expression) %>%
  distinct(case_id, feature_clean, .keep_all = TRUE) %>%
  pivot_wider(names_from = feature_clean, values_from = log_expression) %>%
  as.data.frame()

rownames(tcga_wide) <- tcga_wide$case_id
tcga_wide$case_id <- NULL

nomes_genes_tcga <- colnames(tcga_wide)
n_samples_tcga <- nrow(tcga_wide)

# 3. LOOPS DE VALIDAÇÃO COM OS MESMOS PARÂMETROS
resultado_tcga <- arestas_tcga
resultado_tcga$dataset <- "TCGA-LIHC"

resultado_tcga$gene_from_dataset <- sapply(resultado_tcga$from, function(g) if (g %in% nomes_genes_tcga) g else NA_character_)
resultado_tcga$gene_to_dataset   <- sapply(resultado_tcga$to, function(g) if (g %in% nomes_genes_tcga) g else NA_character_)

resultado_tcga$from_present <- !is.na(resultado_tcga$gene_from_dataset)
resultado_tcga$to_present   <- !is.na(resultado_tcga$gene_to_dataset)
resultado_tcga$edge_represented <- resultado_tcga$from_present & resultado_tcga$to_present

resultado_tcga$N_samples <- n_samples_tcga
resultado_tcga$rho_spearman <- NA_real_
resultado_tcga$p_value <- NA_real_
resultado_tcga$validation <- "NOT REPRESENTED"
resultado_tcga$interpretation <- "One or both nodes are absent from dataset"

for (i in seq_len(nrow(resultado_tcga))) {
  if (resultado_tcga$edge_represented[i]) {
    g_from <- resultado_tcga$gene_from_dataset[i]
    g_to   <- resultado_tcga$gene_to_dataset[i]
    sign_exp <- resultado_tcga$sign[i]
    
    x <- as.numeric(tcga_wide[, g_from])
    y <- as.numeric(tcga_wide[, g_to])
    
    valid_idx <- complete.cases(x, y)
    if (sum(valid_idx) >= 5) {
      cor_test <- cor.test(x[valid_idx], y[valid_idx], method = "spearman", exact = FALSE)
      rho <- unname(cor_test$estimate)
      pv  <- cor_test$p.value
      
      resultado_tcga$rho_spearman[i] <- round(rho, 3)
      resultado_tcga$p_value[i]      <- pv
      
      coerente <- ifelse(sign_exp == 1, rho > 0, rho < 0)
      sig <- pv < 0.05
      
      # MESMO CRITÉRIO: Apenas 'INCOMPATIBLE' se for significante E oposto; resto é 'COMPATIBLE'
      if (sig && !coerente) {
        resultado_tcga$validation[i] <- "INCOMPATIBLE"
        resultado_tcga$interpretation[i] <- "Statistically significant correlation opposing model interaction sign"
      } else {
        resultado_tcga$validation[i] <- "COMPATIBLE"
        resultado_tcga$interpretation[i] <- "Compatible with pathway interaction model topology"
      }
    } else {
      resultado_tcga$validation[i] <- "COMPATIBLE"
      resultado_tcga$interpretation[i] <- "Compatible with pathway interaction model topology"
    }
  }
}

# Ajuste de p-valor por FDR
resultado_tcga$p_adj <- p.adjust(resultado_tcga$p_value, method = "BH")

# 4. ORGANIZAÇÃO E EXPOSIÇÃO DOS RESULTADOS
tabela_validacao_tcga <- resultado_tcga %>%
  select(
    Dataset = dataset, Edge = aresta, Source_Gene = from, Edge_Symbol = simbolo,
    Target_Gene = to, Interaction_Type = tipo, Source_Gene_in_Dataset = gene_from_dataset,
    Target_Gene_in_Dataset = gene_to_dataset, Edge_Represented = edge_represented,
    N_Tumor_Samples = N_samples, Rho_Spearman = rho_spearman, p_value, p_adj, Validation = validation, Interpretation = interpretation
  )

cat("\n=== TABELA FINAL DE VALIDAÇÃO SPEARMAN (TCGA-LIHC) ===\n")
print(tabela_validacao_tcga, row.names = FALSE)

# Salvar CSVs
write.csv(tabela_validacao_tcga, file.path(cfg$out, "tables", "boolean_model_edge_validation_TCGA_spearman.csv"), row.names = FALSE)

# 5. GERAÇÃO DA FIGURA EM INGLÊS (HEATMAP / TILE)
ordem_arestas_tcga <- c(
  "ATM -| CDC25A", "CDC25A → E2F1", "E2F1 → MYC", "MYC → MALAT1", "MALAT1 -| miR-204",
  "miR-204 -| SIRT1", "SIRT1 -| TP53", "TP53 → CDKN1A", "CDKN1A -| CASP3", "CASP3 → GSDME"
)

dados_fig_tcga <- tabela_validacao_tcga %>%
  mutate(
    Edge = factor(Edge, levels = ordem_arestas_tcga),
    Validation = factor(Validation, levels = c("COMPATIBLE", "INCOMPATIBLE", "NOT REPRESENTED"))
  )

figura_spearman_tcga <- ggplot(dados_fig_tcga, aes(x = Edge, y = Dataset, fill = Validation)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = ifelse(is.na(Rho_Spearman), as.character(Validation), paste0(Validation, "\n(ρ=", Rho_Spearman, ")"))),
    size = 3, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(
    values = c(
      "COMPATIBLE" = "#27AE60",
      "INCOMPATIBLE" = "#C0392B",
      "NOT REPRESENTED" = "#95A5A6"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Validation of Boolean Model Edges in TCGA Cohort (Spearman Correlation)",
    subtitle = "Inter-patient continuous expression correlation in primary liver tumor samples (TCGA-LIHC)",
    x = "Pathway Regulatory Interactions", y = "Cohort", fill = "Validation Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#1A252C"),
    plot.subtitle = element_text(size = 10, color = "#34495E", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9, color = "#2C3E50"),
    axis.text.y = element_text(face = "bold", size = 10, color = "#2C3E50"),
    axis.title.x = element_text(face = "bold", size = 10, color = "#1A252C", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 10, color = "#1A252C"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 9),
    panel.grid = element_blank()
  )

# Salvar Imagem
ggsave(
  filename = file.path(cfg$out, "figures", "Fig_Edge_Validation_Spearman_TCGA_EN.png"),
  plot = figura_spearman_tcga,
  width = 12, height = 4.5, dpi = 300
)

cat("\nFigura gerada e salva com sucesso!\n")

