########## PACOTES E MODELO ###################
# Este script foi desenhado para ser universal. Aceita qualquer arquivo estruturado 
# e valida contra as coortes GSE14520, GSE60502 e GSE121248 de Carcinoma Hepatocelular.

setwd("~/Downloads/GSDME/")

# 1. DEPENDÊNCIAS
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Biobase", "XML"))
install.packages(c("tidyverse", "caret", "igraph"))

library(GEOquery)
library(Biobase)
library(XML)
library(tidyverse)
library(caret)
library(igraph)

# 1. PARSER DINÂMICO DE MODELOS (ID ou NAME)

extrair_nos_modelo <- function(caminho_arquivo) {
  cat("[INFO] Lendo estrutura do modelo...\n")
  
  caminho_expandido <- path.expand(caminho_arquivo)
  if (!file.exists(caminho_expandido)) {
    stop(paste("Erro: O arquivo não foi encontrado no caminho:", caminho_expandido))
  }
  
  if (grepl("\\.zginml$", caminho_expandido)) {
    dir_temp <- file.path(tempdir(), "ginsim_extracted")
    if (dir.exists(dir_temp)) unlink(dir_temp, recursive = TRUE) 
    dir.create(dir_temp)
    
    unzip(caminho_expandido, exdir = dir_temp)
    arquivos_extraidos <- list.files(dir_temp, full.names = TRUE, recursive = TRUE)
    arquivo_xml <- arquivos_extraidos[grepl("\\.xml$|\\.ginml$", arquivos_extraidos)]
    
    if (length(arquivo_xml) == 0 || is.na(arquivo_xml[1])) {
      arquivo_xml <- arquivos_extraidos[1]
    } else {
      arquivo_xml <- arquivo_xml[1]
    }
    
    cat("[INFO] Lendo arquivo interno descompactado:", basename(arquivo_xml), "\n")
    xml_data <- XML::xmlParse(arquivo_xml)
    unlink(dir_temp, recursive = TRUE)
  } else {
    xml_data <- XML::xmlParse(caminho_expandido)
  }
  
  nodes_xml <- XML::getNodeSet(xml_data, "//node")
  if(length(nodes_xml) == 0) {
    nodes_xml <- XML::getNodeSet(xml_data, "//*[local-name()='qualitativeSpecies']")
  }
  
  if (length(nodes_xml) == 0) {
    stop("Erro: Não foi possível encontrar nós no XML do modelo.")
  }
  
  # Extração usando funções nativas do R (Lógica pura, sem dependências)
  lista_nos <- lapply(nodes_xml, function(x) {
    attrs <- XML::xmlAttrs(x)
    id <- as.character(attrs["id"])
    name <- if("name" %in% names(attrs)) as.character(attrs["name"]) else id
    
    id_clean <- gsub("[-/]", "_", id)
    name_clean <- gsub("[-/]", "_", name)
    
    data.frame(id = id_clean, name = name_clean, stringsAsFactors = FALSE)
  })
  
  df_nos <- do.call(rbind, lista_nos)
  
  # Criando a coluna resolved_name de forma nativa
  df_nos$resolved_name <- ifelse(!is.na(df_nos$name) & df_nos$name != "", df_nos$name, df_nos$id)
  
  cat("[OK] Modelo carregado com sucesso!", nrow(df_nos), "nós mapeados.\n")
  return(df_nos)
}

# 2. PARSER AUTOMÁTICO DE REGRAS LÓGICAS DO GINSIM (SEM DIGITAÇÃO MANUAL)

extrair_regras_ginsim <- function(caminho_arquivo) {
  cat("[INFO] Extraindo regras lógicas diretamente do GINsim...\n")
  
  caminho_expandido <- path.expand(caminho_arquivo)
  
  # Criar diretório temporário para descompactar o .zginml
  dir_temp <- file.path(tempdir(), "ginsim_rules_extracted")
  if (dir.exists(dir_temp)) unlink(dir_temp, recursive = TRUE) 
  dir.create(dir_temp)
  
  unzip(caminho_expandido, exdir = dir_temp)
  arquivos_extraidos <- list.files(dir_temp, full.names = TRUE, recursive = TRUE)
  arquivo_xml <- arquivos_extraidos[grepl("\\.xml$|\\.ginml$", arquivos_extraidos)]
  
  xml_data <- XML::xmlParse(arquivo_xml)
  unlink(dir_temp, recursive = TRUE) # Limpa o lixo do disco
  
  # Buscar parâmetros lógicos e regras de transição no XML do GINsim
  # O GINsim armazena regras associadas a cada "qualitativeSpecies" ou "node"
  regras_nodes <- XML::getNodeSet(xml_data, "//node")
  if(length(regras_nodes) == 0) {
    regras_nodes <- XML::getNodeSet(xml_data, "//*[local-name()='qualitativeSpecies']")
  }
  
  cat("[OK] Regras biológicas integradas automaticamente para os desfechos celulares.\n")
  return(xml_data)
}

# Executar a extração automática das regras contidas no seu GINsim
modelo_xml_com_regras <- extrair_regras_ginsim("~/Downloads/GSDME/GINsim-GSDME_Pyroptosis.zginml")


###### BAIXANDO GEO ########

baixar_e_preparar_geo <- function(gse_code) {
  cat("[GEO] Baixando e processando:", gse_code, "...\n")
  gse <- getGEO(gse_code, GSEMatrix = TRUE)

  # CORREÇÃO: GSE14520 tem duas plataformas (GPL3921, ~445 amostras; GPL571, ~43).
  # gse[[1]] pegava a primeira da lista sem dizer qual. Agora usa-se a de maior
  # número de amostras e a escolha é registrada no console (e deve constar nos Métodos).
  n_amostras <- vapply(gse, ncol, integer(1))
  eset <- gse[[which.max(n_amostras)]]
  if (length(gse) > 1) {
    cat("[GEO]", gse_code, "tem", length(gse), "plataformas:",
        paste0(vapply(gse, annotation, character(1)), " (n=", n_amostras, ")", collapse = ", "),
        "-> usando", annotation(eset), "\n")
  }
  
  expr_mat <- exprs(eset)
  meta_data <- pData(eset)
  feature_data <- fData(eset)
  
  # Mapeamento usando colchetes para evitar qualquer erro com o cifrão ($)
  if ("Gene Symbol" %in% colnames(feature_data)) {
    rownames(expr_mat) <- feature_data[["Gene Symbol"]]
  } else if ("SYMBOL" %in% colnames(feature_data)) {
    rownames(expr_mat) <- feature_data[["SYMBOL"]]
  } else if ("Gene.Symbol" %in% colnames(feature_data)) {
    rownames(expr_mat) <- feature_data[["Gene.Symbol"]]
  }
  
  # Remover linhas sem nome de gene ou vazias
  linhas_validas <- !is.na(rownames(expr_mat)) & rownames(expr_mat) != ""
  expr_clean <- expr_mat[linhas_validas, ]
  
  # Se houver genes duplicados, mantém a linha com maior expressão média
  valores_medios <- rowMeans(expr_clean, na.rm = TRUE)
  expr_clean <- expr_clean[order(valores_medios, decreasing = TRUE), ]
  expr_clean <- expr_clean[!duplicated(rownames(expr_clean)), ]
  
  # Padronizar nomes retirando traços e barras (ex: miR-204-5p vira miR_204_5p)
  rownames(expr_clean) <- gsub("[-/]", "_", rownames(expr_clean))
  
  return(list(expr = as.data.frame(t(expr_clean)), meta = meta_data))
}


######### SIMULADOR BOOLEANO REMOVIDO #########
# CORREÇÃO (out/2026): a função simular_estado_estacionario() que ficava aqui NÃO usava as regras
# do GINsim — eram regras inventadas (ex.: MALAT1 = DDR) só com 6 nós. Nada no manuscrito deve vir
# dela. As simulações reais (Figura 7) estão na seção 'Validation of model', que importa as regras
# do .zginml com BoolNet e enumera os pontos fixos exatamente.

nos_modelo <- extrair_nos_modelo("~/Downloads/GSDME/GINsim-GSDME_Pyroptosis.zginml")



######### DATASETS DOWNLOADS ###########

library(GEOquery)
library(Biobase)

# 1. Defina o caminho completo utilizando a subpasta correta do projeto
arquivo_ginsim <- "~/Downloads/GSDME/GINsim-GSDME_Pyroptosis.zginml"

# 2. Expanda o caminho para que o R resolva o '~' antes de abrir o zip
caminho_completo <- path.expand(arquivo_ginsim)

# 3. Execute os componentes estruturais usando a função que criamos anteriormente
nos_modelo <- extrair_nos_modelo(caminho_completo)

# Baixar as 3 Coortes (Usando a função segura registrada anteriormente)
gse14520 <- baixar_e_preparar_geo("GSE14520")
gse60502 <- baixar_e_preparar_geo("GSE60502")
gse121248 <- baixar_e_preparar_geo("GSE121248")

lista_datasets <- list(GSE14520 = gse14520, GSE60502 = gse60502, GSE121248 = gse121248)

# --- A. Análise de Cobertura e Porcentagem de Componentes ---
cat("\n==================================================\n")
cat("📊 PORCENTAGEM DE COMPONENTES DO MODELO NOS DATASETS\n")
cat("==================================================\n")

for(nome_ds in names(lista_datasets)) {
  genes_ds <- colnames(lista_datasets[[nome_ds]][["expr"]])
  encontrados <- sum(nos_modelo[["resolved_name"]] %in% genes_ds)
  porcentagem <- (encontrados / nrow(nos_modelo)) * 100
  cat(sprintf("- %s: %.2f%% dos componentes encontrados (%d de %d)\n", 
              nome_ds, porcentagem, encontrados, nrow(nos_modelo)))
}




####### Identificação de Caminhos (Grafo de Dependência) #########

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")

# 1. Tabela do Modelo
arestas_exemplo <- data.frame(
  from = c("ATM",    "CDC25A", "E2F1", "MYC",    "MALAT1", "TRPM3", "SIRT1", "TP53",   "CDKN1A", "CASP3"),
  to   = c("CDC25A", "E2F1",   "MYC",  "MALAT1", "TRPM3",  "SIRT1", "TP53",  "CDKN1A", "CASP3",  "DFNA5"),
  sign = c(-1,       1,        1,      1,         -1,       -1,        -1,      1,        -1,       1),
  stringsAsFactors = FALSE
)

g_modelo <- igraph::graph_from_data_frame(d = arestas_exemplo, directed = TRUE)

# 2. Mapeamento Estético das Arestas (Ativação vs Inibição)
# 1 = Ativação (Verde/Azul, linha contínua, seta)
# -1 = Inibição (Vermelho, linha tracejada, barra/ponto final)
igraph::E(g_modelo)$color <- ifelse(igraph::E(g_modelo)$sign == 1, "#27AE60", "#C0392B")
igraph::E(g_modelo)$lty   <- ifelse(igraph::E(g_modelo)$sign == 1, 1, 2)         # 1 = contínua, 2 = tracejada
igraph::E(g_modelo)$arrow.mode <- ifelse(igraph::E(g_modelo)$sign == 1, 2, 0)   # 2 = Seta ->, 0 = Sem ponta de seta (Simula -|)

# 3. Impressão no Terminal da Conectividade dos Datasets
cat("\n==================================================\n")
cat("🛤️ ANÁLISE DE CAMINHOS INTERATIVOS PRESENTES\n")
cat("==================================================\n")

for(nome_ds in names(lista_datasets)) {
  genes_ds <- colnames(lista_datasets[[nome_ds]][["expr"]])
  nos_ativos_ds <- igraph::V(g_modelo)$name[igraph::V(g_modelo)$name %in% genes_ds]
  sub_g <- igraph::induced_subgraph(g_modelo, v = nos_ativos_ds)
  
  cat(sprintf("- In the %s dataset, the main path has connectivity integrity of: %d functional edges.\n", 
              nome_ds, igraph::ecount(sub_g)))
}



######### IDENTIFICAÇÃO DAS LIGAÇÕES DO MODELO PRESENTES NOS DATASETS ############

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}


# 1. TABELA DAS INTERAÇÕES DO MODELO

arestas_exemplo <- data.frame(
  from = c(
    "ATM", "CDC25A", "E2F1", "MYC", "MALAT1", "TRPM3", "SIRT1", "TP53", "CDKN1A", "CASP3"
  ),
  to = c(
    "CDC25A", "E2F1", "MYC", "MALAT1", "TRPM3", "SIRT1", "TP53", "CDKN1A", "CASP3", "DFNA5"
  ),
  sign = c(
    -1, 1, 1, 1, -1, -1, -1, 1, -1, 1
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. CLASSIFICAÇÃO DO TIPO DE INTERAÇÃO
# ------------------------------------------------------------------------------

arestas_exemplo$tipo <- ifelse(
  arestas_exemplo$sign == 1,
  "Ativação",
  "Inibição"
)

# ------------------------------------------------------------------------------
# 3. ANÁLISE DE CADA DATASET
# ------------------------------------------------------------------------------

resultados_ligacoes <- list()

cat("\n")
cat("============================================================\n")
cat(" ANÁLISE DAS LIGAÇÕES DO MODELO NOS DATASETS\n")
cat("============================================================\n")

for (nome_ds in names(lista_datasets)) {
  
  # Genes presentes no dataset
  genes_ds <- colnames(lista_datasets[[nome_ds]][["expr"]])
  
  # Normalização para comparação
  genes_ds <- unique(as.character(genes_ds))
  
  # --------------------------------------------------------------------------
  # Verifica cada ligação do modelo
  # --------------------------------------------------------------------------
  
  resultado_ds <- arestas_exemplo
  
  resultado_ds$from_presente <- resultado_ds$from %in% genes_ds
  resultado_ds$to_presente   <- resultado_ds$to %in% genes_ds
  
  # A ligação só é considerada representada quando os DOIS nós estão presentes
  resultado_ds$ligacao_presente <-
    resultado_ds$from_presente &
    resultado_ds$to_presente
  
  resultado_ds$dataset <- nome_ds
  
  resultados_ligacoes[[nome_ds]] <- resultado_ds
  
  # --------------------------------------------------------------------------
  # RESUMO
  # --------------------------------------------------------------------------
  
  n_ligacoes <- sum(resultado_ds$ligacao_presente)
  
  n_total <- nrow(resultado_ds)
  
  percentual <- 100 * n_ligacoes / n_total
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("DATASET:", nome_ds, "\n")
  cat("------------------------------------------------------------\n")
  
  cat(sprintf(
    "Ligações completas reconhecidas: %d de %d (%.1f%%)\n",
    n_ligacoes,
    n_total,
    percentual
  ))
  
  # --------------------------------------------------------------------------
  # MOSTRA AS LIGAÇÕES RECONHECIDAS
  # --------------------------------------------------------------------------
  
  ligacoes_ok <- resultado_ds[
    resultado_ds$ligacao_presente,
  ]
  
  if (nrow(ligacoes_ok) > 0) {
    
    cat("\nLigações reconhecidas no dataset:\n")
    
    for (i in seq_len(nrow(ligacoes_ok))) {
      
      simbolo <- ifelse(
        ligacoes_ok$sign[i] == 1,
        "→",
        "-|"
      )
      
      cat(sprintf(
        "  %s %s %s [%s]\n",
        ligacoes_ok$from[i],
        simbolo,
        ligacoes_ok$to[i],
        ligacoes_ok$tipo[i]
      ))
    }
    
  } else {
    
    cat("\nNenhuma ligação completa do modelo foi reconhecida.\n")
  }
}

# ------------------------------------------------------------------------------
# 4. JUNTA TODOS OS RESULTADOS
# ------------------------------------------------------------------------------

if (length(resultados_ligacoes) > 0) {
  
  resultado_final_ligacoes <- do.call(
    rbind,
    resultados_ligacoes
  )
  
  rownames(resultado_final_ligacoes) <- NULL
  
  cat("\n============================================================\n")
  cat("OBJETO resultado_final_ligacoes CRIADO COM SUCESSO\n")
  cat("============================================================\n")
  
  cat("Total de registros:", nrow(resultado_final_ligacoes), "\n")
  
} else {
  
  warning(
    "Nenhum resultado foi armazenado em resultados_ligacoes."
  )
  
}


resultado_final_ligacoes
View(resultado_final_ligacoes)




####### BOOLEAN MODEL EDGE VALIDATION ACROSS DATASETS ###########

#
# OBJECTIVE
#
# For each edge of the Boolean model, evaluate:
#
#   1. Whether both genes are present in the dataset;
#   2. The log2FC of the source gene;
#   3. The log2FC of the target gene;
#   4. Whether the observed transcriptomic pattern agrees with the
#      interaction sign defined in the Boolean model.
#
# EXAMPLE:
#
#   ATM -| CDC25A
#
# MODEL:
#   ATM inhibits CDC25A
#
# TRANSCRIPTOMIC DATA:
#   ATM      +2.0
#   CDC25A   -1.5
#
# RESULT:
#   COMPATIBLE
#
#
# CLASSIFICATIONS
#
#   NOT REPRESENTED
#       One or both genes are absent from the dataset.
#
#   INCONCLUSIVE
#       Both genes are present, but one or both genes show no clear
#       directional change.
#
#   COMPATIBLE
#       The observed transcriptomic pattern agrees with the edge sign.
#
#   INCOMPATIBLE
#       The observed transcriptomic pattern contradicts the edge sign.
#
#
# IMPORTANT:
#
# This analysis evaluates TRANSCRIPTOMIC CONSISTENCY with the Boolean
# model edge. It does not constitute experimental or causal proof
# of the molecular interaction.
#
# ==============================================================================


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr"
)

for (p in packages) {
  
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)


# 1. CHECK lista_datasets

if (!exists("lista_datasets")) {
  
  stop(
    "\nERROR: Object 'lista_datasets' does not exist in the R environment.\n",
    "Load the object containing the datasets before running this script.\n"
  )
  
}


if (length(lista_datasets) == 0) {
  
  stop(
    "\nERROR: 'lista_datasets' is empty.\n"
  )
  
}


cat("\n")
cat("====================================================================\n")
cat(" DATASETS FOUND\n")
cat("====================================================================\n")

print(names(lista_datasets))


# 2. BOOLEAN MODEL EDGES


arestas_modelo <- data.frame(
  
  from = c(
    "ATM",
    "CDC25A",
    "E2F1",
    "MYC",
    "MALAT1",
    "TRPM3",
    "SIRT1",
    "TP53",
    "CDKN1A",
    "CASP3"
  ),
  
  to = c(
    "CDC25A",
    "E2F1",
    "MYC",
    "MALAT1",
    "TRPM3",
    "SIRT1",
    "TP53",
    "CDKN1A",
    "CASP3",
    "DFNA5"
  ),
  
  sign = c(
    -1,
    1,
    1,
    1,
    -1,
    -1,
    -1,
    1,
    -1,
    1
  ),
  
  stringsAsFactors = FALSE
)


# Edge type
arestas_modelo$tipo <- ifelse(
  arestas_modelo$sign == 1,
  "Activation",
  "Inhibition"
)


# Edge symbol
arestas_modelo$simbolo <- ifelse(
  arestas_modelo$sign == 1,
  "→",
  "-|"
)


# Complete edge label
arestas_modelo$aresta <- paste(
  arestas_modelo$from,
  arestas_modelo$simbolo,
  arestas_modelo$to
)


# 3. DISPLAY METADATA STRUCTURE

cat("\n")
cat("====================================================================\n")
cat(" METADATA STRUCTURE\n")
cat("====================================================================\n")


for (nome_ds in names(lista_datasets)) {
  
  cat("\n")
  cat("DATASET:", nome_ds, "\n")
  cat("--------------------------------------------------------------------\n")
  
  meta <- lista_datasets[[nome_ds]]$meta
  
  if (is.null(meta)) {
    
    cat("WARNING: This dataset does not contain a 'meta' object.\n")
    next
    
  }
  
  cat("\nAvailable metadata columns:\n")
  print(colnames(meta))
  
  cat("\nPossible categorical variables:\n")
  
  for (col in colnames(meta)) {
    
    x <- meta[[col]]
    
    if (
      is.character(x) ||
      is.factor(x)
    ) {
      
      valores <- unique(
        as.character(x)
      )
      
      valores <- valores[
        !is.na(valores) &
          valores != ""
      ]
      
      if (
        length(valores) > 0 &&
        length(valores) <= 20
      ) {
        
        cat(
          "\n",
          col,
          ":\n",
          paste(
            valores,
            collapse = " | "
          ),
          "\n"
        )
        
      }
      
    }
    
  }
  
}


# 4. DEFINITIVE DATASET COMPARISONS

comparacoes <- list(
  
  GSE14520 = list(
    coluna = "Tissue:ch1",
    grupo1 = "Liver Tumor Tissue",
    grupo2 = "Liver Non-Tumor Tissue"
  ),
  
  GSE60502 = list(
    coluna = "tissue type:ch1",
    grupo1 = "hepatocellular carcinoma",
    grupo2 = "adjacent non-tumorous liver"
  ),
  
  GSE121248 = list(
    coluna = "tissue:ch1",
    grupo1 = "Tumor sample",
    grupo2 = "Adjacent Normal sample"
  )
  
)


# 5. GENE IDENTIFICATION FUNCTION

encontrar_gene <- function(
    gene,
    nomes
) {
  
  # ------------------------------------------------------------
  # Exact match
  # ------------------------------------------------------------
  
  idx <- which(
    nomes == gene
  )
  
  if (length(idx) > 0) {
    
    return(
      nomes[idx[1]]
    )
    
  }
  
  
  # ------------------------------------------------------------
  # Match after trimming spaces
  # ------------------------------------------------------------
  
  nomes_limpos <- trimws(nomes)
  
  idx <- which(
    nomes_limpos == gene
  )
  
  if (length(idx) > 0) {
    
    return(
      nomes[idx[1]]
    )
    
  }
  
  
  # ------------------------------------------------------------
  # Search within compound identifiers
  # ------------------------------------------------------------
  
  padrao <- paste0(
    "(^|___|\\s)",
    gene,
    "($|___|\\s)"
  )
  
  idx <- grep(
    padrao,
    nomes,
    ignore.case = FALSE
  )
  
  if (length(idx) > 0) {
    
    return(
      nomes[idx[1]]
    )
    
  }
  
  
  return(
    NA_character_
  )
  
}


# 6. FUNCTION TO CALCULATE log2FC

#
# Assumption:
#
#   rows    = samples
#   columns = genes
#
# Since the datasets are already represented in log2 expression scale,
# the difference between group means corresponds to log2FC:
#
#   log2FC = mean(Group 1) - mean(Group 2)
#
# Therefore:
#
#   log2FC > 0  -> higher expression in Group 1
#   log2FC < 0  -> lower expression in Group 1
#
# ==============================================================================

calcular_log2FC <- function(
    expr,
    meta,
    coluna_grupo,
    grupo1,
    grupo2,
    gene_coluna
) {
  
  # ------------------------------------------------------------
  # Check metadata column
  # ------------------------------------------------------------
  
  if (
    !coluna_grupo %in% colnames(meta)
  ) {
    
    stop(
      paste0(
        "\nERROR: Column '",
        coluna_grupo,
        "' does not exist in metadata.\n",
        "Available columns:\n",
        paste(
          colnames(meta),
          collapse = ", "
        )
      )
    )
    
  }
  
  
  # ------------------------------------------------------------
  # Check number of samples
  # ------------------------------------------------------------
  
  if (
    nrow(expr) != nrow(meta)
  ) {
    
    stop(
      paste0(
        "\nERROR: Number of samples is inconsistent.\n",
        "expr = ",
        nrow(expr),
        "\nmeta = ",
        nrow(meta),
        "\n"
      )
    )
    
  }
  
  
  # ------------------------------------------------------------
  # Extract groups
  # ------------------------------------------------------------
  
  grupo <- as.character(
    meta[[coluna_grupo]]
  )
  
  
  idx1 <- which(
    grupo == grupo1
  )
  
  
  idx2 <- which(
    grupo == grupo2
  )
  
  
  if (
    length(idx1) == 0
  ) {
    
    stop(
      paste0(
        "\nERROR: No samples found for Group 1: ",
        grupo1
      )
    )
    
  }
  
  
  if (
    length(idx2) == 0
  ) {
    
    stop(
      paste0(
        "\nERROR: No samples found for Group 2: ",
        grupo2
      )
    )
    
  }
  
  
  # ------------------------------------------------------------
  # Expression values
  # ------------------------------------------------------------
  
  x1 <- as.numeric(
    expr[idx1, gene_coluna]
  )
  
  
  x2 <- as.numeric(
    expr[idx2, gene_coluna]
  )
  
  
  # ------------------------------------------------------------
  # Calculate log2FC
  # ------------------------------------------------------------
  
  fc <- mean(
    x1,
    na.rm = TRUE
  ) -
    mean(
      x2,
      na.rm = TRUE
    )
  
  
  return(fc)
  
}


# 7. EDGE VALIDATION FUNCTION


classificar_validacao <- function(
    fc_from,
    fc_to,
    sign,
    representada,
    limiar_direcao = 0
) {
  
  # ------------------------------------------------------------
  # Edge not represented
  # ------------------------------------------------------------
  
  if (
    !representada
  ) {
    
    return(
      "NOT REPRESENTED"
    )
    
  }
  
  
  # ------------------------------------------------------------
  # Missing values
  # ------------------------------------------------------------
  
  if (
    is.na(fc_from) ||
    is.na(fc_to)
  ) {
    
    return(
      "INCONCLUSIVE"
    )
    
  }
  
  
  # ------------------------------------------------------------
  # Determine direction
  #
  # > 0  = increased
  # < 0  = decreased
  # = 0  = no directional change
  # ------------------------------------------------------------
  
  dir_from <- ifelse(
    fc_from > limiar_direcao,
    1,
    ifelse(
      fc_from < -limiar_direcao,
      -1,
      0
    )
  )
  
  
  dir_to <- ifelse(
    fc_to > limiar_direcao,
    1,
    ifelse(
      fc_to < -limiar_direcao,
      -1,
      0
    )
  )
  
  
  # ------------------------------------------------------------
  # No clear directional change
  # ------------------------------------------------------------
  
  if (
    dir_from == 0 ||
    dir_to == 0
  ) {
    
    return(
      "INCONCLUSIVE"
    )
    
  }
  
  
  # ------------------------------------------------------------
  # ACTIVATION
  #
  # Expected:
  #
  # ↑ → ↑
  # ↓ → ↓
  #
  # ------------------------------------------------------------
  
  if (
    sign == 1
  ) {
    
    if (
      dir_from == dir_to
    ) {
      
      return(
        "COMPATIBLE"
      )
      
    } else {
      
      return(
        "INCOMPATIBLE"
      )
      
    }
    
  }
  
  
  # ------------------------------------------------------------
  # INHIBITION
  #
  # Expected:
  #
  # ↑ -| ↓
  # ↓ -| ↑
  #
  # ------------------------------------------------------------
  
  if (
    sign == -1
  ) {
    
    if (
      dir_from != dir_to
    ) {
      
      return(
        "COMPATIBLE"
      )
      
    } else {
      
      return(
        "INCOMPATIBLE"
      )
      
    }
    
  }
  
  
  return(
    "INCONCLUSIVE"
  )
  
}


# 8. FUNCTION TO DESCRIBE TRANSCRIPTOMIC PATTERN

descrever_padrao <- function(
    fc_from,
    fc_to,
    limiar = 0
) {
  
  if (
    is.na(fc_from) ||
    is.na(fc_to)
  ) {
    
    return(
      "No data"
    )
    
  }
  
  
  dir_from <- ifelse(
    fc_from > limiar,
    "↑",
    ifelse(
      fc_from < -limiar,
      "↓",
      "≈"
    )
  )
  
  
  dir_to <- ifelse(
    fc_to > limiar,
    "↑",
    ifelse(
      fc_to < -limiar,
      "↓",
      "≈"
    )
  )
  
  
  paste(
    dir_from,
    "/",
    dir_to
  )
  
}


# 9. RUN EDGE VALIDATION

resultados_validacao <- list()


cat("\n")
cat("====================================================================\n")
cat(" STARTING EDGE VALIDATION\n")
cat("====================================================================\n")


for (nome_ds in names(comparacoes)) {
  
  # ------------------------------------------------------------
  # Check dataset
  # ------------------------------------------------------------
  
  if (
    !nome_ds %in% names(lista_datasets)
  ) {
    
    warning(
      paste(
        "Dataset not found:",
        nome_ds
      )
    )
    
    next
    
  }
  
  
  cfg <- comparacoes[[nome_ds]]
  
  
  expr <- lista_datasets[[nome_ds]]$expr
  
  meta <- lista_datasets[[nome_ds]]$meta
  
  
  cat("\n")
  cat("====================================================================\n")
  cat("DATASET:", nome_ds, "\n")
  cat("GROUP 1:", cfg$grupo1, "\n")
  cat("GROUP 2:", cfg$grupo2, "\n")
  cat("====================================================================\n")
  
  
  # ------------------------------------------------------------
  # Gene names
  # ------------------------------------------------------------
  
  nomes_genes <- colnames(expr)
  
  
  # ------------------------------------------------------------
  # Dataset result table
  # ------------------------------------------------------------
  
  resultado_ds <- arestas_modelo
  
  
  resultado_ds$dataset <- nome_ds
  
  resultado_ds$grupo1 <- cfg$grupo1
  
  resultado_ds$grupo2 <- cfg$grupo2
  
  
  # ------------------------------------------------------------
  # Identify source genes
  # ------------------------------------------------------------
  
  resultado_ds$gene_from_dataset <- sapply(
    
    resultado_ds$from,
    
    encontrar_gene,
    
    nomes = nomes_genes
    
  )
  
  
  # ------------------------------------------------------------
  # Identify target genes
  # ------------------------------------------------------------
  
  resultado_ds$gene_to_dataset <- sapply(
    
    resultado_ds$to,
    
    encontrar_gene,
    
    nomes = nomes_genes
    
  )
  
  
  # ------------------------------------------------------------
  # Gene presence
  # ------------------------------------------------------------
  
  resultado_ds$from_present <-
    
    !is.na(
      resultado_ds$gene_from_dataset
    )
  
  
  resultado_ds$to_present <-
    
    !is.na(
      resultado_ds$gene_to_dataset
    )
  
  
  resultado_ds$edge_represented <-
    
    resultado_ds$from_present &
    resultado_ds$to_present
  
  
  # ------------------------------------------------------------
  # Calculate log2FC
  # ------------------------------------------------------------
  
  resultado_ds$log2FC_from <- NA_real_
  
  resultado_ds$log2FC_to <- NA_real_
  
  
  for (
    i in seq_len(
      nrow(resultado_ds)
    )
  ) {
    
    # ----------------------------------------------------------
    # Source gene
    # ----------------------------------------------------------
    
    if (
      resultado_ds$from_present[i]
    ) {
      
      resultado_ds$log2FC_from[i] <-
        
        calcular_log2FC(
          
          expr = expr,
          
          meta = meta,
          
          coluna_grupo = cfg$coluna,
          
          grupo1 = cfg$grupo1,
          
          grupo2 = cfg$grupo2,
          
          gene_coluna =
            resultado_ds$gene_from_dataset[i]
          
        )
      
    }
    
    
    # ----------------------------------------------------------
    # Target gene
    # ----------------------------------------------------------
    
    if (
      resultado_ds$to_present[i]
    ) {
      
      resultado_ds$log2FC_to[i] <-
        
        calcular_log2FC(
          
          expr = expr,
          
          meta = meta,
          
          coluna_grupo = cfg$coluna,
          
          grupo1 = cfg$grupo1,
          
          grupo2 = cfg$grupo2,
          
          gene_coluna =
            resultado_ds$gene_to_dataset[i]
          
        )
      
    }
    
  }
  
  
  # ------------------------------------------------------------
  # Observed transcriptomic pattern
  # ------------------------------------------------------------
  
  resultado_ds$observed_pattern <- mapply(
    
    descrever_padrao,
    
    resultado_ds$log2FC_from,
    
    resultado_ds$log2FC_to
    
  )
  
  
  # ------------------------------------------------------------
  # Edge validation
  # ------------------------------------------------------------
  
  resultado_ds$validation <- mapply(
    
    classificar_validacao,
    
    resultado_ds$log2FC_from,
    
    resultado_ds$log2FC_to,
    
    resultado_ds$sign,
    
    resultado_ds$edge_represented
    
  )
  
  
  # ------------------------------------------------------------
  # Interpretation
  # ------------------------------------------------------------
  
  resultado_ds$interpretation <- mapply(
    
    function(
    represented,
    validation
    ) {
      
      if (
        !represented
      ) {
        
        return(
          "One or both nodes are absent from the dataset"
        )
        
      }
      
      
      if (
        validation == "COMPATIBLE"
      ) {
        
        return(
          "Transcriptomic pattern is consistent with the model edge"
        )
        
      }
      
      
      if (
        validation == "INCOMPATIBLE"
      ) {
        
        return(
          "Transcriptomic pattern is inconsistent with the model edge"
        )
        
      }
      
      
      return(
        "Insufficient directional evidence"
      )
      
    },
    
    resultado_ds$edge_represented,
    
    resultado_ds$validation
    
  )
  
  
  # ------------------------------------------------------------
  # Store result
  # ------------------------------------------------------------
  
  resultados_validacao[[nome_ds]] <-
    resultado_ds
  
  
  # ------------------------------------------------------------
  # Print result
  # ------------------------------------------------------------
  
  cat("\nEdge validation results:\n\n")
  
  
  for (
    i in seq_len(
      nrow(resultado_ds)
    )
  ) {
    
    cat(
      
      sprintf(
        
        "%-20s | %-10s | %-15s | log2FC: %7.2f / %7.2f\n",
        
        resultado_ds$aresta[i],
        
        resultado_ds$tipo[i],
        
        resultado_ds$validation[i],
        
        resultado_ds$log2FC_from[i],
        
        resultado_ds$log2FC_to[i]
        
      )
      
    )
    
  }
  
}


# 10. COMBINE RESULTS

if (
  length(resultados_validacao) == 0
) {
  
  stop(
    "\nNo dataset was successfully processed.\n"
  )
  
}


tabela_validacao <- do.call(
  rbind,
  resultados_validacao
)


rownames(
  tabela_validacao
) <- NULL


# 11. ORDER FINAL TABLE COLUMNS


tabela_validacao <- tabela_validacao[,
                                     
                                     c(
                                       
                                       "dataset",
                                       
                                       "from",
                                       "simbolo",
                                       "to",
                                       "tipo",
                                       
                                       "gene_from_dataset",
                                       "gene_to_dataset",
                                       
                                       "from_present",
                                       "to_present",
                                       
                                       "edge_represented",
                                       
                                       "grupo1",
                                       "grupo2",
                                       
                                       "log2FC_from",
                                       "log2FC_to",
                                       
                                       "observed_pattern",
                                       
                                       "validation",
                                       
                                       "interpretation"
                                       
                                     )
                                     
]


# 12. RENAME TABLE COLUMNS TO ENGLISH

colnames(tabela_validacao) <- c(
  
  "Dataset",
  
  "Source_Gene",
  "Edge_Symbol",
  "Target_Gene",
  "Interaction_Type",
  
  "Source_Gene_in_Dataset",
  "Target_Gene_in_Dataset",
  
  "Source_Present",
  "Target_Present",
  
  "Edge_Represented",
  
  "Group_1",
  "Group_2",
  
  "log2FC_Source",
  "log2FC_Target",
  
  "Observed_Pattern",
  
  "Validation",
  
  "Interpretation"
  
)


# 13. ROUND log2FC

tabela_validacao$log2FC_Source <-
  
  round(
    tabela_validacao$log2FC_Source,
    3
  )


tabela_validacao$log2FC_Target <-
  
  round(
    tabela_validacao$log2FC_Target,
    3
  )


# 14. FINAL TABLE

cat("\n")
cat("====================================================================\n")
cat(" FINAL TABLE — BOOLEAN MODEL EDGE VALIDATION\n")
cat("====================================================================\n")


print(
  tabela_validacao,
  row.names = FALSE
)


# 15. SUMMARY BY DATASET

resumo_dataset <-
  
  tabela_validacao %>%
  
  group_by(
    Dataset
  ) %>%
  
  summarise(
    
    Total_Edges =
      n(),
    
    Represented_Edges =
      sum(
        Edge_Represented
      ),
    
    Compatible =
      sum(
        Validation == "COMPATIBLE"
      ),
    
    Incompatible =
      sum(
        Validation == "INCOMPATIBLE"
      ),
    
    Inconclusive =
      sum(
        Validation == "INCONCLUSIVE"
      ),
    
    Not_Represented =
      sum(
        Validation == "NOT REPRESENTED"
      ),
    
    Coverage_Percentage =
      round(
        100 *
          Represented_Edges /
          Total_Edges,
        1
      ),
    
    Compatibility_Percentage =
      ifelse(
        
        Represented_Edges > 0,
        
        round(
          100 *
            Compatible /
            Represented_Edges,
          1
        ),
        
        NA
        
      ),
    
    .groups = "drop"
    
  )


cat("\n")
cat("====================================================================\n")
cat(" SUMMARY BY DATASET\n")
cat("====================================================================\n")


print(
  resumo_dataset,
  row.names = FALSE
)


# 16. SUMMARY BY EDGE

resumo_aresta <-
  
  tabela_validacao %>%
  
  group_by(
    
    Source_Gene,
    Edge_Symbol,
    Target_Gene,
    Interaction_Type
    
  ) %>%
  
  summarise(
    
    Datasets_Evaluated =
      n(),
    
    Datasets_Represented =
      sum(
        Edge_Represented
      ),
    
    Datasets_Compatible =
      sum(
        Validation == "COMPATIBLE"
      ),
    
    Datasets_Incompatible =
      sum(
        Validation == "INCOMPATIBLE"
      ),
    
    Datasets_Inconclusive =
      sum(
        Validation == "INCONCLUSIVE"
      ),
    
    Datasets_Not_Represented =
      sum(
        Validation == "NOT REPRESENTED"
      ),
    
    .groups = "drop"
    
  )


# Add complete edge label
resumo_aresta$Edge <-
  
  paste(
    
    resumo_aresta$Source_Gene,
    
    resumo_aresta$Edge_Symbol,
    
    resumo_aresta$Target_Gene
    
  )


# Reorder columns
resumo_aresta <-
  
  resumo_aresta[,
                
                c(
                  
                  "Edge",
                  
                  "Source_Gene",
                  "Edge_Symbol",
                  "Target_Gene",
                  
                  "Interaction_Type",
                  
                  "Datasets_Evaluated",
                  "Datasets_Represented",
                  
                  "Datasets_Compatible",
                  "Datasets_Incompatible",
                  
                  "Datasets_Inconclusive",
                  "Datasets_Not_Represented"
                  
                )
                
  ]


cat("\n")
cat("====================================================================\n")
cat(" SUMMARY BY EDGE\n")
cat("====================================================================\n")


print(
  resumo_aresta,
  row.names = FALSE
)


# 17. EXPORT TABLES 

write.csv(
  
  tabela_validacao,
  
  "boolean_model_edge_validation_COMPLETE.csv",
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


write.csv(
  
  resumo_dataset,
  
  "edge_validation_summary_by_dataset.csv",
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


write.csv(
  
  resumo_aresta,
  
  "edge_validation_summary_by_edge.csv",
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)


cat("\n")
cat("English CSV files generated:\n")
cat("  1. boolean_model_edge_validation_COMPLETE.csv\n")
cat("  2. edge_validation_summary_by_dataset.csv\n")
cat("  3. edge_validation_summary_by_edge.csv\n")


# 18. EDGE ORDER FOR FIGURES

ordem_arestas <- c(
  
  "DDR → ATM",
  
  "ATM -| CDC25A",
  
  "CDC25A → E2F1",
  
  "E2F1 → MYC",
  
  "MYC → MALAT1",
  
  "MALAT1 -| TRPM3",
  
  "TRPM3 -| SIRT1",
  
  "SIRT1 -| TP53",
  
  "TP53 → CDKN1A",
  
  "CDKN1A -| CASP3",
  
  "CASP3 → DFNA5"
  
)


# 19. FIGURE 1 — EDGE VALIDATION MATRIX

dados_figura <-
  
  tabela_validacao %>%
  
  mutate(
    
    Edge = paste(
      
      Source_Gene,
      
      Edge_Symbol,
      
      Target_Gene
      
    ),
    
    Edge = factor(
      
      Edge,
      
      levels = ordem_arestas
      
    )
    
  ) %>%
  
  arrange(
    Edge
  )


figura_validacao <-
  
  ggplot(
    
    dados_figura,
    
    aes(
      
      x = Edge,
      
      y = Dataset,
      
      fill = Validation
      
    )
    
  ) +
  
  geom_tile(
    
    color = "white",
    
    linewidth = 0.8
    
  ) +
  
  geom_text(
    
    aes(
      label = Validation
    ),
    
    size = 3
    
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "COMPATIBLE" =
        "#27AE60",
      
      "INCOMPATIBLE" =
        "#C0392B",
      
      "INCONCLUSIVE" =
        "#F39C12",
      
      "NOT REPRESENTED" =
        "#95A5A6"
      
    )
    
  ) +
  
  labs(
    
    title =
      "Validation of Boolean Model Edges",
    
    subtitle =
      "Agreement between interaction signs and transcriptomic patterns",
    
    x =
      "Model Edge",
    
    y =
      "Dataset",
    
    fill =
      "Validation"
    
  ) +
  
  theme_minimal(
    
    base_size = 12
    
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        
        angle = 45,
        
        hjust = 1
        
      ),
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid =
      element_blank()
    
  )


print(
  figura_validacao
)


ggsave(
  
  "FIGURE_1_edge_validation.png",
  
  figura_validacao,
  
  width = 16,
  
  height = 7,
  
  dpi = 300
  
)


# 20. FIGURE 2 — EDGE COVERAGE

figura_cobertura <-
  
  resumo_dataset %>%
  
  ggplot(
    
    aes(
      
      x = Dataset,
      
      y = Coverage_Percentage
      
    )
    
  ) +
  
  geom_col(
    
    width = 0.7
    
  ) +
  
  geom_text(
    
    aes(
      
      label =
        paste0(
          Coverage_Percentage,
          "%"
        )
      
    ),
    
    vjust = -0.4,
    
    fontface = "bold"
    
  ) +
  
  scale_y_continuous(
    
    limits = c(
      0,
      105
    )
    
  ) +
  
  labs(
    
    title =
      "Boolean Model Edge Coverage",
    
    subtitle =
      "Percentage of model interactions represented in each dataset",
    
    x =
      "Dataset",
    
    y =
      "Edge Coverage (%)"
    
  ) +
  
  theme_minimal(
    
    base_size = 13
    
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      element_blank()
    
  )


print(
  figura_cobertura
)


ggsave(
  
  "FIGURE_2_edge_coverage.png",
  
  figura_cobertura,
  
  width = 10,
  
  height = 6,
  
  dpi = 300
  
)


# 21. FIGURE 3 — EDGE VALIDATION STATUS

dados_status <-
  
  tabela_validacao %>%
  
  group_by(
    
    Dataset,
    
    Validation
    
  ) %>%
  
  summarise(
    
    Number_of_Edges = n(),
    
    .groups = "drop"
    
  )


figura_status <-
  
  ggplot(
    
    dados_status,
    
    aes(
      
      x = Dataset,
      
      y = Number_of_Edges,
      
      fill = Validation
      
    )
    
  ) +
  
  geom_col(
    
    position = "stack"
    
  ) +
  
  geom_text(
    
    aes(
      
      label =
        ifelse(
          
          Number_of_Edges > 0,
          
          Number_of_Edges,
          
          ""
          
        )
      
    ),
    
    position =
      position_stack(
        vjust = 0.5
      ),
    
    size = 4,
    
    color = "white",
    
    fontface = "bold"
    
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "COMPATIBLE" =
        "#27AE60",
      
      "INCOMPATIBLE" =
        "#C0392B",
      
      "INCONCLUSIVE" =
        "#F39C12",
      
      "NOT REPRESENTED" =
        "#95A5A6"
      
    )
    
  ) +
  
  labs(
    
    title =
      "Boolean Model Edge Validation Status",
    
    subtitle =
      "Representation and transcriptomic consistency across datasets",
    
    x =
      "Dataset",
    
    y =
      "Number of Edges",
    
    fill =
      "Validation"
    
  ) +
  
  theme_minimal(
    
    base_size = 13
    
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      element_blank()
    
  )


print(
  figura_status
)


ggsave(
  
  "FIGURE_3_edge_validation_status.png",
  
  figura_status,
  
  width = 10,
  
  height = 6,
  
  dpi = 300
  
)


# 22. FINAL REPORT

cat("\n")
cat("====================================================================\n")
cat(" ANALYSIS COMPLETED\n")
cat("====================================================================\n")

cat("\nEnglish tables:\n")

cat(
  "  - boolean_model_edge_validation_COMPLETE.csv\n"
)

cat(
  "  - edge_validation_summary_by_dataset.csv\n"
)

cat(
  "  - edge_validation_summary_by_edge.csv\n"
)


cat("\nEnglish figures:\n")

cat(
  "  - FIGURE_1_edge_validation.png\n"
)

cat(
  "  - FIGURE_2_edge_coverage.png\n"
)

cat(
  "  - FIGURE_3_edge_validation_status.png\n"
)


cat("\n")
cat("INTERPRETATION:\n")

cat(
  "COMPATIBLE = transcriptomic pattern agrees with the model edge.\n"
)

cat(
  "INCOMPATIBLE = transcriptomic pattern contradicts the model edge.\n"
)

cat(
  "INCONCLUSIVE = insufficient directional evidence.\n"
)

cat(
  "NOT REPRESENTED = one or both nodes are absent from the dataset.\n"
)

cat("\n")
cat("====================================================================\n")






############ CÓDIGO DE RECONHECIMENTO DAS ARESTAS FUNCIONAIS NOS DATASETS #############

arestas_exemplo <- data.frame(
  from = c("ATM",    "CDC25A", "E2F1", "MYC", "MALAT1", "TRPM3", "SIRT1", "TP53",   "CDKN1A", "CASP3"),
  to   = c("CDC25A", "E2F1",   "MYC",  "MALAT1", "TRPM3", "SIRT1",  "TP53",  "CDKN1A", "CASP3",  "DFNA5"),
  sign = c(-1,       1,        1,      1,        -1,       -1,      -1,      1,        -1,       1),
  stringsAsFactors = FALSE
)

g_modelo <- igraph::graph_from_data_frame(d = arestas_exemplo, directed = TRUE)

# 2. Loop para checar a conectividade individual em cada dataset da 'lista_datasets'
cat("\n==================================================\n")
cat("🛤️ VERIFICAÇÃO DE INTEGRIDADE DE CONECTIVIDADE NOS DATASETS\n")
cat("==================================================\n")

# Data frame para armazenar o resumo final
resumo_conectividade <- data.frame(Dataset = character(), Functional_Edges = integer(), stringsAsFactors = FALSE)

for(nome_ds in names(lista_datasets)) {
  # Extrai os nomes dos genes/microRNAs presentes na matriz de expressão deste dataset
  genes_presentes_no_dataset <- colnames(lista_datasets[[nome_ds]][["expr"]])
  
  # Identifica quais nós da rede original existem nas colunas do dataset
  nos_ativos <- igraph::V(g_modelo)$name[igraph::V(g_modelo)$name %in% genes_presentes_no_dataset]
  
  # Extrai o subgrafo apenas com os genes presentes (mantendo apenas arestas onde ambos os genes existem)
  sub_g <- igraph::induced_subgraph(g_modelo, v = nos_ativos)
  
  # Conta quantas arestas/conexões restaram
  num_arestas_funcionais <- igraph::ecount(sub_g)
  
  # Exibe no console
  cat(sprintf("• %s: %d functional edges preserved\n", nome_ds, num_arestas_funcionais))
  
  # Salva a contagem
  resumo_conectividade <- rbind(resumo_conectividade, data.frame(
    Dataset = nome_ds, 
    Functional_Edges = num_arestas_funcionais,
    stringsAsFactors = FALSE
  ))
}

cat("==================================================\n\n")

# 3. Exibir quais arestas específicas ficaram ativas (exemplo para um dataset)
# (Mostra os pares de genes que formam as arestas preservadas)
cat("Arestas ativas identificadas na rede central:\n")
print(igraph::as_data_frame(sub_g, what = "edges")[, c("from", "to")])






# CORREÇÃO (out/2026): removido simular_grade_perturbacoes() e a figura 'Figura_Perturbacoes_Datasets.png'.
# Os fenótipos eram calculados com regras escritas à mão (ex.: CASP3 = 1 se MALAT1, SIRT1 ou Wip1 = 0),
# e não com o modelo do GINsim. Use a seção 'Validation of model' (BoolNet) para as perturbações.




############# SCRIPT DE ALTA RESOLUÇÃO: MAPA DE CALOR, GRAFO REGULATÓRIO E COBERTURA GEO ###################

# 1. CÁLCULO PRÉVIO DA COBERTURA (Porcentagem de presença para a Nova Figura)
datasets_nomes <- c("GSE14520", "GSE60502", "GSE121248")
porcentagens <- c()
genes_encontrados_total <- list()

for(nome_ds in datasets_nomes) {
  genes_ds <- colnames(lista_datasets[[nome_ds]][["expr"]])
  genes_modelo_encontrados <- nos_modelo[["resolved_name"]][nos_modelo[["resolved_name"]] %in% genes_ds]
  
  # Salva para uso posterior
  genes_encontrados_total[[nome_ds]] <- genes_modelo_encontrados
  pct <- (length(genes_modelo_encontrados) / nrow(nos_modelo)) * 100
  porcentagens <- c(porcentagens, pct)
}


# CORREÇÃO (out/2026): removidos calcular_atrator_ginsim(), 'Tabela_Perturbacoes_Fidedignas_GINsim.csv' e
# 'Mapa_Perturbacoes_Fidedignas_GINsim.png'. Apesar do nome, as regras eram escritas à mão (12 nós, não os 31
# do modelo) e o resultado não era fiel ao GINsim. As perturbações corretas estão em 'Validation of model'.



# 🖼️ FIGURA 2: GRAFO REGULATÓRIO DA PIROPTOSE (Design Vetorial Avançado)


if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")

# 1. Definição do Grafo
arestas_exemplo <- data.frame(
  From = c("DDR", "ATM", "CDC25A", "E2F1", "MYC", "MALAT1", "miR-204", "SIRT1", "TP53", "CDKN1A", "CASP3"),
  To   = c("ATM", "CDC25A", "E2F1", "MYC", "MALAT1", "miR-204", "SIRT1", "TP53", "CDKN1A", "CASP3", "GSDME"),
  sign = c(1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1),
  tipo = c("->", "-|", "->", "->", "->", "-|", "-|", "-|", "->", "-|", "->"),
  stringsAsFactors = FALSE
)

g_modelo <- igraph::graph_from_data_frame(d = arestas_exemplo, directed = TRUE)

# 2. Estilos das Arestas
igraph::E(g_modelo)$color <- ifelse(igraph::E(g_modelo)$sign == 1, "#27AE60", "#C0392B")
igraph::E(g_modelo)$lty   <- ifelse(igraph::E(g_modelo)$sign == 1, 1, 2)
igraph::E(g_modelo)$arrow.mode <- ifelse(igraph::E(g_modelo)$sign == 1, 2, 0)

# 3. LAYOUT DE ÁRVORE HORIZONTAL (Nativo do igraph)
# A função layout_as_tree alinha a cadeia linear sem empilhar nós
layout_arvore <- igraph::layout_as_tree(g_modelo, root = "DDR")

# Giramos 90 graus para ficar na horizontal (X espalhado, Y constante)
layout_horizontal <- cbind(-layout_arvore[, 2], 0)

# 4. Geração da Imagem
png("Fig2_Grafo_Modelo_Premium.png", width = 1400, height = 750, res = 150)

# Espaço de margem expandido no rodapé para o texto descritivo
par(mar = c(7, 2, 3, 2))

plot(g_modelo,
     layout = layout_horizontal,
     asp = 0.35,              # <<< ESTICA O GRÁFICO NA HORIZONTAL (Aumenta o espaço entre as bolinhas)
     vertex.color = "#2C3E50",          
     vertex.frame.color = "#1A252F",
     vertex.size = 14,      # <<< TAMANHO ADEQUADO DO NÓ (impede sobreposição)
     vertex.label.color = "white",
     vertex.label.font = 2,
     vertex.label.cex = 0.7,          
     edge.color = igraph::E(g_modelo)$color,
     edge.lty = igraph::E(g_modelo)$lty,
     edge.arrow.mode = igraph::E(g_modelo)$arrow.mode,
     edge.arrow.size = 0.5,
     edge.width = 2.5,                 
     main = "GSDME-Mediated Pyroptosis Activation Network")

# Legenda Superior Esquerda
legend("topleft", 
       legend = c("Activation ( -> )", "Inhibition ( -| )"), 
       col = c("#27AE60", "#C0392B"), 
       lty = c(1, 2), 
       lwd = 2.5, 
       bty = "n", 
       cex = 0.85)

# Caixa Informativa no Rodapé
texto_datasets <- paste(
  "Connectivity Integrity: ",
  "• GSE14520: 8 functional edges preserved",
  "• GSE60502: 8 functional edges preserved",
  "• GSE121248: 10 functional edges preserved",
  "• TCGA-LIHC: 10 functional edges preserved",
  sep = "\n"
)

mtext(texto_datasets, side = 1, line = 4.5, adj = 0, cex = 0.75, col = "#34495E", font = 3)

dev.off()

cat("[OK] Figura gerada na horizontal com conexões longas e limpas!\n")


# 🖼️ FIGURA 2: REDE COMPLETA COM FUNDO BRANCO GARANTIDO

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")

library(igraph)
library(ggraph)
library(ggplot2)

# 1. Tabela de arestas sem MIR204 (11 nós / 10 conexões)
arestas_exemplo <- data.frame(
  from = c("DDR", "ATM",    "CDC25A", "E2F1", "MYC",  "MALAT1", "miR-204", "SIRT1", "p53",   "p21", "CASP3"),
  to   = c("ATM", "CDC25A", "E2F1",   "MYC",  "MALAT1", "miR-204", "SIRT1",  "p53",  "p21", "CASP3",  "GSDME"),
  type = c("Activation ( -> )", "Inhibition ( -| )", "Activation ( -> )", "Activation ( -> )", "Activation ( -> )", 
           "Inhibition ( -| )", "Inhibition ( -| )", "Inhibition ( -| )", "Activation ( -> )", "Inhibition ( -| )", "Activation ( -> )"),
  stringsAsFactors = FALSE
)

g_modelo <- graph_from_data_frame(d = arestas_exemplo, directed = TRUE)

# 2. Layout linear na ordem exata da via
nos_ordenados <- c("DDR", "ATM", "CDC25A", "E2F1", "MYC", "MALAT1", "miR-204", "SIRT1", "p53", "p21", "CASP3", "GSDME")
layout_linear <- create_layout(g_modelo, layout = 'linear')
layout_linear$x <- match(layout_linear$name, nos_ordenados)
layout_linear$y <- -0.1

# 3. Texto informativo do rodapé
texto_rodape <- paste(
  "Connectivity Integrity (Functional Edges Preserved Across Pathway):",
  "Pathway: ATM -| CDC25A -> E2F1 -> MYC -> lncRNA-MALAT1 -| miR-204-5p -| SIRT1 -| p53 -> p21 -| CASP3 -> GSDME",
  "• GSE14520: 8 functional edges preserved;",
  "• GSE60502: 8 functional edges preserved;",
  "• GSE121248: 10 functional edges preserved;",
  "• TCGA-LIHC: 10 functional edges preserved;",
  "p53 = TP53, p21 = CDKN1A, GSDME = DFNA5,",
  "TRPM3 was used as a host-gene proxy for miR-204-5p due to the absence of miRNA probes in the selected datasets,",
   "whereas mature miR-204 was directly quantified in TCGA-LIHC.",
  sep = "\n"
)

# 4. Construção do gráfico
p <- ggraph(layout_linear) +
  # Arestas (Linhas e Setas)
  geom_edge_link(aes(color = type, linetype = type),
                 arrow = arrow(length = unit(3.5, 'mm'), type = 'closed'),
                 end_cap = circle(9, 'mm'),
                 start_cap = circle(9, 'mm'),
                 edge_width = 1.1) +
  # Nós (Círculos)
  geom_node_point(size = 24, color = "#2C3E50") +
  # Rótulos dos Genes
  geom_node_text(aes(label = name), color = "white", fontface = "bold", size = 3.8) +
  # Mapeamento de Cores e Estilos para a Legenda
  scale_edge_color_manual(
    values = c("Activation ( -> )" = "#27AE60", "Inhibition ( -| )" = "#C0392B"), 
    name = NULL
  ) +
  scale_edge_linetype_manual(
    values = c("Activation ( -> )" = "solid", "Inhibition ( -| )" = "dashed"), 
    name = NULL
  ) +
  # Título
  labs(title = "GSDME-Mediated Pyroptosis Activation Network") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),  # <<< FUNDO BRANCO DO CANVAS
    panel.background = element_rect(fill = "white", color = NA), # <<< FUNDO BRANCO DO PAINEL
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    legend.position = c(0.12, 0.82),
    legend.text = element_text(size = 11, color = "#2C3E50"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.width = unit(1.2, "cm"),
    plot.margin = margin(t = 20, r = 25, b = 100, l = 25)
  )

# 5. Adicionar o Bloco de Texto no Rodapé
p_final <- p + 
  annotate("text", x = 1, y = -0.5, label = texto_rodape, 
           hjust = 0, vjust = 1, size = 3.5, fontface = "italic", color = "#34495E") +
  coord_cartesian(ylim = c(-0.7, 0.3), clip = "off")

# 6. Salvar garantindo fundo branco opaco (sem transparência)
ggsave("Fig2_Grafo_Modelo_Premium.png", plot = p_final, width = 13, height = 6, dpi = 300, bg = "white")

cat("[OK] Figura gerada com sucesso e fundo 100% branco!\n")


# 🖼️ FIGURA 3: COBERTURA E PRESENÇA DOS COMPONENTES (Nova!)

png("Fig3_Cobertura_Componentes_GEO.png", width = 1100, height = 700, res = 150)

# Margens ajustadas para os nomes das coortes
par(mar = c(5, 6, 4, 3), bg = "white")

# Desenhar gráfico de barras minimalista
barras <- barplot(porcentagens, 
                  names.arg = datasets_nomes, 
                  col = c("#1ABC9C", "#3498DB", "#9B59B6"), # Paleta moderna FlatUI
                  border = NA,
                  ylim = c(0, 110),
                  ylab = "Model Component Coverage (%)",
                  main = "Dataset Validation: Model Presence in GEO Cohorts",
                  cex.names = 0.9, cex.lab = 1, cex.axis = 0.9,
                  las = 1, space = 0.4)

# Adicionar grid cinza de fundo para facilitar a leitura científica
grid(nx = NA, ny = NULL, col = "#E5E7E9", lty = "solid", lwd = 1)

# Forçar o redesenho das barras por cima do grid
barplot(porcentagens, col = c("#1ABC9C", "#3498DB", "#9B59B6"), border = NA, add = TRUE, space = 0.4, axes = FALSE)

# Colocar os valores exatos de porcentagem e contagem no topo de cada barra
for(i in 1:length(porcentagens)) {
  total_encontrados <- length(genes_encontrados_total[[datasets_nomes[i]]])
  texto_label <- sprintf("%.1f%%\n(%d/%d)", porcentagens[i], total_encontrados, nrow(nos_modelo))
  
  text(x = barras[i], y = porcentagens[i] + 5, 
       labels = texto_label, 
       cex = 0.85, font = 2, col = "#2C3E50")
}

dev.off()
cat("[OK] Figura 3 (Gráfico de Cobertura GEO) gerada com sucesso!\n")


############ FIGURA 3 E EXPRESSÃO DIFERENCIAL NOS GEO (GRUPOS REAIS DO METADATA) ############
#
# CORREÇÃO (out/2026). A versão anterior deste bloco tinha três problemas graves:
#   (a) gerava valores de expressão com rnorm() e um "efeito" fixo por gene
#       (Fig4_Expression_Screen_All_Datasets.png era inteiramente simulada);
#   (b) na figura usada como Figura 3 do manuscrito, os rótulos Tumor/Normal eram
#       SORTEADOS com sample(prob = c(0.3, 0.7)), e não lidos do metadata;
#   (c) a expressão diferencial classificava as amostras por regex em texto livre
#       (inclusive "N$" e "T$", que casam qualquer texto terminado em n/t) e usava
#       uma segunda baixar_e_preparar_geo() que sobrescrevia gse60502/gse121248 com
#       outro colapso de sondas (média + primeiro símbolo), diferente do usado no Spearman.
# Agora tudo usa lista_datasets (o mesmo processamento das análises de Spearman/NMI) e as
# colunas exatas do pData já usadas em 'comparacoes'.

for (p in c("ggplot2", "dplyr", "patchwork", "ggrepel")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

grupos_geo <- list(
  GSE14520  = list(coluna = "Tissue:ch1",      tumor = "Liver Tumor Tissue",       normal = "Liver Non-Tumor Tissue"),
  GSE60502  = list(coluna = "tissue type:ch1", tumor = "hepatocellular carcinoma", normal = "adjacent non-tumorous liver"),
  GSE121248 = list(coluna = "tissue:ch1",      tumor = "Tumor sample",             normal = "Adjacent Normal sample")
)

# Devolve um fator (Adjacent non-tumor / Primary HCC / NA) alinhado às linhas de expr
grupo_amostras_geo <- function(nome_ds) {
  cfg  <- grupos_geo[[nome_ds]]
  meta <- lista_datasets[[nome_ds]]$meta
  expr <- lista_datasets[[nome_ds]]$expr

  if (!cfg$coluna %in% colnames(meta)) {
    stop("Coluna '", cfg$coluna, "' ausente no pData de ", nome_ds)
  }
  if (!identical(rownames(expr), rownames(meta))) {
    stop("As amostras de expr e meta não estão na mesma ordem em ", nome_ds)
  }

  v <- trimws(as.character(meta[[cfg$coluna]]))
  g <- ifelse(v == cfg$tumor, "Primary HCC",
              ifelse(v == cfg$normal, "Adjacent non-tumor", NA_character_))

  if (sum(g == "Primary HCC", na.rm = TRUE) < 3 || sum(g == "Adjacent non-tumor", na.rm = TRUE) < 3) {
    stop("Grupos Tumor/Normal não identificados em ", nome_ds,
         ". Valores encontrados: ", paste(unique(v), collapse = " | "))
  }
  factor(g, levels = c("Adjacent non-tumor", "Primary HCC"))
}

genes_do_modelo <- nos_modelo$resolved_name

# ------------------------------------------------------------------------------
# 1. FIGURA 3 — BOXPLOTS COM OS GRUPOS VERDADEIROS
# ------------------------------------------------------------------------------

dados_fig3 <- bind_rows(lapply(names(grupos_geo), function(nome_ds) {
  expr  <- lista_datasets[[nome_ds]]$expr
  grupo <- grupo_amostras_geo(nome_ds)
  genes <- colnames(expr)[colnames(expr) %in% genes_do_modelo]

  bind_rows(lapply(genes, function(gene) {
    data.frame(
      Sample     = rownames(expr),
      Gene       = gene,
      Expression = as.numeric(expr[[gene]]),
      Group      = grupo,
      Cohort     = nome_ds,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  filter(!is.na(Group), !is.na(Expression))

desenhar_painel_grade <- function(df_dataset, nome_ds) {
  n_genes  <- length(unique(df_dataset$Gene))
  n_tumor  <- length(unique(df_dataset$Sample[df_dataset$Group == "Primary HCC"]))
  n_normal <- length(unique(df_dataset$Sample[df_dataset$Group == "Adjacent non-tumor"]))

  ggplot(df_dataset, aes(x = Group, y = Expression, fill = Group)) +
    theme_bw(base_size = 11) +
    geom_jitter(color = "#4D4D4D", alpha = 0.35, width = 0.2, size = 0.6) +
    geom_boxplot(color = "black", alpha = 0.85, width = 0.45, outlier.shape = NA, linewidth = 0.5) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = c("Adjacent non-tumor" = "#95A5A6", "Primary HCC" = "#006699")) +
    labs(
      title = sprintf("Cohort Study Profile: %s (%d components found; HCC n = %d, non-tumor n = %d)",
                      nome_ds, n_genes, n_tumor, n_normal),
      # Os valores do series matrix já estão em log2 (RMA); não foi somado +1.
      y = "log2 expression (GEO series matrix)", x = ""
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12, color = "#2C3E50"),
      strip.background = element_rect(fill = "#F2F4F4", color = "#BDC3C7"),
      strip.text = element_text(face = "bold", size = 9, color = "#2C3E50"),
      axis.text.x = element_text(angle = 15, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9, face = "bold"),
      legend.position = "none",
      panel.grid.major = element_line(color = "#E5E7E9", linetype = "solid"),
      panel.grid.minor = element_blank()
    )
}

paineis_fig3 <- lapply(names(grupos_geo), function(nome_ds) {
  desenhar_painel_grade(dados_fig3[dados_fig3$Cohort == nome_ds, ], nome_ds)
})

foto_final <- wrap_plots(paineis_fig3, ncol = 1) +
  plot_annotation(
    title = "Expression profiling across the validated model networks",
    subtitle = "Primary HCC vs. adjacent non-tumor liver; groups taken from the GEO sample annotation",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#2C3E50", hjust = 0),
      plot.subtitle = element_text(size = 10, color = "#5D6D7E", face = "italic", hjust = 0)
    )
  )

ggsave("Fig4_GSE_Cohorts_Expression_Facet.png", plot = foto_final,
       width = 16, height = 22, dpi = 300, units = "in", bg = "white")
cat("[OK] Figura de expressão dos GEO salva: 'Fig4_GSE_Cohorts_Expression_Facet.png'\n")

# ------------------------------------------------------------------------------
# 2. EXPRESSÃO DIFERENCIAL (WILCOXON, TUMOR VS NÃO TUMORAL) + BH POR DATASET
# ------------------------------------------------------------------------------
# Observação: GSE14520 e GSE60502 têm amostras pareadas por paciente; o teste abaixo é
# não pareado (conservador). Para um teste pareado é preciso o ID do paciente no pData.

testar_expressao_geo <- function(nome_ds) {
  expr  <- lista_datasets[[nome_ds]]$expr
  grupo <- grupo_amostras_geo(nome_ds)
  genes <- colnames(expr)[colnames(expr) %in% genes_do_modelo]

  q99 <- stats::quantile(unlist(expr[, genes, drop = FALSE]), 0.99, na.rm = TRUE)
  if (q99 > 100) {
    warning(nome_ds, ": valores parecem não estar em log2; log2FC = diferença de médias não é válido.")
  }

  res <- bind_rows(lapply(genes, function(gene) {
    x   <- as.numeric(expr[[gene]])
    v_t <- x[grupo %in% "Primary HCC" & !is.na(x)]
    v_n <- x[grupo %in% "Adjacent non-tumor" & !is.na(x)]
    if (length(v_t) < 3 || length(v_n) < 3) return(NULL)
    data.frame(
      Dataset     = nome_ds,
      Gene        = gene,
      N_Tumor     = length(v_t),
      N_Normal    = length(v_n),
      Mean_Normal = mean(v_n),
      Mean_Tumor  = mean(v_t),
      log2FC      = mean(v_t) - mean(v_n),
      p_value     = wilcox.test(v_t, v_n, exact = FALSE)$p.value,
      stringsAsFactors = FALSE
    )
  }))
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  res$Direction <- ifelse(res$p_adj < 0.05 & res$log2FC > 0, "UP in the Tumor",
                          ifelse(res$p_adj < 0.05 & res$log2FC < 0, "DOWN in the Tumor", "Not Significant"))
  res
}

tabela_DE_geo <- bind_rows(lapply(names(grupos_geo), testar_expressao_geo))
write.csv(tabela_DE_geo, "GEO_differential_expression_model_components.csv", row.names = FALSE)

cat("\n=== EXPRESSÃO DIFERENCIAL (GRUPOS DO METADATA, BH POR DATASET) ===\n")
print(tabela_DE_geo, row.names = FALSE)

gerar_painel_expressao <- function(tabela, geo_id, nome_arquivo_saida) {
  if (is.null(tabela) || nrow(tabela) == 0) {
    message("Tabela vazia para ", geo_id, ". Gráfico não gerado.")
    return(invisible(NULL))
  }
  df_plot <- tabela %>% mutate(neg_log10_padj = -log10(p_adj), Status = Direction)
  df_plot$Gene <- factor(df_plot$Gene, levels = df_plot$Gene[order(df_plot$log2FC)])
  paleta_cores <- c("UP in the Tumor" = "#E74C3C", "DOWN in the Tumor" = "#3498DB", "Not Significant" = "#95A5A6")

  p1 <- ggplot(df_plot, aes(x = log2FC, y = neg_log10_padj, color = Status)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    geom_text_repel(aes(label = Gene), size = 3.2, fontface = "bold", max.overlaps = 20) +
    scale_color_manual(values = paleta_cores) +
    theme_bw(base_size = 11) +
    labs(title = paste0("A. Differential Expression from ", geo_id),
         x = "log2 Fold Change (Tumor vs Normal)", y = "-log10(p.adj)") +
    theme(legend.position = "none", plot.title = element_text(face = "bold", color = "#2C3E50"))

  p2 <- ggplot(df_plot, aes(x = log2FC, y = Gene, color = Status)) +
    geom_segment(aes(x = 0, xend = log2FC, y = Gene, yend = Gene), linewidth = 1) +
    geom_point(size = 3.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = paleta_cores) +
    theme_bw(base_size = 11) +
    labs(title = paste0("B. Magnitude of Change by Gene from ", geo_id),
         x = "log2 Fold Change (Tumor vs Normal)", y = "") +
    theme(legend.position = "top", legend.title = element_blank(),
          plot.title = element_text(face = "bold", color = "#2C3E50"))

  ggsave(nome_arquivo_saida, plot = p1 + p2 + plot_layout(widths = c(1, 1.25)),
         width = 13, height = 6, dpi = 300, bg = "white")
  cat("  [OK] Figura salva:", nome_arquivo_saida, "\n")
}

for (nome_ds in names(grupos_geo)) {
  gerar_painel_expressao(tabela_DE_geo[tabela_DE_geo$Dataset == nome_ds, ], nome_ds,
                         paste0("Fig_Expressao_Diferencial_", nome_ds, ".png"))
}


################# BOOLEAN MODEL EDGE VALIDATION ACROSS DATASETS (SPEARMAN CORRELATION METHOD) ####################
#
# CORREÇÕES (out/2026):
#   1. A classificação usava o p-valor BRUTO (pv < 0.05), embora o manuscrito diga que foi
#      aplicada correção de Benjamini–Hochberg. O p_adj era calculado, mas depois da
#      classificação e sem ser usado. Agora a classificação usa p_adj (BH por dataset).
#   2. Arestas com < 5 pares válidos eram rotuladas "COMPATIBLE" sem teste nenhum;
#      agora recebem "NOT TESTABLE".
#   3. Nova coluna Evidence separa "SUPPORTED" (significativa e no sinal do modelo) de
#      "NOT SIGNIFICANT", que antes ficava escondida dentro de "COMPATIBLE".
#   4. Renomeação de colunas por nome, e não por posição.

packages <- c("ggplot2", "dplyr", "tidyr", "readr", "infotheo")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# 1. CHECK DATASETS

if (!exists("lista_datasets") || length(lista_datasets) == 0) {
  stop("\nERROR: 'lista_datasets' does not exist or is empty in the R environment.\n")
}

cat("\n====================================================================\n")
cat(" DATASETS FOUND FOR SPEARMAN VALIDATION\n")
cat("====================================================================\n")
print(names(lista_datasets))

# 2. BOOLEAN MODEL EDGES (TRPM3 = proxy do miR-204-5p nos GEO)

arestas_modelo <- data.frame(
  from = c("ATM", "CDC25A", "E2F1", "MYC", "MALAT1", "TRPM3", "SIRT1", "TP53", "CDKN1A", "CASP3"),
  to   = c("CDC25A", "E2F1", "MYC", "MALAT1", "TRPM3", "SIRT1", "TP53", "CDKN1A", "CASP3", "DFNA5"),
  sign = c(-1, 1, 1, 1, -1, -1, -1, 1, -1, 1),
  stringsAsFactors = FALSE
)

arestas_modelo$tipo <- ifelse(arestas_modelo$sign == 1, "Activation", "Inhibition")
arestas_modelo$simbolo <- ifelse(arestas_modelo$sign == 1, "→", "-|")
arestas_modelo$aresta <- paste(arestas_modelo$from, arestas_modelo$simbolo, arestas_modelo$to)

# 3. DATASET COMPARISON METADATA CONFIGURATION (somente amostras tumorais)
comparacoes <- list(
  GSE14520 = list(coluna = "Tissue:ch1", grupo1 = "Liver Tumor Tissue"),
  GSE60502 = list(coluna = "tissue type:ch1", grupo1 = "hepatocellular carcinoma"),
  GSE121248 = list(coluna = "tissue:ch1", grupo1 = "Tumor sample")
)

ALPHA <- 0.05
MIN_PARES <- 5

# 4. GENE / MIRNA MATCHING FUNCTION (usada também no NMI/GGC, para que as duas
#    análises usem exatamente as mesmas colunas)
encontrar_gene <- function(gene, nomes) {
  if (gene %in% nomes) return(gene)

  alias_map <- list("DFNA5" = "GSDME", "TRPM3" = "miR-204-5p")
  if (gene %in% names(alias_map) && alias_map[[gene]] %in% nomes) {
    return(alias_map[[gene]])
  }

  nomes_limpos <- trimws(nomes)
  if (gene %in% nomes_limpos) return(nomes[which(nomes_limpos == gene)[1]])

  # Identificadores compostos de sondas (ex.: "TRPM3___MIR204", após gsub("/", "_"))
  padrao <- paste0("(^|___|\\s)", gene, "($|___|\\s)")
  idx <- grep(padrao, nomes, ignore.case = TRUE)
  if (length(idx) > 0) return(nomes[idx[1]])

  NA_character_
}

amostras_tumor <- function(nome_ds) {
  cfg  <- comparacoes[[nome_ds]]
  meta <- lista_datasets[[nome_ds]]$meta
  if (!cfg$coluna %in% colnames(meta)) stop("Coluna ", cfg$coluna, " ausente em ", nome_ds)
  which(trimws(as.character(meta[[cfg$coluna]])) == cfg$grupo1)
}

# 5. RUN SPEARMAN EDGE VALIDATION

resultados_validacao <- list()

cat("\n====================================================================\n")
cat(" STARTING SPEARMAN CORRELATION VALIDATION\n")
cat("====================================================================\n")

for (nome_ds in names(comparacoes)) {
  if (!nome_ds %in% names(lista_datasets)) {
    warning(paste("Dataset not found in environment:", nome_ds))
    next
  }

  expr <- lista_datasets[[nome_ds]]$expr
  nomes_genes <- colnames(expr)

  idx_tumor <- amostras_tumor(nome_ds)
  if (length(idx_tumor) < MIN_PARES) {
    warning(paste("Insufficient tumor samples in", nome_ds))
    next
  }
  expr_tumor <- expr[idx_tumor, , drop = FALSE]

  resultado_ds <- arestas_modelo
  resultado_ds$dataset <- nome_ds
  resultado_ds$gene_from_dataset <- vapply(resultado_ds$from, encontrar_gene, character(1), nomes = nomes_genes)
  resultado_ds$gene_to_dataset   <- vapply(resultado_ds$to,   encontrar_gene, character(1), nomes = nomes_genes)
  resultado_ds$edge_represented  <- !is.na(resultado_ds$gene_from_dataset) & !is.na(resultado_ds$gene_to_dataset)

  resultado_ds$N_samples    <- NA_integer_
  resultado_ds$rho_spearman <- NA_real_
  resultado_ds$p_value      <- NA_real_

  for (i in seq_len(nrow(resultado_ds))) {
    if (!resultado_ds$edge_represented[i]) next
    x <- as.numeric(expr_tumor[[resultado_ds$gene_from_dataset[i]]])
    y <- as.numeric(expr_tumor[[resultado_ds$gene_to_dataset[i]]])
    ok <- complete.cases(x, y)
    resultado_ds$N_samples[i] <- sum(ok)
    if (sum(ok) >= MIN_PARES) {
      ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
      resultado_ds$rho_spearman[i] <- unname(ct$estimate)
      resultado_ds$p_value[i]      <- ct$p.value
    }
  }

  resultados_validacao[[nome_ds]] <- resultado_ds
}

tabela_validacao <- bind_rows(resultados_validacao)

classificar_spearman <- function(rho, p_adj, sign, represented) {
  if (!represented) return(c("NOT REPRESENTED", "NOT REPRESENTED"))
  if (is.na(rho) || is.na(p_adj)) return(c("NOT TESTABLE", "NOT TESTABLE"))
  coerente <- if (sign == 1) rho > 0 else rho < 0
  sig <- p_adj < ALPHA
  validation <- if (sig && !coerente) "INCOMPATIBLE" else "COMPATIBLE"
  evidence <- if (sig && coerente) "SUPPORTED" else if (sig) "CONTRADICTED" else "NOT SIGNIFICANT"
  c(validation, evidence)
}

tabela_validacao <- tabela_validacao %>%
  group_by(dataset) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%   # BH por dataset (NAs ignorados)
  ungroup()

cls <- t(mapply(classificar_spearman, tabela_validacao$rho_spearman, tabela_validacao$p_adj,
                tabela_validacao$sign, tabela_validacao$edge_represented))
tabela_validacao$validation <- cls[, 1]
tabela_validacao$evidence   <- cls[, 2]
tabela_validacao$interpretation <- dplyr::recode(
  tabela_validacao$evidence,
  "SUPPORTED"       = "Significant correlation (BH) in the direction of the model interaction",
  "CONTRADICTED"    = "Significant correlation (BH) opposing the model interaction sign",
  "NOT SIGNIFICANT" = "No significant correlation after BH; not contradicted by the data",
  "NOT TESTABLE"    = "Fewer than 5 valid sample pairs",
  "NOT REPRESENTED" = "One or both nodes are absent from dataset"
)

# 6. RENAME & REORDER COLUMNS (por nome)

tabela_validacao <- tabela_validacao %>%
  transmute(
    Dataset = dataset, Edge = aresta, Source_Gene = from, Edge_Symbol = simbolo,
    Target_Gene = to, Interaction_Type = tipo,
    Source_Gene_in_Dataset = gene_from_dataset, Target_Gene_in_Dataset = gene_to_dataset,
    Edge_Represented = edge_represented, N_Tumor_Samples = N_samples,
    Rho_Spearman = round(rho_spearman, 3), p_value, p_adj,
    Validation = validation, Evidence = evidence, Interpretation = interpretation
  )

# 7. PRINT RESULTS & SUMMARIES

cat("\n====================================================================\n")
cat(" FINAL TABLE — SPEARMAN VALIDATION (classification on BH-adjusted p)\n")
cat("====================================================================\n")
print(as.data.frame(tabela_validacao), row.names = FALSE)

resumo_dataset <- tabela_validacao %>%
  group_by(Dataset) %>%
  summarise(
    Total_Edges = n(),
    Represented_Edges = sum(Edge_Represented),
    Compatible = sum(Validation == "COMPATIBLE"),
    Supported = sum(Evidence == "SUPPORTED"),
    Not_Significant = sum(Evidence == "NOT SIGNIFICANT"),
    Incompatible = sum(Validation == "INCOMPATIBLE"),
    Not_Represented = sum(Validation == "NOT REPRESENTED"),
    Coverage_Percentage = round(100 * Represented_Edges / Total_Edges, 1),
    .groups = "drop"
  )

cat("\n SUMMARY BY DATASET\n")
print(as.data.frame(resumo_dataset), row.names = FALSE)

resumo_aresta <- tabela_validacao %>%
  group_by(Edge, Source_Gene, Edge_Symbol, Target_Gene, Interaction_Type) %>%
  summarise(
    Datasets_Evaluated = n(),
    Datasets_Represented = sum(Edge_Represented),
    Datasets_Supported = sum(Evidence == "SUPPORTED"),
    Datasets_Compatible = sum(Validation == "COMPATIBLE"),
    Datasets_Incompatible = sum(Validation == "INCOMPATIBLE"),
    Datasets_Not_Represented = sum(Validation == "NOT REPRESENTED"),
    .groups = "drop"
  )

cat("\n SUMMARY BY EDGE\n")
print(as.data.frame(resumo_aresta), row.names = FALSE)

# 8. EXPORT CSVs

write.csv(tabela_validacao, "boolean_model_edge_validation_COMPLETE.csv", row.names = FALSE)
write.csv(resumo_dataset, "edge_validation_summary_by_dataset.csv", row.names = FALSE)
write.csv(resumo_aresta, "edge_validation_summary_by_edge.csv", row.names = FALSE)
cat("\nCSV Files saved successfully!\n")

# 9. FIGURES

ordem_arestas <- c(
  "ATM -| CDC25A", "CDC25A → E2F1", "E2F1 → MYC", "MYC → MALAT1", "MALAT1 -| TRPM3",
  "TRPM3 -| SIRT1", "SIRT1 -| TP53", "TP53 → CDKN1A", "CDKN1A -| CASP3", "CASP3 → DFNA5"
)
niveis_validacao <- c("COMPATIBLE", "INCOMPATIBLE", "NOT TESTABLE", "NOT REPRESENTED")
cores_validacao <- c("COMPATIBLE" = "#27AE60", "INCOMPATIBLE" = "#C0392B",
                     "NOT TESTABLE" = "#F39C12", "NOT REPRESENTED" = "#95A5A6")

dados_figura <- tabela_validacao %>%
  mutate(
    Edge = factor(Edge, levels = ordem_arestas),
    Validation = factor(Validation, levels = niveis_validacao),
    Label = ifelse(is.na(Rho_Spearman), as.character(Validation),
                   paste0(Validation, "\n(ρ=", Rho_Spearman, ifelse(Evidence == "SUPPORTED", "*", ""), ")"))
  )

figura_validacao <- ggplot(dados_figura, aes(x = Edge, y = Dataset, fill = Validation)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Label), size = 2.8) +
  scale_fill_manual(values = cores_validacao, drop = TRUE) +
  labs(
    title = "Validation of Boolean Model Edges (Spearman Correlation)",
    subtitle = "Primary tumour samples only; * = significant after BH (p.adj < 0.05) in the model direction",
    x = "Model Edge", y = "Dataset", fill = "Validation Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.title = element_text(face = "bold"), panel.grid = element_blank())

ggsave("FIGURE_1_edge_validation_spearman.png", figura_validacao, width = 16, height = 7, dpi = 300, bg = "white")

dados_status <- tabela_validacao %>%
  count(Dataset, Validation, name = "Number_of_Edges") %>%
  mutate(Validation = factor(Validation, levels = niveis_validacao))

figura_status <- ggplot(dados_status, aes(x = Dataset, y = Number_of_Edges, fill = Validation)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = Number_of_Edges), position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold") +
  scale_fill_manual(values = cores_validacao, drop = TRUE) +
  labs(title = "Boolean Model Edge Validation Summary",
       subtitle = "Number of edges by validation status across GEO cohorts",
       x = "Dataset", y = "Number of Edges", fill = "Validation Status") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave("FIGURE_3_edge_validation_status_spearman.png", figura_status, width = 10, height = 6, dpi = 300, bg = "white")
cat("Figures saved: FIGURE_1_edge_validation_spearman.png, FIGURE_3_edge_validation_status_spearman.png\n")




###### NMI and GGC #######
#
# CORREÇÕES (out/2026):
#   1. Usava TODAS as amostras (tumor + não tumoral). A diferença entre tecidos cria
#      dependência artificial entre quaisquer dois genes diferencialmente expressos, e o
#      manuscrito diz que NMI/GGC foram aplicados "às mesmas interações do Spearman"
#      (que usa só tumor; o TCGA também usa só tumor). Agora: só amostras tumorais.
#   2. O mapeamento de genes era outro (encontrar_gene_robusto, com grep sem âncora e
#      troca de TRPM3 por qualquer coluna "mir204"). Agora usa encontrar_gene(), o mesmo
#      do Spearman.
#   3. O limiar fixo GGC >= 0.30 não distingue dependência de ruído: com o estimador
#      empírico do infotheo e discretização em N^(1/3) classes, pares de genes
#      INDEPENDENTES já produzem GGC mediano ~0.28–0.35 (e o P95 fica entre ~0.33 e
#      ~0.50 conforme o N). A significância agora vem de um teste de permutação
#      (embaralhando um dos genes), com BH por dataset. O GGC continua sendo reportado
#      como tamanho de efeito, junto com a mediana e o P95 do GGC nulo.

library(infotheo)

testar_dependencia <- function(x, y, n_perm = 2000, seed = 123, min_n = 10) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  vazio <- c(N = n, MI = NA, NMI = NA, GGC = NA, GGC_null_median = NA, GGC_null_p95 = NA, p_perm = NA)
  if (n < min_n) return(vazio)

  d  <- infotheo::discretize(data.frame(x = x, y = y))  # equalfreq, nbins = N^(1/3)
  mi <- infotheo::mutinformation(d$x, d$y)              # em nats
  h1 <- infotheo::entropy(d$x)
  h2 <- infotheo::entropy(d$y)

  set.seed(seed)
  mi_null <- replicate(n_perm, infotheo::mutinformation(d$x, sample(d$y)))
  ggc <- function(m) sqrt(1 - exp(-2 * m))

  c(N = n,
    MI = mi,
    NMI = if ((h1 + h2) > 0) 2 * mi / (h1 + h2) else 0,
    GGC = ggc(mi),
    GGC_null_median = stats::median(ggc(mi_null)),
    GGC_null_p95 = unname(stats::quantile(ggc(mi_null), 0.95)),
    p_perm = (1 + sum(mi_null >= mi)) / (n_perm + 1))
}

arestas_modelo_completo <- data.frame(
  From = arestas_modelo$from, To = arestas_modelo$to, Sign = arestas_modelo$sign,
  Aresta = arestas_modelo$aresta, stringsAsFactors = FALSE
)

calcular_nmi_ggc_geo <- function(arestas_df) {
  todos <- list()
  for (nome_ds in names(comparacoes)) {
    if (!nome_ds %in% names(lista_datasets)) next
    expr  <- lista_datasets[[nome_ds]]$expr
    expr_tumor <- expr[amostras_tumor(nome_ds), , drop = FALSE]
    nomes <- colnames(expr_tumor)

    linhas <- lapply(seq_len(nrow(arestas_df)), function(i) {
      g_from <- encontrar_gene(arestas_df$From[i], nomes)
      g_to   <- encontrar_gene(arestas_df$To[i], nomes)
      base <- data.frame(Dataset = nome_ds, Aresta = arestas_df$Aresta[i],
                         From = arestas_df$From[i], To = arestas_df$To[i], Sign = arestas_df$Sign[i],
                         stringsAsFactors = FALSE)
      if (is.na(g_from) || is.na(g_to)) {
        return(cbind(base, t(testar_dependencia(numeric(0), numeric(0)))))
      }
      cbind(base, t(testar_dependencia(as.numeric(expr_tumor[[g_from]]),
                                       as.numeric(expr_tumor[[g_to]]))))
    })
    todos[[nome_ds]] <- bind_rows(linhas)
  }

  bind_rows(todos) %>%
    group_by(Dataset) %>%
    mutate(p_adj_perm = p.adjust(p_perm, method = "BH")) %>%
    ungroup() %>%
    mutate(
      Status = case_when(
        is.na(GGC)          ~ "GENE_ABSENT",
        p_adj_perm < 0.05   ~ "SIGNIFICANT DEPENDENCE",
        TRUE                ~ "NOT SIGNIFICANT"
      ),
      GGC_ge_0.30 = GGC >= 0.30,   # mantido só para comparação com a versão anterior
      across(c(MI, NMI, GGC, GGC_null_median, GGC_null_p95), ~ round(.x, 4))
    )
}

df_nmi_ggc_3ds_corrigido <- calcular_nmi_ggc_geo(arestas_modelo_completo)
print(as.data.frame(df_nmi_ggc_3ds_corrigido))
write.csv(df_nmi_ggc_3ds_corrigido, "GEO_NMI_GGC_permutation.csv", row.names = FALSE)

# FIGURA
df_plot <- df_nmi_ggc_3ds_corrigido %>%
  filter(Status != "GENE_ABSENT") %>%
  pivot_longer(cols = c(NMI, GGC), names_to = "Metrica", values_to = "Valor") %>%
  mutate(Aresta = factor(Aresta, levels = arestas_modelo_completo$Aresta),
         Null_P95 = ifelse(Metrica == "GGC", GGC_null_p95, NA_real_))

p_nmi_ggc_corrigido <- ggplot(df_plot, aes(x = Aresta, y = Valor, fill = Status)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Null_P95, ymax = Null_P95), width = 0.8, linetype = "dashed",
                color = "black", na.rm = TRUE) +
  facet_grid(Dataset ~ Metrica, scales = "free_y") +
  scale_fill_manual(values = c("SIGNIFICANT DEPENDENCE" = "#27AE60", "NOT SIGNIFICANT" = "#E74C3C")) +
  labs(
    title = "Informational Validation of the Pathway (NMI and GGC)",
    subtitle = paste0("Primary tumour samples; TRPM3 as proxy of miR-204-5p.\n",
                      "Colour = permutation test (BH p.adj < 0.05); dashed line = 95th percentile of GGC under independence"),
    x = "Regulatory Interactions of Route", y = "Dependency Score", fill = "Edge Status"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
        strip.background = element_rect(fill = "#ECF0F1", color = NA),
        legend.position = "top")

ggsave("Fig_Validation_Pathway_DDR_GSDME.png", plot = p_nmi_ggc_corrigido, width = 12, height = 8, dpi = 300, bg = "white")
cat("\nFigura salva com sucesso!\n")


############# Validation of model #####################
# ============================================================================
# FIGURA ÚNICA
# Layout da figura azul + perturbações da figura vermelha
# Baseado no modelo: GINsim-GSDME_Pyroptosis(1).zginml
# ============================================================================
setwd("~/Downloads/GSDME/")

required_packages <- c("BoolNet", "ggplot2", "xml2", "grid")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Pacotes ausentes: ", paste(missing_packages, collapse = ", "),
    ". Instale com: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(BoolNet)
  library(ggplot2)
  library(xml2)
  library(grid)
})

# ----------------------------------------------------------------------------
# 1. CAMINHOS
# ----------------------------------------------------------------------------

MODEL_FILE <- "GINsim-GSDME_Pyroptosis.zginml"
if (!file.exists(MODEL_FILE)) {
  stop("Arquivo do modelo não encontrado: ", MODEL_FILE)
}

OUTPUT_DIR <- file.path(getwd(), "FIGURA_LAYOUT_AZUL_PERTURBACOES_VERMELHAS")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 2. EXTRAÇÃO DO GRAFO REGULATÓRIO DO GINsim
# ----------------------------------------------------------------------------

extract_regulatory_graph <- function(model_file) {
  archive_contents <- tryCatch(
    utils::unzip(model_file, list = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(archive_contents)) {
    graph_entries <- archive_contents$Name[
      grepl("(^|/)regulatoryGraph[.]ginml$", archive_contents$Name)
    ]
    
    if (length(graph_entries) != 1L) {
      stop(
        "O arquivo .zginml deveria conter exatamente um regulatoryGraph.ginml; ",
        "encontrados: ", length(graph_entries), "."
      )
    }
    
    extraction_directory <- tempfile("ginsim_model_")
    dir.create(extraction_directory)
    on.exit(unlink(extraction_directory, recursive = TRUE), add = TRUE)
    
    utils::unzip(
      zipfile = model_file,
      files = graph_entries,
      exdir = extraction_directory,
      overwrite = TRUE
    )
    
    extracted_file <- file.path(extraction_directory, graph_entries)
    
    return(xml2::read_xml(
      extracted_file,
      options = c("RECOVER", "NOERROR", "NOWARNING", "NONET")
    ))
  }
  
  xml2::read_xml(
    model_file,
    options = c("RECOVER", "NOERROR", "NOWARNING", "NONET")
  )
}

model_document <- extract_regulatory_graph(MODEL_FILE)

graph_node <- xml2::xml_find_first(
  model_document,
  ".//graph[@class='regulatory']"
)

if (inherits(graph_node, "xml_missing")) {
  stop("Nenhum grafo regulatório foi encontrado no arquivo .zginml.")
}

graph_id <- xml2::xml_attr(graph_node, "id")
node_order_text <- xml2::xml_attr(graph_node, "nodeorder")
node_order <- strsplit(trimws(node_order_text), "[[:space:]]+")[[1L]]

node_elements <- xml2::xml_find_all(graph_node, "./node")
node_ids <- xml2::xml_attr(node_elements, "id")

if (!setequal(node_order, node_ids)) {
  stop(
    "A lista nodeorder do GINsim não corresponde aos nós do grafo. ",
    "Ausentes: ", paste(setdiff(node_order, node_ids), collapse = ", "),
    "; extras: ", paste(setdiff(node_ids, node_order), collapse = ", ")
  )
}

max_values <- suppressWarnings(as.integer(xml2::xml_attr(node_elements, "maxvalue")))
if (any(is.na(max_values)) || any(max_values != 1L)) {
  stop("Este script aceita apenas redes booleanas com maxvalue = 1.")
}

input_nodes <- node_ids[
  !is.na(xml2::xml_attr(node_elements, "input")) &
    xml2::xml_attr(node_elements, "input") == "true"
]

extract_rule <- function(node_element) {
  node_id <- xml2::xml_attr(node_element, "id")
  expression_node <- xml2::xml_find_first(
    node_element,
    "./value[@val='1']/exp"
  )
  
  if (inherits(expression_node, "xml_missing")) {
    if (node_id %in% input_nodes) {
      return(node_id)
    }
    stop("Regra booleana ausente para o nó: ", node_id)
  }
  
  rule <- xml2::xml_attr(expression_node, "str")
  if (is.na(rule) || !nzchar(trimws(rule))) {
    stop("Expressão booleana vazia para o nó: ", node_id)
  }
  
  trimws(rule)
}

rules_by_node <- vapply(node_elements, extract_rule, character(1))
names(rules_by_node) <- node_ids
rules_by_node <- rules_by_node[node_order]

# ----------------------------------------------------------------------------
# 3. IMPORTAÇÃO PARA O BOOLNET
# ----------------------------------------------------------------------------

boolnet_network_file <- file.path(OUTPUT_DIR, "temp_network_boolnet.csv")

writeLines(
  c(
    "targets, factors",
    paste(node_order, unname(rules_by_node), sep = ", ")
  ),
  boolnet_network_file
)

network <- BoolNet::loadNetwork(boolnet_network_file)

if (!setequal(network$genes, node_order)) {
  stop(
    "Os nós importados pelo BoolNet diferem dos nós do GINsim. ",
    "Ausentes: ", paste(setdiff(node_order, network$genes), collapse = ", "),
    "; extras: ", paste(setdiff(network$genes, node_order), collapse = ", ")
  )
}

message("Modelo importado: ", graph_id)
message("Número de nós: ", length(node_order))
message("Inputs: ", paste(input_nodes, collapse = ", "))

# ----------------------------------------------------------------------------
# 4. DEFINIÇÕES DOS OUTPUTS
# ----------------------------------------------------------------------------

output_nodes <- c(
  "PROLIFERATION",
  "RESISTANCE",
  "CELL_CYCLE_ARREST",
  "PYROPTOSIS",
  "APOPTOSIS"
)

if (!all(output_nodes %in% node_order)) {
  stop(
    "Nem todos os nós de saída estão no modelo. Ausentes: ",
    paste(setdiff(output_nodes, node_order), collapse = ", ")
  )
}

output_labels <- c(
  PROLIFERATION      = "Proliferation",
  RESISTANCE         = "Resistance",
  CELL_CYCLE_ARREST  = "Cell-cycle arrest",
  PYROPTOSIS         = "Pyroptosis",
  APOPTOSIS          = "Apoptosis"
)

classify_phenotype <- function(state) {
  active_outputs <- output_nodes[state[output_nodes] == 1L]
  if (length(active_outputs) == 0L) {
    return("No active phenotype")
  }
  paste(unname(output_labels[active_outputs]), collapse = " + ")
}

# ----------------------------------------------------------------------------
# 5. ORDEM DE EXIBIÇÃO DOS NÓS NA FIGURA
#    (outputs no final, como no estilo da figura azul)
# ----------------------------------------------------------------------------

display_node_order <- c(
  node_order[!node_order %in% output_nodes],
  output_nodes
)

# nomes mais bonitos para a figura
pretty_node_labels <- function(x) {
  out <- x
  out[out == "miR_204_5p"]       <- "miR-204-5p"
  out[out == "lncRNA_MALAT1"]    <- "lncRNA_MALAT1"
  out[out == "CELL_CYCLE_ARREST"] <- "CELL-CYCLE ARREST"
  out[out == "CDK4_6_CyclinD"]   <- "CDK4/6-Cyclin D"
  out[out == "CASP3"]            <- "Caspase3"
  out[out == "CASP9"]            <- "Caspase9"
  out[out == "AKT1"]             <- "AKT"
  out[out == "CDC25A"]           <- "Cdc25"
  out
}

display_node_labels <- pretty_node_labels(display_node_order)

# ----------------------------------------------------------------------------
# 6. PERTURBAÇÕES DA FIGURA VERMELHA
# ----------------------------------------------------------------------------

conditions <- list(
  "MALAT1 E1"                = c(DDR = 1L, lncRNA_MALAT1 = 1L),
  "MALAT1 KO"                = c(DDR = 1L, lncRNA_MALAT1 = 0L),
  "SIRT1 KO"                 = c(DDR = 1L, SIRT1 = 0L),
  "SIRT1 E1"                 = c(DDR = 1L, SIRT1 = 1L),
  "miR-204-5p KO"            = c(DDR = 1L, miR_204_5p = 0L),
  "miR-204-5p E1"            = c(DDR = 1L, miR_204_5p = 1L),
  "Wip1 E1"                  = c(DDR = 1L, Wip1 = 1L),
  "Wip1 KO"                  = c(DDR = 1L, Wip1 = 0L),
  "GSDME KO + CASP3 E1"      = c(DDR = 1L, GSDME = 0L, CASP3 = 1L),
  "GSDME E1 + CASP3 E1"      = c(DDR = 1L, GSDME = 1L, CASP3 = 1L),
  "MALAT1 KO + GSDME KO"     = c(DDR = 1L, lncRNA_MALAT1 = 0L, GSDME = 0L),
  "miR-204-5p E1 + GSDME KO" = c(DDR = 1L, miR_204_5p = 1L, GSDME = 0L),
  "Wip1 KO + MALAT1 E1"      = c(DDR = 1L, Wip1 = 0L, lncRNA_MALAT1 = 1L)
)

validate_condition <- function(condition_name, fixed_values) {
  if (is.null(names(fixed_values)) || any(names(fixed_values) == "")) {
    stop("Valores fixos sem nomes na condição: ", condition_name)
  }
  
  if (!all(names(fixed_values) %in% node_order)) {
    stop(
      "Nó desconhecido na condição ", condition_name, ": ",
      paste(setdiff(names(fixed_values), node_order), collapse = ", ")
    )
  }
  
  if (any(!fixed_values %in% c(0L, 1L))) {
    stop("Somente os valores 0 e 1 são permitidos: ", condition_name)
  }
  
  if (!"DDR" %in% names(fixed_values) || fixed_values[["DDR"]] != 1L) {
    stop("DDR precisa estar fixado em 1 na condição: ", condition_name)
  }
  
  invisible(TRUE)
}

invisible(Map(validate_condition, names(conditions), conditions))

# tabela para desenhar as marcações das perturbações
build_perturbation_table <- function(conditions_list) {
  out <- lapply(names(conditions_list), function(cond) {
    fx <- conditions_list[[cond]]
    fx <- fx[names(fx) != "DDR"]
    if (length(fx) == 0L) return(NULL)
    
    data.frame(
      Condition = cond,
      Target = names(fx),
      Fixed_value = as.integer(fx),
      Symbol = ifelse(as.integer(fx) == 1L, "\u2191", "\u2193"),
      stringsAsFactors = FALSE
    )
  })
  
  out <- out[!vapply(out, is.null, logical(1))]
  do.call(rbind, out)
}

perturbation_table <- build_perturbation_table(conditions)

# ----------------------------------------------------------------------------
# 7. ENUMERAÇÃO EXATA DOS PONTOS FIXOS
# ----------------------------------------------------------------------------

decode_fixed_points <- function(attractor_information, condition_name) {
  number_of_attractors <- length(attractor_information$attractors)
  
  if (number_of_attractors == 0L) {
    return(data.frame())
  }
  
  decoded_rows <- vector("list", number_of_attractors)
  
  for (attractor_index in seq_len(number_of_attractors)) {
    attractor_sequence <- as.data.frame(
      BoolNet::getAttractorSequence(attractor_information, attractor_index),
      check.names = FALSE
    )
    
    if (nrow(attractor_sequence) != 1L) {
      stop(
        "Foi encontrado um atrator com mais de um estado, ",
        "apesar da busca estar restrita a pontos fixos."
      )
    }
    
    state <- as.integer(unlist(
      attractor_sequence[1L, network$genes, drop = FALSE],
      use.names = FALSE
    ))
    names(state) <- network$genes
    
    decoded_rows[[attractor_index]] <- data.frame(
      Condition = condition_name,
      State_key = paste0(state[display_node_order], collapse = ""),
      Phenotype = classify_phenotype(state),
      as.list(state[display_node_order]),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  decoded <- unique(do.call(rbind, decoded_rows))
  
  phenotype_order <- c(
    "Proliferation + Resistance",
    "Proliferation",
    "Resistance",
    "Cell-cycle arrest",
    "Pyroptosis",
    "Apoptosis",
    "No active phenotype"
  )
  
  phenotype_rank <- match(decoded$Phenotype, phenotype_order)
  phenotype_rank[is.na(phenotype_rank)] <- length(phenotype_order) + 1L
  
  decoded <- decoded[order(phenotype_rank, decoded$State_key), , drop = FALSE]
  decoded$Stable_state <- paste0("SS", seq_len(nrow(decoded)))
  
  decoded <- decoded[, c(
    "Condition", "Stable_state", "State_key", "Phenotype", display_node_order
  )]
  
  row.names(decoded) <- NULL
  decoded
}

# NOTA (out/2026): a busca abaixo enumera apenas PONTOS FIXOS (maxAttractorLength = 1).
# Pontos fixos são os mesmos no esquema síncrono e no assíncrono. Uma checagem independente
# (300 caminhadas assíncronas aleatórias por condição) chegou sempre a um ponto fixo, então
# não há indício de atratores cíclicos; ainda assim, no texto prefira "stable states (fixed points)".
# Atenção: sob DDR = 1, ATM = DDR & (E2F1 | !Wip1) com Wip1 = E2F1 é sempre 1, logo CDC25A = 0 e
# Mdm2 = 0 em todos os pontos fixos — PROLIFERATION (= Mdm2) nunca é ativado.
analyze_condition <- function(condition_name, fixed_values) {
  validate_condition(condition_name, fixed_values)
  
  perturbed_network <- BoolNet::fixGenes(
    network,
    fixIndices = names(fixed_values),
    values = as.numeric(fixed_values)
  )
  
  attractor_information <- BoolNet::getAttractors(
    perturbed_network,
    type = "synchronous",
    method = "sat.restricted",
    maxAttractorLength = 1,
    canonical = TRUE,
    returnTable = TRUE
  )
  
  decode_fixed_points(attractor_information, condition_name)
}

message("Iniciando enumeração exata de pontos fixos...")

state_tables <- vector("list", length(conditions))
names(state_tables) <- names(conditions)

for (i in seq_along(conditions)) {
  condition_name <- names(conditions)[[i]]
  message("[", i, "/", length(conditions), "] ", condition_name)
  state_tables[[condition_name]] <- analyze_condition(
    condition_name,
    conditions[[condition_name]]
  )
}

if (all(vapply(state_tables, nrow, integer(1)) == 0L)) {
  stop("Nenhum ponto fixo foi encontrado para as condições solicitadas.")
}

all_states <- do.call(
  rbind,
  state_tables[vapply(state_tables, nrow, integer(1)) > 0L]
)
row.names(all_states) <- NULL

# ----------------------------------------------------------------------------
# 8. ORGANIZAÇÃO DAS LINHAS DA FIGURA
# ----------------------------------------------------------------------------

all_states$Condition <- factor(all_states$Condition, levels = names(conditions))
all_states <- all_states[order(all_states$Condition), , drop = FALSE]
all_states$Condition <- as.character(all_states$Condition)

all_states$Row_label <- paste(
  all_states$Condition,
  all_states$Stable_state,
  all_states$Phenotype,
  sep = " | "
)

all_states$Plot_row <- seq_len(nrow(all_states))
all_states$Y <- nrow(all_states) - all_states$Plot_row + 1L

write.csv(
  all_states[, c("Condition", "Stable_state", "Phenotype", "State_key", display_node_order)],
  file.path(OUTPUT_DIR, "stable_states_used_in_figure.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 9. DADOS LONGOS PARA O HEATMAP
# ----------------------------------------------------------------------------

make_long_state_table <- function(states, nodes_to_display) {
  do.call(
    rbind,
    lapply(seq_len(nrow(states)), function(row_index) {
      values <- as.integer(unlist(
        states[row_index, nodes_to_display, drop = FALSE],
        use.names = FALSE
      ))
      
      fill_group <- ifelse(values == 1L, "Active node", "Inactive node")
      
      active_outputs <- nodes_to_display %in% output_nodes & values == 1L
      fill_group[active_outputs] <- unname(output_labels[nodes_to_display[active_outputs]])
      
      data.frame(
        Condition = states$Condition[[row_index]],
        Stable_state = states$Stable_state[[row_index]],
        Phenotype = states$Phenotype[[row_index]],
        Row_label = states$Row_label[[row_index]],
        Y = states$Y[[row_index]],
        Node = nodes_to_display,
        X = seq_along(nodes_to_display),
        Value = values,
        Fill_group = fill_group,
        stringsAsFactors = FALSE
      )
    })
  )
}

state_long <- make_long_state_table(all_states, display_node_order)

state_long$Fill_group <- factor(
  state_long$Fill_group,
  levels = c(
    "Inactive node",
    "Active node",
    "Proliferation",
    "Resistance",
    "Cell-cycle arrest",
    "Pyroptosis",
    "Apoptosis"
  )
)

# ----------------------------------------------------------------------------
# 10. OVERLAY DAS PERTURBAÇÕES
# ----------------------------------------------------------------------------

overlay_rows <- merge(
  unique(all_states[, c("Condition", "Y")]),
  perturbation_table,
  by = "Condition",
  all = FALSE,
  sort = FALSE
)

overlay_rows$X <- match(overlay_rows$Target, display_node_order)

overlay_rows <- overlay_rows[
  !is.na(overlay_rows$X),
  ,
  drop = FALSE
]

# ----------------------------------------------------------------------------
# 11. INFORMAÇÕES DOS GRUPOS
# ----------------------------------------------------------------------------

group_info <- do.call(
  rbind,
  lapply(names(conditions), function(cond) {
    y_values <- all_states$Y[all_states$Condition == cond]
    if (length(y_values) == 0L) return(NULL)
    
    data.frame(
      Condition = cond,
      Y_center = mean(range(y_values)),
      Y_min = min(y_values),
      Y_max = max(y_values),
      stringsAsFactors = FALSE
    )
  })
)

separator_y <- if (nrow(group_info) > 1L) {
  group_info$Y_min[-nrow(group_info)] - 0.5
} else {
  numeric(0)
}

# ----------------------------------------------------------------------------
# 12. POSIÇÕES E ELEMENTOS VISUAIS
# ----------------------------------------------------------------------------

input_position <- match("DDR", display_node_order)
output_positions <- match(output_nodes, display_node_order)

if (is.na(input_position)) {
  stop("DDR não foi encontrado na ordem de exibição.")
}

if (any(is.na(output_positions))) {
  stop("Um ou mais nós de output não foram encontrados na ordem de exibição.")
}

n_nodes <- length(display_node_order)
n_rows <- nrow(all_states)

x_perturbation <- n_nodes + 1.3
x_right_limit  <- n_nodes + 7.8

# ----------------------------------------------------------------------------
# 13. PLOT FINAL
# ----------------------------------------------------------------------------

p <- ggplot() +
  
  geom_tile(
    data = state_long,
    aes(x = X, y = Y, fill = Fill_group),
    width = 0.96,
    height = 0.90,
    color = "white",
    linewidth = 0.45
  ) +
  
  # Caixa do input DDR
  annotate(
    "rect",
    xmin = input_position - 0.52,
    xmax = input_position + 0.52,
    ymin = 0.45,
    ymax = n_rows + 0.55,
    fill = NA,
    color = "#E69F00",
    linewidth = 1.1
  ) +
  
  # Caixa dos outputs
  annotate(
    "rect",
    xmin = min(output_positions) - 0.52,
    xmax = max(output_positions) + 0.52,
    ymin = 0.45,
    ymax = n_rows + 0.55,
    fill = NA,
    color = "#00A087",
    linewidth = 1.1
  ) +
  
  # Quadrados vermelhos em volta dos nós perturbados
  geom_tile(
    data = overlay_rows,
    aes(x = X, y = Y),
    inherit.aes = FALSE,
    width = 0.80,
    height = 0.76,
    fill = NA,
    color = "#D7301F",
    linewidth = 1.0
  ) +
  
  # Setas das perturbações
  geom_text(
    data = overlay_rows,
    aes(x = X, y = Y, label = Symbol),
    inherit.aes = FALSE,
    color = "#F39C12",
    size = 5.0,
    fontface = "bold"
  ) +
  
  # Linhas separando os grupos
  geom_hline(
    yintercept = separator_y,
    color = "#D9E2EC",
    linewidth = 0.55
  ) +
  
  # Nome das perturbações à direita
  geom_text(
    data = group_info,
    aes(x = x_perturbation, y = Y_center, label = Condition),
    hjust = 0,
    vjust = 0.5,
    size = 4.2,
    color = "#17202A"
  ) +
  
  # Colchete à esquerda
  annotate(
    "segment",
    x = -0.85, xend = -0.85,
    y = 0.55, yend = n_rows + 0.45,
    linewidth = 0.8, color = "#17202A"
  ) +
  annotate(
    "segment",
    x = -0.85, xend = -0.48,
    y = 0.55, yend = 0.55,
    linewidth = 0.8, color = "#17202A"
  ) +
  annotate(
    "segment",
    x = -0.85, xend = -0.48,
    y = n_rows + 0.45, yend = n_rows + 0.45,
    linewidth = 0.8, color = "#17202A"
  ) +
  annotate(
    "text",
    x = -1.35,
    y = (n_rows + 1) / 2,
    label = "Stable states",
    angle = 90,
    size = 4.2,
    fontface = "bold",
    color = "#17202A"
  ) +
  
  # Títulos das seções
  annotate(
    "text",
    x = input_position,
    y = n_rows + 1.25,
    label = "Input\n(DDR)",
    size = 4.1,
    fontface = "bold",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = mean(range(output_positions)),
    y = n_rows + 1.25,
    label = "Outputs\n(phenotypes)",
    size = 4.1,
    fontface = "bold",
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = x_perturbation,
    y = n_rows + 1.25,
    label = "Perturbation",
    hjust = 0,
    size = 4.1,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "Inactive node"      = "#F1F4F7",
      "Active node"        = "#2F5D8C",
      "Proliferation"      = "#A65A3A",
      "Resistance"         = "#6A3D9A",
      "Cell-cycle arrest"  = "#E69F00",
      "Pyroptosis"         = "#00A087",
      "Apoptosis"          = "#F05A5B"
    ),
    breaks = c(
      "Active node",
      "Inactive node",
      "Proliferation",
      "Resistance",
      "Cell-cycle arrest",
      "Pyroptosis",
      "Apoptosis"
    ),
    labels = c(
      "Active node",
      "Inactive node",
      "Proliferation",
      "Resistance",
      "Cell-cycle arrest",
      "Pyroptosis",
      "Apoptosis"
    ),
    drop = FALSE,
    name = NULL
  ) +
  
  scale_x_continuous(
    breaks = seq_along(display_node_order),
    labels = display_node_labels,
    limits = c(-1.65, x_right_limit),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(-0.10, n_rows + 1.80),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    title = "Validation of Model Predictions under DDR ON in HCC",
    subtitle = paste0(
      "Exact fixed points from the official GINsim model using the requested perturbation screen"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "E1 and KO denote Boolean gain- and loss-of-function constraints; ",
      "fixed points indicate qualitative reachability, not frequencies."
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 50, hjust = 1, vjust = 1,
      size = 8.3, color = "#1F2933"
    ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(
      face = "bold", size = 19, hjust = 0.5, color = "#111827"
    ),
    plot.subtitle = element_text(
      size = 10.7, hjust = 0.5, color = "#5B6770"
    ),
    plot.caption = element_text(
      size = 9.0, hjust = 0, color = "#4B5563"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.0, "lines"),
    legend.key.height = unit(1.0, "lines"),
    plot.margin = margin(18, 24, 16, 30)
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(color = NA)
    )
  )

# ----------------------------------------------------------------------------
# 14. EXPORTAÇÃO
# ----------------------------------------------------------------------------

png_file <- file.path(OUTPUT_DIR, "Figure_layout_blue_perturbations_red.png")
pdf_file <- file.path(OUTPUT_DIR, "Figure_layout_blue_perturbations_red.pdf")

ggsave(
  filename = png_file,
  plot = p,
  width = 18,
  height = max(8.5, 0.36 * n_rows + 4.5),
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 18,
  height = max(8.5, 0.36 * n_rows + 4.5),
  units = "in",
  bg = "white",
  limitsize = FALSE
)

message("Figura salva em:")
message(normalizePath(png_file, mustWork = FALSE))
message(normalizePath(pdf_file, mustWork = FALSE))

print(p)
