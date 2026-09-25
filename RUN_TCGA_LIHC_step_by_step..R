# Abra este arquivo no RStudio e execute os blocos UM POR VEZ (Ctrl+Enter).
# Nao use source() neste arquivo inteiro: ele e um roteiro de comandos.

# 1. Instale apenas uma vez, se necessario.
install.packages(c('httr', 'jsonlite', 'xml2', 'ggplot2'))

# 2. Selecione GSDME_TCGA_LIHC.R na pasta extraida do ZIP.
arquivo_script <- file.choose()
source(arquivo_script, encoding = 'UTF-8')

# 3. Selecione o SEU modelo atual (.zginml ou .ginml).
arquivo_modelo <- file.choose()
cfg <- tcga_config(
  model = arquivo_modelo,
  out = file.path(dirname(arquivo_modelo), 'resultados_TCGA_LIHC')
)

# 4. Consulte a coorte e confira os nomes/amostras no manifesto.
job <- tcga_prepare(cfg)
View(job$mapping)
View(job$manifest)

# 5. Baixe os arquivos. Pode demorar; nao feche o R.
# Se a internet cair, execute esta mesma linha novamente.
tcga_download(job)

# 6. Extraia a expressao dos genes e do miR-204-5p maduro.
expressao <- tcga_extract(job)

# 7. Gere tabelas e figuras (PDF e PNG 600 dpi; textos em ingles).
resultados <- tcga_analyze(job, expressao)

# 8. Confira as comparacoes pareadas e as correlacoes.
View(resultados$paired_tests)
View(resultados$correlations)
cat('Resultados em:', cfg$out, '\n')

# PARA RETOMAR EM OUTRA SESSAO, sem consultar tudo novamente:
# source(file.choose(), encoding = 'UTF-8')  # escolha GSDME_TCGA_LIHC.R
# job <- readRDS(file.choose())             # escolha resultados_TCGA_LIHC/TCGA_job.rds
# tcga_download(job)
# expressao <- tcga_extract(job)
# resultados <- tcga_analyze(job, expressao)

# SE MUDAR O MODELO: execute novamente os passos 3, 4, 6 e 7.
# Os arquivos GDC ja baixados permanecem no cache.
# Se quiser atualizar a consulta GDC: job <- tcga_prepare(cfg, refresh = TRUE)

# ATUALIZAR SOMENTE AS FIGURAS DE EXPRESSAO (dados ja salvos):
# source(file.choose(), encoding = 'UTF-8') # escolha o novo modulo
# job <- readRDS(file.choose())            # escolha TCGA_job.rds
# tcga_plot_expression(job)
