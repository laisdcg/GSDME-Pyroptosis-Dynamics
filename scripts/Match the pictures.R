install.packages("magick")
library(magick)

juntar <- function(cima, baixo, saida, largura = 4000) {
  a <- image_scale(image_read(cima),  as.character(largura))
  b <- image_scale(image_read(baixo), as.character(largura))
  a <- image_annotate(image_border(a, "white", "60x0"), "A", size = 90, weight = 700, location = "+10+10")
  b <- image_annotate(image_border(b, "white", "60x0"), "B", size = 90, weight = 700, location = "+10+10")
  fig <- image_append(c(a, image_blank(largura + 120, 60, "white"), b), stack = TRUE)
  image_write(fig, saida, density = 300)
}

# Fig. 5 — Spearman
juntar("FIGURE_1_edge_validation_spearman.png",
       "resultados_TCGA_LIHC/figures/Fig_Edge_Validation_Spearman_TCGA_EN.png",
       "Figure5_Spearman_GEO_TCGA.png")

# Fig. 6 — NMI/GGC
juntar("Fig_Validation_Pathway_DDR_GSDME.png",
       "resultados_TCGA_LIHC/figures/Fig_Pathway_Validation_DDR_GSDME_TCGA_EN.png",
       "Figure6_NMI_GGC_GEO_TCGA.png")