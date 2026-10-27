if (!requireNamespace("magick", quietly = TRUE)) install.packages("magick")
library(magick)

# Junta o painel GEO (A, em cima) e o painel TCGA-LIHC (B, embaixo) numa figura só.
juntar <- function(cima, baixo, saida, largura = 4000) {
  a <- image_scale(image_read(cima),  as.character(largura))
  b <- image_scale(image_read(baixo), as.character(largura))
  a <- image_annotate(image_border(a, "white", "60x0"), "A", size = 90, weight = 700, location = "+10+10")
  b <- image_annotate(image_border(b, "white", "60x0"), "B", size = 90, weight = 700, location = "+10+10")
  fig <- image_append(c(a, image_blank(largura + 120, 60, "white"), b), stack = TRUE)
  image_write(fig, saida, density = 300)
}

TCGA_FIG <- "resultados_TCGA_LIHC/figures"

# Fig. 6 — Spearman, via principal
juntar("FIGURE_1_edge_validation_spearman.png",
       file.path(TCGA_FIG, "Fig_Edge_Validation_Spearman_TCGA_EN.png"),
       "Figure6_Spearman_GEO_TCGA.png")

# Fig. 7 — Spearman, ramo pró-morte (piroptótico)
juntar("FIGURE_1b_pyroptotic_branch_spearman.png",
       file.path(TCGA_FIG, "Fig_Pyroptotic_Branch_Spearman_TCGA_EN.png"),
       "Figure7_Spearman_PyroptoticBranch_GEO_TCGA.png")

# Fig. 8 — NMI/GGC, via principal
juntar("Fig_Validation_Pathway_DDR_GSDME.png",
       file.path(TCGA_FIG, "Fig_Pathway_Validation_DDR_GSDME_TCGA_EN.png"),
       "Figure8_NMI_GGC_GEO_TCGA.png")

# Fig. 9 — NMI/GGC, ramo pró-morte (piroptótico)
juntar("Fig_Validation_Pyroptotic_Branch_GEO.png",
       file.path(TCGA_FIG, "Fig_Pyroptotic_Branch_Validation_TCGA_EN.png"),
       "Figure9_NMI_GGC_PyroptoticBranch_GEO_TCGA.png")
