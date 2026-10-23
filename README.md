# GSDME-Pyroptosis-Dynamics

Code, model and results for the manuscript:

> **Dynamic Modeling and Transcriptomic Analysis Reveal GSDME as the Execution Switch between Apoptosis and Pyroptosis in Hepatocellular Carcinoma**
> Laís de Carvalho Gonçalves, Djorkaeff Oliveira Fontinele, João Paulo Matos Santos Lima, Shantanu Gupta
> Bioinformatics Multidisciplinary Environment (BioME), Digital Metropolis Institute (IMD), Universidade Federal do Rio Grande do Norte (UFRN), Natal, Brazil

This repository contains a Boolean regulatory model of the MALAT1 / miR-204-5p / SIRT1 / p53 / CASP3 / GSDME axis in hepatocellular carcinoma (HCC), the in silico perturbation analysis, and the transcriptomic validation of the model edges in three GEO cohorts and in TCGA-LIHC.

---

## Overview

The study asks which regulatory conditions allow CASP3 activation to end in GSDME-dependent pyroptosis rather than apoptosis or tumour-cell survival. The workflow has three steps:

1. **Model construction.** A logical model built in GINsim from prior knowledge.
2. **Model perturbation.** Stable states under persistent DNA damage (`DDR = 1`), followed by single and combined gain-of-function (E1, node = 1) and loss-of-function (KO, node = 0) perturbations.
3. **Two-pronged validation.**
   - **(3A) Experimental:** the predicted phenotypes are compared with phenotypes reported in independent GoF/LoF experimental studies.
   - **(3B) Database:** the model edges are assessed in patient data (GEO and TCGA-LIHC) using Spearman's correlation (direction) and NMI/GGC (dependence).

**Main result:** upstream regulators (MALAT1, miR-204-5p, SIRT1, Wip1) determine whether the network reaches death-associated states. GSDME availability then decides the execution phenotype: CASP3 with functional GSDME leads to pyroptosis, and CASP3 without GSDME leads to apoptosis. In patient data, CDC25A–E2F1 and CASP3–GSDME are the most consistently supported interactions.

---

## Repository structure

```
GSDME-Pyroptosis-Dynamics/
├── model/
│   ├── GINsim-GSDME_Pyroptosis.zginml        # Boolean model (miR-204-5p node named TRPM3; used for GEO)
│   └── GINsim-GSDME_Pyroptosis_TCGA.zginml   # Same model (miR-204-5p node named MIR204; used for TCGA)
├── scripts/
│   ├── bulkrna-seq-gsdme.R      # GEO analyses + Boolean perturbation analysis (BoolNet)
│   ├── tcga-lihc-gsdme.R        # TCGA-LIHC download, expression, Spearman and NMI/GGC
│   └── Match the pictures.R     # Combines GEO (A) and TCGA (B) panels into Figures 5 and 6
├── results/                     # Tables (.csv) and figures (.png/.pdf) produced by the scripts
├── LICENSE
└── README.md
```

The two `.zginml` files contain **identical** logical rules (31 nodes, 68 interactions). They differ only in the name given to the `miR_204_5p` node. In the GEO microarrays, the mature miRNA is not measured, so its host gene **TRPM3** is used as a proxy. In TCGA-LIHC, mature **miR-204-5p** (MIMAT0000265) is quantified directly from miRNA-seq.

---

## Boolean model

- **31 nodes and 68 regulatory interactions**, including 5 phenotypic outputs: `PROLIFERATION`, `RESISTANCE`, `CELL_CYCLE_ARREST`, `APOPTOSIS`, `PYROPTOSIS`.
- `DDR` is an input node kept active (`DDR = 1`) in all simulations, representing persistent DNA-damage signalling.
- Key rules:
  - `GSDME = CASP3 & p53_K`
  - `PYROPTOSIS = CASP3 & GSDME`
  - `APOPTOSIS = CASP3 & !GSDME`
- Stable states are obtained by exact enumeration of fixed points with **BoolNet**, using the rules imported directly from the `.zginml` file.

Perturbations analysed (all with `DDR = 1`):

| Single | Combined |
|---|---|
| MALAT1 E1 / KO | GSDME KO + CASP3 E1 |
| SIRT1 E1 / KO | GSDME E1 + CASP3 E1 |
| miR-204-5p E1 / KO | MALAT1 KO + GSDME KO |
| Wip1 E1 / KO | miR-204-5p E1 + GSDME KO |
| | Wip1 KO + MALAT1 E1 |

---

## Transcriptomic validation

### Cohorts

| Cohort | Platform / assay | Samples used |
|---|---|---|
| GSE14520 | Affymetrix GPL3921 | 222 HCC, 212 adjacent non-tumour |
| GSE60502 | Affymetrix GPL96 | 18 HCC, 18 adjacent non-tumour (paired) |
| GSE121248 | Affymetrix GPL570 | 70 HCC, 37 adjacent non-tumour |
| TCGA-LIHC | RNA-seq (STAR – Counts, TPM) and miRNA-seq (isoform quantification) | 359 primary tumours; 49 tumour–normal pairs |

Tumour and non-tumour groups are taken from the sample annotation of each GEO series and from the GDC sample type (`Primary Tumor`, `Solid Tissue Normal`).

### Pathway edges evaluated

```text
ATM ┤ CDC25A → E2F1 → MYC → MALAT1 ┤ miR-204-5p ┤ SIRT1 ┤ TP53 → CDKN1A ┤ CASP3 → GSDME
```

(TRPM3 replaces miR-204-5p in the GEO cohorts.)

### Statistics

- **Differential expression**
  - GEO: two-sided Wilcoxon rank-sum test on log2 series-matrix values, with Benjamini–Hochberg (BH) correction within each cohort.
  - TCGA-LIHC: paired Wilcoxon signed-rank test on log2(TPM + 1) or log2(RPM + 1), with BH correction.
- **Spearman correlation**
  - Primary tumour samples only, at least 5 valid pairs per edge, with BH correction within each cohort.
  - An edge is *incompatible* only when the correlation is significant (adjusted p < 0.05) in the direction opposite to the model sign. An asterisk marks edges that are significant in the predicted direction.
- **NMI / GGC**
  - Primary tumour samples only. Data are discretised with `infotheo::discretize` (equal frequency, N^(1/3) bins), and mutual information is estimated in nats.
  - NMI = 2·MI / (H₁ + H₂); GGC = √(1 − e^(−2·MI)).
  - Significance comes from a permutation test (2,000 permutations of the target gene) with BH correction. GGC is reported as the effect size, together with the 95th percentile of its null distribution.
  - A fixed GGC threshold is **not** used: with these sample sizes, independent genes already reach GGC ≈ 0.30.

---

## How to reproduce

### Requirements

R (≥ 4.2) and the following packages:

```r
# CRAN
install.packages(c("tidyverse", "ggplot2", "dplyr", "tidyr", "readr", "patchwork",
                   "ggrepel", "igraph", "infotheo", "BoolNet", "xml2", "XML",
                   "httr", "jsonlite", "magick", "caret"))
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Biobase"))
```

### 1. GEO analyses and Boolean perturbations: `scripts/bulkrna-seq-gsdme.R`

1. Put both `.zginml` files in a working folder (the script uses `~/Downloads/GSDME/` by default; edit `setwd()` at the top).
2. Run the script from top to bottom. It downloads GSE14520, GSE60502 and GSE121248 with GEOquery and then runs, in order:
   - component coverage
   - expression profiling and differential expression (Figure 3)
   - Spearman validation (Figure 5A)
   - NMI/GGC validation (Figure 6A)
   - Boolean perturbation analysis with BoolNet (Figure 7)

### 2. TCGA-LIHC analyses: `scripts/tcga-lihc-gsdme.R`

1. Run the script. When `file.choose()` opens, select **`GINsim-GSDME_Pyroptosis_TCGA.zginml`** (the version with the `MIR204` node).
2. The script queries the NCI GDC API and downloads the open-access STAR–Counts and miRNA isoform files. Each file is checked against its MD5 checksum. If the download is interrupted, run `tcga_download(job)` again: files that are already valid are skipped.
3. Outputs are written to `resultados_TCGA_LIHC/tables` and `resultados_TCGA_LIHC/figures`, next to the model file.

### 3. Combine panels: `scripts/Match the pictures.R`

This script stacks the GEO (A) and TCGA (B) panels into Figures 5 and 6.

---

## Figures and corresponding outputs

| Manuscript figure | Produced by | Main output files in `results/` |
|---|---|---|
| Fig. 1 – Workflow | Schematic | — |
| Fig. 2 – Boolean model | GINsim | `model/*.zginml` |
| Fig. 3 – GEO expression | `bulkrna-seq-gsdme.R` | `Fig4_GSE_Cohorts_Expression_Facet.png`, `GEO_differential_expression_model_components.csv`, `Fig_Expressao_Diferencial_GSE*.png` |
| Fig. 4 – TCGA expression | `tcga-lihc-gsdme.R` | `Figure_TCGA_Facet_Components_Style.png`, `06_paired_expression_tests.csv`, `04_expression_long.csv` |
| Fig. 5 – Spearman | both scripts + `Match the pictures.R` | `Figure5_Spearman_GEO_TCGA.png`, `boolean_model_edge_validation_COMPLETE.csv` (GEO), `boolean_model_edge_validation_TCGA_spearman.csv` |
| Fig. 6 – NMI/GGC | both scripts + `Match the pictures.R` | `Figure6_NMI_GGC_GEO_TCGA.png`, `GEO_NMI_GGC_permutation.csv`, `11_validacao_nmi_ggc_tcga.csv` |
| Fig. 7 – In silico perturbations | `bulkrna-seq-gsdme.R` (BoolNet) | `Figure_layout_blue_perturbations_red.png/.pdf`, `stable_states_used_in_figure.csv` |

TCGA selection tables: `01_model_name_mapping.csv` (model node → gene symbol), `02_all_files_selection_audit.csv` (every GDC file and why it was selected or excluded) and `03_download_manifest.csv` (files used).

---

## Data availability

All transcriptomic data are public:

- **NCBI GEO:** [GSE14520](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520), [GSE60502](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60502), [GSE121248](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121248)
- **NCI Genomic Data Commons:** [TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC) (open-access files only)

---

## Citation

If you use this model or code, please cite the manuscript above (citation details will be added upon publication).

## License

Code is released under the [MIT License](LICENSE). The GINsim model files (`model/*.zginml`) are provided to allow reproduction of the companion manuscript.

## Contact

Corresponding author: Shantanu Gupta — shantanu.gupta@imd.ufrn.br
