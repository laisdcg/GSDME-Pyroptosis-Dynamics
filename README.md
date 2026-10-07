# SIRT1–GSDME cell-death switch in hepatocellular carcinoma (HCC)

Boolean-network modeling (GINsim), Monte Carlo seed consensus, in silico
epistasis/rescue, and independent TCGA-LIHC expression validation for a
MALAT1 / miR-204-5p / SIRT1 / p53 / CASP3 / GSDME logical model of
pyroptotic competence in liver cancer.

> 📄 Manuscript: *Boolean modeling and transcriptomic integration of the MALAT1/miR-204-5p/SIRT1 axis in the regulation of GSDME-dependent pyroptosis in hepatocellular carcinoma)*.

---

## Repository structure

```
GSDME-pyroptosis-HCC/
├── model/
│   └── GINsim-miR_204_GSDME_Pyroptosis.zginml   # Boolean network (GINsim)
├── scripts/
│   ├── bulkrna-seq-gsdme.R          
│   ├── tcga-lihc-gsdme.R              
├── results/                
├── LICENSE
└── README.md
```
## Overview

This repository contains the computational workflow used to investigate how upstream regulatory mechanisms influence the balance between proliferation, resistance, cell-cycle arrest, apoptosis, and pyroptosis in HCC.

The study integrates:

- a **GINsim Boolean regulatory model**;
- *in silico* perturbation analysis;
- three independent GEO cohorts: **GSE14520, GSE60502, and GSE121248**;
- **TCGA-LIHC** RNA-seq and miRNA-seq data;
- regulatory-edge evaluation using **Spearman correlation, NMI, and GGC**.

## Requeriments
```
GEOquery
Biobase
XML
xml2
httr
jsonlite
dplyr
tidyr
ggplot2
patchwork
igraph
infotheo
```

## Data availability

Public datasets:
- NCBI GEO: GSE14520, GSE60502, GSE121248
- NCI GDC: TCGA-LIHC

## Boolean model

The final GINsim network contains:

- **31 nodes**
- **68 regulatory interactions**
- **5 phenotypic outputs**:
  - PROLIFERATION
  - RESISTANCE
  - CELL-CYCLE ARREST
  - APOPTOSIS
  - PYROPTOSIS

DDR was maintained active (`DDR = 1`) during the simulations, representing persistent DNA-damage signaling.

The main regulatory axis investigated was:

```text
ATM ┤ CDC25A → E2F1 → MYC → MALAT1 ┤ miR-204-5p ┤ SIRT1 ┤ TP53 → CDKN1A ┤ CASP3 → GSDME

```
---

## License

Code is released under the [MIT License](LICENSE). The GINsim model file
(`model/*.zginml`) is provided for reproducibility of the companion
manuscript;
