# SIRT1–GSDME cell-death switch in hepatocellular carcinoma (HCC)

Boolean-network modeling (GINsim), Monte Carlo seed consensus, in silico
epistasis/rescue, and independent TCGA-LIHC expression validation for a
MALAT1 / miR-204-5p / SIRT1 / p53 / CASP3 / GSDME logical model of
pyroptotic competence in liver cancer.

> 📄 Manuscript: *Single-cell-informed Boolean epistasis predicts
> upstream control of GSDME-dependent pyroptotic competence in hepatocellular carcinoma (HCC)*.

---

## Repository structure

```
GSDME-pyroptosis-HCC/
├── model/
│   └── GINsim-miR_204_GSDME_Pyroptosis.zginml   # Boolean network (GINsim)
├── scripts/
│   ├── 01_GSDME_systems_oncology_pipeline.R           # per-seed Boolean simulation + GEO validation
│   ├── 02_GSDME_consensus_across_seeds.R              # combines =5 completed seed runs
│   ├── 03_GSDME_in_silico_epistasis.R                 # matched single/double perturbations + rescue
│   ├── 04_GSDME_TCGA_LIHC.R                            # independent TCGA-LIHC expression module
│   └── RUN_04_TCGA_LIHC_step_by_step.R                 # RStudio "run one block at a time" driver
├── results/                
│   ├── seeds/
|       ├── seed_101/
|       ├── seed_204/
|       ├── seed_307/
|       ├── seed_509/
|       ├── seed_811/
│   ├── consensus/
│   ├── epistasis/
│   └── tcga_lihc/
├── LICENSE
└── README.md
```

## Requirements

- R ≥ 4.2
- Packages: `xml2`, `igraph`, `ggplot2`, `scales`, `patchwork`, `Matrix`,
  `data.table`, `httr`, `jsonlite`
  (each script auto-installs what's missing unless run with `--no-install`)
- Internet access for the Boolean pipeline (`--geo ...` downloads from GEO/NCBI)
  and for the TCGA module (queries the GDC API)
- No network needed to re-run `tcga_extract()` / `tcga_analyze()` once GDC
  files are already cached locally

---

## Pipeline order

| # | Script | What it does | Depends on |
|---|--------|---------------|------------|
| 1 | `01_GSDME_systems_oncology_pipeline.R` | Loads the `.zginml`, computes reference attractors, screens literature-grounded perturbations, searches minimal driver sets, runs single-cell-informed simulations (demo + GEO), fits patient multi-omic digital twins, trains a Q-learning agent. Run **once per seed**. | model file |
| 2 | `02_GSDME_consensus_across_seeds.R` | Quality-controls and merges **≥ 5** completed (non-`--quick`) seed folders from step 1 into consensus tables/figures with Monte-Carlo 95% CIs. | ≥ 5 folders from step 1, same pipeline version & model |
| 3 | `03_GSDME_in_silico_epistasis.R` | Matched single/double Boolean perturbations along the MALAT1→miR-204-5p→SIRT1→p53→CASP3→GSDME axis; rescue fractions, GSDME terminal-gate test, GSE125449-informed epistasis. | model file (independent of steps 1–2) |
| 4 | `04_GSDME_TCGA_LIHC.R` | Independent, read-only validation: pulls matched tumor/adjacent-normal TCGA-LIHC RNA-seq (STAR counts) + mature miR-204-5p isomiRs from GDC, paired Wilcoxon tests, Spearman correlations, publication figures. Never touches the Boolean simulation. | model file only (does not read steps 1–3 outputs) |

Steps 1→2 must use the **same model file and pipeline version** — the
consensus script refuses to merge mismatched runs. Step 3 and step 4 are
independent of steps 1–2 and of each other.

---

## Passo a passo completo (execução no R)

> **DDR = 1 (ON) em toda esta versão.** Todas as análises do modelo — estados
> estáveis, simulações assíncronas, consenso e epistasia — usam DDR fixado em
> 1. DDR não foi medido nos dados single-cell: para HCC, os scripts comparam
> o RNA observado com as saídas simuladas **nesse contexto fixo**.

O bloco abaixo assume que você está numa pasta de trabalho com todos os
arquivos `.R` e o `.zginml` juntos.

```r
setwd("/home/usuario/Downloads/GSDME-Pyroptosis-Dynamics")
rscript <- file.path(R.home("bin"), "Rscript")
model_path <- "GINsim-miR_204_GSDME_Pyroptosis.zginml"
stopifnot(file.exists(model_path))

run_seed <- function(seed) {
  status <- system2(rscript, args = c(
    "GSDME_systems_oncology_pipeline.R",
    "--model", shQuote(model_path),
    "--geo", "GSE125449,GSE189903",
    "--geo-platform", "droplet",
    "--geo-max-cells", "2000",
    "--geo-dir", "GEO_scRNA_data",
    "--out", paste0("resultados_finais_GSDME_seed_", seed),
    "--no-demo", "--seed", as.character(seed)
  ))
  if (status != 0L) stop("Falha na semente ", seed, ": status ", status)
  status
}
```

Execute cada linha separadamente e confira o retorno `0`:

```r
status_101 <- run_seed(101)
status_204 <- run_seed(204)
status_307 <- run_seed(307)
status_509 <- run_seed(509)
status_811 <- run_seed(811)
```

Depois, calcule o consenso a partir das cinco pastas geradas:

```r
status_consenso <- system2(rscript, "GSDME_consensus_across_seeds.R")
stopifnot(status_consenso == 0L)
```

Por fim, rode a epistasia separadamente (independente dos passos 1–2):

```r
status_epistasia <- system2(rscript, args = c(
  "GSDME_in_silico_epistasis.R",
  "--model", shQuote(model_path),
  "--seeds", "101,204,307,509,811",
  "--trajectories", "500",
  "--geo-dir", "GEO_scRNA_data",
  "--geo-max-cells", "2000",
  "--out", "resultados_epistasia_GSDME"
))
stopifnot(status_epistasia == 0L)
```

A validação TCGA-LIHC (passo 4) é independente dos passos 1–3 e roda em
outro fluxo, aberto no RStudio bloco por bloco — veja
[`scripts/RUN_04_TCGA_LIHC_step_by_step.R`](scripts/RUN_04_TCGA_LIHC_step_by_step.R).

### Figuras geradas

- O script das sementes (passo 1) gera `Figure_02c_GSDME_dependency_DDR_ON` e
  `Figure_08c_GSDME_dependency_exact_DDR_ON`.
- O consenso (passo 2) gera `Figure_02c_GSDME_dependency_five_phenotypes_consensus`
  e `Figure_02d_GSDME_dependency_contrast_consensus`.
- A epistasia (passo 3) gera `Figure_10_GSDME_dependency_DDR_ON_five_phenotypes`
  e, quando existem trajetórias elegíveis, `Figure_11_GSDME_KO_matched_trajectory_fates`.

Todas são salvas em PNG e PDF dentro da pasta `figures` correspondente a
cada `--out`.

---

## Interpretation guardrails (read before citing a number)

- **Seeds are computational replicates, not biological replicates.** Consensus
  CIs quantify Monte Carlo variability only.
- **GEO/TCGA support is observational**, never causal validation of a
  knockout/overexpression — no dataset here received the modeled perturbations
  experimentally.
- **DDR is fixed ON** in every simulation and stable-state analysis in this
  version; it was not measured in the single-cell/TCGA data.
- No 50/50 apoptosis/pyroptosis split is imposed anywhere: `PYROPTOSIS = GSDME`
  and `GSDME = CASP3` in the attached model, and with CASP3 E1 and endogenous
  GSDME the model decides the fate frequencies from its own logic. Forcing
  `DFNA5 = 1` with `GSDME E1` is a logical control, not a measurement of
  cleavage.
- The GEO datasets used did **not** receive these perturbations
  experimentally. `PPARGC1A` is not part of the supplied model and was not
  added to the scripts.
- TCGA figures never plot a model node that has no measurable transcript
  (e.g. purely conceptual states, or `CDK4_6_CyclinD`, which is expanded into
  its individual components `CDK4`, `CDK6`, `CCND1`, `CCND2`, `CCND3` instead
  of being shown as one panel).
- `--quick` runs are for debugging only and must never be reported or mixed
  into the consensus.

---

## License

Code is released under the [MIT License](LICENSE). The GINsim model file
(`model/*.zginml`) is provided for reproducibility of the companion
manuscript;
