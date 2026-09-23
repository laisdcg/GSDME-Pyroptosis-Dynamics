# GSDME systems-oncology project

## Complete reproducibility guide for the Boolean network, GEO single-cell integration, seed consensus and in silico epistasis

This README is the definitive documentation for the GSDME systems-oncology
project. It supersedes the separate pipeline and epistasis notes and describes
the complete workflow used to generate the final computational results.

The project integrates:

1. a literature-informed logical model created in GINsim;
2. asynchronous Boolean simulations under persistent DNA-damage response;
3. sustained perturbation screening;
4. structural control and minimum driver-set search;
5. probabilistic integration of liver-cancer single-cell RNA-seq;
6. tabular Q-learning for sequential intervention discovery;
7. consensus analysis across five computational random seeds;
8. matched in silico epistasis and rescue simulations.

All publication figures are exported as vector PDF and 600-dpi PNG.

---

## 1. Scientific question

The central model-based question is:

> Does the MALAT1/miR-204-5p/SIRT1 axis control access to the
> p53–mitochondrial–CASP3 cascade, while GSDME availability determines whether
> the terminal cell-death outcome is GSDME-mediated pyroptosis rather than
> apoptosis or resistance?

The proposed regulatory hierarchy is:

`MALAT1 -| miR-204-5p -| SIRT1 -| active p53 -> BAX -> CASP3 -> GSDME -> pyroptosis`

The project tests the internal dynamic consequences of this hypothesis and its
robustness across initial states, random seeds and single-cell-informed states.
It does not substitute for experimental perturbation and rescue in cultured
cells or in vivo models.

---

## 2. Correct scientific scope

The following language is appropriate:

- model-predicted regulatory hierarchy;
- in silico perturbation;
- in silico epistasis and rescue;
- model-projected cell-fate potential;
- observational single-cell concordance;
- candidate driver node or candidate control intervention;
- hypothesis-generating sequential intervention;
- computational robustness across random seeds.

The following claims are not supported by this computational study alone:

- biological proof of causality;
- clinical validation;
- experimental confirmation of knockout or overexpression;
- direct single-cell measurement of GSDME cleavage or caspase activity;
- universal control of every possible state of the Boolean network;
- therapeutic recommendation for patients.

GSDME is interpreted primarily as a **terminal execution gate** for the mode of
cell death. The analysis does not establish GSDME as an upstream master
regulator of the complete network.

---

## 3. Scripts and versions

| File | Version | Function |
|---|---:|---|
| `GSDME_systems_oncology_pipeline.R` | 1.2.1 | Executes the complete pipeline for one random seed |
| `GSDME_consensus_across_seeds.R` | 1.1.0 | Combines five complete seed folders after quality control |
| `GSDME_in_silico_epistasis.R` | 1.0.0 | Performs matched single/double perturbations and rescue analysis |
| `modelo_GINsim_GSDME_available.zginml` | final model | GINsim logical-network source |
| `README.md` | current | Complete execution and interpretation guide |

Do not rename node IDs inside the GINsim model unless every logical rule,
mapping and intervention that references those IDs is also updated.

---

## 4. Final logical model

The imported GINsim model contains:

- 27 nodes;
- 39 regulatory interactions;
- one formal input: `GSDME_availability`;
- one persistent experimental condition: `DDR_fixed_ON = 1`.

`DDR_fixed_ON` is not treated as an additional variable input in the final
interpretation. It is a fixed background condition applied to all therapeutic
perturbations.

### Logical rules

| Node ID | Logical rule or status |
|---|---|
| `GSDME_availability` | unique input |
| `GSDME_N` | `CASP3 & GSDME_availability` |
| `CASP3` | `CASP9` |
| `CASP9` | `CYTOCHROME_C` |
| `CYTOCHROME_C` | `MOMP` |
| `MOMP` | `BAX` |
| `BAX` | `MITO_DAMAGE & PUMA & !BCL2` |
| `BCL2` | `!PUMA` |
| `PUMA` | `p53_ACTIVE` |
| `p53_ACTIVE` | `DDR_fixed_ON & !SIRT1` |
| `miR_204_5p` | `!lncRNA_MALAT1` |
| `lncRNA_MALAT1` | `!p53_ACTIVE` |
| `ROS` | `DDR_fixed_ON` |
| `MITO_DAMAGE` | `ROS & !PGC1A` |
| `PGC1A` | `SIRT1` |
| `p21` | `p53_ACTIVE` |
| `CyclinD_CDK46` | `!p21` |
| `RB1` | `!CyclinD_CDK46` |
| `E2F1` | `!RB1` |
| `DDR_fixed_ON` | fixed to `1` during therapeutic analyses |
| `PYROPTOSIS_GSDME` | `GSDME_N` |
| `SURVIVAL` | `SIRT1 & !CASP3` |
| `RESISTANCE` | `DDR_fixed_ON & SURVIVAL` |
| `PROLIFERATION` | `E2F1 & SURVIVAL & !DDR_fixed_ON` |
| `APOPTOSIS` | `CASP3 & !GSDME_N` |
| `CELL_CYCLE_ARREST` | `p21 & RB1 & !E2F1 & !CASP3 & !GSDME_N` |
| `SIRT1` | `!miR_204_5p` |

These rules are model assumptions supported to different degrees by the
literature. Simulation results that follow directly from an encoded rule must
be interpreted as consequences of the model, not as independent experimental
verification of that rule.

---

## 5. Required software

- R 4.2 or later;
- RStudio is recommended but not required;
- internet connection for the first package installation and first GEO
  download;
- sufficient memory to read the compressed GSE140228 Droplet matrix.

### Required R packages

```r
install.packages(c(
  "xml2",
  "igraph",
  "ggplot2",
  "Matrix",
  "data.table",
  "patchwork"
))
```

The scripts install missing packages automatically unless `--no-install` is
provided.

---

## 6. Recommended project organization

```text
GSDME_Oncology_Pipeline/
├── GSDME_systems_oncology_pipeline.R
├── GSDME_consensus_across_seeds.R
├── GSDME_in_silico_epistasis.R
├── modelo_GINsim_GSDME_available.zginml
├── README.md
├── GEO_scRNA_data/
│   ├── GSE125449/
│   └── GSE140228/
├── resultados_finais_GSDME_seed_101/
├── resultados_finais_GSDME_seed_204/
├── resultados_finais_GSDME_seed_307/
├── resultados_finais_GSDME_seed_509/
├── resultados_finais_GSDME_seed_811/
├── resultados_consenso_GSDME/
└── resultados_epistasia_GSDME/
```

Each output directory contains `figures`, `tables` and `logs` subdirectories.
Always use a new output directory for a new analysis to avoid mixing files from
different versions or quick tests.

---

## 7. Set the working directory in RStudio

Open RStudio and set the working directory to the folder containing the three
scripts and the GINsim model:

```r
setwd("/home/usuario/Downloads/Artigo com Shantanu/GSDME_Oncology_Pipeline")
```

Change the path to the real location on the computer. Confirm it with:

```r
getwd()
list.files()
```

Create the path to the `Rscript` executable:

```r
rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
```

---

## 8. Optional diagnostic run

Use `--quick` only to confirm that installation, model parsing, GEO reading and
figure generation work:

```r
status_test <- system2(
  rscript,
  args = c(
    "GSDME_systems_oncology_pipeline.R",
    "--model", "modelo_GINsim_GSDME_available.zginml",
    "--geo", "GSE125449,GSE140228",
    "--geo-platform", "droplet",
    "--geo-max-cells", "500",
    "--geo-dir", "GEO_scRNA_data",
    "--out", "teste_GEO_GSDME",
    "--no-demo",
    "--quick",
    "--seed", "204"
  )
)

status_test
```

`0` indicates successful completion. Never report quick-mode numerical values
in the manuscript.

---

## 9. Final pipeline executions: five seeds

The main pipeline processes one seed per execution. The final analysis uses the
same model, datasets, cell limit and parameters for all five seeds.

Define a helper function in the RStudio Console:

```r
run_final_seed <- function(seed) {
  output_directory <- paste0("resultados_finais_GSDME_seed_", seed)

  system2(
    rscript,
    args = c(
      "GSDME_systems_oncology_pipeline.R",
      "--model", "modelo_GINsim_GSDME_available.zginml",
      "--geo", "GSE125449,GSE140228",
      "--geo-platform", "droplet",
      "--geo-max-cells", "2000",
      "--geo-dir", "GEO_scRNA_data",
      "--out", output_directory,
      "--no-demo",
      "--seed", as.character(seed)
    )
  )
}
```

Execute and verify each seed separately:

```r
status_101 <- run_final_seed(101)
status_101
```

```r
status_204 <- run_final_seed(204)
status_204
```

```r
status_307 <- run_final_seed(307)
status_307
```

```r
status_509 <- run_final_seed(509)
status_509
```

```r
status_811 <- run_final_seed(811)
status_811
```

All five statuses must be `0`.

Confirm together:

```r
c(
  seed_101 = status_101,
  seed_204 = status_204,
  seed_307 = status_307,
  seed_509 = status_509,
  seed_811 = status_811
)
```

The commands intentionally omit `--quick`.

---

## 10. Consensus across seeds

The consensus script rejects incomplete runs, possible quick runs, different
pipeline versions, model mismatches and folder/manifest seed mismatches.

Run:

```r
status_consensus <- system2(
  rscript,
  args = c(
    "GSDME_consensus_across_seeds.R",
    "--parent", ".",
    "--folders",
    paste(
      c(
        "resultados_finais_GSDME_seed_101",
        "resultados_finais_GSDME_seed_204",
        "resultados_finais_GSDME_seed_307",
        "resultados_finais_GSDME_seed_509",
        "resultados_finais_GSDME_seed_811"
      ),
      collapse = ","
    ),
    "--out", "resultados_consenso_GSDME"
  )
)

status_consensus
```

The status must be `0`. Then open:

```text
resultados_consenso_GSDME/tables/01_run_quality_control.csv
```

Every row must have `quality_control_pass = TRUE`.

### Meaning of the five seeds

The seeds control computational randomness in:

- sampled initial network states;
- asynchronous node-update order;
- cell subsampling where applicable;
- stochastic Q-learning exploration.

They measure algorithmic sensitivity and reproducibility. They are not
patients, cell lines, biological samples or biological replicates.

The final manuscript must report consensus means and variability across seeds,
not select the most visually favorable seed.

---

## 11. Final in silico epistasis and rescue analysis

The epistasis script evaluates the model-predicted order of the axis with
matched single and double interventions. Unlike the main pipeline, this script
runs all five seeds internally.

### Diagnostic epistasis test

```r
status_epistasis_test <- system2(
  rscript,
  args = c(
    "GSDME_in_silico_epistasis.R",
    "--model", "modelo_GINsim_GSDME_available.zginml",
    "--seeds", "101,204,307,509,811",
    "--trajectories", "60",
    "--geo-dir", "GEO_scRNA_data",
    "--geo-max-cells", "200",
    "--out", "teste_epistasia_GSDME",
    "--quick"
  )
)

status_epistasis_test
```

### Final epistasis analysis

```r
status_epistasis <- system2(
  rscript,
  args = c(
    "GSDME_in_silico_epistasis.R",
    "--model", "modelo_GINsim_GSDME_available.zginml",
    "--seeds", "101,204,307,509,811",
    "--trajectories", "500",
    "--max-steps", "350",
    "--geo-dir", "GEO_scRNA_data",
    "--geo-max-cells", "2000",
    "--geo-cell-draws", "1",
    "--geo-subset-seed", "204",
    "--out", "resultados_epistasia_GSDME"
  )
)

status_epistasis
```

The final status must be `0`, and
`resultados_epistasia_GSDME/tables/17_analysis_manifest.csv` must show:

- `quick_mode = FALSE`;
- all five seeds;
- 500 trajectories per condition and initial-state stratum;
- the MD5 checksum of the final GINsim model;
- GEO enabled, if the single-cell module was used.

---

## 12. Epistasis design

The analysis includes 19 reference, single and double intervention conditions.
The central rescue comparisons are:

| Upstream perturbation | Downstream perturbation | Double perturbation | Model question |
|---|---|---|---|
| MALAT1 OFF | miR-204-5p OFF | MALAT1 OFF + miR-204-5p OFF | Does miR-204-5p inhibition rescue MALAT1 inhibition? |
| miR-204-5p ON | SIRT1 ON | miR-204-5p ON + SIRT1 ON | Does SIRT1 restoration rescue miR-204-5p activation? |
| SIRT1 OFF | active p53 OFF | SIRT1 OFF + active p53 OFF | Does p53 inhibition block the SIRT1-loss phenotype? |
| active p53 ON | CASP3 OFF | active p53 ON + CASP3 OFF | Is CASP3 required downstream of p53? |
| CASP3 ON | GSDME unavailable | CASP3 ON + GSDME unavailable | Does GSDME determine the terminal mode of death? |
| MALAT1 OFF | GSDME unavailable | MALAT1 OFF + GSDME unavailable | Is MALAT1-driven pyroptosis conditional on GSDME? |

### Initial-state strata

The network analysis is conducted in two state spaces:

1. `Global random`: Bernoulli initial states sampling the broad network space.
2. `Resistant-local`: states generated around the resistant reference attractor
   with an 8% node-flip probability.

The same initial states and reproducible asynchronous random streams are used
for matched single/double comparisons within each seed. This reduces
differences caused only by unmatched random initialization.

### Rescue fraction

For pyroptosis probability `P`:

```text
rescue fraction =
  [P(upstream) − P(double)] /
  [P(upstream) − P(reference)]
```

Interpretation:

- approximately 1: complete model-based rescue;
- between 0.5 and 1: partial rescue;
- approximately 0: no rescue;
- below 0: the combined perturbation intensified the phenotype;
- missing: the upstream perturbation had no measurable effect, making the
  denominator zero.

### Downstream dominance

The double perturbation is compared with both single perturbations. A high
downstream-dominance score means that the double-perturbation phenotype is
closer to the downstream single perturbation than to the upstream single
perturbation. This is model-based epistatic ordering, not experimental genetic
epistasis.

### GSDME terminal-gate metrics

The script separately calculates:

```text
pyroptosis drop = P(pyroptosis | GSDME ON) − P(pyroptosis | GSDME OFF)

apoptosis gain = P(apoptosis | GSDME OFF) − P(apoptosis | GSDME ON)

resistance gain = P(resistance | GSDME OFF) − P(resistance | GSDME ON)
```

These quantities distinguish suppression of pyroptosis from redirection toward
apoptosis or resistance.

---

## 13. Main pipeline modules

### 13.1 GINsim import and validation

The `.zginml` archive is opened, the regulatory graph is located and node IDs,
labels, rules, input status and signed edges are imported. The parser checks
duplicated or missing IDs, unknown nodes in edges and unknown variables in
rules.

### 13.2 Asynchronous Boolean dynamics

At every step, the target value of every node is calculated from the current
state. One unstable, non-clamped node is selected for update. The simulation
stops at a fixed point, a repeated state or the maximum step limit.

Endpoint frequencies therefore represent the fraction of sampled stochastic
asynchronous trajectories ending with each node or fate active.

### 13.3 Perturbation screen

All curated perturbations are applied under `DDR_fixed_ON = 1`:

| Perturbation | Model purpose |
|---|---|
| GSDME loss | Test pyroptosis-to-apoptosis/resistance redirection |
| CASP3 inhibition | Test dependence of GSDME cleavage on CASP3 |
| miR-204-5p activation/inhibition | Test regulation of SIRT1 |
| MALAT1 activation/inhibition | Test upstream non-coding-RNA control |
| SIRT1 activation/inhibition | Test survival/resistance versus mitochondrial death |
| BCL2 inhibition | Test removal of BAX inhibition |

The literature evidence and its scope are kept separate from simulated
frequencies in `04_scientific_perturbation_evidence.csv`.

### 13.4 Structural control and minimum driver sets

The pipeline calculates in-degree, out-degree, betweenness, PageRank and
strongly connected components. Candidate sustained interventions are evaluated
from the resistant state and locally perturbed states.

An intervention is considered robust in the tested search when it reaches:

```text
PYROPTOSIS_GSDME = 1
RESISTANCE = 0
PROLIFERATION = 0
```

in at least 95% of sampled trajectories.

This is a minimum set within the searched intervention space, not an exhaustive
proof over all `2^27 = 134,217,728` possible Boolean states.

### 13.5 Single-cell binarization

Counts are library-size normalized and log-transformed. A two-component
Gaussian mixture estimates the posterior probability that each mapped gene is
in a high-expression state. Posterior probabilities are sampled to create
cell-specific initial Boolean states.

This is a probabilistic transcriptional initialization. It does not measure:

- active p53 protein;
- mitochondrial outer-membrane permeabilization;
- enzymatically active caspases;
- GSDME-N protein;
- membrane pore formation.

### 13.6 Observational concordance

For each GEO experiment, the pipeline compares observed single-cell node
probabilities with simulated perturbation endpoints:

```text
weighted concordance =
  1 − sum[weight × |observed probability − simulated frequency|] / sum(weights)
```

Weights are:

- 4 for a directly perturbed node;
- 2 for first-order regulatory targets;
- 1 for other measurable network nodes.

A high score means state similarity. It does not mean that the perturbation was
performed in the GEO study.

### 13.7 Reinforcement learning

The Q-learning agent may select:

- no intervention;
- MALAT1 inhibition;
- miR-204-5p activation;
- SIRT1 inhibition;
- PGC-1alpha inhibition;
- BCL2 inhibition.

The largest reward is assigned only when `PYROPTOSIS_GSDME = 1`. Intermediate
rewards are assigned to GSDME-N, CASP3 and BAX, with penalties for resistance,
survival, proliferation, apoptosis without pyroptosis and intervention cost.

The learned sequence is a temporal model hypothesis and not a treatment
recommendation.

---

## 14. GEO datasets and permitted roles

| Dataset | Analysed content | Role | Limitation |
|---|---|---|---|
| GSE125449 | Cells annotated as malignant; Set1 and Set2 kept separate | Tumour-state initialization and observational concordance | No matching perturbation arms; scRNA-seq does not measure active proteins or GSDME cleavage |
| GSE140228 | Sorted CD45+ immune cells from tumour and other tissues | Immune-microenvironment context | Not a direct tumour-cell validation dataset |

GSE125449 and GSE140228 are never pooled as if they were replicates of the same
cellular compartment.

The epistasis script uses only malignant cells from GSE125449 for cell-level
response projection. GSE140228 remains in the main pipeline as immune context.

### Recognized gene aliases

| Logical node | Gene symbols or aliases |
|---|---|
| `GSDME_availability` | GSDME, DFNA5 |
| `CYTOCHROME_C` | CYCS, CYTOCHROME_C |
| `PUMA` | BBC3, PUMA |
| `p53_ACTIVE` | TP53, P53 |
| `miR_204_5p` | MIR204, MIR204-5P, MIR_204_5P |
| `PGC1A` | PPARGC1A, PGC1A |
| `p21` | CDKN1A, P21 |
| `CyclinD_CDK46` | CCND1, CDK4, CDK6 |

Mature miR-204-5p is generally not captured reliably by standard scRNA-seq.
Missing molecular activities are not presented as directly observed evidence.

---

## 15. Study overview

> We developed a literature-informed Boolean network to examine how the
> MALAT1/miR-204-5p/SIRT1 axis may regulate access to a mitochondrial
> CASP3-dependent death cascade under persistent DNA-damage response. The model
> was evaluated through asynchronous perturbation screening, attractor-control
> analysis, single-cell-informed state initialization, seed-consensus analysis,
> reinforcement learning and matched in silico epistasis/rescue simulations.
> GSE125449 provided malignant-cell observational states, whereas GSE140228 was
> restricted to immune-microenvironmental context. The resulting predictions
> prioritize a GSDME-dependent terminal-death gate for future experimental
> validation and do not constitute biological or clinical proof.

---

## 16. Key biological and methodological references

1. Jiang G, Wen L, Zheng H, Jian Z, Deng W. miR-204-5p targeting SIRT1
   regulates hepatocellular carcinoma progression. *Cell Biochem Funct.*
   2016;34:505–510. PMID: 27748572. DOI: 10.1002/cbf.3223.
2. Hou Z, Xu X, Zhou L, et al. The long non-coding RNA MALAT1 promotes the
   migration and invasion of hepatocellular carcinoma by sponging miR-204 and
   releasing SIRT1. *Tumour Biol.* 2017;39:1010428317718135. PMID: 28720061.
   DOI: 10.1177/1010428317718135.
3. Wang Y, Gao W, Shi X, et al. Chemotherapy drugs induce pyroptosis through
   caspase-3 cleavage of a gasdermin. *Nature.* 2017;547:99–103. PMID: 28459430.
   DOI: 10.1038/nature22393.
4. Liu D, et al. SIRT1 inhibition-induced mitochondrial damage promotes
   GSDME-dependent pyroptosis in hepatocellular carcinoma cells. *Molecular
   Biotechnology.* 2024. PMID: 38044396. DOI: 10.1007/s12033-023-00964-z.
5. Sun X, et al. Germacrone induces caspase-3/GSDME activation and enhances ROS
   production, causing HepG2 pyroptosis. *Exp Ther Med.* 2022;24:456.
   PMID: 35747157. DOI: 10.3892/etm.2022.11383.
6. Ma L, Hernandez MO, Zhao Y, et al. Tumor cell biodiversity drives
   microenvironmental reprogramming in liver cancer. *Cancer Cell.*
   2019;36:418–430.e6. PMID: 31588021. DOI: 10.1016/j.ccell.2019.08.007.
   GEO: GSE125449.
7. Zhang Q, He Y, Luo N, et al. Landscape and dynamics of single immune cells
   in hepatocellular carcinoma. *Cell.* 2019;179:829–845.e20.
   PMID: 31675496. DOI: 10.1016/j.cell.2019.10.003. GEO: GSE140228.
8. Hou W, Tamura T, Ching WK, Akutsu T. Finding and analyzing the minimum set
   of driver nodes in control of Boolean networks. *Advances in Complex
   Systems.* 2016;19:1650006. DOI: 10.1142/S0219525916500065.
9. Yang G, Gómez Tejeda Zañudo J, Albert R. Target control in logical models
   using the domain of influence of nodes. *Front Physiol.* 2018;9:454.
   DOI: 10.3389/fphys.2018.00454.
10. Magaña-López G, et al. scBoolSeq: linking scRNA-seq statistics and Boolean
    dynamics. *PLoS Comput Biol.* 2024;20:e1011620.
    DOI: 10.1371/journal.pcbi.1011620.
