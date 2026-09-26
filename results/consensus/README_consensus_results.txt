GSDME CONSENSUS ACROSS COMPUTATIONAL SEEDS
==========================================

Seeds included: 101, 204, 307, 509, 811
Number of seeds: 5
Pipeline version: 1.6.1-HCC-GSDME-switch-DDR-ON

Interpretation:
- Means, SDs and 95% t intervals quantify Monte Carlo/seed variability.
- Seeds are computational replicates, not biological replicates.
- GEO concordance is observational support and is not causal validation.
- A driver/intervention should not be called robust from a single seed.
- Biological causality requires experimental perturbation and rescue.

Primary tables:
01_run_quality_control.csv
02_perturbation_consensus.csv
04_driver_set_consensus.csv
05_driver_action_frequency.csv
07_GEO_perturbation_concordance_consensus.csv
08_GEO_fate_consensus.csv
10_RL_target_consensus.csv
11_RL_action_stability.csv
13_seed_level_phenotype_endpoints.csv
15_GSDME_dependency_contrasts_by_seed.csv
16_GSDME_dependency_contrasts_consensus.csv
17_exact_stable_states_consensus_DDR_ON.csv
18_convergence_consensus_DDR_ON.csv
19_feedback_loop_consensus_DDR_ON.csv
20_stable_states_and_convergence_DDR_ON.csv

Figures:
- Figure 01 is copied after verifying an identical network model.
- Figures 02-07 summarize all included computational seeds.
- Supplementary Figure S1 shows pyroptosis endpoint uncertainty.

Publication note:
Report the exact seeds, full-run parameters, pipeline version, and the
distinction between computational robustness and biological evidence.
