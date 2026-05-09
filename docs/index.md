# HapSelect

HapSelect is an R package for haplotype based genomic selection. It partitions the genome into linkage disequilibrium blocks (haploblocks), estimates the breeding value contribution of each block per individual (localGEBV or true haplotype effect), uses a genetic algorithm (GA) to select the set of parents that maximises coverage of high-value haplotype alleles, and performs basic simulation to compare GA parents to truncation selection (TS) parents performance over time.

The package is designed for breeding programs running genomic selection, where identifying and targeting the genomic regions with the most exploitable diversity can improve the efficiency of parent choice.

## Workflow at a glance

| Stage | What it does | Key function |
|-------|-------------|--------------|
| Pairwise LD | Compute r² between all marker pairs | `pairwise_ld()` / `plink_pairwise_ld()` |
| Haploblocking | Partition genome into LD-based haploblocks | `def_blocks()` |
| Genomic Prediction | Compute marker effects and prediction accuracy | `create_marker_effects_file()`, `n_fold_cross_validation()`, `cross_validation()` |
| LocalGEBV | Estimate per-block breeding value per individual | `compute_local_GEBV()` |
| Visualisation | Explore haploblock structure, localGEBV patterns, and more | `plot_haploblocks()`, `plot_ld_decay()`, … |
| Parent selection | Optimise founder set using a genetic algorithm | `genetic_algorithm()` |
| Basic Simulation | Compare GA and TS parent performance over time and explore how diversity is captured | `GA_vs_TS_simulation()` |

See the [Overview](overview.md) for a full explanation of the method, or jump to [Getting Started](getting-started.md) to run the example workflow.

## Installation (quick for experienced users)

Install with DevTools:

```r
install.packages("devtools")
devtools::install("path/to/HapSelect")
```

Install with base R command:

```r
install.packages("path/to/HapSelect", type = "source", repos = NULL)
```
## Authors

Will Shaffer, Zane Carter, Victor Papin — The University of Queensland
