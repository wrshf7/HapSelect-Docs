# HapSelect

HapSelect is an R package for haplotype-block-based genomic selection. It partitions the genome into linkage disequilibrium blocks, estimates the breeding value contribution of each block per individual (localGEBV), and uses a genetic algorithm to select the set of parents that maximises coverage of high-value haplotype alleles.

The package is designed for plant breeding programs running genomic selection, where identifying and targeting the genomic regions with the most exploitable diversity can improve the efficiency of parent choice.

## Workflow at a glance

| Stage | What it does | Key function |
|-------|-------------|--------------|
| Pairwise LD | Compute r² between all marker pairs | `pairwise_ld()` / `plink_pairwise_ld()` |
| Block definition | Partition genome into LD-based haploblocks | `def_blocks()` |
| Local GEBV | Estimate per-block breeding value per individual | `compute_local_GEBV()` |
| Visualisation | Explore haploblock structure and localGEBV patterns | `plot_haploblocks()`, `plot_ld_decay()`, … |
| Parent selection | Optimise founder set using a genetic algorithm | `genetic_algorithm()` |

See the [Overview](overview.md) for a full explanation of the method, or jump to [Getting Started](getting-started.md) to run the example workflow.

## Installation

```r
install.packages("devtools")
devtools::install("path/to/HapSelect")
```

## Authors

Will Shaffer, Zane Carter, Victor Papin — The University of Queensland
