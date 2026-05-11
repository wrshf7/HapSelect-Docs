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
| Basic Simulation | Compare GA and TS parent performance over time and explore how diversity is captured using the [genomicSimulation](https://github.com/vllrs/genomicSimulation) R package. | `GA_vs_TS_simulation()` |

See the [Overview](overview.md) for a full explanation of the method, or jump to [Getting Started](getting-started.md) to run the example workflow.

## Installation

Install with DevTools:

```r
install.packages("devtools")
devtools::install("path/to/HapSelect")
```

Install with base R command:

```r
install.packages("path/to/HapSelect", type = "source", repos = NULL)
```

!!! warning
    [genomicSimulation](https://github.com/vllrs/genomicSimulation), [PLINK 1.9](https://www.cog-genomics.org/plink/), and [RTools 4.5](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html) (for Windows users) dependencies must be installed separately. Alternatively run the helper installation scripts to install all needed software. The package is usable without PLINK 1.9, but related functions that call PLINK will fail with an error. See [Installation](./installation.md) for more details.

## Important Papers
<strong>localGEBV Method and Haploblock Formation:</strong><br>
[Shaffer et al. 2025. Local genomic estimates provide a powerful framework for haplotype discovery. bioRxiv (under review).](https://doi.org/10.1101/2025.08.28.672830)

<strong>Origin of the localGEBV Method and Parent Optimization with a Genetic Algorithm:</strong><br>
[Kemper et al. 2012. Long-term selection strategies for complex traits using high-density genetic markers. J Dairy Sci.](https://doi.org/10.3168/jds.2011-5289)

<strong>The First Implementation of Haploblocking with localGEBV:</strong><br>
[Voss-Fels et al. 2019. Breeding improves wheat productivity under contrasting agrochemical input levels. Nat Plants.](https://doi.org/10.1038/s41477-019-0445-5)

<strong>The Concept of the Ultimate Genotype:</strong><br>
[Hays et al. 2024. Potential approaches to create ultimate genotypes in crops and livestock. Nat Genet.](https://doi.org/10.1038/s41588-024-01942-0)

<strong>genomicSimulation Package Description and Simulation Comparing the Genetic Algorithm to Truncation Selection:</strong><br>
[Villiers et al. 2024. Evolutionary computing to assemble standing genetic diversity and achieve long-term genetic gain. Plant Genome.](https://doi.org/10.1002/tpg2.20467)

## Authors

Will Shaffer, Zane Carter, Victor Papin — The University of Queensland
