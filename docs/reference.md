# Function Reference

Quick reference for all exported HapSelect functions.

## Map & LD

| Function | Description |
|----------|-------------|
| `order_map(map)` | Sort a marker map by chromosome and position |
| `pairwise_ld(geno, parallelize)` | Compute pairwise r² using the internal C++ engine |
| `plink_pairwise_ld(prefix)` | Compute pairwise r² using PLINK v1.9 |

## Haploblock Definition

| Function | Description |
|----------|-------------|
| `def_blocks(ld, map, method, threshold, tolerance, tol_reset, start, parallel)` | Partition genome into haplotype blocks |
| `block_obj_to_df(block_obj, map)` | Convert block list to a tidy data frame |

## Local GEBV

| Function | Description |
|----------|-------------|
| `compute_local_GEBV(geno, marker_effects, haploblocks_df, set_missing_NA, center)` | Estimate per-block breeding values for each individual |

## Visualizations

| Function | Description |
|----------|-------------|
| `marker_effects_plot(marker_effects, chr, pos)` | Manhattan plot of marker effects |
| `unique_haplo_effects_plot(haplo_obj)` | Distribution of unique haplotype effects per block |
| `block_var_funnel_plot(haplo_obj, mean_line)` | Block variance vs. block size |
| `plot_haploblocks(haploblock_df)` | Haploblock boundaries across chromosomes |
| `plot_marker_density(map_df, bin_size)` | Marker density histogram |
| `plot_ld_decay(map, ld, max_kb, span, k, method)` | Smoothed LD decay curve |

## Parent Selection

| Function | Description |
|----------|-------------|
| `genetic_algorithm(localGEBV, n_founders, popSize, maxiter, run, selfing, pmutation, pcrossover, pelite)` | GA-based founder selection |

## Bundled Datasets

| Dataset | Description |
|---------|-------------|
| `map` | Example marker map |
| `geno` | Example genotype matrix |
| `pairwise_ld` | Precomputed LD for example data |
| `marker_effects` | Example marker effect estimates |
