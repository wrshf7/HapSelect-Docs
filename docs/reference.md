# Function Reference

Quick reference for all exported HapSelect functions.

## Map & LD

| Function | Description |
|----------|-------------|
| `order_map(map)` | Sort a marker map by chromosome and position |
| `pairwise_ld(genotype_matrix, parallelize)` | Compute pairwise r² using the internal C++ engine |
| `plink_pairwise_ld(prefix, ld_window, ld_window_kb, ld_window_r2, extra_args)` | Compute pairwise r² from a PLINK binary fileset |
| `plink_pairwise_ld_geno(geno, ld_window, ld_window_kb, ld_window_r2, extra_args)` | Compute pairwise r² from a genotype data frame via PLINK |
| `format_plink_ld(ld_path, bim)` | Convert a PLINK `.ld` output file into HapSelect LD format |

## Haploblock Definition

| Function | Description |
|----------|-------------|
| `def_blocks(ld, map, method, tolerance, tol_reset, threshold, start, parallel)` | Partition genome into haplotype blocks |
| `block_obj_to_df(block_obj, map)` | Convert block list to a tidy data frame |
| `block_summary(block_df)` | Compute summary statistics across all blocks |

## Genomic Prediction

| Function | Description |
|----------|-------------|
| `create_marker_effects_file(geno, BLUE, h2_method, ploidy)` | Fit marker effects and return an effects data frame (wraps `solve_marker_effects`) |
| `solve_marker_effects(geno, BLUE, h2_method, ploidy)` | Solve marker effects via RR-BLUP and report variance components |
| `compute_prediction_accuracy(geno, marker_effects, BLUE)` | Correlate predicted GEBVs with observed phenotypes |
| `n_fold_cross_validation(geno, BLUE, nfold, h2_method, ploidy)` | Estimate prediction accuracy by n-fold cross-validation |
| `cross_validation(geno, BLUE, train_prop, fold, h2_method, ploidy)` | Estimate prediction accuracy by random train/validation splits |

## Local GEBV

| Function | Description |
|----------|-------------|
| `compute_local_GEBV(geno, marker_effects, haploblocks_df, marker_pecov, set_missing_NA, mean_adjust, parallel, chunk_size)` | Estimate per-block breeding values for each individual |
| `center_genotypes(geno)` | Center genotype matrix by subtracting per-marker mean dosage |

## Haploblock Variance Test

| Function | Description |
|----------|-------------|
| `haploblock_var_test(haploblock_obj, geno, gen_var, threshold)` | Chi-square test of block variance against a null expectation; adds p-values and FDR to the haploblock table |

## Visualisations

| Function | Description |
|----------|-------------|
| `marker_effects_plot(marker_effects, chr, pos, colors)` | Manhattan plot of marker effects |
| `unique_haplo_effects_plot(haplo_obj, colors, pos_type)` | Distribution of unique haplotype effects per block |
| `block_var_funnel_plot(haplo_obj, mean_line, scale_colors)` | Block variance vs. block size |
| `plot_haploblocks(haploblock_df, block_fill, chrom_fill, height, single_width_bp)` | Haploblock boundaries across chromosomes |
| `plot_marker_density(map, bin_size, height, chrom_fill, col_low, col_mid, col_high)` | Marker density histogram |
| `plot_ld_decay(map, ld, max_kb, point_color, curve_color, alpha, span, method, k)` | Smoothed LD decay curve |

## Parent Selection & Basic Simulation

| Function | Description |
|----------|-------------|
| `select_top_blocks(haploblock_obj, n, perc_total, perc_of_total_var)` | Subset haploblock object to top blocks by variance |
| `genetic_algorithm(localGEBV, n_founders, popSize, maxiter, run, selfing, pmutation, pcrossover, pelite)` | GA-based founder selection |
| `GA_vs_TS_simulation(GA_output, geno, marker_effects, map, genetic_map_position, num_gen, num_sim_reps, num_cross_per_gen, num_TS_parents, mean_adjust, max_cM_chr, PCA, colors, alpha)` | Forward simulation comparing GA vs. truncation-selected parents |

## Bundled Datasets

| Dataset | Description |
|---------|-------------|
| `map` | Example marker map |
| `geno` | Example genotype matrix |
| `ld_pairs` | Precomputed LD for example data |
| `marker_effects` | Example marker effect estimates |
| `marker_pecov` | Example marker prediction error covariance matrix |
