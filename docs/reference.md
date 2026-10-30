# Function Reference

Quick reference for all exported HapSelect functions.

## Bundled Datasets

| Dataset | Description |
|---------|-------------|
| `HapSelect::map` | Example marker map |
| `HapSelect::geno` | Example genotype matrix |
| `HapSelect::ld_pairs` | Pre-computed LD for example data using PLINK 1.9 |
| `HapSelect::geno_phased` | Phased genotype matrix (geno) using BEAGLE implementation |
| `HapSelect::marker_effects` | Example of marker effect estimates |
| `HapSelect::BLUE` | Example file of BLUE inputs to calculate marker effects |

## Environment Settings
| Function | Description |
|----------|-------------|
| `set_beagle_path()` | Sets the folder location of BEAGLE if not in `PATH` or the `/user/bin` folder for the entire session |
| `set_plink_path()` | Sets the folder location of PLINK if not in `PATH` or the `/user/bin` folder for the entire session |

## File Processing, Manipulation, and Import
| Function | Description |
|---------|-------------|
| `read_vcf_geno(path, prefer_ds, ploidy)` | Basic VCF reader to import a VCF file to dosage genotypes in HapSelect format |
| `read_vcf_phased(path, ploidy)` | Basic VCF reader to import a VCF file to phased alleles in HapSelect format |
| `beagle_impute(vcf_in, out_prefix, ref, map, extra_args)` | Wrapper function to call BEAGLE to read a VCF file, impute missing genotypes, and write out a vcf file |
| `beagle_impute_geno(geno, ref, map, extra_args, prefer_ds)` | Wrapper function to call BEAGLE and impute missing dosage genotypes |
| `beagle_phase_geno(geno, ref, map, extra_args)` | Wrapper function to call BEAGLE and impute dosage genotypes and phase to HapSelect phased format |
| `center_genotypes(geno)` | Center genotype/phased genotype matrix by subtracting per-marker mean dosage |
| `order_map(map, verbose)` | Function to order and format the map file |
| `order_geno(geno, map, verbose)` | Function to order the genotype table the same as the map and subset to mapped markers |

## LD Calculation

| Function | Description |
|----------|-------------|
| `pairwise_ld(genotype_matrix, parallelize, window, min_r2, min_obs)` | Compute pairwise r² using the internal C++ engine |
| `plink_pairwise_ld(prefix, ld_window, ld_window_kb, ld_window_r2, extra_args)` | Compute pairwise r² from a PLINK binary fileset |
| `plink_pairwise_ld_geno(geno, ld_window, ld_window_kb, ld_window_r2, extra_args)` | Compute pairwise r² from a genotype data frame via PLINK |
| `format_plink_ld(ld_path, bim)` | Convert a PLINK `.ld` output file into HapSelect LD format |

## Haploblock Definition

| Function | Description |
|----------|-------------|
| `def_blocks(map, strategy)` | Partition genome into haplotype blocks using one of three strategies |
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

## LocalGEBV and Haplotype Effects

| Function | Description |
|----------|-------------|
| `compute_local_GEBV(geno, marker_effects, haploblocks_df, set_missing_NA, mean_adjust, parallel, chunk_size)` | Estimate per-block breeding values for each individual |
| `compute_haplotype_effects(geno, marker_effects, haploblocks_df, set_missing_NA, mean_adjust, parallel, chunk_size)` | Estimate per-block breeding values for each individual |

## Visualisations

| Function | Description |
|----------|-------------|
| `marker_effects_plot(marker_effects, chr, pos, colors)` | Manhattan plot of marker effects |
| `unique_haplo_effects_plot(haplo_obj, colors, pos_type)` / `unique_localGEBV_effects_plot(haplo_obj, colors, pos_type)` | Distribution of unique haplotype/localGEBV effects per block |
| `local_gebv_block_var_funnel_plot(haplo_obj, mean_line, scale_colors)` / `haplo_block_var_funnel_plot(haplo_obj, mean_line, scale_colors)`| Block variance vs. block size |
| `plot_haploblocks(haploblock_df, block_fill, chrom_fill, height, single_width_bp)` | Haploblock boundaries across chromosomes |
| `plot_marker_density(map, bin_size, height, chrom_fill, col_low, col_mid, col_high)` | Marker density histogram |
| `plot_ld_decay(map, ld, max_kb, point_color, curve_color, alpha, span, method, k)` | Smoothed LD decay curve |
| `block_variance_manhattan_plot(haplo_obj, colors, pos_type)` | Manhattan-style plot of haploblock variance |
| `GA_progress_plot(parent_selection_object, max_color, mean_color, ribbon_color, ribbon_alpha, max_linewidth, mean_linewidth)` | Plot convergence and optimization curves of a GA run |
| `plot_haploblock_variance(haploblock_obj, variance, scale_viridis_option, chrom_fill, height, single_width_bp)` | Plot the haploblock boundaries and color according to variance |



## Parent Selection & Basic Simulation

| Function | Description |
|----------|-------------|
| `select_top_blocks(haploblock_obj, n, perc_total, perc_of_total_var)` | Subset haploblock object to top blocks by variance |
| `local_gebv_parent_selection(haploblock_obj, strategy, n_founders, popSize, maxiter, run, selfing, pmutation, pcrossover, maximize, monitor)` | GA-based founder selection for localGEBV |
| `haplotype_parent_selection(haploblock_obj, n_founders, popSize, maxiter, run, selfing, pmutation, pcrossover, maximize, monitor)` | GA-based founder selection for haplotypes |
| `localGEBV_vs_TS_simulation(GA_output, geno, marker_effects, map, genetic_map_position, num_gen, num_sim_reps, num_cross_per_gen, num_TS_parents, mean_adjust, max_cM_chr, PCA, colors, alpha)` | Forward simulation comparing GA vs. truncation-selected parents |
| `Haplotype_vs_TS_simulation(GA_output, geno_phased, marker_effects, map, genetic_map_position, num_gen, num_sim_reps, num_cross_per_gen, num_TS_parents, mean_adjust, max_cM_chr, PCA, colors, alpha)` | Forward simulation comparing GA vs. truncation-selected parents |

