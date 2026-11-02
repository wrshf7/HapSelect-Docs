# Function Reference

Quick reference for all exported HapSelect functions. Follow the link in each section heading for full parameter documentation and examples.

## Bundled Datasets

| Dataset | Description |
|---------|-------------|
| `HapSelect::map` | Example marker map (`SNP`, `Chromosome`, `Position`) |
| `HapSelect::geno` | Example dosage genotype matrix (localGEBV) |
| `HapSelect::geno_phased` | Phased genotype matrix of `geno` using the BEAGLE implementation (haplotypes) |
| `HapSelect::ld_pairs` | Pre-computed LD for example data using PLINK 1.9 |
| `HapSelect::marker_effects` | Example of marker effect estimates (`SNP`, `Effect`) |
| `HapSelect::BLUE` | Example file of BLUE inputs to calculate marker effects |

## [Configuration](configuration.md)

| Function | Description |
|----------|-------------|
| `set_plink_path(path)` | Use a specific PLINK 1.9 executable (or the folder containing it) for the rest of the session. `NULL` clears the setting |
| `set_beagle_path(path)` | Use a specific BEAGLE 5.5 `.jar` file (or the folder containing it) for the rest of the session. `NULL` clears the setting |
| `set_cpu_cores(cores)` | Limit the number of CPU cores HapSelect may use for the rest of the session. `NULL` clears the limit |
| `cpu_cores()` | Report the number of CPU cores HapSelect will use |

See also [Progress Reporting](progress.md) for the progress bar and `HAPSELECT_MACHINE_LOGGING` settings.

## [File Processing, Manipulation, and Import](workflow/data-import.md)

| Function | Description |
|---------|-------------|
| `read_vcf_geno(path, prefer_ds, ploidy)` | Basic VCF reader to import a VCF file to dosage genotypes in HapSelect format |
| `read_vcf_phased(path, ploidy)` | Basic VCF reader to import a phased VCF file to phased alleles in HapSelect format |
| `beagle_impute(vcf_in, out_prefix, ref, map, extra_args)` | Wrapper function to call BEAGLE to read a VCF file, impute and phase missing genotypes, and write out a VCF file |
| `beagle_impute_geno(geno, ref, map, extra_args, prefer_ds)` | Wrapper function to call BEAGLE and impute missing dosage genotypes |
| `beagle_phase_geno(geno, ref, map, extra_args)` | Wrapper function to call BEAGLE to impute and phase dosage genotypes into HapSelect phased format |
| `center_genotypes(geno)` | Center genotype/phased genotype matrix by subtracting per-marker mean dosage |
| `order_map(map, verbose)` | Order the map file by chromosome and position, and convert chromosome labels to numbers ([Pairwise LD](workflow/pairwise-ld.md)) |
| `order_geno(geno, map, verbose)` | Order the genotype table the same as the map and subset to markers in both ([Pairwise LD](workflow/pairwise-ld.md)) |

## [LD Calculation](workflow/pairwise-ld.md)

| Function | Description |
|----------|-------------|
| `pairwise_ld(genotype_matrix, parallelize, window, min_r2, min_obs)` | Compute pairwise r² using the internal C++ engine |
| `plink_pairwise_ld(prefix, ld_window, ld_window_kb, ld_window_r2, extra_args)` | Compute pairwise r² from a PLINK binary fileset |
| `plink_pairwise_ld_geno(geno, ld_window, ld_window_kb, ld_window_r2, extra_args)` | Compute pairwise r² from a genotype data frame via PLINK |
| `format_plink_ld(ld_path, bim)` | Convert a PLINK `.ld` output file into HapSelect LD format |

## [Haploblock Definition](workflow/haploblocks.md)

| Function | Description |
|----------|-------------|
| `def_blocks(strategy, map, geno)` | Partition genome into haplotype blocks using one of three strategies. `geno` is only used by `graph_strategy()` |
| `ld_strategy(ld, method, tolerance, tol_reset, threshold, start, parallel)` | LD-based blocking strategy: grow blocks from seed marker pairs while LD stays high |
| `window_strategy(window, method)` | Window-based blocking strategy: fixed number of markers or fixed map distance per block |
| `graph_strategy(theta_core, theta_core_by_chr, theta_extend, theta_bridge, theta_refill, ld_min_r2, window_ld, window_core, window_extend, min_links, max_gap_snps, max_gap_markers, max_gap_position, min_block_snps, parallel)` | Graph-based blocking strategy: blocks from connected components of a local LD graph |
| `block_obj_to_df(block_obj, map)` | Convert block list to a tidy data frame |
| `block_summary(block_df)` | Compute summary statistics across all blocks |

## [Genomic Prediction](workflow/marker-effects.md)

| Function | Description |
|----------|-------------|
| `create_marker_effects_file(geno, BLUE, h2_method, ploidy, mean_impute)` | Fit marker effects and return an effects data frame (wraps `solve_marker_effects()`) |
| `solve_marker_effects(geno, BLUE, h2_method, ploidy, mean_impute)` | Solve marker effects via RR-BLUP and report variance components |
| `compute_prediction_accuracy(geno, marker_effects, BLUE, mean_impute)` | Correlate predicted GEBVs with observed phenotypes, using the model from `solve_marker_effects()` |
| `n_fold_cross_validation(geno, BLUE, nfold, h2_method, ploidy, mean_impute)` | Estimate prediction accuracy by n-fold cross-validation |
| `cross_validation(geno, BLUE, train_prop, fold, h2_method, ploidy, mean_impute)` | Estimate prediction accuracy by random train/validation splits |

## [LocalGEBV and Haplotype Effects](workflow/local-gebv.md)

| Function | Description |
|----------|-------------|
| `compute_local_GEBV(geno, marker_effects, haploblocks_df, set_missing_NA, mean_adjust, parallel, chunk_size)` | Estimate per-block breeding values (localGEBV) for each individual from a dosage genotype matrix |
| `compute_haplotype_effects(geno, marker_effects, haploblocks_df, set_missing_NA, mean_adjust, parallel, chunk_size)` | Estimate per-block haplotype effects for each chromosome from a phased genotype matrix |

## [Visualisations](workflow/visualizations.md)

| Function | Description |
|----------|-------------|
| `marker_effects_plot(marker_effects, chr, pos, colors)` | Manhattan-style plot of marker effects |
| `unique_haplo_effects_plot(haplo_obj, colors, pos_type)` / `unique_localGEBV_effects_plot(haplo_obj, colors, pos_type)` | Manhattan-style plot of the unique haplotype/localGEBV effects in each block |
| `block_variance_manhattan_plot(haplo_obj, colors, pos_type)` | Manhattan-style plot of haploblock variance |
| `local_gebv_block_var_funnel_plot(haplo_obj, mean_line, scale_colors)` / `haplo_block_var_funnel_plot(haplo_obj, mean_line, scale_colors)` | Funnel plot of localGEBV/haplotype effects vs. scaled block variance |
| `block_var_funnel_plot(haplo_obj, mean_line, scale_colors)` | Generic funnel plot, maintained for backwards compatibility |
| `plot_haploblocks(haploblock_df, block_fill, chrom_fill, height, single_width_bp)` | Haploblock boundaries across chromosomes |
| `plot_haploblocks_variance(haploblock_obj, variance, scale_viridis_option, chrom_fill, height, single_width_bp)` | Haploblock boundaries across chromosomes, coloured by block variance |
| `plot_marker_density(map, bin_size, height, chrom_fill, col_low, col_mid, col_high)` | Marker density along chromosomes |
| `plot_ld_decay(map, ld, max_kb, point_color, curve_color, alpha, span, method, k)` | Smoothed LD decay curve |

## [Parent Selection](workflow/parent-selection.md) & [Basic Simulation](workflow/basic-simulation.md)

| Function | Description |
|----------|-------------|
| `select_top_blocks(haploblock_obj, n, perc_total, perc_of_total_var)` | Subset haploblock object to top blocks by variance |
| `local_gebv_parent_selection(haploblock_obj, strategy, n_founders, popSize, maxiter, run, pmutation, pcrossover, maximize, monitor)` | GA-based founder selection for localGEBV |
| `haplotype_parent_selection(haploblock_obj, strategy, n_founders, popSize, maxiter, run, pmutation, pcrossover, maximize, monitor)` | GA-based founder selection for haplotypes |
| `localGEBV_vs_TS_simulation(GA_output, geno, marker_effects, map, genetic_map_position, num_gen, num_sim_reps, num_cross_per_gen, num_TS_parents, mean_adjust, maximize, max_cM_chr, PCA, colors, alpha)` | Forward simulation comparing GA vs. truncation-selected parents |
| `Haplotype_vs_TS_simulation(GA_output, geno_phased, marker_effects, map, genetic_map_position, num_gen, num_sim_reps, num_cross_per_gen, num_TS_parents, mean_adjust, maximize, max_cM_chr, PCA, colors, alpha)` | Forward simulation comparing GA vs. truncation-selected parents |
