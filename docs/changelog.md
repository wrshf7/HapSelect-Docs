# Changelog

All notable changes to HapSelect are documented here. This project follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

## [2.0.0] - 2026-10-07 

### Added
- graphing haploblocking strategy added to `def_blocks()` for alternative haploblocking strategies
- window-based haploblocking (physical/genetic distance based on the map or number of SNP) added to `def_blocks()` for alternative haploblocking strategies
- `strategy` argument with method-specific container function call added to `def_blocks()` to differentiate methodology used for haploblocking
- HapSelect R package text banner added on load
- Rudimentary imputation and phasing using a BEAGLE backbone added to simplify workflows: `beagle_impute_geno()` and `beagle_phase_geno()` (imputes and phases)
- VCF read and write functions added for file importation and exportation simplicity (Note: does not currently handle aneuploidy or variable ploidy)
- `write_vcf_geno()`, `read_vcf_geno()`, `read_vcf_phased()`
- Visualization for haploblock variance on physical/genetic chromosome location added to explore regional variance and block size: `plot_haploblocks_variance()`
- Visualization for GA optimization and convergence `GA_progress_plot()`
- Limitations to CPU handlers and machine-readable step reporting added: `set_cpu_cores()`
- Custom C++ LD computation added to remove reliance on PLINK backbone: `pairwise_ld()` with C++ backbone
- Functions to set PLINK executable/binary and BEAGLE .jar file locations added: `set_beagle_path()` and `set_plink_path()`
- Added more file and input checks, more informative error messages
- Actual wheat example input files added
- More extensive testthat coverage to ensure future updates produce viable code

### Fixed
- Fixed bug that produced `NA` prediction accuracy if any BLUE/adj phenotype/de-regressed BLUP inputs were `NA`
- Small reporting or file error bugs fixed (note nothing that would affect any results users see)
- Some conflicts between example input files and old input files that were not updated are now updated and correct

### Changed
- Step Reporting and progress reporting updated
- `def_blocks()` function syntax changed to support evolving haploblocking strategies
- detection of PLINK improved
- updated graphical progress bar calls to reduce RStudio GUI performance impacts
- `order_map()` will now coerce chromosomes to numerics and report
- many issues with files will now report inputs causing issues (e.g., markers missing from the map file)

## [1.0.1] - 2026-07-15

### Added
- `block_variance_manhattan_plot()` for visualizing block variance across chromosomes.

### Fixed
- `unique_localGEBV_effects_plot()` now correctly plots the localGEBV effect (`Haplotype_Effect`) instead of reusing the haplotype effect column.
- Clarified the PLINK-not-found error message on Windows to make it explicit that PLINK must be installed or the executable placed in the expected directory.

### Changed
- Updated README documentation and refreshed the overview diagram image.

## [1.0.0] - 2026-07-03

Initial release of HapSelect.