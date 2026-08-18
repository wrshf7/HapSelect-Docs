# Changelog

All notable changes to HapSelect are documented here. This project follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

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