# Installation

## System Dependencies

Some HapSelect functions (PLINK-backed LD) require external tools. Install them with the provided scripts before using the package.

=== "Linux (Ubuntu/Debian)"
    ```bash
    ./inst/scripts/install/install_linux.sh
    ```

=== "macOS"
    ```bash
    ./inst/scripts/install/install_mac.sh
    ```

=== "Windows"
    ```powershell
    Powershell -ExecutionPolicy bypass -File inst/scripts/install/install_windows.ps1
    ```

## R Package

Install directly from the package source using `devtools`:

```r
install.packages("devtools")
devtools::install("path/to/HapSelect")
```

Or from a downloaded `.zip` archive:

```r
install.packages("/path/to/HapSelect", repos = NULL)
```

## R Dependencies

HapSelect imports the following packages, which will be installed automatically:

| Package | Purpose |
|---------|---------|
| `rrBLUP` | Marker effect estimation |
| `GA` | Genetic algorithm for parent selection |
| `ggplot2` | Visualizations |
| `cowplot` | Plot composition |
| `dplyr` / `purrr` | Data manipulation |
| `furrr` / `future` | Parallel computation |
| `Rcpp` | C++ extensions for LD calculation |
| `progressr` | Progress reporting |
