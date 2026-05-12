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

### Important Notes 
- The locations described are realtive to the root folder in the R package directory. Please download and extract the `.zip` file and set your current directory to the root folder. This can be accomplished with the `cd` command from `bash` on a Linux/Unix or MacOS command line, for example. Example command: `cd /User/Will/HapSelect`.
- System dependencies included in the install scripts include [PLINK 1.9](https://www.cog-genomics.org/plink/) and, for Windows users, [RTools 4.5](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html). PLINK 1.9 is optional and the package will install without it, but RTools 4.5 is required for Windows users. PLINK LD function calls will not run, however, without PLINK 1.9 isntalled and available to the system `PATH` variable.
    - You can manually install the software using the links above if you do not wish to use the command line or Windows PowerShell to run the install scripts or lack administrator access to do so.
- Other systems (Unix/Linux and MacOS) must have a C++ compiler compatible with the `Rcpp` R package. This is normally detected and referenced by default during R installation.

## R Dependencies - CRAN

HapSelect imports the following packages, which will be installed automatically. However, in the case of installation issues with non-zero exit status, please install these individually with `install.packages("package_name")`.

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

## R Dependencies - GitHub

HapSelect imports the [genomicSimulation](https://github.com/vllrs/genomicSimulation) R package, which must be installed separately. Please see the linked GitHub for installation instructions.

## R Package - HapSelect

1. Download the `.zip` file [HapSelect](https://github.com/wrshf7/HapSelect).
2. Extract the `.zip` file.
3. Set the working directory (in R) to the extracted file (e.g., `setwd("C:/Users/Will/Downloads/HapSelect")`).
4. Use one of the methods below to install from the root of the directory:


**Install directly from the package source using `devtools`:**

```r
install.packages("devtools")
library(devtools)

#if setwd("path/to/HapSelect") was not run:
devtools::install("path/to/HapSelect")

#if the working directory is in the extracted directory:
devtools::install("./")

```

**Or with the base R command:**

```r
#if setwd("path/to/HapSelect") was not run:
install.packages("/path/to/HapSelect", repos = NULL, type = "source")

#if the working directory is in the extracted directory:
install.packages("./", repos = NULL, type = "source")
```

