# Installation

## System Dependencies

Some HapSelect functions (PLINK-backed LD calculation and BEAGLE imputation and phasing) require external tools. See below for more information.


## Quick Installer Tools (all but Beagle)

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
- The locations described are relative to the root folder in the R package directory. Please download and extract the `.zip` file and set your current directory to the root folder. This can be accomplished with the `cd` command from `bash` on a Linux/Unix or MacOS command line, for example. Example command: `cd /User/Will/HapSelect`.
- System dependencies included in the install scripts include [PLINK 1.9](https://www.cog-genomics.org/plink/) and, for Windows users, [RTools 4.5](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html). PLINK 1.9 is optional and the package will install without it, but RTools 4.5 is required for Windows users. PLINK LD function calls will not run, however, without PLINK 1.9 installed and available to the system `PATH` variable.
    - You can manually install the software using the links above if you do not wish to use the command line or Windows PowerShell to run the install scripts or lack administrator access to do so.
- Other systems (Unix/Linux and MacOS) must have a C++ compiler compatible with the `Rcpp` R package. This is normally detected and referenced by default during R installation.

### Install RTools Instructions for Windows Users
- Visit [RTools website](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html)
- Click the [Rtools45 installer hyperlink](https://cran.r-project.org/bin/windows/Rtools/rtools45/files/rtools45-6768-6492.exe)
- Download and save the installer: `rtools45-6768-6492.exe`
- Run the installation wizard to install it (requires [R](https://cran.r-project.org/bin/windows/base/) to be installed)
- It is our recommendation to use the recommended installation settings unless you have a specific reason not to (and know the consequences)

!!! warning
    The other installer link (64-bit ARM RTools45 Installer) is for a VERY specific CPU architecture! It will NOT apply to the vast majority that have an Intel or AMD CPU!

## R Dependencies - CRAN

HapSelect imports the following packages, which will be installed automatically. However, in the case of installation issues with non-zero exit status, please install these individually with `install.packages("package_name")`.


| Package | Purpose |
|---------|---------|
| `rrBLUP` | Marker effect estimation |
| `GA` | Genetic algorithm for parent selection |
| `ggplot2` | Visualisations |
| `cowplot` | Plot composition |
| `dplyr` / `purrr` | Data manipulation |
| `furrr` / `future` | Parallel computation |
| `Rcpp` | C++ extensions for LD calculation |
| `progressr` | Progress reporting |

For convenience, please run the following function to install all R packages:

```r
install.packages(c("rrBLUP", "GA", "ggplot2", "cowplot", "dplyr", "purrr", "furrr", "future", "Rcpp", "progressr"))
```

## R Dependencies - GitHub

HapSelect imports the [genomicSimulation](https://github.com/vllrs/genomicSimulation) R package, which is installed automatically if you choose to use the provided installation scripts. Otherwise it must be installed separately. Please see the linked GitHub for installation instructions.

Brief Instructions from the [genomicSimulation Manual](https://vllrs.github.io/genomicSimulation/doc/gSvignette.html):
1. Download the appropriate compiled file (NOT source!): e.g., `genomicSimulation_0.3.zip` for Windows users
2. In R, run `install.packages("path/to/downloaded/file/genomicSimulation_0.3.zip", repos = NULL)` for Windows users or `install.packages("path/to/downloaded/file/genomicSimulation_0.3.tar.gz", repos = NULL)`

## R Package - HapSelect

1. Download the [HapSelect](https://github.com/wrshf7/HapSelect/releases/) `.zip` file.
2. Extract the `.zip` file.
3. Set the working directory (in R) to the extracted file (e.g., `setwd("C:/Users/Will/Downloads/HapSelect-2.0.0/HapSelect-2.0.0/")`) and run `install.packages("HapSelect-2.0.0/", repos = NULL, type = "source")`
4. **OR** run `install.packages("path/to/HapSelect-2.0.0/HapSelect-2.0.0/", type = "source", repos = NULL)`

## PLINK and BEAGLE (Optional)
PLINK1.9 executable (.exe)/binary download: [PLINK1.9](https://www.cog-genomics.org/plink/)
BEAGLE .jar download: [BEAGLE5.5](https://faculty.washington.edu/browning/beagle/beagle.html)
To use various functions, such as: `plink_pairwise_ld_geno()`, `beagle_impute()`, and `beagle_phase()`, PLINK and BEAGLE will need to be available.

The simplest way is to place the `plink.exe` and `beagle.27Feb25.75f.jar` files into the /user/bin folder (Windows)
- For example, `C:/Users/Will/bin` folder: `C:/Users/Will/bin/plink.exe` and `C:/Users/Will/bin/beagle.27Feb25.75f.jar` for Windows
- For Linux and Unix: create the `/home/username/bin` (same as `~/bin` but replace `username` with your username) folder manually ad place `plink` and `beagle.27Feb25.75f.jar` in this folder (make sure they are executable! Run `chmod 755 plink` and `chmod 755 beagle.27Feb25.75f.jar` in the terminal to make them executable)
- For Mac: create `/Users/username/bin` and replicate the process for Linux and Unix. Ensure you replace `username` with your username

Alternative method 1: place the executable files/binaries anywhere in the system `PATH` variable, which R will pull.

Alternative method 2: modify `PATH` temporarily in R to include folders where the files are located

!!! warning
    Do not attempt to modify the system `PATH` variable unless you know what you are doing! Modifying the `PATH` within R is temporary, but the system `PATH` is something many programs rely on!

Alternative method 3: address the default location directly within HapSelect using the `set_beagle_path()` and `set_plink_path()` functions