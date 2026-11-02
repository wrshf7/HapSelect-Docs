# Configuration

HapSelect has a small number of session-wide settings: where to find the optional external tools (PLINK and BEAGLE) and how many CPU cores it may use. Every setting lasts for the current R session only. See [Making settings permanent](#making-settings-permanent) to apply them every time R starts.

| Function | Description |
|:---|:---|
| `set_plink_path(path)` | Use a specific PLINK 1.9 executable |
| `set_beagle_path(path)` | Use a specific BEAGLE 5.5 `.jar` file |
| `set_cpu_cores(cores)` | Limit the number of CPU cores HapSelect may use |
| `cpu_cores()` | Report the number of CPU cores HapSelect will use |

---

## External Tool Locations

### How PLINK and BEAGLE are found

PLINK (for `plink_pairwise_ld()` and `plink_pairwise_ld_geno()`) and BEAGLE (for the `beagle_*()` functions) are optional. When one is needed, HapSelect searches the following locations in order and uses the first match:

1. The location given to `set_plink_path()` / `set_beagle_path()`
2. The `bin` folder in your home directory used by the install scripts (`~/bin` on Linux/macOS, `C:/Users/<username>/bin` on Windows)
3. Every folder on the system `PATH`

Each candidate is checked by actually running it. PLINK must report version 1.9 and BEAGLE must report version 5.5 - anything else is skipped (or, for a location you set yourself, is an error). BEAGLE also needs a Java runtime available on the `PATH`. See [Installation](installation.md) for how to install both tools.

### `set_plink_path()`

```r
set_plink_path(path = "C:/Users/<username>/tools/plink.exe")  #the executable
set_plink_path(path = "C:/Users/<username>/tools")            #or the folder containing it

set_plink_path(path = NULL)                                   #clear the setting and search automatically again
```

### `set_beagle_path()`

```r
set_beagle_path(path = "C:/Users/<username>/tools/beagle.27Feb25.75f.jar")  #the .jar file
set_beagle_path(path = "C:/Users/<username>/tools")                         #or the folder containing it

set_beagle_path(path = NULL)                                                #clear the setting and search automatically again
```

| Parameter | Default | Description |
|:---|:---|:---|
| `path` | `NULL` | Path to the executable / `.jar` file, or the folder containing it. `NULL` clears the setting |

Both functions check the tool straight away and stop with an error if it cannot be found or is the wrong version, in which case the previous setting is kept. On success, they invisibly return the full path that will be used. The path is checked again every time the tool is run, in case it has been moved or upgraded since.

!!! tip
    On Windows, the PLINK check also catches a common mix-up: the PuTTY suite ships its own program called `plink.exe`, which is not PLINK. If `set_plink_path()` reports *"This may not be PLINK, do you have PuTTY installed?"*, point it at the PLINK 1.9 executable directly.

---

## CPU Cores

### `set_cpu_cores()`

By default, HapSelect may use every core available to the R session. `set_cpu_cores()` sets a ceiling on that for the rest of the session, which is useful on shared machines or when running other jobs alongside HapSelect.

```r
set_cpu_cores(cores = 4)     #use at most 4 cores
set_cpu_cores(cores = NULL)  #clear the limit and use every available core
```

| Parameter | Default | Description |
|:---|:---|:---|
| `cores` | `NULL` | Maximum number of cores, a number ≥ 1 (fractions are rounded down). `NULL` clears the limit |

The limit is a ceiling, not a target: HapSelect never uses more cores than are actually available, so setting a limit higher than the machine has does nothing.

!!! tip
    "Available" cores come from `future::availableCores()`, which already respects limits set by job schedulers such as SLURM or PBS (e.g., `SLURM_CPUS_PER_TASK`) and R's `mc.cores` option. On a cluster you will usually not need `set_cpu_cores()` at all.

### `cpu_cores()`

Returns the number of cores HapSelect will use, after applying the limit.

```r
cpu_cores()
#> [1] 4
```

### What the limit applies to

| Function | Parallel option | Notes |
|:---|:---|:---|
| `pairwise_ld()` | `parallelize = TRUE` | One worker per core, split by chromosome |
| `ld_strategy()` / `graph_strategy()` | `parallel = TRUE` | Split by chromosome, so never more workers than chromosomes |
| `compute_local_GEBV()` / `compute_haplotype_effects()` | `parallel = TRUE` | Blocks are split into chunks across workers |
| `plink_pairwise_ld()` / `plink_pairwise_ld_geno()` | always | Passed to PLINK as `--threads`, unless `--threads` is given in `extra_args` |
| `beagle_impute()` / `beagle_impute_geno()` / `beagle_phase_geno()` | always | Passed to BEAGLE as `nthreads=`, unless `nthreads=` is given in `extra_args` |

!!! warning
    Each parallel worker is a separate R process with its own copy of the data it is given. On large datasets, memory rather than the number of cores is often the limiting factor - if R runs out of memory in parallel mode, lower `set_cpu_cores()` (or, for localGEBV, `chunk_size`).

---

## Making settings permanent

Every setting above is stored as an R option, so you can apply them every time R starts by adding the options to your `.Rprofile` file (e.g., with `usethis::edit_r_profile()`):

```r
options(
  HapSelect.plink_path  = "C:/Users/<username>/tools/plink.exe",
  HapSelect.beagle_path = "C:/Users/<username>/tools/beagle.27Feb25.75f.jar",
  HapSelect.cpu_cores   = 4
)
```

Options set this way are checked the first time they are used, rather than when they are set.

---

## Startup Banner

HapSelect prints a banner with the installed version and a link to the matching version of this documentation when it is attached with `library(HapSelect)`. To load it silently (e.g., in scripts):

```r
suppressPackageStartupMessages(library(HapSelect))
```

!!! tip
    Features and functions differ between versions. Always use the documentation for your installed version - the version selector at the top of this site, or the link in the banner, will take you there. `packageVersion("HapSelect")` reports the installed version.
