# Getting Started

This page walks through the complete HapSelect workflow using the built-in example data.

## Load the Package

```r
library(HapSelect)
```

## Required Data

HapSelect expects two primary inputs:

- **Map file** — marker positions with columns for marker ID, chromosome, and position
- **Genotype matrix** — individuals × markers, typically coded as allele dosage (0/1/2)
- **Marker effects** — estimated SNP effects from a genomic prediction model

Example datasets are bundled with the package:

```r
data("map", package = "HapSelect")
data("geno", package = "HapSelect")
data("marker_effects", package = "HapSelect")
```

## Step 1 — Order the Map

Markers must be sorted by chromosome and position before LD computation:

```r
map <- order_map(map = map)
```

## Step 2 — Compute Pairwise LD

```r
# Internal R implementation (slow for large datasets)
ld_pairs <- pairwise_ld(geno, parallelize = FALSE)

# Recommended: PLINK-backed implementation
ld_pairs <- plink_pairwise_ld("path/to/plink_prefix")
```

!!! tip
    For large marker panels, use `plink_pairwise_ld()` with a PLINK v1.9 binary fileset. See [Pairwise LD](workflow/pairwise-ld.md) for details.

## Step 3 — Define Haplotype Blocks

```r
haploblocks <- def_blocks(
  ld        = ld_pairs,
  map       = map,
  method    = "flanking",
  threshold = 0.2,
  tolerance = 4,
  tol_reset = TRUE,
  start     = "LD",
  parallel  = FALSE
)

haploblocks <- block_obj_to_df(haploblocks, map)
```

## Step 4 — Compute Local GEBV

```r
haploblock_obj <- compute_local_GEBV(
  geno           = geno,
  marker_effects = marker_effects,
  haploblocks_df = haploblocks,
  set_missing_NA = TRUE,
  center         = TRUE
)
```

## Step 5 — Visualize

```r
marker_effects_plot(marker_effects = marker_effects$Effect,
                    chr = map$Chromosome, pos = map$Position)

unique_haplo_effects_plot(haplo_obj = haploblock_obj)

plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks)

plot_ld_decay(map = map, ld = ld_pairs, max_kb = 500)
```

## Step 6 — Select Parents

```r
GA_output <- genetic_algorithm(
  localGEBV  = localGEBV,
  n_founders = 20,
  popSize    = 10,
  maxiter    = 300,
  run        = 150,
  selfing    = FALSE
)

GA_output$One_Solution
```

See the individual workflow pages for full parameter documentation.
