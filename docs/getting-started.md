# Getting Started

This page walks through the complete HapSelect workflow using the built-in example data.

## Load the Package

```r
library(HapSelect)
```

## Required Data

HapSelect expects two primary inputs (map file and genotype file). The rest can be computed in the R package or externally as additional inputs:

- **Map file** — marker positions with columns for marker ID, chromosome, and position.
    - First column should be named `SNP` and should be a character vector, not numeric or factorized.
    - Second column should be named `Chromosome` and should be numeric for proper sorting.
    - Third column should be named `Position` and should be numeric. It can be a physical position or genetic map position.

- **Genotype matrix** — individuals × markers, coded as allele dosage (0/1/2 for diploid).
    - First three columns should be identical to the map file.
    - Columns `4:ncol(genotype_file)` should be individuals and their genotype for each marker (rows are markers, columns are individuals).
    - Names of columns 4 onwards should be the individuals' identifiers.
    - Genotypes should be dosage format: i.e., number of copies of the alternative allele (integer counts only for consistency).
    - `NA` values are allowed and are handled differently at each sage via options in the relevant functions.

- **LD file** - pairwise LD between each marker within a chromosome.
    - Can be computed internally using either the in-built function or the PLINK 1.9 wrapper function.
    - We HIGHLY recommend using the PLINK 1.9 wrapper function if PLINK 1.9 is installed. because R is not built for large, iterative computions require to compute LD pairs.
    - Pairs not present (i.e., missing) in the data frame object are allowed and are handled in the haploblocking function.
    - Columns:
        - `Chrom`: the chromosome each SNP pair belongs to.
        - `Locus1`: numerical integer for the first marker in the marker pair. This should correspond to the order of the marker in the **ordered map file** (see below for more details).
        - `Locus2`: similar to `Locus1`, this is the numerical integer of the second marker in the marker pair.
        - `Name1`: Character name of the first marker in the pair as seen in the genotype, map, and marker effeccts file.
        - `Name2`: Similar to `Name1`, this corresponds to the name of the second marker in the pair as seen in the genotype, map, and marker effects file.
        - `LD`: The numerical LD value computed. This is typically an $R^{2}$ value.

- **Marker effects** — estimated SNP effects from a genomic prediction model
    - Can be computed for basic cases in the package with BLUE/dergressed BLUP/singular adjusted phenotype or provided externally.
    - First column: `SNP` corresponding to the `SNP` column in the map and genotype file. It should be formatted as a character vector.
    - Second column: `Effect` the allele subsititution (marker) effects corresponding to the SNP in column one.

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
