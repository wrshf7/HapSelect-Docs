# Pairwise LD

Linkage disequilibrium (LD) between all marker pairs is the foundation of LD-based haploblock definition (`ld_strategy()`, see [Haplotype Blocks](haploblocks.md)). Before computing LD, the map and genotype files should be put into the same, sorted order. If your genotypes are in a VCF file, see [Data Import, Imputation, and Phasing](data-import.md) first.

!!! tip
    `window_strategy()` does not use LD and `graph_strategy()` computes the LD it needs internally, so a pairwise LD table is only required for `ld_strategy()` and `plot_ld_decay()`.

## Preparing the Map and Genotypes

### `order_map()`

Sorts a marker map by chromosome and position. Required before LD computation and haploblocking.

```r
map <- order_map(map = map, verbose = FALSE)
```

| Parameter | Default | Description |
|:---|:---|:---|
| `map` | | Map file with the marker ID, chromosome, and position in the first three columns |
| `verbose` | `FALSE` | If `TRUE`, print how chromosome labels were converted to numbers |

The first three columns are renamed `SNP`, `Chromosome`, and `Position`. Markers with a missing chromosome are sorted last.

**Non-numeric chromosomes**

HapSelect requires numeric chromosomes. If the `Chromosome` column is not numeric (e.g., `"1A"`, `"chr2"`, `"chrX"`), every distinct label is sorted into natural order and numbered from 1. Natural order means `"chr2"` comes before `"chr10"`:

```r
map <- order_map(map = map, verbose = TRUE)
#> Chromosome labels were mapped to integers:
#>   chr2  -> 1
#>   chr10 -> 2
#>   chrX  -> 3
```

!!! warning
    Labels that differ only in case or punctuation (e.g., `"Chr1"` and `"chr1"`) are numbered as **separate** chromosomes, and labels that look like missing values (e.g., `"NA"`, `"."`, `"unknown"`) are numbered as real chromosomes. Both produce a warning - fix these labels yourself before ordering the map!

    Non-character marker IDs and non-numeric positions are also coerced with a warning. Check the output is what you expect.

!!! tip
    Run `order_map()` once and use the resulting map everywhere downstream. The chromosome numbers it assigns are what `theta_core_by_chr` in `graph_strategy()` and the plotting functions will refer to.

### `order_geno()`

Puts a genotype matrix into the same order as an ordered map and keeps only the markers in both. `order_map()` only reorders the three map columns, so use `order_geno()` to bring the genotypes along with it.

```r
map  <- order_map(map = map)
geno <- order_geno(geno = geno, map = map)
```

| Parameter | Default | Description |
|:---|:---|:---|
| `geno` | | Genotype matrix (dosage or phased), with the marker ID, chromosome, and position in the first three columns |
| `map` | `NULL` | Ordered map from `order_map()`. If `NULL`, `geno` is ordered by its own first three columns instead |
| `verbose` | `FALSE` | Passed to the chromosome numbering when `map = NULL` |

Markers are matched by **name**, not by row, so `geno` can arrive in any order. The chromosome and position columns of the result are taken from `map`.

- Markers in `geno` but not `map`, or in `map` but not `geno`, are dropped with a warning listing (up to five of) them.
- Duplicated marker IDs in either file are an error, because a marker could not be matched to one row.
- `map` must already have numeric chromosomes - run `order_map()` on it first.

!!! tip
    Because only shared markers are kept, a smaller map can be used to extract a subset of markers from a large genotype matrix, e.g., `order_geno(geno, map[map$Chromosome == 1, ])` for chromosome 1 only.

---

## Functions

### `pairwise_ld()`

Computes pairwise r² values between markers within each chromosome using HapSelect's internal C++ implementation. Markers on different chromosomes are never compared.

```r
ld_pairs <- pairwise_ld(
  genotype_matrix = geno,
  parallelize     = FALSE,
  window          = NULL,
  min_r2          = NULL,
  min_obs         = 2L
)
```

| Parameter | Default | Description |
|:---|:---|:---|
| `genotype_matrix` | | Dosage genotype matrix. See `HapSelect::geno` for formatting examples or [Getting Started](../getting-started.md) |
| `parallelize` | `FALSE` | Process chromosomes in parallel, using up to the number of cores allowed by `set_cpu_cores()` (see [Configuration](../configuration.md)) |
| `window` | `NULL` | Only compare markers at most this many positions apart. `NULL` compares every pair on a chromosome |
| `min_r2` | `NULL` | Drop pairs with r² below this value. `NULL` keeps every pair |
| `min_obs` | `2L` | Number of individuals that must be observed (non-missing) at both markers for an r² to be computed |

**`window`**

Without a window, the number of pairs grows with the square of the number of markers on a chromosome, so dense panels produce very large tables. A window keeps the number of pairs proportional to the number of markers. Choose a window larger than the longest block you expect to find with `ld_strategy()`, as pairs outside the window are missing from the table and are handled as missing by the haploblocking function.

**`min_obs`**

Two individuals always lie on a straight line, so a pair observed in only two individuals gives r² = 1 whatever the genotypes were. Raising `min_obs` to `3L` (or higher) refuses those pairs, at the cost of losing pairs between sparsely genotyped markers.

!!! tip
    A pair whose r² is undefined (one marker is monomorphic, or fewer than `min_obs` individuals are observed at both markers) is left out of the table rather than returned as `NA`.

**`parallelize`**

Work is split by chromosome. If `pairwise_ld()` is run serially without a window on a dataset with more than ~1,000 markers on a chromosome, it will print a message suggesting `parallelize = TRUE`. Set `parallelize = FALSE` explicitly to silence it.

### `plink_pairwise_ld_geno()` and `plink_pairwise_ld()`

Wraps PLINK v1.9 for LD computation. This can be computed from `geno` or from a PLINK binary fileset (`.bed`/`.bim`/`.fam`).

```r
#from scratch with a genotype file:
ld_pairs <- plink_pairwise_ld_geno(geno = geno, ld_window = 999999,
                                   ld_window_kb = 1e6, ld_window_r2 = 0)

#with the PLINK fileset already made:
ld_pairs <- plink_pairwise_ld("path/to/plink_prefix")
```

!!! warning
    PLINK 1.9 MUST be installed in order to use these functions. If it cannot be found you will get an error message. See [Installation](../installation.md) for more details, or use `set_plink_path()` to point HapSelect at a specific executable (see [Configuration](../configuration.md)).

#### Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `geno` | | Genotype matrix used throughout HapSelect. See `HapSelect::geno` for formatting examples or [Getting Started](../getting-started.md) |
| `prefix` | | PLINK binary prefix for .bed .bim and .fam files (e.g., prefix.bed) |
| `ld_window` | `999999` | Maximum number of SNP pairs ahead to compute LD |
| `ld_window_kb` | `1000000` | Maximum window size in kb to compute LD |
| `ld_window_r2` | `0` | Minimum r² for a pair to be reported |
| `extra_args` | `character()` | Additional PLINK arguments, e.g., `c("--maf", "0.05")` |

!!! warning
    PLINK reads hard genotype calls, so dosages in `geno` must be whole numbers (0/1/2). Fractional dosages (e.g., from `read_vcf_geno(prefer_ds = TRUE)`) cannot be written to PLINK format - use `pairwise_ld()` for those instead.

!!! tip
    PLINK is limited to the number of cores allowed by `set_cpu_cores()`. To override this for a single call, pass `--threads` in `extra_args`, e.g., `extra_args = c("--threads", "4")`.

### `format_plink_ld()`

If you would rather run PLINK yourself, `format_plink_ld()` converts a PLINK `--r2` output file into the HapSelect LD format.

```r
ld_pairs <- format_plink_ld(
  ld_path = "path/to/plink_output.ld",     #output of plink --r2
  bim     = "path/to/plink_prefix.bim"     #.bim file of the fileset LD was computed from
)
```

`Locus1` and `Locus2` are taken from the order of the markers in the `.bim` file, restarting at 1 on each chromosome, so the `.bim` file must be in the same order as the ordered map.

!!! warning
    PLINK's own defaults (`--ld-window 10 --ld-window-kb 1000 --ld-window-r2 0.2`) are much stricter than HapSelect's, and only report pairs within 10 markers with r² ≥ 0.2. On chromosome 1 of the example data this returns 490 pairs instead of 45,451! To match `plink_pairwise_ld()`, run PLINK with `--r2 --ld-window 999999 --ld-window-kb 1000000 --ld-window-r2 0`.

---

## Loading Example LD Data

```r
ld_pairs <- HapSelect::ld_pairs
```

The LD dataframe must contain the following columns:

```r
c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD")
```

### Other Information
Pairs not present (i.e., missing) in the data frame object are allowed and are handled in the haploblocking function.

Columns:

- `Chrom`: the chromosome each SNP pair belongs to (numeric).
- `Locus1`: numerical integer for the first marker in the marker pair. This should correspond to the order of the marker in the **ordered map file** above, counting from 1 on each chromosome.
- `Locus2`: similar to `Locus1`, this is the numerical integer of the second marker in the marker pair.
- `Name1`: Character name of the first marker in the pair as seen in the genotype, map, and marker effects file. Names must match exactly
- `Name2`: Similar to `Name1`, this corresponds to the name of the second marker in the pair as seen in the genotype, map, and marker effects file.
- `LD`: The numerical LD value computed. This is typically an $\text{r}^{2}$ value, and must be between 0 and 1.

```r
head(HapSelect::ld_pairs)
```
