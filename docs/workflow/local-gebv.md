# LocalGEBV

Local genomic estimated breeding values (localGEBV) quantify the breeding value contribution of each haploblock for each individual.

## `compute_local_GEBV()`

```r
haploblock_obj <- compute_local_GEBV(
  geno           = geno,
  marker_effects = marker_effects,
  haploblocks_df = haploblocks,
  set_missing_NA = TRUE,
  center         = TRUE
)
```

| Parameter | Description |
|-----------|-------------|
| `geno` | Genotype matrix (individuals × markers) |
| `marker_effects` | Data frame with columns for SNP ID and estimated effect |
| `haploblocks_df` | Haploblock data frame from `block_obj_to_df()` |
| `set_missing_NA` | If `TRUE`, replaces localGEBV with any missing genotypes with `NA` |
| `mean_adjust` | Centers the genotypes, required in most cases |

### `set_missing_NA`

-`TRUE` if a haplotype/marker configuration for an individual in a haploblock is missing $>=$ 1 genotype call, the localGEBV value for that individual's block is set to `NA`.

-`FALSE` localGEBV are calculated with all non-missing genotypes. The missing genotypes are imputed to the mean and thus do not affect the localGEBV value.

### `mean_adjust`

-`TRUE` centers the genotypes according to [Van Raden, 2008](https://doi.org/10.3168/jds.2007-0980). This is a necessary step in the vast majority of cases to prevent biasing localGEBV and GEBV. It does not affect haploblock variance, but it may affect how parents are chosen in the GA.

-`FALSE` does not center the genotypes. 

!!! warning
    If marker effects are computed using an additive genetic model like GBLUP backsolve, rrBLUP, BayesC, etc., this **WILL** lead to biased localGEBV and GEBV!!! Only set this to `FALSE` if you are confident in what you are doing!

## Output Object

`compute_local_GEBV()` returns a list with:

- **`$Haploblocks`** — haploblock dataframe from [Haplotype Blocks](workflow/haploblocks.md) extended with `Block_Var` (population variance of localGEBV across individuals)
- **`$Haplotype_ID_Matrix`** - blocks × individuals matrix of unique haplotype IDs listed in `$Haplotypes` dataframe
- **`$Haplotype_Effect_Matrix`** — blocks × individuals matrix of localGEBV values corresponding to `$Haplotype_ID_Matrix` IDs and `$Haplotypes` dataframe IDs and values.
- **`$Haplotypes`** - dataframe containing the unique haplotype/genotype configurations in each haploblock as well as their localGEBV estimates.

## Marker Effects Input

The `marker_effects` data frame must have:

- Column 1: SNP ID matching the map file
- Column 2: Estimated marker effect (e.g., from `rrBLUP::mixed.solve()`)

See [Marker Effects](workflow/marker-effects.md) for more information.

```r
head(HapSelect::marker_effects)
```
