# Local GEBV

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
| `set_missing_NA` | Replace missing genotypes with `NA` rather than imputing |
| `center` | Center localGEBV values within each block |

## Output Object

`compute_local_GEBV()` returns a list with:

- **`$Haploblocks`** — haploblock data frame extended with `Block_Var` (variance of localGEBV across individuals)
- **`$Haplotype_Effect_Matrix`** — blocks × individuals matrix of localGEBV values

## Marker Effects Input

The `marker_effects` data frame must have:

- Column 1: SNP ID matching the map file
- Column 2: Estimated marker effect (e.g., from `rrBLUP::mixed.solve()`)

```r
data("marker_effects", package = "HapSelect")
head(marker_effects)
```
