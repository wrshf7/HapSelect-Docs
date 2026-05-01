# Haplotype Blocks

Haplotype blocks (haploblocks) are contiguous genomic regions with high internal LD, defined here using a flanking-marker method.

## `def_blocks()`

Partitions the genome into haploblocks from a pairwise LD matrix.

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
```

| Parameter | Description |
|-----------|-------------|
| `ld` | Pairwise LD output from `pairwise_ld()` or `plink_pairwise_ld()` |
| `map` | Ordered marker map |
| `method` | Block definition method (`"flanking"`) |
| `threshold` | Minimum r² to consider markers in the same block |
| `tolerance` | Number of consecutive markers allowed to fall below threshold before closing a block |
| `tol_reset` | Reset tolerance counter when a strong LD pair is encountered |
| `start` | Seed strategy for block boundaries (`"LD"`) |
| `parallel` | Use parallel processing via `furrr` |

## `block_obj_to_df()`

Converts the block list object to a tidy data frame for downstream use.

```r
haploblocks <- block_obj_to_df(haploblocks, map)
```

The resulting data frame has one row per haploblock with columns for block ID, chromosome, start/end position, and constituent markers.
