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
| `method` | Block definition method (`"flanking"` or `"average"`) |
| `threshold` | Minimum r² value to consider adding a marker to the block |
| `tolerance` | Integer number of markers allowed to fall below threshold before terminating a block |
| `tol_reset` | If `TRUE`, reset tolerance counter when a marker is successfully added to a block |
| `start` | Seed strategy for block boundaries (`"LD"` or `"Beginning"`) |
| `parallel` | Use parallel processing via `furrr` |

**`method`**
- `"flanking"` → compares the LD of the adjacent marker only (i.e., compares current first/last marker in the block to the next marker)
- `"average"` → compares average LD across block (i.e., averages the next marker's LD to all markers currently in the block)

**`tolerance`**
- If the next marker to compare does not meet or exceed `threshold`, the counter is incremented by 1. Once the counter is greater than the `tolerance` value the block is terminated. Any markers between a successfully added marker and a the block that did not meet the threshold will also be added to the block. This paramemter helps to accomodate for reference alignment error, genotyping error, structural variation, and other systematic errors

**`tol_reset`**
- `TRUE` the `tolerance` counter is reset to 0 when a marker meets the threshold.
- `FALSE` the counter is not reset and will keep incrementing even if a marker is successfully added.

**`start`**
- `LD` -> starts from the highest LD pair and the block will be extended both to the left and to the right according to the distance coordinate. The tolerance counter is reset between left extension and right extension.
- `Beginning` -> starts at the beginning of the chromosome, so blocks are only extended to the right.

## `block_obj_to_df()`

Converts the block list object to a tidy data frame for downstream use.

```r
haploblocks <- block_obj_to_df(haploblocks, map)
```

The resulting data frame has one row per haploblock with columns for block ID, chromosome, start/end position, and constituent markers. The constituent markers are `;` delimited.

!!! warning Do not use any kind of whitespace as a delimeter, this will not be considered a delimeter and will be part of the marker name!
