# Pairwise LD

Linkage disequilibrium (LD) between all marker pairs is the foundation of haploblock definition.

## Functions

### `order_map()`

Sorts a marker map by chromosome and position. Required before LD computation.

```r
map <- order_map(map = map)
```

### `pairwise_ld()`

Computes pairwise r² values using HapSelect's internal C++ implementation.

```r
ld_pairs <- pairwise_ld(geno, parallelize = FALSE)
```

!!! warning
    This is slow for moderate to large marker panels. Use `plink_pairwise_ld()` for production datasets.

### `plink_pairwise_ld()`

Wraps PLINK v1.9 for fast LD computation. This can be computed from `geno` or from a PLINK binary fileset (`.bed`/`.bim`/`.fam`).

```r
#from scratch with a genotype file:
ld_pairs <- plink_pairwise_ld_geno(geno = geno, ld_window = 999999, 
                                ld_window_kb = 1e6, ld_window_r2 = 0)

#with the PLINK fileset already made:
ld_pairs <- plink_pairwise_ld("path/to/plink_prefix")
```

#### Parameters

| Parameter | Description |
|:---|:---|
| `prefix` | PLINK binary prefix |
| `ld_window` | Maximum number of SNP pairs ahead to compute LD |
| `ld_window_kb` | LD window size in kb |
| `ld_window_r2` | Minimum LD threshold before termination of LD calculation |
| `extra_args` | Additional PLINK arguments |

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
- `Locus1`: numerical integer for the first marker in the marker pair. This should correspond to the order of the marker in the **ordered map file** above.
- `Locus2`: similar to `Locus1`, this is the numerical integer of the second marker in the marker pair.
- `Name1`: Character name of the first marker in the pair as seen in the genotype, map, and marker effects file.
- `Name2`: Similar to `Name1`, this corresponds to the name of the second marker in the pair as seen in the genotype, map, and marker effects file.
- `LD`: The numerical LD value computed. This is typically an $\text{r}^{2}$ value.

```r
head(HapSelect::ld_pairs)
```

## Benchmarking

HapSelect includes a benchmarking script to compare both implementations on a synthetic dataset:

```bash
Rscript inst/scripts/benchmarks/benchmark_ld.R \
  --n_markers=500 \
  --n_individuals=200 \
  --n_chr=5 \
  --missing_rate=0.02 \
  --seed=1
```

The output reports wall-clock runtimes for both methods on identical data.
