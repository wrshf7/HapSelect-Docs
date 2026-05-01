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
    This is slow for large marker panels. Use `plink_pairwise_ld()` for production datasets.

### `plink_pairwise_ld()`

Wraps PLINK v1.9 for fast LD computation. Requires a PLINK binary fileset (`.bed`/`.bim`/`.fam`).

```r
ld_pairs <- plink_pairwise_ld("path/to/plink_prefix")
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
