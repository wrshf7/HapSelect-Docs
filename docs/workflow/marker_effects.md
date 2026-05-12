# Marker Effects

## `create_marker_effects_file()`
This function was created to create the marker effects file using the `geno` input file and a new input called `BLUE`

`BLUE` is a dataframe containing two columns:
1. `Individual`, which gives the individual's identifier found in the genotype file, `geno`. Each individual **MUST** only appear **once**!
2. `BLUE`, which contains a singular, numeric value. This value can be a Best Linear Unbiased Estimate (BLUE), de-regressed Best Linear Unbiased Prediction (BLUP) extracted from a random effect with an identity covariance structure (i.e., broad-sense heritability/repeatability effect with no additive genetic effect fit), or an adjusted phenotype.

### Marker Effect Parameters

#### `h2_method`

Controls scaling of marker variance to compute teh additive genetic variance and thus narrow-sense heritability.

Options:

- `"VanRaden"` (recommended)
- `"marker_num"` (only generally used if the genotype matrix was scaled; not applicable 99% of the time)

The `"VanRaden` method is most commonly utilized in animal and plant breeding. It follows method 1 of [Van Raden, 2008](https://doi.org/10.3168/jds.2007-0980). The marker variance is multiplied by $2*\sum{p(1-p)}$ such that each marker is scaled equivalently and is assumed to contribute equivalent variance (infinitesimal model like GBLUP or BayesC).

`"marker_num"` should only be utilzied when the genotypes have been pre-scaled to unit variance (i.e., not dosage values).

We plan to offer Yang's method in the future.

---

#### `ploidy`

Integer ploidy level.

Must be provided as an integer. Example:

```r
2L
```

NOT:

```r
2
```

---

## `n_fold_cross_validation()`
This function performs the marker effects calculation, but is utilized to evaluate model performance.

```r
CV <- n_fold_cross_validation(
  geno = geno,
  BLUE = BLUE,
  nfold = 5L,
  h2_method = "VanRaden",
  ploidy = 2L
)
```

The data is split into `nfold` groups, where each group is utilzied once for validation while training occurs on the other groups. That is, `nfold` - 1 groups are utilzied for training and 1 group is utilzied for validation. This is repeated until each group is utilized for validation.

`nfold` **MUST** be an integer (i.e., `5L` and not `5`)

---

## `cross_validation()`

```r
CV <- cross_validation(
  geno = geno,
  BLUE = BLUE,
  train_prop = 0.9,
  fold = 5L,
  h2_method = "VanRaden",
  ploidy = 2L
)
```

A more generalized version of `n_fold_cross_validation()`. This function splits the data into `train_prop` and 1 - `train_prop` proportions randomly each iteration.

`train_prop` must be a value between 0 and 1 (i.e., the proportion). This is taken as `train_prop` * nubmer of rows rounded down. The rest of the data is utilized for validation.

`fold` is an integer value (e.g., `5L` and not `5`) that is utilized to repeat the random sampling with replacement `fold` times for `fold` validations.

## Important Information
!!! warning
    You should not try to utilize the validation functions on a small population size as it is highly inadvisable to train marker effects on small population sizes. What is considered "small" is relative to the number of markers, heritability, effective population size, and population structure. For crops, less than 200 individuals would be inadvisable. For livestock, this number is likely much greater.
