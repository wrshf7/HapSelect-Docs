# Data Import, Imputation, and Phasing

HapSelect can read genotypes straight from a VCF file, and impute and phase them with [BEAGLE 5.5](https://faculty.washington.edu/browning/beagle/beagle.html), so that the genotype matrix (`geno`, localGEBV) and phased genotype matrix (`geno_phased`, haplotypes) can be produced without leaving R.

| Function | Input | Output |
|:---|:---|:---|
| `read_vcf_geno()` | VCF file | Dosage genotype matrix (`geno`) |
| `read_vcf_phased()` | Phased VCF file | Phased genotype matrix (`geno_phased`) |
| `beagle_impute_geno()` | `geno` with missing values | Imputed `geno` |
| `beagle_phase_geno()` | `geno` | Imputed and phased `geno_phased` |
| `beagle_impute()` | VCF file | Imputed and phased VCF file |
| `center_genotypes()` | `geno` or `geno_phased` | Centered genotypes |

See [Getting Started](../getting-started.md) for the structure of the genotype and phased genotype matrices.

!!! warning
    The VCF functions do not currently handle aneuploidy or variable ploidy - every genotype call in a file must have the same number of alleles. BEAGLE, and therefore every `beagle_*()` function, only supports diploid genotypes.

---

## Reading VCF Files

### `read_vcf_geno()`

Reads a VCF file (`.vcf` or gzipped `.vcf.gz`) into a dosage genotype matrix: the number of alternative (ALT) alleles at each marker for each individual.

```r
geno <- read_vcf_geno(
  path      = "path/to/genotypes.vcf.gz",
  prefer_ds = FALSE,              #if TRUE, read imputed dosages from the DS field instead of GT
  ploidy    = 2L                  #number of alleles per genotype call
)
```

| Parameter | Default | Description |
|:---|:---|:---|
| `path` | | Path to a `.vcf` or `.vcf.gz` file |
| `prefer_ds` | `FALSE` | Read dosages from the `DS` field instead of `GT` |
| `ploidy` | `2L` | Number of alleles in every genotype call. Dosages run from 0 to `ploidy` |

The result has the `SNP` (from the VCF `ID` column), `Chromosome`, and `Position` columns followed by one column per sample.

- Phased (`0|1`) and unphased (`0/1`) calls are both read.
- A call with any missing allele (`./.`, `0/.`) becomes `NA`.
- A call with a second ALT allele (e.g., `1/2`) is an error, because a single ALT count cannot express it. Split multi-allelic sites into biallelic markers before reading (e.g., with `bcftools norm -m -`).
- A call with a different number of alleles to `ploidy` is an error.

**`prefer_ds`**

- `FALSE` -> dosages are the hard calls from `GT` (whole numbers).
- `TRUE` -> dosages are read from the `DS` field, which imputation software (including BEAGLE with a reference panel) writes alongside `GT`. These are fractional (e.g., `0.83` rather than `1`), keeping the imputation's uncertainty. If **any** record lacks `DS`, the whole file is read from `GT` instead.

!!! warning
    Fractional dosages from `prefer_ds = TRUE` work with `pairwise_ld()` and all haploblocking strategies, but cannot be written to PLINK format, so they will not work with `plink_pairwise_ld_geno()`!

!!! tip
    Chromosomes are read as they are written in the VCF (e.g., `"1A"`, `"chr2"`). Run `order_map()` and `order_geno()` on the result to number the chromosomes and sort the markers before using it anywhere else:

    ```r
    geno <- read_vcf_geno("path/to/genotypes.vcf.gz")
    map  <- order_map(geno[, 1:3], verbose = TRUE)
    geno <- order_geno(geno, map)
    ```

### `read_vcf_phased()`

Reads a **phased** VCF file into the phased genotype matrix used for haplotype effects. Each sample's call is split into one column per chromosome copy (`sample_1`, `sample_2`, ...) of 0 (REF) / 1 (ALT) alleles.

```r
geno_phased <- read_vcf_phased(
  path   = "path/to/phased.vcf.gz",
  ploidy = 2L
)
```

| Parameter | Default | Description |
|:---|:---|:---|
| `path` | | Path to a phased `.vcf` or `.vcf.gz` file |
| `ploidy` | `2L` | Number of alleles in every genotype call |

!!! warning
    Every call must be fully phased (`|` between every allele) with no missing alleles, because each allele has to be placed on a chromosome copy. An unphased (`0/1`) or missing (`.|.`) call is an error naming the sample - phase and impute the genotypes first (e.g., with `beagle_phase_geno()`), or use `read_vcf_geno()` to read the file as dosages.

---

## Imputation and Phasing with BEAGLE

The `beagle_*()` functions run [BEAGLE 5.5](https://faculty.washington.edu/browning/beagle/beagle.html) for you. BEAGLE always phases, and phasing fills in missing genotypes as it goes, so a single BEAGLE run both phases and imputes. The functions only differ in what they read and return.

!!! warning
    BEAGLE 5.5 and a Java runtime (Java 8 or later) MUST be installed to use these functions. HapSelect looks for a `beagle*.jar` file in the install scripts' `~/bin` folder, then in every folder on the system `PATH`, and only accepts BEAGLE 5.5. Use `set_beagle_path()` to point HapSelect at a specific `.jar` file (see [Configuration](../configuration.md)). See [Installation](../installation.md) for more details.

All three functions share the following optional parameters:

| Parameter | Default | Description |
|:---|:---|:---|
| `ref` | `NULL` | Path to a reference panel VCF for reference-based imputation |
| `map` | `NULL` | Path to a PLINK-format genetic map, passed to BEAGLE's `map=` argument. This is **not** the HapSelect map data frame! |
| `extra_args` | `character()` | Additional BEAGLE `key=value` arguments, e.g., `c("ne=100", "window=40")` |

!!! tip
    Without `ref`, BEAGLE imputes missing genotypes from the samples in your own data. A reference panel is only needed to impute markers that are absent from your data, or to improve accuracy with a larger, related population. Without `map`, BEAGLE assumes a constant recombination rate of 1 cM per Mb.

!!! tip
    BEAGLE is limited to the number of cores allowed by `set_cpu_cores()`. To override this for a single call, pass `nthreads=` in `extra_args`, e.g., `extra_args = c("nthreads=4")`.

### `beagle_phase_geno()`

Imputes and phases a dosage genotype matrix, returning the phased genotype matrix (`geno_phased`) that `compute_haplotype_effects()` expects.

```r
geno_phased <- beagle_phase_geno(
  geno       = geno,
  ref        = NULL,
  map        = NULL,
  extra_args = character()
)
```

### `beagle_impute_geno()`

Imputes missing genotypes in a dosage genotype matrix and returns a dosage genotype matrix in the same layout, sample order, and sample names as `geno`.

```r
geno_imputed <- beagle_impute_geno(
  geno       = geno,
  ref        = NULL,
  map        = NULL,
  extra_args = character(),
  prefer_ds  = FALSE
)
```

| Parameter | Default | Description |
|:---|:---|:---|
| `prefer_ds` | `FALSE` | If `TRUE`, return BEAGLE's fractional imputed dosages (`DS`) instead of hard calls. See `read_vcf_geno()` above |

!!! warning
    BEAGLE only writes `DS` when imputing against a reference panel. `prefer_ds = TRUE` without `ref` has no effect and produces a warning.

!!! tip
    `beagle_impute_geno()` discards the phasing BEAGLE produced, because a dosage cannot express which chromosome an allele sits on. If you need both a dosage matrix (localGEBV) and a phased matrix (haplotypes) from the same imputation, run `beagle_phase_geno()` and sum each individual's two chromosome columns to get dosages.

For both `beagle_phase_geno()` and `beagle_impute_geno()`, `geno` must be diploid (dosages of 0, 1, 2, or `NA`).

### `beagle_impute()`

Runs BEAGLE on a VCF file directly and writes the imputed and phased result to `<out_prefix>.vcf.gz`. Use this when your data is already in VCF format or too large to hold comfortably in R.

```r
out_vcf <- beagle_impute(
  vcf_in     = "path/to/genotypes.vcf.gz",
  out_prefix = "path/to/imputed",          #BEAGLE writes path/to/imputed.vcf.gz and path/to/imputed.log
  ref        = NULL,
  map        = NULL,
  extra_args = character()
)

#read the result back in
geno_phased <- read_vcf_phased(out_vcf)
```

`beagle_impute()` returns the path to the output VCF.

---

## Other Genotype Preparation

### `center_genotypes()`

Centers a genotype (or phased genotype) matrix by subtracting each marker's mean from every individual. Missing values are ignored when calculating the mean and remain `NA`.

```r
geno_centered <- center_genotypes(geno = geno)
```

!!! tip
    You do not need to center genotypes yourself before computing localGEBV or haplotype effects - `compute_local_GEBV()` and `compute_haplotype_effects()` do this internally with `mean_adjust = TRUE`. `center_genotypes()` is provided for custom analyses outside the HapSelect workflow.
