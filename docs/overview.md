# Process Overview

HapSelect implements a four-stage pipeline that moves from raw genotype data to an optimised set of founder parents for a genomic selection program.

![HapSelect process overview](assets/overview-diagram.png){ .overview-diagram }

---

## A — Block Formation

The first stage partitions the genome into haplotype blocks using linkage disequilibrium (LD) between markers. We highly recommend using [PLINK 1.9](https://www.cog-genomics.org/plink/) to compute marker effects. If PLINK 1.9 is installed and available in the `PATH` variable, we offer a wrapper function to compute LD and format the output appropriately with the `plink_pairwise_ld_geno()` function. See [Pairwise LD](workflow/pairwise-ld.md) for more information.

Starting from a pared-down VCF structured genotype matrix of **N individuals × M markers** (see [Pairwise LD](workflow/pairwise-ld.md) for more details), HapSelect:

1. Computes pairwise r² between all marker pairs within each chromosome
2. Uses the r² matrix to define contiguous regions of high internal LD — **haploblocks**

The result is a set of **J haploblocks** that tile the genome, each capturing a segment of co-inherited variation. 
Haploblock IDs follow the convention of **Chromosome:Block**, identifying blocks nested within chromosome.

A custom haploblock dataframe can be provided, but it must follow the output structure produced from `block_obj_to_df()`.

```r
#compute pairwise LD
ld_pairs   <- plink_pairwise_ld_geno(geno = geno, ld_window = 999999,
              ld_window_kb = 1e6, ld_window_r2 = 0)

#make a haploblock list object
haploblocks <- def_blocks(ld = ld_pairs, map = map, method = "flanking",
                          threshold = 0.2, tolerance = 4)

#turn the list object into a singular data frame
haploblocks <- block_obj_to_df(haploblocks, map)
```
!!! tip
    There are many parameters that affect the behaior of `def_blocks()` with various defaults. Be sure to check out all options and how they can influence the result. See [Haplotype Blocks](workflow/haploblocks.md) for more details!
---

## B — SNP Effect Estimation

Marker effects are estimated independently using a genomic prediction model — GBLUP with backsolve, BayesR, rrBLUP, other Bayesian methods, or any equivalent method that produces per-SNP effect estimates (\(\hat{u}_m\)).

HapSelect offers a basic genomic prediction model by integrating the [rrBLUP R Package](https://cran.r-project.org/web/packages/rrBLUP/index.html). However, the implementation is limited to utilizing singular BLUE, adjusted phenotype, or deregressed BLUP for each individual and no other effects may be included in the model. For any users needing to employ more advanced modeling, please consult the packages below. 

Any model that returns a vector of marker effect estimates aligned to the same map can be used and it is up to the user to ensure marker effect estimates are correctly generated. Please see [localGEBV](workflow/local-gebv.md) for file structure details.

### Example R Packages for Genomic Prediction:
#### SNP Based Models:
[rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/index.html)
[BGLR](http://cran.r-project.org/web/packages/BGLR/index.html)
[Sommer](http://cran.r-project.org/web/packages/BGLR/index.html)
[ASREML-R](https://asreml.kb.vsni.co.uk/)

#### GBLUP Based Models:
[rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/index.html)
[Sommer](http://cran.r-project.org/web/packages/BGLR/index.html)
[ASREML-R](https://asreml.kb.vsni.co.uk/)

```r
# Example using rrBLUP integration into package
marker_effects            <- create_marker_effects_file(geno = geno, BLUE = BLUE, h2_method = "VanRaden", ploidy = 2L)
```

---

## C — Local GEBV

With haploblocks defined and marker effects estimated, HapSelect calculates a **local genomic estimated breeding value (localGEBV)** for each individual in each haploblock:

$$
\text{localGEBV}_{ij} = \sum_{m \in \text{block}_j} z_{im} \cdot \hat{u}_m
$$

where \(z_{im}\) is the **centered** (see [Van Raden, 2008](https://doi.org/10.3168/jds.2007-0980) allele dosage of individual \(i\) at marker \(m\) in haploblock \(j\), and \(\hat{u}_m\) is the estimated marker effect.

This produces an **\text{N × J}** test matrix matrix of localGEBV values — one value per individual per haploblock for \(J\) haploblocks — which captures the breeding value contribution of each genomic region separately.

```r
haploblock_obj <- compute_local_GEBV(
  geno           = geno,
  marker_effects = marker_effects,
  haploblocks_df = haploblocks,
  center         = TRUE
)
```

---

## D — Haploblock Variance

Not all haploblocks are equally informative. HapSelect ranks blocks by the **variance of localGEBV across individuals**:

$$
\text{var}(\text{localGEBV}_j) = \frac{\sum_i \left( \text{localGEBV}_{ij} - \overline{\text{localGEBV}_j} \right)^2}{N}
$$

where \text{var}(\text{localGEBV}_j) is the haploblock variance.

High-variance blocks are those where individuals differ substantially in their haplotype (localGEBV) effects — these are the genomic regions where parent choice will have the greatest impact on offspring breeding value.

The funnel plot (`block_var_funnel_plot`) visualises this across all blocks, with localGEBV effects on the x-axis and block variance on the y-axis scaled using a 0 to 1 min-max scaling procedure. Similarly, localGEBV effects can be visualized in a Manhattan-style plot using the `unique_haplo_effects_plot()` function (for more details and information, see [Visualizations](workflow/visualizations.md)).

---

## Parent Selection

The top-ranked haploblocks (by variance) are used as targets for the genetic algorithm. `genetic_algorithm()` searches for a set of **n_founders** parents that collectively carry the highest-value haplotype alleles across the selected blocks, subject to constraints on crossing scheme and population size.

```r
GA_output <- genetic_algorithm(
  localGEBV  = localGEBV,
  n_founders = 20,
  maxiter    = 300
)
GA_output$One_Solution
```

See [Parent Selection](workflow/parent-selection.md) for full parameter details.
