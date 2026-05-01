# Process Overview

HapSelect implements a four-stage pipeline that moves from raw genotype data to an optimised set of founder parents for a genomic selection program.

![HapSelect process overview](assets/overview-diagram.png){ .overview-diagram }

---

## A — Block Formation

The first stage partitions the genome into haplotype blocks using linkage disequilibrium (LD) between markers.

Starting from a genotype matrix of **N individuals × M markers**, HapSelect:

1. Computes pairwise r² between all marker pairs within each chromosome
2. Uses the r² matrix to define contiguous regions of high internal LD — **haploblocks**

The result is a set of **J haploblocks** that tile the genome, each capturing a segment of co-inherited variation.

```r
ld_pairs   <- plink_pairwise_ld("path/to/plink_prefix")
haploblocks <- def_blocks(ld = ld_pairs, map = map, method = "flanking",
                          threshold = 0.2, tolerance = 4)
haploblocks <- block_obj_to_df(haploblocks, map)
```

---

## B — SNP Effect Estimation

Marker effects are estimated independently using a genomic prediction model — GBLUP, BayesR, or any equivalent method that produces per-SNP effect estimates (\(\hat{u}_m\)).

HapSelect does not perform this step itself; it consumes the output. Any model that returns a vector of marker effects aligned to the same map can be used.

```r
# Example using rrBLUP
fit            <- rrBLUP::mixed.solve(y = phenotypes, Z = geno)
marker_effects <- data.frame(SNP = colnames(geno), Effect = fit$u)
```

---

## C — Local GEBV

With haploblocks defined and marker effects estimated, HapSelect calculates a **local genomic estimated breeding value** for each individual in each haploblock:

$$
\text{localGEBV}_{ij} = \sum_{m \in \text{block}_j} z_{im} \cdot \hat{u}_m
$$

where \(z_{im}\) is the allele dosage of individual \(i\) at marker \(m\), and \(\hat{u}_m\) is the estimated marker effect.

This produces an **N × J matrix** of localGEBV values — one value per individual per haploblock — which captures the breeding value contribution of each genomic region separately.

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

High-variance blocks are those where individuals differ substantially in their haplotype effects — these are the genomic regions where parent choice will have the greatest impact on offspring breeding value.

The funnel plot (`block_var_funnel_plot`) visualises this across all blocks, with haplotype effect on the x-axis and scaled block variance on the y-axis.

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
