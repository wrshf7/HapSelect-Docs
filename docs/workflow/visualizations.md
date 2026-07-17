# Visualisations

HapSelect provides a collection of publication-quality `ggplot2` visualisation functions for exploring marker effects, haplotype architecture, haploblock structure, linkage disequilibrium, and optimization progress.

All plotting functions return standard `ggplot2` objects, allowing further customization using normal `ggplot2` syntax.

```r
p <- marker_effects_plot(...)

p +
  labs(title = "Marker Effects") +
  theme_minimal()
```

As a special case:

We use default posotion values because of the ability to choose between genetic map position (cM) and physical map position (bp), the following is an example convert the positions on some graphs to Mb:

```r
haploblock_plot <- haploblock_plot + scale_x_continuous(breaks = seq(start, end, by), labels = as.character(seq(start, end, by) / 1e6), limits = c(start, end)) + labs(x = "Position, Mb)
```

where `breaks = seq(start, end, by / 1e6)` controls where labels should go and `seq` is a function in R that takes a start position (where the x axis starts), an end position (where the x axis ends), and a by argument (where should labels be put on the x axis). The `labels = as.character(seq(start, end, by) / 1e6)` controls the names of the labels set by `breaks` and you will notice the entire sequence vector is divided by `1e6`, or 1 million. This converts the labels into Mb. Similarly, 100 Kb increments can be defined by replacing `1e6` with `1e5` and Gb can be defined by replacing `1e6` with `1e9`.


---

# Genome-wide Marker Visualisations

## Marker Effects Plot

Creates a Manhattan-style scatter plot displaying estimated marker effects across the genome.

```r
marker_effects_plot(
  marker_effects = marker_effects$Effect,
  chr            = map$Chromosome,
  pos            = map$Position
)
```

### Required Inputs

| Parameter | Description |
|:--|:--|
| `marker_effects` | Numeric vector containing estimated marker effects |
| `chr` | Numeric chromosome identifier for every marker |
| `pos` | Genomic position corresponding to every marker |

### Optional Parameters

| Parameter | Default | Description |
|:--|:--|:--|
| `colors` | `c("#A01FF0","#A7A8AA")` | Two alternating chromosome colors |

### Details

* Missing values are automatically removed.
* Marker effects, chromosome vector, and position vector must all have identical lengths.
* Chromosomes are ordered numerically.
* Returns a `ggplot` object.

!!! tip
    If you have non-numeric chromosomes you must convert them into numerics. You can then re-map the numbered chromosomes into their non-numeric values using ggplot syntax. AI tools are incredibly helpful to help with this - either provide the source code or mention it is a Manhattan-plot style ggplot and you are trying to replace the chromosome labels with another set of labels.

---

# Haplotype Effect Visualisations

These functions visualize estimated effects assigned to individual haplotypes or local genomic breeding values in a Manhattan-style plot.

## Unique Haplotype Effects

Plots the estimated effect for every unique haplotype.

```r
unique_haplo_effects_plot(
    haplo_obj = haploblock_obj
)
```

---

## Unique localGEBV Effects

Plots estimated local genomic breeding values for every unique haplotype.

```r
unique_localGEBV_effects_plot(
    haplo_obj = haploblock_obj
)
```

---

## Block Variance Manhattan Plot

Displays estimated haploblock variance along the genome in a Manhattan-style plot.

```r
block_variance_manhattan_plot(
    haplo_obj = haploblock_obj
)
```

### Required Inputs

| Parameter | Description |
|:--|:--|
| `haplo_obj` | Object returned by `compute_local_GEBV()` |

### Shared Optional Parameters

| Parameter | Default | Description |
|:--|:--|:--|
| `colors` | `c("#A01FF0","#A7A8AA")` | Alternating chromosome colours |
| `pos_type` | `"midpoint"` | Position used to represent each haploblock |

### Position Options

| Value | Description |
|:--|:--|
| `"midpoint"` | Uses the midpoint of each haploblock |
| `"start"` | Uses the first marker position |

### Returns

A `ggplot2` object.

---

# Haploblock Variance Visualisations

These plots examine the relationship between haplotype effects and estimated haploblock variance.

Internally, block variance is transformed using

```
log10(Block_Var)
```

before being min-max scaled to improve visualization.

---

## Haploblock Variance Funnel Plot

Plots haplotype effects against scaled haploblock variance.

```r
haplo_block_var_funnel_plot(
    haplo_obj = haploblock_obj
)
```

---

## localGEBV Variance Funnel Plot

Plots localGEBV values against scaled haploblock variance.

```r
local_gebv_block_var_funnel_plot(
    haplo_obj = haploblock_obj
)
```

---

## Generic Block Variance Funnel Plot

Maintained for backwards compatibility.

```r
block_var_funnel_plot(
    haplo_obj = haploblock_obj
)
```

### Required Inputs

| Parameter | Description |
|:--|:--|
| `haplo_obj` | Output from `compute_local_GEBV()` |

### Optional Parameters

| Parameter | Default | Description |
|:--|:--|:--|
| `mean_line` | `TRUE` | Draw dashed mean variance line, change to `FALSE` to remove |
| `scale_colors` | `c("blue","purple","red")` | Gradient colors for effect size on the x axis |

### Notes

* Larger values indicate higher estimated block variance.
* Color represents effect magnitude.

---

# Genome Structure Visualisations

## Haploblock Structure Plot

Displays haploblock boundaries across chromosomes.

```r
plot_haploblocks(
    haploblock_df = haploblock_obj$Haploblocks
)
```

### Required Inputs

| Parameter | Description |
|:--|:--|
| `haploblock_df` | Haploblock dataframe |

### Optional Parameters

| Parameter | Default | Description |
|:--|:--|:--|
| `block_fill` | `"#A01FF0"` | Haploblock fill colour |
| `chrom_fill` | `NA` | Chromosome background colour |
| `height` | `0.30` | Chromosome track thickness |
| `single_width_bp` | `NULL` | Width used to display single-marker blocks - mostly deprecated |

### Notes

* Transparent chromosome backgrounds are produced with `chrom_fill = NA`.
* Single-marker blocks are automatically given a visible width if `single_width_bp = NULL`.

---

## Marker Density Plot

Plots marker density along chromosomes using fixed genomic bins.

```r
plot_marker_density(
    map = map
)
```

### Required Inputs

| Parameter | Description |
|:--|:--|
| `map` | Ordered marker map containing `Chromosome` and `Position` |

### Optional Parameters

| Parameter | Default | Description |
|:--|:--|:--|
| `bin_size` | `500000` | Bin width used for counting markers |
| `height` | `0.30` | Chromosome thickness |
| `chrom_fill` | `NA` | Chromosome fill colour |
| `col_low` | `"white"` | Low-density colour |
| `col_mid` | `"purple"` | Mid-density colour |
| `col_high` | `"red"` | High-density colour |

### Notes

* Smaller bins increase resolution.
* Larger bins produce smoother density estimates.
* Final bins are automatically truncated at chromosome ends.

---

# Linkage Disequilibrium

## LD Decay Plot

Plots pairwise linkage disequilibrium (R²) as a function of physical distance.

```r
plot_ld_decay(
    map = map,
    ld = ld_pairs
)
```

### Required Inputs

| Parameter | Description |
|:--|:--|
| `map` | Marker map containing SNP names, chromosomes, and positions |
| `ld` | Pairwise LD dataframe |

### Optional Parameters

| Parameter | Default | Description |
|:--|:--|:--|
| `max_kb` | `500` | Maximum pairwise distance |
| `method` | `"gam_cr"` | Curve fitting method |
| `k` | `50` | Number of basis functions for GAM |
| `span` | `0.30` | LOESS smoothing span |
| `point_color` | `"#A7A8AA"` | Scatter colour |
| `curve_color` | `"#A01FF0"` | Fitted curve colour |
| `alpha` | `0.20` | Point transparency |

### Curve Options

| Method | Description |
|:--|:--|
| `"gam_tp"` | Thin-plate regression spline GAM |
| `"gam_cr"` | Cubic regression spline GAM |
| `"exp"` | Exponential decay model |
| `"loess"` | LOESS smoother |

### Notes

* `"exp"` produces a monotonic decay curve.
* `span` only affects LOESS.
* `k` only affects GAM fitting.

---

# Genetic Algorithm Diagnostics

## GA Progress Plot

Visualizes convergence of the genetic algorithm during parent selection.

```r
GA_progress_plot(
    parent_selection_object
)
```

### Required Inputs

| Parameter | Description |
|:--|:--|
| `parent_selection_object` | Object returned from `select_parents()` |

### Optional Parameters

| Parameter | Default | Description |
|:--|:--|:--|
| `max_color` | `"#A01FF0"` | Maximum fitness line colour |
| `mean_color` | `"#A01FF0"` | Mean fitness line colour |
| `ribbon_color` | `"#A01FF0"` | Interquartile ribbon colour |
| `ribbon_alpha` | `0.35` | Ribbon transparency |
| `max_linewidth` | `1.2` | Maximum fitness line width |
| `mean_linewidth` | `0.8` | Mean fitness line width |

### Details

The figure displays:

* Maximum objective value by generation.
* Mean objective value by generation.
* Interquartile range across the population.
* Final maximum value labelled on the plot.

Returns a `ggplot2` object.
