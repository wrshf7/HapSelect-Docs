# Visualizations

HapSelect provides several `ggplot2`-based visualization functions for exploring:

- marker effects
- haploblock structure
- localGEBV architecture
- haploblock variance
- marker density
- LD decay

All visualization functions return `ggplot2` objects and can be modified using standard `ggplot2` syntax.

Example:

```r
p <- marker_effects_plot(...)

p + labs(x = "Position (cM)")
```

---

## Marker Effects Plot

Manhattan-style plot of estimated marker effects across the genome.

```r
marker_effects_plot(
  marker_effects = marker_effects$Effect,
  chr            = map$Chromosome,
  pos            = map$Position
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `marker_effects` | Numeric vector of marker effects |
| `chr` | Chromosome vector corresponding to each marker |
| `pos` | Marker position vector |

### Optional Parameters

| Parameter | Description |
|:---|:---|
| `colors` | Vector of length 2 defining alternating chromosome colors |

### Notes

- Missing values (`NA`) are automatically removed
- Chromosomes should generally be numeric for proper ordering

---

## Unique Haplotype Effects Plot

Displays the distribution of unique haplotype effects within each haploblock.

```r
unique_haplo_effects_plot(
  haplo_obj = haploblock_obj
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `haplo_obj` | Object returned from `compute_local_GEBV()` |

### Optional Parameters

| Parameter | Description |
|:---|:---|
| `colors` | Vector of length 2 defining alternating chromosome colors |
| `pos_type` | Haploblock positioning method |

### `pos_type` Options

| Option | Description |
|:---|:---|
| `"midpoint"` | Plot effects at the midpoint of the haploblock |
| `"start"` | Plot effects at the first marker in the haploblock |

---

## Block Variance Funnel Plot

Plots haploblock variance (`Block_Var`) against haploblock size.

```r
block_var_funnel_plot(
  haplo_obj = haploblock_obj,
  mean_line = FALSE
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `haplo_obj` | Object returned from `compute_local_GEBV()` |

### Optional Parameters

| Parameter | Description |
|:---|:---|
| `mean_line` | Adds dashed mean variance line |
| `scale_colors` | Vector of length 3 defining low/mid/high variance colors |

### Notes

- Haploblock variance is internally log-scaled for visualization
- Useful for identifying high-variance haploblocks

---

## Haploblock Structure Plot

Visualizes haploblock boundaries across chromosomes.

```r
plot_haploblocks(
  haploblock_df = haploblock_obj$Haploblocks
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `haploblock_df` | Haploblock dataframe |

### Optional Parameters

| Parameter | Description |
|:---|:---|
| `block_fill` | Haploblock fill color |
| `chrom_fill` | Chromosome background color |
| `height` | Chromosome track thickness |
| `single_width_bp` | Width of single-marker blocks (mostly deprecated) |

### Notes

- `chrom_fill = NA` produces transparent chromosome backgrounds
- Smaller `height` values create thinner chromosome tracks

---

## Marker Density Plot

Displays marker density across chromosomes using fixed genomic bins.

```r
plot_marker_density(
  map_df = map,
  bin_size = 500000
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `map_df` | Ordered map dataframe |

### Optional Parameters

| Parameter | Description |
|:---|:---|
| `bin_size` | Genomic bin size used for marker counting |
| `height` | Chromosome track thickness |
| `chrom_fill` | Chromosome background color |
| `col_low` | Low-density color |
| `col_mid` | Mid-density color |
| `col_high` | High-density color |

### Notes

- Smaller `bin_size` values increase resolution but also increase noise
- Larger bins produce smoother density patterns
- Final bins are automatically capped at chromosome ends

---

## LD Decay Plot

Plots linkage disequilibrium (`r²`) against physical distance.

```r
plot_ld_decay(
  map    = map,
  ld     = ld_pairs,
  max_kb = 500,
  span   = 0.3,
  k      = 10,
  method = "gam_cr"
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `map` | Ordered map dataframe |
| `ld` | Pairwise LD dataframe |

### Optional Parameters

| Parameter | Description |
|:---|:---|
| `max_kb` | Maximum pairwise marker distance in kb |
| `method` | Curve fitting method |
| `span` | LOESS smoothing parameter |
| `k` | Integer number of GAM basis functions |
| `point_color` | LD point color |
| `curve_color` | Fitted curve color |
| `alpha` | Point transparency between 0 and 1 |

### `method` Options

| Option | Description |
|:---|:---|
| `"gam_tp"` | Thin-plate regression spline GAM |
| `"gam_cr"` | Cubic regression spline GAM |
| `"exp"` | Exponential decay model |
| `"loess"` | LOESS smoothing |

### Notes

- `"exp"` enforces monotonic decay
- `span` is only used by `"loess"`
- `k` is only used by `"gam_tp"` and `"gam_cr"`
- Smaller `span` or larger `k` values may increase overfitting
