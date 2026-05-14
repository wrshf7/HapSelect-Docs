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
  pos            = map$Position,
  colors         = c("#A01FF0", "#A7A8AA")
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `marker_effects` | Numeric vector of marker effects |
| `chr` | Chromosome vector corresponding to each marker |
| `pos` | Marker position vector |

### Optional Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `colors` | `c("#A01FF0", "#A7A8AA")` | Vector of length 2 defining alternating chromosome colors |

### Notes

- Missing values (`NA`) are automatically removed
- Chromosomes should generally be numeric for proper ordering

---

## Unique Haplotype Effects Plot

Displays the distribution of unique haplotype effects within each haploblock.

```r
unique_haplo_effects_plot(
  haplo_obj = haploblock_obj,
  colors    = c("#A01FF0", "#A7A8AA"),
  pos_type  = "midpoint"
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `haplo_obj` | Object returned from `compute_local_GEBV()` |

### Optional Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `colors` | `c("#A01FF0", "#A7A8AA")` | Vector of length 2 defining alternating chromosome colors |
| `pos_type` | `"midpoint"` | Haploblock positioning method |

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
  haplo_obj    = haploblock_obj,
  mean_line    = FALSE,
  scale_colors = c("blue", "purple", "red")
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `haplo_obj` | Object returned from `compute_local_GEBV()` |

### Optional Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `mean_line` | `FALSE` | Adds dashed mean variance line |
| `scale_colors` | `c("blue", "purple", "red")` | Vector of length 3 defining low/mid/high variance colors |

### Notes

- Haploblock variance is internally log-scaled for visualization
- Useful for identifying high-variance haploblocks

---

## Haploblock Structure Plot

Visualizes haploblock boundaries across chromosomes.

```r
plot_haploblocks(
  haploblock_df = haploblock_obj$Haploblocks,
  block_fill    = "#A01FF0",
  chrom_fill    = NA,
  height        = 0.30,
  single_width_bp = NULL
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `haploblock_df` | Haploblock dataframe |

### Optional Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `block_fill` | `"#A01FF0"` | Haploblock fill color |
| `chrom_fill` | `NA` | Chromosome background color |
| `height` | `0.30` | Chromosome track thickness |
| `single_width_bp` | `NULL` | Width of single-marker blocks (mostly deprecated) |

### Notes

- `chrom_fill = NA` produces transparent chromosome backgrounds
- Smaller `height` values create thinner chromosome tracks

---

## Marker Density Plot

Displays marker density across chromosomes using fixed genomic bins.

```r
plot_marker_density(
  map_df     = map,
  bin_size   = 500000,
  height     = 0.30,
  chrom_fill = NA,
  col_low    = "white",
  col_mid    = "purple",
  col_high   = "red"
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `map_df` | Ordered map dataframe |

### Optional Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `bin_size` | `500000` | Genomic bin size used for marker counting |
| `height` | `0.30` | Chromosome track thickness |
| `chrom_fill` | `NA` | Chromosome background color |
| `col_low` | `"white"` | Low-density color |
| `col_mid` | `"purple"` | Mid-density color |
| `col_high` | `"red"` | High-density color |

### Notes

- Smaller `bin_size` values increase resolution but also increase noise
- Larger bins produce smoother density patterns
- Final bins are automatically capped at chromosome ends

---

## LD Decay Plot

Plots linkage disequilibrium (`r²`) against physical distance.

```r
plot_ld_decay(
  map         = map,
  ld          = ld_pairs,
  max_kb      = 500,
  span        = 0.3,
  k           = 10,
  method      = "gam_cr",
  point_color = "#A7A8AA",
  curve_color = "#A01FF0",
  alpha       = 0.2
)
```

### Required Inputs

| Parameter | Description |
|:---|:---|
| `map` | Ordered map dataframe |
| `ld` | Pairwise LD dataframe |

### Optional Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `max_kb` | `500` | Maximum pairwise marker distance in kb |
| `method` | `"gam_cr"` | Curve fitting method |
| `span` | `0.3` | LOESS smoothing parameter |
| `k` | `10` | Integer number of GAM basis functions |
| `point_color` | `"#A7A8AA"` | LD point color |
| `curve_color` | `"#A01FF0"` | Fitted curve color |
| `alpha` | `0.2` | Point transparency between 0 and 1 |

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
