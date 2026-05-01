# Visualizations

HapSelect provides six `ggplot2`-based plots for exploring marker, haploblock, and localGEBV data.

## Marker Effects

```r
marker_effects_plot(
  marker_effects = marker_effects$Effect,
  chr            = map$Chromosome,
  pos            = map$Position
)
```

Manhattan-style plot of estimated marker effects across the genome.

## Unique Haplotype Effects

```r
unique_haplo_effects_plot(haplo_obj = haploblock_obj)
```

Displays the distribution of unique haplotype effects within each block.

## Block Variance Funnel

```r
block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE)
```

Plots `Block_Var` (localGEBV variance) against block size — useful for identifying high-value, high-diversity blocks.

## Haploblock Structure

```r
plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks)
```

Visualizes haploblock boundaries across chromosomes.

## Marker Density

```r
plot_marker_density(map_df = map, bin_size = 500000)
```

Histogram of marker density in 500 kb bins across the genome.

## LD Decay

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

Smoothed LD decay curve (r² vs. physical distance). Supports `"gam_cr"` and loess smoothing methods.
