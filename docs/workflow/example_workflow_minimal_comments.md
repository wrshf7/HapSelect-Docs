# HapSelect Workflow with Minimal Comments

Minimal workflow demonstrating:

1. LD calculation
2. Haploblock construction
3. Marker effect estimation
4. localGEBV computation
5. Visualisation
6. GA parent selection
7. Recurrent truncation selection simulation on GA selected parents and TS selected parents

---

## Load Package and Example Data

```r
library(HapSelect)

map         <- HapSelect::map
geno        <- HapSelect::geno
geno_phased <- HapSelect::geno_phased
```

---

## Order Map File

```r
# simulate unordered map
map2 <- map[sample(1:nrow(map), nrow(map)), ]

# reorder
map2 <- order_map(map = map2)
```

---
## LD Calculation
### Pairwise LD (Native Function)

```r
# very slow; demonstration only
# ld_pairs <- pairwise_ld(geno, parallelize = FALSE)
```

---

### Pairwise LD with PLINK (Recommended)

```r
ld_pairs <- plink_pairwise_ld_geno(
  geno = geno,
  ld_window = 999999,
  ld_window_kb = 1e6,
  ld_window_r2 = 0
)
```

Requires PLINK v1.9 installed and available in the system `PATH`.


## Load Example LD Data

```r
ld_pairs <- HapSelect::ld_pairs
```

---

## Construct Haploblocks

```r
haploblocks <- def_blocks(
  ld = ld_pairs,
  map = map,
  method = "flanking",
  threshold = 0.4,
  tolerance = 3,
  tol_reset = TRUE,
  start = "beginning",
  parallel = FALSE
)
```

---

## Convert to Dataframe

```r
haploblocks <- block_obj_to_df(haploblocks, map)
```

---

### Summarise Haploblocks

```r
block_summary(block_df = haploblocks)
```

---

## Marker Effect Estimation
### Estimate marker effects if not provided
#### Important Note:
- One BLUE, deregressed BLUP, or adjusted phenotype per individual only
- rrBLUP implementation
- no other effects in the model, see [Documentation Overview](https://wrshf7.github.io/HapSelect-Docs/overview/) for more details and alternative R packages for independent modelling.

```r
BLUE <- HapSelect::BLUE

marker_effects <- create_marker_effects_file(
  geno = geno,
  BLUE = BLUE,
  h2_method = "VanRaden",
  ploidy = 2L
)
```

---

### Cross Validation

#### N-fold

```r
CV <- n_fold_cross_validation(
  geno = geno,
  BLUE = BLUE,
  nfold = 5L,
  h2_method = "VanRaden",
  ploidy = 2L
)
```

#### Random Train/Test

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

---

## Compute localGEBV

```r
marker_effects <- HapSelect::marker_effects

#localGEBV
haploblock_obj <- compute_local_GEBV(
  geno = geno,
  marker_effects = marker_effects,
  haploblocks_df = haploblocks,
  set_missing_NA = TRUE,
  mean_adjust = TRUE,
  parallel = TRUE
)

#haplotypes
haploblock_obj <- compute_haplotype_effects(
  geno = geno_phased,
  marker_effects = marker_effects,
  haploblocks_df = haploblocks,
  set_missing_NA = TRUE,
  mean_adjust = TRUE,
  parallel = TRUE
)
```

---

## Visualisations

### Editing Plots Example

```r
haploblock_plot <- haploblock_plot + labs(x = "cM")

haploblock_plot <- haploblock_plot + scale_x_continuous(breaks = seq(start, end, by), labels = as.character(seq(start, end, by) / 1e6), limits = c(start, end)) + labs(x = "Position, Mb)
```


### Marker Effects Plot

```r
marker_plot <- marker_effects_plot(
  marker_effects = marker_effects$Effect,
  chr = map$Chromosome,
  pos = map$Position,
  colors = c("#A01FF0", "#A7A8AA")
)

marker_plot
```

---

### Unique localGEBV/Haplotype Effects Plot

```r
localGEBV_plot <- unique_localGEBV_effects_plot(
  haplo_obj = haploblock_obj,
  colors = c("#A01FF0", "#A7A8AA"),
  pos_type = "midpoint"
)

localGEBV_plot

haplo_eff_plot <- unique_haplo_effects_plot(
  haplo_obj = haploblock_obj,
  colors = c("#A01FF0", "#A7A8AA"),
  pos_type = "midpoint"
)

haplo_eff_plot
```

---

### Funnel Plot

```r
#localGEBV
funnel_plot <- local_gebv_block_var_funnel_plot(
  haplo_obj = haploblock_obj,
  mean_line = FALSE,
  scale_colors = c("blue", "purple", "red")
)

#haplotype
funnel_plot <- haplo_block_var_funnel_plot(
  haplo_obj = haploblock_obj,
  mean_line = FALSE,
  scale_colors = c("blue", "purple", "red")
)

funnel_plot
```

---

### Haploblock Position Plot

```r
haploblock_plot <- plot_haploblocks(
  haploblock_df = haploblock_obj$Haploblocks, #or haploblocks (haplobock dataframe)
  block_fill = "#A01FF0",
  chrom_fill = NA,
  height = 0.30,
  single_width_bp = NULL
)

haploblock_plot
```

---

### Marker Density Plot

```r
marker_density_plot <- plot_marker_density(
  map = map,
  bin_size = 500e3, #500 kb
  height = 0.3,
  chrom_fill = NA,
  col_low = "white",
  col_mid = "purple",
  col_high = "red"
)

marker_density_plot
```

---

### LD Decay Plot

```r
ld_decay_plot <- plot_ld_decay(
  map = map,
  ld = ld_pairs,
  max_kb = 500,
  span = 0.3,
  k = 10L,
  method = "gam_cr",
  point_color = "#A7A8AA",
  curve_color = "#A01FF0",
  alpha = 0.2
)

ld_decay_plot
```

### Block Variance Manhattan Plot

```r
variance_plot = block_variance_manhattan_plot(
  haplo_obj = haploblock_obj,
  colors = c("#A01FF0", "#A7A8AA"),
  pos_type = "midpoint"
)
variance_plot

```

---

## Select Haploblocks for the GA

### Top N Blocks

```r
haploblock_obj <- select_top_blocks(
  haploblock_obj = haploblock_obj,
  n = 15
)

nrow(haploblock_obj$Haploblocks_GA)
```

### Top Percentage of Blocks

```r
haploblock_obj <- select_top_blocks(
  haploblock_obj = haploblock_obj,
  perc_total = 0.5
)

nrow(haploblock_obj$Haploblocks_GA)
```

### Variance Explained Threshold

```r
haploblock_obj <- select_top_blocks(
  haploblock_obj = haploblock_obj,
  perc_of_total_var = 0.9
)

nrow(haploblock_obj$Haploblocks_GA)
```

---

## Genetic Algorithm

```r
localGEBV_parent_obj <- local_gebv_parent_selection(
  haploblock_obj = haploblock_obj,
  n_founders = 20,
  popSize = 10,
  maxiter = 300,
  run = 150,
  strategy = "no_selfing",
  pmutation = 0.6,
  pcrossover = 0.6,
  maximize = TRUE,
  monitor = TRUE
)

haplotype_parent_obj <- haplotype_parent_selection(
  haploblock_obj = haploblock_obj,
  n_founders = 20,
  popSize = 10,
  maxiter = 300,
  run = 150,
  strategy = "OHS",
  pmutation = 0.6,
  pcrossover = 0.6,
  maximize = TRUE,
  monitor = TRUE
)
```

---

### Inspect One Solution

```r
localGEBV_parent_obj$selected_founders #substitute haplotype_parent_obj for haplotypes
```

---

### Optimisation Progress Plot

```r
optimisation_progress_plot <- GA_progress_plot(
  parent_selection_object = localGEBV_parent_obj, #substitute haplotype_parent_obj for haplotypes
  max_color = "#A01FF0",
  mean_color = "#A01FF0",
  ribbon_color = "#A01FF0",
  ribbon_alpha = 0.35,
  max_linewidth = 1.2,
  mean_linewidth = 0.8
)
```

---

## Simulate Recurrent Selection

```r
localGEBV_Sim <- localGEBV_vs_TS_simulation(
  GA_output = localGEBV_parent_obj,
  geno = geno,
  marker_effects = marker_effects,
  map = map,
  genetic_map_position = NULL,
  num_gen = 50,
  num_sim_reps = 30,
  num_cross_per_gen = 100,
  num_TS_parents = NULL,
  mean_adjust = TRUE,
  maximize = TRUE,
  max_cM_chr = 100,
  PCA = TRUE,
  colors = c("green", "#d95f02", "#A01FF0", "gray80"),
  alpha = c(1,1,1,0.5)
)

Haplotype_Sim <- Haplotype_vs_TS_simulation(
  GA_output = localGEBV_parent_obj,
  geno_phased = geno,
  marker_effects = marker_effects,
  map = map,
  genetic_map_position = NULL,
  num_gen = 50,
  num_sim_reps = 30,
  num_cross_per_gen = 100,
  num_TS_parents = NULL,
  mean_adjust = TRUE,
  maximize = TRUE,
  max_cM_chr = 100,
  PCA = TRUE,
  colors = c("green", "#d95f02", "#A01FF0", "gray80"),
  alpha = c(1,1,1,0.5)
)
```

---

### Display Simulation Results

```r
localGEBV_Sim$Simulation_Plot

localGEBV_Sim$PCA_Plot

#substitute localGEBV_Sim with Haplotype_Sim for haplotypes
```

---

