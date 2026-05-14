# Basic GA vs TS Simulation

The `GA_vs_TS_simulation()` function performs a simplified recurrent genomic selection simulation comparing:

- parents selected using the HapSelect genetic algorithm (GA)
- parents selected using standard genomic truncation selection (TS)

The simulation evaluates long-term breeding value trajectories across generations and can optionally visualize population structure using PCA.

---

# Overview

The simulation performs the following steps:

1. Extracts GA-selected parents from the `genetic_algorithm()` output
2. Selects TS parents based on highest GEBV
3. Converts physical marker positions to genetic positions (if needed)
4. Simulates recurrent crossing and selection across generations
5. Repeats the simulation multiple times to quantify Monte Carlo variability
6. Summarizes breeding value trajectories
7. Generates trajectory and PCA visualizations

The simulation is intentionally simple and designed primarily for:

- comparing parent selection strategies
- evaluating long-term gain potential
- exploring diversity retention
- demonstrating GA optimization behavior

---

# Running the Simulation

```r
parent_sln_obj <- GA_vs_TS_simulation(
  GA_output           = GA_output,
  geno                = geno,
  marker_effects      = marker_effects,
  map                 = map,
  genetic_map_position = NULL,
  num_gen             = 50,
  num_sim_reps        = 30,
  num_cross_per_gen   = 1000,
  num_TS_parents      = NULL,
  mean_adjust         = TRUE,
  max_cM_chr          = 100,
  PCA                 = TRUE,
  colors              = c("green", "#d95f02", "#A01FF0", "gray80"),
  alpha               = c(1,1,1,0.5)
)
```

---

# Required Inputs

| Parameter | Description |
|:---|:---|
| `GA_output` | Output object returned from `genetic_algorithm()` |
| `geno` | Genotype dataframe used throughout the HapSelect workflow |
| `marker_effects` | Marker effects dataframe |
| `map` | Ordered map dataframe |

---

# Optional Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `genetic_map_position` | `NULL` | Optional vector of genetic map positions in cM |
| `num_gen` | `50` | Number of recurrent selection generations |
| `num_sim_reps` | `30` | Number of simulation replicates |
| `num_cross_per_gen` | `1000` | Number of progeny generated per generation |
| `num_TS_parents` | `NULL` | Number of truncation-selected parents |
| `mean_adjust` | `TRUE` | Whether to internally center markers |
| `max_cM_chr` | `100` | Chromosome genetic length assumption if no genetic map is supplied |
| `PCA` | `TRUE` | Whether to compute PCA visualization |
| `colors` | `c("green", "#d95f02", "#A01FF0", "gray80")` | Plot colors |
| `alpha` | `c(1,1,1,0.5)` | PCA transparency values |

---

# Simulation Workflow

## 1. Parent Extraction

The function extracts:

- GA-selected parents from `GA_output`
- TS-selected parents based on highest GEBV

If `num_TS_parents = NULL`, the number of TS parents automatically matches the number of GA parents.

---

## 2. Marker Compatibility Checks

The function internally verifies that:

- marker IDs match between `geno`, `marker_effects`, and `map`
- marker IDs contain no duplicates
- marker IDs contain no missing values

The simulation will stop immediately if incompatibilities are detected.

---

## 3. GEBV Calculation

The function computes genomic estimated breeding values (GEBV) using:

```r
GEBV = Zu
```

where:

- `Z` is the genotype matrix
- `u` is the marker effect vector

If `mean_adjust = TRUE`, markers are internally centered before GEBV calculation.

This should remain `TRUE` in almost all analyses.

---

## 4. TS Parent Selection

Truncation-selected parents are chosen as the individuals with the highest GEBV values.

Example:

```r
num_TS_parents = 20
```

selects the top 20 individuals ranked by GEBV.

---

## 5. Genetic Map Construction

The simulation requires marker positions in centiMorgans (cM).

Two approaches are supported.

---

## Using a True Genetic Map (Recommended)

```r
genetic_map_position = map$cM
```

The vector must:

- be the same length as `map`
- contain no missing values
- be in the same marker order as `map`

This is the most biologically realistic option.

---

## Inferring Genetic Positions from Physical Distance

If:

```r
genetic_map_position = NULL
```

then marker positions are inferred proportionally from physical distance.

Each chromosome is assumed to span:

```r
max_cM_chr = 100
```

by default.

Example:

```r
max_cM_chr = 150
```

would assume each chromosome spans 150 cM.

---

# Recurrent Selection Simulation

The simulation uses the `genomicSimulation` package internally.

For each generation:

1. Mean breeding value is recorded
2. Top individuals are selected by GEBV
3. Random crosses generate the next generation
4. The process repeats for `num_gen` generations

This procedure is independently repeated:

```r
num_sim_reps
```

times to quantify stochastic simulation variability.

---

# Understanding Simulation Parameters

## `num_gen`

Controls the number of recurrent selection generations.

### Larger Values

- better evaluate long-term gain
- reveal selection plateaus
- better assess diversity preservation

### Smaller Values

- faster simulations
- emphasize short-term gain

---

## `num_sim_reps`

Controls the number of independent simulation replicates.

Replicates differ because of:

- recombination randomness
- stochastic inheritance
- random mating patterns

### Larger Values

Advantages:

- smoother trajectories
- reduced stochastic noise
- more stable estimates

Disadvantages:

- slower runtime
- increased memory usage

---

## `num_cross_per_gen`

Controls the number of progeny generated each generation.

### Larger Values

Advantages:

- stronger selection intensity
- greater opportunity for favorable recombination
- smoother trajectories

Disadvantages:

- slower simulations
- increased memory usage

---

## `mean_adjust`

Controls genotype centering before GEBV calculation.

Recommended setting:

```r
mean_adjust = TRUE
```

This should only be disabled if marker effects were estimated from an uncentered genotype matrix.

---

## `PCA`

If `TRUE`, PCA is performed on the genotype matrix.

The PCA visualization highlights:

- GA-selected parents
- TS-selected parents
- overlapping parents
- non-selected individuals

This helps visualize:

- diversity retention
- population structure
- overlap between strategies

---

# PCA Plot Colors

The `colors` argument must contain exactly four valid R colors in the following order:

```r
colors = c(
  "green",
  "#d95f02",
  "#A01FF0",
  "gray80"
)
```

| Position | Meaning |
|:---|:---|
| 1 | GA-selected parents |
| 2 | TS-selected parents |
| 3 | Individuals selected by both methods |
| 4 | Non-selected individuals |

!!! tip
    Even if `PCA = FALSE`, this should still be specified with 4 colors as the first two are utilized simultaneously for the `Simulation_Plot`

---

# PCA Transparency

The `alpha` argument controls point transparency on the PCA plot.

```r
alpha = c(1,1,1,0.5)
```

| Position | Meaning |
|:---|:---|
| 1 | GA-selected parents |
| 2 | TS-selected parents |
| 3 | Overlapping parents |
| 4 | Non-selected individuals |

All values must be between:

```r
0 and 1
```

---

# Simulation Outputs

If:

```r
PCA = TRUE
```

the function returns:

```r
list(
  Simulation_Plot,
  PCA_Plot,
  PCA_df
)
```

---

## `Simulation_Plot`

A `ggplot2` trajectory plot showing:

- mean breeding value across generations
- GA-selected parent trajectories
- TS-selected parent trajectories
- standard error ribbons across replicates

Display:

```r
parent_sln_obj$Simulation_Plot
```

---

## `PCA_Plot`

A PCA visualization of the genotype matrix showing selected parents.

Display:

```r
parent_sln_obj$PCA_Plot
```

---

## `PCA_df`

A dataframe containing:

- accession names
- PC coordinates
- parent selection groupings

Access:

```r
head(parent_sln_obj$PCA_df)
```

---

# Example Workflow

```r
# Run GA
GA_output <- genetic_algorithm(
  localGEBV = haploblock_obj$Haplotype_Effect_Matrix_GA,
  n_founders = 20
)

# Run simulation
parent_sln_obj <- GA_vs_TS_simulation(
  GA_output = GA_output,
  geno = geno,
  marker_effects = marker_effects,
  map = map,
  num_gen = 50,
  num_sim_reps = 30
)

# Display plots
parent_sln_obj$Simulation_Plot

parent_sln_obj$PCA_Plot
```

---

# Notes

- Simulations utilize random mating
- Selection occurs entirely on GEBV
- Missing genotype values are internally replaced prior to simulation
- Temporary files are automatically created and deleted internally
- Simulations are stochastic and results will vary slightly between runs
