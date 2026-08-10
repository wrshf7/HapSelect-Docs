# Basic GA vs TS Simulation

The `localGEBV_vs_TS_simulation()` and `Haplotype_vs_TS_simulation` functions perform a simplified recurrent genomic selection simulation comparing:

- parents selected using the HapSelect genetic algorithm (GA) with localGEBV or haplotype effects
- parents selected using standard genomic truncation selection (TS)

The simulation evaluates long-term breeding value trajectories across generations and can optionally visualise population structure using PCA to demonstrate where parents are being sourced.

## Advantages of GA vs TS

In general, if a population has been under selection for a given trait, the best parents identified by whole-genome GEBV via TS tend to be related. In other words, they're good because they share the same good genomic segements; but conversely, they also share the same bad genomic segments. The GA with the default fitness functions, in comparison, will attempt to match indviduals by potential contribution. This necessarily means individuals tend to be less related if the population has been under selection for the given trait. This also means genetic variance tends to be greater in these individuals, even if the mean GEBV of the GA parents is initially worse than the TS parents. Thus, long-term selection potential is usually greater. This simulation attempts to demonstrate that diversity potential.

In the PCA, GA, TS, and overlapping parents are marked for visualisation. If the TS-selected parents tend to cluster to a small region of the PCA plot, it usually means that trait has been under selection in that population and such a scenario is where GA-selected parents shine. However, if TS are dispersed and thus unrelated, then realised differences between the GA and TS-selected parents tend to be minimal or even unfavourable.

---

## Overview

The simulation performs the following steps:

1. Extracts GA-selected parents from the GA output object (`haplotype_parent_selection()` or `local_gebv_parent_selection()` functions
2. Selects TS parents based on highest GEBV
3. Converts physical marker positions to proportional genetic positions in cM (if needed)
4. Simulates recurrent random crossing and selection across generations
5. Repeats the simulation multiple times to quantify Monte Carlo variability
6. Summarises breeding value trajectories
7. Generates trajectory and PCA visualisations of selected parents

The simulation is intentionally simple and designed primarily for:

- comparing parent selection strategies
- evaluating long-term gain potential
- exploring diversity retention
- demonstrating GA optimisation behaviour

---

## Running the Simulation

```r
#localGEBV
localGEBV_Sim  <- localGEBV_vs_TS_simulation(
  GA_output            = localGEBV_parent_obj,
  geno                 = geno,
  marker_effects       = marker_effects,
  map                  = map,
  genetic_map_position = NULL,
  num_gen              = 50,
  num_sim_reps         = 30,
  num_cross_per_gen    = 1000,
  num_TS_parents       = NULL,
  mean_adjust          = TRUE,
  maximize             = TRUE,
  max_cM_chr           = 100,
  PCA                  = TRUE,
  colors               = c("green", "#d95f02", "#A01FF0", "gray80"),
  alpha                = c(1,1,1,0.5)
)

#Haplotype
Haplotype_Sim  <- Haplotype_vs_TS_simulation(
  GA_output            = haplotype_parent_obj,
  geno_phased          = geno,
  marker_effects       = marker_effects,
  map                  = map,
  genetic_map_position = NULL,
  num_gen              = 50,
  num_sim_reps         = 30,
  num_cross_per_gen    = 100,
  num_TS_parents       = NULL,
  mean_adjust          = TRUE,
  maximize             = TRUE,
  max_cM_chr           = 100,
  PCA                  = TRUE,
  colors               = c("green", "#d95f02", "#A01FF0", "gray80"),
  alpha                = c(1,1,1,0.5)
)

```

---

## Required Inputs

| Parameter | Description |
|:---|:---|
| `GA_output` | Output object returned from `local_gebv_parent_selection()` or `haplotype_parent_selection()` |
| `geno` | Genotype/Haplotype dataframe used throughout the HapSelect workflow |
| `marker_effects` | Marker effects dataframe |
| `map` | Ordered map dataframe |

---

## Optional Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `genetic_map_position` | `NULL` | Optional vector of genetic map positions in cM - defined proportional to physical position if `NULL` |
| `num_gen` | `50` | Number of recurrent selection generations |
| `num_sim_reps` | `30` | Number of simulation replicates |
| `num_cross_per_gen` | `1000` | Number of progeny generated per generation by random mating |
| `num_TS_parents` | `NULL` | Number of truncation-selected parents, default is the same as supplied GA parents |
| `mean_adjust` | `TRUE` | Whether to internally center markers for TS GEBV |
| `maximize` | `TRUE` | If `TRUE`, assumes larger GEBV is better. If `FALSE`, tries to minimize GEBV |
| `max_cM_chr` | `100` | Per-chromosome genetic length assumption if no genetic map is supplied |
| `PCA` | `TRUE` | Whether to compute PCA visualisation |
| `colors` | `c("green", "#d95f02", "#A01FF0", "gray80")` | Plot colors for GA, TS, Overlap, and Not Selected individuals |
| `alpha` | `c(1,1,1,0.5)` | PCA transparency values for GA, TS, Overlap, and Not Selected individuals |

---

## Simulation Workflow

### 1. Marker Compatibility Checks

The function internally verifies that:

- marker IDs match between `geno`, `marker_effects`, and `map`
- marker IDs contain no duplicates
- marker IDs contain no missing values

The simulation will stop immediately if incompatibilities are detected.

---

### 2. GEBV Calculation

The function computes genomic estimated breeding values (GEBV) using:

```r
GEBV = Zu
```

where:

- `Z` is the genotype matrix
- `u` is the marker effect vector

If `mean_adjust = TRUE`, markers are internally centered before GEBV calculation.

This should remain `TRUE` in almost all analyses, but in practice will not affect the simulation because it will not change GEBV ranking.

---

### 3. TS Parent Selection

Truncation-selected parents are chosen as the individuals with the highest (or lowest when `maximize = FALSE` GEBV values.


Example:

```r
num_TS_parents = 20
```

selects the top 20 individuals ranked by GEBV.

If `num_TS_parents = NULL`, the number of TS parents automatically matches the number of GA parents.


---

### 5. Genetic Map Construction

The simulation requires marker positions in centiMorgans (cM).

Two approaches are supported.

---

#### Using a True Genetic Map (Recommended)

```r
#example of how to supply the genetic map positions
genetic_map_position = map$cM
```

The vector must:

- be the same length as `map`
- contain no missing values
- be in the same marker order as `map`

This is the most biologically realistic option.

---

### Inferring Genetic Positions from Physical Distance

If:

```r
genetic_map_position = NULL
```

then marker positions are inferred proportionally from physical distance relative to the first and last marker positions (i.e., the first marker is 0 cM and the last is 100 cM with the rest based on proportion of distance between the first and last marker).

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

## Recurrent Selection Simulation

The simulation uses the `genomicSimulation` package internally.

For each generation:

1. Mean breeding value is recorded
2. Top individuals are selected by GEBV
3. Random crosses generate the next generation
4. The process repeats for `num_gen` generations

This procedure is independently repeated `num_sim_reps` times to quantify stochastic simulation variability. Variability is derived from recombination differences, random mating sampling, Mendelian sampling (random inheritance), and, in the case of unphased data, random phasing of heterozygotes initially.

---

## Understanding Simulation Parameters

### `num_gen`

Controls the number of recurrent selection generations.

#### Larger Values

- better evaluate long-term gain
- reveal selection plateaus
- better assess diversity preservation

#### Smaller Values

- faster simulations
- emphasize short-term gain

---

### `num_sim_reps`

Controls the number of independent simulation replicates.

Replicates differ because of:

- recombination randomness
- stochastic inheritance
- random mating patterns
- random phasing (if using dosage heterozygotes)

#### Larger Values

Advantages:

- smoother trajectories
- reduced stochastic noise
- more stable estimates

Disadvantages:

- slower runtime
- increased memory usage

---

### `num_cross_per_gen`

Controls the number of progeny generated each generation.

#### Larger Values

Advantages:

- stronger selection intensity
- greater opportunity for favourable recombination
- smoother trajectories
- more stable representation of genomic variability (i.e., less drift)

Disadvantages:

- slower simulations
- increased memory usage

---

### `mean_adjust`

Controls genotype centering before GEBV calculation.

Recommended setting:

```r
mean_adjust = TRUE
```

This should only be disabled if you know what you are doing! However, in practice, it affects very little in this step.

!!! Warning
    Centering is currently not conducted in GenomicSimulation simulations - therefore, BV may be mean shifted!

---

### `PCA`

If `TRUE`, PCA is performed on the genotype matrix.

The PCA visualisation highlights:

- GA-selected parents
- TS-selected parents
- overlapping parents
- non-selected individuals

This helps visualise:

- diversity retention
- population structure
- overlap between strategies
- presence of selection pressure

---

## PCA Plot Colors

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
    Even if `PCA = FALSE`, this should still be specified with 4 colors as the first two are utilised simultaneously for the `Simulation_Plot`

---

## PCA Transparency

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

## Simulation Outputs

If:

```r
PCA = TRUE
```

the function returns:

```r
list(
  Simulation_Plot,
  PCA_Plot,
  Simulation_Summary,
  PCA_df
)
```

---

### `Simulation_Plot`

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

### `PCA_Plot`

A PCA visualisation of the genotype matrix showing selected parents.

Display:

```r
parent_sln_obj$PCA_Plot
```

---

### `Simulation_Summary`

A data frame in long-format giving the mean and standard deviation of simulation replicates at each generation for GA and TS.

Access:

```r
parent_sln_obj$Simulation_Summary
```

---

### `PCA_df`

A dataframe containing:

- accession names
- PC coordinates
- parent selection groupings

Access:

```r
head(parent_sln_obj$PCA_df)
```

---

## Notes

- Simulations utilise random mating
- Selection each generation occurs entirely on GEBV (TS)
- Missing genotype values are internally replaced prior to simulation (currently as heterozygous and randomly phased by GenomicSimulation)
- Temporary files are automatically created and deleted internally
- Simulations are stochastic and results will vary slightly between runs
- Imputing before beginning HapSelect and providing a phased input for simulation (haplotype genotype matrix and function call, even if parent selection was with localGEBV) will reduce some of the stochastic behavior between runs
