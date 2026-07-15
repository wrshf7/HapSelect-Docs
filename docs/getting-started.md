# Getting Started

This page walks through the complete HapSelect using localGEBV or haplotypes using the built-in example datasets.

## Load the Package

```r
library(HapSelect)
```

## Required Data

HapSelect expects two primary inputs (map file and genotype file). The rest can be computed in the R package or externally as additional inputs:

- **Map file** — marker positions with columns for marker ID, chromosome, and position.
    - First column should be named `SNP` and should be a character vector, not numeric or factorized.
    - Second column should be named `Chromosome` and should be numeric for proper sorting.
    - Third column should be named `Position` and should be numeric. It can be a physical position or genetic map position.

- **Genotype matrix (localGEBV)** — markers × individuals, coded as allele dosage (0/1/2 for diploid, 0/1/.../n for polyploid) at a biallelic marker.
    - First three columns should be identical to the map file.
    - Columns `4:ncol(genotype_file)` should be individuals and their genotype for each marker (rows are markers, columns are individuals).
    - Names of columns 4 onwards should be the individuals' identifiers.
    - Genotypes should be dosage format: i.e., number of copies of the alternative allele (integer counts only for consistency).
    - `NA` values are allowed and are handled differently at each stage via options in the relevant functions.

- **Phased genotype matrix (haplotypes)** — markers × chromosomes coded as the allele (0 for reference, 1 for alternative) at a biallelic marker.
    - First three columns should be identical to the map file.
    - Names of columns 4 onwards should be the individuals' chromosomal identifiers. The nomencalture adopted in `HapSelect` is `individualID_chromosomeNum` as the identifier, where the last `_` in the string is used to identify the chromosome.
    - Alleles should be dosage format: i.e., 0 for the reference allele and 1 for the alternative allele.
    - `NA` values are allowed and are handled differently at each stage via options in the relevant functions.

!!! warning
    Computing LD, haploblocks, and haplotype effects are supported for polyploids in localGEBV and haplotype methods, but parent selection algorithms and simulation currently only support diploid! We plan to update the localGEBV method to support polyploidy in the near future and the haplotype method at an indeterminate date.

- **LD file** - pairwise LD between each marker within a chromosome.
    - Can be computed internally using either the in-built function or the PLINK 1.9 wrapper function.
    - Can be computed externally and read in (see file structure below)
    - We HIGHLY recommend using the PLINK 1.9 wrapper function if PLINK 1.9 is installed, because R is not built for large, iterative computations required to compute LD pairs. Even small to moderate sized SNP panels will take too long in R.
    - Pairs not present (i.e., missing) in the data frame object are allowed and are handled in the haploblocking function call.
    - Columns:
        - `Chrom`: the chromosome each SNP pair belongs to.
        - `Locus1`: numerical integer for the first marker in the marker pair. This should correspond to the order of the marker in the **ordered map file** (see below for more details).
        - `Locus2`: similar to `Locus1`, this is the numerical integer of the second marker in the marker pair.
        - `Name1`: Character name of the first marker in the pair as seen in the genotype, map, and marker effects file.
        - `Name2`: Similar to `Name1`, this corresponds to the name of the second marker in the pair as seen in the genotype, map, and marker effects file.
        - `LD`: The numerical LD value computed. This is typically an $R^{2}$ value.

- **Marker effects file** — estimated SNP effects from a genomic prediction model
    - Can be computed for basic cases in the package (localGEBV genotype dosage matrix only!) with BLUE/dergressed BLUP/singular adjusted phenotype or provided externally.
    - First column: `SNP` corresponding to the `SNP` column in the map and genotype file. It should be formatted as a character vector.
    - Second column: `Effect` the allele substitution (marker) effects corresponding to the SNP in column one. This is the allele subsitution effect when subsituting one reference (0) allele with one alternative (1) allele.

!!! warning
    Ensure that the genotype coding in the genotype matrix is the same that was utilised to estimate marker effects. Marker effects computed within the package will be the same as long as the genotype matrix was not changed post-estimation.

- **Haploblock file** - dataframe produced as a standard part of the workflow.
    - May be provided if custom haploblocking is desired (e.g., blocking based on physical/map distance, number of markers per block, etc.)
        - We plan to provide haploblocking by distance/number of markers functionality in a future update.
    - Columns:
        - `Block`: SNP names within a block as identified in the map, genotype, and marker effects files. They should be separated by a ";" (no whitespaces as this will break downstream computations!). It must be a character!
        - `Block_ID`: Unique block identifier. We utilise the format "chromosome:block_nested_within_chromosome". As an example, "3:132" corresponds to the 132nd block on chromosome 3. It must be a character!
        - `Num_SNP`: The number of SNP in the block (number of markers in `Block`, must be numeric).
        - `First_SNP`: The first marker (based on physical/map position) in the block (must be a character).
        - `Last_SNP`: The last marker (based on physical/map position) in the block (must be a character).
        - `Chrom`: Numerical value of the chromosome the block resides on.
        - `Start_Pos`: The numerical  position of `First_SNP` (usually physical or genetic map position).
        - `End_Pos`: The numerical position of `Last_SNP` (usually physical or genetic map position).

Example datasets are bundled with the package along with function examples to compute other files:

```r
#Two required, primary inputs:
map  <- HapSelect::map
geno <- HapSelect::geno #localGEBV
geno <- HapSelect::geno_phased #haplotype method

#Optional inputs that can also be created in the R package workflow:
ld_pairs       <- HapSelect::ld_pairs
marker_effects <- HapSelect::marker_effects
BLUE           <- HapSelect::BLUE #example file to compute marker effects

#Haploblock dataframe can be computed by running the def_blocks() and block_obj_to_df() functions on the example data

```

## Step 1 — Order the Map

Markers must be sorted by chromosome and position before LD computation:

```r
map <- order_map(map = map)
```

## Step 2 — Compute Pairwise LD

```r
# Internal R implementation (slow for large datasets)
ld_pairs <- pairwise_ld(geno, parallelize = FALSE)

# Recommended: PLINK-backed implementation (recommended)
ld_pairs <- plink_pairwise_ld_geno(geno = geno)
```

!!! tip
    For large marker panels, use `plink_pairwise_ld_geno()`. See [Pairwise LD](workflow/pairwise-ld.md) for details.

## Step 3 — Define Haplotype Blocks

```r
haploblocks <- def_blocks(
  ld        = ld_pairs,
  map       = map,
  method    = "flanking", #compares new marker to adjacent marker only
  threshold = 0.2,        #LD threshold to include a marker
  tolerance = 4,          #how many markers can fail the check for inclusion
  tol_reset = TRUE,       #if a marker is added, reset the tolerance counter
  start     = "beginning",       #blocking starts at the highest LD pair
  parallel  = FALSE       #parallelise for large datasets (~25k+ markers)
)

haploblocks <- block_obj_to_df(haploblocks, map)

#Report basic block statistics:
block_summary(block_df = haploblocks)

#Note: Mean_Block_Size does not include singleton blocks (size 0) in the computation. It is the mean size of multi-SNP blocks.
```

## Step 4 — Compute localGEBV/haplotype effects

```r
#localGEBV
haploblock_obj <- compute_local_GEBV(
  geno           = geno_phased,
  marker_effects = marker_effects,
  haploblocks_df = haploblocks,
  set_missing_NA = TRUE,            #if TRUE, any missing genotypes in a block will return NA for the localGEBV
  mean_adjust    = TRUE             #centeres the genotype matrix - this should be TRUE in almost every case! Setting to FALSE will lead to biased localGEBV!
)

#haplotypes
haploblock_obj <- compute_haplotype_effects(
  geno           = geno,
  marker_effects = marker_effects,
  haploblocks_df = haploblocks,
  set_missing_NA = TRUE,            #if TRUE, any missing alleles in a block will return NA for the haplotype effect
  mean_adjust    = TRUE             #centeres the allele matrix - this should be TRUE in almost every case! Setting to FALSE will lead to biased haplotype effects!
)

```

!!! tip
    Check out the full workflow if you need a basic model to compute marker effects (localGEBV genotype matrix only!).

## Step 5 — Visualise

```r
#Manhattan-style plot of allele subsitution effects
marker_effects_plot(marker_effects = marker_effects$Effect,
                    chr = map$Chromosome, pos = map$Position)

#Manhattan-style plots of localGEBV and haplotype effects
unique_haplo_effects_plot(haplo_obj = haploblock_obj)     #haplotypes
unique_localGEBV_effects_plot(haplo_obj = haploblock_obj) #localGEBV

#Visualization of block sizes and block positions
plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks)

#Tornado plots demonstrating the spread of localGEBV and haplotype effects with blocks ordered by variance
haplo_block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE) #haplotypes
local_gebv_block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE)

#Plot to visualize LD decay
plot_ld_decay(map = map, ld = ld_pairs, max_kb = 500) #warning: may be slow for large dataset! Use method = "exp" for faster computation and monotonic decay. For a more flexible and still relatively quick approach, explore, method = "gam_tp".
```

!!! tip
    All plots return a ggplot object. Positions for localGEBV/haplotype effects, haploblocks, and marker densities are in absolute units to accomodate both physical and genetic maps. All ggplots can be modified to change axis scale, legend name, axis name, etc. See [Visualisations](workflow/visualizations.md) for examples.

## Step 6 — Select Top Blocks and Parents

!!! tip
    Generally, a small subset of blocks explains the vast majority of total variance in localGEBV. While the whole genome can be used for optimisation, results are often completley identical and it is much more computationally efficient to use a subset of blocks. In past examples, 9-13% of blocks have explaind at least 90% of genetic variance.

```r
There are three methods we have implemented to select haploblocks:
#1: select the top n blocks based on Block_Var
#2: select the top x% of blocks (round up)
#3: select the top blocks explaining x% of the total sum of block variance (round up)


#1: select top 15 haploblocks (arbitrary)

haploblock_obj = select_top_blocks(haploblock_obj = haploblock_obj, n = 15) #n denotes the number of blocks
#Objects are added to the haploblock_obj for the GA


#2 select top 50% of haploblocks
haploblock_obj = select_top_blocks(haploblock_obj = haploblock_obj, perc_total = 0.5) #perc_total is the percentage of blocks to select (rounded up)

#3 select the top blocks explaining at least 90% of the total block variance
haploblock_obj = select_top_blocks(haploblock_obj = haploblock_obj, perc_of_total_var = 0.9) #perc_of_total_var is the proportion of blocks exlaining at lest a given threshold of the total block variation (rounded up)

#Parent Selection (localGEBV)
GA_output = local_gebv_parent_selection(
    haploblock_obj = haploblock_obj,    #haploblock object from select_top_blocks()
    n_founders = 20,                    #number of parents to select
    popSize = 10,                       #parameter to aid optimization - higher uses more memory and each iteration is slower, but leads to faster convergence
    maxiter = 300,                      #max number of iterations
    run = 150,                          #max number of iterations of no change before termination
    pmutation = 0.6,                    #probability of swapping one idividual randomly in a population each iteration - aids in parameter space exploration
    pcrossover = 0.6,                   #probability of exchanging members between two populations in each iteration - aids in parameter space exploration
    maximize = TRUE,                    #if TRUE, maximise potential GEBV, if FALSE, minimise
    monitor = TRUE,                     #if TRUE, print out optimisation statisitcs each iteration
    strategy = "no_selfing"             #"no_selfing" does not allow selfing in optimisation (i.e., two distinct parents per block) whereas "selfing" does allow the same parent to be utilized twice per block
)

#Parent Selection (haplotypes - diploid only!)
GA_output = ohs_parent_selection(
    haploblock_obj = haploblock_obj,    #haploblock object from select_top_blocks()
    n_founders = 20,                    #number of parents to select
    popSize = 10,                       #parameter to aid optimization - higher uses more memory and each iteration is slower, but leads to faster convergence
    maxiter = 300,                      #max number of iterations
    run = 150,                          #max number of iterations of no change before termination
    pmutation = 0.6,                    #probability of swapping one idividual randomly in a population each iteration - aids in parameter space exploration
    pcrossover = 0.6,                   #probability of exchanging members between two populations in each iteration - aids in parameter space exploration
    maximize = TRUE,                    #if TRUE, maximize potential GEBV, if FALSE, minimize
    monitor = TRUE,                     #if TRUE, print out optimisation statisitcs each iteration
    strategy = "OHS"                    #"OHS" enforces haplotypes must come from two different parents, "OPV" allows the same haplotype to be utilised twice, "Haploid_OHS" allows the same parent to contribute two haplotypes, but not from the same chromosome
)

#one unique set of parents - other solutions may exist, see the GA_output$GA@solution object for other potential sets of parents with the same fitness (note: many are the same combination in different permutations)
GA_output$selected_founders

```

!!! warning
    Be sure to use `select_top_blocks()` to regenerate `haploblock_obj`. If you would like to use the whole genome, specify `n = nrow(haploblock_obj$Haploblocks)`, `perc_total = 1`, or `perc_of_total_var = 1`


## Step 7 - Basic Simulation and Diversity Analysis

```r

#Basic Recurrent Truncation Selection Simulation and Diversity Analysis (haplotypes)
parent_sln_obj = OHS_vs_TS_simulation(
    GA_output = GA_output,
    geno_phased = geno_phased,                             #phased genotype matrix from previous steps
    marker_effects = marker_effects,                       #marker effects dataframe previously used
    map = map,                                             #map dataframe previously used
    genetic_map_position = NULL,                           #vector of genetic map positions for each markers (in cM) - if NULL, map positions are assumed to be proportionally distributed to chromosomal physical position (0 to max_cM_chr)
    num_gen = 50,                                          #number of generations to simulate forward
    num_sim_reps = 10,                                     #number of simulation replicates to compute rate of genetic gain mean and standard deviation
    num_cross_per_gen = 100,                               #for each generation, the number of progreny should be produced from random mating
    num_TS_parents = NULL,                                 #number of truncation selected parents to use for comparison - if NULL, it is the same as the number of genetic algorithim selected parents
    mean_adjust = TRUE,                                    #should GEBV be centered? This is important in most cases to prevent a mean shift (bias) in simulation
    max_cM_chr = 100,                                      #if genetic_map_position = NULL, the max cM per chromosome used to assign marker map positions in cM proportional to physical distance
    PCA = TRUE,                                            #if TRUE, perform PCA and highlight GA vs TS selected parents
    maximize = TRUE,                                       #if TRUE, recurrent selection selects greater GEBV individuals each generation. If FALSE, it selects lower GEBV individuals
    colors = c("green", "#d95f02", "#A01FF0", "gray80"),   #the four colors assigned to GA-selected parents, TS-selected parents, overlapping parents (PCA only), and individuals not selected (PCA only)
    alpha = c(1,1,1,0.5)                                   #transparency values assigned to the same four groups on the PCA plot
)

#Basic Recurrent Truncation Selection Simulation and Diversity Analysis (localGEBV)
parent_sln_obj = localGEBV_vs_TS_simulation(
    GA_output = GA_output,
    geno = geno,                                           #genotype matrix from previous steps
    marker_effects = marker_effects,                       #marker effects dataframe previously used
    map = map,                                             #map dataframe previously used
    genetic_map_position = NULL,                           #vector of genetic map positions for each markers (in cM) - if NULL, map positions are assumed to be proportionally distributed to chromosomal physical position (0 to max_cM_chr)
    num_gen = 50,                                          #number of generations to simulate forward
    num_sim_reps = 10,                                     #number of simulation replicates to compute rate of genetic gain mean and standard deviation
    num_cross_per_gen = 100,                               #for each generation, the number of progreny should be produced from random mating
    num_TS_parents = NULL,                                 #number of truncation selected parents to use for comparison - if NULL, it is the same as the number of genetic algorithim selected parents
    mean_adjust = TRUE,                                    #should GEBV be centered? This is important in most cases to prevent a mean shift (bias) in simulation
    max_cM_chr = 100,                                      #if genetic_map_position = NULL, the max cM per chromosome used to assign marker map positions in cM proportional to physical distance
    PCA = TRUE,                                            #if TRUE, perform PCA and highlight GA vs TS selected parents
    maximize = TRUE,                                       #if TRUE, recurrent selection selects greater GEBV individuals each generation. If FALSE, it selects lower GEBV individuals
    colors = c("green", "#d95f02", "#A01FF0", "gray80"),   #the four colors assigned to GA-selected parents, TS-selected parents, overlapping parents (PCA only), and individuals not selected (PCA only)
    alpha = c(1,1,1,0.5)                                   #transparency values assigned to the same four groups on the PCA plot
)

```

See the individual workflow pages for full parameter documentation.
