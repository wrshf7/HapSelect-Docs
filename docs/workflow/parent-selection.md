# Parent Selection

## How the Genetic Algorithm Works

The HapSelect genetic algorithm (GA) attempts to identify a founder set that maximizes the potential to recover favorable haplotypes across the selected haploblocks.

Rather than optimizing overall GEBV directly, the GA optimizes the ability of the selected founder pool to produce highly favorable offspring combinations across genomic regions.

The GA core algorithm relies on the following R package while the fitness, mutation, and crossover functions are custom designed in HapSelect:

[GA: Genetic Algorithms](https://cran.r-project.org/web/packages/GA/index.html)

---

## The Optimization Objective

For each haploblock:

1. All possible pairwise crosses among the selected founders are evaluated
2. The expected offspring localGEBV for each cross is calculated
3. The highest-scoring cross for that haploblock is retained
4. The process repeats across all haploblocks
5. The fitness values are summed across haploblocks

The GA therefore attempts to maximize:

- favorable haplotype complementarity
- genomic coverage of elite haplotypes
- the best achievable offspring configuration across blocks

rather than simply selecting the individuals with the highest total GEBV.

---

## The Fitness Function

Conceptually, the fitness function is:

$$
\mathrm{Fitness} = \sum_{j=1}^{J} \max_{(i,k)} \left( \frac{localGEBV_{ij} + localGEBV_{kj}}{2} \right)
$$
where:

| Symbol | Meaning |
|:---|:---|
| *J* | Number of selected haploblocks |
| $localGEBV_{ij}$ | localGEBV of individual `i` at haploblock `j` |
| *(i,k)* | Pairwise founder combinations |

For each block, the GA identifies the founder pair with the highest expected offspring value (EPD) and sums these optimal values across all blocks.

---

## Why This Differs From Truncation Selection

Traditional truncation selection (TS):

- selects the individuals with the highest overall GEBV
- tends to repeatedly favor the same highly elite individuals
    - these individuals are often related (clustered on the PCA) and therefore share similar "good" and "bad" genomic segments. This increases inbreeding more quickly as well as unfavorable LD.
- may rapidly reduce diversity

In contrast, the HapSelect GA:

- searches for complementary founder combinations
- rewards founder sets that collectively cover favorable haplotypes
- allows different founders to contribute to different genomic regions
- can retain valuable rare haplotypes ignored by standard TS

This means an individual with only moderate total GEBV may still be highly valuable if it contributes an elite haplotype at a specific high-variance block.

---

## Pairwise Crossing Strategy

The GA assumes that:

- favorable haplotypes can be combined through recombination
- different founder pairs may be optimal for different haploblocks
- no single founder pair is necessarily optimal genome-wide

For each haploblock, the algorithm evaluates:

```r
combn(founders, 2)
```

to test all possible pairwise combinations among the selected founders.

If:

```r
selfing = TRUE
```

then self-crosses are also evaluated.

---

## Evolutionary Search Procedure

The GA evolves founder sets over multiple iterations using:

| Operation | Purpose |
|:---|:---|
| Population initialization | Generate random founder sets |
| Fitness evaluation | Score founder sets using haploblock complementarity |
| Mutation | Randomly replace founders |
| Crossover | Swap part of two founder sets |
| Elite sampling | Bias replacement toward high-performing solutions (greatest GEBV) during Mutation and if Crossover yields overlapping individuals |

The search attempts to balance:

- exploration of new founder combinations
- exploitation of high-performing founder sets

---

## Mutation

Mutation randomly replaces one founder within a solution:

```r
#expressed in probability: i.e., between 0 and 1
pmutation = 0.1
```

Higher mutation rates:

- increase exploration
- reduce risk of local optima

but may:

- destabilize convergence
- slow optimization if many "good" founder sets are systematically disrupted every generation

---

## Crossover

Crossover exchanges founders between two high-performing solutions.

The algorithm:

1. Retains part of each founder set
2. Combines non-overlapping founders
3. Fills missing founders from elite solutions based on `pelite`

This helps preserve useful founder combinations while still exploring new combinations.

---

## Elite Founder Sampling

The:

```r
pelite
```

parameter controls how strongly crossover favors founders from high-performing solutions.

Example:

```r
#expressed in probability: i.e., between 0 and 1
pelite = 0.2
```

means replacement founders are preferentially sampled from the top 20% of solutions ranked by overall fitness (GEBV).

Smaller values:

- increase selection pressure
- may accelerate convergence if elite individuals contain most of the best haplotypes
    - this may also cause convergence to a local optima if the value is too high

Larger values:

- maintain greater diversity
- improve exploration

---

## Interpretation of the Final Founder Set

The final GA-selected founder set should (generally) be interpreted as:

- a complementary breeding population
- a set of parents with strong collective haplotype coverage
- a founder pool optimized for long-term recombination potential

rather than simply the top individuals ranked by overall GEBV.

## Preparing Input - `select_top_blocks()`

Select top haploblocks by `Block_Var`. The `select_top_blocks()` function can select blocks in three different ways and returns a modified `haploblock_obj` object by adding on two additional dataframes to the list structure.

```r
#select the top n blocks using the n argument
haploblock_obj <- select_top_blocks(
  haploblock_obj = haploblock_obj,
  n = 15
)

#select the top x% of blocks using the perc_total argument
haploblock_obj <- select_top_blocks(
  haploblock_obj = haploblock_obj,
  perc_total = 0.5
)

#select the nubmer of blocks explaining at least x% of the total block variance utilizing the perc_total argument
haploblock_obj <- select_top_blocks(
  haploblock_obj = haploblock_obj,
  perc_of_total_var = 0.9
)
```

- `n = number` requires an integer number of blocks, denoted by `number` to select with the greatest haploblock variance
- `perc_total` denotes the proportion (between 0 and 1) of blocks to select (rounded up) from the total number of blocks with the greatest variance
- `perc_of_total_var` denotes selecting the proportion (between 0 and 1) of blocks (rounded up) that explain at least that proportion of the total block variance and have the greatest variance (i.e., the minimum number of blocks required to achieve the target proportion). This is the recommended option to use because it:
    - Is usually the most biologically meaningful approach.
    - Dynamically adapts to the architecture of the trait.
    - Retains more blocks for highly polygenic traits.
    - Retains fewer blocks when major-effect haploblocks dominate.

### Comparison of Selection Strategies

| Method | Strengths | Weaknesses |
| :--- | :--- | :--- |
| Top `n` blocks | Simple and interpretable | Arbitrary cutoff |
| Top percentage | Scales with dataset size | May retain weak blocks or discard high-variance blocks |
| Variance explained | Biologically adaptive | Number of retained blocks varies between traits and architecture |

The output is a modified `haploblock_obj` that contains subsetted matrices and dataframes needed for the GA.

---

### Computational Considerations

Increasing the number of retained haploblocks:

- increases optimization dimensionality
- increases GA runtime
- increases memory usage
- may slow convergence substantially

However, retaining too few blocks may:

- miss favorable rare haplotypes
- oversimplify trait architecture
- reduce long-term genetic gain potential
- may ignore how the contribution small, cumulative effects have to overall fitness and optimization

Users are encouraged to experiment with multiple selection thresholds depending on breeding goals and computational resources. The [Basic Simulation](workflow/basic-simulation.md) functions will be useful for interpreting the GA output.


## `genetic_algorithm()`

!!! warning
    These individuals are selected via a heuristic search optimization and are thus never guaranteed to be the best set of individuals! The heuristic search optimization is necessary to make most problems computationally feasible. Generally, unless stuck in a very pre-mature local optima, the results are the best or close to the best solution. Furthermore, more than one unique set of parents with the same overall fitness may exist in smaller scenarios. The `GA_output` object contains all solutions.


```r
GA_output <- genetic_algorithm(
  localGEBV  = localGEBV,
  n_founders = 20,
  popSize    = 10,
  maxiter    = 300,
  run        = 150,
  selfing    = FALSE,
  pmutation  = 0.2,
  pcrossover = 0.8,
  pelite     = 0.5
)
```

| Parameter | Description |
|---|---|
| `n_founders` | Integer number of parents to choose |
| `popSize` | Integer number of parental sets per simulation iteration |
| `maxiter` | Maximum iterations before termination |
| `run` | Iterations without improvement before stopping |
| `selfing` | Allow selfing (i.e., for the fitness function allow the same parent at a block). If `FALSE` requires two different parents per chosen block |
| `pmutation` | Mutation probability - swaps out one random individual for another from the total population |
| `pcrossover` | Crossover probability - swaps half of each population; if there is overlap, non-overlapping parents are chosen randomly from the total population |
| `pelite` | Elite proportion (between 0 and 1) - constrains choosing individuals from the total population to the highest `pelite` proportion based on GEBV when needed to find non-overlapping parents for `pcrossover` |

<table>
  <colgroup>
    <col style="width: 20%; white-space: nowrap;">
    <col style="width: 80%;">
  </colgroup>

  <thead>
    <tr>
      <th>Parameter</th>
      <th>Description</th>
    </tr>
  </thead>

  <tbody>
    <tr>
      <td><code>n_founders</code></td>
      <td>Integer number of parents to choose</td>
    </tr>

    <tr>
      <td><code>popSize</code></td>
      <td>Integer number of parental sets per simulation iteration</td>
    </tr>

    <tr>
      <td><code>maxiter</code></td>
      <td>Maximum iterations before termination</td>
    </tr>

    <tr>
      <td><code>run</code></td>
      <td>Iterations without improvement before stopping</td>
    </tr>

    <tr>
      <td><code>selfing</code></td>
      <td>
        Allow selfing (i.e., for the fitness function allow the same
        parent at a block). If <code>FALSE</code>, requires two different
        parents per chosen block.
      </td>
    </tr>

    <tr>
      <td><code>pmutation</code></td>
      <td>
        Mutation probability — swaps out one random individual for another
        from the total population.
      </td>
    </tr>

    <tr>
      <td><code>pcrossover</code></td>
      <td>
        Crossover probability — swaps half of each population; if there is
        overlap, non-overlapping parents are chosen randomly from the total
        population.
      </td>
    </tr>

    <tr>
      <td><code>pelite</code></td>
      <td>
        Elite proportion (between 0 and 1) — constrains choosing
        individuals from the total population to the highest
        <code>pelite</code> proportion based on GEBV when needed to find
        non-overlapping parents for <code>pcrossover</code>.
      </td>
    </tr>
  </tbody>
</table>

The genetic algorithm (GA) balances:

1. **Exploration**  
   Searching broadly across possible parental combinations

2. **Exploitation**  
   Refining highly fit parental combinations already discovered

Improper parameter tuning can lead to:

- premature convergence
- failure to converge
- excessive runtime
- oscillation around suboptimal solutions

---

### `popSize`

Controls the number of candidate parental sets evaluated per iteration.

#### Trade-offs

| Smaller `popSize` | Larger `popSize` |
|:---|:---|
| Faster individual iterations | Slower individual iterations |
| Less memory usage | Higher memory usage |
| Late convergence | Better exploration and earlier convergence |
| Higher risk of local optima | Lower risk of local optima |


Increasing `popSize` generally reduces the number of iterations needed for convergence, but each iteration becomes more computationally expensive.

---

### `maxiter`

Maximum number of iterations allowed.

#### Notes

- Too small → GA may terminate before convergence
- Too large → unnecessary runtime after convergence

Generally:

- small problems converge quickly
- highly polygenic architectures with many parents may require many iterations

---

### `run`

Number of iterations allowed without improvement before stopping.

#### Trade-offs

| Smaller `run` | Larger `run` |
|:---|:---|
| Faster termination | More exhaustive search |
| May stop too early | Longer runtime |
| Risk missing optimum | Better convergence stability |


---

### `pmutation`

Mutation probability.

Mutation randomly substitutes one individual within populations from the total population.

#### Trade-off

| Low mutation | High mutation |
|:---|:---|
| Stable convergence | More exploration |
| Risk local optima | Risk instability and "overshooting" |

##### Important Notes

Overly large mutation probabilities can prevent convergence entirely because high-performing parental sets are continuously disrupted.

---

### `pcrossover`

Probability populations exchange parental subsets (half of each pair swapped). If there are overlapping individuals after swapping, then the duplicates are dropped and non-duplicate individuals are randomly sampled from the total population.

#### Trade-offs

| Low crossover | High crossover |
|:---|:---|
| Less swapping | Greater swapping |
| Slower exploration | Faster exploration |
| More stable solutions | Greater instability |

Very high crossover rates may cause the GA to overshoot promising solutions and continuously disrupt near-optimal parental combinations.


---

### `pelite`

Restricts replacement individuals to the proportion of the population ranked by GEBV.

#### Trade-offs

| Smaller `pelite` | Larger `pelite` |
|:---|:---|
| Faster convergence | Greater diversity |
| Stronger selection pressure | Better exploration |
| Risk premature convergence | Slower convergence |

Very aggressive elite selection may reduce genetic diversity within the GA search process and increase the likelihood of local optima.

---

### If convergence is unstable:

- decrease `pmutation`
- decrease `pcrossover`
- increase `pelite`
- increase `run`

### If convergence is too slow:

- increase `popSize`
- slightly increase `pmutation`
- slightly increase `pcrossover`

### If solutions appear trapped in local optima:

- increase `popSize`
- increase `pmutation`
- increase `pelite`
- increase `run`

---

## Output

```r
# One optimal set of selected parents
GA_output$One_Solution
```

The output contains the selected parent IDs and names. More information is available in the internal `GA` object about GA performance, change in the best localGEBV and meanGEBV of the parental sets over iterations, etc. There may be more than one set of parents that give rise to the same optimum. Additional sets of parents are contained in the `GA` object and only the first set is presented in `$One_Solution`
