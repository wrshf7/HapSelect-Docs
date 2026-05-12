# Parent Selection

HapSelect uses a genetic algorithm (GA) to select an optimal set of founders that maximizes coverage of high-value haplotype alleles across target haploblocks.

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

- `n = number` denotes the integer number of blocks to select with the greatest haploblock variance
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

Users are encouraged to experiment with multiple selection thresholds depending on breeding goals and computational resources. The [Basic Simulation](workflow/basic-simulation.md) functions will be useful for interpreting the GA output.


## `genetic_algorithm()`

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
