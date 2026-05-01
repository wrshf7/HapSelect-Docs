# Parent Selection

HapSelect uses a genetic algorithm (GA) to select an optimal set of founders that maximizes coverage of high-value haplotype alleles across target haploblocks.

## Preparing Input

Select top haploblocks by `Block_Var` and extract the corresponding localGEBV matrix:

```r
haploblock_effects <- haploblock_obj$Haploblocks
haploblock_effects <- haploblock_effects[order(haploblock_effects$Block_Var, decreasing = TRUE), ]
haploblock_top     <- haploblock_effects[1:15, ]  # top 15 blocks

localGEBV <- haploblock_obj$Haplotype_Effect_Matrix
localGEBV <- localGEBV[row.names(localGEBV) %in% haploblock_top$Block_ID, ]
localGEBV <- as.data.frame(t(as.matrix(localGEBV)))
```

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
|-----------|-------------|
| `localGEBV` | Individuals × blocks matrix of localGEBV values |
| `n_founders` | Number of parents to select |
| `popSize` | GA population size |
| `maxiter` | Maximum number of generations |
| `run` | Stop early if fitness does not improve for this many generations |
| `selfing` | Allow self-crosses in the candidate set |
| `pmutation` | Mutation probability |
| `pcrossover` | Crossover probability |
| `pelite` | Proportion of elite solutions carried forward |

## Output

```r
# One optimal set of selected parents
GA_output$One_Solution
```

The output contains the selected parent IDs and the GA fitness trace.
