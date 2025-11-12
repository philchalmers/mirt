# Extract various elements from estimated model objects

A generic function to extract the internal objects from estimated
models.

## Usage

``` r
extract.mirt(x, what, item = 1, ...)
```

## Arguments

- x:

  mirt model of class 'SingleGroupClass', 'MultipleGroupClass',
  'MixedClass' or 'DiscreteGroupClass'

- what:

  a string indicating what to extract

- item:

  if necessary, which item to extract information from. Defaults to 1 if
  not specified

- ...:

  additional arguments to pass to
  [`summary()`](https://rdrr.io/r/base/summary.html) (e.g., for
  different factor rotations)

## Details

Objects which can be extracted from mirt objects include:

- logLik:

  observed log-likelihood

- logPrior:

  log term contributed by prior parameter distributions

- G2:

  goodness of fit statistic

- df:

  degrees of freedom

- p:

  p-value for G2 statistic

- RMSEA:

  root mean-square error of approximation based on G2

- CFI:

  CFI fit statistic

- TLI:

  TLI fit statistic

- AIC:

  AIC

- BIC:

  BIC

- SABIC:

  sample size adjusted BIC

- HQ:

  HQ

- LLhistory:

  EM log-likelihood history

- tabdata:

  a tabular version of the raw response data input. Frequencies are
  stored in `freq`

- freq:

  frequencies associated with `tabdata`

- K:

  an integer vector indicating the number of unique elements for each
  item

- mins:

  an integer vector indicating the lowest category found in the input
  `data`

- model:

  input model syntax

- method:

  estimation method used

- itemtype:

  a vector of item types for each respective item (e.g., 'graded',
  '2PL', etc)

- itemnames:

  a vector of item names from the input data

- factorNames:

  a vector of factor names from the model definition

- rowID:

  an integer vector indicating all valid row numbers used in the model
  estimation (when all cases are used this will be `1:nrow(data)`)

- data:

  raw input data of item responses

- covdata:

  raw input data of data used as covariates

- tabdatalong:

  similar to `tabdata`, however the responses have been transformed into
  dummy coded variables

- fulldatalong:

  analogous to `tabdatafull`, but for the raw input data instead of the
  tabulated frequencies

- EMhistory:

  if saved, extract the EM iteration history

- F:

  unrotated factor loadings matrix

- Frot:

  rotated factor loadings matrix

- exp_resp:

  expected probability of the unique response patterns

- survey.weights:

  if supplied, the vector of survey weights used during estimation (NULL
  if missing)

- converged:

  a logical value indicating whether the model terminated within the
  convergence criteria

- iterations:

  number of iterations it took to reach the convergence criteria

- nest:

  number of freely estimated parameters

- parvec:

  vector containing uniquely estimated parameters

- vcov:

  parameter covariance matrix (associated with parvec)

- condnum:

  the condition number of the Hessian (if computed). Otherwise NA

- constrain:

  a list of item parameter constraints to indicate which item parameters
  were equal during estimation

- Prior:

  prior density distribution for the latent traits

- thetaPosterior:

  posterior distribution for latent traits when using EM algorithm

- key:

  if supplied, the data scoring key

- nfact:

  number of latent traits/factors

- nitems:

  number of items

- ngroups:

  number of groups

- groupNames:

  character vector of unique group names

- group:

  a character vector indicating the group membership

- invariance:

  a character vector indicating `invariance` input from
  [`multipleGroup`](https://philchalmers.github.io/mirt/reference/multipleGroup.md)

- secondordertest:

  a logical indicating whether the model passed the second-order test
  based on the Hessian matrix. Indicates whether model is a potential
  local maximum solution

- SEMconv:

  logical; check whether the supplemented EM information matrix
  converged. Will be `NA` if not applicable

- time:

  estimation time, broken into different sections

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`extract.group`](https://philchalmers.github.io/mirt/reference/extract.group.md),
[`extract.item`](https://philchalmers.github.io/mirt/reference/extract.item.md),
[`mod2values`](https://philchalmers.github.io/mirt/reference/mod2values.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
mod <- mirt(Science, 1)

extract.mirt(mod, 'logLik')
#> [1] -1608.87
extract.mirt(mod, 'K')  # unique categories for each item
#> [1] 4 4 4 4

# multiple group model
grp <- rep(c('G1', 'G2'), each = nrow(Science)/2)
mod2 <- multipleGroup(Science, 1, grp)

grp1 <- extract.group(mod2, 1) #extract single group model
extract.mirt(mod2, 'parvec')
#>  [1]  0.8316800  4.8913029  2.5597789 -1.3592138  0.8289135  2.4772334
#>  [7]  0.7025447 -2.1135580  2.8809926  5.4820292  2.2126281 -2.5451325
#> [13]  0.6920437  3.0133748  0.8608953 -1.5516666  1.1955161  4.8900152
#> [19]  2.6929286 -1.5545256  1.7621713  3.6416308  1.1724290 -2.5499919
#> [25]  2.4324588  6.4188352  2.6336468 -1.8297283  1.4875546  3.8081350
#> [31]  1.1473714 -1.8343958
extract.mirt(grp1, 'parvec')
#>  [1]  0.8316800  4.8913029  2.5597789 -1.3592138  0.8289135  2.4772334
#>  [7]  0.7025447 -2.1135580  2.8809926  5.4820292  2.2126281 -2.5451325
#> [13]  0.6920437  3.0133748  0.8608953 -1.5516666

# }
```
