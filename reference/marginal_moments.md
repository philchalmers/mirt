# Function to calculate the marginal moments for items and bundles fit via IRT models

Given an estimated model and a prior density function, compute the
marginal moments for either and item or a bundle of items. Function
returns the first found moments implied by the model and select density
function (MEAN, VAR, SKEW, and KURT). Currently limited to
unidimensional IRT models.

## Usage

``` r
marginal_moments(
  mod,
  which.items = NULL,
  group = NULL,
  bundle = TRUE,
  density = NULL,
  Theta_lim = c(-6, 6),
  quadpts = 121,
  ...
)
```

## Arguments

- mod:

  an object of class `'SingleGroupClass'` or `'MultipleGroupClass'`

- which.items:

  vector indicating which items to use in the computation of the
  expected values. Default (`NULL`) uses all available items

- group:

  optional indicator to return only specific group information for
  multiple group models. Default compute moments for each group,
  returning a `list`

- bundle:

  logical; given `which.items`, should the composite of the response
  functions be used as a collective bundle? If `TRUE` then
  [`expected.test`](https://philchalmers.github.io/mirt/reference/expected.test.md)
  will be used, otherwise
  [`expected.item`](https://philchalmers.github.io/mirt/reference/expected.item.md)
  will be used

- density:

  a density function to use for integration. Default assumes the latent
  traits are from a normal (Gaussian) distribution. Function definition
  must be of the form `function(quadrature, mean, sd)` as the values of
  the mean/variance are extracted and passed from the supplied model

- Theta_lim:

  range of integration grid to use when forming expectations

- quadpts:

  number of discrete quadrature to use in the computations

- ...:

  additional arguments passed to the density function

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`expected.item`](https://philchalmers.github.io/mirt/reference/expected.item.md),
[`expected.test`](https://philchalmers.github.io/mirt/reference/expected.test.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

# single group
dat <- expand.table(deAyala)
mod <- mirt(dat)
TS <- rowSums(dat)

# expected moments of total scores given model
marginal_moments(mod)
#>    MEAN   VAR   SKEW  KURT
#> 1 2.917 1.297 -0.253 2.032
c(mean=mean(TS), var=var(TS))  # first two moments of data
#>     mean      var 
#> 2.911790 2.056759 

# same, however focusing on individual items
marginal_moments(mod, bundle = FALSE)
#>         MEAN   VAR   SKEW  KURT
#> Item.1 0.887 0.015 -2.121 8.454
#> Item.2 0.646 0.089 -0.597 2.053
#> Item.3 0.568 0.075 -0.260 1.874
#> Item.4 0.428 0.075  0.275 1.888
#> Item.5 0.388 0.040  0.394 2.367
cbind(mean=colMeans(dat), var=apply(dat, 2, var)) # first two moments of data
#>             mean        var
#> Item.1 0.8874547 0.09988393
#> Item.2 0.6440488 0.22926165
#> Item.3 0.5659915 0.24565765
#> Item.4 0.4269680 0.24467881
#> Item.5 0.3873272 0.23731694

############################################
## same as above, however with multiple group model

set.seed(1234)
group <- sample(c('G1', 'G2'), nrow(dat), replace=TRUE)
modMG <- multipleGroup(dat, group=group,
 invariance=c(colnames(dat), 'free_mean', 'free_var'))
coef(modMG, simplify=TRUE)
#> $G1
#> $items
#>           a1      d g u
#> Item.1 1.244  2.588 0 1
#> Item.2 2.021  1.000 0 1
#> Item.3 1.572  0.397 0 1
#> Item.4 1.566 -0.413 0 1
#> Item.5 0.995 -0.548 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> 
#> $G2
#> $items
#>           a1      d g u
#> Item.1 1.244  2.588 0 1
#> Item.2 2.021  1.000 0 1
#> Item.3 1.572  0.397 0 1
#> Item.4 1.566 -0.413 0 1
#> Item.5 0.995 -0.548 0 1
#> 
#> $means
#>     F1 
#> -0.004 
#> 
#> $cov
#>       F1
#> F1 1.003
#> 
#> 

# expected moments of total scores given model
marginal_moments(modMG)
#> $G1
#>   MEAN   VAR   SKEW  KURT
#> 1 2.92 1.295 -0.255 2.034
#> 
#> $G2
#>    MEAN   VAR   SKEW KURT
#> 1 2.915 1.299 -0.251 2.03
#> 
marginal_moments(modMG, group = 'G1') # specific group only
#>   MEAN   VAR   SKEW  KURT
#> 1 2.92 1.295 -0.255 2.034

# same, however focusing on individual items
marginal_moments(modMG, bundle = FALSE)
#> $G1
#>         MEAN   VAR   SKEW  KURT
#> Item.1 0.887 0.015 -2.124 8.476
#> Item.2 0.647 0.089 -0.600 2.058
#> Item.3 0.568 0.075 -0.262 1.876
#> Item.4 0.429 0.075  0.273 1.887
#> Item.5 0.388 0.040  0.393 2.366
#> 
#> $G2
#>         MEAN   VAR   SKEW  KURT
#> Item.1 0.887 0.015 -2.118 8.432
#> Item.2 0.645 0.089 -0.594 2.048
#> Item.3 0.567 0.075 -0.257 1.872
#> Item.4 0.428 0.075  0.277 1.889
#> Item.5 0.387 0.040  0.396 2.368
#> 


# }
```
