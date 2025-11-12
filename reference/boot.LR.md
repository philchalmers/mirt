# Parametric bootstrap likelihood-ratio test

Given two fitted models, compute a parametric bootstrap test to
determine whether the less restrictive models fits significantly better
than the more restricted model. Note that this hypothesis test also
works when prior parameter distributions are included for either model.
Function can be run in parallel after using a suitable
[`mirtCluster`](https://philchalmers.github.io/mirt/reference/mirtCluster.md)
definition.

## Usage

``` r
boot.LR(mod, mod2, R = 1000, verbose = interactive())
```

## Arguments

- mod:

  an estimated model object, more constrained than `mod2`

- mod2:

  an estimated model object

- R:

  number of parametric bootstraps to use.

- verbose:

  logical; include additional information in the console?

## Value

a p-value evaluating whether the more restrictive model fits
significantly worse than the less restrictive model

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

# standard
dat <- expand.table(LSAT7)
mod1 <- mirt(dat, 1)
mod2 <- mirt(dat, 1, '3PL')

# standard LR test
anova(mod1, mod2)
#>          AIC    SABIC       HQ      BIC    logLik  X2 df     p
#> mod1 5337.61 5354.927 5356.263 5386.688 -2658.805             
#> mod2 5346.11 5372.085 5374.089 5419.726 -2658.055 1.5  5 0.913

# bootstrap LR test (run in parallel to save time)
if(interactive()) mirtCluster()
boot.LR(mod1, mod2, R=200)
#> [1] 0.3084577

# }
```
