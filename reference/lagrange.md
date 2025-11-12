# Lagrange test for freeing parameters

Lagrange (i.e., score) test to test whether parameters should be freed
from a more constrained baseline model.

## Usage

``` r
lagrange(mod, parnum, SE.type = "Oakes", type = "Richardson", ...)
```

## Arguments

- mod:

  an estimated model

- parnum:

  a vector, or list of vectors, containing one or more parameter
  locations/sets of locations to be tested. See objects returned from
  [`mod2values`](https://philchalmers.github.io/mirt/reference/mod2values.md)
  for the locations

- SE.type:

  type of information matrix estimator to use. See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  further details

- type:

  type of numerical algorithm passed to
  [`numerical_deriv`](https://philchalmers.github.io/mirt/reference/numerical_deriv.md)
  to obtain the gradient terms

- ...:

  additional arguments to pass to
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`wald`](https://philchalmers.github.io/mirt/reference/wald.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
dat <- expand.table(LSAT7)
mod <- mirt(dat, 1, 'Rasch')
(values <- mod2values(mod))
#>    group   item     class   name parnum value lbound ubound   est const nconst
#> 1    all Item.1      dich     a1      1 1.000   -Inf    Inf FALSE  none   none
#> 2    all Item.1      dich      d      2 1.868   -Inf    Inf  TRUE  none   none
#> 3    all Item.1      dich      g      3 0.000      0      1 FALSE  none   none
#> 4    all Item.1      dich      u      4 1.000      0      1 FALSE  none   none
#> 5    all Item.2      dich     a1      5 1.000   -Inf    Inf FALSE  none   none
#> 6    all Item.2      dich      d      6 0.791   -Inf    Inf  TRUE  none   none
#> 7    all Item.2      dich      g      7 0.000      0      1 FALSE  none   none
#> 8    all Item.2      dich      u      8 1.000      0      1 FALSE  none   none
#> 9    all Item.3      dich     a1      9 1.000   -Inf    Inf FALSE  none   none
#> 10   all Item.3      dich      d     10 1.461   -Inf    Inf  TRUE  none   none
#> 11   all Item.3      dich      g     11 0.000      0      1 FALSE  none   none
#> 12   all Item.3      dich      u     12 1.000      0      1 FALSE  none   none
#> 13   all Item.4      dich     a1     13 1.000   -Inf    Inf FALSE  none   none
#> 14   all Item.4      dich      d     14 0.521   -Inf    Inf  TRUE  none   none
#> 15   all Item.4      dich      g     15 0.000      0      1 FALSE  none   none
#> 16   all Item.4      dich      u     16 1.000      0      1 FALSE  none   none
#> 17   all Item.5      dich     a1     17 1.000   -Inf    Inf FALSE  none   none
#> 18   all Item.5      dich      d     18 1.993   -Inf    Inf  TRUE  none   none
#> 19   all Item.5      dich      g     19 0.000      0      1 FALSE  none   none
#> 20   all Item.5      dich      u     20 1.000      0      1 FALSE  none   none
#> 21   all  GROUP GroupPars MEAN_1     21 0.000   -Inf    Inf FALSE  none   none
#> 22   all  GROUP GroupPars COV_11     22 1.022      0    Inf  TRUE  none   none
#>    prior.type prior_1 prior_2
#> 1        none     NaN     NaN
#> 2        none     NaN     NaN
#> 3        none     NaN     NaN
#> 4        none     NaN     NaN
#> 5        none     NaN     NaN
#> 6        none     NaN     NaN
#> 7        none     NaN     NaN
#> 8        none     NaN     NaN
#> 9        none     NaN     NaN
#> 10       none     NaN     NaN
#> 11       none     NaN     NaN
#> 12       none     NaN     NaN
#> 13       none     NaN     NaN
#> 14       none     NaN     NaN
#> 15       none     NaN     NaN
#> 16       none     NaN     NaN
#> 17       none     NaN     NaN
#> 18       none     NaN     NaN
#> 19       none     NaN     NaN
#> 20       none     NaN     NaN
#> 21       none     NaN     NaN
#> 22       none     NaN     NaN

# test all fixed slopes individually
parnum <- values$parnum[values$name == 'a1']
lagrange(mod, parnum)
#>            X2 df          p
#> 1  0.36714267  1 0.54456589
#> 5  0.04789243  1 0.82677204
#> 9  4.69717570  1 0.03021223
#> 13 1.44960175  1 0.22859187
#> 17 2.55408012  1 0.11000983

# compare to LR test for first two slopes
mod2 <- mirt(dat, 'F = 1-5
                   FREE = (1, a1)', 'Rasch')
coef(mod2, simplify=TRUE)$items
#>              a1         d g u
#> Item.1 1.157956 1.9393358 0 1
#> Item.2 1.000000 0.7850452 0 1
#> Item.3 1.000000 1.4501952 0 1
#> Item.4 1.000000 0.5175858 0 1
#> Item.5 1.000000 1.9787464 0 1
anova(mod, mod2)
#>           AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> mod  5341.802 5352.192 5352.994 5371.248 -2664.901               
#> mod2 5343.264 5355.386 5356.321 5377.618 -2664.632 0.538  1 0.463

mod2 <- mirt(dat, 'F = 1-5
                   FREE = (2, a1)', 'Rasch')
coef(mod2, simplify=TRUE)$items
#>               a1         d g u
#> Item.1 1.0000000 1.8746626 0 1
#> Item.2 0.9464702 0.7810629 0 1
#> Item.3 1.0000000 1.4661198 0 1
#> Item.4 1.0000000 0.5234068 0 1
#> Item.5 1.0000000 1.9997347 0 1
anova(mod, mod2)
#>           AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> mod  5341.802 5352.192 5352.994 5371.248 -2664.901               
#> mod2 5343.720 5355.842 5356.777 5378.075 -2664.860 0.081  1 0.775

mod2 <- mirt(dat, 'F = 1-5
                   FREE = (3, a1)', 'Rasch')
coef(mod2, simplify=TRUE)$items
#>              a1         d g u
#> Item.1 1.000000 1.8101732 0 1
#> Item.2 1.000000 0.7649559 0 1
#> Item.3 1.848252 1.7774916 0 1
#> Item.4 1.000000 0.5042538 0 1
#> Item.5 1.000000 1.9316244 0 1
anova(mod, mod2)
#>           AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> mod  5341.802 5352.192 5352.994 5371.248 -2664.901               
#> mod2 5335.130 5347.252 5348.187 5369.485 -2660.565 8.671  1 0.003

# test slopes first two slopes and last three slopes jointly
lagrange(mod, list(parnum[1:2], parnum[3:5]))
#>                X2 df          p
#> 1.5     0.4591775  2 0.79486042
#> 9.13.17 9.1527189  3 0.02732783

# test all 5 slopes and first + last jointly
lagrange(mod, list(parnum[1:5], parnum[c(1, 5)]))
#>                   X2 df          p
#> 1.5.9.13.17 9.861713  5 0.07924974
#> 1.17        2.898233  2 0.23477762

# }
```
