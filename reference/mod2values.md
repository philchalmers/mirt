# Convert an estimated mirt model to a data.frame

Given an estimated model from any of mirt's model fitting functions this
function will convert the model parameters into the design data frame of
starting values and other parameter characteristics (similar to using
the `pars = 'values'` for obtaining starting values).

## Usage

``` r
mod2values(x)
```

## Arguments

- x:

  an estimated model x from the mirt package

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`extract.mirt`](https://philchalmers.github.io/mirt/reference/extract.mirt.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
dat <- expand.table(LSAT7)
mod <- mirt(dat, "F=1-5
                  CONSTRAIN=(1-5, a1)")
values <- mod2values(mod)
values
#>    group   item     class   name parnum value lbound ubound   est const nconst
#> 1    all Item.1      dich     a1      1 1.011   -Inf    Inf  TRUE     1   none
#> 2    all Item.1      dich      d      2 1.868   -Inf    Inf  TRUE  none   none
#> 3    all Item.1      dich      g      3 0.000      0      1 FALSE  none   none
#> 4    all Item.1      dich      u      4 1.000      0      1 FALSE  none   none
#> 5    all Item.2      dich     a1      5 1.011   -Inf    Inf  TRUE     1   none
#> 6    all Item.2      dich      d      6 0.791   -Inf    Inf  TRUE  none   none
#> 7    all Item.2      dich      g      7 0.000      0      1 FALSE  none   none
#> 8    all Item.2      dich      u      8 1.000      0      1 FALSE  none   none
#> 9    all Item.3      dich     a1      9 1.011   -Inf    Inf  TRUE     1   none
#> 10   all Item.3      dich      d     10 1.461   -Inf    Inf  TRUE  none   none
#> 11   all Item.3      dich      g     11 0.000      0      1 FALSE  none   none
#> 12   all Item.3      dich      u     12 1.000      0      1 FALSE  none   none
#> 13   all Item.4      dich     a1     13 1.011   -Inf    Inf  TRUE     1   none
#> 14   all Item.4      dich      d     14 0.521   -Inf    Inf  TRUE  none   none
#> 15   all Item.4      dich      g     15 0.000      0      1 FALSE  none   none
#> 16   all Item.4      dich      u     16 1.000      0      1 FALSE  none   none
#> 17   all Item.5      dich     a1     17 1.011   -Inf    Inf  TRUE     1   none
#> 18   all Item.5      dich      d     18 1.993   -Inf    Inf  TRUE  none   none
#> 19   all Item.5      dich      g     19 0.000      0      1 FALSE  none   none
#> 20   all Item.5      dich      u     20 1.000      0      1 FALSE  none   none
#> 21   all  GROUP GroupPars MEAN_1     21 0.000   -Inf    Inf FALSE  none   none
#> 22   all  GROUP GroupPars COV_11     22 1.000      0    Inf FALSE  none   none
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

# use the converted values as starting values in a new model, and reduce TOL
mod2 <- mirt(dat, 1, pars = values, TOL=1e-5)
coef(mod2, simplify=TRUE)
#> $items
#>           a1     d g u
#> Item.1 1.011 1.868 0 1
#> Item.2 1.011 0.791 0 1
#> Item.3 1.011 1.461 0 1
#> Item.4 1.011 0.521 0 1
#> Item.5 1.011 1.993 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 

# use parameters on different dataset
mod3 <- mirt(expand.table(LSAT6), pars=values)
coef(mod3, simplify=TRUE)
#> $items
#>           a1     d g u
#> Item_1 0.755 2.730 0 1
#> Item_2 0.755 0.999 0 1
#> Item_3 0.755 0.240 0 1
#> Item_4 0.755 1.307 0 1
#> Item_5 0.755 2.100 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 

# supports differing itemtypes on second model
sv <- mirt(Science, itemtype=c('graded', rep('gpcm', 3)), pars='values')
mod3 <- mirt(Science, pars = sv)  # itemtype omitted
coef(mod3, simplify=TRUE)$items
#>                a1       d1       d2         d3 ak0 ak1 ak2 ak3 d0
#> Comfort 1.0278666 4.842467 2.629427 -1.4610757  NA  NA  NA  NA NA
#> Work    0.8622915 1.735970 2.608453  0.8615093   0   1   2   3  0
#> Future  2.1388465 4.498623 6.612162  4.8108956   0   1   2   3  0
#> Benefit 0.7208997 2.096130 2.894803  1.7189648   0   1   2   3  0
extract.mirt(mod3, 'itemtype')
#> [1] "graded" "gpcm"   "gpcm"   "gpcm"  


# }
```
