# Function to calculate the empirical (marginal) reliability

Given secondary latent trait estimates and their associated standard
errors returned from
[`fscores`](https://philchalmers.github.io/mirt/reference/fscores.md),
compute the empirical reliability.

## Usage

``` r
empirical_rxx(Theta_SE, T_as_X = FALSE)
```

## Arguments

- Theta_SE:

  a matrix of latent trait estimates returned from
  [`fscores`](https://philchalmers.github.io/mirt/reference/fscores.md)
  with the options `full.scores = TRUE` and `full.scores.SE = TRUE`

- T_as_X:

  logical; should the observed variance be equal to
  `var(X) = var(T) + E(E^2)` or `var(X) = var(T)` when computing
  empirical reliability estimates? Default (`FALSE`) uses the former

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`fscores`](https://philchalmers.github.io/mirt/reference/fscores.md),
[`marginal_rxx`](https://philchalmers.github.io/mirt/reference/marginal_rxx.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

dat <- expand.table(deAyala)
itemstats(dat)
#> $overall
#>      N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  19601            2.912          1.434 0.233 0.074 0.608     0.898
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1 19601 2 0.887 0.316   0.447         0.246       0.605
#> Item.2 19601 2 0.644 0.479   0.688         0.439       0.510
#> Item.3 19601 2 0.566 0.496   0.680         0.416       0.523
#> Item.4 19601 2 0.427 0.495   0.673         0.405       0.529
#> Item.5 19601 2 0.387 0.487   0.602         0.312       0.581
#> 
#> $proportions
#>            0     1
#> Item.1 0.113 0.887
#> Item.2 0.356 0.644
#> Item.3 0.434 0.566
#> Item.4 0.573 0.427
#> Item.5 0.613 0.387
#> 
mod <- mirt(dat)

theta_se <- fscores(mod, full.scores.SE = TRUE)
empirical_rxx(theta_se)
#>        F1 
#> 0.6200703 

theta_se <- fscores(mod, full.scores.SE = TRUE, method = 'ML')
empirical_rxx(theta_se)
#>        F1 
#> 0.5636644 
empirical_rxx(theta_se, T_as_X = TRUE)
#>        F1 
#> 0.2258948 

# }
```
