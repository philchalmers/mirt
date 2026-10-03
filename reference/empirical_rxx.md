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

Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
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
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found
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
