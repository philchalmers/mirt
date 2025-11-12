# Compute profiled-likelihood (or posterior) confidence intervals

Computes profiled-likelihood based confidence intervals. Supports the
inclusion of equality constraints. Object returns the confidence
intervals and whether the respective interval could be found.

## Usage

``` r
PLCI.mirt(
  mod,
  parnum = NULL,
  alpha = 0.05,
  search_bound = TRUE,
  step = 0.5,
  lower = TRUE,
  upper = TRUE,
  inf2val = 30,
  NealeMiller = FALSE,
  verbose = interactive(),
  ...
)
```

## Arguments

- mod:

  a converged mirt model

- parnum:

  a numeric vector indicating which parameters to estimate. Use
  [`mod2values`](https://philchalmers.github.io/mirt/reference/mod2values.md)
  to determine parameter numbers. If `NULL`, all possible parameters are
  used

- alpha:

  two-tailed alpha critical level

- search_bound:

  logical; use a fixed grid of values around the ML estimate to
  determine more suitable optimization bounds? Using this has much
  better behaviour than setting fixed upper/lower bound values and
  searching from more extreme ends

- step:

  magnitude of steps used when `search_bound` is `TRUE`. Smaller values
  create more points to search a suitable bound for (up to the lower
  bound value visible with
  [`mod2values`](https://philchalmers.github.io/mirt/reference/mod2values.md)).
  When upper/lower bounds are detected this value will be adjusted
  accordingly

- lower:

  logical; search for the lower CI?

- upper:

  logical; search for the upper CI?

- inf2val:

  a numeric used to change parameter bounds which are infinity to a
  finite number. Decreasing this too much may not allow a suitable bound
  to be located. Default is 30

- NealeMiller:

  logical; use the Neale and Miller (1997) approximation? Default is
  `FALSE`

- verbose:

  logical; include additional information in the console?

- ...:

  additional arguments to pass to the estimation functions

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Chalmers, R. P., Pek, J., & Liu, Y. (2017). Profile-likelihood
Confidence Intervals in Item Response Theory Models. *Multivariate
Behavioral Research, 52*, 533-550.
[doi:10.1080/00273171.2017.1329082](https://doi.org/10.1080/00273171.2017.1329082)

Neale, M. C. & Miller, M. B. (1997). The use of likelihood-based
confidence intervals in genetic models. *Behavior Genetics, 27*,
113-120.

## See also

[`boot.mirt`](https://philchalmers.github.io/mirt/reference/boot.mirt.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
if(interactive()) mirtCluster() #use all available cores to estimate CI's in parallel
dat <- expand.table(LSAT7)
mod <- mirt(dat, 1)

result <- PLCI.mirt(mod)
result
#>      Item class parnam parnum     value lower_2.5 upper_97.5 lower_conv
#> 1  Item.1  dich     a1      1 0.9879254 0.6705382  1.3819761       TRUE
#> 2  Item.1  dich      d      2 1.8560605 1.6203325  2.1474211       TRUE
#> 3  Item.2  dich     a1      5 1.0808847 0.7816287  1.4614714       TRUE
#> 4  Item.2  dich      d      6 0.8079786 0.6373152  0.9994390       TRUE
#> 5  Item.3  dich     a1      9 1.7058006 1.1965700  2.6063138       TRUE
#> 6  Item.3  dich      d     10 1.8042187 1.4754050  2.3765517       TRUE
#> 7  Item.4  dich     a1     13 0.7651853 0.5211915  1.0554405       TRUE
#> 8  Item.4  dich      d     14 0.4859966 0.3417041  0.6365806       TRUE
#> 9  Item.5  dich     a1     17 0.7357980 0.4551545  1.0555079       TRUE
#> 10 Item.5  dich      d     18 1.8545127 1.6438380  2.0976038       TRUE
#>    upper_conv
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE
#> 6        TRUE
#> 7        TRUE
#> 8        TRUE
#> 9        TRUE
#> 10       TRUE


mod2 <- mirt(Science, 1)
result2 <- PLCI.mirt(mod2)
result2
#>       Item  class parnam parnum      value  lower_2.5 upper_97.5 lower_conv
#> 1  Comfort graded     a1      1  1.0417547  0.7008509   1.453043       TRUE
#> 2  Comfort graded     d1      2  4.8641542  4.0112194   5.966782       TRUE
#> 3  Comfort graded     d2      3  2.6399417  2.2332766   3.115799       TRUE
#> 4  Comfort graded     d3      4 -1.4660135 -1.7996417  -1.171326       TRUE
#> 5     Work graded     a1      5  1.2259618  0.8942261   1.620993       TRUE
#> 6     Work graded     d1      6  2.9240027  2.4851419   3.430756       TRUE
#> 7     Work graded     d2      7  0.9011651  0.6307977   1.195238       TRUE
#> 8     Work graded     d3      8 -2.2665647 -2.6975967  -1.893557       TRUE
#> 9   Future graded     a1      9  2.2933717  1.5687257   3.986967       TRUE
#> 10  Future graded     d1     10  5.2339928  4.1279881   7.822343       TRUE
#> 11  Future graded     d2     11  2.2137728  1.6589420   3.416513       TRUE
#> 12  Future graded     d3     12 -1.9637062 -3.0256934  -1.453916       TRUE
#> 13 Benefit graded     a1     13  1.0949151  0.7659052   1.501997       TRUE
#> 14 Benefit graded     d1     14  3.3479196  2.8453841   3.940600       TRUE
#> 15 Benefit graded     d2     15  0.9916289  0.7267311   1.282508       TRUE
#> 16 Benefit graded     d3     16 -1.6882599 -2.0443917  -1.375841       TRUE
#>    upper_conv
#> 1        TRUE
#> 2        TRUE
#> 3        TRUE
#> 4        TRUE
#> 5        TRUE
#> 6        TRUE
#> 7        TRUE
#> 8        TRUE
#> 9        TRUE
#> 10       TRUE
#> 11       TRUE
#> 12       TRUE
#> 13       TRUE
#> 14       TRUE
#> 15       TRUE
#> 16       TRUE

# only estimate CI's slopes
sv <- mod2values(mod2)
parnum <- sv$parnum[sv$name == 'a1']
result3 <- PLCI.mirt(mod2, parnum)
result3
#>      Item  class parnam parnum    value lower_2.5 upper_97.5 lower_conv
#> 1 Comfort graded     a1      1 1.041755 0.7008509   1.453043       TRUE
#> 2    Work graded     a1      5 1.225962 0.8942261   1.620993       TRUE
#> 3  Future graded     a1      9 2.293372 1.5687257   3.986967       TRUE
#> 4 Benefit graded     a1     13 1.094915 0.7659052   1.501997       TRUE
#>   upper_conv
#> 1       TRUE
#> 2       TRUE
#> 3       TRUE
#> 4       TRUE

# }
```
