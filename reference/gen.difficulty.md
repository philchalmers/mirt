# Generalized item difficulty summaries

Function provides the four generalized item difficulty representations
for polytomous response models described by Ali, Chang, and Anderson
(2015). These estimates are used to gauge how difficult a polytomous
item may be.

## Usage

``` r
gen.difficulty(mod, type = "IRF", interval = c(-30, 30), ...)
```

## Arguments

- mod:

  a single factor model estimated by
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)

- type:

  type of generalized difficulty parameter to report. Can be `'IRF'` to
  use the item response function (default), `'mean'` to find the average
  of the difficulty estimates, `'median'` the median of the difficulty
  estimates, and `'trimmed'` to find the trimmed mean after removing the
  first and last difficulty estimates

- interval:

  interval range to search for `'IRF'` type

- ...:

  additional arguments to pass to
  [`uniroot`](https://rdrr.io/r/stats/uniroot.html)

## References

Ali, U. S., Chang, H.-H., & Anderson, C. J. (2015). *Location indices
for ordinal polytomous items based on item response theory* (Research
Report No. RR-15-20). Princeton, NJ: Educational Testing Service.
http://dx.doi.org/10.1002/ets2.12065

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

mod <- mirt(Science, 1)
coef(mod, simplify=TRUE, IRTpars = TRUE)$items
#>                a        b1         b2        b3
#> Comfort 1.041755 -4.669193 -2.5341299 1.4072541
#> Work    1.225962 -2.385068 -0.7350678 1.8488053
#> Future  2.293372 -2.282226 -0.9652918 0.8562529
#> Benefit 1.094915 -3.057698 -0.9056673 1.5419094

gen.difficulty(mod)
#>    Comfort       Work     Future    Benefit 
#> -2.3089094 -0.5741303 -0.9207845 -0.8530161 
gen.difficulty(mod, type = 'mean')
#>    Comfort       Work     Future    Benefit 
#> -1.9320231 -0.4237770 -0.7970883 -0.8071519 

# also works for dichotomous items (though this is unnecessary)
dat <- expand.table(LSAT7)
mod <- mirt(dat, 1)
coef(mod, simplify=TRUE, IRTpars = TRUE)$items
#>                a          b g u
#> Item.1 0.9879254 -1.8787456 0 1
#> Item.2 1.0808847 -0.7475160 0 1
#> Item.3 1.7058006 -1.0576962 0 1
#> Item.4 0.7651853 -0.6351358 0 1
#> Item.5 0.7357980 -2.5204102 0 1

gen.difficulty(mod)
#>     Item.1     Item.2     Item.3     Item.4     Item.5 
#> -1.8787448 -0.7475182 -1.0576961 -0.6351601 -2.5204127 
gen.difficulty(mod, type = 'mean')
#>     Item.1     Item.2     Item.3     Item.4     Item.5 
#> -1.8787456 -0.7475160 -1.0576962 -0.6351358 -2.5204102 

# }
```
