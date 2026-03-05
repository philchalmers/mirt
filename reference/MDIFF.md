# Compute multidimensional difficulty index

Returns a matrix containing the MDIFF values (Reckase, 2009). Only
supported for items of class 'dich' and 'graded'.

## Usage

``` r
MDIFF(x, which.items = NULL, group = NULL)
```

## Arguments

- x:

  an object of class 'SingleGroupClass', or an object of class
  'MultipleGroupClass' if a suitable `group` input were supplied

- which.items:

  a vector indicating which items to select. If NULL is used (the
  default) then MDISC will be computed for all items

- group:

  group argument to pass to
  [`extract.group`](https://philchalmers.github.io/mirt/reference/extract.group.md)
  function. Required when the input object is a multiple-group model

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Reckase, M. D. (2009). Multidimensional Item Response Theory. Springer.

## See also

[`extract.group`](https://philchalmers.github.io/mirt/reference/extract.group.md),
[`MDISC`](https://philchalmers.github.io/mirt/reference/MDISC.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

mod <- mirt(Science, 2)
MDIFF(mod)
#>           MDIFF_1    MDIFF_2   MDIFF_3
#> Comfort -3.892829 -2.1408175 1.1974522
#> Work    -1.806316 -0.5623950 1.4163829
#> Future  -2.486824 -1.0433836 0.9256248
#> Benefit -2.316311 -0.6940864 1.1869606

mod <- mirt(expand.table(LSAT7), 2)
MDIFF(mod)
#>           MDIFF_1
#> Item.1 -1.2103317
#> Item.2 -0.7903550
#> Item.3 -0.8774211
#> Item.4 -0.6408189
#> Item.5 -2.4615861

# }
```
