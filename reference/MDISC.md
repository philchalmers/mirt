# Compute multidimensional discrimination index

Returns a vector containing the MDISC values for each item in the model
input object (Reckase, 2009).

## Usage

``` r
MDISC(x, group = NULL)
```

## Arguments

- x:

  an object of class 'SingleGroupClass', or an object of class
  'MultipleGroupClass' if a suitable `group` input were supplied

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

[`extract.group`](https://philchalmers.github.io/mirt/reference/extract.group.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

mod <- mirt(Science, 2)
MDISC(mod)
#>  Comfort     Work   Future  Benefit 
#> 1.338530 2.050473 1.875270 1.722043 

# }
```
