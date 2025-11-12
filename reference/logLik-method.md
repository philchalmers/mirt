# Extract log-likelihood

Extract the observed-data log-likelihood.

## Usage

``` r
# S4 method for class 'SingleGroupClass'
logLik(object)
```

## Arguments

- object:

  an object of class `SingleGroupClass`, `MultipleGroupClass`, or
  `MixedClass`

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Examples

``` r
# \donttest{
x <- mirt(Science, 1)
logLik(x)
#> [1] -1608.87

# }
```
