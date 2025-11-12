# Create all possible combinations of vector input

This function constructs all possible k-way combinations of an input
vector. It is primarily useful when used in conjunction with the
[`mdirt`](https://philchalmers.github.io/mirt/reference/mdirt.md)
function, though users may have other uses for it as well. See
[`expand.grid`](https://rdrr.io/r/base/expand.grid.html) for more
flexible combination formats.

## Usage

``` r
thetaComb(theta, nfact, intercept = FALSE)
```

## Arguments

- theta:

  the vector from which all possible combinations should be obtained

- nfact:

  the number of observations (and therefore the number of columns to
  return in the matrix of combinations)

- intercept:

  logical; should a vector of 1's be appended to the first column of the
  result to include an intercept design component? Default is `FALSE`

## Value

a matrix with all possible combinations

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# all possible joint combinations for the vector -4 to 4
thetaComb(-4:4, 2)
#>       [,1] [,2]
#>  [1,]   -4   -4
#>  [2,]   -3   -4
#>  [3,]   -2   -4
#>  [4,]   -1   -4
#>  [5,]    0   -4
#>  [6,]    1   -4
#>  [7,]    2   -4
#>  [8,]    3   -4
#>  [9,]    4   -4
#> [10,]   -4   -3
#> [11,]   -3   -3
#> [12,]   -2   -3
#> [13,]   -1   -3
#> [14,]    0   -3
#> [15,]    1   -3
#> [16,]    2   -3
#> [17,]    3   -3
#> [18,]    4   -3
#> [19,]   -4   -2
#> [20,]   -3   -2
#> [21,]   -2   -2
#> [22,]   -1   -2
#> [23,]    0   -2
#> [24,]    1   -2
#> [25,]    2   -2
#> [26,]    3   -2
#> [27,]    4   -2
#> [28,]   -4   -1
#> [29,]   -3   -1
#> [30,]   -2   -1
#> [31,]   -1   -1
#> [32,]    0   -1
#> [33,]    1   -1
#> [34,]    2   -1
#> [35,]    3   -1
#> [36,]    4   -1
#> [37,]   -4    0
#> [38,]   -3    0
#> [39,]   -2    0
#> [40,]   -1    0
#> [41,]    0    0
#> [42,]    1    0
#> [43,]    2    0
#> [44,]    3    0
#> [45,]    4    0
#> [46,]   -4    1
#> [47,]   -3    1
#> [48,]   -2    1
#> [49,]   -1    1
#> [50,]    0    1
#> [51,]    1    1
#> [52,]    2    1
#> [53,]    3    1
#> [54,]    4    1
#> [55,]   -4    2
#> [56,]   -3    2
#> [57,]   -2    2
#> [58,]   -1    2
#> [59,]    0    2
#> [60,]    1    2
#> [61,]    2    2
#> [62,]    3    2
#> [63,]    4    2
#> [64,]   -4    3
#> [65,]   -3    3
#> [66,]   -2    3
#> [67,]   -1    3
#> [68,]    0    3
#> [69,]    1    3
#> [70,]    2    3
#> [71,]    3    3
#> [72,]    4    3
#> [73,]   -4    4
#> [74,]   -3    4
#> [75,]   -2    4
#> [76,]   -1    4
#> [77,]    0    4
#> [78,]    1    4
#> [79,]    2    4
#> [80,]    3    4
#> [81,]    4    4

# all possible binary combinations for four observations
thetaComb(c(0,1), 4)
#>       [,1] [,2] [,3] [,4]
#>  [1,]    0    0    0    0
#>  [2,]    1    0    0    0
#>  [3,]    0    1    0    0
#>  [4,]    1    1    0    0
#>  [5,]    0    0    1    0
#>  [6,]    1    0    1    0
#>  [7,]    0    1    1    0
#>  [8,]    1    1    1    0
#>  [9,]    0    0    0    1
#> [10,]    1    0    0    1
#> [11,]    0    1    0    1
#> [12,]    1    1    0    1
#> [13,]    0    0    1    1
#> [14,]    1    0    1    1
#> [15,]    0    1    1    1
#> [16,]    1    1    1    1

# all possible binary combinations for four observations (with intercept)
thetaComb(c(0,1), 4, intercept=TRUE)
#>       [,1] [,2] [,3] [,4] [,5]
#>  [1,]    1    0    0    0    0
#>  [2,]    1    1    0    0    0
#>  [3,]    1    0    1    0    0
#>  [4,]    1    1    1    0    0
#>  [5,]    1    0    0    1    0
#>  [6,]    1    1    0    1    0
#>  [7,]    1    0    1    1    0
#>  [8,]    1    1    1    1    0
#>  [9,]    1    0    0    0    1
#> [10,]    1    1    0    0    1
#> [11,]    1    0    1    0    1
#> [12,]    1    1    1    0    1
#> [13,]    1    0    0    1    1
#> [14,]    1    1    0    1    1
#> [15,]    1    0    1    1    1
#> [16,]    1    1    1    1    1
```
