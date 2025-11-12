# Change polytomous items to dichotomous item format

Transforms a matrix of items into a new matrix where the select
polytomous items have been converted into comparable dichotomous items
with the same information.

## Usage

``` r
poly2dich(data, which.items = 1:ncol(data), sep = "_cat.")
```

## Arguments

- data:

  an object of class `data.frame` or `matrix`

- which.items:

  a vector indicating which items should be transformed into the
  dichotomous form. Default uses all input items

- sep:

  character vector pattern to append to each item name in `data`

## Value

Returns an integer matrix

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
data(Science)

head(Science)
#>   Comfort Work Future Benefit
#> 1       4    4      3       2
#> 2       3    3      3       3
#> 3       3    2      2       3
#> 4       3    2      2       3
#> 5       3    4      4       1
#> 6       4    4      3       3
newScience <- poly2dich(Science)
head(newScience)
#>      Comfort_cat.1 Comfort_cat.2 Comfort_cat.3 Comfort_cat.4 Work_cat.1
#> [1,]             0             0             0             1          0
#> [2,]             0             0             1             0          0
#> [3,]             0             0             1             0          0
#> [4,]             0             0             1             0          0
#> [5,]             0             0             1             0          0
#> [6,]             0             0             0             1          0
#>      Work_cat.2 Work_cat.3 Work_cat.4 Future_cat.1 Future_cat.2 Future_cat.3
#> [1,]          0          0          1            0            0            1
#> [2,]          0          1          0            0            0            1
#> [3,]          1          0          0            0            1            0
#> [4,]          1          0          0            0            1            0
#> [5,]          0          0          1            0            0            0
#> [6,]          0          0          1            0            0            1
#>      Future_cat.4 Benefit_cat.1 Benefit_cat.2 Benefit_cat.3 Benefit_cat.4
#> [1,]            0             0             1             0             0
#> [2,]            0             0             0             1             0
#> [3,]            0             0             0             1             0
#> [4,]            0             0             0             1             0
#> [5,]            1             1             0             0             0
#> [6,]            0             0             0             1             0

newScience2 <- poly2dich(Science, which.items = 2)
head(newScience2)
#>   Comfort Work_cat.1 Work_cat.2 Work_cat.3 Work_cat.4 Future Benefit
#> 1       4          0          0          0          1      3       2
#> 2       3          0          0          1          0      3       3
#> 3       3          0          1          0          0      2       3
#> 4       3          0          1          0          0      2       3
#> 5       3          0          0          0          1      4       1
#> 6       4          0          0          0          1      3       3

# }
```
