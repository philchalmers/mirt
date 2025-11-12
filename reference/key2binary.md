# Score a test by converting response patterns to binary data

The `key2binary` function will convert response pattern data to a
dichotomous format, given a response key.

## Usage

``` r
key2binary(fulldata, key, score_missing = FALSE)
```

## Arguments

- fulldata:

  an object of class `data.frame`, `matrix`, or `table` with the
  response patterns

- key:

  a vector or matrix consisting of the 'correct' response to the items.
  Each value/row corresponds to each column in `fulldata`. If the input
  is a matrix, multiple scoring keys can be supplied for each item. NA
  values are used to indicate no scoring key (or in the case of a matrix
  input, no additional scoring keys)

- score_missing:

  logical; should missing data elements be returned as incorrect (i.e.,
  0)? If `FALSE`, all missing data terms will be kept as missing

## Value

Returns a numeric matrix with all the response patterns in dichotomous
format

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
data(SAT12)
head(SAT12)
#>   Item.1 Item.2 Item.3 Item.4 Item.5 Item.6 Item.7 Item.8 Item.9 Item.10
#> 1      1      4      5      2      3      1      2      1      3       1
#> 2      3      4      2      8      3      3      2      8      3       1
#> 3      1      4      5      4      3      2      2      3      3       2
#> 4      2      4      4      2      3      3      2      4      3       2
#> 5      2      4      5      2      3      2      2      1      1       2
#> 6      1      4      3      1      3      2      2      3      3       1
#>   Item.11 Item.12 Item.13 Item.14 Item.15 Item.16 Item.17 Item.18 Item.19
#> 1       2       4       2       1       5       3       4       4       1
#> 2       2       8       2       1       5       2       4       1       1
#> 3       2       1       3       1       5       5       4       1       3
#> 4       2       4       2       1       5       2       4       1       3
#> 5       2       4       2       1       5       4       4       5       1
#> 6       2       3       2       1       5       5       4       4       1
#>   Item.20 Item.21 Item.22 Item.23 Item.24 Item.25 Item.26 Item.27 Item.28
#> 1       4       3       3       4       1       3       5       1       3
#> 2       4       3       3       8       1       8       4       1       4
#> 3       4       3       3       1       1       3       4       1       3
#> 4       4       3       1       5       2       5       4       1       3
#> 5       4       3       3       3       1       1       5       1       3
#> 6       4       3       3       4       1       1       4       1       4
#>   Item.29 Item.30 Item.31 Item.32
#> 1       1       5       4       5
#> 2       5       8       4       8
#> 3       4       4       4       1
#> 4       4       2       4       2
#> 5       1       2       4       1
#> 6       2       3       4       3
key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)

dicho.SAT12 <- key2binary(SAT12, key)
head(dicho.SAT12)
#>      Item.1 Item.2 Item.3 Item.4 Item.5 Item.6 Item.7 Item.8 Item.9 Item.10
#> [1,]      1      1      1      1      1      1      1      1      1       1
#> [2,]      0      1      0      0      1      0      1      0      1       1
#> [3,]      1      1      1      0      1      0      1      0      1       0
#> [4,]      0      1      0      1      1      0      1      0      1       0
#> [5,]      0      1      1      1      1      0      1      1      0       0
#> [6,]      1      1      0      0      1      0      1      0      1       1
#>      Item.11 Item.12 Item.13 Item.14 Item.15 Item.16 Item.17 Item.18 Item.19
#> [1,]       1       1       1       1       1       1       1       1       1
#> [2,]       1       0       1       1       1       0       1       0       1
#> [3,]       1       0       0       1       1       0       1       0       0
#> [4,]       1       1       1       1       1       0       1       0       0
#> [5,]       1       1       1       1       1       0       1       0       1
#> [6,]       1       0       1       1       1       0       1       1       1
#>      Item.20 Item.21 Item.22 Item.23 Item.24 Item.25 Item.26 Item.27 Item.28
#> [1,]       1       1       1       1       1       1       1       1       1
#> [2,]       1       1       1       0       1       0       0       1       0
#> [3,]       1       1       1       0       1       1       0       1       1
#> [4,]       1       1       0       0       0       0       0       1       1
#> [5,]       1       1       1       0       1       0       1       1       1
#> [6,]       1       1       1       1       1       0       0       1       0
#>      Item.29 Item.30 Item.31 Item.32
#> [1,]       1       1       1       1
#> [2,]       0       0       1       0
#> [3,]       0       0       1       0
#> [4,]       0       0       1       0
#> [5,]       1       0       1       0
#> [6,]       0       0       1       0

# multiple scoring keys
key2 <- cbind(c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5),
              c(2,3,NA,1,rep(NA, 28)))
dicho.SAT12 <- key2binary(SAT12, key2)

# keys from raw character responses
resp <- as.data.frame(matrix(c(
  "B","B","D","D","E",
  "B","A","D","D","E",
  "B","A","D","C","E",
  "D","D","D","C","E",
  "B","C","A","D","A"), ncol=5, byrow=TRUE))

key <- c("B", "D", "D", "C", "E")

d01 <- key2binary(resp, key)
head(d01)
#>      V1 V2 V3 V4 V5
#> [1,]  1  0  1  0  1
#> [2,]  1  0  1  0  1
#> [3,]  1  0  1  1  1
#> [4,]  0  1  1  1  1
#> [5,]  1  0  0  0  0

# score/don't score missing values
resp[1,1] <- NA
d01NA <- key2binary(resp, key) # without scoring
d01NA
#>      V1 V2 V3 V4 V5
#> [1,] NA  0  1  0  1
#> [2,]  1  0  1  0  1
#> [3,]  1  0  1  1  1
#> [4,]  0  1  1  1  1
#> [5,]  1  0  0  0  0

d01 <- key2binary(resp, key, score_missing = TRUE) # with scoring
d01
#>      V1 V2 V3 V4 V5
#> [1,]  0  0  1  0  1
#> [2,]  1  0  1  0  1
#> [3,]  1  0  1  1  1
#> [4,]  0  1  1  1  1
#> [5,]  1  0  0  0  0

```
