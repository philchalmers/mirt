# Expand summary table of patterns and frequencies

The `expand.table` function expands a summary table of unique response
patterns to a full sized data-set. By default the response frequencies
are assumed to be on rightmost column of the input data, though this can
be modified.

## Usage

``` r
expand.table(tabdata, freq = colnames(tabdata)[ncol(tabdata)], sample = FALSE)
```

## Arguments

- tabdata:

  An object of class `data.frame` or `matrix` with the unique response
  patterns and the number of frequencies in the rightmost column (though
  see `freq` for details on how to omit this column)

- freq:

  either a character vector specifying the column in `tabdata` to be
  used as the frequency count indicator for each response pattern
  (defaults to the right-most column) or a integer vector of length
  `nrow(tabdata)` specifying the frequency counts. When using the latter
  approach the `tabdata` input should not include any information
  regarding the counts, and instead should only include the unique
  response patterns themselves

- sample:

  logical; randomly switch the rows in the expanded table? This does not
  change the expanded data, only the row locations

## Value

Returns a numeric matrix with all the response patterns.

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
data(LSAT7)
head(LSAT7) # frequency in right-most column
#>   Item.1 Item.2 Item.3 Item.4 Item.5 freq
#> 1      0      0      0      0      0   12
#> 2      0      0      0      0      1   19
#> 3      0      0      0      1      0    1
#> 4      0      0      0      1      1    7
#> 5      0      0      1      0      0    3
#> 6      0      0      1      0      1   19
LSAT7full <- expand.table(LSAT7)
head(LSAT7full)
#>   Item.1 Item.2 Item.3 Item.4 Item.5
#> 1      0      0      0      0      0
#> 2      0      0      0      0      0
#> 3      0      0      0      0      0
#> 4      0      0      0      0      0
#> 5      0      0      0      0      0
#> 6      0      0      0      0      0
dim(LSAT7full)
#> [1] 1000    5

# randomly switch rows in the expanded response table
LSAT7samp <- expand.table(LSAT7, sample = TRUE)
head(LSAT7samp)
#>   Item.1 Item.2 Item.3 Item.4 Item.5
#> 1      1      1      1      1      1
#> 2      1      0      1      1      1
#> 3      1      1      0      1      1
#> 4      1      0      1      1      1
#> 5      1      1      1      1      1
#> 6      1      1      1      1      1
colMeans(LSAT7full)
#> Item.1 Item.2 Item.3 Item.4 Item.5 
#>  0.828  0.658  0.772  0.606  0.843 
colMeans(LSAT7samp) #equal
#> Item.1 Item.2 Item.3 Item.4 Item.5 
#>  0.828  0.658  0.772  0.606  0.843 

#--------

# \donttest{
# Generate data from separate response pattern matrix and freq vector
# The following uses Table 2.1 from de Ayala (2009)
f <- c(691,2280,242,235,158,184,1685,1053,134,462,92,65,571,79,87,41,1682,702,
       370,63,626,412,166,52,28,15,2095,1219,500,187,40,3385)

pat <- matrix(c(
   0, 0, 0, 0, 0,
   1, 0, 0, 0, 0,
   0, 1, 0, 0, 0,
   0, 0, 1, 0, 0,
   0, 0, 0, 1, 0,
   0, 0, 0, 0, 1,
   1, 1, 0, 0, 0,
   1, 0, 1, 0, 0,
   0, 1, 1, 0, 0,
   1, 0, 0, 1, 0,
   0, 1, 0, 1, 0,
   0, 0, 1, 1, 0,
   1, 0, 0, 0, 1,
   0, 1, 0, 0, 1,
   0, 0, 1, 0, 1,
   0, 0, 0, 1, 1,
   1, 1, 1, 0, 0,
   1, 1, 0, 1, 0,
   1, 0, 1, 1, 0,
   0, 1, 1, 1, 0,
   1, 1, 0, 0, 1,
   1, 0, 1, 0, 1,
   1, 0, 0, 1, 1,
   0, 1, 1, 0, 1,
   0, 1, 0, 1, 1,
   0, 0, 1, 1, 1,
   1, 1, 1, 1, 0,
   1, 1, 1, 0, 1,
   1, 1, 0, 1, 1,
   1, 0, 1, 1, 1,
   0, 1, 1, 1, 1,
   1, 1, 1, 1, 1), ncol=5, byrow=TRUE)

colnames(pat) <- paste0('Item.', 1:5)
head(pat)
#>      Item.1 Item.2 Item.3 Item.4 Item.5
#> [1,]      0      0      0      0      0
#> [2,]      1      0      0      0      0
#> [3,]      0      1      0      0      0
#> [4,]      0      0      1      0      0
#> [5,]      0      0      0      1      0
#> [6,]      0      0      0      0      1

table2.1 <- expand.table(pat, freq = f)
dim(table2.1)
#> [1] 19601     5

# }

```
