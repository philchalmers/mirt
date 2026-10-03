# Generic item summary statistics

Function to compute generic item summary statistics that do not require
prior fitting of IRT models. Contains information about sample sizes
(`N`), number of observed categories (`K`), (standardized) coefficient
alpha (and alpha if an item is removed; (`alpha_if_rm`)), mean/SD and
frequency of total scores, reduced item-total correlations
(`cor_if_rm`), average/sd of the correlation between items, squared
multiple correlation (`smc`), response frequencies, and conditional
mean/sd information given the unweighted sum scores. Summary information
involving the total scores only included for responses with no missing
data to ensure the metric is meaningful, however standardized statistics
(e.g., correlations) utilize all possible response information.

## Usage

``` r
itemstats(
  data,
  group = NULL,
  use_ts = TRUE,
  itemfreq = "proportions",
  ts.tables = FALSE
)
```

## Arguments

- data:

  An object of class `data.frame` or `matrix` with the response patterns

- group:

  optional grouping variable to condition on when computing summary
  information

- use_ts:

  logical; include information that is conditional on a meaningful total
  score?

- itemfreq:

  character vector indicting whether to include item response
  `"proportions"` or `"counts"` for each item. If set to `'none'` then
  this will be omitted

- ts.tables:

  logical; include mean/sd summary information pertaining to the
  unweighted total score?

## Value

Returns a list containing the summary statistics

## References

Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`empirical_plot`](https://philchalmers.github.io/mirt/reference/empirical_plot.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# dichotomous data example
LSAT7full <- expand.table(LSAT7)
head(LSAT7full)
#>   Item.1 Item.2 Item.3 Item.4 Item.5
#> 1      0      0      0      0      0
#> 2      0      0      0      0      0
#> 3      0      0      0      0      0
#> 4      0      0      0      0      0
#> 5      0      0      0      0      0
#> 6      0      0      0      0      0
itemstats(LSAT7full)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found
itemstats(LSAT7full, itemfreq='counts')
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

# behaviour with missing data
LSAT7full[1:5,1] <- NA
itemstats(LSAT7full)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

# data with no meaningful total score
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
itemstats(SAT12, use_ts=FALSE)
#> $overall
#>     N
#> 1 600
#> 
#> $itemstats
#>           N K  mean    sd
#> Item.1  600 6 2.497 1.188
#> Item.2  600 6 3.385 1.356
#> Item.3  600 6 3.212 1.534
#> Item.4  600 6 2.762 1.370
#> Item.5  600 6 2.868 0.911
#> Item.6  600 5 2.358 1.135
#> Item.7  600 6 2.422 0.908
#> Item.8  600 6 2.925 1.370
#> Item.9  600 5 2.907 0.567
#> Item.10 600 6 2.320 1.490
#> Item.11 600 5 2.017 0.199
#> Item.12 600 6 3.642 1.184
#> Item.13 600 5 2.317 0.956
#> Item.14 600 6 1.798 1.432
#> Item.15 600 6 4.535 1.087
#> Item.16 600 6 3.368 1.135
#> Item.17 600 5 3.968 0.343
#> Item.18 600 6 3.020 1.514
#> Item.19 600 5 1.900 1.053
#> Item.20 600 6 3.870 0.483
#> Item.21 600 6 2.937 0.554
#> Item.22 600 5 2.985 0.442
#> Item.23 600 6 2.755 1.437
#> Item.24 600 6 1.502 1.037
#> Item.25 600 6 2.740 1.380
#> Item.26 600 6 3.923 1.265
#> Item.27 600 6 1.240 0.766
#> Item.28 600 6 3.262 0.937
#> Item.29 600 6 2.285 1.306
#> Item.30 600 6 3.703 1.553
#> Item.31 600 6 3.788 0.899
#> Item.32 600 6 3.023 1.303
#> 
#> $proportions
#>             1     2     3     4     5     8
#> Item.1  0.283 0.203 0.267 0.232 0.013 0.002
#> Item.2  0.212 0.022 0.070 0.568 0.127 0.002
#> Item.3  0.165 0.183 0.260 0.098 0.280 0.013
#> Item.4  0.165 0.378 0.148 0.172 0.128 0.008
#> Item.5  0.093 0.143 0.620 0.093 0.048 0.002
#> Item.6  0.160 0.582 0.107 0.043 0.108    NA
#> Item.7  0.025 0.760 0.007 0.190 0.017 0.002
#> Item.8  0.202 0.205 0.207 0.250 0.133 0.003
#> Item.9  0.065 0.010 0.885 0.033 0.007    NA
#> Item.10 0.422 0.215 0.165 0.028 0.167 0.003
#> Item.11 0.003 0.983 0.008 0.003 0.002    NA
#> Item.12 0.072 0.082 0.218 0.415 0.205 0.008
#> Item.13 0.110 0.662 0.070 0.118 0.040    NA
#> Item.14 0.723 0.027 0.108 0.022 0.117 0.003
#> Item.15 0.035 0.062 0.060 0.025 0.817 0.002
#> Item.16 0.070 0.105 0.413 0.215 0.195 0.002
#> Item.17 0.008 0.005 0.010 0.963 0.013    NA
#> Item.18 0.303 0.033 0.165 0.352 0.142 0.005
#> Item.19 0.548 0.053 0.358 0.030 0.010    NA
#> Item.20 0.012 0.002 0.105 0.873 0.007 0.002
#> Item.21 0.050 0.008 0.915 0.013 0.012 0.002
#> Item.22 0.028 0.005 0.935 0.017 0.015    NA
#> Item.23 0.290 0.177 0.128 0.313 0.087 0.005
#> Item.24 0.728 0.162 0.042 0.022 0.045 0.002
#> Item.25 0.240 0.170 0.375 0.065 0.142 0.008
#> Item.26 0.020 0.227 0.030 0.262 0.460 0.002
#> Item.27 0.862 0.093 0.012 0.020 0.010 0.003
#> Item.28 0.082 0.010 0.530 0.337 0.037 0.005
#> Item.29 0.340 0.295 0.205 0.085 0.067 0.008
#> Item.30 0.150 0.110 0.107 0.183 0.440 0.010
#> Item.31 0.075 0.020 0.012 0.833 0.058 0.002
#> Item.32 0.125 0.183 0.443 0.075 0.162 0.012
#> 

# extra total scores tables
dat <- key2binary(SAT12,
                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,
                           5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
itemstats(dat, ts.tables=TRUE)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

# grouping information
group <- gl(2, 300, labels=c('G1', 'G2'))
itemstats(dat, group=group)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found


#####
# polytomous data example
itemstats(Science)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

# polytomous data with missing
newScience <- Science
newScience[1:5,1] <- NA
itemstats(newScience)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

# unequal categories
newScience[,1] <- ifelse(Science[,1] == 1, NA, Science[,1])
itemstats(newScience)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

merged <- data.frame(LSAT7full[1:392,], Science)
itemstats(merged)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found
```
