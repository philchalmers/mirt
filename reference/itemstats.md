# Generic item summary statistics

Function to compute generic item summary statistics that do not require
prior fitting of IRT models. Contains information about sample sizes
(`N`), number of observed categories (`K`), coefficient alpha (and alpha
if an item is deleted), mean/SD and frequency of total scores, reduced
item-total correlations, average/sd of the correlation between items,
response frequencies, and conditional mean/sd information given the
unweighted sum scores. Summary information involving the total scores
only included for responses with no missing data to ensure the metric is
meaningful, however standardized statistics (e.g., correlations) utilize
all possible response information.

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

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
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
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  1000            3.707          1.199 0.143 0.052 0.453     0.886
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1 1000 2 0.828 0.378   0.530         0.246       0.396
#> Item.2 1000 2 0.658 0.475   0.600         0.247       0.394
#> Item.3 1000 2 0.772 0.420   0.611         0.313       0.345
#> Item.4 1000 2 0.606 0.489   0.592         0.223       0.415
#> Item.5 1000 2 0.843 0.364   0.461         0.175       0.438
#> 
#> $proportions
#>            0     1
#> Item.1 0.172 0.828
#> Item.2 0.342 0.658
#> Item.3 0.228 0.772
#> Item.4 0.394 0.606
#> Item.5 0.157 0.843
#> 
itemstats(LSAT7full, itemfreq='counts')
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  1000            3.707          1.199 0.143 0.052 0.453     0.886
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1 1000 2 0.828 0.378   0.530         0.246       0.396
#> Item.2 1000 2 0.658 0.475   0.600         0.247       0.394
#> Item.3 1000 2 0.772 0.420   0.611         0.313       0.345
#> Item.4 1000 2 0.606 0.489   0.592         0.223       0.415
#> Item.5 1000 2 0.843 0.364   0.461         0.175       0.438
#> 
#> $counts
#>          0   1
#> Item.1 172 828
#> Item.2 342 658
#> Item.3 228 772
#> Item.4 394 606
#> Item.5 157 843
#> 

# behaviour with missing data
LSAT7full[1:5,1] <- NA
itemstats(LSAT7full)
#> $overall
#>  N.complete    N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>         995 1000            3.726          1.172 0.137 0.052 0.426     0.888
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1  995 2 0.832 0.374   0.515         0.222       0.396
#> Item.2 1000 2 0.658 0.475   0.600         0.247       0.364
#> Item.3 1000 2 0.772 0.420   0.611         0.313       0.316
#> Item.4 1000 2 0.606 0.489   0.592         0.223       0.384
#> Item.5 1000 2 0.843 0.364   0.461         0.175       0.418
#> 
#> $proportions
#>            0     1  <NA>
#> Item.1 0.167 0.828 0.005
#> Item.2 0.342 0.658    NA
#> Item.3 0.228 0.772    NA
#> Item.4 0.394 0.606    NA
#> Item.5 0.157 0.843    NA
#> 

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
#> $overall
#>    N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  600           18.202          5.054 0.108 0.075 0.798     2.272
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1  600 2 0.283 0.451   0.380         0.300       0.793
#> Item.2  600 2 0.568 0.496   0.539         0.464       0.785
#> Item.3  600 2 0.280 0.449   0.446         0.371       0.789
#> Item.4  600 2 0.378 0.485   0.325         0.235       0.796
#> Item.5  600 2 0.620 0.486   0.424         0.340       0.791
#> Item.6  600 2 0.160 0.367   0.414         0.351       0.791
#> Item.7  600 2 0.760 0.427   0.366         0.289       0.793
#> Item.8  600 2 0.202 0.402   0.307         0.233       0.795
#> Item.9  600 2 0.885 0.319   0.189         0.127       0.798
#> Item.10 600 2 0.422 0.494   0.465         0.383       0.789
#> Item.11 600 2 0.983 0.128   0.181         0.156       0.797
#> Item.12 600 2 0.415 0.493   0.173         0.076       0.803
#> Item.13 600 2 0.662 0.474   0.438         0.358       0.790
#> Item.14 600 2 0.723 0.448   0.411         0.333       0.791
#> Item.15 600 2 0.817 0.387   0.393         0.325       0.792
#> Item.16 600 2 0.413 0.493   0.367         0.278       0.794
#> Item.17 600 2 0.963 0.188   0.238         0.202       0.796
#> Item.18 600 2 0.352 0.478   0.576         0.508       0.783
#> Item.19 600 2 0.548 0.498   0.401         0.314       0.792
#> Item.20 600 2 0.873 0.333   0.376         0.318       0.792
#> Item.21 600 2 0.915 0.279   0.190         0.136       0.798
#> Item.22 600 2 0.935 0.247   0.284         0.238       0.795
#> Item.23 600 2 0.313 0.464   0.338         0.253       0.795
#> Item.24 600 2 0.728 0.445   0.422         0.346       0.791
#> Item.25 600 2 0.375 0.485   0.383         0.297       0.793
#> Item.26 600 2 0.460 0.499   0.562         0.489       0.783
#> Item.27 600 2 0.862 0.346   0.425         0.367       0.791
#> Item.28 600 2 0.530 0.500   0.465         0.383       0.789
#> Item.29 600 2 0.340 0.474   0.407         0.324       0.791
#> Item.30 600 2 0.440 0.497   0.255         0.159       0.799
#> Item.31 600 2 0.833 0.373   0.479         0.419       0.788
#> Item.32 600 2 0.162 0.368   0.110         0.037       0.802
#> 
#> $proportions
#>             0     1
#> Item.1  0.717 0.283
#> Item.2  0.432 0.568
#> Item.3  0.720 0.280
#> Item.4  0.622 0.378
#> Item.5  0.380 0.620
#> Item.6  0.840 0.160
#> Item.7  0.240 0.760
#> Item.8  0.798 0.202
#> Item.9  0.115 0.885
#> Item.10 0.578 0.422
#> Item.11 0.017 0.983
#> Item.12 0.585 0.415
#> Item.13 0.338 0.662
#> Item.14 0.277 0.723
#> Item.15 0.183 0.817
#> Item.16 0.587 0.413
#> Item.17 0.037 0.963
#> Item.18 0.648 0.352
#> Item.19 0.452 0.548
#> Item.20 0.127 0.873
#> Item.21 0.085 0.915
#> Item.22 0.065 0.935
#> Item.23 0.687 0.313
#> Item.24 0.272 0.728
#> Item.25 0.625 0.375
#> Item.26 0.540 0.460
#> Item.27 0.138 0.862
#> Item.28 0.470 0.530
#> Item.29 0.660 0.340
#> Item.30 0.560 0.440
#> Item.31 0.167 0.833
#> Item.32 0.838 0.162
#> 
#> $total.score_frequency
#>      4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
#> Freq 1 2 2 2 5 7 14 14 17 35 51 45 41 46 50 44 44 20 35 31 18 18 18 19  7  7  3
#>      31 32
#> Freq  1  3
#> 
#> $total.score_means
#>                0        1
#> Item.1  16.99535 21.25294
#> Item.2  15.07722 20.57478
#> Item.3  16.79630 21.81548
#> Item.4  16.92225 20.30396
#> Item.5  15.46930 19.87634
#> Item.6  17.28968 22.98958
#> Item.7  14.91667 19.23904
#> Item.8  17.42171 21.28926
#> Item.9  15.55072 18.54614
#> Item.10 16.19597 20.95257
#> Item.11 11.20000 18.32034
#> Item.12 17.46724 19.23695
#> Item.13 15.10837 19.78338
#> Item.14 14.84940 19.48387
#> Item.15 14.01818 19.14082
#> Item.16 16.64773 20.40726
#> Item.17 12.04545 18.43599
#> Item.18 16.05913 22.15166
#> Item.19 15.97048 20.03951
#> Item.20 13.21053 18.92557
#> Item.21 15.05882 18.49362
#> Item.22 12.76923 18.57932
#> Item.23 17.04854 20.72872
#> Item.24 14.71166 19.50343
#> Item.25 16.70400 20.69778
#> Item.26 15.58025 21.27899
#> Item.27 12.84337 19.06190
#> Item.28 15.70567 20.41509
#> Item.29 16.72727 21.06373
#> Item.30 17.06250 19.65152
#> Item.31 12.79000 19.28400
#> Item.32 17.95825 19.46392
#> 
#> $total.score_sds
#>                0        1
#> Item.1  4.495009 5.115311
#> Item.2  3.791287 4.583007
#> Item.3  4.322840 5.013323
#> Item.4  4.333771 5.444014
#> Item.5  4.280262 4.756698
#> Item.6  4.448894 5.353788
#> Item.7  4.313744 4.825059
#> Item.8  4.575452 5.661917
#> Item.9  5.007454 4.961288
#> Item.10 4.278821 4.736481
#> Item.11 4.184628 4.985966
#> Item.12 4.861326 5.147431
#> Item.13 4.274965 4.679466
#> Item.14 4.008502 4.822098
#> Item.15 4.212219 4.744448
#> Item.16 4.361290 5.155818
#> Item.17 4.613410 4.923352
#> Item.18 3.955174 4.446051
#> Item.19 4.349442 4.854746
#> Item.20 3.714174 4.809188
#> Item.21 5.092786 4.954396
#> Item.22 3.923355 4.906759
#> Item.23 4.505479 5.276912
#> Item.24 3.896380 4.816172
#> Item.25 4.379749 5.124107
#> Item.26 3.761839 4.627003
#> Item.27 3.775666 4.692898
#> Item.28 4.340725 4.593648
#> Item.29 4.241013 5.281332
#> Item.30 4.817161 4.984363
#> Item.31 3.364386 4.622779
#> Item.32 4.901850 5.638524
#> 

# grouping information
group <- gl(2, 300, labels=c('G1', 'G2'))
itemstats(dat, group=group)
#> $G1
#> $G1$overall
#>    N mean_total.score sd_total.score ave.r sd.r alpha SEM.alpha
#>  300           17.987          5.051 0.107 0.09 0.796     2.282
#> 
#> $G1$itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1  300 2 0.290 0.455   0.410         0.331       0.789
#> Item.2  300 2 0.573 0.495   0.534         0.458       0.783
#> Item.3  300 2 0.290 0.455   0.458         0.382       0.787
#> Item.4  300 2 0.353 0.479   0.342         0.255       0.793
#> Item.5  300 2 0.630 0.484   0.480         0.401       0.786
#> Item.6  300 2 0.130 0.337   0.398         0.340       0.789
#> Item.7  300 2 0.727 0.446   0.383         0.303       0.790
#> Item.8  300 2 0.213 0.410   0.352         0.277       0.791
#> Item.9  300 2 0.897 0.305   0.203         0.144       0.796
#> Item.10 300 2 0.403 0.491   0.433         0.349       0.788
#> Item.11 300 2 0.993 0.082   0.154         0.138       0.796
#> Item.12 300 2 0.413 0.493   0.123         0.026       0.803
#> Item.13 300 2 0.647 0.479   0.441         0.359       0.788
#> Item.14 300 2 0.727 0.446   0.394         0.316       0.790
#> Item.15 300 2 0.793 0.406   0.407         0.337       0.789
#> Item.16 300 2 0.377 0.485   0.338         0.249       0.793
#> Item.17 300 2 0.957 0.204   0.220         0.181       0.795
#> Item.18 300 2 0.360 0.481   0.567         0.497       0.781
#> Item.19 300 2 0.537 0.499   0.432         0.347       0.788
#> Item.20 300 2 0.863 0.344   0.355         0.293       0.791
#> Item.21 300 2 0.910 0.287   0.202         0.147       0.795
#> Item.22 300 2 0.927 0.261   0.260         0.211       0.794
#> Item.23 300 2 0.260 0.439   0.323         0.242       0.793
#> Item.24 300 2 0.707 0.456   0.421         0.342       0.789
#> Item.25 300 2 0.370 0.484   0.406         0.321       0.790
#> Item.26 300 2 0.473 0.500   0.539         0.463       0.783
#> Item.27 300 2 0.857 0.351   0.474         0.418       0.787
#> Item.28 300 2 0.537 0.499   0.435         0.350       0.788
#> Item.29 300 2 0.343 0.476   0.378         0.293       0.791
#> Item.30 300 2 0.440 0.497   0.245         0.149       0.798
#> Item.31 300 2 0.810 0.393   0.518         0.457       0.784
#> Item.32 300 2 0.180 0.385   0.036        -0.041       0.803
#> 
#> $G1$proportions
#>             0     1
#> Item.1  0.710 0.290
#> Item.2  0.427 0.573
#> Item.3  0.710 0.290
#> Item.4  0.647 0.353
#> Item.5  0.370 0.630
#> Item.6  0.870 0.130
#> Item.7  0.273 0.727
#> Item.8  0.787 0.213
#> Item.9  0.103 0.897
#> Item.10 0.597 0.403
#> Item.11 0.007 0.993
#> Item.12 0.587 0.413
#> Item.13 0.353 0.647
#> Item.14 0.273 0.727
#> Item.15 0.207 0.793
#> Item.16 0.623 0.377
#> Item.17 0.043 0.957
#> Item.18 0.640 0.360
#> Item.19 0.463 0.537
#> Item.20 0.137 0.863
#> Item.21 0.090 0.910
#> Item.22 0.073 0.927
#> Item.23 0.740 0.260
#> Item.24 0.293 0.707
#> Item.25 0.630 0.370
#> Item.26 0.527 0.473
#> Item.27 0.143 0.857
#> Item.28 0.463 0.537
#> Item.29 0.657 0.343
#> Item.30 0.560 0.440
#> Item.31 0.190 0.810
#> Item.32 0.820 0.180
#> 
#> 
#> $G2
#> $G2$overall
#>    N mean_total.score sd_total.score ave.r sd.r alpha SEM.alpha
#>  300           18.417          5.056  0.11 0.08   0.8     2.262
#> 
#> $G2$itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1  300 2 0.277 0.448   0.352         0.271       0.796
#> Item.2  300 2 0.563 0.497   0.547         0.472       0.787
#> Item.3  300 2 0.270 0.445   0.438         0.363       0.792
#> Item.4  300 2 0.403 0.491   0.305         0.213       0.799
#> Item.5  300 2 0.610 0.489   0.371         0.283       0.796
#> Item.6  300 2 0.190 0.393   0.426         0.360       0.792
#> Item.7  300 2 0.793 0.406   0.344         0.270       0.796
#> Item.8  300 2 0.190 0.393   0.265         0.190       0.799
#> Item.9  300 2 0.873 0.333   0.180         0.116       0.801
#> Item.10 300 2 0.440 0.497   0.495         0.415       0.789
#> Item.11 300 2 0.973 0.161   0.215         0.184       0.799
#> Item.12 300 2 0.417 0.494   0.222         0.127       0.803
#> Item.13 300 2 0.677 0.469   0.434         0.354       0.792
#> Item.14 300 2 0.720 0.450   0.428         0.351       0.792
#> Item.15 300 2 0.840 0.367   0.375         0.310       0.794
#> Item.16 300 2 0.450 0.498   0.391         0.303       0.795
#> Item.17 300 2 0.970 0.171   0.258         0.226       0.798
#> Item.18 300 2 0.343 0.476   0.588         0.522       0.784
#> Item.19 300 2 0.560 0.497   0.369         0.279       0.796
#> Item.20 300 2 0.883 0.322   0.398         0.343       0.794
#> Item.21 300 2 0.920 0.272   0.175         0.123       0.800
#> Item.22 300 2 0.943 0.232   0.309         0.266       0.797
#> Item.23 300 2 0.367 0.483   0.348         0.260       0.797
#> Item.24 300 2 0.750 0.434   0.421         0.347       0.793
#> Item.25 300 2 0.380 0.486   0.360         0.272       0.796
#> Item.26 300 2 0.447 0.498   0.590         0.520       0.784
#> Item.27 300 2 0.867 0.341   0.374         0.314       0.794
#> Item.28 300 2 0.523 0.500   0.498         0.418       0.789
#> Item.29 300 2 0.337 0.473   0.437         0.357       0.792
#> Item.30 300 2 0.440 0.497   0.265         0.170       0.801
#> Item.31 300 2 0.857 0.351   0.435         0.376       0.792
#> Item.32 300 2 0.143 0.351   0.196         0.128       0.800
#> 
#> $G2$proportions
#>             0     1
#> Item.1  0.723 0.277
#> Item.2  0.437 0.563
#> Item.3  0.730 0.270
#> Item.4  0.597 0.403
#> Item.5  0.390 0.610
#> Item.6  0.810 0.190
#> Item.7  0.207 0.793
#> Item.8  0.810 0.190
#> Item.9  0.127 0.873
#> Item.10 0.560 0.440
#> Item.11 0.027 0.973
#> Item.12 0.583 0.417
#> Item.13 0.323 0.677
#> Item.14 0.280 0.720
#> Item.15 0.160 0.840
#> Item.16 0.550 0.450
#> Item.17 0.030 0.970
#> Item.18 0.657 0.343
#> Item.19 0.440 0.560
#> Item.20 0.117 0.883
#> Item.21 0.080 0.920
#> Item.22 0.057 0.943
#> Item.23 0.633 0.367
#> Item.24 0.250 0.750
#> Item.25 0.620 0.380
#> Item.26 0.553 0.447
#> Item.27 0.133 0.867
#> Item.28 0.477 0.523
#> Item.29 0.663 0.337
#> Item.30 0.560 0.440
#> Item.31 0.143 0.857
#> Item.32 0.857 0.143
#> 
#> 


#####
# polytomous data example
itemstats(Science)
#> $overall
#>    N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  392           11.668          2.003 0.275 0.098 0.598      1.27
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Comfort 392 4 3.120 0.588   0.596         0.352       0.552
#> Work    392 4 2.722 0.807   0.666         0.332       0.567
#> Future  392 4 2.990 0.757   0.748         0.488       0.437
#> Benefit 392 4 2.837 0.802   0.684         0.363       0.541
#> 
#> $proportions
#>             1     2     3     4
#> Comfort 0.013 0.082 0.679 0.227
#> Work    0.084 0.250 0.526 0.140
#> Future  0.036 0.184 0.536 0.245
#> Benefit 0.054 0.255 0.492 0.199
#> 

# polytomous data with missing
newScience <- Science
newScience[1:5,1] <- NA
itemstats(newScience)
#> $overall
#>  N.complete   N mean_total.score sd_total.score ave.r sd.r alpha SEM.alpha
#>         387 392           11.672          2.011 0.276  0.1 0.605     1.264
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Comfort 387 4 3.119 0.590   0.596         0.352       0.552
#> Work    392 4 2.722 0.807   0.646         0.311       0.576
#> Future  392 4 2.990 0.757   0.740         0.481       0.449
#> Benefit 392 4 2.837 0.802   0.684         0.370       0.538
#> 
#> $proportions
#>             1     2     3     4  <NA>
#> Comfort 0.013 0.082 0.668 0.224 0.013
#> Work    0.084 0.250 0.526 0.140    NA
#> Future  0.036 0.184 0.536 0.245    NA
#> Benefit 0.054 0.255 0.492 0.199    NA
#> 

# unequal categories
newScience[,1] <- ifelse(Science[,1] == 1, NA, Science[,1])
itemstats(newScience)
#> $overall
#>  N.complete   N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>         387 392           11.731          1.917  0.26 0.092 0.572     1.254
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Comfort 387 3 3.147 0.540   0.556         0.314       0.552
#> Work    392 4 2.722 0.807   0.656         0.325       0.517
#> Future  392 4 2.990 0.757   0.746         0.490       0.409
#> Benefit 392 4 2.837 0.802   0.684         0.370       0.525
#> 
#> $proportions
#>             1     2     3     4  <NA>
#> Comfort    NA 0.082 0.679 0.227 0.013
#> Work    0.084 0.250 0.526 0.140    NA
#> Future  0.036 0.184 0.536 0.245    NA
#> Benefit 0.054 0.255 0.492 0.199    NA
#> 

merged <- data.frame(LSAT7full[1:392,], Science)
itemstats(merged)
#> $overall
#>  N.complete   N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>         387 392           14.331          2.231 0.037 0.167 0.379     1.759
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1  387 2 0.568 0.496   0.193        -0.030       0.417
#> Item.2  392 2 0.232 0.423   0.041        -0.146       0.443
#> Item.3  392 2 0.605 0.490   0.231         0.014       0.405
#> Item.4  392 2 0.467 0.500   0.272         0.051       0.392
#> Item.5  392 2 0.760 0.428   0.201         0.011       0.402
#> Comfort 392 4 3.120 0.588   0.519         0.288       0.286
#> Work    392 4 2.722 0.807   0.608         0.299       0.251
#> Future  392 4 2.990 0.757   0.676         0.418       0.185
#> Benefit 392 4 2.837 0.802   0.600         0.291       0.261
#> 
#> $proportions
#>             0     1     2     3     4  <NA>
#> Item.1  0.426 0.561    NA    NA    NA 0.013
#> Item.2  0.768 0.232    NA    NA    NA    NA
#> Item.3  0.395 0.605    NA    NA    NA    NA
#> Item.4  0.533 0.467    NA    NA    NA    NA
#> Item.5  0.240 0.760    NA    NA    NA    NA
#> Comfort    NA 0.013 0.082 0.679 0.227    NA
#> Work       NA 0.084 0.250 0.526 0.140    NA
#> Future     NA 0.036 0.184 0.536 0.245    NA
#> Benefit    NA 0.054 0.255 0.492 0.199    NA
#> 
```
