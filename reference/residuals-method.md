# Compute model residuals

Return model implied residuals for linear dependencies between items or
at the person level. If the latent trait density was approximated (e.g.,
Davidian curves, Empirical histograms, etc) then passing
`use_dentype_estimate = TRUE` will use the internally saved quadrature
and density components (where applicable).

## Usage

``` r
# S4 method for class 'SingleGroupClass'
residuals(
  object,
  type = "LD",
  p.adjust = "none",
  df.p = FALSE,
  approx.z = FALSE,
  full.scores = FALSE,
  QMC = FALSE,
  printvalue = NULL,
  tables = FALSE,
  verbose = TRUE,
  Theta = NULL,
  suppress = NA,
  theta_lim = c(-6, 6),
  quadpts = NULL,
  fold = TRUE,
  upper = TRUE,
  technical = list(),
  ...
)
```

## Arguments

- object:

  an object of class `SingleGroupClass` or `MultipleGroupClass`.
  Bifactor models are automatically detected and utilized for better
  accuracy

- type:

  type of residuals to be displayed. Can be either `'LD'` or `'LDG2'`
  for a local dependence matrix based on the X2 or G2 statistics (Chen &
  Thissen, 1997), `'Q3'` for the statistic proposed by Yen (1984),
  `'JSI'` for the jack-knife statistic proposed Edwards et al. (2018),
  `'exp'` for the expected values for the frequencies of every response
  pattern, and `'expfull'` for the expected values for every
  theoretically observable response pattern. For the 'LD' and 'LDG2'
  types, the upper diagonal elements represent the standardized
  residuals in the form of signed Cramers V coefficients

- p.adjust:

  method to use for adjusting all p-values (see
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for available
  options). Default is `'none'`

- df.p:

  logical; print the degrees of freedom and p-values?

- approx.z:

  logical; transform \\\chi^2(df)\\ information from LD tests into
  approximate z-ratios instead using the transformation \\z=\sqrt{2 \*
  \chi^2} - \sqrt{2 \* df - 1}\\?

- full.scores:

  logical; compute relevant statistics for each subject in the original
  data?

- QMC:

  logical; use quasi-Monte Carlo integration? If `quadpts` is omitted
  the default number of nodes is 5000

- printvalue:

  a numeric value to be specified when using the `res='exp'` option.
  Only prints patterns that have standardized residuals greater than
  `abs(printvalue)`. The default (NULL) prints all response patterns

- tables:

  logical; for LD type, return the observed, expected, and standardized
  residual tables for each item combination?

- verbose:

  logical; allow information to be printed to the console?

- Theta:

  a matrix of factor scores used for statistics that require empirical
  estimates (i.e., Q3). If supplied, arguments typically passed to
  [`fscores()`](https://philchalmers.github.io/mirt/reference/fscores.md)
  will be ignored and these values will be used instead

- suppress:

  a numeric value indicating which parameter local dependency
  combinations to flag as being too high (for LD, LDG2, and Q3 the
  standardize correlations are used; for JSI, the z-ratios are used).
  Absolute values for the standardized estimates greater than this value
  will be returned, while all values less than this value will be set to
  missing

- theta_lim:

  range for the integration grid

- quadpts:

  number of quadrature nodes to use. The default is extracted from model
  (if available) or generated automatically if not available

- fold:

  logical; apply the sum 'folding' described by Edwards et al. (2018)
  for the JSI statistic?

- upper:

  logical; which portion of the matrix (upper versus lower triangle)
  should the `suppress` argument be applied to?

- technical:

  list of technical arguments when models are re-estimated (see
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  details)

- ...:

  additional arguments to be passed to
  [`fscores()`](https://philchalmers.github.io/mirt/reference/fscores.md)

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Chen, W. H. & Thissen, D. (1997). Local dependence indices for item
pairs using item response theory. *Journal of Educational and Behavioral
Statistics, 22*, 265-289.

Edwards, M. C., Houts, C. R. & Cai, L. (2018). A Diagnostic Procedure to
Detect Departures From Local Independence in Item Response Theory
Models. *Psychological Methods, 23*, 138-149.

Yen, W. (1984). Effects of local item dependence on the fit and equating
performance of the three parameter logistic model. *Applied
Psychological Measurement, 8*, 125-145.

## Examples

``` r
# \donttest{

x <- mirt(Science, 1)
residuals(x)
#> LD matrix (lower triangle) and standardized residual correlations (upper triangle)
#> 
#> Upper triangle summary:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.147  -0.136  -0.111  -0.045   0.041   0.152 
#> 
#>         Comfort   Work Future Benefit
#> Comfort         -0.147 -0.101   0.152
#> Work     25.512         0.088  -0.141
#> Future   12.002  9.208         -0.122
#> Benefit  27.321 23.235 17.461        
residuals(x, tables = TRUE)
#> $Comfort_Work
#> $Comfort_Work$Obs
#>    
#>       1   2   3   4
#>   1   2   0   1   2
#>   2   2  11  15   4
#>   3  24  71 148  23
#>   4   5  16  42  26
#> 
#> $Comfort_Work$Exp
#>           [,1]      [,2]       [,3]      [,4]
#> [1,]  1.046123  1.787854   1.922925  0.262908
#> [2,]  5.722535 11.386923  13.625623  1.972172
#> [3,] 23.322258 69.348472 140.029013 31.718143
#> [4,]  3.487224 14.109434  50.763258 21.495135
#> 
#> $Comfort_Work$std_res
#>    
#>              1          2          3          4
#>   1  0.9326115 -1.3371066 -0.6655567  3.3878247
#>   2 -1.5561252 -0.1146625  0.3723298  1.4439717
#>   3  0.1403392  0.1983204  0.6736015 -1.5479971
#>   4  0.8100928  0.5033119 -1.2299596  0.9716541
#> 
#> 
#> $Comfort_Future
#> $Comfort_Future$Obs
#>    
#>       1   2   3   4
#>   1   2   2   1   0
#>   2   3  11  11   7
#>   3   8  51 155  52
#>   4   1   8  43  37
#> 
#> $Comfort_Future$Exp
#>           [,1]      [,2]       [,3]       [,4]
#> [1,] 0.7680477  1.812936   2.091572  0.3472546
#> [2,] 3.6749826 11.241583  15.084907  2.7057791
#> [3,] 9.1976778 52.384929 148.579175 54.2561060
#> [4,] 0.7527106  6.828614  42.541147 39.7325783
#> 
#> $Comfort_Future$std_res
#>    
#>               1           2           3           4
#>   1  1.40572316  0.13893088 -0.75477217 -0.58928310
#>   2 -0.35209910 -0.07205318 -1.05174604  2.61058722
#>   3 -0.39491253 -0.19134814  0.52675890 -0.30629167
#>   4  0.28503062  0.44826372  0.07035076 -0.43351011
#> 
#> 
#> $Comfort_Benefit
#> $Comfort_Benefit$Obs
#>    
#>       1   2   3   4
#>   1   4   0   1   0
#>   2   5  10  13   4
#>   3  11  81 133  41
#>   4   1   9  46  33
#> 
#> $Comfort_Benefit$Exp
#>            [,1]     [,2]       [,3]       [,4]
#> [1,]  0.6710457  1.95217   1.959779  0.4368153
#> [2,]  3.6241191 12.16983  13.682082  3.2312240
#> [3,] 14.8454892 71.28246 131.343241 46.9466978
#> [4,]  2.3173267 14.55596  44.803655 28.1781067
#> 
#> $Comfort_Benefit$std_res
#>    
#>              1          2          3          4
#>   1  4.0637950 -1.3972006 -0.6855953 -0.6609200
#>   2  0.7227359 -0.6219892 -0.1844000  0.4276774
#>   3 -0.9980547  1.1509727  0.1445624 -0.8679073
#>   4 -0.8653661 -1.4562597  0.1787310  0.9083677
#> 
#> 
#> $Work_Future
#> $Work_Future$Obs
#>    
#>       1   2   3   4
#>   1   7  10  14   2
#>   2   3  28  57  10
#>   3   3  31 122  50
#>   4   1   3  17  34
#> 
#> $Work_Future$Exp
#>           [,1]      [,2]      [,3]      [,4]
#> [1,] 4.6004355 12.483087  14.34426  2.150355
#> [2,] 5.8020446 27.786418  52.26417 10.780046
#> [3,] 3.6903551 28.936581 118.01439 55.699488
#> [4,] 0.3005835  3.061977  23.67397 28.411829
#> 
#> $Work_Future$std_res
#>    
#>               1           2           3           4
#>   1  1.11874979 -0.70279865 -0.09089724 -0.10253275
#>   2 -1.16328067  0.04051800  0.65507895 -0.23757996
#>   3 -0.35936723  0.38358699  0.36688237 -0.76367799
#>   4  1.27571400 -0.03541829 -1.37166696  1.04838334
#> 
#> 
#> $Work_Benefit
#> $Work_Benefit$Obs
#>    
#>       1   2   3   4
#>   1   4   8  12   9
#>   2   6  34  47  11
#>   3   8  52 111  35
#>   4   3   6  23  23
#> 
#> $Work_Benefit$Exp
#>          [,1]      [,2]      [,3]      [,4]
#> [1,] 4.233251 13.170723  13.31610  2.858063
#> [2,] 7.660483 32.196149  44.96112 11.814929
#> [3,] 8.417022 46.931424 106.55195 44.440425
#> [4,] 1.147224  7.662122  26.95958 19.679427
#> 
#> $Work_Benefit$std_res
#>    
#>              1          2          3          4
#>   1 -0.1133671 -1.4247756 -0.3606627  3.6330344
#>   2 -0.5999380  0.3179060  0.3040694 -0.2370851
#>   3 -0.1437407  0.7398678  0.4309125 -1.4161278
#>   4  1.7298113 -0.6004660 -0.7625934  0.7485258
#> 
#> 
#> $Future_Benefit
#> $Future_Benefit$Obs
#>    
#>       1   2   3   4
#>   1   5   1   6   2
#>   2   5  32  30   5
#>   3   8  53 118  31
#>   4   3  14  39  40
#> 
#> $Future_Benefit$Exp
#>          [,1]      [,2]       [,3]       [,4]
#> [1,] 2.960508  6.706374   4.143629  0.5829072
#> [2,] 7.760422 29.142832  29.887302  5.4775063
#> [3,] 9.227453 52.779344 109.948223 36.3417805
#> [4,] 1.509597 11.331867  47.809604 36.3906499
#> 
#> $Future_Benefit$std_res
#>    
#>               1           2           3           4
#>   1  1.18532911 -2.20351681  0.91195681  1.85608815
#>   2 -0.99090687  0.52926095  0.02061456 -0.20402700
#>   3 -0.40407692  0.03037266  0.76788760 -0.88610042
#>   4  1.21303423  0.79260501 -1.27408623  0.59832081
#> 
#> 
residuals(x, type = 'exp')
#>    Comfort Work Future Benefit freq    exp std.res
#> 1        1    1      1       1    2  0.124   5.324
#> 2        1    3      2       1    1  0.067   3.605
#> 3        1    4      2       3    1  0.019   7.046
#> 4        1    4      3       1    1  0.006  12.642
#> 5        2    1      1       1    1  0.460   0.796
#> 6        2    1      2       4    1  0.095   2.930
#> 7        2    2      1       1    1  0.351   1.095
#> 8        2    2      2       2    4  2.147   1.264
#> 9        2    2      2       3    2  1.616   0.302
#> 10       2    2      3       1    1  0.377   1.015
#> 11       2    2      3       2    1  1.716  -0.547
#> 12       2    2      3       3    1  2.228  -0.823
#> 13       2    2      4       3    1  0.251   1.497
#> 14       2    3      1       3    1  0.213   1.707
#> 15       2    3      2       2    2  1.545   0.366
#> 16       2    3      2       3    1  1.489  -0.401
#> 17       2    3      3       2    3  2.248   0.502
#> 18       2    3      3       3    3  3.910  -0.460
#> 19       2    3      3       4    2  1.035   0.948
#> 20       2    3      4       1    1  0.041   4.763
#> 21       2    3      4       3    2  0.866   1.218
#> 22       2    4      2       1    1  0.029   5.690
#> 23       2    4      4       3    2  0.259   3.418
#> 24       2    4      4       4    1  0.184   1.902
#> 25       3    1      1       1    1  0.644   0.444
#> 26       3    1      1       3    2  0.638   1.705
#> 27       3    1      2       2    2  3.923  -0.971
#> 28       3    1      2       3    4  3.077   0.526
#> 29       3    1      3       2    5  3.567   0.759
#> 30       3    1      3       3    5  5.097  -0.043
#> 31       3    1      3       4    3  1.160   1.709
#> 32       3    1      4       3    1  0.764   0.269
#> 33       3    1      4       4    1  0.320   1.203
#> 34       3    2      1       2    1  1.781  -0.585
#> 35       3    2      1       4    1  0.154   2.157
#> 36       3    2      2       1    1  2.250  -0.833
#> 37       3    2      2       2   10  8.471   0.525
#> 38       3    2      2       3    7  8.067  -0.376
#> 39       3    2      3       1    1  2.221  -0.819
#> 40       3    2      3       2   16 11.743   1.242
#> 41       3    2      3       3   22 19.584   0.546
#> 42       3    2      3       4    5  4.950   0.023
#> 43       3    2      4       1    1  0.193   1.835
#> 44       3    2      4       2    1  1.293  -0.258
#> 45       3    2      4       3    3  3.778  -0.400
#> 46       3    2      4       4    2  1.716   0.217
#> 47       3    3      1       3    2  0.971   1.044
#> 48       3    3      2       1    1  1.764  -0.575
#> 49       3    3      2       2   13  7.852   1.837
#> 50       3    3      2       3   10  9.780   0.070
#> 51       3    3      2       4    2  1.975   0.018
#> 52       3    3      3       1    5  3.363   0.893
#> 53       3    3      3       2   23 20.447   0.565
#> 54       3    3      3       3   52 45.243   1.005
#> 55       3    3      3       4    8 14.757  -1.759
#> 56       3    3      4       2    7  4.389   1.246
#> 57       3    3      4       3   13 16.979  -0.966
#> 58       3    3      4       4   12 10.323   0.522
#> 59       3    4      1       3    1  0.090   3.030
#> 60       3    4      2       3    1  1.098  -0.093
#> 61       3    4      3       2    2  3.043  -0.598
#> 62       3    4      3       3    4  8.562  -1.559
#> 63       3    4      3       4    4  3.637   0.191
#> 64       3    4      4       1    1  0.160   2.100
#> 65       3    4      4       2    1  1.287  -0.253
#> 66       3    4      4       3    6  6.466  -0.183
#> 67       3    4      4       4    3  5.659  -1.118
#> 68       4    1      1       4    1  0.006  12.587
#> 69       4    1      2       2    1  0.325   1.183
#> 70       4    1      2       4    2  0.057   8.158
#> 71       4    1      3       4    1  0.300   1.278
#> 72       4    2      2       1    1  0.193   1.836
#> 73       4    2      2       3    3  1.017   1.967
#> 74       4    2      3       3    8  4.462   1.675
#> 75       4    2      3       4    2  1.448   0.459
#> 76       4    2      4       2    1  0.433   0.862
#> 77       4    2      4       4    1  1.080  -0.077
#> 78       4    3      2       3    1  1.665  -0.516
#> 79       4    3      3       2    2  4.874  -1.302
#> 80       4    3      3       3   21 14.002   1.870
#> 81       4    3      3       4    3  5.965  -1.214
#> 82       4    3      4       2    2  2.097  -0.067
#> 83       4    3      4       3    5 10.419  -1.679
#> 84       4    3      4       4    8  8.830  -0.279
#> 85       4    4      3       2    1  0.966   0.034
#> 86       4    4      3       3    2  3.571  -0.831
#> 87       4    4      3       4    3  2.063   0.652
#> 88       4    4      4       2    2  0.906   1.149
#> 89       4    4      4       3    6  5.780   0.092
#> 90       4    4      4       4   12  7.470   1.657
residuals(x, suppress = .15)
#> LD matrix (lower triangle) and standardized residual correlations (upper triangle)
#> 
#> Upper triangle summary:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.147  -0.136  -0.111  -0.045   0.041   0.152 
#> 
#>         Comfort Work Future Benefit
#> Comfort                       0.152
#> Work                               
#> Future                             
#> Benefit  27.321                    
residuals(x, df.p = TRUE)
#> Degrees of freedom (lower triangle) and p-values:
#> 
#>         Comfort  Work Future Benefit
#> Comfort         0.002  0.213   0.001
#> Work          9        0.418   0.006
#> Future        9 9.000          0.042
#> Benefit       9 9.000  9.000        
#> 
#> LD matrix (lower triangle) and standardized residual correlations (upper triangle)
#> 
#> Upper triangle summary:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.147  -0.136  -0.111  -0.045   0.041   0.152 
#> 
#>         Comfort   Work Future Benefit
#> Comfort         -0.147 -0.101   0.152
#> Work     25.512         0.088  -0.141
#> Future   12.002  9.208         -0.122
#> Benefit  27.321 23.235 17.461        
residuals(x, df.p = TRUE, p.adjust = 'fdr') # apply FWE control
#> Degrees of freedom (lower triangle) and p-values:
#> 
#>         Comfort  Work Future Benefit
#> Comfort         0.007  0.256   0.007
#> Work          9        0.418   0.011
#> Future        9 9.000          0.063
#> Benefit       9 9.000  9.000        
#> 
#> LD matrix (lower triangle) and standardized residual correlations (upper triangle)
#> 
#> Upper triangle summary:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.147  -0.136  -0.111  -0.045   0.041   0.152 
#> 
#>         Comfort   Work Future Benefit
#> Comfort         -0.147 -0.101   0.152
#> Work     25.512         0.088  -0.141
#> Future   12.002  9.208         -0.122
#> Benefit  27.321 23.235 17.461        

# Pearson's X2 estimate for goodness-of-fit
full_table <- residuals(x, type = 'expfull')
head(full_table)
#>   Comfort Work Future Benefit freq   exp    res
#> 1       1    1      1       1    2 0.124  5.324
#> 2       1    1      1       2    0 0.160 -0.400
#> 3       1    1      1       3    0 0.054 -0.233
#> 4       1    1      1       4    0 0.005 -0.073
#> 5       1    1      2       1    0 0.092 -0.303
#> 6       1    1      2       2    0 0.219 -0.468
X2 <- with(full_table, sum((freq - exp)^2 / exp))
df <- nrow(full_table) - extract.mirt(x, 'nest') - 1
p <- pchisq(X2, df = df, lower.tail=FALSE)
data.frame(X2, df, p, row.names='Pearson-X2')
#>                  X2  df            p
#> Pearson-X2 689.3347 239 2.942933e-45

# above FOG test as a function
PearsonX2 <- function(x){
   full_table <- residuals(x, type = 'expfull')
   X2 <- with(full_table, sum((freq - exp)^2 / exp))
   df <- nrow(full_table) - extract.mirt(x, 'nest') - 1
   p <- pchisq(X2, df = df, lower.tail=FALSE)
   data.frame(X2, df, p, row.names='Pearson-X2')
}
PearsonX2(x)
#>                  X2  df            p
#> Pearson-X2 689.3347 239 2.942933e-45


# extract results manually
out <- residuals(x, df.p = TRUE, verbose=FALSE)
str(out)
#> List of 2
#>  $ df.p: 'mirt_matrix' num [1:4, 1:4] NA 9 9 9 0.00245 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "Comfort" "Work" "Future" "Benefit"
#>   .. ..$ : chr [1:4] "Comfort" "Work" "Future" "Benefit"
#>  $ LD  : 'mirt_matrix' num [1:4, 1:4] NA 25.512 12.002 27.321 -0.147 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "Comfort" "Work" "Future" "Benefit"
#>   .. ..$ : chr [1:4] "Comfort" "Work" "Future" "Benefit"
out$df.p[1,2]
#> [1] 0.002454207

# with and without supplied factor scores
Theta <- fscores(x)
residuals(x, type = 'Q3', Theta=Theta)
#> Q3 summary statistics:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.320  -0.249  -0.225  -0.190  -0.205   0.085 
#> 
#>         Comfort   Work Future Benefit
#> Comfort   1.000 -0.203 -0.252   0.085
#> Work     -0.203  1.000 -0.208  -0.242
#> Future   -0.252 -0.208  1.000  -0.320
#> Benefit   0.085 -0.242 -0.320   1.000
residuals(x, type = 'Q3', method = 'ML')
#> Q3 summary statistics:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.514  -0.426  -0.357  -0.311  -0.270   0.053 
#> 
#>         Comfort   Work Future Benefit
#> Comfort   1.000 -0.262 -0.419   0.053
#> Work     -0.262  1.000 -0.428  -0.295
#> Future   -0.419 -0.428  1.000  -0.514
#> Benefit   0.053 -0.295 -0.514   1.000

# Edwards et al. (2018) JSI statistic
N <- 250
a <- rnorm(10, 1.7, 0.3)
d <- rnorm(10)
dat <- simdata(a, d, N=250, itemtype = '2PL')

mod <- mirt(dat, 1)
residuals(mod, type = 'JSI')
#> JSI summary statistics:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.751  -0.239  -0.057  -0.016   0.243   0.754 
#> 
#>         Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> Item_1         -0.049 -0.222 -0.081 -0.751  0.108  0.355 -0.327 -0.018   0.703
#> Item_2  -0.049         0.245  0.754 -0.057  0.146 -0.531 -0.357 -0.339   0.370
#> Item_3  -0.222  0.245        -0.696  0.243 -0.059 -0.505  0.115  0.535  -0.027
#> Item_4  -0.081  0.754 -0.696         0.127  0.252 -0.088  0.318 -0.239  -0.249
#> Item_5  -0.751 -0.057  0.243  0.127         0.461 -0.156 -0.107  0.252  -0.249
#> Item_6   0.108  0.146 -0.059  0.252  0.461        -0.152 -0.356 -0.006  -0.310
#> Item_7   0.355 -0.531 -0.505 -0.088 -0.156 -0.152         0.234  0.188  -0.236
#> Item_8  -0.327 -0.357  0.115  0.318 -0.107 -0.356  0.234        -0.151   0.389
#> Item_9  -0.018 -0.339  0.535 -0.239  0.252 -0.006  0.188 -0.151         -0.215
#> Item_10  0.703  0.370 -0.027 -0.249 -0.249 -0.310 -0.236  0.389 -0.215        
residuals(mod, type = 'JSI', fold=FALSE) # unfolded
#> JSI summary statistics:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.480  -0.121  -0.029  -0.008   0.107   0.419 
#> 
#>         Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> Item_1         -0.058 -0.075 -0.038 -0.480  0.044  0.234 -0.168 -0.022   0.371
#> Item_2   0.009         0.166  0.419 -0.004  0.078 -0.250 -0.201 -0.196   0.197
#> Item_3  -0.148  0.079        -0.377  0.180 -0.054 -0.232  0.045  0.275  -0.017
#> Item_4  -0.043  0.335 -0.319         0.070  0.114 -0.035  0.152 -0.122  -0.115
#> Item_5  -0.272 -0.053  0.063  0.056         0.142 -0.089 -0.059  0.081  -0.093
#> Item_6   0.064  0.068 -0.005  0.139  0.319        -0.048 -0.195 -0.008  -0.175
#> Item_7   0.120 -0.281 -0.273 -0.052 -0.066 -0.104         0.109  0.099  -0.137
#> Item_8  -0.159 -0.155  0.071  0.166 -0.049 -0.161  0.125        -0.073   0.198
#> Item_9   0.004 -0.143  0.260 -0.117  0.171  0.003  0.089 -0.078         -0.098
#> Item_10  0.332  0.172 -0.010 -0.134 -0.155 -0.136 -0.099  0.191 -0.117        

# LD between items 1-2
aLD <- numeric(10)
aLD[1:2] <- rnorm(2, 2.55, 0.15)
a2 <- cbind(a, aLD)
dat <- simdata(a2, d, N=250, itemtype = '2PL')

mod <- mirt(dat, 1)

# JSI executed in parallel over multiple cores
if(interactive()) mirtCluster()
residuals(mod, type = 'JSI')
#> JSI summary statistics:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -1.142  -0.427  -0.071   0.032   0.437   2.777 
#> 
#>         Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> Item_1          2.777 -0.427 -0.734 -0.783 -0.090 -0.312 -0.741 -0.860   0.740
#> Item_2   2.777        -0.121 -0.338 -0.624 -0.882 -0.071 -0.600 -0.698   0.418
#> Item_3  -0.427 -0.121         1.148  0.192 -1.142  0.589 -0.247  0.166   0.145
#> Item_4  -0.734 -0.338  1.148         0.517 -0.257 -0.010  0.038  0.496  -0.328
#> Item_5  -0.783 -0.624  0.192  0.517         0.576  0.281  0.092  0.482  -0.226
#> Item_6  -0.090 -0.882 -1.142 -0.257  0.576        -0.158  1.277  0.437   0.982
#> Item_7  -0.312 -0.071  0.589 -0.010  0.281 -0.158         0.398  0.118  -0.637
#> Item_8  -0.741 -0.600 -0.247  0.038  0.092  1.277  0.398         0.780  -0.274
#> Item_9  -0.860 -0.698  0.166  0.496  0.482  0.437  0.118  0.780         -0.672
#> Item_10  0.740  0.418  0.145 -0.328 -0.226  0.982 -0.637 -0.274 -0.672        

# }
```
