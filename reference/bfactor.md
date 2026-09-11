# Full-Information Item Bi-factor and Two-Tier Analysis

`bfactor` fits a confirmatory maximum likelihood
two-tier/bifactor/testlet model to dichotomous and polytomous data under
the item response theory paradigm. The IRT models are fit using a
dimensional reduction EM algorithm so that regardless of the number of
specific factors estimated the model only uses the number of factors in
the second-tier structure plus 1. For the bifactor model the maximum
number of dimensions is only 2 since the second-tier only consists of a
ubiquitous unidimensional factor. See
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
appropriate methods to be used on the objects returned from the
estimation.

## Usage

``` r
bfactor(
  data,
  model,
  model2 = paste0("G = 1-", ncol(data)),
  group = NULL,
  quadpts = NULL,
  invariance = "",
  ...
)
```

## Arguments

- data:

  a `matrix` or `data.frame` that consists of numerically ordered data,
  organized in the form of integers, with missing data coded as `NA`

- model:

  a numeric vector specifying which factor loads on which item. For
  example, if for a 4 item test with two specific factors, the first
  specific factor loads on the first two items and the second specific
  factor on the last two, then the vector is `c(1,1,2,2)`. For items
  that should only load on the second-tier factors (have no specific
  component) `NA` values may be used as place-holders. These numbers
  will be translated into a format suitable for
  [`mirt.model()`](https://philchalmers.github.io/mirt/reference/mirt.model.md),
  combined with the definition in `model2`, with the letter 'S' added to
  the respective factor number

  Alternatively, input can be specified using the
  [`mirt.model`](https://philchalmers.github.io/mirt/reference/mirt.model.md)
  syntax with the restriction that each item must load on exactly one
  specific factor (or no specific factors, if it is only predicted by
  the general factor specified in `model2`)

  IMPORTANT: Additional model information (e.g., keywords such as
  `CONSTRAIN`, `MEAN`, etc) must be specified in the `model2` input
  instead of this specification

- model2:

  a two-tier model specification object defined by
  [`mirt.model()`](https://philchalmers.github.io/mirt/reference/mirt.model.md)
  or a string to be passed to
  [`mirt.model`](https://philchalmers.github.io/mirt/reference/mirt.model.md).
  By default the model will fit a unidimensional model in the
  second-tier, and therefore be equivalent to the bifactor model

- group:

  a factor variable indicating group membership used for multiple group
  analyses

- quadpts:

  number of quadrature nodes to use after accounting for the reduced
  number of dimensions. Scheme is the same as the one used in
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
  however it is in regards to the reduced dimensions (e.g., a bifactor
  model has 2 dimensions to be integrated)

- invariance:

  see
  [`multipleGroup`](https://philchalmers.github.io/mirt/reference/multipleGroup.md)
  for details, however, the specific factor variances and means will be
  constrained according to the dimensional reduction algorithm

- ...:

  additional arguments to be passed to the estimation engine. See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  more details and examples

## Value

function returns an object of class `SingleGroupClass`
([SingleGroupClass-class](https://philchalmers.github.io/mirt/reference/SingleGroupClass-class.md))
or
`MultipleGroupClass`([MultipleGroupClass-class](https://philchalmers.github.io/mirt/reference/MultipleGroupClass-class.md)).

## Details

`bfactor` follows the item factor analysis strategy explicated by
Gibbons and Hedeker (1992), Gibbons et al. (2007), and Cai (2010).
Nested models may be compared via an approximate chi-squared difference
test or by a reduction in AIC or BIC (accessible via
[`anova`](https://rdrr.io/r/stats/anova.html)). See
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for more
details regarding the IRT estimation approach used in this package.

The two-tier model has a specific block diagonal covariance structure
between the primary and secondary latent traits. Namely, the secondary
latent traits are assumed to be orthogonal to all traits and have a
fixed variance of 1, while the primary traits can be organized to vary
and covary with other primary traits in the model.

\$\$\Sigma\_{two-tier} = \left(\begin{array}{cc} G & 0 \\ 0 & diag(S)
\end{array} \right)\$\$

The bifactor model is a special case of the two-tier model when \\G\\
above is a 1x1 matrix, and therefore only 1 primary factor is being
modeled. Evaluation of the numerical integrals for the two-tier model
requires only \\ncol(G) + 1\\ dimensions for integration since the \\S\\
second order (or 'specific') factors require only 1 integration grid due
to the dimension reduction technique.

Note: for multiple group two-tier analyses only the second-tier means
and variances should be freed since the specific factors are not treated
independently due to the dimension reduction technique.

## References

Cai, L. (2010). A two-tier full-information item factor analysis model
with applications. *Psychometrika, 75*, 581-612.

Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Bradlow, E.T., Wainer, H., & Wang, X. (1999). A Bayesian random effects
model for testlets. *Psychometrika, 64*, 153-168.

Gibbons, R. D., & Hedeker, D. R. (1992). Full-information Item Bi-Factor
Analysis. *Psychometrika, 57*, 423-436.

Gibbons, R. D., Darrell, R. B., Hedeker, D., Weiss, D. J., Segawa, E.,
Bhaumik, D. K., Kupfer, D. J., Frank, E., Grochocinski, V. J., & Stover,
A. (2007). Full-Information item bifactor analysis of graded response
data. *Applied Psychological Measurement, 31*, 4-19.

Wainer, H., Bradlow, E.T., & Wang, X. (2007). Testlet response theory
and its applications. New York, NY: Cambridge University Press.

## See also

[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# \donttest{

### load SAT12 and compute bifactor model with 3 specific factors
data(SAT12)
data <- key2binary(SAT12,
  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
mod1 <- bfactor(data, specific)
#> 
summary(mod1)
#>             G     S1     S2     S3    h2
#> Item.1  0.408         0.227        0.218
#> Item.2  0.620                0.339 0.499
#> Item.3  0.557        -0.074        0.316
#> Item.4  0.281                0.310 0.175
#> Item.5  0.478                0.254 0.293
#> Item.6  0.534         0.272        0.359
#> Item.7  0.474  0.421               0.402
#> Item.8  0.354         0.273        0.200
#> Item.9  0.218  0.532               0.331
#> Item.10 0.485  0.379               0.379
#> Item.11 0.644  0.332               0.525
#> Item.12 0.070                0.159 0.030
#> Item.13 0.522  0.274               0.348
#> Item.14 0.479                0.456 0.438
#> Item.15 0.599  0.240               0.416
#> Item.16 0.389         0.205        0.193
#> Item.17 0.664  0.118               0.454
#> Item.18 0.716  0.078               0.519
#> Item.19 0.452                0.020 0.205
#> Item.20 0.658                0.179 0.465
#> Item.21 0.281  0.345               0.198
#> Item.22 0.702 -0.028               0.494
#> Item.23 0.324                0.266 0.175
#> Item.24 0.585  0.103               0.353
#> Item.25 0.373                0.330 0.248
#> Item.26 0.643                0.212 0.459
#> Item.27 0.737  0.155               0.568
#> Item.28 0.526                0.076 0.282
#> Item.29 0.419         0.707        0.675
#> Item.30 0.245               -0.096 0.069
#> Item.31 0.833 -0.087               0.702
#> Item.32 0.078         0.016        0.006
#> 
#> SS loadings:  8.435 1.03 0.748 0.781 
#> Proportion Var:  0.264 0.032 0.023 0.024 
#> 
#> Factor correlations: 
#> 
#>    G S1 S2 S3
#> G  1         
#> S1 0  1      
#> S2 0  0  1   
#> S3 0  0  0  1
itemplot(mod1, 18, drop.zeros = TRUE) #drop the zero slopes to allow plotting


# complete factor score predictions (general + specific factors)
eaps <- fscores(mod1)
#> Warning: High-dimensional models factor scores should use quasi-Monte Carlo integration. Pass QMC=TRUE

# factor score predictions for general factors only (more accurate due
# to lower dimensional integration)
eaps_gen <- fscores(mod1, method = 'EAP_general')
head(cbind(eaps_gen, NA, eaps))
#>                G              G         S1         S2          S3
#> [1,]  2.52983691 NA  2.52984165  0.1930511  0.6978441  0.57291912
#> [2,] -0.06966393 NA -0.06876657  0.7643308 -0.6155095 -0.04731736
#> [3,]  0.07791027 NA  0.07794124 -0.2428963 -0.4726210  0.32041577
#> [4,] -0.62501109 NA -0.62339611  0.3677513 -0.4792017  0.86868085
#> [5,]  0.69852130 NA  0.69876943 -0.8731009  0.4935891  0.45513250
#> [6,]  0.44291463 NA  0.44307267  0.7012492 -0.4925127 -0.04974689

# similar EAP estimates, but with respect to sum-scores
fscores(mod1, method = 'EAPsum_2.0', full.scores=FALSE)
#>    Sum.Scores      G  SE_G observed expected std.res
#> 0           0 -3.153 0.553        0    0.014   0.120
#> 1           1 -2.974 0.549        0    0.085   0.291
#> 2           2 -2.787 0.539        0    0.282   0.531
#> 3           3 -2.593 0.524        0    0.699   0.836
#> 4           4 -2.400 0.507        1    1.424   0.355
#> 5           5 -2.212 0.491        2    2.520   0.328
#> 6           6 -2.029 0.478        2    4.017   1.006
#> 7           7 -1.851 0.468        2    5.918   1.611
#> 8           8 -1.679 0.460        5    8.223   1.124
#> 9           9 -1.509 0.456        7   10.929   1.189
#> 10         10 -1.343 0.454       14   14.040   0.011
#> 11         11 -1.179 0.453       14   17.545   0.846
#> 12         12 -1.017 0.454       17   21.403   0.952
#> 13         13 -0.857 0.456       35   25.520   1.877
#> 14         14 -0.699 0.458       51   29.738   3.899
#> 15         15 -0.542 0.462       45   33.838   1.919
#> 16         16 -0.385 0.466       41   37.559   0.561
#> 17         17 -0.228 0.471       46   40.635   0.842
#> 18         18 -0.070 0.477       50   42.820   1.097
#> 19         19  0.092 0.484       44   43.911   0.013
#> 20         20  0.257 0.492       44   43.760   0.036
#> 21         21  0.427 0.502       20   42.278   3.426
#> 22         22  0.601 0.512       35   39.445   0.708
#> 23         23  0.781 0.524       31   35.336   0.729
#> 24         24  0.966 0.536       18   30.151   2.213
#> 25         25  1.156 0.549       18   24.247   1.269
#> 26         26  1.350 0.564       18   18.123   0.029
#> 27         27  1.551 0.580       19   12.351   1.892
#> 28         28  1.758 0.597        7    7.456   0.167
#> 29         29  1.969 0.615        7    3.799   1.642
#> 30         30  2.180 0.632        3    1.500   1.225
#> 31         31  2.377 0.645        1    0.388   0.982
#> 32         32  2.530 0.650        3    0.044  14.137

# alternative model definition via ?mirt.model syntax
specific2 <- "S1 = 7,9,10,11,13,15,17,18,21,22,24,27,31
              S2 = 1,3,6,8,16,29,32
              S3 = 2,4,5,12,14,19,20,23,25,26,28,30"
mod2 <- bfactor(data, specific2)
#> 
anova(mod1, mod2) # same
#>          AIC    SABIC       HQ      BIC    logLik X2 df   p
#> mod1 19062.1 19179.44 19226.42 19484.21 -9435.052          
#> mod2 19062.1 19179.44 19226.42 19484.21 -9435.052  0  0 NaN

# also equivalent using item names instead (not run)
specific3 <- "S1 = Item.7, Item.9, Item.10, Item.11, Item.13, Item.15,
                Item.17, Item.18, Item.21, Item.22, Item.24, Item.27, Item.31
              S2 = Item.1, Item.3, Item.6, Item.8, Item.16, Item.29, Item.32
              S3 = Item.2, Item.4, Item.5, Item.12, Item.14, Item.19,
                Item.20, Item.23, Item.25, Item.26, Item.28, Item.30"
# mod3 <- bfactor(data, specific3)
# anova(mod1, mod2, mod3)  # all same

### Try with fixed guessing parameters added
guess <- rep(.1,32)
mod2 <- bfactor(data, specific, guess = guess)
#> 
coef(mod2)
#> $Item.1
#>        a1 a2    a3 a4      d   g u
#> par 1.225  0 0.624  0 -1.822 0.1 1
#> 
#> $Item.2
#>        a1 a2 a3    a4     d   g u
#> par 1.721  0  0 0.954 0.171 0.1 1
#> 
#> $Item.3
#>        a1 a2     a3 a4      d   g u
#> par 2.415  0 -0.459  0 -2.602 0.1 1
#> 
#> $Item.4
#>        a1 a2 a3    a4      d   g u
#> par 0.745  0  0 0.695 -0.989 0.1 1
#> 
#> $Item.5
#>        a1 a2 a3    a4     d   g u
#> par 1.048  0  0 0.603 0.419 0.1 1
#> 
#> $Item.6
#>       a1 a2    a3 a4      d   g u
#> par 3.06  0 0.501  0 -5.002 0.1 1
#> 
#> $Item.7
#>        a1    a2 a3 a4     d   g u
#> par 1.121 0.839  0  0 1.373 0.1 1
#> 
#> $Item.8
#>        a1 a2    a3 a4      d   g u
#> par 1.956  0 1.443  0 -3.772 0.1 1
#> 
#> $Item.9
#>        a1    a2 a3 a4     d   g u
#> par 0.512 1.236  0  0 2.484 0.1 1
#> 
#> $Item.10
#>       a1    a2 a3 a4      d   g u
#> par 1.68 1.506  0  0 -1.031 0.1 1
#> 
#> $Item.11
#>        a1    a2 a3 a4     d   g u
#> par 1.655 0.842  0  0 5.441 0.1 1
#> 
#> $Item.12
#>        a1 a2 a3    a4      d   g u
#> par 0.129  0  0 0.364 -0.641 0.1 1
#> 
#> $Item.13
#>        a1    a2 a3 a4     d   g u
#> par 1.183 0.477  0  0 0.679 0.1 1
#> 
#> $Item.14
#>        a1 a2 a3    a4     d   g u
#> par 1.125  0  0 1.058 1.164 0.1 1
#> 
#> $Item.15
#>        a1    a2 a3 a4     d   g u
#> par 1.435 0.317  0  0 1.863 0.1 1
#> 
#> $Item.16
#>       a1 a2    a3 a4      d   g u
#> par 0.95  0 0.573  0 -0.783 0.1 1
#> 
#> $Item.17
#>        a1    a2 a3 a4     d   g u
#> par 1.547 0.059  0  0 4.112 0.1 1
#> 
#> $Item.18
#>        a1    a2 a3 a4      d   g u
#> par 2.731 0.094  0  0 -1.808 0.1 1
#> 
#> $Item.19
#>        a1 a2 a3    a4      d   g u
#> par 0.918  0  0 0.101 -0.001 0.1 1
#> 
#> $Item.20
#>        a1 a2 a3    a4     d   g u
#> par 1.456  0  0 0.593 2.501 0.1 1
#> 
#> $Item.21
#>        a1    a2 a3 a4    d   g u
#> par 0.596 0.493  0  0 2.49 0.1 1
#> 
#> $Item.22
#>        a1     a2 a3 a4     d   g u
#> par 1.554 -0.242  0  0 3.428 0.1 1
#> 
#> $Item.23
#>        a1 a2 a3    a4      d   g u
#> par 0.908  0  0 0.766 -1.488 0.1 1
#> 
#> $Item.24
#>        a1    a2 a3 a4     d   g u
#> par 1.379 0.001  0  0 1.132 0.1 1
#> 
#> $Item.25
#>       a1 a2 a3    a4      d   g u
#> par 1.03  0  0 1.094 -1.164 0.1 1
#> 
#> $Item.26
#>        a1 a2 a3    a4      d   g u
#> par 1.985  0  0 0.747 -0.663 0.1 1
#> 
#> $Item.27
#>        a1    a2 a3 a4     d   g u
#> par 1.909 0.348  0  0 2.642 0.1 1
#> 
#> $Item.28
#>        a1 a2 a3    a4      d   g u
#> par 1.213  0  0 0.142 -0.097 0.1 1
#> 
#> $Item.29
#>        a1 a2    a3 a4      d   g u
#> par 1.938  0 2.339  0 -2.209 0.1 1
#> 
#> $Item.30
#>        a1 a2 a3     a4      d   g u
#> par 0.479  0  0 -0.128 -0.527 0.1 1
#> 
#> $Item.31
#>        a1    a2 a3 a4     d   g u
#> par 3.173 -0.82  0  0 3.316 0.1 1
#> 
#> $Item.32
#>        a1 a2     a3 a4      d   g u
#> par 0.534  0 -0.053  0 -2.786 0.1 1
#> 
#> $GroupPars
#>     MEAN_1 MEAN_2 MEAN_3 MEAN_4 COV_11 COV_21 COV_31 COV_41 COV_22 COV_32
#> par      0      0      0      0      1      0      0      0      1      0
#>     COV_42 COV_33 COV_43 COV_44
#> par      0      1      0      1
#> 
anova(mod1, mod2)
#>           AIC    SABIC       HQ      BIC    logLik     X2 df   p
#> mod1 19062.10 19179.44 19226.42 19484.21 -9435.052              
#> mod2 19009.55 19126.88 19173.87 19431.65 -9408.775 52.553  0 NaN

## don't estimate specific factor for item 32
specific[32] <- NA
mod3 <- bfactor(data, specific)
#> 
anova(mod3, mod1)
#>           AIC    SABIC       HQ      BIC    logLik   X2 df     p
#> mod3 19060.12 19176.23 19222.73 19477.83 -9435.062              
#> mod1 19062.10 19179.44 19226.42 19484.21 -9435.052 0.02  1 0.886

# same, but with syntax (not run)
specific3 <- "S1 = 7,9,10,11,13,15,17,18,21,22,24,27,31
              S2 = 1,3,6,8,16,29
              S3 = 2,4,5,12,14,19,20,23,25,26,28,30"
# mod3b <- bfactor(data, specific3)
# anova(mod3b)


#########
# mixed itemtype example

# simulate data
a <- matrix(c(
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5),ncol=3,byrow=TRUE)

d <- matrix(c(
-1.0,NA,NA,
-1.5,NA,NA,
 1.5,NA,NA,
 0.0,NA,NA,
2.5,1.0,-1,
3.0,2.0,-0.5,
3.0,2.0,-0.5,
3.0,2.0,-0.5,
2.5,1.0,-1,
2.0,0.0,NA,
-1.0,NA,NA,
-1.5,NA,NA,
 1.5,NA,NA,
 0.0,NA,NA),ncol=3,byrow=TRUE)
items <- rep('2PL', 14)
items[5:10] <- 'graded'

sigma <- diag(3)
dataset <- simdata(a,d,5000,itemtype=items,sigma=sigma)
itemstats(dataset)
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  5000           15.116          4.495 0.175 0.034 0.733     2.323
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  5000 2 0.311 0.463   0.434         0.345       0.720
#> Item_2  5000 2 0.229 0.420   0.386         0.302       0.724
#> Item_3  5000 2 0.770 0.421   0.378         0.293       0.724
#> Item_4  5000 2 0.493 0.500   0.455         0.360       0.718
#> Item_5  5000 4 1.893 0.975   0.585         0.414       0.711
#> Item_6  5000 4 2.150 0.883   0.537         0.374       0.716
#> Item_7  5000 4 2.168 0.881   0.538         0.376       0.715
#> Item_8  5000 4 2.138 0.889   0.549         0.388       0.714
#> Item_9  5000 4 1.857 0.978   0.586         0.413       0.711
#> Item_10 5000 3 1.327 0.746   0.523         0.387       0.713
#> Item_11 5000 2 0.304 0.460   0.427         0.338       0.721
#> Item_12 5000 2 0.234 0.424   0.403         0.320       0.723
#> Item_13 5000 2 0.752 0.432   0.417         0.333       0.722
#> Item_14 5000 2 0.491 0.500   0.457         0.362       0.718
#> 
#> $proportions
#>             0     1     2     3
#> Item_1  0.689 0.311    NA    NA
#> Item_2  0.771 0.229    NA    NA
#> Item_3  0.230 0.770    NA    NA
#> Item_4  0.507 0.493    NA    NA
#> Item_5  0.114 0.194 0.378 0.314
#> Item_6  0.079 0.089 0.434 0.397
#> Item_7  0.079 0.083 0.431 0.408
#> Item_8  0.081 0.093 0.433 0.393
#> Item_9  0.115 0.214 0.370 0.301
#> Item_10 0.168 0.337 0.495    NA
#> Item_11 0.696 0.304    NA    NA
#> Item_12 0.766 0.234    NA    NA
#> Item_13 0.248 0.752    NA    NA
#> Item_14 0.509 0.491    NA    NA
#> 

specific <- "S1 = 1-7
             S2 = 8-14"
simmod <- bfactor(dataset, specific)
#> 
coef(simmod, simplify=TRUE)
#> $items
#>            a1    a2    a3      d  g  u    d1     d2     d3
#> Item_1  1.032 0.505 0.000 -1.003  0  1    NA     NA     NA
#> Item_2  0.953 0.474 0.000 -1.477  0  1    NA     NA     NA
#> Item_3  0.868 0.494 0.000  1.442  0  1    NA     NA     NA
#> Item_4  1.002 0.406 0.000 -0.037  0  1    NA     NA     NA
#> Item_5  1.040 0.578 0.000     NA NA NA 2.536  1.029 -0.998
#> Item_6  0.963 0.465 0.000     NA NA NA 2.899  1.928 -0.514
#> Item_7  0.969 0.504 0.000     NA NA NA 2.929  2.009 -0.462
#> Item_8  1.011 0.000 0.680     NA NA NA 2.993  1.974 -0.563
#> Item_9  1.023 0.000 0.469     NA NA NA 2.471  0.888 -1.054
#> Item_10 1.002 0.000 0.461     NA NA NA 1.951 -0.023     NA
#> Item_11 1.001 0.000 0.543 -1.041  0  1    NA     NA     NA
#> Item_12 1.039 0.000 0.603 -1.502  0  1    NA     NA     NA
#> Item_13 1.033 0.000 0.429  1.377  0  1    NA     NA     NA
#> Item_14 1.017 0.000 0.401 -0.046  0  1    NA     NA     NA
#> 
#> $means
#>  G S1 S2 
#>  0  0  0 
#> 
#> $cov
#>    G S1 S2
#> G  1  0  0
#> S1 0  1  0
#> S2 0  0  1
#> 


#########
# General testlet response model (Wainer, 2007)

# simulate data
set.seed(1234)
a <- matrix(0, 12, 4)
a[,1] <- rlnorm(12, .2, .3)
ind <- 1
for(i in 1:3){
   a[ind:(ind+3),i+1] <- a[ind:(ind+3),1]
   ind <- ind+4
}
print(a)
#>            [,1]      [,2]     [,3]      [,4]
#>  [1,] 0.8503394 0.8503394 0.000000 0.0000000
#>  [2,] 1.3274088 1.3274088 0.000000 0.0000000
#>  [3,] 1.6910208 1.6910208 0.000000 0.0000000
#>  [4,] 0.6042850 0.6042850 0.000000 0.0000000
#>  [5,] 1.3892130 0.0000000 1.389213 0.0000000
#>  [6,] 1.4216480 0.0000000 1.421648 0.0000000
#>  [7,] 1.0279618 0.0000000 1.027962 0.0000000
#>  [8,] 1.0366667 0.0000000 1.036667 0.0000000
#>  [9,] 1.0311394 0.0000000 0.000000 1.0311394
#> [10,] 0.9351846 0.0000000 0.000000 0.9351846
#> [11,] 1.0584888 0.0000000 0.000000 1.0584888
#> [12,] 0.9052755 0.0000000 0.000000 0.9052755
d <- rnorm(12, 0, .5)
sigma <- diag(c(1, .5, 1, .5))
dataset <- simdata(a,d,2000,itemtype=rep('2PL', 12),sigma=sigma)
itemstats(dataset)
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  2000                6          2.929 0.175 0.068 0.717     1.558
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  2000 2 0.426 0.495   0.438         0.287       0.708
#> Item_2  2000 2 0.502 0.500   0.560         0.425       0.689
#> Item_3  2000 2 0.571 0.495   0.575         0.445       0.686
#> Item_4  2000 2 0.502 0.500   0.383         0.224       0.716
#> Item_5  2000 2 0.464 0.499   0.549         0.413       0.690
#> Item_6  2000 2 0.436 0.496   0.561         0.428       0.688
#> Item_7  2000 2 0.440 0.497   0.500         0.356       0.698
#> Item_8  2000 2 0.693 0.462   0.474         0.339       0.701
#> Item_9  2000 2 0.511 0.500   0.481         0.334       0.701
#> Item_10 2000 2 0.456 0.498   0.465         0.316       0.704
#> Item_11 2000 2 0.458 0.498   0.459         0.309       0.705
#> Item_12 2000 2 0.540 0.498   0.475         0.327       0.702
#> 
#> $proportions
#>             0     1
#> Item_1  0.575 0.426
#> Item_2  0.498 0.502
#> Item_3  0.430 0.571
#> Item_4  0.498 0.502
#> Item_5  0.536 0.464
#> Item_6  0.564 0.436
#> Item_7  0.559 0.440
#> Item_8  0.308 0.693
#> Item_9  0.488 0.511
#> Item_10 0.543 0.456
#> Item_11 0.541 0.458
#> Item_12 0.460 0.540
#> 

# estimate by applying constraints and freeing the latent variances
specific <- "S1 = 1-4
             S2 = 5-8
             S3 = 9-12"
model <- "G = 1-12
          CONSTRAIN = (1, a1, a2), (2, a1, a2), (3, a1, a2), (4, a1, a2),
            (5, a1, a3), (6, a1, a3), (7, a1, a3), (8, a1, a3),
            (9, a1, a4), (10, a1, a4), (11, a1, a4), (12, a1, a4)
          COV = S1*S1, S2*S2, S3*S3"

simmod <- bfactor(dataset, specific, model)
#> 
coef(simmod, simplify=TRUE)
#> $items
#>            a1    a2    a3    a4      d g u
#> Item_1  0.794 0.794 0.000 0.000 -0.359 0 1
#> Item_2  1.544 1.544 0.000 0.000  0.011 0 1
#> Item_3  1.762 1.762 0.000 0.000  0.479 0 1
#> Item_4  0.544 0.544 0.000 0.000  0.011 0 1
#> Item_5  1.386 0.000 1.386 0.000 -0.244 0 1
#> Item_6  1.497 0.000 1.497 0.000 -0.449 0 1
#> Item_7  0.853 0.000 0.853 0.000 -0.312 0 1
#> Item_8  0.953 0.000 0.953 0.000  1.101 0 1
#> Item_9  0.981 0.000 0.000 0.981  0.058 0 1
#> Item_10 0.913 0.000 0.000 0.913 -0.217 0 1
#> Item_11 0.868 0.000 0.000 0.868 -0.204 0 1
#> Item_12 0.966 0.000 0.000 0.966  0.206 0 1
#> 
#> $means
#>  G S1 S2 S3 
#>  0  0  0  0 
#> 
#> $cov
#>    G    S1    S2    S3
#> G  1 0.000 0.000 0.000
#> S1 0 0.452 0.000 0.000
#> S2 0 0.000 1.135 0.000
#> S3 0 0.000 0.000 0.432
#> 

# Constrained testlet model (Bradlow, 1999)
model2 <- "G = 1-12
          CONSTRAIN = (1, a1, a2), (2, a1, a2), (3, a1, a2), (4, a1, a2),
            (5, a1, a3), (6, a1, a3), (7, a1, a3), (8, a1, a3),
            (9, a1, a4), (10, a1, a4), (11, a1, a4), (12, a1, a4),
            (GROUP, COV_22, COV_33, COV_44)
          COV = S1*S1, S2*S2, S3*S3"

simmod2 <- bfactor(dataset, specific, model2)
#> 
coef(simmod2, simplify=TRUE)
#> $items
#>            a1    a2    a3    a4      d g u
#> Item_1  0.744 0.744 0.000 0.000 -0.360 0 1
#> Item_2  1.453 1.453 0.000 0.000  0.010 0 1
#> Item_3  1.664 1.664 0.000 0.000  0.482 0 1
#> Item_4  0.509 0.509 0.000 0.000  0.011 0 1
#> Item_5  1.541 0.000 1.541 0.000 -0.241 0 1
#> Item_6  1.670 0.000 1.670 0.000 -0.445 0 1
#> Item_7  0.968 0.000 0.968 0.000 -0.313 0 1
#> Item_8  1.075 0.000 1.075 0.000  1.098 0 1
#> Item_9  0.927 0.000 0.000 0.927  0.059 0 1
#> Item_10 0.854 0.000 0.000 0.854 -0.218 0 1
#> Item_11 0.813 0.000 0.000 0.813 -0.205 0 1
#> Item_12 0.908 0.000 0.000 0.908  0.207 0 1
#> 
#> $means
#>  G S1 S2 S3 
#>  0  0  0  0 
#> 
#> $cov
#>    G    S1    S2    S3
#> G  1 0.000 0.000 0.000
#> S1 0 0.667 0.000 0.000
#> S2 0 0.000 0.667 0.000
#> S3 0 0.000 0.000 0.667
#> 
anova(simmod2, simmod)
#>              AIC    SABIC       HQ      BIC   logLik     X2 df     p
#> simmod2 30256.59 30317.19 30308.00 30396.61 -15103.3                
#> simmod  30248.79 30314.24 30304.32 30400.02 -15097.4 11.795  2 0.003


#########
# Two-tier model

# simulate data
set.seed(1234)
a <- matrix(c(
  1,0,0.5,NA,NA,
  1,0,0.5,NA,NA,
  1,0,0.5,NA,NA,
  1,0,0.5,NA,NA,
  1,0,0.5,NA,NA,
  1,0,NA,0.5,NA,
  1,0,NA,0.5,NA,
  1,0,NA,0.5,NA,
  0,1,NA,0.5,NA,
  0,1,NA,0.5,NA,
  0,1,NA,0.5,NA,
  0,1,NA,NA,0.5,
  0,1,NA,NA,0.5,
  0,1,NA,NA,0.5,
  0,1,NA,NA,0.5,
  0,1,NA,NA,0.5),ncol=5,byrow=TRUE)

d <- matrix(rnorm(16))
items <- rep('2PL', 16)

sigma <- diag(5)
sigma[1,2] <- sigma[2,1] <- .4
dataset <- simdata(a,d,2000,itemtype=items,sigma=sigma)
itemstats(dataset)
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  2000            7.049          3.089  0.11 0.055 0.666     1.786
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  2000 2 0.280 0.449   0.400         0.268       0.651
#> Item_2  2000 2 0.565 0.496   0.399         0.252       0.653
#> Item_3  2000 2 0.697 0.460   0.378         0.241       0.655
#> Item_4  2000 2 0.123 0.329   0.283         0.181       0.661
#> Item_5  2000 2 0.591 0.492   0.395         0.249       0.654
#> Item_6  2000 2 0.582 0.493   0.424         0.280       0.650
#> Item_7  2000 2 0.366 0.482   0.430         0.291       0.648
#> Item_8  2000 2 0.390 0.488   0.403         0.258       0.653
#> Item_9  2000 2 0.394 0.489   0.407         0.263       0.652
#> Item_10 2000 2 0.316 0.465   0.413         0.277       0.650
#> Item_11 2000 2 0.399 0.490   0.442         0.302       0.647
#> Item_12 2000 2 0.320 0.467   0.420         0.284       0.649
#> Item_13 2000 2 0.363 0.481   0.430         0.291       0.648
#> Item_14 2000 2 0.501 0.500   0.428         0.282       0.649
#> Item_15 2000 2 0.672 0.470   0.405         0.266       0.651
#> Item_16 2000 2 0.491 0.500   0.445         0.302       0.647
#> 
#> $proportions
#>             0     1
#> Item_1  0.720 0.280
#> Item_2  0.435 0.565
#> Item_3  0.303 0.697
#> Item_4  0.877 0.123
#> Item_5  0.410 0.591
#> Item_6  0.418 0.582
#> Item_7  0.634 0.366
#> Item_8  0.610 0.390
#> Item_9  0.606 0.394
#> Item_10 0.684 0.316
#> Item_11 0.601 0.399
#> Item_12 0.680 0.320
#> Item_13 0.637 0.363
#> Item_14 0.499 0.501
#> Item_15 0.328 0.672
#> Item_16 0.509 0.491
#> 

specific <- "S1 = 1-5
             S2 = 6-11
             S3 = 12-16"
model <- '
    G1 = 1-8
    G2 = 9-16
    COV = G1*G2'

# quadpts dropped for faster estimation, but not as precise
simmod <- bfactor(dataset, specific, model, quadpts = 15, TOL = 1e-3)
#> 
coef(simmod, simplify=TRUE)
#> $items
#>            a1    a2    a3    a4    a5      d g u
#> Item_1  1.109 0.000 1.040 0.000 0.000 -1.330 0 1
#> Item_2  0.990 0.000 0.128 0.000 0.000  0.316 0 1
#> Item_3  0.906 0.000 0.362 0.000 0.000  0.994 0 1
#> Item_4  0.995 0.000 0.352 0.000 0.000 -2.341 0 1
#> Item_5  0.915 0.000 0.494 0.000 0.000  0.447 0 1
#> Item_6  1.011 0.000 0.000 0.597 0.000  0.419 0 1
#> Item_7  1.088 0.000 0.000 0.589 0.000 -0.714 0 1
#> Item_8  0.926 0.000 0.000 0.340 0.000 -0.539 0 1
#> Item_9  0.000 0.936 0.000 0.354 0.000 -0.521 0 1
#> Item_10 0.000 0.908 0.000 0.696 0.000 -0.971 0 1
#> Item_11 0.000 1.016 0.000 0.342 0.000 -0.507 0 1
#> Item_12 0.000 1.058 0.000 0.000 0.490 -0.955 0 1
#> Item_13 0.000 1.147 0.000 0.000 0.046 -0.711 0 1
#> Item_14 0.000 1.070 0.000 0.000 0.095  0.003 0 1
#> Item_15 0.000 1.047 0.000 0.000 0.536  0.910 0 1
#> Item_16 0.000 1.136 0.000 0.000 0.840 -0.052 0 1
#> 
#> $means
#> G1 G2 S1 S2 S3 
#>  0  0  0  0  0 
#> 
#> $cov
#>       G1    G2 S1 S2 S3
#> G1 1.000 0.407  0  0  0
#> G2 0.407 1.000  0  0  0
#> S1 0.000 0.000  1  0  0
#> S2 0.000 0.000  0  1  0
#> S3 0.000 0.000  0  0  1
#> 
summary(simmod)
#>            G1    G2    S1    S2    S3    h2
#> Item_1  0.486       0.456             0.444
#> Item_2  0.502       0.065             0.256
#> Item_3  0.462       0.185             0.247
#> Item_4  0.497       0.176             0.278
#> Item_5  0.459       0.248             0.272
#> Item_6  0.489             0.289       0.322
#> Item_7  0.517             0.280       0.346
#> Item_8  0.471             0.173       0.252
#> Item_9        0.474       0.179       0.257
#> Item_10       0.443       0.339       0.311
#> Item_11       0.505       0.170       0.284
#> Item_12       0.513             0.237 0.319
#> Item_13       0.559             0.023 0.313
#> Item_14       0.531             0.047 0.285
#> Item_15       0.506             0.259 0.323
#> Item_16       0.514             0.380 0.408
#> 
#> SS loadings:  1.886 2.053 0.338 0.368 0.27 
#> Proportion Var:  0.118 0.128 0.021 0.023 0.017 
#> 
#> Factor correlations: 
#> 
#>       G1 G2 S1 S2 S3
#> G1 1.000            
#> G2 0.407  1         
#> S1 0.000  0  1      
#> S2 0.000  0  0  1   
#> S3 0.000  0  0  0  1
itemfit(simmod, QMC=TRUE)
#>       item   S_X2 df.S_X2 RMSEA.S_X2 p.S_X2
#> 1   Item_1 14.270       9      0.017  0.113
#> 2   Item_2  7.004      10      0.000  0.725
#> 3   Item_3 13.697       9      0.016  0.134
#> 4   Item_4 11.930      10      0.010  0.290
#> 5   Item_5 22.962      10      0.025  0.011
#> 6   Item_6 16.035      10      0.017  0.099
#> 7   Item_7  7.401       9      0.000  0.595
#> 8   Item_8 10.845      10      0.007  0.370
#> 9   Item_9 18.119      10      0.020  0.053
#> 10 Item_10 13.087      10      0.012  0.219
#> 11 Item_11 13.155      10      0.013  0.215
#> 12 Item_12 15.821       9      0.019  0.071
#> 13 Item_13 19.170      10      0.021  0.038
#> 14 Item_14 15.028      10      0.016  0.131
#> 15 Item_15 22.121       9      0.027  0.009
#> 16 Item_16  6.090      10      0.000  0.808
M2(simmod, QMC=TRUE)
#>           M2 df     p RMSEA RMSEA_5 RMSEA_95 SRMSR   TLI CFI
#> stats 73.523 87 0.848     0       0    0.007 0.016 1.005   1
residuals(simmod, QMC=TRUE)
#> LD matrix (lower triangle) and standardized residual correlations (upper triangle)
#> 
#> Upper triangle summary:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -0.048  -0.010   0.001   0.000   0.010   0.041 
#> 
#>         Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> Item_1         -0.003  0.003 -0.003  0.000 -0.002 -0.016  0.003  0.004   0.019
#> Item_2   0.017        -0.015  0.015  0.008  0.012 -0.010  0.015 -0.030   0.020
#> Item_3   0.018  0.473         0.004 -0.002  0.029  0.006 -0.029  0.010   0.005
#> Item_4   0.013  0.478  0.030        -0.004 -0.027  0.009  0.027 -0.014   0.009
#> Item_5   0.000  0.135  0.008  0.026        -0.029  0.017  0.006 -0.039  -0.001
#> Item_6   0.005  0.267  1.703  1.510  1.686         0.007 -0.001  0.002   0.005
#> Item_7   0.486  0.203  0.082  0.151  0.607  0.097        -0.010 -0.014   0.004
#> Item_8   0.014  0.435  1.666  1.446  0.072  0.002  0.203        -0.021   0.017
#> Item_9   0.031  1.841  0.203  0.382  3.044  0.005  0.368  0.853         -0.001
#> Item_10  0.708  0.800  0.043  0.148  0.002  0.056  0.026  0.602  0.002        
#> Item_11  3.398  0.026  0.001  0.032  0.279  0.166  0.838  0.071  0.359   0.137
#> Item_12  0.051  1.006  0.000  0.391  0.249  0.087  0.477  0.286  0.004   0.032
#> Item_13  0.927  2.865  0.198  0.024  0.218  0.002  3.025  0.924  0.585   0.103
#> Item_14  0.158  0.180  0.088  1.854  0.181  0.440  0.457  0.229  0.206   0.065
#> Item_15  0.275  0.701  0.004  4.564  1.125  0.122  3.620  0.122  0.172   0.006
#> Item_16  0.101  0.001  0.500  0.426  0.766  0.018  0.024  0.848  0.774   0.255
#>         Item_11 Item_12 Item_13 Item_14 Item_15 Item_16
#> Item_1    0.041   0.005   0.022   0.009   0.012  -0.007
#> Item_2    0.004   0.022  -0.038  -0.009  -0.019   0.001
#> Item_3   -0.001   0.000   0.010   0.007   0.001   0.016
#> Item_4    0.004  -0.014  -0.003  -0.030  -0.048  -0.015
#> Item_5    0.012   0.011  -0.010  -0.010   0.024   0.020
#> Item_6   -0.009   0.007  -0.001   0.015  -0.008   0.003
#> Item_7    0.020   0.015   0.039  -0.015  -0.043   0.003
#> Item_8    0.006  -0.012  -0.021  -0.011   0.008   0.021
#> Item_9    0.013  -0.001  -0.017   0.010  -0.009   0.020
#> Item_10  -0.008  -0.004   0.007  -0.006   0.002  -0.011
#> Item_11          -0.008  -0.003  -0.024   0.029  -0.014
#> Item_12   0.124          -0.014   0.012  -0.003  -0.001
#> Item_13   0.015   0.381           0.019  -0.007   0.006
#> Item_14   1.200   0.284   0.718           0.001  -0.010
#> Item_15   1.664   0.014   0.090   0.003           0.002
#> Item_16   0.398   0.004   0.080   0.203   0.005        

# EAP predictions for all factors (high dimensional)
eaps_all <- fscores(simmod, QMC=TRUE, quadpts=50000)
head(eaps_all)
#>              G1          G2         S1           S2          S3
#> [1,]  0.5191147 -0.05812303 -0.2298246  0.421019053 -0.02374122
#> [2,] -0.7348744 -0.95615920 -0.2423637 -0.101862792 -0.41703218
#> [3,] -0.6147126  1.19712303 -0.4783205 -0.472646037  0.44121560
#> [4,]  0.9284838  0.59274153 -0.2548899  0.986259541 -0.25352275
#> [5,]  0.6460876  1.29373519  0.6411221 -0.005344923  0.38031801
#> [6,] -0.7212222 -1.02053207 -0.4440952  0.189161622 -0.47899054
maps <- fscores(simmod, method = 'MAP')
head(maps)
#>              G1          G2         S1           S2          S3
#> [1,]  0.5135992 -0.03912492 -0.2056159  0.419322965 -0.02778409
#> [2,] -0.6774244 -0.87671079 -0.2194985 -0.077330672 -0.39920301
#> [3,] -0.5799759  1.15757307 -0.4439183 -0.448083235  0.41194310
#> [4,]  0.8983374  0.57134257 -0.2333886  0.967336777 -0.26136242
#> [5,]  0.6156947  1.23660200  0.6409160 -0.005859276  0.35081023
#> [6,] -0.6658814 -0.93907518 -0.4207128  0.209132981 -0.45956922

######
# Multiple-group bifactor example from Cai, Yang, Hansen (2011)

# Table 1 info from Cai, Yang, and Hansen (2011)
intercept <- c(1, .25, -.25, -1, 1, .25, -.25, -1,1,
               .25, -.25, -1,1, .25, -.25, -1)
theta0 <- c(1, 1.4,1.7,2,
            1.4,1.7,2,1,
            1.7,2,1,1.4,
            2, 1, 1.4, 1.7)
thetan <- c(.8,1.5,1.2,1,
            1,.8,1.5,1.2,
            1.2,1,.8,1.5,
            1.5,1.2,1,.8)
thetaN <- matrix(0, 16, 4)
thetaN[1:4, 1] <- thetan[1:4]
thetaN[1:4+4, 2] <- thetan[1:4+4]
thetaN[1:4+8, 3] <- thetan[1:4+8]
thetaN[1:4+12, 4] <- thetan[1:4+12]

as <- cbind(theta0, thetaN)
as
#>       theta0                
#>  [1,]    1.0 0.8 0.0 0.0 0.0
#>  [2,]    1.4 1.5 0.0 0.0 0.0
#>  [3,]    1.7 1.2 0.0 0.0 0.0
#>  [4,]    2.0 1.0 0.0 0.0 0.0
#>  [5,]    1.4 0.0 1.0 0.0 0.0
#>  [6,]    1.7 0.0 0.8 0.0 0.0
#>  [7,]    2.0 0.0 1.5 0.0 0.0
#>  [8,]    1.0 0.0 1.2 0.0 0.0
#>  [9,]    1.7 0.0 0.0 1.2 0.0
#> [10,]    2.0 0.0 0.0 1.0 0.0
#> [11,]    1.0 0.0 0.0 0.8 0.0
#> [12,]    1.4 0.0 0.0 1.5 0.0
#> [13,]    2.0 0.0 0.0 0.0 1.5
#> [14,]    1.0 0.0 0.0 0.0 1.2
#> [15,]    1.4 0.0 0.0 0.0 1.0
#> [16,]    1.7 0.0 0.0 0.0 0.8

# data generation for focal group in publication does not have response
# for items 13-16. However, full-data approach presented first
N <- 1000
itemtype <- '2PL'
gmeans <- c(1, -.5, 0, .5, 0)
sigma <- diag(c(.8, 1.2, 1.5, 1, 1))

datG1 <- simdata(as, intercept, N=N, itemtype='2PL')
datG2 <- simdata(as, intercept, N=N, itemtype='2PL',
  mu = gmeans, sigma = sigma)

dat <- rbind(datG1, datG2)
group <- rep(c('G1', 'G2'), each=N)

specific <- "S1 = 1-4
             S2 = 5-8
             S3 = 9-12
             S4 = 13-16"

mod <- bfactor(dat, specific, group=group, SE=TRUE,
  invariance=c('free_means', 'free_vars', colnames(dat)))
#> 
#> 
#> Calculating information matrix...
coef(mod, simplify=TRUE)
#> $G1
#> $items
#>            a1    a2    a3    a4    a5      d g u
#> Item_1  1.053 0.609 0.000 0.000 0.000  0.793 0 1
#> Item_2  1.524 1.392 0.000 0.000 0.000  0.237 0 1
#> Item_3  1.991 1.262 0.000 0.000 0.000 -0.170 0 1
#> Item_4  2.398 1.070 0.000 0.000 0.000 -1.038 0 1
#> Item_5  1.447 0.000 0.779 0.000 0.000  1.053 0 1
#> Item_6  1.760 0.000 0.702 0.000 0.000  0.228 0 1
#> Item_7  1.951 0.000 1.416 0.000 0.000 -0.310 0 1
#> Item_8  1.001 0.000 1.103 0.000 0.000 -0.975 0 1
#> Item_9  1.750 0.000 0.000 1.053 0.000  0.847 0 1
#> Item_10 1.990 0.000 0.000 0.928 0.000  0.147 0 1
#> Item_11 1.077 0.000 0.000 0.965 0.000 -0.319 0 1
#> Item_12 1.759 0.000 0.000 2.115 0.000 -1.305 0 1
#> Item_13 2.104 0.000 0.000 0.000 1.606  1.129 0 1
#> Item_14 1.122 0.000 0.000 0.000 1.068  0.333 0 1
#> Item_15 1.408 0.000 0.000 0.000 1.041 -0.142 0 1
#> Item_16 1.752 0.000 0.000 0.000 0.985 -1.015 0 1
#> 
#> $means
#>  G S1 S2 S3 S4 
#>  0  0  0  0  0 
#> 
#> $cov
#>    G S1 S2 S3 S4
#> G  1  0  0  0  0
#> S1 0  1  0  0  0
#> S2 0  0  1  0  0
#> S3 0  0  0  1  0
#> S4 0  0  0  0  1
#> 
#> 
#> $G2
#> $items
#>            a1    a2    a3    a4    a5      d g u
#> Item_1  1.053 0.609 0.000 0.000 0.000  0.793 0 1
#> Item_2  1.524 1.392 0.000 0.000 0.000  0.237 0 1
#> Item_3  1.991 1.262 0.000 0.000 0.000 -0.170 0 1
#> Item_4  2.398 1.070 0.000 0.000 0.000 -1.038 0 1
#> Item_5  1.447 0.000 0.779 0.000 0.000  1.053 0 1
#> Item_6  1.760 0.000 0.702 0.000 0.000  0.228 0 1
#> Item_7  1.951 0.000 1.416 0.000 0.000 -0.310 0 1
#> Item_8  1.001 0.000 1.103 0.000 0.000 -0.975 0 1
#> Item_9  1.750 0.000 0.000 1.053 0.000  0.847 0 1
#> Item_10 1.990 0.000 0.000 0.928 0.000  0.147 0 1
#> Item_11 1.077 0.000 0.000 0.965 0.000 -0.319 0 1
#> Item_12 1.759 0.000 0.000 2.115 0.000 -1.305 0 1
#> Item_13 2.104 0.000 0.000 0.000 1.606  1.129 0 1
#> Item_14 1.122 0.000 0.000 0.000 1.068  0.333 0 1
#> Item_15 1.408 0.000 0.000 0.000 1.041 -0.142 0 1
#> Item_16 1.752 0.000 0.000 0.000 0.985 -1.015 0 1
#> 
#> $means
#>      G     S1     S2     S3     S4 
#>  0.979 -0.766 -0.065  0.414 -0.179 
#> 
#> $cov
#>        G    S1    S2    S3    S4
#> G  0.684 0.000 0.000 0.000 0.000
#> S1 0.000 1.268 0.000 0.000 0.000
#> S2 0.000 0.000 2.116 0.000 0.000
#> S3 0.000 0.000 0.000 0.996 0.000
#> S4 0.000 0.000 0.000 0.000 1.005
#> 
#> 

# DIF testing for general dimension only with likelihood ratio tests (not run)
if(FALSE){
  mirtCluster()  # allocate computing cores
  DIF(mod, which.par=c('a1', 'd'), scheme='drop')
}

## same analysis, however items 13:16 do not exist in the focal group
datG2[,13:16] <- NA    # use this to match publication
dat <- rbind(datG1, datG2)
head(dat)
#>      Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> [1,]      1      1      1      1      1      1      0      1      1       1
#> [2,]      1      0      1      0      0      1      0      0      0       1
#> [3,]      1      1      1      1      1      1      1      0      1       1
#> [4,]      1      1      1      1      1      1      1      1      1       1
#> [5,]      0      0      0      0      1      0      0      0      1       0
#> [6,]      1      1      1      1      1      1      1      1      1       1
#>      Item_11 Item_12 Item_13 Item_14 Item_15 Item_16
#> [1,]       1       1       1       1       1       1
#> [2,]       0       0       1       0       1       1
#> [3,]       1       1       1       1       1       1
#> [4,]       1       0       0       0       0       0
#> [5,]       0       0       0       0       0       0
#> [6,]       0       0       1       1       1       0
tail(dat)
#>         Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> [1995,]      1      1      1      1      1      1      1      1      1       1
#> [1996,]      0      1      1      1      1      1      0      0      0       1
#> [1997,]      0      1      1      1      1      1      1      0      1       1
#> [1998,]      0      1      1      1      1      1      1      1      1       1
#> [1999,]      1      1      0      1      1      1      0      0      1       1
#> [2000,]      1      0      1      1      1      1      1      1      1       1
#>         Item_11 Item_12 Item_13 Item_14 Item_15 Item_16
#> [1995,]       1       1      NA      NA      NA      NA
#> [1996,]       0       0      NA      NA      NA      NA
#> [1997,]       1       1      NA      NA      NA      NA
#> [1998,]       1       1      NA      NA      NA      NA
#> [1999,]       1       1      NA      NA      NA      NA
#> [2000,]       1       1      NA      NA      NA      NA

# specify mean/cov structure explicitly
model2 <- "G = 1-16
           MEAN [G2] = G, S1, S2, S3
           COV [G2] = G*G, S1*S1, S2*S2, S3*S3"

mod2 <- bfactor(dat, specific, model2, group=group, invariance=colnames(dat))
#> 
coef(mod2, simplify=TRUE)
#> $G1
#> $items
#>            a1    a2    a3    a4    a5      d g u
#> Item_1  1.043 0.615 0.000 0.000 0.000  0.799 0 1
#> Item_2  1.532 1.337 0.000 0.000 0.000  0.233 0 1
#> Item_3  2.000 1.262 0.000 0.000 0.000 -0.168 0 1
#> Item_4  2.393 1.102 0.000 0.000 0.000 -1.041 0 1
#> Item_5  1.423 0.000 0.800 0.000 0.000  1.050 0 1
#> Item_6  1.722 0.000 0.715 0.000 0.000  0.214 0 1
#> Item_7  1.951 0.000 1.447 0.000 0.000 -0.310 0 1
#> Item_8  1.014 0.000 1.064 0.000 0.000 -0.969 0 1
#> Item_9  1.732 0.000 0.000 1.039 0.000  0.834 0 1
#> Item_10 1.987 0.000 0.000 0.897 0.000  0.124 0 1
#> Item_11 1.088 0.000 0.000 0.926 0.000 -0.317 0 1
#> Item_12 1.845 0.000 0.000 2.279 0.000 -1.349 0 1
#> Item_13 1.839 0.000 0.000 0.000 1.280  1.007 0 1
#> Item_14 1.175 0.000 0.000 0.000 1.068  0.307 0 1
#> Item_15 1.440 0.000 0.000 0.000 1.194 -0.150 0 1
#> Item_16 1.840 0.000 0.000 0.000 1.010 -1.001 0 1
#> 
#> $means
#>  G S1 S2 S3 S4 
#>  0  0  0  0  0 
#> 
#> $cov
#>    G S1 S2 S3 S4
#> G  1  0  0  0  0
#> S1 0  1  0  0  0
#> S2 0  0  1  0  0
#> S3 0  0  0  1  0
#> S4 0  0  0  0  1
#> 
#> 
#> $G2
#> $items
#>            a1    a2    a3    a4    a5      d g u
#> Item_1  1.043 0.615 0.000 0.000 0.000  0.799 0 1
#> Item_2  1.532 1.337 0.000 0.000 0.000  0.233 0 1
#> Item_3  2.000 1.262 0.000 0.000 0.000 -0.168 0 1
#> Item_4  2.393 1.102 0.000 0.000 0.000 -1.041 0 1
#> Item_5  1.423 0.000 0.800 0.000 0.000  1.050 0 1
#> Item_6  1.722 0.000 0.715 0.000 0.000  0.214 0 1
#> Item_7  1.951 0.000 1.447 0.000 0.000 -0.310 0 1
#> Item_8  1.014 0.000 1.064 0.000 0.000 -0.969 0 1
#> Item_9  1.732 0.000 0.000 1.039 0.000  0.834 0 1
#> Item_10 1.987 0.000 0.000 0.897 0.000  0.124 0 1
#> Item_11 1.088 0.000 0.000 0.926 0.000 -0.317 0 1
#> Item_12 1.845 0.000 0.000 2.279 0.000 -1.349 0 1
#> Item_13 1.839 0.000 0.000 0.000 1.280  1.007 0 1
#> Item_14 1.175 0.000 0.000 0.000 1.068  0.307 0 1
#> Item_15 1.440 0.000 0.000 0.000 1.194 -0.150 0 1
#> Item_16 1.840 0.000 0.000 0.000 1.010 -1.001 0 1
#> 
#> $means
#>      G     S1     S2     S3     S4 
#>  1.049 -0.887 -0.150  0.329  0.000 
#> 
#> $cov
#>       G   S1    S2    S3 S4
#> G  0.69 0.00 0.000 0.000  0
#> S1 0.00 1.23 0.000 0.000  0
#> S2 0.00 0.00 2.095 0.000  0
#> S3 0.00 0.00 0.000 0.989  0
#> S4 0.00 0.00 0.000 0.000  1
#> 
#> 

# }
```
