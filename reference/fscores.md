# Compute factor score estimates (a.k.a, ability estimates, latent trait estimates, etc)

Computes MAP, EAP, ML (Embretson & Reise, 2000), EAP for sum-scores
(Thissen et al., 1995), or WLE (Warm, 1989) factor scores with a
multivariate normal prior distribution using equally spaced quadrature.
EAP scores for models with more than three factors are generally not
recommended since the integration grid becomes very large, resulting in
slower estimation and less precision if the `quadpts` are too low.
Therefore, MAP scores should be used instead of EAP scores for higher
dimensional models. Multiple imputation variants are possible for each
estimator if a parameter information matrix was computed, which are
useful if the sample size/number of items were small. As well, if the
model contained latent regression predictors this information will be
used in computing MAP and EAP estimates (for these models,
`full.scores=TRUE` will always be used). Finally, plausible value
imputation is also available, and will also account for latent
regression predictor effects.

## Usage

``` r
fscores(
  object,
  method = "EAP",
  full.scores = TRUE,
  rotate = "oblimin",
  Target = NULL,
  response.pattern = NULL,
  append_response.pattern = FALSE,
  na.rm = FALSE,
  plausible.draws = 0,
  plausible.type = "normal",
  quadpts = NULL,
  item_weights = rep(1, extract.mirt(object, "nitems")),
  returnER = FALSE,
  T_as_X = FALSE,
  EAPsum.scores = FALSE,
  return.acov = FALSE,
  mean = NULL,
  cov = NULL,
  covdata = NULL,
  verbose = interactive(),
  full.scores.SE = FALSE,
  theta_lim = c(-6, 6),
  MI = 0,
  use_dentype_estimate = FALSE,
  QMC = FALSE,
  custom_den = NULL,
  custom_theta = NULL,
  expected.info = FALSE,
  min_expected = 1,
  max_theta = 20,
  start = NULL,
  ...
)
```

## Arguments

- object:

  a computed model object of class `SingleGroupClass`,
  `MultipleGroupClass`, or `DiscreteClass`

- method:

  type of factor score estimation method. Can be:

  - `"EAP"` for the expected a-posteriori (default). For models fit
    using
    [`mdirt`](https://philchalmers.github.io/mirt/reference/mdirt.md)
    this will return the posterior classification probabilities

  - `"MAP"` for the maximum a-posteriori (i.e, Bayes modal)

  - `"ML"` for maximum likelihood

  - `"WLE"` or `"WML"` for weighted maximum-likelihood estimation

  - `"EAPsum"` for the expected a-posteriori for each sum score

  - `"plausible"` for a single plausible value imputation for each case.
    This is equivalent to setting `plausible.draws = 1`

  - `"classify"` for the posteriori classification probabilities (only
    applicable when the input model was of class `MixtureClass`)

- full.scores:

  if `FALSE` then a summary table with factor scores for each unique
  pattern is displayed as a formatted `matrix` object. Otherwise, a
  matrix of factor scores for each response pattern in the data is
  returned (default)

- rotate:

  prior rotation to be used when estimating the factor scores. See
  [`summary-method`](https://philchalmers.github.io/mirt/reference/summary-method.md)
  for details. If the object is not an exploratory model then this
  argument is ignored

- Target:

  target rotation; see
  [`summary-method`](https://philchalmers.github.io/mirt/reference/summary-method.md)
  for details

- response.pattern:

  an optional argument used to calculate the factor scores and standard
  errors for a given response vector or matrix/data.frame

- append_response.pattern:

  logical; should the inputs from `response.pattern` also be appended to
  the factor score output?

- na.rm:

  logical; remove rows with any missing values? This is generally not
  required due to the nature of computing factors scores, however for
  the "EAPsum" method this may be necessary to ensure that the
  sum-scores correspond to the same composite score

- plausible.draws:

  number of plausible values to draw for future researchers to perform
  secondary analyses of the latent trait scores. Typically used in
  conjunction with latent regression predictors (see
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  details), but can also be generated when no predictor variables were
  modelled. If `plausible.draws` is greater than 0 a list of plausible
  values will be returned

- plausible.type:

  type of plausible values to obtain. Can be either `'normal'` (default)
  to use a normal approximation based on the ACOV matrix, or `'MH'` to
  obtain Metropolis-Hastings samples from the posterior (silently passes
  object to
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
  therefore arguments like `technical` can be supplied to increase the
  number of burn-in draws and discarded samples)

- quadpts:

  number of quadrature to use per dimension. If not specified, a
  suitable one will be created which decreases as the number of
  dimensions increases (and therefore for estimates such as EAP, will be
  less accurate). This is determined from the switch statement
  `quadpts <- switch(as.character(nfact), '1'=121, '2'=61, '3'=31, '4'=19, '5'=11, '6'=7, 5)`

- item_weights:

  a user-defined weight vector used in the likelihood expressions to add
  more/less weight for a given observed response. Default is a vector of
  1's, indicating that all the items receive the same weight

- returnER:

  logical; return empirical reliability (also known as marginal
  reliability) estimates as a numeric values?

- T_as_X:

  logical; should the observed variance be equal to
  `var(X) = var(T) + E(E^2)` or `var(X) = var(T)` when computing
  empirical reliability estimates? Default (`FALSE`) uses the former

- EAPsum.scores:

  logical; include the model-implied expected values and variance for
  the item and total scores when using `method = 'EAPsum'` with
  `full.scores=FALSE`? This information is included in the hidden
  `'fit'` attribute which can be extracted via `attr(., 'fit')` for
  later use

- return.acov:

  logical; return a list containing covariance matrices instead of
  factors scores? `impute = TRUE` not supported with this option

- mean:

  a vector for custom latent variable means. If NULL, the default for
  'group' values from the computed mirt object will be used

- cov:

  a custom matrix of the latent variable covariance matrix. If NULL, the
  default for 'group' values from the computed mirt object will be used

- covdata:

  when latent regression model has been fitted, and the
  `response.pattern` input is used to score individuals, then this
  argument is used to include the latent regression covariate terms for
  each row vector supplied to `response.pattern`

- verbose:

  logical; print verbose output messages?

- full.scores.SE:

  logical; when `full.scores == TRUE`, also return the standard errors
  associated with each respondent? Default is `FALSE`

- theta_lim:

  lower and upper range to evaluate latent trait integral for each
  dimension. If omitted, a range will be generated automatically based
  on the number of dimensions

- MI:

  a number indicating how many multiple imputation draws to perform.
  Default is 0, indicating that no MI draws will be performed

- use_dentype_estimate:

  logical; if the density of the latent trait was estimated in the model
  (e.g., via Davidian curves or empirical histograms), should this
  information be used to compute the latent trait estimates? Only
  applicable for EAP-based estimates (EAP, EAPsum, and plausible)

- QMC:

  logical; use quasi-Monte Carlo integration? If `quadpts` is omitted
  the default number of nodes is 5000

- custom_den:

  a function used to define the integration density (if required). The
  NULL default assumes that the multivariate normal distribution with
  the 'GroupPars' hyper-parameters are used. At the minimum must be of
  the form:

  `function(Theta, ...)`

  where Theta is a matrix of latent trait values (will be a grid of
  values if `method == 'EAPsum'` or `method == 'EAP'`, otherwise Theta
  will have only 1 row). Additional arguments may included and are
  caught through the `fscores(...)` input. The function *must* return a
  numeric vector of density weights (one for each row in Theta)

- custom_theta:

  a matrix of custom integration nodes to use instead of the default,
  where each column corresponds to the respective dimension in the model

- expected.info:

  logical; instead of using the observed information when
  `method = 'ML'` or `method = 'WLE'` use the expected information
  computed via
  [`testinfo`](https://philchalmers.github.io/mirt/reference/testinfo.md)?
  Only currently supported for unidimensional models

- min_expected:

  when computing goodness of fit tests when `method = 'EAPsum'`, this
  value is used to collapse across the conditioned total scores until
  the expected values are greater than this value. Note that this only
  affect the goodness of fit tests and not the returned EAP for sum
  scores table

- max_theta:

  the maximum/minimum value any given factor score estimate will achieve
  using any modal estimator method (e.g., MAP, WLE, ML)

- start:

  a matrix of starting values to use for iterative estimation methods.
  Default will start at a vector of 0's for each response pattern, or
  will start at the EAP estimates (unidimensional models only). Must be
  in the form that matches `full.scores = FALSE` (mostly used in the
  `mirtCAT` package)

- ...:

  additional arguments to be passed to `nlm`

## Details

The function will return either a table with the computed scores and
standard errors, the original data matrix with scores appended to the
rightmost column, or the scores only. By default the latent means and
covariances are determined from the estimated object, though these can
be overwritten. Iterative estimation methods can be estimated in
parallel to decrease estimation times if a
[`mirtCluster`](https://philchalmers.github.io/mirt/reference/mirtCluster.md)
object is available.

If the input object is a discrete latent class object estimated from
[`mdirt`](https://philchalmers.github.io/mirt/reference/mdirt.md) then
the returned results will be with respect to the posterior
classification for each individual. The method inputs for
`'DiscreteClass'` objects may only be `'EAP'`, for posterior
classification of each response pattern, or `'EAPsum'` for posterior
classification based on the raw sum-score. For more information on these
algorithms refer to the `mirtCAT` package and the associated JSS paper
(Chalmers, 2016).

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Chalmers, R. P. (2016). Generating Adaptive and Non-Adaptive Test
Interfaces for Multidimensional Item Response Theory Applications.
*Journal of Statistical Software, 71*(5), 1-39.
[doi:10.18637/jss.v071.i05](https://doi.org/10.18637/jss.v071.i05)

Embretson, S. E. & Reise, S. P. (2000). Item Response Theory for
Psychologists. Erlbaum.

Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995).
Item Response Theory for Scores on Tests Including Polytomous Items with
Ordered Responses. *Applied Psychological Measurement, 19*, 39-49.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*, 427-450.

## See also

[`averageMI`](https://philchalmers.github.io/mirt/reference/averageMI.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
mod <- mirt(Science)
tabscores <- fscores(mod, full.scores = FALSE)
head(tabscores)
#>      Comfort Work Future Benefit         F1     SE_F1
#> [1,]       1    1      1       1 -2.7492669 0.6293525
#> [2,]       1    3      2       1 -1.4198318 0.5772364
#> [3,]       1    4      2       3 -0.7141976 0.6200139
#> [4,]       1    4      3       1 -0.4469265 0.6509531
#> [5,]       2    1      1       1 -2.5437807 0.5909114
#> [6,]       2    1      2       4 -1.2478570 0.5840105

# convert scores into expected total score information with 95% CIs
E.total <- expected.test(mod, Theta=tabscores[,'F1'])
E.total_2.5 <- expected.test(mod, Theta=tabscores[,'F1'] +
                                       tabscores[,'SE_F1'] * qnorm(.05/2))
E.total_97.5 <- expected.test(mod, Theta=tabscores[,'F1'] +
                                       tabscores[,'SE_F1'] * qnorm(1-.05/2))

data.frame(Total_score=rowSums(tabscores[,1:4]),
           E.total, E.total_2.5, E.total_97.5) |> head()
#>   Total_score   E.total E.total_2.5 E.total_97.5
#> 1           4  6.791606    5.321810     9.084082
#> 2           7  9.266018    7.128071    11.296189
#> 3          10 10.584682    8.296461    12.504975
#> 4           9 11.041648    8.691107    13.034195
#> 5           5  7.141179    5.576233     9.330947
#> 6           9  9.592533    7.415339    11.582060

# \donttest{
fullscores <- fscores(mod)
fullscores_with_SE <- fscores(mod, full.scores.SE=TRUE)
head(fullscores)
#>              F1
#> [1,]  0.4015613
#> [2,]  0.0520324
#> [3,] -0.8906436
#> [4,] -0.8906436
#> [5,]  0.7653806
#> [6,]  0.6695350
head(fullscores_with_SE)
#>              F1     SE_F1
#> [1,]  0.4015613 0.5978747
#> [2,]  0.0520324 0.5549554
#> [3,] -0.8906436 0.5421855
#> [4,] -0.8906436 0.5421855
#> [5,]  0.7653806 0.6385998
#> [6,]  0.6695350 0.5761860

# convert scores into expected total score information with 95% CIs
E.total <- expected.test(mod, Theta=fullscores[,'F1'])
E.total_2.5 <- expected.test(mod, Theta=fullscores_with_SE[,'F1'] +
                                 fullscores_with_SE[,'SE_F1'] * qnorm(.05/2))
E.total_97.5 <- expected.test(mod, Theta=fullscores_with_SE[,'F1'] +
                               fullscores_with_SE[,'SE_F1'] * qnorm(1-.05/2))

data.frame(Total_score=rowSums(Science),
           E.total, E.total_2.5, E.total_97.5) |> head()
#>   Total_score  E.total E.total_2.5 E.total_97.5
#> 1          13 12.34882   10.484432     14.15621
#> 2          12 11.81643    9.994252     13.53226
#> 3          10 10.26478    8.250620     11.99721
#> 4          10 10.26478    8.250620     11.99721
#> 5          12 12.93068   10.976668     14.67563
#> 6          14 12.77499   11.020512     14.43522

# change method argument to use MAP estimates
fullscores <- fscores(mod, method='MAP')
head(fullscores)
#>               F1
#> [1,]  0.42140793
#> [2,]  0.05866659
#> [3,] -0.91936052
#> [4,] -0.91936052
#> [5,]  0.79013883
#> [6,]  0.68302194

# calculate MAP for a given response vector
fscores(mod, method='MAP', response.pattern = c(1,2,3,4))
#>         F1 SE_F1
#> [1,] -0.37 0.605
# or matrix
fscores(mod, method='MAP', response.pattern = rbind(c(1,2,3,4), c(2,2,1,3)))
#>          F1 SE_F1
#> [1,] -0.370 0.605
#> [2,] -1.694 0.577

# return only the scores and their SEs
fscores(mod, method='MAP', response.pattern = c(1,2,3,4))
#>         F1 SE_F1
#> [1,] -0.37 0.605

# use custom latent variable properties (diffuse prior for MAP is very close to ML)
fscores(mod, method='MAP', cov = matrix(1000), full.scores = FALSE)
#>       Comfort Work Future Benefit     F1  SE_F1
#>  [1,]       1    1      1       1 -9.340  9.636
#>  [2,]       1    3      2       1 -2.104  0.689
#>  [3,]       1    4      2       3 -1.151  0.705
#>  [4,]       1    4      3       1 -0.781  0.807
#>  [5,]       2    1      1       1 -4.394  1.179
#>  [6,]       2    1      2       4 -1.863  0.668
#>  [7,]       2    2      1       1 -3.248  0.812
#>  [8,]       2    2      2       2 -1.911  0.588
#>  [9,]       2    2      2       3 -1.606  0.608
#> [10,]       2    2      3       1 -1.468  0.714
#> [11,]       2    2      3       2 -1.168  0.635
#> [12,]       2    2      3       3 -0.797  0.639
#> [13,]       2    2      4       3  0.240  0.789
#> [14,]       2    3      1       3 -2.142  0.702
#> [15,]       2    3      2       2 -1.609  0.626
#> [16,]       2    3      2       3 -1.254  0.630
#> [17,]       2    3      3       2 -0.756  0.656
#> [18,]       2    3      3       3 -0.332  0.698
#> [19,]       2    3      3       4  0.083  0.761
#> [20,]       2    3      4       1  0.189  0.878
#> [21,]       2    3      4       3  0.763  0.681
#> [22,]       2    4      2       1 -1.822  0.694
#> [23,]       2    4      4       3  1.307  0.753
#> [24,]       2    4      4       4  2.095  1.039
#> [25,]       3    1      1       1 -3.538  1.001
#> [26,]       3    1      1       3 -2.510  0.697
#> [27,]       3    1      2       2 -1.926  0.615
#> [28,]       3    1      2       3 -1.587  0.643
#> [29,]       3    1      3       2 -1.088  0.674
#> [30,]       3    1      3       3 -0.658  0.696
#> [31,]       3    1      3       4 -0.281  0.807
#> [32,]       3    1      4       3  0.542  0.734
#> [33,]       3    1      4       4  1.038  0.736
#> [34,]       3    2      1       2 -2.348  0.625
#> [35,]       3    2      1       4 -1.913  0.709
#> [36,]       3    2      2       1 -1.861  0.616
#> [37,]       3    2      2       2 -1.577  0.594
#> [38,]       3    2      2       3 -1.255  0.602
#> [39,]       3    2      3       1 -1.025  0.663
#> [40,]       3    2      3       2 -0.798  0.629
#> [41,]       3    2      3       3 -0.409  0.668
#> [42,]       3    2      3       4 -0.037  0.744
#> [43,]       3    2      4       1  0.015  0.900
#> [44,]       3    2      4       2  0.205  0.789
#> [45,]       3    2      4       3  0.643  0.684
#> [46,]       3    2      4       4  1.102  0.718
#> [47,]       3    3      1       3 -1.629  0.780
#> [48,]       3    3      2       1 -1.518  0.659
#> [49,]       3    3      2       2 -1.239  0.617
#> [50,]       3    3      2       3 -0.882  0.630
#> [51,]       3    3      2       4 -0.621  0.708
#> [52,]       3    3      3       1 -0.549  0.711
#> [53,]       3    3      3       2 -0.347  0.687
#> [54,]       3    3      3       3  0.086  0.687
#> [55,]       3    3      3       4  0.490  0.676
#> [56,]       3    3      4       2  0.732  0.681
#> [57,]       3    3      4       3  1.048  0.650
#> [58,]       3    3      4       4  1.525  0.746
#> [59,]       3    4      1       3 -1.381  0.943
#> [60,]       3    4      2       3 -0.608  0.726
#> [61,]       3    4      3       2  0.099  0.768
#> [62,]       3    4      3       3  0.538  0.677
#> [63,]       3    4      3       4  0.961  0.675
#> [64,]       3    4      4       1  1.218  0.758
#> [65,]       3    4      4       2  1.272  0.749
#> [66,]       3    4      4       3  1.599  0.766
#> [67,]       3    4      4       4  2.422  1.040
#> [68,]       4    1      1       4 -2.106  0.772
#> [69,]       4    1      2       2 -1.638  0.654
#> [70,]       4    1      2       4 -1.014  0.724
#> [71,]       4    1      3       4  0.361  0.778
#> [72,]       4    2      2       1 -1.568  0.651
#> [73,]       4    2      2       3 -0.935  0.628
#> [74,]       4    2      3       3  0.054  0.711
#> [75,]       4    2      3       4  0.492  0.710
#> [76,]       4    2      4       2  0.762  0.723
#> [77,]       4    2      4       4  1.760  0.931
#> [78,]       4    3      2       3 -0.485  0.704
#> [79,]       4    3      3       2  0.139  0.718
#> [80,]       4    3      3       3  0.533  0.652
#> [81,]       4    3      3       4  0.929  0.655
#> [82,]       4    3      4       2  1.216  0.719
#> [83,]       4    3      4       3  1.527  0.734
#> [84,]       4    3      4       4  2.258  0.968
#> [85,]       4    4      3       2  0.640  0.707
#> [86,]       4    4      3       3  0.981  0.661
#> [87,]       4    4      3       4  1.465  0.749
#> [88,]       4    4      4       2  2.028  1.010
#> [89,]       4    4      4       3  2.390  1.015
#> [90,]       4    4      4       4  7.112 10.616
fscores(mod, method='ML', full.scores = FALSE)
#>       Comfort Work Future Benefit     F1 SE_F1
#>  [1,]       1    1      1       1   -Inf    NA
#>  [2,]       1    3      2       1 -2.105 0.689
#>  [3,]       1    4      2       3 -1.151 0.705
#>  [4,]       1    4      3       1 -0.781 0.807
#>  [5,]       2    1      1       1 -4.400 1.182
#>  [6,]       2    1      2       4 -1.864 0.668
#>  [7,]       2    2      1       1 -3.250 0.813
#>  [8,]       2    2      2       2 -1.912 0.588
#>  [9,]       2    2      2       3 -1.606 0.609
#> [10,]       2    2      3       1 -1.469 0.714
#> [11,]       2    2      3       2 -1.169 0.635
#> [12,]       2    2      3       3 -0.797 0.639
#> [13,]       2    2      4       3  0.240 0.789
#> [14,]       2    3      1       3 -2.143 0.702
#> [15,]       2    3      2       2 -1.609 0.626
#> [16,]       2    3      2       3 -1.254 0.630
#> [17,]       2    3      3       2 -0.756 0.656
#> [18,]       2    3      3       3 -0.332 0.698
#> [19,]       2    3      3       4  0.083 0.762
#> [20,]       2    3      4       1  0.189 0.878
#> [21,]       2    3      4       3  0.763 0.681
#> [22,]       2    4      2       1 -1.823 0.694
#> [23,]       2    4      4       3  1.308 0.754
#> [24,]       2    4      4       4  2.097 1.041
#> [25,]       3    1      1       1 -3.542 1.004
#> [26,]       3    1      1       3 -2.511 0.698
#> [27,]       3    1      2       2 -1.926 0.615
#> [28,]       3    1      2       3 -1.588 0.644
#> [29,]       3    1      3       2 -1.088 0.674
#> [30,]       3    1      3       3 -0.658 0.696
#> [31,]       3    1      3       4 -0.281 0.808
#> [32,]       3    1      4       3  0.543 0.734
#> [33,]       3    1      4       4  1.038 0.736
#> [34,]       3    2      1       2 -2.349 0.625
#> [35,]       3    2      1       4 -1.914 0.709
#> [36,]       3    2      2       1 -1.862 0.616
#> [37,]       3    2      2       2 -1.577 0.594
#> [38,]       3    2      2       3 -1.255 0.602
#> [39,]       3    2      3       1 -1.026 0.663
#> [40,]       3    2      3       2 -0.798 0.629
#> [41,]       3    2      3       3 -0.409 0.668
#> [42,]       3    2      3       4 -0.037 0.744
#> [43,]       3    2      4       1  0.015 0.901
#> [44,]       3    2      4       2  0.206 0.789
#> [45,]       3    2      4       3  0.643 0.684
#> [46,]       3    2      4       4  1.102 0.719
#> [47,]       3    3      1       3 -1.630 0.780
#> [48,]       3    3      2       1 -1.518 0.659
#> [49,]       3    3      2       2 -1.240 0.617
#> [50,]       3    3      2       3 -0.883 0.630
#> [51,]       3    3      2       4 -0.621 0.708
#> [52,]       3    3      3       1 -0.549 0.711
#> [53,]       3    3      3       2 -0.348 0.688
#> [54,]       3    3      3       3  0.086 0.687
#> [55,]       3    3      3       4  0.490 0.676
#> [56,]       3    3      4       2  0.733 0.681
#> [57,]       3    3      4       3  1.048 0.650
#> [58,]       3    3      4       4  1.526 0.746
#> [59,]       3    4      1       3 -1.383 0.943
#> [60,]       3    4      2       3 -0.609 0.726
#> [61,]       3    4      3       2  0.099 0.768
#> [62,]       3    4      3       3  0.538 0.678
#> [63,]       3    4      3       4  0.962 0.675
#> [64,]       3    4      4       1  1.219 0.759
#> [65,]       3    4      4       2  1.273 0.749
#> [66,]       3    4      4       3  1.600 0.767
#> [67,]       3    4      4       4  2.425 1.042
#> [68,]       4    1      1       4 -2.107 0.772
#> [69,]       4    1      2       2 -1.639 0.654
#> [70,]       4    1      2       4 -1.015 0.724
#> [71,]       4    1      3       4  0.361 0.778
#> [72,]       4    2      2       1 -1.568 0.651
#> [73,]       4    2      2       3 -0.936 0.628
#> [74,]       4    2      3       3  0.054 0.711
#> [75,]       4    2      3       4  0.492 0.710
#> [76,]       4    2      4       2  0.763 0.723
#> [77,]       4    2      4       4  1.761 0.932
#> [78,]       4    3      2       3 -0.485 0.704
#> [79,]       4    3      3       2  0.139 0.718
#> [80,]       4    3      3       3  0.533 0.652
#> [81,]       4    3      3       4  0.929 0.656
#> [82,]       4    3      4       2  1.216 0.720
#> [83,]       4    3      4       3  1.527 0.735
#> [84,]       4    3      4       4  2.261 0.969
#> [85,]       4    4      3       2  0.640 0.707
#> [86,]       4    4      3       3  0.982 0.662
#> [87,]       4    4      3       4  1.465 0.749
#> [88,]       4    4      4       2  2.030 1.011
#> [89,]       4    4      4       3  2.392 1.016
#> [90,]       4    4      4       4    Inf    NA

# EAPsum table of values based on total scores
(fs <- fscores(mod, method = 'EAPsum', full.scores = FALSE))
#>    Sum.Scores          F1     SE_F1 observed   expected   std.res
#> 4           4 -2.74926694 0.6293525        2  0.1241448 5.3239625
#> 5           5 -2.43096701 0.6167666        1  0.7902670 0.2359282
#> 6           6 -2.08095804 0.6101826        2  2.7598787 0.4574033
#> 7           7 -1.71797996 0.6019326        1  7.2515944 2.3215286
#> 8           8 -1.36354348 0.5979810       11 15.8925644 1.2272684
#> 9           9 -1.01186764 0.6038596       32 29.6740255 0.4269890
#> 10         10 -0.64926712 0.6101301       58 48.3632492 1.3857117
#> 11         11 -0.28688279 0.6048405       70 68.4845962 0.1831184
#> 12         12  0.08240446 0.5996636       91 80.5530281 1.1639907
#> 13         13  0.48709473 0.6133534       56 65.6680929 1.1930636
#> 14         14  0.93406958 0.6170256       36 42.6372130 1.0164625
#> 15         15  1.38421709 0.6218261       20 22.3313442 0.4933430
#> 16         16  1.85351389 0.6544097       12  7.4700015 1.6574396

# convert expected counts back into marginal probability distribution
within(fs,
   `P(y)` <- expected / sum(observed))
#>    Sum.Scores          F1     SE_F1 observed   expected   std.res        P(y)
#> 4           4 -2.74926694 0.6293525        2  0.1241448 5.3239625 0.000316696
#> 5           5 -2.43096701 0.6167666        1  0.7902670 0.2359282 0.002015987
#> 6           6 -2.08095804 0.6101826        2  2.7598787 0.4574033 0.007040507
#> 7           7 -1.71797996 0.6019326        1  7.2515944 2.3215286 0.018498965
#> 8           8 -1.36354348 0.5979810       11 15.8925644 1.2272684 0.040542256
#> 9           9 -1.01186764 0.6038596       32 29.6740255 0.4269890 0.075699045
#> 10         10 -0.64926712 0.6101301       58 48.3632492 1.3857117 0.123375636
#> 11         11 -0.28688279 0.6048405       70 68.4845962 0.1831184 0.174705603
#> 12         12  0.08240446 0.5996636       91 80.5530281 1.1639907 0.205492419
#> 13         13  0.48709473 0.6133534       56 65.6680929 1.1930636 0.167520645
#> 14         14  0.93406958 0.6170256       36 42.6372130 1.0164625 0.108768401
#> 15         15  1.38421709 0.6218261       20 22.3313442 0.4933430 0.056967715
#> 16         16  1.85351389 0.6544097       12  7.4700015 1.6574396 0.019056126

# list of error VCOV matrices for EAPsum (works for other estimators as well)
acovs <- fscores(mod, method = 'EAPsum', full.scores = FALSE, return.acov = TRUE)
acovs
#> $`0`
#>           [,1]
#> [1,] 0.3960846
#> 
#> $`1`
#>           [,1]
#> [1,] 0.3804011
#> 
#> $`2`
#>           [,1]
#> [1,] 0.3723228
#> 
#> $`3`
#>           [,1]
#> [1,] 0.3623228
#> 
#> $`4`
#>           [,1]
#> [1,] 0.3575813
#> 
#> $`5`
#>           [,1]
#> [1,] 0.3646464
#> 
#> $`6`
#>           [,1]
#> [1,] 0.3722588
#> 
#> $`7`
#>          [,1]
#> [1,] 0.365832
#> 
#> $`8`
#>           [,1]
#> [1,] 0.3595964
#> 
#> $`9`
#>           [,1]
#> [1,] 0.3762024
#> 
#> $`10`
#>           [,1]
#> [1,] 0.3807206
#> 
#> $`11`
#>           [,1]
#> [1,] 0.3866677
#> 
#> $`12`
#>           [,1]
#> [1,] 0.4282521
#> 

# WLE estimation, run in parallel using available cores
if(interactive()) mirtCluster()
head(fscores(mod, method='WLE', full.scores = FALSE))
#>      Comfort Work Future Benefit         F1     SE_F1
#> [1,]       1    1      1       1 -5.6980733 1.5782658
#> [2,]       1    3      2       1 -2.1191038 0.6332801
#> [3,]       1    4      2       3 -1.1387624 0.6557002
#> [4,]       1    4      3       1 -0.8489387 0.7000117
#> [5,]       2    1      1       1 -4.0112458 1.1423934
#> [6,]       2    1      2       4 -1.8957020 0.6698434

# multiple imputation using 30 draws for EAP scores. Requires information matrix
mod <- mirt(Science, 1, SE=TRUE)
fs <- fscores(mod, MI = 30)
head(fs)
#>               F1
#> [1,]  0.40337517
#> [2,]  0.05309148
#> [3,] -0.87541372
#> [4,] -0.87541372
#> [5,]  0.73035366
#> [6,]  0.68055257

# plausible values for future work
pv <- fscores(mod, plausible.draws = 5)
lapply(pv, function(x) c(mean=mean(x), var=var(x), min=min(x), max=max(x)))
#> [[1]]
#>         mean          var          min          max 
#>  0.008078887  0.954639097 -2.990804153  2.714129271 
#> 
#> [[2]]
#>        mean         var         min         max 
#>  0.01798039  0.93077370 -3.15510927  3.24355926 
#> 
#> [[3]]
#>        mean         var         min         max 
#>  0.05404685  0.98508826 -3.67842686  2.78694047 
#> 
#> [[4]]
#>         mean          var          min          max 
#> -0.009914966  0.962301671 -3.063918906  2.859407646 
#> 
#> [[5]]
#>        mean         var         min         max 
#> -0.01718763  0.88897331 -3.45397245  2.68883001 
#> 

## define a custom_den function (*must* return a numeric vector).
#  EAP with a uniform prior between -3 and 3
fun <- function(Theta, ...) as.numeric(dunif(Theta, min = -3, max = 3))
head(fscores(mod, custom_den = fun))
#>               F1
#> [1,]  0.62740797
#> [2,]  0.07362014
#> [3,] -1.23412036
#> [4,] -1.23412036
#> [5,]  1.25196942
#> [6,]  1.00691502

# compare EAP estimators with same modified prior
fun <- function(Theta, ...) as.numeric(dnorm(Theta, mean=.5))
head(fscores(mod, custom_den = fun))
#>              F1
#> [1,]  0.5801169
#> [2,]  0.2057691
#> [3,] -0.7411082
#> [4,] -0.7411082
#> [5,]  0.9677165
#> [6,]  0.8355336
head(fscores(mod, method = 'EAP', mean=.5))
#>              F1
#> [1,]  0.5801169
#> [2,]  0.2057691
#> [3,] -0.7411082
#> [4,] -0.7411082
#> [5,]  0.9677165
#> [6,]  0.8355336

if (FALSE) { # \dontrun{
# custom MAP prior: standard truncated normal between 5 and -2
library(msm)
# need the :: scope for parallel to see the function (not require if no mirtCluster() defined)
fun <- function(Theta, ...) msm::dtnorm(Theta, mean = 0, sd = 1, lower = -2, upper = 5)
head(fscores(mod, custom_den = fun, method = 'MAP', full.scores = FALSE))
} # }


####################
# scoring via response.pattern input (with latent regression structure)
# simulate data
set.seed(1234)
N <- 1000

# covariates
X1 <- rnorm(N); X2 <- rnorm(N)
covdata <- data.frame(X1, X2)
Theta <- matrix(0.5 * X1 + -1 * X2 + rnorm(N, sd = 0.5))

# items and response data
a <- matrix(1, 20); d <- matrix(rnorm(20))
dat <- simdata(a, d, 1000, itemtype = '2PL', Theta=Theta)

# conditional model using X1 and X2 as predictors of Theta
mod <- mirt(dat, 1, 'Rasch', covdata=covdata, formula = ~ X1 + X2)
coef(mod, simplify=TRUE)
#> $items
#>         a1      d g u
#> Item_1   1 -0.409 0 1
#> Item_2   1  0.491 0 1
#> Item_3   1  0.313 0 1
#> Item_4   1  1.965 0 1
#> Item_5   1  1.753 0 1
#> Item_6   1 -0.246 0 1
#> Item_7   1 -1.077 0 1
#> Item_8   1  0.533 0 1
#> Item_9   1 -1.232 0 1
#> Item_10  1  0.603 0 1
#> Item_11  1 -0.404 0 1
#> Item_12  1  1.238 0 1
#> Item_13  1  1.033 0 1
#> Item_14  1  1.524 0 1
#> Item_15  1 -0.548 0 1
#> Item_16  1  2.075 0 1
#> Item_17  1 -0.695 0 1
#> Item_18  1 -1.200 0 1
#> Item_19  1  0.121 0 1
#> Item_20  1  0.523 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>       F1
#> F1 0.215
#> 
#> $lr.betas
#>                 F1
#> (Intercept)  0.000
#> X1           0.527
#> X2          -1.036
#> 

# all EAP estimates that include latent regression information
fs <- fscores(mod, full.scores.SE=TRUE)
head(fs)
#>              F1     SE_F1
#> [1,]  0.3085492 0.3442152
#> [2,] -0.3474178 0.3408490
#> [3,]  2.0484395 0.3903304
#> [4,] -1.9418498 0.3670495
#> [5,] -0.7631517 0.3432688
#> [6,]  2.1019471 0.3922103

# score only two response patterns
rp <- dat[1:2, ]
cd <- covdata[1:2, ]

fscores(mod, response.pattern=rp, covdata=cd)
#>          F1 SE_F1
#> [1,]  0.309 0.344
#> [2,] -0.347 0.341
fscores(mod, response.pattern=rp[2,], covdata=cd[2,]) # just one pattern
#>          F1 SE_F1
#> [1,] -0.347 0.341

# }
```
