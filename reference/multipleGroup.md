# Multiple Group Estimation

`multipleGroup` performs a full-information maximum-likelihood multiple
group analysis for any combination of dichotomous and polytomous data
under the item response theory paradigm using either Cai's (2010)
Metropolis-Hastings Robbins-Monro (MHRM) algorithm or with an EM
algorithm approach. This function may be used for detecting differential
item functioning (DIF), thought the
[`DIF`](https://philchalmers.github.io/mirt/reference/DIF.md) function
may provide a more convenient approach. If the grouping variable is not
specified then the `dentype` input can be modified to fit mixture models
to estimate any latent group components.

## Usage

``` r
multipleGroup(
  data,
  model = 1,
  group,
  itemtype = NULL,
  invariance = "",
  method = "EM",
  dentype = "Gaussian",
  itemdesign = NULL,
  item.formula = NULL,
  nruns = 1,
  return_max = TRUE,
  GenRandomPars = FALSE,
  verbose = interactive(),
  ...
)
```

## Arguments

- data:

  a `matrix` or `data.frame` that consists of numerically ordered data,
  organized in the form of integers, with missing data coded as `NA`

- model:

  string to be passed to, or a model object returned from,
  [`mirt.model`](https://philchalmers.github.io/mirt/reference/mirt.model.md)
  declaring how the global model is to be estimated (useful to apply
  constraints here)

- group:

  a `character` or `factor` vector indicating group membership. If a
  `character` vector is supplied this will be automatically transformed
  into a [`factor`](https://rdrr.io/r/base/factor.html) variable. As
  well, the first level of the (factorized) grouping variable will be
  treated as the "reference" group

- itemtype:

  can be same type of input as is documented in
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
  however may also be a `ngroups` by `nitems` matrix specifying the type
  of IRT models for each group, respectively. Rows of this input
  correspond to the levels of the `group` input. For mixture models the
  rows correspond to the respective mixture grouping variables to be
  constructed, and the IRT models should be within these mixtures

- invariance:

  a character vector containing the following possible options:

  `'free_mean'` or `'free_means'`

  :   freely estimate all latent means in all focal groups (reference
      group constrained to a vector of 0's)

  `'free_var'`, `'free_vars'`, `'free_variance'`, or `'free_variances'`

  :   freely estimate all latent variances in focal groups (reference
      group variances all constrained to 1)

  `'slopes'`

  :   to constrain all the slopes to be equal across all groups

  `'intercepts'`

  :   to constrain all the intercepts to be equal across all groups,
      note for nominal models this also includes the category specific
      slope parameters

  Additionally, specifying specific item name bundles (from
  `colnames(data)`) will constrain all freely estimated parameters in
  each item to be equal across groups. This is useful for selecting
  'anchor' items for vertical and horizontal scaling, and for detecting
  differential item functioning (DIF) across groups

- method:

  a character object that is either `'EM'`, `'QMCEM'`, or `'MHRM'`
  (default is `'EM'`). See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  details

- dentype:

  type of density form to use for the latent trait parameters. Current
  options include all of the methods described in
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md), as
  well as

  - `'mixture-#'` estimates mixtures of Gaussian distributions, where
    the `#` placeholder represents the number of potential grouping
    variables (e.g., `'mixture-3'` will estimate 3 underlying classes).
    Each class is assigned the group name `MIXTURE_#`, where `#` is the
    class number.

    Note that internally the mixture coefficients are stored as log
    values where the first mixture group coefficient is fixed at 0.
    Additionally, it is recommended to use the `nruns` argument as
    mixture IRT models are known to contain local maximums

- itemdesign:

  see [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)
  for details

- item.formula:

  see [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)
  for details

- nruns:

  a numeric value indicating how many times the model should be fit to
  the data when using random starting values, which is particularly
  useful when evaluating mixture IRT Models. If greater than 1,
  `GenRandomPars` is set to `TRUE` by default. Using this returns a list
  of fitted model objects, where the model with the highest
  log-likelihood should generally be selected as the model best
  associated with the MLE (this is done automatically if
  `return_max = TRUE`). Note that if a
  [`mirtCluster`](https://philchalmers.github.io/mirt/reference/mirtCluster.md)
  was defined earlier then the runs will be run in parallel

- return_max:

  logical; when `nruns > 1`, return the model that has the most optimal
  maximum likelihood criteria? If FALSE, returns a list of all the
  estimated objects

- GenRandomPars:

  see [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)
  for details

- verbose:

  see [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)
  for details

- ...:

  additional arguments to be passed to the estimation engine. See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  details and examples

## Value

function returns an object of class `MultipleGroupClass`
([MultipleGroupClass-class](https://philchalmers.github.io/mirt/reference/MultipleGroupClass-class.md)).

## Details

By default the estimation in `multipleGroup` assumes that the models are
maximally independent, and therefore could initially be performed by
sub-setting the data and running identical models with `mirt` and
aggregating the results (e.g., log-likelihood). However, constrains may
be automatically imposed across groups by invoking various `invariance`
keywords. Users may also supply a list of parameter equality constraints
to by `constrain` argument, of define equality constraints using the
[`mirt.model`](https://philchalmers.github.io/mirt/reference/mirt.model.md)
syntax (recommended).

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Magnus, B. E. and Garnier-Villarreal (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, 27, 261-279.

Wall, M., M., Park, J., Y., and Moustaki I. (2015). IRT modeling in the
presence of zero-inflation with application to psychiatric disorder
severity. *Applied Psychological Measurement* 39: 583-597.

## See also

[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
[`DIF`](https://philchalmers.github.io/mirt/reference/DIF.md),
[`extract.group`](https://philchalmers.github.io/mirt/reference/extract.group.md),
[`DRF`](https://philchalmers.github.io/mirt/reference/DRF.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

# single factor
set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
itemtype <- rep('2PL', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))

# marginal information
itemstats(dat)
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  2000            7.888          3.555 0.188 0.053 0.777     1.678
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  2000 2 0.609 0.488   0.525         0.414       0.762
#> Item_2  2000 2 0.392 0.488   0.540         0.431       0.761
#> Item_3  2000 2 0.456 0.498   0.526         0.413       0.762
#> Item_4  2000 2 0.684 0.465   0.461         0.349       0.768
#> Item_5  2000 2 0.538 0.499   0.528         0.415       0.762
#> Item_6  2000 2 0.649 0.477   0.334         0.207       0.779
#> Item_7  2000 2 0.681 0.466   0.524         0.419       0.762
#> Item_8  2000 2 0.432 0.495   0.504         0.389       0.764
#> Item_9  2000 2 0.302 0.459   0.446         0.334       0.769
#> Item_10 2000 2 0.274 0.446   0.388         0.273       0.774
#> Item_11 2000 2 0.734 0.442   0.457         0.351       0.768
#> Item_12 2000 2 0.470 0.499   0.585         0.481       0.756
#> Item_13 2000 2 0.585 0.493   0.561         0.455       0.759
#> Item_14 2000 2 0.592 0.492   0.538         0.428       0.761
#> Item_15 2000 2 0.490 0.500   0.458         0.336       0.769
#> 
#> $proportions
#>             0     1
#> Item_1  0.392 0.609
#> Item_2  0.608 0.392
#> Item_3  0.544 0.456
#> Item_4  0.316 0.684
#> Item_5  0.462 0.538
#> Item_6  0.351 0.649
#> Item_7  0.318 0.681
#> Item_8  0.569 0.432
#> Item_9  0.698 0.302
#> Item_10 0.726 0.274
#> Item_11 0.266 0.734
#> Item_12 0.530 0.470
#> Item_13 0.416 0.585
#> Item_14 0.408 0.592
#> Item_15 0.509 0.490
#> 

# conditional information
itemstats(dat, group=group)
#> $D1
#> $D1$overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  1000             7.82          3.346 0.159 0.047  0.74     1.705
#> 
#> $D1$itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  1000 2 0.605 0.489   0.484         0.361       0.725
#> Item_2  1000 2 0.368 0.483   0.507         0.388       0.722
#> Item_3  1000 2 0.461 0.499   0.471         0.343       0.727
#> Item_4  1000 2 0.673 0.469   0.439         0.315       0.729
#> Item_5  1000 2 0.526 0.500   0.500         0.376       0.723
#> Item_6  1000 2 0.654 0.476   0.355         0.222       0.739
#> Item_7  1000 2 0.683 0.466   0.508         0.394       0.722
#> Item_8  1000 2 0.431 0.495   0.462         0.333       0.728
#> Item_9  1000 2 0.287 0.453   0.419         0.298       0.731
#> Item_10 1000 2 0.274 0.446   0.372         0.249       0.735
#> Item_11 1000 2 0.739 0.439   0.401         0.282       0.732
#> Item_12 1000 2 0.456 0.498   0.563         0.448       0.715
#> Item_13 1000 2 0.584 0.493   0.535         0.416       0.719
#> Item_14 1000 2 0.592 0.492   0.483         0.358       0.725
#> Item_15 1000 2 0.487 0.500   0.451         0.321       0.729
#> 
#> $D1$proportions
#>             0     1
#> Item_1  0.395 0.605
#> Item_2  0.632 0.368
#> Item_3  0.539 0.461
#> Item_4  0.327 0.673
#> Item_5  0.474 0.526
#> Item_6  0.346 0.654
#> Item_7  0.317 0.683
#> Item_8  0.569 0.431
#> Item_9  0.713 0.287
#> Item_10 0.726 0.274
#> Item_11 0.261 0.739
#> Item_12 0.544 0.456
#> Item_13 0.416 0.584
#> Item_14 0.408 0.592
#> Item_15 0.513 0.487
#> 
#> 
#> $D2
#> $D2$overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  1000            7.955          3.754 0.217 0.065 0.807      1.65
#> 
#> $D2$itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  1000 2 0.612 0.488   0.563         0.464       0.792
#> Item_2  1000 2 0.417 0.493   0.569         0.469       0.792
#> Item_3  1000 2 0.450 0.498   0.578         0.479       0.791
#> Item_4  1000 2 0.695 0.461   0.484         0.381       0.798
#> Item_5  1000 2 0.551 0.498   0.553         0.451       0.793
#> Item_6  1000 2 0.644 0.479   0.316         0.195       0.811
#> Item_7  1000 2 0.680 0.467   0.540         0.443       0.794
#> Item_8  1000 2 0.432 0.496   0.544         0.440       0.794
#> Item_9  1000 2 0.317 0.466   0.470         0.365       0.799
#> Item_10 1000 2 0.275 0.447   0.403         0.297       0.804
#> Item_11 1000 2 0.729 0.445   0.508         0.413       0.796
#> Item_12 1000 2 0.483 0.500   0.606         0.511       0.788
#> Item_13 1000 2 0.585 0.493   0.587         0.491       0.790
#> Item_14 1000 2 0.591 0.492   0.589         0.493       0.790
#> Item_15 1000 2 0.494 0.500   0.466         0.351       0.801
#> 
#> $D2$proportions
#>             0     1
#> Item_1  0.388 0.612
#> Item_2  0.583 0.417
#> Item_3  0.550 0.450
#> Item_4  0.305 0.695
#> Item_5  0.449 0.551
#> Item_6  0.356 0.644
#> Item_7  0.320 0.680
#> Item_8  0.568 0.432
#> Item_9  0.683 0.317
#> Item_10 0.725 0.275
#> Item_11 0.271 0.729
#> Item_12 0.517 0.483
#> Item_13 0.415 0.585
#> Item_14 0.409 0.591
#> Item_15 0.506 0.494
#> 
#> 

mod_configural <- multipleGroup(dat, 1, group = group) #completely separate analyses
# limited information fit statistics
M2(mod_configural)
#>            M2  df     p RMSEA RMSEA_5 RMSEA_95 SRMSR.D1 SRMSR.D2   TLI CFI
#> stats 142.987 180 0.981     0       0        0    0.024    0.019 1.005   1

mod_metric <- multipleGroup(dat, 1, group = group, invariance=c('slopes')) #equal slopes
# equal intercepts, free variance and means
mod_scalar2 <- multipleGroup(dat, 1, group = group,
                             invariance=c('slopes', 'intercepts', 'free_var','free_means'))
mod_scalar1 <- multipleGroup(dat, 1, group = group,  #fixed means
                             invariance=c('slopes', 'intercepts', 'free_var'))
mod_fullconstrain <- multipleGroup(dat, 1, group = group,
                             invariance=c('slopes', 'intercepts'))
extract.mirt(mod_fullconstrain, 'time') #time of estimation components
#> TOTAL:   Data  Estep  Mstep     SE   Post 
#>  0.282  0.054  0.062  0.148  0.000  0.000 

# optionally use Newton-Raphson for (generally) faster convergence in the
#  M-step's, though occasionally less stable
mod_fullconstrain <- multipleGroup(dat, 1, group = group, optimizer = 'NR',
                             invariance=c('slopes', 'intercepts'))
extract.mirt(mod_fullconstrain, 'time') #time of estimation components
#> TOTAL:   Data  Estep  Mstep     SE   Post 
#>  0.183  0.053  0.070  0.043  0.000  0.000 

summary(mod_scalar2)
#> 
#> ----------
#> GROUP: D1 
#>          F1    h2
#>  [1,] 0.544 0.296
#>  [2,] 0.577 0.332
#>  [3,] 0.529 0.280
#>  [4,] 0.460 0.212
#>  [5,] 0.536 0.288
#>  [6,] 0.253 0.064
#>  [7,] 0.570 0.325
#>  [8,] 0.498 0.248
#>  [9,] 0.459 0.210
#> [10,] 0.373 0.139
#> [11,] 0.487 0.237
#> [12,] 0.632 0.400
#> [13,] 0.596 0.355
#> [14,] 0.561 0.315
#> [15,] 0.414 0.172
#> 
#> SS loadings:  3.872 
#> Proportion Var:  0.258 
#> 
#> Factor correlations: 
#> 
#>    F1
#> F1  1
#> 
#> ----------
#> GROUP: D2 
#>          F1    h2
#>  [1,] 0.634 0.402
#>  [2,] 0.665 0.443
#>  [3,] 0.619 0.383
#>  [4,] 0.548 0.300
#>  [5,] 0.626 0.392
#>  [6,] 0.313 0.098
#>  [7,] 0.659 0.434
#>  [8,] 0.587 0.345
#>  [9,] 0.546 0.298
#> [10,] 0.453 0.205
#> [11,] 0.575 0.331
#> [12,] 0.718 0.515
#> [13,] 0.683 0.467
#> [14,] 0.650 0.423
#> [15,] 0.498 0.248
#> 
#> SS loadings:  5.284 
#> Proportion Var:  0.352 
#> 
#> Factor correlations: 
#> 
#>    F1
#> F1  1
coef(mod_scalar2, simplify=TRUE)
#> $D1
#> $items
#>            a1      d g u
#> Item_1  1.104  0.538 0 1
#> Item_2  1.201 -0.630 0 1
#> Item_3  1.061 -0.265 0 1
#> Item_4  0.882  0.900 0 1
#> Item_5  1.082  0.164 0 1
#> Item_6  0.445  0.636 0 1
#> Item_7  1.180  0.976 0 1
#> Item_8  0.977 -0.377 0 1
#> Item_9  0.879 -1.035 0 1
#> Item_10 0.685 -1.118 0 1
#> Item_11 0.948  1.213 0 1
#> Item_12 1.389 -0.224 0 1
#> Item_13 1.262  0.429 0 1
#> Item_14 1.153  0.453 0 1
#> Item_15 0.774 -0.070 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> 
#> $D2
#> $items
#>            a1      d g u
#> Item_1  1.104  0.538 0 1
#> Item_2  1.201 -0.630 0 1
#> Item_3  1.061 -0.265 0 1
#> Item_4  0.882  0.900 0 1
#> Item_5  1.082  0.164 0 1
#> Item_6  0.445  0.636 0 1
#> Item_7  1.180  0.976 0 1
#> Item_8  0.977 -0.377 0 1
#> Item_9  0.879 -1.035 0 1
#> Item_10 0.685 -1.118 0 1
#> Item_11 0.948  1.213 0 1
#> Item_12 1.389 -0.224 0 1
#> Item_13 1.262  0.429 0 1
#> Item_14 1.153  0.453 0 1
#> Item_15 0.774 -0.070 0 1
#> 
#> $means
#>    F1 
#> 0.066 
#> 
#> $cov
#>       F1
#> F1 1.595
#> 
#> 
residuals(mod_scalar2)
#> $D1
#>         Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> Item_1      NA -0.036 -0.056  0.036 -0.014  0.014  0.017 -0.020 -0.008  -0.066
#> Item_2   1.305     NA -0.052 -0.032 -0.031  0.033 -0.028  0.034  0.038   0.042
#> Item_3   3.144  2.659     NA -0.047  0.057  0.071 -0.037 -0.050  0.033  -0.047
#> Item_4   1.279  1.040  2.173     NA  0.029  0.035  0.032 -0.039 -0.038   0.032
#> Item_5   0.198  0.932  3.291  0.862     NA  0.029  0.048 -0.028  0.027  -0.027
#> Item_6   0.204  1.059  5.032  1.197  0.830     NA  0.020  0.034  0.023   0.041
#> Item_7   0.277  0.804  1.362  1.018  2.337  0.392     NA -0.023 -0.040   0.038
#> Item_8   0.402  1.143  2.484  1.545  0.774  1.125  0.518     NA -0.027   0.025
#> Item_9   0.067  1.478  1.062  1.471  0.717  0.523  1.631  0.726     NA   0.052
#> Item_10  4.303  1.756  2.228  1.054  0.744  1.650  1.453  0.610  2.692      NA
#> Item_11  0.361  3.951  0.988  1.197  0.910  1.320  0.139  2.072  2.295   1.368
#> Item_12  1.718  0.586  1.038  0.725  0.218  0.454  1.323  0.473  0.322   0.604
#> Item_13  0.194  3.315  1.368  1.158  0.797  0.517  0.045  0.395  0.200   2.690
#> Item_14  0.074  0.804  3.798  2.946  0.648  0.228  0.223  0.464  0.693   0.378
#> Item_15  0.164  1.266  0.869  6.382  4.884  1.983  0.583  1.563  1.130   2.406
#>         Item_11 Item_12 Item_13 Item_14 Item_15
#> Item_1   -0.019   0.041  -0.014  -0.009   0.013
#> Item_2   -0.063   0.024   0.058   0.028   0.036
#> Item_3   -0.031   0.032  -0.037  -0.062   0.029
#> Item_4    0.035   0.027  -0.034  -0.054   0.080
#> Item_5   -0.030   0.015   0.028  -0.025  -0.070
#> Item_6    0.036   0.021   0.023  -0.015   0.045
#> Item_7    0.012   0.036  -0.007   0.015   0.024
#> Item_8   -0.046   0.022   0.020  -0.022   0.040
#> Item_9   -0.048  -0.018  -0.014   0.026   0.034
#> Item_10   0.037   0.025   0.052  -0.019  -0.049
#> Item_11      NA   0.010  -0.024  -0.053   0.021
#> Item_12   0.094      NA  -0.016  -0.016   0.012
#> Item_13   0.555   0.254      NA  -0.013   0.022
#> Item_14   2.862   0.255   0.159      NA   0.009
#> Item_15   0.424   0.148   0.486   0.073      NA
#> 
#> $D2
#>         Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> Item_1      NA -0.029  0.031 -0.026  0.021 -0.042 -0.006  0.039  0.010  -0.029
#> Item_2   0.861     NA  0.047 -0.060  0.049 -0.061 -0.036 -0.036  0.043   0.035
#> Item_3   0.973  2.251     NA  0.064  0.039 -0.029  0.033  0.045  0.042   0.031
#> Item_4   0.673  3.546  4.082     NA  0.034  0.045  0.029  0.033 -0.031  -0.053
#> Item_5   0.455  2.426  1.534  1.151     NA -0.042 -0.043 -0.030  0.034  -0.026
#> Item_6   1.750  3.773  0.861  2.016  1.752     NA -0.038  0.023 -0.060   0.022
#> Item_7   0.040  1.293  1.085  0.831  1.854  1.475     NA  0.023 -0.021   0.025
#> Item_8   1.530  1.310  1.985  1.122  0.915  0.545  0.539     NA -0.028  -0.026
#> Item_9   0.106  1.832  1.780  0.934  1.137  3.654  0.453  0.760     NA   0.037
#> Item_10  0.815  1.230  0.958  2.779  0.685  0.466  0.622  0.690  1.393      NA
#> Item_11  0.822  0.826  0.800  1.861  0.336  0.218  0.062  2.847  0.285   0.407
#> Item_12  0.087  1.536  1.097  0.639  0.309  0.454  0.096  0.556  0.398   0.520
#> Item_13  0.582  1.458  0.844  1.093  0.620  0.403  0.401  0.619  1.349   0.478
#> Item_14  0.403  1.550  0.914  1.292  0.544  0.799  0.534  3.153  0.137   0.344
#> Item_15  0.338  1.306  1.241  0.979  2.379  1.262  0.615  0.292  0.276   0.325
#>         Item_11 Item_12 Item_13 Item_14 Item_15
#> Item_1    0.029   0.009   0.024   0.020   0.018
#> Item_2    0.029  -0.039   0.038   0.039  -0.036
#> Item_3    0.028   0.033   0.029   0.030   0.035
#> Item_4    0.043  -0.025  -0.033   0.036  -0.031
#> Item_5    0.018  -0.018  -0.025   0.023  -0.049
#> Item_6    0.015  -0.021   0.020  -0.028  -0.036
#> Item_7   -0.008  -0.010  -0.020   0.023  -0.025
#> Item_8    0.053  -0.024  -0.025   0.056   0.017
#> Item_9   -0.017  -0.020  -0.037   0.012  -0.017
#> Item_10   0.020   0.023  -0.022   0.019   0.018
#> Item_11      NA   0.032   0.017   0.018   0.024
#> Item_12   1.010      NA   0.014  -0.022  -0.011
#> Item_13   0.280   0.185      NA   0.022  -0.012
#> Item_14   0.333   0.480   0.479      NA  -0.007
#> Item_15   0.580   0.131   0.145   0.048      NA
#> 
plot(mod_configural)

plot(mod_configural, type = 'info')

plot(mod_configural, type = 'trace')

plot(mod_configural, type = 'trace', which.items = 1:4)

itemplot(mod_configural, 2)

itemplot(mod_configural, 2, type = 'RE')


anova(mod_metric, mod_configural) #equal slopes only
#>                     AIC    SABIC       HQ      BIC    logLik    X2 df p
#> mod_metric     35937.61 36046.68 36030.15 36189.65 -17923.80           
#> mod_configural 35927.53 36072.96 36050.92 36263.58 -17903.76 40.08 15 0
anova(mod_scalar2, mod_metric) #equal intercepts, free variance and mean
#>                  AIC    SABIC       HQ      BIC    logLik      X2 df   p
#> mod_scalar2 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> mod_metric  35937.61 36046.68 36030.15 36189.65 -17923.80 -16.947 13 NaN
anova(mod_scalar1, mod_scalar2) #fix mean
#>                  AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> mod_scalar1 35893.96 35969.10 35957.71 36067.58 -17915.98               
#> mod_scalar2 35894.66 35972.22 35960.47 36073.89 -17915.33 1.296  1 0.255
anova(mod_fullconstrain, mod_scalar1) #fix variance
#>                        AIC    SABIC       HQ      BIC    logLik     X2 df p
#> mod_fullconstrain 35917.51 35990.22 35979.20 36085.53 -17928.75            
#> mod_scalar1       35893.96 35969.10 35957.71 36067.58 -17915.98 25.552  1 0

# compared all at once (in order of most constrained to least)
anova(mod_fullconstrain, mod_scalar2, mod_configural)
#>                        AIC    SABIC       HQ      BIC    logLik     X2 df     p
#> mod_fullconstrain 35917.51 35990.22 35979.20 36085.53 -17928.75                
#> mod_scalar2       35894.66 35972.22 35960.47 36073.89 -17915.33 26.848  2     0
#> mod_configural    35927.53 36072.96 36050.92 36263.58 -17903.76 23.133 28 0.726


# test whether first 6 slopes should be equal across groups
values <- multipleGroup(dat, 1, group = group, pars = 'values')
values
#>     group    item     class   name parnum  value lbound ubound   est const
#> 1      D1  Item_1      dich     a1      1  0.851   -Inf    Inf  TRUE  none
#> 2      D1  Item_1      dich      d      2  0.541   -Inf    Inf  TRUE  none
#> 3      D1  Item_1      dich      g      3  0.000      0      1 FALSE  none
#> 4      D1  Item_1      dich      u      4  1.000      0      1 FALSE  none
#> 5      D1  Item_2      dich     a1      5  0.851   -Inf    Inf  TRUE  none
#> 6      D1  Item_2      dich      d      6 -0.536   -Inf    Inf  TRUE  none
#> 7      D1  Item_2      dich      g      7  0.000      0      1 FALSE  none
#> 8      D1  Item_2      dich      u      8  1.000      0      1 FALSE  none
#> 9      D1  Item_3      dich     a1      9  0.851   -Inf    Inf  TRUE  none
#> 10     D1  Item_3      dich      d     10 -0.220   -Inf    Inf  TRUE  none
#> 11     D1  Item_3      dich      g     11  0.000      0      1 FALSE  none
#> 12     D1  Item_3      dich      u     12  1.000      0      1 FALSE  none
#> 13     D1  Item_4      dich     a1     13  0.851   -Inf    Inf  TRUE  none
#> 14     D1  Item_4      dich      d     14  0.941   -Inf    Inf  TRUE  none
#> 15     D1  Item_4      dich      g     15  0.000      0      1 FALSE  none
#> 16     D1  Item_4      dich      u     16  1.000      0      1 FALSE  none
#> 17     D1  Item_5      dich     a1     17  0.851   -Inf    Inf  TRUE  none
#> 18     D1  Item_5      dich      d     18  0.190   -Inf    Inf  TRUE  none
#> 19     D1  Item_5      dich      g     19  0.000      0      1 FALSE  none
#> 20     D1  Item_5      dich      u     20  1.000      0      1 FALSE  none
#> 21     D1  Item_6      dich     a1     21  0.851   -Inf    Inf  TRUE  none
#> 22     D1  Item_6      dich      d     22  0.752   -Inf    Inf  TRUE  none
#> 23     D1  Item_6      dich      g     23  0.000      0      1 FALSE  none
#> 24     D1  Item_6      dich      u     24  1.000      0      1 FALSE  none
#> 25     D1  Item_7      dich     a1     25  0.851   -Inf    Inf  TRUE  none
#> 26     D1  Item_7      dich      d     26  0.927   -Inf    Inf  TRUE  none
#> 27     D1  Item_7      dich      g     27  0.000      0      1 FALSE  none
#> 28     D1  Item_7      dich      u     28  1.000      0      1 FALSE  none
#> 29     D1  Item_8      dich     a1     29  0.851   -Inf    Inf  TRUE  none
#> 30     D1  Item_8      dich      d     30 -0.339   -Inf    Inf  TRUE  none
#> 31     D1  Item_8      dich      g     31  0.000      0      1 FALSE  none
#> 32     D1  Item_8      dich      u     32  1.000      0      1 FALSE  none
#> 33     D1  Item_9      dich     a1     33  0.851   -Inf    Inf  TRUE  none
#> 34     D1  Item_9      dich      d     34 -1.019   -Inf    Inf  TRUE  none
#> 35     D1  Item_9      dich      g     35  0.000      0      1 FALSE  none
#> 36     D1  Item_9      dich      u     36  1.000      0      1 FALSE  none
#> 37     D1 Item_10      dich     a1     37  0.851   -Inf    Inf  TRUE  none
#> 38     D1 Item_10      dich      d     38 -1.178   -Inf    Inf  TRUE  none
#> 39     D1 Item_10      dich      g     39  0.000      0      1 FALSE  none
#> 40     D1 Item_10      dich      u     40  1.000      0      1 FALSE  none
#> 41     D1 Item_11      dich     a1     41  0.851   -Inf    Inf  TRUE  none
#> 42     D1 Item_11      dich      d     42  1.228   -Inf    Inf  TRUE  none
#> 43     D1 Item_11      dich      g     43  0.000      0      1 FALSE  none
#> 44     D1 Item_11      dich      u     44  1.000      0      1 FALSE  none
#> 45     D1 Item_12      dich     a1     45  0.851   -Inf    Inf  TRUE  none
#> 46     D1 Item_12      dich      d     46 -0.150   -Inf    Inf  TRUE  none
#> 47     D1 Item_12      dich      g     47  0.000      0      1 FALSE  none
#> 48     D1 Item_12      dich      u     48  1.000      0      1 FALSE  none
#> 49     D1 Item_13      dich     a1     49  0.851   -Inf    Inf  TRUE  none
#> 50     D1 Item_13      dich      d     50  0.419   -Inf    Inf  TRUE  none
#> 51     D1 Item_13      dich      g     51  0.000      0      1 FALSE  none
#> 52     D1 Item_13      dich      u     52  1.000      0      1 FALSE  none
#> 53     D1 Item_14      dich     a1     53  0.851   -Inf    Inf  TRUE  none
#> 54     D1 Item_14      dich      d     54  0.455   -Inf    Inf  TRUE  none
#> 55     D1 Item_14      dich      g     55  0.000      0      1 FALSE  none
#> 56     D1 Item_14      dich      u     56  1.000      0      1 FALSE  none
#> 57     D1 Item_15      dich     a1     57  0.851   -Inf    Inf  TRUE  none
#> 58     D1 Item_15      dich      d     58 -0.047   -Inf    Inf  TRUE  none
#> 59     D1 Item_15      dich      g     59  0.000      0      1 FALSE  none
#> 60     D1 Item_15      dich      u     60  1.000      0      1 FALSE  none
#> 61     D1   GROUP GroupPars MEAN_1     61  0.000   -Inf    Inf FALSE  none
#> 62     D1   GROUP GroupPars COV_11     62  1.000      0    Inf FALSE  none
#> 63     D2  Item_1      dich     a1     63  0.851   -Inf    Inf  TRUE  none
#> 64     D2  Item_1      dich      d     64  0.541   -Inf    Inf  TRUE  none
#> 65     D2  Item_1      dich      g     65  0.000      0      1 FALSE  none
#> 66     D2  Item_1      dich      u     66  1.000      0      1 FALSE  none
#> 67     D2  Item_2      dich     a1     67  0.851   -Inf    Inf  TRUE  none
#> 68     D2  Item_2      dich      d     68 -0.536   -Inf    Inf  TRUE  none
#> 69     D2  Item_2      dich      g     69  0.000      0      1 FALSE  none
#> 70     D2  Item_2      dich      u     70  1.000      0      1 FALSE  none
#> 71     D2  Item_3      dich     a1     71  0.851   -Inf    Inf  TRUE  none
#> 72     D2  Item_3      dich      d     72 -0.220   -Inf    Inf  TRUE  none
#> 73     D2  Item_3      dich      g     73  0.000      0      1 FALSE  none
#> 74     D2  Item_3      dich      u     74  1.000      0      1 FALSE  none
#> 75     D2  Item_4      dich     a1     75  0.851   -Inf    Inf  TRUE  none
#> 76     D2  Item_4      dich      d     76  0.941   -Inf    Inf  TRUE  none
#> 77     D2  Item_4      dich      g     77  0.000      0      1 FALSE  none
#> 78     D2  Item_4      dich      u     78  1.000      0      1 FALSE  none
#> 79     D2  Item_5      dich     a1     79  0.851   -Inf    Inf  TRUE  none
#> 80     D2  Item_5      dich      d     80  0.190   -Inf    Inf  TRUE  none
#> 81     D2  Item_5      dich      g     81  0.000      0      1 FALSE  none
#> 82     D2  Item_5      dich      u     82  1.000      0      1 FALSE  none
#> 83     D2  Item_6      dich     a1     83  0.851   -Inf    Inf  TRUE  none
#> 84     D2  Item_6      dich      d     84  0.752   -Inf    Inf  TRUE  none
#> 85     D2  Item_6      dich      g     85  0.000      0      1 FALSE  none
#> 86     D2  Item_6      dich      u     86  1.000      0      1 FALSE  none
#> 87     D2  Item_7      dich     a1     87  0.851   -Inf    Inf  TRUE  none
#> 88     D2  Item_7      dich      d     88  0.927   -Inf    Inf  TRUE  none
#> 89     D2  Item_7      dich      g     89  0.000      0      1 FALSE  none
#> 90     D2  Item_7      dich      u     90  1.000      0      1 FALSE  none
#> 91     D2  Item_8      dich     a1     91  0.851   -Inf    Inf  TRUE  none
#> 92     D2  Item_8      dich      d     92 -0.339   -Inf    Inf  TRUE  none
#> 93     D2  Item_8      dich      g     93  0.000      0      1 FALSE  none
#> 94     D2  Item_8      dich      u     94  1.000      0      1 FALSE  none
#> 95     D2  Item_9      dich     a1     95  0.851   -Inf    Inf  TRUE  none
#> 96     D2  Item_9      dich      d     96 -1.019   -Inf    Inf  TRUE  none
#> 97     D2  Item_9      dich      g     97  0.000      0      1 FALSE  none
#> 98     D2  Item_9      dich      u     98  1.000      0      1 FALSE  none
#> 99     D2 Item_10      dich     a1     99  0.851   -Inf    Inf  TRUE  none
#> 100    D2 Item_10      dich      d    100 -1.178   -Inf    Inf  TRUE  none
#> 101    D2 Item_10      dich      g    101  0.000      0      1 FALSE  none
#> 102    D2 Item_10      dich      u    102  1.000      0      1 FALSE  none
#> 103    D2 Item_11      dich     a1    103  0.851   -Inf    Inf  TRUE  none
#> 104    D2 Item_11      dich      d    104  1.228   -Inf    Inf  TRUE  none
#> 105    D2 Item_11      dich      g    105  0.000      0      1 FALSE  none
#> 106    D2 Item_11      dich      u    106  1.000      0      1 FALSE  none
#> 107    D2 Item_12      dich     a1    107  0.851   -Inf    Inf  TRUE  none
#> 108    D2 Item_12      dich      d    108 -0.150   -Inf    Inf  TRUE  none
#> 109    D2 Item_12      dich      g    109  0.000      0      1 FALSE  none
#> 110    D2 Item_12      dich      u    110  1.000      0      1 FALSE  none
#> 111    D2 Item_13      dich     a1    111  0.851   -Inf    Inf  TRUE  none
#> 112    D2 Item_13      dich      d    112  0.419   -Inf    Inf  TRUE  none
#> 113    D2 Item_13      dich      g    113  0.000      0      1 FALSE  none
#> 114    D2 Item_13      dich      u    114  1.000      0      1 FALSE  none
#> 115    D2 Item_14      dich     a1    115  0.851   -Inf    Inf  TRUE  none
#> 116    D2 Item_14      dich      d    116  0.455   -Inf    Inf  TRUE  none
#> 117    D2 Item_14      dich      g    117  0.000      0      1 FALSE  none
#> 118    D2 Item_14      dich      u    118  1.000      0      1 FALSE  none
#> 119    D2 Item_15      dich     a1    119  0.851   -Inf    Inf  TRUE  none
#> 120    D2 Item_15      dich      d    120 -0.047   -Inf    Inf  TRUE  none
#> 121    D2 Item_15      dich      g    121  0.000      0      1 FALSE  none
#> 122    D2 Item_15      dich      u    122  1.000      0      1 FALSE  none
#> 123    D2   GROUP GroupPars MEAN_1    123  0.000   -Inf    Inf FALSE  none
#> 124    D2   GROUP GroupPars COV_11    124  1.000      0    Inf FALSE  none
#>     nconst prior.type prior_1 prior_2
#> 1     none       none     NaN     NaN
#> 2     none       none     NaN     NaN
#> 3     none       none     NaN     NaN
#> 4     none       none     NaN     NaN
#> 5     none       none     NaN     NaN
#> 6     none       none     NaN     NaN
#> 7     none       none     NaN     NaN
#> 8     none       none     NaN     NaN
#> 9     none       none     NaN     NaN
#> 10    none       none     NaN     NaN
#> 11    none       none     NaN     NaN
#> 12    none       none     NaN     NaN
#> 13    none       none     NaN     NaN
#> 14    none       none     NaN     NaN
#> 15    none       none     NaN     NaN
#> 16    none       none     NaN     NaN
#> 17    none       none     NaN     NaN
#> 18    none       none     NaN     NaN
#> 19    none       none     NaN     NaN
#> 20    none       none     NaN     NaN
#> 21    none       none     NaN     NaN
#> 22    none       none     NaN     NaN
#> 23    none       none     NaN     NaN
#> 24    none       none     NaN     NaN
#> 25    none       none     NaN     NaN
#> 26    none       none     NaN     NaN
#> 27    none       none     NaN     NaN
#> 28    none       none     NaN     NaN
#> 29    none       none     NaN     NaN
#> 30    none       none     NaN     NaN
#> 31    none       none     NaN     NaN
#> 32    none       none     NaN     NaN
#> 33    none       none     NaN     NaN
#> 34    none       none     NaN     NaN
#> 35    none       none     NaN     NaN
#> 36    none       none     NaN     NaN
#> 37    none       none     NaN     NaN
#> 38    none       none     NaN     NaN
#> 39    none       none     NaN     NaN
#> 40    none       none     NaN     NaN
#> 41    none       none     NaN     NaN
#> 42    none       none     NaN     NaN
#> 43    none       none     NaN     NaN
#> 44    none       none     NaN     NaN
#> 45    none       none     NaN     NaN
#> 46    none       none     NaN     NaN
#> 47    none       none     NaN     NaN
#> 48    none       none     NaN     NaN
#> 49    none       none     NaN     NaN
#> 50    none       none     NaN     NaN
#> 51    none       none     NaN     NaN
#> 52    none       none     NaN     NaN
#> 53    none       none     NaN     NaN
#> 54    none       none     NaN     NaN
#> 55    none       none     NaN     NaN
#> 56    none       none     NaN     NaN
#> 57    none       none     NaN     NaN
#> 58    none       none     NaN     NaN
#> 59    none       none     NaN     NaN
#> 60    none       none     NaN     NaN
#> 61    none       none     NaN     NaN
#> 62    none       none     NaN     NaN
#> 63    none       none     NaN     NaN
#> 64    none       none     NaN     NaN
#> 65    none       none     NaN     NaN
#> 66    none       none     NaN     NaN
#> 67    none       none     NaN     NaN
#> 68    none       none     NaN     NaN
#> 69    none       none     NaN     NaN
#> 70    none       none     NaN     NaN
#> 71    none       none     NaN     NaN
#> 72    none       none     NaN     NaN
#> 73    none       none     NaN     NaN
#> 74    none       none     NaN     NaN
#> 75    none       none     NaN     NaN
#> 76    none       none     NaN     NaN
#> 77    none       none     NaN     NaN
#> 78    none       none     NaN     NaN
#> 79    none       none     NaN     NaN
#> 80    none       none     NaN     NaN
#> 81    none       none     NaN     NaN
#> 82    none       none     NaN     NaN
#> 83    none       none     NaN     NaN
#> 84    none       none     NaN     NaN
#> 85    none       none     NaN     NaN
#> 86    none       none     NaN     NaN
#> 87    none       none     NaN     NaN
#> 88    none       none     NaN     NaN
#> 89    none       none     NaN     NaN
#> 90    none       none     NaN     NaN
#> 91    none       none     NaN     NaN
#> 92    none       none     NaN     NaN
#> 93    none       none     NaN     NaN
#> 94    none       none     NaN     NaN
#> 95    none       none     NaN     NaN
#> 96    none       none     NaN     NaN
#> 97    none       none     NaN     NaN
#> 98    none       none     NaN     NaN
#> 99    none       none     NaN     NaN
#> 100   none       none     NaN     NaN
#> 101   none       none     NaN     NaN
#> 102   none       none     NaN     NaN
#> 103   none       none     NaN     NaN
#> 104   none       none     NaN     NaN
#> 105   none       none     NaN     NaN
#> 106   none       none     NaN     NaN
#> 107   none       none     NaN     NaN
#> 108   none       none     NaN     NaN
#> 109   none       none     NaN     NaN
#> 110   none       none     NaN     NaN
#> 111   none       none     NaN     NaN
#> 112   none       none     NaN     NaN
#> 113   none       none     NaN     NaN
#> 114   none       none     NaN     NaN
#> 115   none       none     NaN     NaN
#> 116   none       none     NaN     NaN
#> 117   none       none     NaN     NaN
#> 118   none       none     NaN     NaN
#> 119   none       none     NaN     NaN
#> 120   none       none     NaN     NaN
#> 121   none       none     NaN     NaN
#> 122   none       none     NaN     NaN
#> 123   none       none     NaN     NaN
#> 124   none       none     NaN     NaN
constrain <- list(c(1, 63), c(5,67), c(9,71), c(13,75), c(17,79), c(21,83))
equalslopes <- multipleGroup(dat, 1, group = group, constrain = constrain)
anova(equalslopes, mod_configural)
#>                     AIC    SABIC       HQ      BIC    logLik     X2 df     p
#> equalslopes    35935.51 36066.40 36046.56 36237.96 -17913.76                
#> mod_configural 35927.53 36072.96 36050.92 36263.58 -17903.76 19.983  6 0.003

# same as above, but using mirt.model syntax
newmodel <- '
    F = 1-15
    CONSTRAINB = (1-6, a1)'
equalslopes <- multipleGroup(dat, newmodel, group = group)
coef(equalslopes, simplify=TRUE)
#> $D1
#> $items
#>            a1      d g u
#> Item_1  1.246  0.546 0 1
#> Item_2  1.356 -0.720 0 1
#> Item_3  1.199 -0.201 0 1
#> Item_4  1.006  0.861 0 1
#> Item_5  1.224  0.131 0 1
#> Item_6  0.515  0.673 0 1
#> Item_7  1.305  0.999 0 1
#> Item_8  0.943 -0.328 0 1
#> Item_9  0.916 -1.058 0 1
#> Item_10 0.731 -1.079 0 1
#> Item_11 0.851  1.186 0 1
#> Item_12 1.515 -0.251 0 1
#> Item_13 1.322  0.444 0 1
#> Item_14 1.058  0.451 0 1
#> Item_15 0.885 -0.061 0 1
#> 
#> $means
#> F 
#> 0 
#> 
#> $cov
#>   F
#> F 1
#> 
#> 
#> $D2
#> $items
#>            a1      d g u
#> Item_1  1.246  0.599 0 1
#> Item_2  1.356 -0.462 0 1
#> Item_3  1.199 -0.264 0 1
#> Item_4  1.006  1.001 0 1
#> Item_5  1.224  0.266 0 1
#> Item_6  0.515  0.631 0 1
#> Item_7  1.374  1.031 0 1
#> Item_8  1.268 -0.367 0 1
#> Item_9  1.061 -0.952 0 1
#> Item_10 0.826 -1.115 0 1
#> Item_11 1.282  1.308 0 1
#> Item_12 1.625 -0.107 0 1
#> Item_13 1.523  0.493 0 1
#> Item_14 1.545  0.532 0 1
#> Item_15 0.885 -0.030 0 1
#> 
#> $means
#> F 
#> 0 
#> 
#> $cov
#>   F
#> F 1
#> 
#> 

############
# vertical scaling (i.e., equating when groups answer items others do not)
dat2 <- dat
dat2[group == 'D1', 1:2] <- dat2[group != 'D1', 14:15] <- NA
head(dat2)
#>      Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> [1,]     NA     NA      1      1      1      0      1      1      0       0
#> [2,]     NA     NA      1      1      1      1      1      0      1       0
#> [3,]     NA     NA      0      1      1      1      1      1      0       0
#> [4,]     NA     NA      1      1      1      1      1      1      0       0
#> [5,]     NA     NA      1      0      1      1      1      1      1       0
#> [6,]     NA     NA      1      1      1      1      1      0      0       0
#>      Item_11 Item_12 Item_13 Item_14 Item_15
#> [1,]       1       1       0       1       1
#> [2,]       0       1       1       1       1
#> [3,]       1       0       1       1       1
#> [4,]       1       1       1       1       1
#> [5,]       1       1       1       1       1
#> [6,]       1       1       1       1       0
tail(dat2)
#>         Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> [1995,]      1      1      0      0      0      1      1      1      0       0
#> [1996,]      1      0      1      1      1      0      1      1      1       1
#> [1997,]      0      1      0      0      0      1      0      1      0       0
#> [1998,]      0      1      0      0      1      0      0      0      0       0
#> [1999,]      1      1      0      1      1      1      1      1      0       0
#> [2000,]      0      0      0      0      0      0      1      0      0       0
#>         Item_11 Item_12 Item_13 Item_14 Item_15
#> [1995,]       1       0       1      NA      NA
#> [1996,]       1       1       0      NA      NA
#> [1997,]       1       0       1      NA      NA
#> [1998,]       0       1       0      NA      NA
#> [1999,]       1       1       1      NA      NA
#> [2000,]       0       0       0      NA      NA

# items with missing responses need to be constrained across groups for identification
nms <- colnames(dat2)
mod <- multipleGroup(dat2, 1, group, invariance = nms[c(1:2, 14:15)])

# this will throw an error without proper constraints (SEs cannot be computed either)
# mod <- multipleGroup(dat2, 1, group)

# model still does not have anchors, therefore need to add a few (here use items 3-5)
mod_anchor <- multipleGroup(dat2, 1, group,
                            invariance = c(nms[c(1:5, 14:15)], 'free_means', 'free_var'))
coef(mod_anchor, simplify=TRUE)
#> $D1
#> $items
#>            a1      d g u
#> Item_1  1.108  0.542 0 1
#> Item_2  1.160 -0.564 0 1
#> Item_3  1.073 -0.272 0 1
#> Item_4  0.871  0.895 0 1
#> Item_5  1.089  0.160 0 1
#> Item_6  0.582  0.685 0 1
#> Item_7  1.292  1.009 0 1
#> Item_8  0.906 -0.328 0 1
#> Item_9  0.855 -1.049 0 1
#> Item_10 0.746 -1.089 0 1
#> Item_11 0.866  1.200 0 1
#> Item_12 1.432 -0.247 0 1
#> Item_13 1.244  0.440 0 1
#> Item_14 1.000  0.449 0 1
#> Item_15 0.852 -0.061 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> 
#> $D2
#> $items
#>            a1      d g u
#> Item_1  1.108  0.542 0 1
#> Item_2  1.160 -0.564 0 1
#> Item_3  1.073 -0.272 0 1
#> Item_4  0.871  0.895 0 1
#> Item_5  1.089  0.160 0 1
#> Item_6  0.374  0.596 0 1
#> Item_7  1.090  0.942 0 1
#> Item_8  0.988 -0.437 0 1
#> Item_9  0.854 -1.018 0 1
#> Item_10 0.657 -1.164 0 1
#> Item_11 1.023  1.228 0 1
#> Item_12 1.324 -0.207 0 1
#> Item_13 1.212  0.399 0 1
#> Item_14 1.000  0.449 0 1
#> Item_15 0.852 -0.061 0 1
#> 
#> $means
#>    F1 
#> 0.077 
#> 
#> $cov
#>       F1
#> F1 1.658
#> 
#> 

# check if identified by computing information matrix
mod_anchor <- multipleGroup(dat2, 1, group, pars = mod2values(mod_anchor), TOL=NaN, SE=TRUE,
                            invariance = c(nms[c(1:5, 14:15)], 'free_means', 'free_var'))
mod_anchor
#> 
#> Call:
#> multipleGroup(data = dat2, model = 1, group = group, invariance = c(nms[c(1:5, 
#>     14:15)], "free_means", "free_var"), pars = mod2values(mod_anchor), 
#>     TOL = NaN, SE = TRUE)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within NaN tolerance after 1 EM iterations.
#> mirt version: 1.45.6 
#> M-step optimizer: nlminb 
#> EM acceleration: Ramsay 
#> Number of rectangular quadrature: 61
#> Latent density type: Gaussian 
#> 
#> Information matrix estimated with method: Oakes
#> Second-order test: model is a possible local maximum
#> Condition number of information matrix =  87.67822
#> 
#> Log-likelihood = -15563.53
#> Estimated parameters: 48 
#> AIC = 31223.06
#> BIC = 31491.91; SABIC = 31339.41
#> 
coef(mod_anchor)
#> $D1
#> $Item_1
#>            a1     d  g  u
#> par     1.108 0.542  0  1
#> CI_2.5  0.842 0.330 NA NA
#> CI_97.5 1.375 0.753 NA NA
#> 
#> $Item_2
#>            a1      d  g  u
#> par     1.160 -0.564  0  1
#> CI_2.5  0.883 -0.785 NA NA
#> CI_97.5 1.436 -0.342 NA NA
#> 
#> $Item_3
#>            a1      d  g  u
#> par     1.073 -0.272  0  1
#> CI_2.5  0.896 -0.408 NA NA
#> CI_97.5 1.250 -0.137 NA NA
#> 
#> $Item_4
#>            a1     d  g  u
#> par     0.871 0.895  0  1
#> CI_2.5  0.710 0.765 NA NA
#> CI_97.5 1.031 1.026 NA NA
#> 
#> $Item_5
#>            a1     d  g  u
#> par     1.089 0.160  0  1
#> CI_2.5  0.906 0.024 NA NA
#> CI_97.5 1.272 0.295 NA NA
#> 
#> $Item_6
#>            a1     d  g  u
#> par     0.582 0.685  0  1
#> CI_2.5  0.405 0.543 NA NA
#> CI_97.5 0.760 0.828 NA NA
#> 
#> $Item_7
#>            a1     d  g  u
#> par     1.292 1.009  0  1
#> CI_2.5  1.027 0.818 NA NA
#> CI_97.5 1.558 1.199 NA NA
#> 
#> $Item_8
#>            a1      d  g  u
#> par     0.906 -0.328  0  1
#> CI_2.5  0.702 -0.476 NA NA
#> CI_97.5 1.111 -0.179 NA NA
#> 
#> $Item_9
#>            a1      d  g  u
#> par     0.855 -1.049  0  1
#> CI_2.5  0.642 -1.216 NA NA
#> CI_97.5 1.069 -0.882 NA NA
#> 
#> $Item_10
#>            a1      d  g  u
#> par     0.746 -1.089  0  1
#> CI_2.5  0.542 -1.253 NA NA
#> CI_97.5 0.950 -0.926 NA NA
#> 
#> $Item_11
#>            a1     d  g  u
#> par     0.866 1.200  0  1
#> CI_2.5  0.651 1.025 NA NA
#> CI_97.5 1.081 1.375 NA NA
#> 
#> $Item_12
#>            a1      d  g  u
#> par     1.432 -0.247  0  1
#> CI_2.5  1.152 -0.421 NA NA
#> CI_97.5 1.713 -0.074 NA NA
#> 
#> $Item_13
#>            a1     d  g  u
#> par     1.244 0.440  0  1
#> CI_2.5  0.994 0.273 NA NA
#> CI_97.5 1.494 0.607 NA NA
#> 
#> $Item_14
#>            a1     d  g  u
#> par     1.000 0.449  0  1
#> CI_2.5  0.785 0.294 NA NA
#> CI_97.5 1.215 0.603 NA NA
#> 
#> $Item_15
#>            a1      d  g  u
#> par     0.852 -0.061  0  1
#> CI_2.5  0.655 -0.205 NA NA
#> CI_97.5 1.049  0.083 NA NA
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0      1
#> CI_2.5      NA     NA
#> CI_97.5     NA     NA
#> 
#> 
#> $D2
#> $Item_1
#>            a1     d  g  u
#> par     1.108 0.542  0  1
#> CI_2.5  0.842 0.330 NA NA
#> CI_97.5 1.375 0.753 NA NA
#> 
#> $Item_2
#>            a1      d  g  u
#> par     1.160 -0.564  0  1
#> CI_2.5  0.883 -0.785 NA NA
#> CI_97.5 1.436 -0.342 NA NA
#> 
#> $Item_3
#>            a1      d  g  u
#> par     1.073 -0.272  0  1
#> CI_2.5  0.896 -0.408 NA NA
#> CI_97.5 1.250 -0.137 NA NA
#> 
#> $Item_4
#>            a1     d  g  u
#> par     0.871 0.895  0  1
#> CI_2.5  0.710 0.765 NA NA
#> CI_97.5 1.031 1.026 NA NA
#> 
#> $Item_5
#>            a1     d  g  u
#> par     1.089 0.160  0  1
#> CI_2.5  0.906 0.024 NA NA
#> CI_97.5 1.272 0.295 NA NA
#> 
#> $Item_6
#>            a1     d  g  u
#> par     0.374 0.596  0  1
#> CI_2.5  0.236 0.454 NA NA
#> CI_97.5 0.512 0.738 NA NA
#> 
#> $Item_7
#>            a1     d  g  u
#> par     1.090 0.942  0  1
#> CI_2.5  0.823 0.722 NA NA
#> CI_97.5 1.357 1.163 NA NA
#> 
#> $Item_8
#>            a1      d  g  u
#> par     0.988 -0.437  0  1
#> CI_2.5  0.747 -0.636 NA NA
#> CI_97.5 1.228 -0.239 NA NA
#> 
#> $Item_9
#>            a1      d  g  u
#> par     0.854 -1.018  0  1
#> CI_2.5  0.635 -1.219 NA NA
#> CI_97.5 1.072 -0.817 NA NA
#> 
#> $Item_10
#>            a1      d  g  u
#> par     0.657 -1.164  0  1
#> CI_2.5  0.469 -1.350 NA NA
#> CI_97.5 0.845 -0.978 NA NA
#> 
#> $Item_11
#>            a1     d  g  u
#> par     1.023 1.228  0  1
#> CI_2.5  0.766 1.003 NA NA
#> CI_97.5 1.280 1.452 NA NA
#> 
#> $Item_12
#>            a1      d  g  u
#> par     1.324 -0.207  0  1
#> CI_2.5  1.010 -0.442 NA NA
#> CI_97.5 1.637  0.027 NA NA
#> 
#> $Item_13
#>            a1     d  g  u
#> par     1.212 0.399  0  1
#> CI_2.5  0.923 0.178 NA NA
#> CI_97.5 1.502 0.620 NA NA
#> 
#> $Item_14
#>            a1     d  g  u
#> par     1.000 0.449  0  1
#> CI_2.5  0.785 0.294 NA NA
#> CI_97.5 1.215 0.603 NA NA
#> 
#> $Item_15
#>            a1      d  g  u
#> par     0.852 -0.061  0  1
#> CI_2.5  0.655 -0.205 NA NA
#> CI_97.5 1.049  0.083 NA NA
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par      0.008  1.658
#> CI_2.5  -0.148  1.093
#> CI_97.5  0.164  2.222
#> 
#> 
coef(mod_anchor, printSE=TRUE)
#> $D1
#> $Item_1
#>        a1     d logit(g) logit(u)
#> par 1.108 0.542     -999      999
#> SE  0.136 0.108       NA       NA
#> 
#> $Item_2
#>        a1      d logit(g) logit(u)
#> par 1.160 -0.564     -999      999
#> SE  0.141  0.113       NA       NA
#> 
#> $Item_3
#>        a1      d logit(g) logit(u)
#> par 1.073 -0.272     -999      999
#> SE  0.090  0.069       NA       NA
#> 
#> $Item_4
#>        a1     d logit(g) logit(u)
#> par 0.871 0.895     -999      999
#> SE  0.082 0.066       NA       NA
#> 
#> $Item_5
#>        a1     d logit(g) logit(u)
#> par 1.089 0.160     -999      999
#> SE  0.093 0.069       NA       NA
#> 
#> $Item_6
#>        a1     d logit(g) logit(u)
#> par 0.582 0.685     -999      999
#> SE  0.091 0.073       NA       NA
#> 
#> $Item_7
#>        a1     d logit(g) logit(u)
#> par 1.292 1.009     -999      999
#> SE  0.135 0.097       NA       NA
#> 
#> $Item_8
#>        a1      d logit(g) logit(u)
#> par 0.906 -0.328     -999      999
#> SE  0.104  0.076       NA       NA
#> 
#> $Item_9
#>        a1      d logit(g) logit(u)
#> par 0.855 -1.049     -999      999
#> SE  0.109  0.085       NA       NA
#> 
#> $Item_10
#>        a1      d logit(g) logit(u)
#> par 0.746 -1.089     -999      999
#> SE  0.104  0.083       NA       NA
#> 
#> $Item_11
#>        a1     d logit(g) logit(u)
#> par 0.866 1.200     -999      999
#> SE  0.110 0.089       NA       NA
#> 
#> $Item_12
#>        a1      d logit(g) logit(u)
#> par 1.432 -0.247     -999      999
#> SE  0.143  0.089       NA       NA
#> 
#> $Item_13
#>        a1     d logit(g) logit(u)
#> par 1.244 0.440     -999      999
#> SE  0.127 0.085       NA       NA
#> 
#> $Item_14
#>       a1     d logit(g) logit(u)
#> par 1.00 0.449     -999      999
#> SE  0.11 0.079       NA       NA
#> 
#> $Item_15
#>        a1      d logit(g) logit(u)
#> par 0.852 -0.061     -999      999
#> SE  0.100  0.073       NA       NA
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> SE      NA     NA
#> 
#> 
#> $D2
#> $Item_1
#>        a1     d logit(g) logit(u)
#> par 1.108 0.542     -999      999
#> SE  0.136 0.108       NA       NA
#> 
#> $Item_2
#>        a1      d logit(g) logit(u)
#> par 1.160 -0.564     -999      999
#> SE  0.141  0.113       NA       NA
#> 
#> $Item_3
#>        a1      d logit(g) logit(u)
#> par 1.073 -0.272     -999      999
#> SE  0.090  0.069       NA       NA
#> 
#> $Item_4
#>        a1     d logit(g) logit(u)
#> par 0.871 0.895     -999      999
#> SE  0.082 0.066       NA       NA
#> 
#> $Item_5
#>        a1     d logit(g) logit(u)
#> par 1.089 0.160     -999      999
#> SE  0.093 0.069       NA       NA
#> 
#> $Item_6
#>        a1     d logit(g) logit(u)
#> par 0.374 0.596     -999      999
#> SE  0.071 0.073       NA       NA
#> 
#> $Item_7
#>        a1     d logit(g) logit(u)
#> par 1.090 0.942     -999      999
#> SE  0.136 0.112       NA       NA
#> 
#> $Item_8
#>        a1      d logit(g) logit(u)
#> par 0.988 -0.437     -999      999
#> SE  0.123  0.101       NA       NA
#> 
#> $Item_9
#>        a1      d logit(g) logit(u)
#> par 0.854 -1.018     -999      999
#> SE  0.111  0.102       NA       NA
#> 
#> $Item_10
#>        a1      d logit(g) logit(u)
#> par 0.657 -1.164     -999      999
#> SE  0.096  0.095       NA       NA
#> 
#> $Item_11
#>        a1     d logit(g) logit(u)
#> par 1.023 1.228     -999      999
#> SE  0.131 0.114       NA       NA
#> 
#> $Item_12
#>        a1      d logit(g) logit(u)
#> par 1.324 -0.207     -999      999
#> SE  0.160  0.120       NA       NA
#> 
#> $Item_13
#>        a1     d logit(g) logit(u)
#> par 1.212 0.399     -999      999
#> SE  0.148 0.113       NA       NA
#> 
#> $Item_14
#>       a1     d logit(g) logit(u)
#> par 1.00 0.449     -999      999
#> SE  0.11 0.079       NA       NA
#> 
#> $Item_15
#>        a1      d logit(g) logit(u)
#> par 0.852 -0.061     -999      999
#> SE  0.100  0.073       NA       NA
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par  0.008  1.658
#> SE   0.080  0.288
#> 
#> 


#############
# DIF test for each item (using all other items as anchors)
itemnames <- colnames(dat)
refmodel <- multipleGroup(dat, 1, group = group, SE=TRUE,
                          invariance=c('free_means', 'free_var', itemnames))

# loop over items (in practice, run in parallel to increase speed). May be better to use ?DIF
estmodels <- vector('list', ncol(dat))
for(i in 1:ncol(dat))
    estmodels[[i]] <- multipleGroup(dat, 1, group = group, verbose = FALSE,
                             invariance=c('free_means', 'free_var', itemnames[-i]))
anova(refmodel, estmodels[[1]])
#>                     AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel       35894.66 35972.22 35960.47 36073.89 -17915.33               
#> estmodels[[1]] 35898.45 35980.86 35968.37 36088.88 -17915.22 0.213  2 0.899
(anovas <- lapply(estmodels, function(x, refmodel) anova(refmodel, x),
   refmodel=refmodel))
#> [[1]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35898.45 35980.86 35968.37 36088.88 -17915.22 0.213  2 0.899
#> 
#> [[2]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.81 35979.22 35966.73 36087.24 -17914.41 1.847  2 0.397
#> 
#> [[3]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35893.66 35976.07 35963.58 36084.09 -17912.83 5.001  2 0.082
#> 
#> [[4]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35897.07 35979.48 35967.00 36087.50 -17914.54 1.586  2 0.453
#> 
#> [[5]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35897.97 35980.38 35967.89 36088.40 -17914.99 0.688  2 0.709
#> 
#> [[6]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df    p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33              
#> x        35894.87 35977.28 35964.79 36085.30 -17913.43 3.793  2 0.15
#> 
#> [[7]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35897.53 35979.94 35967.45 36087.96 -17914.76 1.131  2 0.568
#> 
#> [[8]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35897.20 35979.61 35967.12 36087.63 -17914.60 1.462  2 0.481
#> 
#> [[9]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35898.46 35980.87 35968.38 36088.89 -17915.23 0.197  2 0.906
#> 
#> [[10]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35897.67 35980.08 35967.59 36088.10 -17914.83 0.993  2 0.609
#> 
#> [[11]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.07 35978.48 35965.99 36086.50 -17914.03 2.593  2 0.273
#> 
#> [[12]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35897.57 35979.98 35967.49 36088.00 -17914.78 1.092  2 0.579
#> 
#> [[13]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35898.47 35980.88 35968.39 36088.90 -17915.23 0.192  2 0.908
#> 
#> [[14]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.15 35978.57 35966.08 36086.58 -17914.08 2.505  2 0.286
#> 
#> [[15]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.96 35979.37 35966.88 36087.39 -17914.48 1.699  2 0.428
#> 

# family-wise error control
p <- do.call(rbind, lapply(anovas, function(x) x[2, 'p']))
p.adjust(p, method = 'BH')
#>  [1] 0.9084069 0.8299980 0.8299980 0.8299980 0.8862364 0.8299980 0.8299980
#>  [8] 0.8299980 0.9084069 0.8299980 0.8299980 0.8299980 0.9084069 0.8299980
#> [15] 0.8299980

# same as above, except only test if slopes vary (1 df)
# constrain all intercepts
estmodels <- vector('list', ncol(dat))
for(i in 1:ncol(dat))
    estmodels[[i]] <- multipleGroup(dat, 1, group = group, verbose = FALSE,
                             invariance=c('free_means', 'free_var', 'intercepts',
                             itemnames[-i]))

(anovas <- lapply(estmodels, function(x, refmodel) anova(refmodel, x),
   refmodel=refmodel))
#> [[1]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.52 35976.50 35964.38 36081.35 -17915.26 0.143  1 0.705
#> 
#> [[2]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.60 35976.58 35964.46 36081.43 -17915.30 0.061  1 0.804
#> 
#> [[3]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35894.66 35974.65 35962.53 36079.49 -17914.33 1.997  1 0.158
#> 
#> [[4]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.39 35976.38 35964.26 36081.22 -17915.20 0.268  1 0.605
#> 
#> [[5]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.56 35976.54 35964.42 36081.39 -17915.28 0.103  1 0.748
#> 
#> [[6]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35893.62 35973.60 35961.48 36078.45 -17913.81 3.042  1 0.081
#> 
#> [[7]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35895.67 35975.66 35963.54 36080.50 -17914.84 0.985  1 0.321
#> 
#> [[8]]
#>               AIC    SABIC       HQ      BIC    logLik   X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33              
#> x        35896.27 35976.26 35964.13 36081.10 -17915.13 0.39  1 0.532
#> 
#> [[9]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.66 35976.64 35964.52 36081.49 -17915.33 0.001  1 0.979
#> 
#> [[10]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35896.12 35976.11 35963.99 36080.95 -17915.06 0.538  1 0.463
#> 
#> [[11]]
#>               AIC    SABIC       HQ      BIC    logLik   X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33              
#> x        35894.22 35974.21 35962.08 36079.05 -17914.11 2.44  1 0.118
#> 
#> [[12]]
#>               AIC    SABIC       HQ      BIC    logLik   X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33              
#> x        35895.87 35975.86 35963.74 36080.70 -17914.94 0.79  1 0.374
#> 
#> [[13]]
#>               AIC    SABIC       HQ      BIC    logLik   X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33              
#> x        35896.58 35976.57 35964.44 36081.41 -17915.29 0.08  1 0.778
#> 
#> [[14]]
#>               AIC    SABIC       HQ      BIC    logLik  X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33             
#> x        35894.16 35974.15 35962.02 36078.99 -17914.08 2.5  1 0.114
#> 
#> [[15]]
#>               AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> refmodel 35894.66 35972.22 35960.47 36073.89 -17915.33               
#> x        35894.99 35974.97 35962.85 36079.82 -17914.49 1.672  1 0.196
#> 

# quickly test with Wald test using DIF()
mod_configural2 <- multipleGroup(dat, 1, group = group, SE=TRUE)
DIF(mod_configural2, which.par = c('a1', 'd'), Wald=TRUE, p.adjust = 'fdr')
#> NOTE: No hyper-parameters were estimated in the DIF model. 
#>       For effective DIF testing, freeing the focal group hyper-parameters is recommended.
#>         groups      W df     p adj_p
#> Item_1   D1,D2  4.636  2 0.098 0.246
#> Item_2   D1,D2  7.265  2 0.026 0.099
#> Item_3   D1,D2 10.375  2 0.006 0.055
#> Item_4   D1,D2  3.210  2 0.201 0.335
#> Item_5   D1,D2  3.618  2 0.164 0.307
#> Item_6   D1,D2  0.820  2 0.664 0.804
#> Item_7   D1,D2  0.575  2  0.75 0.804
#> Item_8   D1,D2  5.782  2 0.056 0.167
#> Item_9   D1,D2  3.802  2 0.149 0.307
#> Item_10  D1,D2  0.722  2 0.697 0.804
#> Item_11  D1,D2  8.922  2 0.012 0.058
#> Item_12  D1,D2  2.340  2  0.31 0.463
#> Item_13  D1,D2  2.162  2 0.339 0.463
#> Item_14  D1,D2  9.848  2 0.007 0.055
#> Item_15  D1,D2  0.204  2 0.903 0.903



#############
# Three group model where the latent variable parameters are constrained to
# be equal in the focal groups

set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
itemtype <- rep('2PL', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
dataset3 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2, dataset3)
group <- rep(c('D1', 'D2', 'D3'), each=N)

# marginal information
itemstats(dat)
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  3000             7.89          3.567  0.19 0.054 0.779     1.676
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  3000 2 0.605 0.489   0.525         0.415       0.764
#> Item_2  3000 2 0.403 0.491   0.563         0.458       0.761
#> Item_3  3000 2 0.457 0.498   0.504         0.388       0.767
#> Item_4  3000 2 0.676 0.468   0.458         0.345       0.770
#> Item_5  3000 2 0.545 0.498   0.539         0.429       0.763
#> Item_6  3000 2 0.645 0.478   0.334         0.207       0.782
#> Item_7  3000 2 0.688 0.464   0.529         0.426       0.764
#> Item_8  3000 2 0.427 0.495   0.498         0.382       0.767
#> Item_9  3000 2 0.292 0.455   0.454         0.344       0.770
#> Item_10 3000 2 0.284 0.451   0.393         0.279       0.775
#> Item_11 3000 2 0.738 0.440   0.451         0.345       0.770
#> Item_12 3000 2 0.460 0.499   0.598         0.496       0.757
#> Item_13 3000 2 0.591 0.492   0.556         0.449       0.761
#> Item_14 3000 2 0.591 0.492   0.548         0.441       0.762
#> Item_15 3000 2 0.488 0.500   0.452         0.330       0.772
#> 
#> $proportions
#>             0     1
#> Item_1  0.395 0.605
#> Item_2  0.597 0.403
#> Item_3  0.543 0.457
#> Item_4  0.324 0.676
#> Item_5  0.455 0.545
#> Item_6  0.355 0.645
#> Item_7  0.312 0.688
#> Item_8  0.573 0.427
#> Item_9  0.708 0.292
#> Item_10 0.716 0.284
#> Item_11 0.262 0.738
#> Item_12 0.540 0.460
#> Item_13 0.409 0.591
#> Item_14 0.409 0.591
#> Item_15 0.512 0.488
#> 

# conditional information
itemstats(dat, group=group)
#> $D1
#> $D1$overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  1000             7.82          3.346 0.159 0.047  0.74     1.705
#> 
#> $D1$itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  1000 2 0.605 0.489   0.484         0.361       0.725
#> Item_2  1000 2 0.368 0.483   0.507         0.388       0.722
#> Item_3  1000 2 0.461 0.499   0.471         0.343       0.727
#> Item_4  1000 2 0.673 0.469   0.439         0.315       0.729
#> Item_5  1000 2 0.526 0.500   0.500         0.376       0.723
#> Item_6  1000 2 0.654 0.476   0.355         0.222       0.739
#> Item_7  1000 2 0.683 0.466   0.508         0.394       0.722
#> Item_8  1000 2 0.431 0.495   0.462         0.333       0.728
#> Item_9  1000 2 0.287 0.453   0.419         0.298       0.731
#> Item_10 1000 2 0.274 0.446   0.372         0.249       0.735
#> Item_11 1000 2 0.739 0.439   0.401         0.282       0.732
#> Item_12 1000 2 0.456 0.498   0.563         0.448       0.715
#> Item_13 1000 2 0.584 0.493   0.535         0.416       0.719
#> Item_14 1000 2 0.592 0.492   0.483         0.358       0.725
#> Item_15 1000 2 0.487 0.500   0.451         0.321       0.729
#> 
#> $D1$proportions
#>             0     1
#> Item_1  0.395 0.605
#> Item_2  0.632 0.368
#> Item_3  0.539 0.461
#> Item_4  0.327 0.673
#> Item_5  0.474 0.526
#> Item_6  0.346 0.654
#> Item_7  0.317 0.683
#> Item_8  0.569 0.431
#> Item_9  0.713 0.287
#> Item_10 0.726 0.274
#> Item_11 0.261 0.739
#> Item_12 0.544 0.456
#> Item_13 0.416 0.584
#> Item_14 0.408 0.592
#> Item_15 0.513 0.487
#> 
#> 
#> $D2
#> $D2$overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  1000            7.955          3.754 0.217 0.065 0.807      1.65
#> 
#> $D2$itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  1000 2 0.612 0.488   0.563         0.464       0.792
#> Item_2  1000 2 0.417 0.493   0.569         0.469       0.792
#> Item_3  1000 2 0.450 0.498   0.578         0.479       0.791
#> Item_4  1000 2 0.695 0.461   0.484         0.381       0.798
#> Item_5  1000 2 0.551 0.498   0.553         0.451       0.793
#> Item_6  1000 2 0.644 0.479   0.316         0.195       0.811
#> Item_7  1000 2 0.680 0.467   0.540         0.443       0.794
#> Item_8  1000 2 0.432 0.496   0.544         0.440       0.794
#> Item_9  1000 2 0.317 0.466   0.470         0.365       0.799
#> Item_10 1000 2 0.275 0.447   0.403         0.297       0.804
#> Item_11 1000 2 0.729 0.445   0.508         0.413       0.796
#> Item_12 1000 2 0.483 0.500   0.606         0.511       0.788
#> Item_13 1000 2 0.585 0.493   0.587         0.491       0.790
#> Item_14 1000 2 0.591 0.492   0.589         0.493       0.790
#> Item_15 1000 2 0.494 0.500   0.466         0.351       0.801
#> 
#> $D2$proportions
#>             0     1
#> Item_1  0.388 0.612
#> Item_2  0.583 0.417
#> Item_3  0.550 0.450
#> Item_4  0.305 0.695
#> Item_5  0.449 0.551
#> Item_6  0.356 0.644
#> Item_7  0.320 0.680
#> Item_8  0.568 0.432
#> Item_9  0.683 0.317
#> Item_10 0.725 0.275
#> Item_11 0.271 0.729
#> Item_12 0.517 0.483
#> Item_13 0.415 0.585
#> Item_14 0.409 0.591
#> Item_15 0.506 0.494
#> 
#> 
#> $D3
#> $D3$overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  1000            7.894          3.592 0.194 0.064 0.783     1.672
#> 
#> $D3$itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  1000 2 0.599 0.490   0.526         0.416       0.769
#> Item_2  1000 2 0.424 0.494   0.610         0.512       0.761
#> Item_3  1000 2 0.459 0.499   0.460         0.340       0.776
#> Item_4  1000 2 0.659 0.474   0.452         0.338       0.775
#> Item_5  1000 2 0.557 0.497   0.562         0.456       0.766
#> Item_6  1000 2 0.638 0.481   0.334         0.208       0.786
#> Item_7  1000 2 0.700 0.458   0.539         0.439       0.767
#> Item_8  1000 2 0.419 0.494   0.485         0.369       0.773
#> Item_9  1000 2 0.272 0.445   0.470         0.365       0.773
#> Item_10 1000 2 0.303 0.460   0.405         0.290       0.779
#> Item_11 1000 2 0.747 0.435   0.438         0.333       0.776
#> Item_12 1000 2 0.442 0.497   0.622         0.526       0.759
#> Item_13 1000 2 0.604 0.489   0.545         0.438       0.767
#> Item_14 1000 2 0.589 0.492   0.569         0.465       0.765
#> Item_15 1000 2 0.482 0.500   0.439         0.317       0.778
#> 
#> $D3$proportions
#>             0     1
#> Item_1  0.401 0.599
#> Item_2  0.576 0.424
#> Item_3  0.541 0.459
#> Item_4  0.341 0.659
#> Item_5  0.443 0.557
#> Item_6  0.362 0.638
#> Item_7  0.300 0.700
#> Item_8  0.581 0.419
#> Item_9  0.728 0.272
#> Item_10 0.697 0.303
#> Item_11 0.253 0.747
#> Item_12 0.558 0.442
#> Item_13 0.396 0.604
#> Item_14 0.411 0.589
#> Item_15 0.518 0.482
#> 
#> 

model <- 'F1 = 1-15
          FREE[D2, D3] = (GROUP, MEAN_1), (GROUP, COV_11)
          CONSTRAINB[D2,D3] = (GROUP, MEAN_1), (GROUP, COV_11)'

mod <- multipleGroup(dat, model, group = group, invariance = colnames(dat))
coef(mod, simplify=TRUE)
#> $D1
#> $items
#>            a1      d g u
#> Item_1  1.089  0.517 0 1
#> Item_2  1.302 -0.601 0 1
#> Item_3  0.950 -0.251 0 1
#> Item_4  0.857  0.847 0 1
#> Item_5  1.119  0.195 0 1
#> Item_6  0.447  0.618 0 1
#> Item_7  1.209  1.023 0 1
#> Item_8  0.939 -0.396 0 1
#> Item_9  0.909 -1.110 0 1
#> Item_10 0.695 -1.073 0 1
#> Item_11 0.929  1.232 0 1
#> Item_12 1.454 -0.290 0 1
#> Item_13 1.216  0.457 0 1
#> Item_14 1.199  0.453 0 1
#> Item_15 0.751 -0.085 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> 
#> $D2
#> $items
#>            a1      d g u
#> Item_1  1.089  0.517 0 1
#> Item_2  1.302 -0.601 0 1
#> Item_3  0.950 -0.251 0 1
#> Item_4  0.857  0.847 0 1
#> Item_5  1.119  0.195 0 1
#> Item_6  0.447  0.618 0 1
#> Item_7  1.209  1.023 0 1
#> Item_8  0.939 -0.396 0 1
#> Item_9  0.909 -1.110 0 1
#> Item_10 0.695 -1.073 0 1
#> Item_11 0.929  1.232 0 1
#> Item_12 1.454 -0.290 0 1
#> Item_13 1.216  0.457 0 1
#> Item_14 1.199  0.453 0 1
#> Item_15 0.751 -0.085 0 1
#> 
#> $means
#>    F1 
#> 0.055 
#> 
#> $cov
#>       F1
#> F1 1.467
#> 
#> 
#> $D3
#> $items
#>            a1      d g u
#> Item_1  1.089  0.517 0 1
#> Item_2  1.302 -0.601 0 1
#> Item_3  0.950 -0.251 0 1
#> Item_4  0.857  0.847 0 1
#> Item_5  1.119  0.195 0 1
#> Item_6  0.447  0.618 0 1
#> Item_7  1.209  1.023 0 1
#> Item_8  0.939 -0.396 0 1
#> Item_9  0.909 -1.110 0 1
#> Item_10 0.695 -1.073 0 1
#> Item_11 0.929  1.232 0 1
#> Item_12 1.454 -0.290 0 1
#> Item_13 1.216  0.457 0 1
#> Item_14 1.199  0.453 0 1
#> Item_15 0.751 -0.085 0 1
#> 
#> $means
#>    F1 
#> 0.055 
#> 
#> $cov
#>       F1
#> F1 1.467
#> 
#> 

#############
# Testing main effects in multiple independent grouping variables
set.seed(1234)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
itemtype <- rep('2PL', nrow(a))
N <- 500

# generated data have interaction effect for latent means, as well as a
# main effect across D but no main effect across G
d11 <- simdata(a, d, N, itemtype, mu = 0)
d12 <- simdata(a, d, N, itemtype, mu = 0)
d13 <- simdata(a, d, N, itemtype, mu = 0)
d21 <- simdata(a, d, N, itemtype, mu = 1/2)
d22 <- simdata(a, d, N, itemtype, mu = 1/2)
d23 <- simdata(a, d, N, itemtype, mu = -1)
dat <- do.call(rbind, list(d11, d12, d13, d21, d22, d23))
group <- rep(c('G1.D1', 'G1.D2', 'G1.D3', 'G2.D1', 'G2.D2', 'G2.D3'), each=N)
table(group)
#> group
#> G1.D1 G1.D2 G1.D3 G2.D1 G2.D2 G2.D3 
#>   500   500   500   500   500   500 

if (FALSE) { # \dontrun{
# in practice, group would be organized in a data.frame as follows
df <- data.frame(group)
dfw <- tidyr::separate_wider_delim(df, group, delim='.', names = c('G', 'D'))
head(dfw)

# for use with multipleGroup() combine into a single long group
group <- with(dfw, factor(G):factor(D))

# conditional information
itemstats(dat, group=group)

mod <- multipleGroup(dat, group = group, SE=TRUE,
                     invariance = c(colnames(dat), 'free_mean', 'free_var'))
coef(mod, simplify=TRUE)
sapply(coef(mod, simplify=TRUE), \(x) unname(x$means)) # mean estimates
wald(mod) # view parameter names for later testing

# test for main effect over G group (manually compute marginal mean)
wald(mod, "0 + MEAN_1.123 + MEAN_1.185 = MEAN_1.247 + MEAN_1.309 + MEAN_1.371")

# test for main effect over D group  (manually compute marginal means)
wald(mod, c("0 + MEAN_1.247 = MEAN_1.123 + MEAN_1.309",
            "0 + MEAN_1.247 = MEAN_1.185 + MEAN_1.371"))

# post-hoc tests (better practice would include p.adjust() )
wald(mod, "0 + MEAN_1.247 = MEAN_1.123 + MEAN_1.309") # D1 vs D2
wald(mod, "0 + MEAN_1.247 = MEAN_1.185 + MEAN_1.371") # D1 vs D3
wald(mod, "MEAN_1.123 + MEAN_1.309 = MEAN_1.185 + MEAN_1.371") # D2 vs D3
} # }

#############
# multiple factors

a <- matrix(c(abs(rnorm(5,1,.3)), rep(0,15),abs(rnorm(5,1,.3)),
     rep(0,15),abs(rnorm(5,1,.3))), 15, 3)
d <- matrix(rnorm(15,0,.7),ncol=1)
mu <- c(-.4, -.7, .1)
sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)
itemtype <- rep('2PL', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = mu, sigma = sigma)
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))

# group models
model <- '
   F1 = 1-5
   F2 = 6-10
   F3 = 11-15'

# define mirt cluster to use parallel architecture
if(interactive()) mirtCluster()

# EM approach (not as accurate with 3 factors, but generally good for quick model comparisons)
mod_configural <- multipleGroup(dat, model, group = group) #completely separate analyses
mod_metric <- multipleGroup(dat, model, group = group, invariance=c('slopes')) #equal slopes
mod_fullconstrain <- multipleGroup(dat, model, group = group, #equal means, slopes, intercepts
                             invariance=c('slopes', 'intercepts'))

anova(mod_metric, mod_configural)
#>                     AIC    SABIC       HQ      BIC    logLik     X2 df     p
#> mod_metric     36921.29 37030.37 37013.84 37173.33 -18415.65                
#> mod_configural 36916.10 37061.53 37039.49 37252.15 -18398.05 35.197 15 0.002
anova(mod_fullconstrain, mod_metric)
#>                        AIC    SABIC       HQ      BIC    logLik      X2 df p
#> mod_fullconstrain 37020.58 37093.29 37082.27 37188.60 -18480.29             
#> mod_metric        36921.29 37030.37 37013.84 37173.33 -18415.65 129.284 15 0

# same as above, but with MHRM (generally  more accurate with 3+ factors, but slower)
mod_configural <- multipleGroup(dat, model, group = group, method = 'MHRM')
mod_metric <- multipleGroup(dat, model, group = group, invariance=c('slopes'), method = 'MHRM')
mod_fullconstrain <- multipleGroup(dat, model, group = group, method = 'MHRM',
                             invariance=c('slopes', 'intercepts'))

anova(mod_metric, mod_configural)
#>                     AIC    SABIC       HQ      BIC    logLik     X2 df p
#> mod_metric     36927.79 37036.86 37020.33 37179.83 -18418.89            
#> mod_configural 36918.04 37063.47 37041.43 37254.10 -18399.02 39.745 15 0
anova(mod_fullconstrain, mod_metric)
#>                        AIC    SABIC       HQ      BIC    logLik     X2 df p
#> mod_fullconstrain 37023.82 37096.54 37085.52 37191.85 -18481.91            
#> mod_metric        36927.79 37036.86 37020.33 37179.83 -18418.89 126.04 15 0

############
# polytomous item example
set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
d <- cbind(d, d-1, d-2)
itemtype <- rep('graded', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))
model <- 'F1 = 1-15'

mod_configural <- multipleGroup(dat, model, group = group)
plot(mod_configural)

plot(mod_configural, type = 'SE')

plot(mod_configural, type = 'gen.difficulty')

itemplot(mod_configural, 1)

itemplot(mod_configural, 1, type = 'info')

plot(mod_configural, type = 'trace') # messy, score function typically better

plot(mod_configural, type = 'itemscore')


fs <- fscores(mod_configural, full.scores = FALSE)
head(fs[["D1"]])
#>      Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> [1,]      0      0      0      0      0      0      0      0      0       0
#> [2,]      0      0      0      0      0      0      0      0      0       0
#> [3,]      0      0      0      0      0      0      0      0      0       0
#> [4,]      0      0      0      0      0      0      0      0      0       0
#> [5,]      0      0      0      0      0      0      0      0      0       0
#> [6,]      0      0      0      0      0      0      0      0      0       0
#>      Item_11 Item_12 Item_13 Item_14 Item_15        F1     SE_F1
#> [1,]       0       0       0       0       0 -2.137056 0.6262714
#> [2,]       0       0       0       1       0 -1.785692 0.5741235
#> [3,]       0       0       0       2       0 -1.726621 0.5800312
#> [4,]       0       0       0       3       0 -1.688716 0.5856268
#> [5,]       0       0       0       3       1 -1.439230 0.5535499
#> [6,]       0       2       0       0       0 -1.582715 0.5660598
fscores(mod_configural, method = 'EAPsum', full.scores = FALSE)
#> $D1
#>    Sum.Scores           F1     SE_F1 observed    expected    std.res
#> 0           0 -2.137056441 0.6262714        6  5.73109332 0.11232666
#> 1           1 -1.842264885 0.5888524        6 10.08249657 1.28570635
#> 2           2 -1.645473280 0.5750881       19 15.02644053 1.02506560
#> 3           3 -1.490721889 0.5681733       22 21.11696057 0.19216062
#> 4           4 -1.322711281 0.5477826       31 26.16731435 0.94473275
#> 5           5 -1.172604830 0.5343878       21 30.81114992 1.76752660
#> 6           6 -1.033369339 0.5229023       30 34.97351880 0.84099633
#> 7           7 -0.898167173 0.5101253       37 38.28642385 0.20790344
#> 8           8 -0.770379457 0.4995061       47 40.92147882 0.95021586
#> 9           9 -0.648276432 0.4899687       40 42.89556180 0.44210625
#> 10         10 -0.530510563 0.4810482       54 44.20825137 1.47268067
#> 11         11 -0.417091373 0.4732401       52 44.94497192 1.05234512
#> 12         12 -0.307324568 0.4662902       50 45.16103718 0.72006273
#> 13         13 -0.200657261 0.4600924       35 44.90907262 1.47865198
#> 14         14 -0.096765731 0.4546839       41 44.25218143 0.48888580
#> 15         15  0.004759979 0.4499792       33 43.24436612 1.55783093
#> 16         16  0.104274150 0.4459305       50 41.93576821 1.24529057
#> 17         17  0.202071875 0.4425175       38 40.37507441 0.37378387
#> 18         18  0.298449688 0.4397012       34 38.60588650 0.74128714
#> 19         19  0.393678038 0.4374552       40 36.66805457 0.55024226
#> 20         20  0.488005095 0.4357619       30 34.59858723 0.78179924
#> 21         21  0.581672874 0.4346027       42 32.43092044 1.68031432
#> 22         22  0.674910103 0.4339646       32 30.19547738 0.32839110
#> 23         23  0.767936430 0.4338392       28 27.92026367 0.01509025
#> 24         24  0.860968350 0.4342199       19 25.63069683 1.30972152
#> 25         25  0.954215436 0.4351028       34 23.34994032 2.20398764
#> 26         26  1.047885949 0.4364871       21 21.09940668 0.02164116
#> 27         27  1.142189172 0.4383746       18 18.89862601 0.20671113
#> 28         28  1.237330487 0.4407665       21 16.76546358 1.03418465
#> 29         29  1.333523036 0.4436684       10 14.71661806 1.22949481
#> 30         30  1.430986181 0.4470894       11 12.76726484 0.49459840
#> 31         31  1.529934577 0.4510309        8 10.93114667 0.88655297
#> 32         32  1.630612066 0.4555070        8  9.22112463 0.40213145
#> 33         33  1.733275013 0.4605369        8  7.64818473 0.12721417
#> 34         34  1.838159660 0.4661076        6  6.22141174 0.08876795
#> 35         35  1.945607770 0.4722588        1  4.94877914 1.77506315
#> 36         36  2.055967060 0.4790331        4  3.83475999 0.08438128
#> 37         37  2.169474066 0.4863361        6  2.88068295 1.83785730
#> 38         38  2.286835090 0.4943335        3  2.08640055 0.63249484
#> 39         39  2.408610539 0.5031395        3  1.44515936 1.29338558
#> 40         40  2.534716475 0.5122126        1  0.94647470 0.05501798
#> 41         41  2.667837385 0.5224881        0  0.57961995 0.76132775
#> 42         42  2.809124175 0.5344045        0  0.32363182 0.56888648
#> 43         43  2.954701793 0.5443149        0  0.15780842 0.39725108
#> 44         44  3.126041758 0.5588274        0  0.06566835 0.25625837
#> 45         45  3.344792190 0.5841085        0  0.01877911 0.13703691
#> 
#> $D2
#>    Sum.Scores          F1     SE_F1 observed   expected    std.res
#> 0           0 -2.05130147 0.5912971        9 11.0355976 0.61276507
#> 1           1 -1.75022208 0.5476302       23 16.1395895 1.70766964
#> 2           2 -1.55419674 0.5300849       22 20.9952745 0.21927375
#> 3           3 -1.41173370 0.5243664       31 27.0499126 0.75949299
#> 4           4 -1.24022266 0.4925574       35 30.4339816 0.82767202
#> 5           5 -1.09460723 0.4745106       33 33.1973941 0.03425959
#> 6           6 -0.96352775 0.4602807       34 35.4757938 0.24777621
#> 7           7 -0.83622876 0.4435206       32 36.8803015 0.80361756
#> 8           8 -0.71831303 0.4305401       32 37.8120865 0.94518484
#> 9           9 -0.60733281 0.4193943       31 38.3526137 1.18725571
#> 10         10 -0.50128231 0.4091029       28 38.5164822 1.69452188
#> 11         11 -0.40019407 0.4004018       46 38.4037717 1.22577588
#> 12         12 -0.30319834 0.3928821       41 38.0627585 0.47609052
#> 13         13 -0.20956698 0.3863159       36 37.5255576 0.24903767
#> 14         14 -0.11890072 0.3807251       33 36.8301475 0.63112266
#> 15         15 -0.03072312 0.3759795       37 36.0028603 0.16618336
#> 16         16  0.05538372 0.3719906       33 35.0638735 0.34854040
#> 17         17  0.13975529 0.3687182       33 34.0317999 0.17686952
#> 18         18  0.22271790 0.3661063       42 32.9212969 1.58228730
#> 19         19  0.30456977 0.3641136       28 31.7441592 0.66454193
#> 20         20  0.38557952 0.3627146       29 30.5107753 0.27351003
#> 21         21  0.46600618 0.3618867       30 29.2298023 0.14245880
#> 22         22  0.54609598 0.3616155       24 27.9084016 0.73982959
#> 23         23  0.62608453 0.3618947       28 26.5529477 0.28081988
#> 24         24  0.70620669 0.3627244       29 25.1690033 0.76362260
#> 25         25  0.78669301 0.3641096       27 23.7612468 0.66442083
#> 26         26  0.86777451 0.3660611       29 22.3340727 1.41051144
#> 27         27  0.94969277 0.3686001       19 20.8916388 0.41385858
#> 28         28  1.03268637 0.3717461       17 19.4374534 0.55286188
#> 29         29  1.11700520 0.3755241       13 17.9753823 1.17351161
#> 30         30  1.20293125 0.3799843       19 16.5097589 0.61287376
#> 31         31  1.29071727 0.3851447       17 15.0440068 0.50429607
#> 32         32  1.38065471 0.3910380       17 13.5834528 0.92700555
#> 33         33  1.47313321 0.3977848       16 12.1351562 1.10945390
#> 34         34  1.56833794 0.4053106       11 10.7031592 0.09073350
#> 35         35  1.66666006 0.4136580        8  9.2985852 0.42585525
#> 36         36  1.76883244 0.4232093       11  7.9356798 1.08778253
#> 37         37  1.87438760 0.4333202        3  6.6167296 1.40602924
#> 38         38  1.98429588 0.4442628        6  5.3697725 0.27196902
#> 39         39  2.10060111 0.4574333        4  4.2190260 0.10663252
#> 40         40  2.21882242 0.4694234        2  3.1515306 0.64865644
#> 41         41  2.34438672 0.4824723        2  2.2364769 0.15812718
#> 42         42  2.48431289 0.5010802        0  1.4847938 1.21852115
#> 43         43  2.61160839 0.5110564        0  0.8421851 0.91770643
#> 44         44  2.76572645 0.5249224        0  0.4383660 0.66209216
#> 45         45  2.98799042 0.5545567        0  0.1853438 0.43051577
#> 

# constrain slopes within each group to be equal (but not across groups)
model2 <- 'F1 = 1-15
           CONSTRAIN = (1-15, a1)'
mod_configural2 <- multipleGroup(dat, model2, group = group)
plot(mod_configural2, type = 'SE')

plot(mod_configural2, type = 'RE')

itemplot(mod_configural2, 10)


############
## empirical histogram example (normal and bimodal groups)
set.seed(1234)
a <- matrix(rlnorm(50, .2, .2))
d <- matrix(rnorm(50))
ThetaNormal <- matrix(rnorm(2000))
ThetaBimodal <- scale(matrix(c(rnorm(1000, -2), rnorm(1000,2)))) #bimodal
Theta <- rbind(ThetaNormal, ThetaBimodal)
dat <- simdata(a, d, 4000, itemtype = '2PL', Theta=Theta)
group <- rep(c('G1', 'G2'), each=2000)

EH <- multipleGroup(dat, 1, group=group, dentype="empiricalhist", invariance = colnames(dat))
coef(EH, simplify=TRUE)
#> $G1
#> $items
#>            a1      d g u
#> Item_1  0.759 -1.830 0 1
#> Item_2  1.062 -0.555 0 1
#> Item_3  1.235 -1.117 0 1
#> Item_4  0.611 -0.953 0 1
#> Item_5  1.088 -0.132 0 1
#> Item_6  1.054  0.502 0 1
#> Item_7  0.849  1.581 0 1
#> Item_8  0.870 -0.779 0 1
#> Item_9  0.902  1.658 0 1
#> Item_10 0.812 -1.188 0 1
#> Item_11 0.903  0.661 0 1
#> Item_12 0.841  2.646 0 1
#> Item_13 0.843 -0.064 0 1
#> Item_14 1.009 -0.600 0 1
#> Item_15 1.167  0.013 0 1
#> Item_16 1.008  1.886 0 1
#> Item_17 0.885 -1.124 0 1
#> Item_18 0.790  1.420 0 1
#> Item_19 0.799  1.316 0 1
#> Item_20 1.565  0.294 0 1
#> Item_21 1.002  0.001 0 1
#> Item_22 0.929 -0.497 0 1
#> Item_23 0.903 -0.383 0 1
#> Item_24 1.120  0.702 0 1
#> Item_25 0.850  2.002 0 1
#> Item_26 0.706 -0.208 0 1
#> Item_27 1.048 -1.279 0 1
#> Item_28 0.826 -0.770 0 1
#> Item_29 1.062  0.221 0 1
#> Item_30 0.737 -0.333 0 1
#> Item_31 1.208 -0.178 0 1
#> Item_32 0.917 -0.202 0 1
#> Item_33 0.859 -1.360 0 1
#> Item_34 0.913 -0.252 0 1
#> Item_35 0.767  0.886 0 1
#> Item_36 0.783  0.763 0 1
#> Item_37 0.630  0.550 0 1
#> Item_38 0.736 -0.405 0 1
#> Item_39 0.917 -0.166 0 1
#> Item_40 0.895 -1.251 0 1
#> Item_41 1.251 -0.040 0 1
#> Item_42 0.882  0.202 0 1
#> Item_43 0.796  1.639 0 1
#> Item_44 1.022  1.038 0 1
#> Item_45 0.812 -0.470 0 1
#> Item_46 0.807  0.336 0 1
#> Item_47 0.815 -1.121 0 1
#> Item_48 0.739  0.892 0 1
#> Item_49 0.904  1.036 0 1
#> Item_50 0.895  2.124 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> 
#> $G2
#> $items
#>            a1      d g u
#> Item_1  0.759 -1.830 0 1
#> Item_2  1.062 -0.555 0 1
#> Item_3  1.235 -1.117 0 1
#> Item_4  0.611 -0.953 0 1
#> Item_5  1.088 -0.132 0 1
#> Item_6  1.054  0.502 0 1
#> Item_7  0.849  1.581 0 1
#> Item_8  0.870 -0.779 0 1
#> Item_9  0.902  1.658 0 1
#> Item_10 0.812 -1.188 0 1
#> Item_11 0.903  0.661 0 1
#> Item_12 0.841  2.646 0 1
#> Item_13 0.843 -0.064 0 1
#> Item_14 1.009 -0.600 0 1
#> Item_15 1.167  0.013 0 1
#> Item_16 1.008  1.886 0 1
#> Item_17 0.885 -1.124 0 1
#> Item_18 0.790  1.420 0 1
#> Item_19 0.799  1.316 0 1
#> Item_20 1.565  0.294 0 1
#> Item_21 1.002  0.001 0 1
#> Item_22 0.929 -0.497 0 1
#> Item_23 0.903 -0.383 0 1
#> Item_24 1.120  0.702 0 1
#> Item_25 0.850  2.002 0 1
#> Item_26 0.706 -0.208 0 1
#> Item_27 1.048 -1.279 0 1
#> Item_28 0.826 -0.770 0 1
#> Item_29 1.062  0.221 0 1
#> Item_30 0.737 -0.333 0 1
#> Item_31 1.208 -0.178 0 1
#> Item_32 0.917 -0.202 0 1
#> Item_33 0.859 -1.360 0 1
#> Item_34 0.913 -0.252 0 1
#> Item_35 0.767  0.886 0 1
#> Item_36 0.783  0.763 0 1
#> Item_37 0.630  0.550 0 1
#> Item_38 0.736 -0.405 0 1
#> Item_39 0.917 -0.166 0 1
#> Item_40 0.895 -1.251 0 1
#> Item_41 1.251 -0.040 0 1
#> Item_42 0.882  0.202 0 1
#> Item_43 0.796  1.639 0 1
#> Item_44 1.022  1.038 0 1
#> Item_45 0.812 -0.470 0 1
#> Item_46 0.807  0.336 0 1
#> Item_47 0.815 -1.121 0 1
#> Item_48 0.739  0.892 0 1
#> Item_49 0.904  1.036 0 1
#> Item_50 0.895  2.124 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> 
plot(EH, type = 'empiricalhist', npts = 60)


# DIF test for item 1
EH1 <- multipleGroup(dat, 1, group=group, dentype="empiricalhist", invariance = colnames(dat)[-1])
anova(EH, EH1)
#>          AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> EH  216599.8 217659.4 217358.4 218739.8 -107959.9               
#> EH1 216602.9 217668.7 217365.9 218755.4 -107959.4 0.916  2 0.632

#--------------------------------
# Mixture model (no prior group variable specified)

set.seed(12345)
nitems <- 20
a1 <- matrix(.75, ncol=1, nrow=nitems)
a2 <- matrix(1.25, ncol=1, nrow=nitems)
d1 <- matrix(rnorm(nitems,0,1),ncol=1)
d2 <- matrix(rnorm(nitems,0,1),ncol=1)
itemtype <- rep('2PL', nrow(a1))
N1 <- 500
N2 <- N1*2 # second class twice as large

dataset1 <- simdata(a1, d1, N1, itemtype)
dataset2 <- simdata(a2, d2, N2, itemtype)
dat <- rbind(dataset1, dataset2)
# group <- c(rep('D1', N1), rep('D2', N2))

# Mixture Rasch model (Rost, 1990)
models <- 'F1 = 1-20
           CONSTRAIN = (1-20, a1)'
mod_mix <- multipleGroup(dat, models, dentype = 'mixture-2', GenRandomPars = TRUE)
coef(mod_mix, simplify=TRUE)
#> $MIXTURE_1
#> $items
#>            a1      d g u
#> Item_1  1.302  0.772 0 1
#> Item_2  1.302  1.476 0 1
#> Item_3  1.302 -0.553 0 1
#> Item_4  1.302 -1.588 0 1
#> Item_5  1.302 -1.545 0 1
#> Item_6  1.302  2.084 0 1
#> Item_7  1.302 -0.510 0 1
#> Item_8  1.302  0.606 0 1
#> Item_9  1.302  0.550 0 1
#> Item_10 1.302 -0.167 0 1
#> Item_11 1.302  0.961 0 1
#> Item_12 1.302  2.156 0 1
#> Item_13 1.302  2.082 0 1
#> Item_14 1.302  1.733 0 1
#> Item_15 1.302  0.138 0 1
#> Item_16 1.302  0.434 0 1
#> Item_17 1.302 -0.368 0 1
#> Item_18 1.302 -1.633 0 1
#> Item_19 1.302  1.799 0 1
#> Item_20 1.302 -0.024 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> $class_proportion
#>     pi
#>  0.655
#> 
#> 
#> $MIXTURE_2
#> $items
#>            a1      d g u
#> Item_1  0.681  0.664 0 1
#> Item_2  0.681  0.805 0 1
#> Item_3  0.681 -0.039 0 1
#> Item_4  0.681 -0.416 0 1
#> Item_5  0.681  0.530 0 1
#> Item_6  0.681 -2.045 0 1
#> Item_7  0.681  0.770 0 1
#> Item_8  0.681 -0.314 0 1
#> Item_9  0.681 -0.030 0 1
#> Item_10 0.681 -0.820 0 1
#> Item_11 0.681 -0.252 0 1
#> Item_12 0.681  1.939 0 1
#> Item_13 0.681  0.642 0 1
#> Item_14 0.681  0.593 0 1
#> Item_15 0.681 -0.749 0 1
#> Item_16 0.681  0.852 0 1
#> Item_17 0.681 -0.950 0 1
#> Item_18 0.681 -0.345 0 1
#> Item_19 0.681  1.230 0 1
#> Item_20 0.681  0.271 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> $class_proportion
#>     pi
#>  0.345
#> 
#> 
summary(mod_mix)
#> 
#> ----------
#> GROUP: MIXTURE_1 
#>          F1    h2
#>  [1,] 0.608 0.369
#>  [2,] 0.608 0.369
#>  [3,] 0.608 0.369
#>  [4,] 0.608 0.369
#>  [5,] 0.608 0.369
#>  [6,] 0.608 0.369
#>  [7,] 0.608 0.369
#>  [8,] 0.608 0.369
#>  [9,] 0.608 0.369
#> [10,] 0.608 0.369
#> [11,] 0.608 0.369
#> [12,] 0.608 0.369
#> [13,] 0.608 0.369
#> [14,] 0.608 0.369
#> [15,] 0.608 0.369
#> [16,] 0.608 0.369
#> [17,] 0.608 0.369
#> [18,] 0.608 0.369
#> [19,] 0.608 0.369
#> [20,] 0.608 0.369
#> 
#> SS loadings:  7.386 
#> Proportion Var:  0.369 
#> 
#> Factor correlations: 
#> 
#>    F1
#> F1  1
#> 
#> Class proportion:  0.655 
#> 
#> ----------
#> GROUP: MIXTURE_2 
#>          F1    h2
#>  [1,] 0.371 0.138
#>  [2,] 0.371 0.138
#>  [3,] 0.371 0.138
#>  [4,] 0.371 0.138
#>  [5,] 0.371 0.138
#>  [6,] 0.371 0.138
#>  [7,] 0.371 0.138
#>  [8,] 0.371 0.138
#>  [9,] 0.371 0.138
#> [10,] 0.371 0.138
#> [11,] 0.371 0.138
#> [12,] 0.371 0.138
#> [13,] 0.371 0.138
#> [14,] 0.371 0.138
#> [15,] 0.371 0.138
#> [16,] 0.371 0.138
#> [17,] 0.371 0.138
#> [18,] 0.371 0.138
#> [19,] 0.371 0.138
#> [20,] 0.371 0.138
#> 
#> SS loadings:  2.758 
#> Proportion Var:  0.138 
#> 
#> Factor correlations: 
#> 
#>    F1
#> F1  1
#> 
#> Class proportion:  0.345 
plot(mod_mix)

plot(mod_mix, type = 'trace')

plot(mod_mix, type = 'gen.difficulty')

itemplot(mod_mix, 1, type = 'info')


head(fscores(mod_mix)) # theta estimates
#>         Class_1
#> [1,] -0.3609178
#> [2,] -1.6950807
#> [3,] -0.8217168
#> [4,]  0.8090483
#> [5,]  0.7936363
#> [6,]  0.2286020
head(fscores(mod_mix, method = 'classify')) # classification probability
#>         CLASS_1   CLASS_2
#> [1,] 0.04826707 0.9517329
#> [2,] 0.74223456 0.2577654
#> [3,] 0.02496631 0.9750337
#> [4,] 0.02929597 0.9707040
#> [5,] 0.07220958 0.9277904
#> [6,] 0.42428417 0.5757158
itemfit(mod_mix)
#>       item   S_X2 df.S_X2 RMSEA.S_X2 p.S_X2
#> 1   Item_1 30.498      12      0.032  0.002
#> 2   Item_2 36.271      12      0.037  0.000
#> 3   Item_3 15.756      12      0.014  0.203
#> 4   Item_4  9.334      11      0.000  0.591
#> 5   Item_5 21.638      12      0.023  0.042
#> 6   Item_6 26.071      13      0.026  0.017
#> 7   Item_7  8.624      12      0.000  0.735
#> 8   Item_8 16.652      12      0.016  0.163
#> 9   Item_9 22.772      12      0.024  0.030
#> 10 Item_10 14.622      12      0.012  0.263
#> 11 Item_11  7.357      12      0.000  0.833
#> 12 Item_12 14.964      12      0.013  0.243
#> 13 Item_13 10.222      12      0.000  0.597
#> 14 Item_14 17.414      13      0.015  0.181
#> 15 Item_15  9.436      12      0.000  0.665
#> 16 Item_16 11.408      12      0.000  0.494
#> 17 Item_17 16.483      12      0.016  0.170
#> 18 Item_18 11.717      11      0.007  0.385
#> 19 Item_19 11.425      12      0.000  0.493
#> 20 Item_20 15.291      12      0.014  0.226

# Above works fine, but its generally a good idea to evaluate models
# with multiple random starting values in case local maximums are an issue.
# To do this use the argument "nruns", which returns the best model
if(interactive()) mirtCluster()
mod_mix <- multipleGroup(dat, models, dentype = 'mixture-2', nruns=5)
mod_mix
#> 
#> Call:
#> multipleGroup(data = dat, model = models, dentype = "mixture-2", 
#>     nruns = 5)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within 1e-04 tolerance after 58 EM iterations.
#> mirt version: 1.45.6 
#> M-step optimizer: BFGS 
#> EM acceleration: Ramsay 
#> Number of rectangular quadrature: 61
#> Latent density type: mixture 
#> 
#> Log-likelihood = -17271.09
#> Estimated parameters: 43 
#> AIC = 34628.19
#> BIC = 34856.66; SABIC = 34720.06
#> G2 (1048532) = 12931.52, p = 1
#> RMSEA = 0, CFI = NaN, TLI = NaN

# For obtaining isolated estimates within each mixture, use extract.group()
#   to construct single-group extractions of the mixtures
mix1 <- extract.group(mod_mix, group = "MIXTURE_1")
mix2 <- extract.group(mod_mix, group = "MIXTURE_2")

# EAP estimates per mixture group, ignoring the original mixture structure.
#   Used to demonstrate the behaviour of how the individuals would have
#   been scored if they (deterministically) belonged to one class
data.frame(EAP_mix1=unname(fscores(mix1)),
           EAP_mix2=unname(fscores(mix2)),
           EAP=unname(fscores(mod_mix))) |> head()
#>      EAP_mix1   EAP_mix2        EAP
#> 1 -0.46823421 -0.3555457 -0.3609861
#> 2 -1.65665303 -1.8056233 -1.6950401
#> 3 -0.83303093 -0.8214973 -0.8217853
#> 4  0.46049686  0.8194937  0.8089769
#> 5  0.46049686  0.8194937  0.7935688
#> 6  0.07565017  0.3412955  0.2286019

############
# Mixture 2PL model
mod_mix2 <- multipleGroup(dat, 1, dentype = 'mixture-2', nruns=5)
anova(mod_mix, mod_mix2)
#>               AIC    SABIC       HQ      BIC    logLik     X2 df    p
#> mod_mix  34628.19 34720.06 34713.30 34856.65 -17271.09               
#> mod_mix2 34655.23 34828.29 34815.56 35085.60 -17246.62 48.957 38 0.11
coef(mod_mix2, simplify=TRUE)
#> $MIXTURE_1
#> $items
#>            a1      d g u
#> Item_1  0.517  0.673 0 1
#> Item_2  0.516  0.826 0 1
#> Item_3  0.533 -0.012 0 1
#> Item_4  0.386 -0.363 0 1
#> Item_5  0.452  0.552 0 1
#> Item_6  0.730 -1.903 0 1
#> Item_7  1.031  0.879 0 1
#> Item_8  0.744 -0.274 0 1
#> Item_9  0.900  0.038 0 1
#> Item_10 0.710 -0.777 0 1
#> Item_11 0.539 -0.189 0 1
#> Item_12 0.916  2.080 0 1
#> Item_13 0.729  0.691 0 1
#> Item_14 1.184  0.723 0 1
#> Item_15 0.812 -0.716 0 1
#> Item_16 0.582  0.843 0 1
#> Item_17 0.792 -0.952 0 1
#> Item_18 0.657 -0.289 0 1
#> Item_19 0.498  1.222 0 1
#> Item_20 0.468  0.304 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> $class_proportion
#>    pi
#>  0.35
#> 
#> 
#> $MIXTURE_2
#> $items
#>            a1      d g u
#> Item_1  1.118  0.710 0 1
#> Item_2  1.216  1.405 0 1
#> Item_3  1.244 -0.566 0 1
#> Item_4  1.684 -1.827 0 1
#> Item_5  1.669 -1.802 0 1
#> Item_6  1.474  2.152 0 1
#> Item_7  1.143 -0.511 0 1
#> Item_8  1.244  0.574 0 1
#> Item_9  1.272  0.506 0 1
#> Item_10 1.188 -0.183 0 1
#> Item_11 1.466  0.979 0 1
#> Item_12 1.403  2.208 0 1
#> Item_13 1.388  2.110 0 1
#> Item_14 1.062  1.599 0 1
#> Item_15 1.311  0.112 0 1
#> Item_16 1.191  0.409 0 1
#> Item_17 1.382 -0.388 0 1
#> Item_18 1.500 -1.799 0 1
#> Item_19 1.625  1.955 0 1
#> Item_20 1.353 -0.060 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 
#> $class_proportion
#>    pi
#>  0.65
#> 
#> 
itemfit(mod_mix2)
#>       item   S_X2 df.S_X2 RMSEA.S_X2 p.S_X2
#> 1   Item_1 21.271      12      0.023  0.047
#> 2   Item_2 44.607      13      0.040  0.000
#> 3   Item_3 14.259      12      0.011  0.284
#> 4   Item_4  6.551      11      0.000  0.834
#> 5   Item_5 20.641      12      0.022  0.056
#> 6   Item_6 27.961      13      0.028  0.009
#> 7   Item_7  8.628      12      0.000  0.734
#> 8   Item_8 16.262      12      0.015  0.180
#> 9   Item_9 21.305      12      0.023  0.046
#> 10 Item_10 11.461      12      0.000  0.490
#> 11 Item_11  7.091      12      0.000  0.852
#> 12 Item_12 13.834      12      0.010  0.311
#> 13 Item_13 10.444      12      0.000  0.577
#> 14 Item_14 15.706      13      0.012  0.265
#> 15 Item_15  9.439      12      0.000  0.665
#> 16 Item_16 10.128      12      0.000  0.605
#> 17 Item_17 17.474      12      0.017  0.133
#> 18 Item_18 11.838      11      0.007  0.376
#> 19 Item_19  8.069      11      0.000  0.707
#> 20 Item_20 14.542      12      0.012  0.267

# Compare to single group
mod <- mirt(dat)
anova(mod, mod_mix2)
#>               AIC    SABIC       HQ      BIC    logLik      X2 df p
#> mod      35276.70 35362.16 35355.88 35489.23 -17598.35             
#> mod_mix2 34655.23 34828.29 34815.56 35085.60 -17246.62 703.474 41 0

########################################
# Zero-inflated 2PL IRT model (Wall, Park, and Moustaki, 2015)

n <- 1000
nitems <- 20

a <- rep(2, nitems)
d <- rep(c(-2,-1,0,1,2), each=nitems/5)
zi_p <- 0.2 # Proportion of people in zero class

theta <- rnorm(n, 0, 1)
zeros <- matrix(0, n*zi_p, nitems)
nonzeros <- simdata(a, d, n*(1-zi_p), itemtype = '2PL',
                   Theta = as.matrix(theta[1:(n*(1-zi_p))]))
data <- rbind(nonzeros, zeros)

# define class with extreme theta but fixed item parameters
zi2PL <- "F = 1-20
          START [MIXTURE_1] = (GROUP, MEAN_1, -100), (GROUP, COV_11, .00001),
                              (1-20, a1, 1.0), (1-20, d, 0)
          FIXED [MIXTURE_1] = (GROUP, MEAN_1), (GROUP, COV_11),
                              (1-20, a1), (1-20, d)"

# define custom Theta integration grid that contains extreme theta + normal grid
technical <- list(customTheta = matrix(c(-100, seq(-6,6,length.out=61))))

# fit ZIM-IRT
zi2PL.fit <- multipleGroup(data, zi2PL, dentype = 'mixture-2', technical=technical)
coef(zi2PL.fit, simplify=TRUE)
#> $MIXTURE_1
#> $items
#>         a1 d g u
#> Item_1   1 0 0 1
#> Item_2   1 0 0 1
#> Item_3   1 0 0 1
#> Item_4   1 0 0 1
#> Item_5   1 0 0 1
#> Item_6   1 0 0 1
#> Item_7   1 0 0 1
#> Item_8   1 0 0 1
#> Item_9   1 0 0 1
#> Item_10  1 0 0 1
#> Item_11  1 0 0 1
#> Item_12  1 0 0 1
#> Item_13  1 0 0 1
#> Item_14  1 0 0 1
#> Item_15  1 0 0 1
#> Item_16  1 0 0 1
#> Item_17  1 0 0 1
#> Item_18  1 0 0 1
#> Item_19  1 0 0 1
#> Item_20  1 0 0 1
#> 
#> $means
#>    F 
#> -100 
#> 
#> $cov
#>   F
#> F 0
#> 
#> $class_proportion
#>   pi
#>  0.2
#> 
#> 
#> $MIXTURE_2
#> $items
#>            a1      d g u
#> Item_1  2.145 -1.865 0 1
#> Item_2  2.070 -2.187 0 1
#> Item_3  1.987 -2.169 0 1
#> Item_4  2.276 -2.282 0 1
#> Item_5  1.831 -0.875 0 1
#> Item_6  1.633 -0.769 0 1
#> Item_7  1.743 -0.861 0 1
#> Item_8  1.715 -0.960 0 1
#> Item_9  1.960  0.039 0 1
#> Item_10 1.860  0.061 0 1
#> Item_11 2.089  0.184 0 1
#> Item_12 2.092 -0.002 0 1
#> Item_13 1.964  1.070 0 1
#> Item_14 1.727  0.979 0 1
#> Item_15 1.842  1.162 0 1
#> Item_16 1.714  0.910 0 1
#> Item_17 2.326  2.527 0 1
#> Item_18 1.846  2.095 0 1
#> Item_19 1.855  1.979 0 1
#> Item_20 1.900  2.331 0 1
#> 
#> $means
#> F 
#> 0 
#> 
#> $cov
#>   F
#> F 1
#> 
#> $class_proportion
#>   pi
#>  0.8
#> 
#> 

# classification estimates
pi_hat <- fscores(zi2PL.fit, method = 'classify')
head(pi_hat)
#>            CLASS_1 CLASS_2
#> [1,] 4.549844e-164       1
#> [2,] 6.529859e-210       1
#> [3,] 2.336737e-223       1
#> [4,] 3.525975e-164       1
#> [5,] 1.362734e-115       1
#> [6,]  9.240858e-29       1
tail(pi_hat)
#>           CLASS_1    CLASS_2
#>  [995,] 0.9188292 0.08117085
#>  [996,] 0.9188292 0.08117085
#>  [997,] 0.9188292 0.08117085
#>  [998,] 0.9188292 0.08117085
#>  [999,] 0.9188292 0.08117085
#> [1000,] 0.9188292 0.08117085

# EAP estimates (not useful for zip class)
fs <- fscores(zi2PL.fit)
head(fs)
#>         Class_1
#> [1,]  0.0512440
#> [2,]  0.5477758
#> [3,]  0.7687951
#> [4,]  0.0973191
#> [5,] -0.3378649
#> [6,] -1.4372167
tail(fs)
#>           Class_1
#>  [995,] -92.05499
#>  [996,] -92.05499
#>  [997,] -92.05499
#>  [998,] -92.05499
#>  [999,] -92.05499
#> [1000,] -92.05499

########################################
# Zero-inflated graded response model (Magnus and Garnier-Villarreal, 2022)

n <- 1000
nitems <- 20

a <- matrix(rlnorm(20,.2,.3))

# for the graded model, ensure that there is enough space between the intercepts,
# otherwise closer categories will not be selected often (minimum distance of 0.3 here)
diffs <- t(apply(matrix(runif(20*4, .3, 1), 20), 1, cumsum))
diffs <- -(diffs - rowMeans(diffs))
d <- diffs + rnorm(20)

zi_p <- 0.2 # Proportion of people in zero/lowest category class

theta <- rnorm(n, 0, 1)
zeros <- matrix(0, n*zi_p, nitems)
nonzeros <- simdata(a, d, n*(1-zi_p), itemtype = 'graded',
                    Theta = as.matrix(theta[1:(n*(1-zi_p))]))
data <- rbind(nonzeros, zeros)

# intercepts will be labelled as d1 through d4
apply(data, 2, table)
#>   Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> 0    688    450    414    359    545    454    334    418    451     346
#> 1     36     44    149     60     48    105     44     67     87     103
#> 2     96     52    118     93     68     68     71     58     49      77
#> 3     32    110    130    156     98    124    107     41    101      85
#> 4    148    344    189    332    241    249    444    416    312     389
#>   Item_11 Item_12 Item_13 Item_14 Item_15 Item_16 Item_17 Item_18 Item_19
#> 0     497     327     437     433     456     416     346     550     715
#> 1     159      83     132      63      67      91      60     142      90
#> 2     124     126      94     112      74      60      62      67      77
#> 3      64     146     142     103      45     114     108      70      52
#> 4     156     318     195     289     358     319     424     171      66
#>   Item_20
#> 0     458
#> 1      85
#> 2     131
#> 3      80
#> 4     246

# ignoring zero inflation (bad idea)
modGRM <- mirt(data)
coef(modGRM, simplify=TRUE)
#> $items
#>            a1     d1     d2     d3     d4
#> Item_1  2.818 -1.826 -2.108 -2.958 -3.294
#> Item_2  2.839 -0.045 -0.390 -0.777 -1.579
#> Item_3  1.763  0.319 -0.583 -1.261 -2.134
#> Item_4  2.340  0.701  0.219 -0.426 -1.439
#> Item_5  2.560 -0.715 -1.049 -1.516 -2.237
#> Item_6  2.134  0.028 -0.655 -1.075 -1.896
#> Item_7  3.303  1.013  0.521 -0.141 -1.018
#> Item_8  2.933  0.194 -0.360 -0.809 -1.113
#> Item_9  2.761 -0.060 -0.714 -1.065 -1.795
#> Item_10 2.751  0.836 -0.077 -0.646 -1.238
#> Item_11 2.031 -0.244 -1.229 -2.072 -2.613
#> Item_12 2.590  1.020  0.247 -0.687 -1.675
#> Item_13 2.278  0.128 -0.763 -1.380 -2.433
#> Item_14 3.058  0.039 -0.489 -1.357 -2.148
#> Item_15 4.557 -0.411 -1.117 -1.869 -2.317
#> Item_16 2.128  0.296 -0.311 -0.687 -1.412
#> Item_17 3.586  0.889  0.234 -0.350 -1.277
#> Item_18 2.403 -0.667 -1.617 -2.110 -2.716
#> Item_19 1.822 -1.505 -2.149 -2.868 -3.596
#> Item_20 2.852 -0.144 -0.796 -1.750 -2.374
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 

# Define class with extreme theta but fixed item parameters
#   For GRM in zero-inflated class the intercept values are arbitrary
#   as the model forces the responses all into the first category (hence,
#   spacing arbitrarily set to 1)
ziGRM <- "F = 1-20
          START [MIXTURE_1] = (GROUP, MEAN_1, -100), (GROUP, COV_11, .00001),
                              (1-20, a1, 1.0),
                              (1-20, d1, 2), (1-20, d2, 1), (1-20, d3, 0), (1-20, d4, -1)
          FIXED [MIXTURE_1] = (GROUP, MEAN_1), (GROUP, COV_11),
                              (1-20, a1),
                              (1-20, d1), (1-20, d2), (1-20, d3), (1-20, d4)"

# define custom Theta integration grid that contains extreme theta + normal grid
technical <- list(customTheta = matrix(c(-100, seq(-6,6,length.out=61))))

# fit zero-inflated GRM
ziGRM.fit <- multipleGroup(data, ziGRM, dentype = 'mixture-2', technical=technical)
coef(ziGRM.fit, simplify=TRUE)
#> $MIXTURE_1
#> $items
#>         a1 d1 d2 d3 d4
#> Item_1   1  2  1  0 -1
#> Item_2   1  2  1  0 -1
#> Item_3   1  2  1  0 -1
#> Item_4   1  2  1  0 -1
#> Item_5   1  2  1  0 -1
#> Item_6   1  2  1  0 -1
#> Item_7   1  2  1  0 -1
#> Item_8   1  2  1  0 -1
#> Item_9   1  2  1  0 -1
#> Item_10  1  2  1  0 -1
#> Item_11  1  2  1  0 -1
#> Item_12  1  2  1  0 -1
#> Item_13  1  2  1  0 -1
#> Item_14  1  2  1  0 -1
#> Item_15  1  2  1  0 -1
#> Item_16  1  2  1  0 -1
#> Item_17  1  2  1  0 -1
#> Item_18  1  2  1  0 -1
#> Item_19  1  2  1  0 -1
#> Item_20  1  2  1  0 -1
#> 
#> $means
#>    F 
#> -100 
#> 
#> $cov
#>   F
#> F 0
#> 
#> $class_proportion
#>     pi
#>  0.199
#> 
#> 
#> $MIXTURE_2
#> $items
#>            a1     d1     d2     d3     d4
#> Item_1  1.624 -0.668 -0.951 -1.802 -2.138
#> Item_2  1.525  1.124  0.780  0.395 -0.398
#> Item_3  0.726  1.114  0.205 -0.456 -1.297
#> Item_4  1.088  1.684  1.190  0.549 -0.432
#> Item_5  1.403  0.353  0.020 -0.445 -1.160
#> Item_6  1.052  0.935  0.255 -0.160 -0.961
#> Item_7  1.722  2.342  1.845  1.187  0.327
#> Item_8  1.553  1.391  0.840  0.396  0.096
#> Item_9  1.479  1.081  0.429  0.081 -0.639
#> Item_10 1.371  1.965  1.046  0.486 -0.092
#> Item_11 1.018  0.628 -0.347 -1.171 -1.698
#> Item_12 1.244  2.106  1.311  0.392 -0.563
#> Item_13 1.162  1.093  0.203 -0.403 -1.431
#> Item_14 1.670  1.301  0.773 -0.091 -0.873
#> Item_15 2.599  1.469  0.755 -0.003 -0.455
#> Item_16 0.991  1.194  0.585  0.214 -0.488
#> Item_17 1.916  2.332  1.678  1.100  0.188
#> Item_18 1.302  0.341 -0.602 -1.090 -1.686
#> Item_19 0.986 -0.718 -1.358 -2.069 -2.788
#> Item_20 1.550  1.036  0.386 -0.559 -1.175
#> 
#> $means
#> F 
#> 0 
#> 
#> $cov
#>   F
#> F 1
#> 
#> $class_proportion
#>     pi
#>  0.801
#> 
#> 

# classification estimates
pi_hat <- fscores(ziGRM.fit, method = 'classify')
head(pi_hat)
#>            CLASS_1 CLASS_2
#> [1,] 4.558633e-300       1
#> [2,] 3.846940e-248       1
#> [3,] 1.372135e-214       1
#> [4,] 1.220535e-224       1
#> [5,] 2.774206e-126       1
#> [6,] 8.219653e-210       1
tail(pi_hat)
#>           CLASS_1     CLASS_2
#>  [995,] 0.9932982 0.006701771
#>  [996,] 0.9932982 0.006701771
#>  [997,] 0.9932982 0.006701771
#>  [998,] 0.9932982 0.006701771
#>  [999,] 0.9932982 0.006701771
#> [1000,] 0.9932982 0.006701771

# EAP estimates (not useful for zip class)
fs <- fscores(ziGRM.fit)
head(fs)
#>          Class_1
#> [1,]  0.03415732
#> [2,]  0.77178287
#> [3,] -0.08031843
#> [4,] -0.29479164
#> [5,] -1.27315493
#> [6,] -0.45721982
tail(fs)
#>           Class_1
#>  [995,] -99.34781
#>  [996,] -99.34781
#>  [997,] -99.34781
#>  [998,] -99.34781
#>  [999,] -99.34781
#> [1000,] -99.34781

# }
```
