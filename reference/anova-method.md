# Compare nested models with likelihood-based statistics

Compare nested models using likelihood ratio test (X2), Akaike
Information Criterion (AIC), Bayesian Information Criterion (BIC),
Sample-Size Adjusted BIC (SABIC), and Hannan-Quinn (HQ) Criterion. When
given a sequence of objects, `anova` tests the models against one
another in the order specified. Note that the `object` inputs should be
ordered in terms of most constrained model to least constrained.

## Usage

``` r
# S4 method for class 'SingleGroupClass'
anova(
  object,
  object2,
  ...,
  bounded = FALSE,
  mix = 0.5,
  frame = 1,
  verbose = FALSE
)
```

## Arguments

- object:

  an object of class `SingleGroupClass`, `MultipleGroupClass`, or
  `MixedClass`, reflecting the most constrained model fitted

- object2:

  a second model estimated from any of the mirt package estimation
  methods

- ...:

  additional less constrained model objects to be compared sequentially
  to the previous model

- bounded:

  logical; are the two models comparing a bounded parameter (e.g.,
  comparing a single 2PL and 3PL model with 1 df)? If `TRUE` then a
  50:50 mix of chi-squared distributions is used to obtain the p-value

- mix:

  proportion of chi-squared mixtures. Default is 0.5

- frame:

  (internal parameter not for standard use)

- verbose:

  (deprecated argument)

## Value

a `data.frame`/`mirt_df` object

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Examples

``` r
# \donttest{
x <- mirt(Science, 1)
x2 <- mirt(Science, 2)
anova(x, x2)
#>         AIC    SABIC       HQ      BIC    logLik     X2 df     p
#> x  3249.739 3262.512 3274.922 3313.279 -1608.870                
#> x2 3241.938 3257.106 3271.843 3317.392 -1601.969 13.801  3 0.003

# compare three models sequentially (X2 not always meaningful)
x3 <- mirt(Science, 1, 'gpcm')
x4 <- mirt(Science, 1, 'nominal')
anova(x, x2, x3, x4)
#>         AIC    SABIC       HQ      BIC    logLik      X2 df     p
#> x  3249.739 3262.512 3274.922 3313.279 -1608.870                 
#> x2 3241.938 3257.106 3271.843 3317.392 -1601.969  13.801  3 0.003
#> x3 3257.366 3270.139 3282.549 3320.906 -1612.683 -21.428 -3   NaN
#> x4 3264.910 3284.069 3302.684 3360.220 -1608.455   8.456  8  0.39

# in isolation
anova(x)
#>        AIC    SABIC       HQ      BIC   logLik
#> x 3249.739 3262.512 3274.922 3313.279 -1608.87

# with priors on first model
model <- "Theta = 1-4
          PRIOR = (1-4, a1, lnorm, 0, 10)"
xp <- mirt(Science, model)
anova(xp, x2)
#>         AIC    SABIC       HQ      BIC    logLik   logPost df
#> xp 3249.829 3262.602 3275.012 3313.369 -1608.914 -1622.881 NA
#> x2 3241.938 3257.106 3271.843 3317.392 -1601.969 -1601.969  3
anova(xp)
#>         AIC    SABIC       HQ      BIC    logLik   logPost
#> xp 3249.829 3262.602 3275.012 3313.369 -1608.914 -1622.881

# bounded parameter
dat <- expand.table(LSAT7)
mod <- mirt(dat, 1)
mod2 <- mirt(dat, 1, itemtype = c(rep('2PL', 4), '3PL'))
anova(mod, mod2) #unbounded test
#>           AIC    SABIC       HQ      BIC    logLik    X2 df    p
#> mod  5337.610 5354.927 5356.263 5386.688 -2658.805              
#> mod2 5339.587 5358.636 5360.106 5393.573 -2658.794 0.023  1 0.88
anova(mod, mod2, bounded = TRUE) #bounded
#>           AIC    SABIC       HQ      BIC    logLik    X2 df    p
#> mod  5337.610 5354.927 5356.263 5386.688 -2658.805              
#> mod2 5339.587 5358.636 5360.106 5393.573 -2658.794 0.023  1 0.44

# priors
model <- 'F = 1-5
          PRIOR = (5, g, norm, -1, 1)'
mod1b <- mirt(dat, model, itemtype = c(rep('2PL', 4), '3PL'))
anova(mod1b)
#>            AIC   SABIC       HQ      BIC    logLik   logPost
#> mod1b 5339.571 5358.62 5360.089 5393.557 -2658.786 -2659.705

model2 <- 'F = 1-5
          PRIOR = (1-5, g, norm, -1, 1)'
mod2b <- mirt(dat, model2, itemtype = '3PL')
anova(mod1b, mod2b)
#>            AIC    SABIC       HQ      BIC    logLik   logPost df
#> mod1b 5339.571 5358.620 5360.089 5393.557 -2658.786 -2659.705 NA
#> mod2b 5348.306 5374.282 5376.286 5421.923 -2659.153 -2664.008  4

# }
```
