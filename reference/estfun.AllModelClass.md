# Extract Empirical Estimating Functions

A function for extracting the empirical estimating functions of a fitted
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
[`multipleGroup`](https://philchalmers.github.io/mirt/reference/multipleGroup.md),
[`bfactor`](https://philchalmers.github.io/mirt/reference/bfactor.md),
or [`mdirt`](https://philchalmers.github.io/mirt/reference/mdirt.md)
model. This is the derivative of the log-likelihood with respect to the
parameter vector, evaluated at the observed (case-wise) data. In other
words, this function returns the case-wise scores, evaluated at the
fitted model parameters. Currently, models fitted via the `EM` or `BL`
method are supported. For the computations, the internal `Theta` grid of
the model is being used which was already used during the estimation of
the model itself along with its matching normalized density.

## Usage

``` r
estfun.AllModelClass(
  x,
  weights = extract.mirt(x, "survey.weights"),
  centering = FALSE
)
```

## Arguments

- x:

  a fitted model object of class `SingleGroupClass`,
  `MultipleGroupClass`, or `DiscreteClass`

- weights:

  by default, the `survey.weights` which were (optionally) specified
  when fitting the model are included to calculate the scores. If
  specified by the user, this should be a numeric vector of length equal
  to the total sample size. Note that if not all cases were weighted
  equally when fitting the model, the weights must be corrected by
  taking their square root if the scores are being used to compute the
  outer product of gradients (OPG) estimate of the variance-covariance
  matrix (see examples below).

- centering:

  a boolean variable that allows the centering of the case-wise scores
  (i.e., setting their expected values to 0). If the case-wise scores
  were obtained from maximum likelihood estimates, this setting does not
  affect the result.

## Value

An n x k matrix corresponding to n observations and k parameters

## See also

[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
[`multipleGroup`](https://philchalmers.github.io/mirt/reference/multipleGroup.md),
[`bfactor`](https://philchalmers.github.io/mirt/reference/bfactor.md),
[`mdirt`](https://philchalmers.github.io/mirt/reference/mdirt.md)

## Author

Lennart Schneider <lennart.sch@web.de> and Phil Chalmers; centering
argument contributed by Rudolf Debelak
(<rudolf.debelak@psychologie.uzh.ch>)

## Examples

``` r
# \donttest{
# fit a 2PL on the LSAT7 data and get the scores
mod1 <- mirt(expand.table(LSAT7), 1, SE = TRUE, SE.type = "crossprod")
sc1 <- estfun.AllModelClass(mod1)
# get the gradient
colSums(sc1)
#>          a1.1           d.2          a1.5           d.6          a1.9 
#> -2.477600e-03 -1.749068e-03  8.292831e-05 -1.521604e-03  1.198404e-02 
#>          d.10         a1.13          d.14         a1.17          d.18 
#>  6.498262e-03 -2.188204e-03  4.417930e-03 -8.004882e-04 -5.566577e-04 
# calculate the OPG estimate of the variance-covariance matrix "by hand"
vc1 <- vcov(mod1)
all.equal(crossprod(sc1), chol2inv(chol(vc1)), check.attributes = FALSE)
#> [1] TRUE

# Discrete group
modd <- mdirt(expand.table(LSAT7), 2, SE = TRUE, SE.type = "crossprod")
sc1 <- estfun.AllModelClass(modd)
# get the gradient
colSums(sc1)
#>        a1.1        a2.2        a1.3        a2.4        a1.5        a2.6 
#> 0.010556101 0.012533394 0.015994245 0.016832573 0.021043233 0.020842514 
#>        a1.7        a2.8        a1.9       a2.10       c1.11 
#> 0.013694253 0.011470847 0.007827056 0.009231532 0.059032788 
# calculate the OPG estimate of the variance-covariance matrix "by hand"
vc1 <- vcov(modd)
all.equal(crossprod(sc1), chol2inv(chol(vc1)), check.attributes = FALSE)
#> [1] TRUE

# fit a multiple group 2PL and do the same as above
group <- rep(c("G1", "G2"), 500)
mod2 <- multipleGroup(expand.table(LSAT7), 1, group, SE = TRUE,
  SE.type = "crossprod")
sc2 <- estfun.AllModelClass(mod2)
colSums(sc2)
#>          a1.1           d.2          a1.5           d.6          a1.9 
#> -3.118728e-04 -2.982410e-03 -2.735485e-04 -6.118069e-04  7.440279e-03 
#>          d.10         a1.13          d.14         a1.17          d.18 
#>  7.274125e-03 -9.423217e-04 -3.379603e-04  6.496281e-04 -1.868471e-03 
#>         a1.23          d.24         a1.27          d.28         a1.31 
#> -1.197098e-04 -4.435649e-04 -2.378475e-05  1.458023e-04  6.403262e-04 
#>          d.32         a1.35          d.36         a1.39          d.40 
#>  1.212925e-03 -1.443583e-04  4.381246e-05  4.279434e-04 -9.488836e-04 
vc2 <- vcov(mod2)
all.equal(crossprod(sc2), chol2inv(chol(vc2)), check.attributes = FALSE)
#> [1] TRUE

# fit a bifactor model with 2 specific factors and do the same as above
mod3 <- bfactor(expand.table(LSAT7), c(2, 2, 1, 1, 2), SE = TRUE,
  SE.type = "crossprod")
#> 
#> 
#> Calculating information matrix...
sc3 <- estfun.AllModelClass(mod3)
colSums(sc3)
#>          a1.1          a3.3           d.4          a1.7          a3.9 
#>  0.0018172801  0.0007522911  0.0008693077  0.0026167974 -0.0065562618 
#>          d.10         a1.13         a2.14          d.16         a1.19 
#>  0.0010184719  0.0002582749 -0.0123565804  0.0003107030 -0.0001820340 
#>         a2.20          d.22         a1.25         a3.27          d.28 
#> -0.0048417996  0.0020931538 -0.0005721500 -0.0060778033 -0.0028105998 
vc3 <- vcov(mod3)
all.equal(crossprod(sc3), chol2inv(chol(vc3)), check.attributes = FALSE)
#> [1] TRUE

# fit a 2PL not weighting all cases equally
survey.weights <- c(rep(2, sum(LSAT7$freq) / 2), rep(1, sum(LSAT7$freq) / 2))
survey.weights <- survey.weights / sum(survey.weights) * sum(LSAT7$freq)
mod4 <- mirt(expand.table(LSAT7), 1, SE = TRUE, SE.type = "crossprod",
  survey.weights = survey.weights)
sc4 <- estfun.AllModelClass(mod4,
  weights = extract.mirt(mod4, "survey.weights"))
# get the gradient
colSums(sc4)
#>          a1.1           d.2          a1.5           d.6          a1.9 
#> -0.0067330450 -0.0001045753  0.0019680853  0.0007675188  0.0088929479 
#>          d.10         a1.13          d.14         a1.17          d.18 
#>  0.0012985459  0.0016433039 -0.0008952128 -0.0071939497 -0.0031136208 
# to calculate the OPG estimate of the variance-covariance matrix "by hand",
# the weights must be adjusted by taking their square root
sc4_crp <- estfun.AllModelClass(mod4,
  weights = sqrt(extract.mirt(mod4, "survey.weights")))
vc4 <- vcov(mod4)
all.equal(crossprod(sc4_crp), chol2inv(chol(vc4)), check.attributes = FALSE)
#> [1] TRUE

# }
```
