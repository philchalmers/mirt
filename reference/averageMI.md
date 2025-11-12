# Collapse values from multiple imputation draws

This function computes updated parameter and standard error estimates
using multiple imputation methodology. Given a set of parameter
estimates and their associated standard errors the function returns the
weighted average of the overall between and within variability due to
the multiple imputations according to Rubin's (1987) methodology.

## Usage

``` r
averageMI(par, SEpar, as.data.frame = TRUE)
```

## Arguments

- par:

  a list containing parameter estimates which were computed the imputed
  datasets

- SEpar:

  a list containing standard errors associated with `par`

- as.data.frame:

  logical; return a data.frame instead of a list? Default is TRUE

## Value

returns a list or data.frame containing the updated averaged parameter
estimates, standard errors, and t-values with the associated degrees of
freedom and two tailed p-values

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Rubin, D.B. (1987) Multiple Imputation for Nonresponse in Surveys. Wiley
& Sons, New York.

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

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

mod1 <- mirt(dat, 1, 'Rasch', covdata=covdata, formula = ~ X1 + X2)
coef(mod1, simplify=TRUE)
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

# draw plausible values for secondary analyses
pv <- fscores(mod1, plausible.draws = 10)
pvmods <- lapply(pv, function(x, covdata) lm(x ~ covdata$X1 + covdata$X2),
                 covdata=covdata)

# compute Rubin's multiple imputation average
so <- lapply(pvmods, summary)
par <- lapply(so, function(x) x$coefficients[, 'Estimate'])
SEpar <- lapply(so, function(x) x$coefficients[, 'Std. Error'])
averageMI(par, SEpar)
#>                par SEpar       t      df     p
#> (Intercept)  0.003 0.016   0.209 198.552 0.209
#> covdata$X1   0.528 0.018  28.554  63.384     0
#> covdata$X2  -1.037 0.021 -49.311  35.705     0

# }
```
