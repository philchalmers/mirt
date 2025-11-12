# Imputing plausible data for missing values

Given an estimated model from any of mirt's model fitting functions and
an estimate of the latent trait, impute plausible missing data values.
Returns the original data in a `data.frame` without any NA values. If a
list of `Theta` values is supplied then a list of complete datasets is
returned instead.

## Usage

``` r
imputeMissing(x, Theta, warn = TRUE, ...)
```

## Arguments

- x:

  an estimated model x from the mirt package

- Theta:

  a matrix containing the estimates of the latent trait scores (e.g.,
  via
  [`fscores`](https://philchalmers.github.io/mirt/reference/fscores.md))

- warn:

  logical; print warning messages?

- ...:

  additional arguments to pass

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
dat <- expand.table(LSAT7)
(original <- mirt(dat, 1))
#> 
#> Call:
#> mirt(data = dat, model = 1)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within 1e-04 tolerance after 28 EM iterations.
#> mirt version: 1.45.5 
#> M-step optimizer: BFGS 
#> EM acceleration: Ramsay 
#> Number of rectangular quadrature: 61
#> Latent density type: Gaussian 
#> 
#> Log-likelihood = -2658.805
#> Estimated parameters: 10 
#> AIC = 5337.61
#> BIC = 5386.688; SABIC = 5354.927
#> G2 (21) = 31.7, p = 0.0628
#> RMSEA = 0.023, CFI = NaN, TLI = NaN
NAperson <- sample(1:nrow(dat), 20, replace = TRUE)
NAitem <- sample(1:ncol(dat), 20, replace = TRUE)
for(i in 1:20)
    dat[NAperson[i], NAitem[i]] <- NA
(mod <- mirt(dat, 1))
#> 
#> Call:
#> mirt(data = dat, model = 1)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within 1e-04 tolerance after 25 EM iterations.
#> mirt version: 1.45.5 
#> M-step optimizer: BFGS 
#> EM acceleration: Ramsay 
#> Number of rectangular quadrature: 61
#> Latent density type: Gaussian 
#> 
#> Log-likelihood = -2652.1
#> Estimated parameters: 10 
#> AIC = 5324.2
#> BIC = 5373.277; SABIC = 5341.517
#> 
scores <- fscores(mod, method = 'MAP')

# re-estimate imputed dataset (good to do this multiple times and average over)
fulldata <- imputeMissing(mod, scores)
(fullmod <- mirt(fulldata, 1))
#> 
#> Call:
#> mirt(data = fulldata, model = 1)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within 1e-04 tolerance after 24 EM iterations.
#> mirt version: 1.45.5 
#> M-step optimizer: BFGS 
#> EM acceleration: Ramsay 
#> Number of rectangular quadrature: 61
#> Latent density type: Gaussian 
#> 
#> Log-likelihood = -2660.331
#> Estimated parameters: 10 
#> AIC = 5340.663
#> BIC = 5389.741; SABIC = 5357.98
#> G2 (21) = 30.61, p = 0.0804
#> RMSEA = 0.021, CFI = NaN, TLI = NaN

# with multipleGroup
set.seed(1)
group <- sample(c('group1', 'group2'), 1000, TRUE)
mod2 <- multipleGroup(dat, 1, group, TOL=1e-2)
fs <- fscores(mod2)
fulldata2 <- imputeMissing(mod2, fs)

# }
```
