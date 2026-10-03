# Description of LSAT7 data

Data from Bock & Lieberman (1970); contains 5 dichotomously scored items
obtained from the Law School Admissions Test, section 7.

Data from

## References

Bock, R. D., & Lieberman, M. (1970). Fitting a response model for *n*
dichotomously scored items. *Psychometrika, 35*(2), 179-197.

Bock, R. D., & Lieberman, M. (1970). Fitting a response model for *n*
dichotomously scored items. *Psychometrika, 35*(2), 179-197.

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# \donttest{
dat <- expand.table(LSAT7)
head(dat)
#>   Item.1 Item.2 Item.3 Item.4 Item.5
#> 1      0      0      0      0      0
#> 2      0      0      0      0      0
#> 3      0      0      0      0      0
#> 4      0      0      0      0      0
#> 5      0      0      0      0      0
#> 6      0      0      0      0      0
itemstats(dat)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

# fit 2PL model for each item
(mod <- mirt(dat))
#> 
#> Call:
#> mirt(data = dat)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within 1e-04 tolerance after 28 EM iterations.
#> mirt version: 1.47.5 
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
coef(mod)
#> $Item.1
#>        a1     d g u
#> par 0.988 1.856 0 1
#> 
#> $Item.2
#>        a1     d g u
#> par 1.081 0.808 0 1
#> 
#> $Item.3
#>        a1     d g u
#> par 1.706 1.804 0 1
#> 
#> $Item.4
#>        a1     d g u
#> par 0.765 0.486 0 1
#> 
#> $Item.5
#>        a1     d g u
#> par 0.736 1.855 0 1
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> 

# monotonic splines models (see Winsberg, Thissen, and Wainer, 1984)
mod_monospline <- mirt(dat, itemtype = 'monospline')
anova(mod, mod_monospline)
#>                    AIC    SABIC       HQ      BIC    logLik   X2 df     p
#> mod            5337.61 5354.927 5356.263 5386.688 -2658.805              
#> mod_monospline 5355.36 5389.994 5392.666 5453.515 -2657.680 2.25 10 0.994
plot(mod_monospline)


# compare item 1 trace-lines
i1 <- extract.item(mod, 1)
i1mono <- extract.item(mod_monospline, 1)
theta <- matrix(seq(-6, 6, length.out=100))
twoPL <- probtrace(i1, theta)[,2]
monospline <- probtrace(i1mono, theta)[,2]

plot(twoPL ~ theta, type = 'l')
lines(monospline ~ theta, col='red')


# }

# \donttest{
dat <- expand.table(LSAT7)
head(dat)
#>   Item.1 Item.2 Item.3 Item.4 Item.5
#> 1      0      0      0      0      0
#> 2      0      0      0      0      0
#> 3      0      0      0      0      0
#> 4      0      0      0      0      0
#> 5      0      0      0      0      0
#> 6      0      0      0      0      0
itemstats(dat)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

(mod <- mirt(dat, 1))
#> 
#> Call:
#> mirt(data = dat, model = 1)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within 1e-04 tolerance after 28 EM iterations.
#> mirt version: 1.47.5 
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
coef(mod)
#> $Item.1
#>        a1     d g u
#> par 0.988 1.856 0 1
#> 
#> $Item.2
#>        a1     d g u
#> par 1.081 0.808 0 1
#> 
#> $Item.3
#>        a1     d g u
#> par 1.706 1.804 0 1
#> 
#> $Item.4
#>        a1     d g u
#> par 0.765 0.486 0 1
#> 
#> $Item.5
#>        a1     d g u
#> par 0.736 1.855 0 1
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> 
# }
```
