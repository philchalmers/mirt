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

# fit 2PL model for each item
(mod <- mirt(dat))
#> 
#> Call:
#> mirt(data = dat)
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

(mod <- mirt(dat, 1))
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
