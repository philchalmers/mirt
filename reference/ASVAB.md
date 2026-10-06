# Description of ASVAB data

Table of counts extracted from Mislvey (1985). Data the 16 possible
response patterns observed for four items from the arithmetic reasoning
test of the Armed Services Vocational Aptitude Battery (ASVAB), Form 8A,
from samples of white males and females and black males and females.

## References

Mislevy, R. J. (1985). Estimation of latent group effects. *Journal of
the American Statistical Association, 80*, 993-997.

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

data(ASVAB)
datWM <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, White_Male)))
datWF <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, White_Female)))
datBM <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, Black_Male)))
datBF <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, Black_Female)))

dat <- rbind(datWM, datWF, datBM, datBF)
sex <- rep(c("Male", "Female", "Male", "Female"),
  times=c(nrow(datWM), nrow(datWF), nrow(datBM), nrow(datBF))) |> factor()
color <- rep(c("White", "Black"),
  times=c(nrow(datWM) + nrow(datWF), nrow(datBM) + nrow(datBF))) |> factor()
group <- sex:color

itemstats(dat, group=group)
#> $`Female:Black`
#> $`Female:Black`$overall
#>    N mean.total sd.total std.alpha alpha SEM.alpha mean.r  sd.r
#>  145      1.462    1.014     0.163 0.176     0.921  0.046 0.087
#> 
#> $`Female:Black`$itemstats
#>          N K  mean    sd cor_if_rm   smc alpha_if_rm
#> Item.1 145 2 0.503 0.502     0.213 0.047      -0.078
#> Item.2 145 2 0.421 0.495     0.074 0.026       0.152
#> Item.3 145 2 0.283 0.452     0.064 0.022       0.164
#> Item.4 145 2 0.255 0.437    -0.011 0.013       0.256
#> 
#> $`Female:Black`$proportions
#>            0     1
#> Item.1 0.497 0.503
#> Item.2 0.579 0.421
#> Item.3 0.717 0.283
#> Item.4 0.745 0.255
#> 
#> 
#> $`Female:White`
#> $`Female:White`$overall
#>    N mean.total sd.total std.alpha alpha SEM.alpha mean.r  sd.r
#>  228      2.118    1.255     0.512 0.512     0.877  0.208 0.037
#> 
#> $`Female:White`$itemstats
#>          N K  mean    sd cor_if_rm   smc alpha_if_rm
#> Item.1 228 2 0.618 0.487     0.277 0.083       0.464
#> Item.2 228 2 0.605 0.490     0.312 0.100       0.432
#> Item.3 228 2 0.487 0.501     0.284 0.085       0.458
#> Item.4 228 2 0.408 0.493     0.339 0.117       0.408
#> 
#> $`Female:White`$proportions
#>            0     1
#> Item.1 0.382 0.618
#> Item.2 0.395 0.605
#> Item.3 0.513 0.487
#> Item.4 0.592 0.408
#> 
#> 
#> $`Male:Black`
#> $`Male:Black`$overall
#>    N mean.total sd.total std.alpha alpha SEM.alpha mean.r  sd.r
#>  140      1.443    1.027     0.176  0.18      0.93  0.051 0.102
#> 
#> $`Male:Black`$itemstats
#>          N K  mean    sd cor_if_rm   smc alpha_if_rm
#> Item.1 140 2 0.443 0.499     0.158 0.035       0.029
#> Item.2 140 2 0.400 0.492     0.133 0.022       0.071
#> Item.3 140 2 0.329 0.471    -0.037 0.017       0.304
#> Item.4 140 2 0.271 0.446     0.100 0.055       0.123
#> 
#> $`Male:Black`$proportions
#>            0     1
#> Item.1 0.557 0.443
#> Item.2 0.600 0.400
#> Item.3 0.671 0.329
#> Item.4 0.729 0.271
#> 
#> 
#> $`Male:White`
#> $`Male:White`$overall
#>    N mean.total sd.total std.alpha alpha SEM.alpha mean.r  sd.r
#>  263      2.475    1.361     0.674 0.673     0.779   0.34 0.075
#> 
#> $`Male:White`$itemstats
#>          N K  mean    sd cor_if_rm   smc alpha_if_rm
#> Item.1 263 2 0.741 0.439     0.475 0.230       0.596
#> Item.2 263 2 0.635 0.482     0.361 0.136       0.667
#> Item.3 263 2 0.593 0.492     0.481 0.248       0.588
#> Item.4 263 2 0.506 0.501     0.507 0.282       0.569
#> 
#> $`Male:White`$proportions
#>            0     1
#> Item.1 0.259 0.741
#> Item.2 0.365 0.635
#> Item.3 0.407 0.593
#> Item.4 0.494 0.506
#> 
#> 
```
