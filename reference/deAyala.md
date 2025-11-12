# Description of deAyala data

Mathematics data from de Ayala (2009; pg. 14); 5 item dataset in table
format.

## References

de Ayala, R. J. (2009). *The theory and practice of item response
theory*. Guilford Press.

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
dat <- expand.table(deAyala)
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
#>      N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  19601            2.912          1.434 0.233 0.074 0.608     0.898
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item.1 19601 2 0.887 0.316   0.447         0.246       0.605
#> Item.2 19601 2 0.644 0.479   0.688         0.439       0.510
#> Item.3 19601 2 0.566 0.496   0.680         0.416       0.523
#> Item.4 19601 2 0.427 0.495   0.673         0.405       0.529
#> Item.5 19601 2 0.387 0.487   0.602         0.312       0.581
#> 
#> $proportions
#>            0     1
#> Item.1 0.113 0.887
#> Item.2 0.356 0.644
#> Item.3 0.434 0.566
#> Item.4 0.573 0.427
#> Item.5 0.613 0.387
#> 

# }
```
