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
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

# }
```
