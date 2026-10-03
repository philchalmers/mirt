# Simulated datasets for PIRT-DIF

Three associated datasets for PIRT-DIF, stored as a list (more
information to come).

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# \donttest{
data(pirt_DIF)

# dataset 1
dat1 <- pirt_DIF$dat1
group <- dat1$group
dat <- dat1[,-1]
itemstats(dat, group=group)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

# }
```
