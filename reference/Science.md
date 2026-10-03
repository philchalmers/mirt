# Description of Science data

A 4-item data set borrowed from `ltm` package in R, first example of the
`grm()` function. See more complete documentation therein, as well as
Karlheinz and Melich (1992).

## References

Karlheinz, R. and Melich, A. (1992). Euro-Barometer 38.1: *Consumer
Protection and Perceptions of Science and Technology*. INRA (Europe),
Brussels. \[computer file\]

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# \donttest{
itemstats(Science)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found

mod <- mirt(Science, 1)
plot(mod, type = 'trace')

# }
```
