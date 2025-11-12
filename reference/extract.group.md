# Extract a group from a multiple group mirt object

Extract a single group from an object defined by
[`multipleGroup`](https://philchalmers.github.io/mirt/reference/multipleGroup.md),
or as a mixture model from
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md).

## Usage

``` r
extract.group(x, group)
```

## Arguments

- x:

  model object of class `'MultipleGroupClass'` or `'MixtureClass'`

- group:

  the name of the group to extract

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`extract.item`](https://philchalmers.github.io/mirt/reference/extract.item.md),
[`extract.mirt`](https://philchalmers.github.io/mirt/reference/extract.mirt.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
itemtype <- rep('2PL', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))
models <- 'F1 = 1-15'

mod_configural <- multipleGroup(dat, models, group = group)
group.1 <- extract.group(mod_configural, 'D1') #extract first group
summary(group.1)
#>            F1    h2
#> Item_1  0.532 0.283
#> Item_2  0.582 0.338
#> Item_3  0.487 0.237
#> Item_4  0.466 0.218
#> Item_5  0.542 0.294
#> Item_6  0.315 0.099
#> Item_7  0.599 0.359
#> Item_8  0.477 0.228
#> Item_9  0.464 0.215
#> Item_10 0.391 0.153
#> Item_11 0.438 0.192
#> Item_12 0.655 0.429
#> Item_13 0.604 0.365
#> Item_14 0.519 0.270
#> Item_15 0.453 0.205
#> 
#> SS loadings:  3.885 
#> Proportion Var:  0.259 
#> 
#> Factor correlations: 
#> 
#>    F1
#> F1  1
plot(group.1)

# }
```
