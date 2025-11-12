# Remap item categories to have integer distances of 1

The mirt package's estimation setup requires that all item responses
have spaces equal to 1 (e.g., a Likert scale scored from 1 through 5).
In the event that categories are missing the categories must be
re-coded. This function is automatically called by the package
estimation functions (e.g.,
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)),
however for convince this function has been extracted for users to
better understand the remapping consequences.

## Usage

``` r
remap.distance(data, message = TRUE)
```

## Arguments

- data:

  the response data to remap as a data.frame or matrix

- message:

  logical; print message information pertaining to which items were
  remapped?

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# category 2 for item 1 missing
dat <- Science
dat[,1] <- ifelse(Science[,1] == 2, 1, Science[,1])
apply(dat, 2, table)
#> $Comfort
#> 
#>   1   3   4 
#>  37 266  89 
#> 
#> $Work
#> 
#>   1   2   3   4 
#>  33  98 206  55 
#> 
#> $Future
#> 
#>   1   2   3   4 
#>  14  72 210  96 
#> 
#> $Benefit
#> 
#>   1   2   3   4 
#>  21 100 193  78 
#> 

# mirt() automatically remaps categories
mod <- mirt(dat, 1)
#> "Comfort" re-mapped to ensure all categories have a distance of 1
coef(mod, simplify=TRUE)
#> $items
#>            a1    d1     d2     d3
#> Comfort 0.995 2.626 -1.447     NA
#> Work    1.231 2.928  0.903 -2.271
#> Future  2.338 5.293  2.245 -1.989
#> Benefit 1.079 3.333  0.988 -1.681
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 

# this is the transformed data used by mirt()
remap_dat <- remap.distance(dat)
#> "Comfort" re-mapped to ensure all categories have a distance of 1
apply(remap_dat, 2, table)
#> $Comfort
#> 
#>   1   2   3 
#>  37 266  89 
#> 
#> $Work
#> 
#>   1   2   3   4 
#>  33  98 206  55 
#> 
#> $Future
#> 
#>   1   2   3   4 
#>  14  72 210  96 
#> 
#> $Benefit
#> 
#>   1   2   3   4 
#>  21 100 193  78 
#> 

```
