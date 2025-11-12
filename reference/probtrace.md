# Function to calculate probability trace lines

Given an internal mirt object extracted from an estimated model, or the
single-group estimated model itself, compute the probability trace lines
for all categories.

## Usage

``` r
probtrace(x, Theta)
```

## Arguments

- x:

  either an extracted internal mirt object containing item information
  (see
  [`extract.item`](https://philchalmers.github.io/mirt/reference/extract.item.md))
  or a model of class `SingleGroupClass` typically returned by the
  function
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) or
  [`bfactor`](https://philchalmers.github.io/mirt/reference/bfactor.md)

- Theta:

  a vector (unidimensional) or matrix (unidimensional/multidimensional)
  of latent trait values

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`extract.item`](https://philchalmers.github.io/mirt/reference/extract.item.md),
[`extract.group`](https://philchalmers.github.io/mirt/reference/extract.group.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
mod <- mirt(Science, 1)

# single item probabilty tracelines for Item 2
extr.2 <- extract.item(mod, 2)
Theta <- matrix(seq(-4,4, by = .1))
traceline <- probtrace(extr.2, Theta)
head(data.frame(traceline, Theta=Theta))
#>         P.1       P.2        P.3          P.4 Theta
#> 1 0.8786646 0.1033965 0.01717048 0.0007684150  -4.0
#> 2 0.8649759 0.1147928 0.01936274 0.0008685506  -3.9
#> 3 0.8500065 0.1271837 0.02182811 0.0009817226  -3.8
#> 4 0.8336966 0.1405950 0.02459876 0.0011096245  -3.7
#> 5 0.8159972 0.1550385 0.02771018 0.0012541689  -3.6
#> 6 0.7968730 0.1705082 0.03120137 0.0014175155  -3.5

# probability tracelines for all items in test
tracelines <- probtrace(mod, Theta)
head(tracelines)
#>      Comfort.P.1 Comfort.P.2 Comfort.P.3 Comfort.P.4  Work.P.1  Work.P.2
#> [1,]   0.3324476   0.4891306   0.1748568 0.003564956 0.8786646 0.1033965
#> [2,]   0.3097452   0.4960477   0.1902523 0.003954823 0.8649759 0.1147928
#> [3,]   0.2879244   0.5010453   0.2066432 0.004387139 0.8500065 0.1271837
#> [4,]   0.2670460   0.5040571   0.2240304 0.004866481 0.8336966 0.1405950
#> [5,]   0.2471562   0.5050429   0.2424030 0.005397913 0.8159972 0.1550385
#> [6,]   0.2282864   0.5039894   0.2617372 0.005987029 0.7968730 0.1705082
#>        Work.P.3     Work.P.4 Future.P.1 Future.P.2   Future.P.3   Future.P.4
#> [1,] 0.01717048 0.0007684150  0.9809133 0.01813820 0.0009339072 1.456036e-05
#> [2,] 0.01936274 0.0008685506  0.9761110 0.02269637 0.0011743453 1.831346e-05
#> [3,] 0.02182811 0.0009817226  0.9701371 0.02836330 0.0014765907 2.303394e-05
#> [4,] 0.02459876 0.0011096245  0.9627263 0.03538820 0.0018564770 2.897114e-05
#> [5,] 0.02771018 0.0012541689  0.9535646 0.04406509 0.0023338620 3.643864e-05
#> [6,] 0.03120137 0.0014175155  0.9422860 0.05473458 0.0029336324 4.583085e-05
#>      Benefit.P.1 Benefit.P.2 Benefit.P.3 Benefit.P.4
#> [1,]   0.7372533   0.2300751  0.03036097 0.002310633
#> [2,]   0.7155002   0.2481850  0.03373746 0.002577309
#> [3,]   0.6926969   0.2669559  0.03747256 0.002874673
#> [4,]   0.6689116   0.2862818  0.04160042 0.003206237
#> [5,]   0.6442308   0.3060358  0.04615750 0.003575906
#> [6,]   0.6187588   0.3260706  0.05118258 0.003988026
```
