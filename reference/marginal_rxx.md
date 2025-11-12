# Function to calculate the marginal reliability

Given an estimated model and a prior density function, compute the
marginal reliability (Thissen and Wainer, 2001). This is only available
for unidimensional tests.

## Usage

``` r
marginal_rxx(mod, density = dnorm, ...)
```

## Arguments

- mod:

  an object of class `'SingleGroupClass'`

- density:

  a density function to use for integration. Default assumes the latent
  traits are from a normal (Gaussian) distribution

- ...:

  additional arguments passed to the density function

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Thissen, D. and Wainer, H. (2001). Test Scoring. Lawrence Erlbaum
Associates.

## See also

[`empirical_rxx`](https://philchalmers.github.io/mirt/reference/empirical_rxx.md),
[`extract.group`](https://philchalmers.github.io/mirt/reference/extract.group.md),
[`testinfo`](https://philchalmers.github.io/mirt/reference/testinfo.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

dat <- expand.table(deAyala)
mod <- mirt(dat)

# marginal estimate treating item parameters as known
marginal_rxx(mod)
#> [1] 0.6092894

# compare to alpha
itemstats(dat)$overall$alpha
#> [1] 0.6077281

# \donttest{

# empirical estimate (assuming the same prior)
fscores(mod, returnER = TRUE)
#>        F1 
#> 0.6200703 

# empirical rxx the alternative way, given theta scores and SEs
fs <- fscores(mod, full.scores.SE=TRUE)
head(fs)
#>             F1     SE_F1
#> [1,] -1.580133 0.6699424
#> [2,] -1.580133 0.6699424
#> [3,] -1.580133 0.6699424
#> [4,] -1.580133 0.6699424
#> [5,] -1.580133 0.6699424
#> [6,] -1.580133 0.6699424
empirical_rxx(fs)
#>        F1 
#> 0.6200703 

#############
# example demonstrating correlation attenuation

theta <- rnorm(1000)
X <- theta + rnorm(1000, sd=2)
cor(X, theta)    # correlation without measurement error (what you want)
#> [1] 0.4504707

# measured with a 10 item GRM test
nitems <- 10
a <- matrix(rlnorm(nitems,.2,.3))
diffs <- t(apply(matrix(runif(nitems*4, .3, 1), nitems), 1, cumsum))
diffs <- -(diffs - rowMeans(diffs))
d <- diffs + rnorm(nitems)
dat <- simdata(a, d, itemtype = 'graded', Theta=matrix(theta))

# correlation with total score (attenuated)
cor(rowSums(dat), X)
#> [1] 0.3875882

# fit single group model
mod <- mirt(dat)

# EAP correlation (also attenuated)
fs <- fscores(mod)
cor(fs, X)
#>         [,1]
#> F1 0.3901645

# correction for attenuation, r_x.theta = r_x.theta.hat / sqrt(rxx_theta.hat)
(rxx <- marginal_rxx(mod))  # alternatively, could use empirical_rxx()
#> [1] 0.7563939
cor(fs, X) / sqrt(rxx)  # correction estimate
#>         [,1]
#> F1 0.4486149
cor(X, theta)           # compare to true correlation
#> [1] 0.4504707

# }
```
