# Compute latent regression fixed effect expected values

Create expected values for fixed effects parameters in latent regression
models.

## Usage

``` r
fixef(x)
```

## Arguments

- x:

  an estimated model object from the
  [`mixedmirt`](https://philchalmers.github.io/mirt/reference/mixedmirt.md)
  or [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)
  function

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Chalmers, R. P. (2015). Extended Mixed-Effects Item Response Models with
the MH-RM Algorithm. *Journal of Educational Measurement, 52*, 200-222.
[doi:10.1111/jedm.12072](https://doi.org/10.1111/jedm.12072)

## See also

[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
[`mixedmirt`](https://philchalmers.github.io/mirt/reference/mixedmirt.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

#simulate data
set.seed(1234)
N <- 1000

# covariates
X1 <- rnorm(N); X2 <- rnorm(N)
covdata <- data.frame(X1, X2)
Theta <- matrix(0.5 * X1 + -1 * X2 + rnorm(N, sd = 0.5))

#items and response data
a <- matrix(1, 20); d <- matrix(rnorm(20))
dat <- simdata(a, d, 1000, itemtype = '2PL', Theta=Theta)

#conditional model using X1 and X2 as predictors of Theta
mod1 <- mirt(dat, 1, 'Rasch', covdata=covdata, formula = ~ X1 + X2)

#latent regression fixed effects (i.e., expected values)
fe <- fixef(mod1)
head(fe)
#>              F1
#> [1,]  0.6128940
#> [2,] -0.1661763
#> [3,]  2.1652373
#> [4,] -1.8932285
#> [5,] -0.5021551
#> [6,]  2.2405531

# with mixedmirt()
mod1b <- mixedmirt(dat, covdata, 1, lr.fixed = ~ X1 + X2, fixed = ~ 0 + items)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1965, Max-Change = 0.1716, Max-Change = 0.1474, Max-Change = 0.1256, Max-Change = 0.1060, Max-Change = 0.0884, Max-Change = 0.0838, Max-Change = 0.0746, Max-Change = 0.0642, Max-Change = 0.0512, Max-Change = 0.0449, Max-Change = 0.0425, Max-Change = 0.0340, Max-Change = 0.0309, Max-Change = 0.0248, Max-Change = 0.0231, Max-Change = 0.0157, Max-Change = 0.0173, Max-Change = 0.0177, Max-Change = 0.0151, Max-Change = 0.0149, Max-Change = 0.0139, Max-Change = 0.0103, Max-Change = 0.0071, Max-Change = 0.0092, Max-Change = 0.0086, Max-Change = 0.0073, Max-Change = 0.0079, Max-Change = 0.0073, Max-Change = 0.0063, Max-Change = 0.0070, Max-Change = 0.0069, Max-Change = 0.0032, Max-Change = 0.0047, Max-Change = 0.0042, Max-Change = 0.0033, Max-Change = 0.0047, Max-Change = 0.0034, Max-Change = 0.0043, Max-Change = 0.0048, Max-Change = 0.0044, Max-Change = 0.0080, Max-Change = 0.0022, Max-Change = 0.0017, Max-Change = 0.0045, Max-Change = 0.0037, Max-Change = 0.0041, Max-Change = 0.0038, Max-Change = 0.0040, Max-Change = 0.0035, Max-Change = 0.0043, Max-Change = 0.0025, Max-Change = 0.0057, Max-Change = 0.0039, Max-Change = 0.0038, Max-Change = 0.0032, Max-Change = 0.0017, Max-Change = 0.0055, Max-Change = 0.0042, Max-Change = 0.0036, Max-Change = 0.0026, Max-Change = 0.0024, Max-Change = 0.0054, Max-Change = 0.0044, Max-Change = 0.0109, Max-Change = 0.0031, Max-Change = 0.0054, Max-Change = 0.0019, Max-Change = 0.0011, Max-Change = 0.0019, Max-Change = 0.0055, Max-Change = 0.0013, Max-Change = 0.0055, Max-Change = 0.0033, Max-Change = 0.0044, Max-Change = 0.0044, Max-Change = 0.0031, Max-Change = 0.0014, Max-Change = 0.0081, Max-Change = 0.0021, Max-Change = 0.0045, Max-Change = 0.0032, Max-Change = 0.0042, Max-Change = 0.0034, Max-Change = 0.0071, Max-Change = 0.0059, Max-Change = 0.0019, Max-Change = 0.0014, Max-Change = 0.0050, Max-Change = 0.0038, Max-Change = 0.0012, Max-Change = 0.0046, Max-Change = 0.0034, Max-Change = 0.0024, Max-Change = 0.0047, Max-Change = 0.0046, Max-Change = 0.0034, Max-Change = 0.0041, Max-Change = 0.0044, Max-Change = 0.0024, Max-Change = 0.0054, Max-Change = 0.0059, Max-Change = 0.0032, Max-Change = 0.0045, Max-Change = 0.0028, Max-Change = 0.0024, Max-Change = 0.0044, Max-Change = 0.0050, Max-Change = 0.0011, Max-Change = 0.0037, Max-Change = 0.0037, Max-Change = 0.0026, Max-Change = 0.0029, Max-Change = 0.0040, Max-Change = 0.0037, Max-Change = 0.0056, Max-Change = 0.0046, Max-Change = 0.0020, Max-Change = 0.0016, Max-Change = 0.0041, Max-Change = 0.0064, Max-Change = 0.0059, Max-Change = 0.0042, Max-Change = 0.0026, Max-Change = 0.0036, Max-Change = 0.0028, Max-Change = 0.0044, Max-Change = 0.0035, Max-Change = 0.0026, Max-Change = 0.0034, Max-Change = 0.0037, Max-Change = 0.0045, Max-Change = 0.0049, Max-Change = 0.0034, Max-Change = 0.0028, Max-Change = 0.0026, Max-Change = 0.0077, Max-Change = 0.0028, Max-Change = 0.0037, Max-Change = 0.0035, Max-Change = 0.0021, Max-Change = 0.0026, Max-Change = 0.0023, Max-Change = 0.0033, Max-Change = 0.0049, Max-Change = 0.0058, Max-Change = 0.0024, Max-Change = 0.0040, Max-Change = 0.0035, Max-Change = 0.0030, Max-Change = 0.0035, Max-Change = 0.0052, Max-Change = 0.0016, Max-Change = 0.0026, Max-Change = 0.0044, Max-Change = 0.0045, Max-Change = 0.0029, Max-Change = 0.0046, Max-Change = 0.0050, Max-Change = 0.0030, Max-Change = 0.0010, Max-Change = 0.0043, Max-Change = 0.0055, Max-Change = 0.0051, Max-Change = 0.0041, Max-Change = 0.0045, Max-Change = 0.0042, Max-Change = 0.0023, Max-Change = 0.0017, Max-Change = 0.0029, Max-Change = 0.0010, Max-Change = 0.0035, Max-Change = 0.0013, Max-Change = 0.0021, Max-Change = 0.0036, Max-Change = 0.0040, Max-Change = 0.0047, Max-Change = 0.0043, Max-Change = 0.0035, Max-Change = 0.0026, Max-Change = 0.0032, Max-Change = 0.0013, Max-Change = 0.0048, Max-Change = 0.0029, Max-Change = 0.0019, Max-Change = 0.0034, Max-Change = 0.0019, Max-Change = 0.0031, Max-Change = 0.0015, Max-Change = 0.0041, Max-Change = 0.0040, Max-Change = 0.0041, Max-Change = 0.0038, Max-Change = 0.0033, Max-Change = 0.0021, Max-Change = 0.0025, Max-Change = 0.0032, Max-Change = 0.0031, Max-Change = 0.0047, Max-Change = 0.0022, Max-Change = 0.0030, Max-Change = 0.0035, Max-Change = 0.0071, Max-Change = 0.0041, Max-Change = 0.0028, Max-Change = 0.0055, Max-Change = 0.0067, Max-Change = 0.0018, Max-Change = 0.0019, Max-Change = 0.0058, Max-Change = 0.0030, Max-Change = 0.0064, Max-Change = 0.0038, Max-Change = 0.0062, Max-Change = 0.0029, Max-Change = 0.0028, Max-Change = 0.0040, Max-Change = 0.0042, Max-Change = 0.0031, Max-Change = 0.0026, Max-Change = 0.0022, Max-Change = 0.0069, Max-Change = 0.0058, Max-Change = 0.0021, Max-Change = 0.0035, Max-Change = 0.0068, Max-Change = 0.0051, Max-Change = 0.0015, Max-Change = 0.0036, Max-Change = 0.0043, Max-Change = 0.0034, Max-Change = 0.0013, Max-Change = 0.0022, Max-Change = 0.0038, Max-Change = 0.0036, Max-Change = 0.0039, Max-Change = 0.0021, Max-Change = 0.0030, Max-Change = 0.0031, Max-Change = 0.0058, Max-Change = 0.0060, Max-Change = 0.0024, Max-Change = 0.0018, Max-Change = 0.0020, Max-Change = 0.0038, Max-Change = 0.0064, Max-Change = 0.0051, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0053, gam = 0.1057, Max-Change = 0.0033, gam = 0.0780, Max-Change = 0.0010, gam = 0.0629, Max-Change = 0.0008, gam = 0.0532, Max-Change = 0.0017, gam = 0.0464, Max-Change = 0.0011, gam = 0.0413, Max-Change = 0.0007, gam = 0.0374, Max-Change = 0.0009, gam = 0.0342, Max-Change = 0.0006
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
fe2 <- fixef(mod1b)
head(fe2)
#>              F1
#> [1,]  0.6165760
#> [2,] -0.1671332
#> [3,]  2.1737979
#> [4,] -1.8995517
#> [5,] -0.5047064
#> [6,]  2.2499237

# }
```
