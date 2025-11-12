# Create a user defined group-level object with correct generic functions

Initializes the proper S4 class and methods necessary for mirt functions
to use in estimation for defining customized group-level functions. To
use the defined objects pass to the `mirt(..., customGroup = OBJECT)`
command, and ensure that the class parameters are properly labelled.

## Usage

``` r
createGroup(
  par,
  est,
  den,
  nfact,
  standardize = FALSE,
  gr = NULL,
  hss = NULL,
  gen = NULL,
  lbound = NULL,
  ubound = NULL,
  derivType = "Richardson"
)
```

## Arguments

- par:

  a named vector of the starting values for the parameters

- est:

  a logical vector indicating which parameters should be freely
  estimated by default

- den:

  the probability density function given the Theta/ability values. First
  input contains a vector of all the defined parameters and the second
  input must be a matrix called `Theta`. Function also must return a
  `numeric` vector object corresponding to the associated densities for
  each row in the `Theta` input

- nfact:

  number of factors required for the model. E.g., for unidimensional
  models with only one dimension of integration `nfact = 1`

- standardize:

  logical; use standardization of the quadrature table method proposed
  by Woods and Thissen (2006)? If TRUE, the logical elements named
  `'MEAN_1'` and `'COV_11'` can be included in the parameter vector, and
  when these values are set to FALSE in the `est` input the E-table will
  be standardized to these fixed values (e.g.,
  `par <- c(a1=1, d=0, MEAN_1=0, COV_11=1)` with
  `est <- c(TRUE, TRUE, FALSE, FALSE)` will standardize the E-table to
  have a 0 mean and unit variance)

- gr:

  gradient function (vector of first derivatives) of the log-likelihood
  used in estimation. The function must be of the form `gr(x, Theta)`,
  where `x` is the object defined by `createGroup()` and `Theta` is a
  matrix of latent trait parameters

- hss:

  Hessian function (matrix of second derivatives) of the log-likelihood
  used in estimation. If not specified a numeric approximation will be
  used. The input is identical to the `gr` argument

- gen:

  a function used when `GenRandomPars = TRUE` is passed to the
  estimation function to generate random starting values. Function must
  be of the form `function(object) ...` and must return a vector with
  properties equivalent to the `par` object. If NULL, parameters will
  remain at the defined starting values by default

- lbound:

  optional vector indicating the lower bounds of the parameters. If not
  specified then the bounds will be set to -Inf

- ubound:

  optional vector indicating the lower bounds of the parameters. If not
  specified then the bounds will be set to Inf

- derivType:

  if the `gr` or `hss` terms are not specified this type will be used to
  obtain them numerically. Default is 'Richardson'

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# normal density example, N(mu, sigma^2)
den <- function(obj, Theta) dnorm(Theta, obj@par[1], sqrt(obj@par[2]))
par <- c(mu = 0, sigma2 = .5)
est <- c(FALSE, TRUE)
lbound <- c(-Inf, 0)
grp <- createGroup(par, est, den, nfact = 1, lbound=lbound)

dat <- expand.table(LSAT6)
mod <- mirt(dat, 1, 'Rasch')
modcustom <- mirt(dat, 1, 'Rasch', customGroup=grp)

coef(mod)
#> $Item_1
#>     a1     d g u
#> par  1 2.731 0 1
#> 
#> $Item_2
#>     a1     d g u
#> par  1 0.999 0 1
#> 
#> $Item_3
#>     a1    d g u
#> par  1 0.24 0 1
#> 
#> $Item_4
#>     a1     d g u
#> par  1 1.307 0 1
#> 
#> $Item_5
#>     a1   d g u
#> par  1 2.1 0 1
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0  0.572
#> 
coef(modcustom)
#> $Item_1
#>     a1     d g u
#> par  1 2.729 0 1
#> 
#> $Item_2
#>     a1     d g u
#> par  1 0.998 0 1
#> 
#> $Item_3
#>     a1    d g u
#> par  1 0.24 0 1
#> 
#> $Item_4
#>     a1     d g u
#> par  1 1.306 0 1
#> 
#> $Item_5
#>     a1     d g u
#> par  1 2.099 0 1
#> 
#> $GroupPars
#>     mu sigma2
#> par  0  0.569
#> 
```
