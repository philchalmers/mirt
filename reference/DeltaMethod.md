# Numerical derivative version of delta method

Delta method using numerical derivatives (via
[`numerical_deriv`](https://philchalmers.github.io/mirt/reference/numerical_deriv.md))
for the provided function. Convenient when target transformation
function is easier to automate programmatically instead of using
explicit formula or math expressions. Can also be useful for checking
analytic results.

## Usage

``` r
DeltaMethod(fn, par, acov, ...)
```

## Arguments

- fn:

  a function specifying the type of transformation to make for each new
  parameter of interest. Must be of the form `fn(par, ...)`, or more
  simply `fn(par)`, and return a numeric vector with one element

- par:

  numerical vector passed to `fn(par)` (typically vector of MLEs)

- acov:

  numeric matrix for the ACOV of the MLEs

- ...:

  additional arguments passed to fn

## Value

returns a list of the transformed parameters, ACOV, and SEs

## Examples

``` r
# Slightly modified example from ?msm::deltamethod
# Multiple linear regression, E(y) = alpha + beta1 x + beta2 g
x <- 1:100
g <- rep(0:1, each=50)
y <- rnorm(100, 4*x, 5)
toy.lm <- lm(y ~ x + g)
estmean <- coef(toy.lm)
estvar <- vcov(toy.lm)

# Estimate of (1 / (b0 + b1)) and (1 / (b0 + b1 + b2))
1 / (estmean[1] + estmean[2])
#> (Intercept) 
#>   0.2355943 
1 / (estmean[1] + estmean[2] + estmean[3])
#> (Intercept) 
#>    0.158912 

if (FALSE) { # \dontrun{
## Approximate standard error
msm::deltamethod (~ 1 / (x1 + x2), estmean, estvar)
msm::deltamethod (~ 1 / (x1 + x2 + x3), estmean, estvar)
} # }

# with DeltaMethod
fn <- function(par) 1 / sum(par[1:2])
fn2 <- function(par) 1 / sum(par[1:3])
DeltaMethod(fn, estmean, estvar)$se
#> [1] 0.0615558
DeltaMethod(fn2, estmean, estvar)$se
#> [1] 0.06784448

# index argument for easier flexibility
fn <- function(par, index) 1 / sum(par[index])
DeltaMethod(fn, estmean, estvar, index=1:2)$se
#> [1] 0.0615558
DeltaMethod(fn, estmean, estvar, index=1:3)$se
#> [1] 0.06784448
```
