# Compute numerical derivatives

Compute numerical derivatives using forward/backward difference, central
difference, or Richardson extrapolation.

## Usage

``` r
numerical_deriv(
  par,
  f,
  ...,
  delta = 1e-05,
  gradient = TRUE,
  type = "Richardson"
)
```

## Arguments

- par:

  a vector of parameters to find partial derivative at

- f:

  the objective function being evaluated

- ...:

  additional arguments to be passed to `f`

- delta:

  the term used to perturb the `f` function. Default is 1e-5

- gradient:

  logical; compute the gradient terms? If FALSE then the Hessian is
  computed instead

- type:

  type of difference to compute. Can be either `'forward'` for the
  forward difference, `'central'` for the central difference, or
  `'Richardson'` for the Richardson extrapolation (default). Backward
  difference is achieved by supplying a negative `delta` value with
  `'forward'`. When `type = 'Richardson'`, the default value of `delta`
  is increased to `delta * 1000` for the Hessian and `delta * 10` for
  the gradient to provide a reasonable perturbation starting location
  (each `delta` is halved at each iteration).

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
f <- function(x) 3*x[1]^3 - 4*x[2]^2
par <- c(3,8)

# grad = 9 * x^2 , -8 * y
(actual <- c(9 * par[1]^2, -8 * par[2]))
#> [1]  81 -64
numerical_deriv(par, f, type = 'forward')
#> [1]  81.00027 -64.00004
numerical_deriv(par, f, type = 'central')
#> [1]  81 -64
numerical_deriv(par, f, type = 'Richardson') # default
#> [1]  81 -64

# Hessian = h11 -> 18 * x, h22 -> -8, h12 -> h21 -> 0
(actual <- matrix(c(18 * par[1], 0, 0, -8), 2, 2))
#>      [,1] [,2]
#> [1,]   54    0
#> [2,]    0   -8
numerical_deriv(par, f, type = 'forward', gradient = FALSE)
#>          [,1]      [,2]
#> [1,] 54.00011  0.000000
#> [2,]  0.00000 -7.999574
numerical_deriv(par, f, type = 'central', gradient = FALSE)
#>          [,1]      [,2]
#> [1,] 54.00004  0.000000
#> [2,]  0.00000 -7.999645
numerical_deriv(par, f, type = 'Richardson', gradient = FALSE) # default
#>      [,1] [,2]
#> [1,]   54    0
#> [2,]    0   -8

# }
```
