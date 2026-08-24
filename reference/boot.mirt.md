# Calculate bootstrapped standard errors for estimated models

Given an internal mirt object estimate the bootstrapped standard errors.
It may be beneficial to run the computations using multi-core
architecture (e.g., the `parallel` package). Parameters are organized
from the freely estimated values in `mod2values(x)` (equality
constraints will also be returned in the bootstrapped estimates).

## Usage

``` r
boot.mirt(x, R = 100, boot.fun = NULL, technical = NULL, ...)
```

## Arguments

- x:

  an estimated model object

- R:

  number of draws to use (passed to the `boot()` function)

- boot.fun:

  a user-defined function used to extract the information from the
  bootstrap fitted models. Must be of the form `boot.fun(x)`, where `x`
  is the bootstrap fitted model under investigation, and the return must
  be a numeric vector. If omitted a default function will be defined
  internally that returns the estimated parameters from the `mod`
  object, resulting in bootstrapped parameter estimate results

- technical:

  technical arguments passed to estimation engine. See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  details

- ...:

  additional arguments to be passed on to `boot(...)` and mirt's
  estimation engine

## References

Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# \donttest{

# standard
mod <- mirt(Science, 1)
booted <- boot.mirt(mod, R=20)
plot(booted)

booted
#> 
#> ORDINARY NONPARAMETRIC BOOTSTRAP
#> 
#> 
#> Call:
#> boot.mirt(x = mod, R = 20)
#> 
#> 
#> Bootstrap Statistics :
#>        original        bias    std. error
#> t1*   1.0417547  0.0549530793   0.2663407
#> t2*   4.8641542  0.1848934251   0.5600447
#> t3*   2.6399417  0.0276654119   0.2829615
#> t4*  -1.4660135 -0.0579383796   0.1825341
#> t5*   1.2259618 -0.0004286391   0.2200705
#> t6*   2.9240027  0.0565515796   0.2563031
#> t7*   0.9011651 -0.0125030274   0.1216343
#> t8*  -2.2665647 -0.0492747101   0.2514210
#> t9*   2.2933717  0.0423347984   0.4848237
#> t10*  5.2339928  0.0120627525   0.6633531
#> t11*  2.2137728  0.0415069121   0.3547517
#> t12* -1.9637062 -0.0630642237   0.3346630
#> t13*  1.0949151  0.0555070132   0.2155457
#> t14*  3.3479196  0.0253459515   0.3575720
#> t15*  0.9916289  0.0148941197   0.1592292
#> t16* -1.6882599 -0.0522167965   0.1642088

if (FALSE) { # \dontrun{
#run in parallel using snow back-end using all available cores
mod <- mirt(Science, 1)
booted <- boot.mirt(mod, parallel = 'snow', ncpus = parallel::detectCores())
booted
} # }

####
# bootstrapped CIs for standardized factor loadings
boot.fun <- function(mod){
  so <- summary(mod, verbose=FALSE)
  as.vector(so$rotF)
}

# test to see if it works before running
boot.fun(mod)
#> [1] 0.5220496 0.5844686 0.8030199 0.5410276

# run
booted.loads <- boot.mirt(mod, boot.fun=boot.fun)
#> Warning: EM cycles terminated after 500 iterations.
booted.loads
#> 
#> ORDINARY NONPARAMETRIC BOOTSTRAP
#> 
#> 
#> Call:
#> boot.mirt(x = mod, boot.fun = boot.fun)
#> 
#> 
#> Bootstrap Statistics :
#>      original        bias    std. error
#> t1* 0.5220496  0.0008750989  0.08420131
#> t2* 0.5844686 -0.0103200946  0.06588461
#> t3* 0.8030199 -0.0007947349  0.06583813
#> t4* 0.5410276  0.0066780769  0.08254764

# }
```
