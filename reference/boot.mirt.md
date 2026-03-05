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

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
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
#>        original       bias    std. error
#> t1*   1.0417547  0.058693159   0.1984664
#> t2*   4.8641542  0.150294832   0.4054832
#> t3*   2.6399417  0.037911356   0.2374887
#> t4*  -1.4660135 -0.041767010   0.1077305
#> t5*   1.2259618  0.029002150   0.2097267
#> t6*   2.9240027  0.050212886   0.2951558
#> t7*   0.9011651 -0.009026432   0.1650976
#> t8*  -2.2665647 -0.042071689   0.2048568
#> t9*   2.2933717  0.016950789   0.5777651
#> t10*  5.2339928  0.048281636   0.9918184
#> t11*  2.2137728 -0.061069666   0.3767135
#> t12* -1.9637062 -0.057113034   0.3837490
#> t13*  1.0949151 -0.018524664   0.1797097
#> t14*  3.3479196 -0.031719694   0.3065423
#> t15*  0.9916289 -0.064023193   0.1268906
#> t16* -1.6882599 -0.014627676   0.1691217

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
#> Warning: EM cycles terminated after 500 iterations.
#> Warning: EM cycles terminated after 500 iterations.
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
#>      original       bias    std. error
#> t1* 0.5220496  0.004550328  0.08465327
#> t2* 0.5844686 -0.005786251  0.05942348
#> t3* 0.8030199 -0.009600268  0.06967140
#> t4* 0.5410276  0.014923169  0.08041202

# }
```
