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
#> Warning: EM cycles terminated after 500 iterations.
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
#> t1*   1.0417547  1.740444e-02   0.1904304
#> t2*   4.8641542  1.343367e-01   0.4229983
#> t3*   2.6399417  5.530817e-02   0.2159828
#> t4*  -1.4660135 -3.491278e-02   0.1223034
#> t5*   1.2259618  2.723756e-02   0.2371705
#> t6*   2.9240027  4.155927e-02   0.3352348
#> t7*   0.9011651 -2.975668e-03   0.1389067
#> t8*  -2.2665647 -4.463624e-02   0.2169337
#> t9*   2.2933717  7.284125e-05   0.5191957
#> t10*  5.2339928  6.762286e-02   0.8326838
#> t11*  2.2137728 -5.567994e-02   0.4026079
#> t12* -1.9637062 -6.846435e-02   0.3752243
#> t13*  1.0949151 -3.698904e-02   0.1979288
#> t14*  3.3479196 -1.755712e-02   0.2532703
#> t15*  0.9916289 -6.915037e-02   0.1036222
#> t16* -1.6882599  1.739899e-02   0.1837708

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
#>      original       bias    std. error
#> t1* 0.5220496  0.004664565  0.08413149
#> t2* 0.5844686 -0.006647398  0.06021513
#> t3* 0.8030199 -0.008141530  0.06939387
#> t4* 0.5410276  0.010047291  0.08117562

# }
```
