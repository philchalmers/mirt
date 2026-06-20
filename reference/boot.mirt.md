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
#>        original       bias    std. error
#> t1*   1.0417547  0.042558430   0.2214854
#> t2*   4.8641542  0.111021533   0.4309400
#> t3*   2.6399417  0.012306432   0.2973376
#> t4*  -1.4660135 -0.044498186   0.1805614
#> t5*   1.2259618 -0.010047684   0.1668678
#> t6*   2.9240027  0.058332101   0.2819103
#> t7*   0.9011651 -0.011563287   0.1472493
#> t8*  -2.2665647  0.016315466   0.2186385
#> t9*   2.2933717 -0.001076186   0.4712102
#> t10*  5.2339928 -0.005911719   0.8133380
#> t11*  2.2137728 -0.016181943   0.4009778
#> t12* -1.9637062 -0.027475491   0.3737008
#> t13*  1.0949151  0.012635282   0.2126180
#> t14*  3.3479196 -0.041300089   0.3142282
#> t15*  0.9916289 -0.020732269   0.1858147
#> t16* -1.6882599 -0.012920563   0.1573070

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
#> t1* 0.5220496  0.005531384  0.08303446
#> t2* 0.5844686 -0.006122249  0.06187720
#> t3* 0.8030199 -0.006558949  0.07143927
#> t4* 0.5410276  0.011365206  0.08398957

# }
```
