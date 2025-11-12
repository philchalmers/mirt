# Print the model objects

Print model object summaries to the console.

## Usage

``` r
# S4 method for class 'SingleGroupClass'
print(x)
```

## Arguments

- x:

  an object of class `SingleGroupClass`, `MultipleGroupClass`, or
  `MixedClass`

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Examples

``` r
# \donttest{
x <- mirt(Science, 1)
print(x)
#> 
#> Call:
#> mirt(data = Science, model = 1)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within 1e-04 tolerance after 36 EM iterations.
#> mirt version: 1.45.6 
#> M-step optimizer: BFGS 
#> EM acceleration: Ramsay 
#> Number of rectangular quadrature: 61
#> Latent density type: Gaussian 
#> 
#> Log-likelihood = -1608.87
#> Estimated parameters: 16 
#> AIC = 3249.739
#> BIC = 3313.279; SABIC = 3262.512
#> G2 (239) = 213.56, p = 0.8804
#> RMSEA = 0, CFI = NaN, TLI = NaN
# }
```
