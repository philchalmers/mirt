# Summary of model object

Transforms coefficients into a standardized factor loading's metric. For
`MixedClass` objects, the fixed and random coefficients are printed.
Note that while the output to the console is rounded to three digits,
the returned list of objects is not. For simulations, use
`output <- summary(mod, verbose = FALSE)` to suppress the console
messages.

## Usage

``` r
# S4 method for class 'SingleGroupClass'
summary(
  object,
  SE = TRUE,
  rotate = "oblimin",
  Target = NULL,
  suppress = 0,
  suppress.cor = 0,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  an object of class `SingleGroupClass`, `MultipleGroupClass`, or
  `MixedClass`

- SE:

  logical; include the standard errors for the standardized loadings?
  Requires the initial model to have included and estimated of the
  asymptotic covariance matrix (via, for instance,
  `mirt(..., SE = TRUE)`). If `TRUE` SEs are computed using the delta
  method

- rotate:

  a string indicating which rotation to use for exploratory models,
  primarily from the `GPArotation` package (see documentation therein).

  Rotations currently supported are: `'promax'`, `'oblimin'`,
  `'varimax'`, `'quartimin'`, `'targetT'`, `'targetQ'`, `'pstT'`,
  `'pstQ'`, `'oblimax'`, `'entropy'`, `'quartimax'`, `'simplimax'`,
  `'bentlerT'`, `'bentlerQ'`, `'tandemI'`, `'tandemII'`, `'geominT'`,
  `'geominQ'`, `'cfT'`, `'cfQ'`, `'infomaxT'`, `'infomaxQ'`,
  `'mccammon'`, `'bifactorT'`, `'bifactorQ'`.

  For models that are not exploratory this input will automatically be
  set to `'none'`

- Target:

  a dummy variable matrix indicting a target rotation pattern. This is
  required for rotations such as `'targetT'`, `'targetQ'`, `'pstT'`, and
  `'pstQ'`

- suppress:

  a numeric value indicating which (possibly rotated) factor loadings
  should be suppressed. Typical values are around .3 in most statistical
  software. Default is 0 for no suppression

- suppress.cor:

  same as `suppress`, but for the correlation matrix output

- verbose:

  logical; allow information to be printed to the console?

- ...:

  additional arguments to be passed

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`coef-method`](https://philchalmers.github.io/mirt/reference/coef-method.md)

## Examples

``` r
# \donttest{
x <- mirt(Science, 2)
summary(x)
#> 
#> Rotation:  oblimin 
#> 
#> Rotated factor loadings: 
#> 
#>             F1     F2    h2
#> Comfort  0.602  0.031 0.382
#> Work    -0.057  0.797 0.592
#> Future   0.330  0.515 0.548
#> Benefit  0.723 -0.024 0.506
#> 
#> Rotated SS loadings:  0.997 0.902 
#> 
#> Factor correlations: 
#> 
#>       F1 F2
#> F1 1.000   
#> F2 0.511  1
summary(x, rotate = 'varimax')
#> 
#> Rotation:  varimax 
#> 
#> Rotated factor loadings: 
#> 
#>            F1    F2    h2
#> Comfort 0.579 0.216 0.382
#> Work    0.121 0.760 0.592
#> Future  0.428 0.605 0.548
#> Benefit 0.683 0.200 0.506
#> 
#> Rotated SS loadings:  0.999 1.03 
#> 
#> Factor correlations: 
#> 
#>    F1 F2
#> F1  1   
#> F2  0  1

# }
```
