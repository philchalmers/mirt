# Extract raw coefs from model object

Return a list (or data.frame) of raw item and group level coefficients.
Note that while the output to the console is rounded to three digits,
the returned list of objects is not. Hence, elements from
`cfs <- coef(mod); cfs[[1]]` will contain the non-rounded results
(useful for simulations).

## Usage

``` r
# S4 method for class 'SingleGroupClass'
coef(
  object,
  CI = 0.95,
  printSE = FALSE,
  rotate = "none",
  Target = NULL,
  IRTpars = FALSE,
  rawug = FALSE,
  as.data.frame = FALSE,
  simplify = FALSE,
  unique = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  an object of class `SingleGroupClass`, `MultipleGroupClass`, or
  `MixedClass`

- CI:

  the amount of converged used to compute confidence intervals; default
  is 95 percent confidence intervals

- printSE:

  logical; print the standard errors instead of the confidence
  intervals? When `IRTpars = TRUE` then the delta method will be used to
  compute the associated standard errors from mirt's default
  slope-intercept form

- rotate:

  see `summary` method for details. The default rotation is `'none'`

- Target:

  a dummy variable matrix indicting a target rotation pattern

- IRTpars:

  logical; convert slope intercept parameters into traditional IRT
  parameters? Only applicable to unidimensional models or models with
  simple structure (i.e., only one non-zero slope). If a suitable ACOV
  estimate was computed in the fitted model, and `printSE = FALSE`, then
  suitable CIs will be included based on the delta method (where
  applicable)

- rawug:

  logical; return the untransformed internal g and u parameters? If
  `FALSE`, g and u's are converted with the original format along with
  delta standard errors

- as.data.frame:

  logical; convert list output to a data.frame instead?

- simplify:

  logical; if all items have the same parameter names (indicating they
  are of the same class) then they are collapsed to a matrix, and a list
  of length 2 is returned containing a matrix of item parameters and
  group-level estimates

- unique:

  return the vector of uniquely estimated parameters

- verbose:

  logical; allow information to be printed to the console?

- ...:

  additional arguments to be passed

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`summary-method`](https://philchalmers.github.io/mirt/reference/summary-method.md)

## Examples

``` r
# \donttest{
dat <- expand.table(LSAT7)
x <- mirt(dat, 1)
coef(x)
#> $Item.1
#>        a1     d g u
#> par 0.988 1.856 0 1
#> 
#> $Item.2
#>        a1     d g u
#> par 1.081 0.808 0 1
#> 
#> $Item.3
#>        a1     d g u
#> par 1.706 1.804 0 1
#> 
#> $Item.4
#>        a1     d g u
#> par 0.765 0.486 0 1
#> 
#> $Item.5
#>        a1     d g u
#> par 0.736 1.855 0 1
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> 
coef(x, IRTpars = TRUE)
#> $Item.1
#>         a      b g u
#> par 0.988 -1.879 0 1
#> 
#> $Item.2
#>         a      b g u
#> par 1.081 -0.748 0 1
#> 
#> $Item.3
#>         a      b g u
#> par 1.706 -1.058 0 1
#> 
#> $Item.4
#>         a      b g u
#> par 0.765 -0.635 0 1
#> 
#> $Item.5
#>         a     b g u
#> par 0.736 -2.52 0 1
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> 
coef(x, simplify = TRUE)
#> $items
#>           a1     d g u
#> Item.1 0.988 1.856 0 1
#> Item.2 1.081 0.808 0 1
#> Item.3 1.706 1.804 0 1
#> Item.4 0.765 0.486 0 1
#> Item.5 0.736 1.855 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 

#with computed information matrix
x <- mirt(dat, 1, SE = TRUE)
coef(x)
#> $Item.1
#>            a1     d  g  u
#> par     0.988 1.856  0  1
#> CI_2.5  0.641 1.598 NA NA
#> CI_97.5 1.335 2.114 NA NA
#> 
#> $Item.2
#>            a1     d  g  u
#> par     1.081 0.808  0  1
#> CI_2.5  0.750 0.629 NA NA
#> CI_97.5 1.412 0.987 NA NA
#> 
#> $Item.3
#>            a1     d  g  u
#> par     1.706 1.804  0  1
#> CI_2.5  1.078 1.404 NA NA
#> CI_97.5 2.334 2.205 NA NA
#> 
#> $Item.4
#>            a1     d  g  u
#> par     0.765 0.486  0  1
#> CI_2.5  0.502 0.339 NA NA
#> CI_97.5 1.028 0.633 NA NA
#> 
#> $Item.5
#>            a1     d  g  u
#> par     0.736 1.855  0  1
#> CI_2.5  0.440 1.630 NA NA
#> CI_97.5 1.032 2.079 NA NA
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0      1
#> CI_2.5      NA     NA
#> CI_97.5     NA     NA
#> 
coef(x, printSE = TRUE)
#> $Item.1
#>        a1     d logit(g) logit(u)
#> par 0.988 1.856     -999      999
#> SE  0.177 0.131       NA       NA
#> 
#> $Item.2
#>        a1     d logit(g) logit(u)
#> par 1.081 0.808     -999      999
#> SE  0.169 0.091       NA       NA
#> 
#> $Item.3
#>        a1     d logit(g) logit(u)
#> par 1.706 1.804     -999      999
#> SE  0.320 0.204       NA       NA
#> 
#> $Item.4
#>        a1     d logit(g) logit(u)
#> par 0.765 0.486     -999      999
#> SE  0.134 0.075       NA       NA
#> 
#> $Item.5
#>        a1     d logit(g) logit(u)
#> par 0.736 1.855     -999      999
#> SE  0.151 0.114       NA       NA
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> SE      NA     NA
#> 
coef(x, as.data.frame = TRUE)
#>                        par    CI_2.5   CI_97.5
#> Item.1.a1        0.9879254 0.6405319 1.3353189
#> Item.1.d         1.8560605 1.5983450 2.1137759
#> Item.1.g         0.0000000        NA        NA
#> Item.1.u         1.0000000        NA        NA
#> Item.2.a1        1.0808847 0.7500334 1.4117360
#> Item.2.d         0.8079786 0.6291264 0.9868309
#> Item.2.g         0.0000000        NA        NA
#> Item.2.u         1.0000000        NA        NA
#> Item.3.a1        1.7058006 1.0778209 2.3337803
#> Item.3.d         1.8042187 1.4035692 2.2048683
#> Item.3.g         0.0000000        NA        NA
#> Item.3.u         1.0000000        NA        NA
#> Item.4.a1        0.7651853 0.5022681 1.0281025
#> Item.4.d         0.4859966 0.3391601 0.6328331
#> Item.4.g         0.0000000        NA        NA
#> Item.4.u         1.0000000        NA        NA
#> Item.5.a1        0.7357980 0.4395386 1.0320574
#> Item.5.d         1.8545127 1.6302516 2.0787739
#> Item.5.g         0.0000000        NA        NA
#> Item.5.u         1.0000000        NA        NA
#> GroupPars.MEAN_1 0.0000000        NA        NA
#> GroupPars.COV_11 1.0000000        NA        NA

#two factors
x2 <- mirt(Science, 2)
coef(x2)
#> $Comfort
#>         a1    a2    d1    d2     d3
#> par -1.335 0.097 5.211 2.866 -1.603
#> 
#> $Work
#>         a1    a2    d1    d2     d3
#> par -0.879 1.853 3.704 1.153 -2.904
#> 
#> $Future
#>        a1    a2    d1    d2     d3
#> par -1.47 1.165 4.663 1.957 -1.736
#> 
#> $Benefit
#>         a1 a2    d1    d2     d3
#> par -1.722  0 3.989 1.195 -2.044
#> 
#> $GroupPars
#>     MEAN_1 MEAN_2 COV_11 COV_21 COV_22
#> par      0      0      1      0      1
#> 
coef(x2, rotate = 'varimax')
#> 
#> Rotation:  varimax 
#> 
#> $Comfort
#>        a1    a2    d1    d2     d3
#> par 1.254 0.468 5.211 2.866 -1.603
#> 
#> $Work
#>        a1    a2    d1    d2     d3
#> par 0.323 2.025 3.704 1.153 -2.904
#> 
#> $Future
#>        a1    a2    d1    d2     d3
#> par 1.083 1.531 4.663 1.957 -1.736
#> 
#> $Benefit
#>        a1    a2    d1    d2     d3
#> par 1.653 0.484 3.989 1.195 -2.044
#> 
#> $GroupPars
#>     MEAN_1 MEAN_2 COV_11 COV_21 COV_22
#> par      0      0      1      0      1
#> 

# }
```
