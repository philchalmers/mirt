# Specify model information

The `mirt.model` function scans/reads user input to specify the
confirmatory model. Item locations must be used in the specifications if
no `itemnames` argument is supplied. This is called implicitly by
estimation functions when a string is passed to the `model` argument.

## Usage

``` r
mirt.model(
  input = NULL,
  itemnames = NULL,
  file = "",
  COV = NULL,
  quiet = TRUE,
  ...
)
```

## Arguments

- input:

  input for writing out the model syntax. Can either be a string
  declaration of class character or the so-called Q-matrix or class
  `matrix` that specifies the model either with integer or logical
  values. If the Q-matrix method is chosen covariances terms can be
  specified with the `COV` input

- itemnames:

  a character vector or factor indicating the item names. If a
  data.frame or matrix object is supplied the names will be extracted
  using `colnames(itemnames)`. Supplying this input allows the syntax to
  be specified with the raw item names rather than item locations

- file:

  a input specifying an external file that declares the input.

- COV:

  a symmetric, logical matrix used to declare which covariance terms are
  estimated

- quiet:

  logical argument passed to
  [`scan()`](https://rdrr.io/r/base/scan.html) to suppress console read
  message

- ...:

  additional arguments for [`scan()`](https://rdrr.io/r/base/scan.html)

## Value

Returns a model specification object to be used in
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
[`bfactor`](https://philchalmers.github.io/mirt/reference/bfactor.md),
[`multipleGroup`](https://philchalmers.github.io/mirt/reference/multipleGroup.md),
or
[`mixedmirt`](https://philchalmers.github.io/mirt/reference/mixedmirt.md)

## Details

Factors are first named and then specify which numerical items they
affect (i.e., where the slope is not equal to 0), separated either by
commas or by - to indicate a range of items. Products between factors
may be specified by enclosing the left hand term within brackets. To
finish the declaration of a model simply enter a blank line with only a
carriage return (i.e., the 'enter' or 'return' key), or instead read in
an input version of the model syntax. The associated slopes throughout
the package label these coefficients as `a1, a2, ..., ak`, where the
associated number is assigned according to the respective order of the
defined factors.

For example, if the syntax were

`"G = 1-10 F = 1-5 A = 6-10"`

then the `G` factor would be assigned the slopes `a1` for each item, `F`
assigned the slopes `a2`, and `A` assigned the slopes `a3`. The same
principle applies to the
[`bfactor`](https://philchalmers.github.io/mirt/reference/bfactor.md)
function whereby the slopes are automatically included for the specific
factors after the general factor structure has been assigned.

There is an optional keyword for specifying the correlation between
relationships between factors called `COV`, and non-linear factor
products can be included by enclosing the product combination on the
left hand side of the declaration (e.g., `(F1*F1)` would create a
quadratic factor for `F1`).

The keywords
`CONSTRAIN, CONSTRAINB, PRIOR, FIXED, FREE, START, UBOUND, LBOUND` can
be applied to specific sub-groups in multiple-group models by included
square brackets before the = sign, where groups are separated by commas.
For example, to apply within-group equality constraints to a group
called "male", then specifying:

`CONSTRAIN [male] = (1-5, a1)`

is appropriate, while specifying the same constraints to the sub-groups
"male" and "female" would appear as

`CONSTRAIN [male, female] = (1-5, a1)`

For all other groups in the multi-group model, these within-group
equality constraints would not appear. Therefore, these bracketed group
specifications are useful when modifying priors, starting values,
between/within group equality constraints, and so on when the
specifications for each sub-group may differ.

Additionally, the use of negations can be used to omit specific groups
in the constraint specifications by prefixing the string with a `-`
operator, such as the following which applies between-group constraints
to all groups except "Group2" and "Group3":

`CONSTRAINB [-Group2, -Group3] = (1-5, a1)`

Finally, the keyword `GROUP` can be used to specify the group-level
hyper-parameter terms, such as the means and variance of the default
Gaussian distribution. For example, to set the starting value of the
variance parameter (`COV_11`) to 1.5:

`START = (GROUP, COV_11, 1.5)`

- COV:

  Specify the relationship between the latent factors. Estimating a
  correlation between factors is declared by joining the two factors
  with an asterisk (e.g., F1\*F2), or with an asterisk between three or
  more factors to estimate all the possible correlations (e.g.,
  F1\*F2\*F3). Specifications with the same factor (e.g., F1\*F1) will
  free the variance of said factor instead

- MEAN:

  A comma separated list specifying which latent factor means to freely
  estimate. E.g., `MEAN = F1, F2` will free the latent means for factors
  F1 and F2

- CONSTRAIN:

  A bracketed, comma separated list specifying equality constrains
  between items. The input format is
  `CONSTRAIN = (items, ..., parameterName(s)), (items, ..., parameterName)`.

  For example, in a single group 10-item dichotomous tests, using the
  default 2PL model, the first and last 5 item slopes (a1) can be
  constrained to be equal by using `CONSTRAIN = (1-5, a1), (6-10, a1)`,
  or some combination such as `CONSTRAIN = (1-3,4,5,a1), (6,7,8-10,a1)`.

  When constraining parameters to be equal across items with different
  parameter names, a balanced bracketed vector must be supplied. E.g.,
  setting the first slope for item 1 equal to the second slope in item 3
  would be `CONSTRAIN = (1, 3, a1, a2)`

- CONSTRAINB:

  A bracketed, comma separate list specifying equality constrains
  between groups. The input format is
  `CONSTRAINB = (items, ..., parameterName), (items, ..., parameterName)`.

  For example, in a two group 10-item dichotomous tests, using the
  default 2PL model, the first 5 item slopes (a1) can be constrained to
  be equal across both groups by using `CONSTRAINB = (1-5, a1)`, or some
  combination such as `CONSTRAINB = (1-3,4,5,a1)`

- PRIOR:

  A bracketed, comma separate list specifying prior parameter
  distributions. The input format is
  `PRIOR = (items, ..., parameterName, priorType, val1, val2), (items, ..., parameterName, priorType, val1, val2)`.
  For example, in a single group 10-item dichotomous tests, using the
  default 2PL model, defining a normal prior of N(0,2) for the first 5
  item intercepts (d) can be defined by `PRIOR = (1-5, d, norm, 0, 2)`

  Currently supported priors are of the form: `(items, norm, mean, sd)`
  for the normal/Gaussian, `(items, lnorm, log_mean, log_sd)` for
  log-normal, `(items, beta, alpha, beta)` for beta, and
  `(items, expbeta, alpha, beta)` for the beta distribution after
  applying the function
  [`plogis`](https://rdrr.io/r/stats/Logistic.html) to the input value
  (note, this is specifically for applying a beta prior to the
  lower-bound parameters in 3/4PL models)

- LBOUND:

  A bracketed, comma separate list specifying lower bounds for estimated
  parameters (used in optimizers such as `L-BFGS-B` and `nlminb`). The
  input format is
  `LBOUND = (items, ..., parameterName, value), (items, ..., parameterName, value)`.

  For example, in a single group 10-item dichotomous tests, using the
  3PL model and setting lower bounds for the 'g' parameters for the
  first 5 items to 0.2 is accomplished with `LBOUND = (1-5, g, 0.2)`

- UBOUND:

  same as LBOUND, but specifying upper bounds in estimated parameters

- START:

  A bracketed, comma separate list specifying the starting values for
  individual parameters. The input is of the form
  `(items, ..., parameterName, value)`. For instance, setting the 10th
  and 12th to 15th item slope parameters (a1) to 1.0 is specified with
  `START = (10, 12-15, a1, 1.0)`

  For more hands on control of the starting values pass the argument
  `pars = 'values'` through whatever estimation function is being used

- FIXED:

  A bracketed, comma separate list specifying which parameters should be
  fixed at their starting values (i.e., not freely estimated). The input
  is of the form `(items, ..., parameterName)`. For instance, fixing the
  10th and 12th to 15th item slope parameters (a1) is accomplished with
  `FIXED = (10, 12-15, a1)`

  For more hands on control of the estimated values pass the argument
  `pars = 'values'` through whatever estimation function is being used

- FREE:

  Equivalent to the `FIXED` input, except that parameters are freely
  estimated instead of fixed at their starting value

- NEXPLORE:

  Number of exploratory factors to extract. Usually this is not required
  because passing a numeric value to the `model` argument in the
  estimation function will generate an exploratory factor analysis
  model, however if different start values, priors, lower and upper
  bounds, etc, are desired then this input can be used

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com> and Alexander Robitzsch

## Examples

``` r
# \donttest{

# interactively through the console (not run)
#model <- mirt.model()
#  F1 = 1,2,3,4-10
#  F2 = 10-20
#  (F1*F2) = 1,2,3,4-10
#  COV = F1*F2


# Or alternatively with a string input
s <- 'F1 = 1,2,3,4-10
      F2 = 10-20
      (F1*F2) = 1,2,3,4-10
      COV = F1*F2'
model <- mirt.model(s)

# strings can also be passed to the estimation functions directly,
#   which silently calls mirt.model(). E.g., using the string above:
# mod <- mirt(data, s)


# Q-matrix specification
Q <- matrix(c(1,1,1,0,0,0,0,0,0,1,1,1), ncol=2, dimnames = list(NULL, c('Factor1', 'Factor2')))
COV <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2)
model <- mirt.model(Q, COV=COV)

## constrain various items slopes and all intercepts in single group model to be equal,
#   and use a log-normal prior for all the slopes
s <- 'F = 1-10
      CONSTRAIN = (1-3, 5, 6, a1), (1-10, d)
      PRIOR = (1-10, a1, lnorm, .2, .2)'
model <- mirt.model(s)


## constrain various items slopes and intercepts across groups for use in multipleGroup(),
#  and constrain first two slopes within 'group1' to be equal
s <- 'F = 1-10
      CONSTRAIN = (1-2, a1)
      CONSTRAINB = (1-3, 5, 6, a1), (1-10, d)'
model <- mirt.model(s)


## specify model using raw item names
data(data.read, package = 'sirt')
dat <- data.read

# syntax with variable names
mirtsyn2 <- "
       F1 = A1,B2,B3,C4
       F2 = A1-A4,C2,C4
       MEAN = F1
       COV = F1*F1, F1*F2
       CONSTRAIN=(A2-A4,a2),(A3,C2,d)
       PRIOR = (C3,A2-A4,a2,lnorm, .2, .2),(B3,d,norm,0,.0001)"
# create a mirt model
mirtmodel <- mirt.model(mirtsyn2, itemnames=dat)
# or equivalently:
# mirtmodel <- mirt.model(mirtsyn2, itemnames=colnames(dat))

# mod <- mirt(dat , mirtmodel)

# using sprintf() to functionally fill in information (useful for long tests
# or more complex specifications)
nitems <- 100
s <- sprintf('F = 1-%i
      CONSTRAIN = (%s, a1)
      CONSTRAINB = (%s, a1), (1-%i, d)',
      nitems, "1,2,4,50,100",
      paste0(1:45, collapse=','),
      nitems)
cat(s)
#> F = 1-100
#>       CONSTRAIN = (1,2,4,50,100, a1)
#>       CONSTRAINB = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45, a1), (1-100, d)
model <- mirt.model(s)

    # }
```
