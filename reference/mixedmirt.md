# Mixed effects modeling for MIRT models

`mixedmirt` fits MIRT models using FIML estimation to dichotomous and
polytomous IRT models conditional on fixed and random effect of person
and item level covariates. This can also be understood as 'explanatory
IRT' if only fixed effects are modeled, or multilevel/mixed IRT if
random and fixed effects are included. The method uses the MH-RM
algorithm exclusively. Additionally, computation of the log-likelihood
can be sped up by using parallel estimation via
[`mirtCluster`](https://philchalmers.github.io/mirt/reference/mirtCluster.md).

## Usage

``` r
mixedmirt(
  data,
  covdata = NULL,
  model = 1,
  fixed = ~1,
  random = NULL,
  itemtype = "Rasch",
  lr.fixed = ~1,
  lr.random = NULL,
  itemdesign = NULL,
  constrain = NULL,
  pars = NULL,
  return.design = FALSE,
  SE = TRUE,
  internal_constraints = TRUE,
  technical = list(SEtol = 1e-04),
  ...
)
```

## Arguments

- data:

  a `matrix` or `data.frame` that consists of numerically ordered data,
  organized in the form of integers, with missing data coded as `NA`

- covdata:

  a `data.frame` that consists of the `nrow(data)` by `K` 'person level'
  fixed and random predictors. If missing data are present in this
  object then the observations from `covdata` and `data` will be removed
  row-wise via the `rowSums(is.na(covdata)) > 0`

- model:

  an object returned from, or a string to be passed to,
  [`mirt.model()`](https://philchalmers.github.io/mirt/reference/mirt.model.md)
  to declare how the IRT model is to be estimated. See
  [`mirt.model`](https://philchalmers.github.io/mirt/reference/mirt.model.md)
  and [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)
  for more detail

- fixed:

  a right sided R formula for specifying the fixed effect (aka
  'explanatory') predictors from `covdata` and `itemdesign`. To estimate
  the intercepts for each item the keyword `items` is reserved and
  automatically added to the `itemdesign` input. If any polytomous items
  are being model the `items` are argument is not valid since all
  intercept parameters are freely estimated and identified with the
  parameterizations found in
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md), and
  the first column in the fixed design matrix (commonly the intercept or
  a reference group) is omitted

- random:

  a right sided formula or list of formulas containing crossed random
  effects of the form `v1 + ... v_n | G`, where `G` is the grouping
  variable and `v_n` are random numeric predictors within each group. If
  no intercept value is specified then by default the correlations
  between the `v`'s and `G` are estimated, but can be suppressed by
  including the `~ -1 + ...` or 0 constant. `G` may contain interaction
  terms, such as `group:items` to include cross or person-level
  interactions effects

- itemtype:

  same as itemtype in
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
  except when the `fixed` or `random` inputs are used does not support
  the following item types:
  `c('PC2PL', 'PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')`

- lr.fixed:

  an R formula (or list of formulas) to specify regression effects in
  the latent variables from the variables in `covdata`. This is used to
  construct models such as the so-called 'latent regression model' to
  explain person-level ability/trait differences. If a named list of
  formulas is supplied (where the names correspond to the latent trait
  names in `model`) then specific regression effects can be estimated
  for each factor. Supplying a single formula will estimate the
  regression parameters for all latent traits by default.

- lr.random:

  a list of random effect terms for modeling variability in the latent
  trait scores, where the syntax uses the same style as in the `random`
  argument. Useful for building so-called 'multilevel IRT' models which
  are non-Rasch (multilevel Rasch models do not technically require
  these because they can be built using the `fixed` and `random` inputs
  alone)

- itemdesign:

  a `data.frame` object used to create a design matrix for the items,
  where each `nrow(itemdesign) == nitems` and the number of columns is
  equal to the number of fixed effect predictors (i.e., item
  intercepts). By default an `items` variable is reserved for modeling
  the item intercept parameters

- constrain:

  a list indicating parameter equality constrains. See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  more detail

- pars:

  used for parameter starting values. See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  more detail

- return.design:

  logical; return the design matrices before they have (potentially)
  been reassigned?

- SE:

  logical; compute the standard errors by approximating the information
  matrix using the MHRM algorithm? Default is TRUE

- internal_constraints:

  logical; use the internally defined constraints for constraining
  effects across persons and items? Default is TRUE. Setting this to
  FALSE runs the risk of under-identification

- technical:

  the technical list passed to the MH-RM estimation engine, with the
  SEtol default increased to .0001. Additionally, the argument
  `RANDSTART` is available to indicate at which iteration (during the
  burn-in stage) the additional random effect variables should begin to
  be approximated (i.e., elements in `lr.random` and `random`). The
  default for `RANDSTART` is to start at iteration 100, and when random
  effects are included the default number of burn-in iterations is
  increased from 150 to 200. See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  further details

- ...:

  additional arguments to be passed to the MH-RM estimation engine. See
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  more details and examples

## Value

function returns an object of class `MixedClass`
([MixedClass-class](https://philchalmers.github.io/mirt/reference/MixedClass-class.md)).

## Details

For dichotomous response models, `mixedmirt` follows the general form

\$\$P(x = 1\|\theta, \psi) = g + \frac{(u - g)}{1 + exp(-1 \* \[\theta
a + X \beta + Z \delta\])}\$\$

where X is a design matrix with associated \\\beta\\ fixed effect
intercept coefficients, and Z is a design matrix with associated
\\\delta\\ random effects for the intercepts. For simplicity and easier
interpretation, the unique item intercept values typically found in \\X
\beta\\ are extracted and reassigned within mirt's 'intercept'
parameters (e.g., `'d'`). To observe how the design matrices are
structured prior to reassignment and estimation pass the argument
`return.design = TRUE`.

Polytomous IRT models follow a similar format except the item intercepts
are automatically estimated internally, rendering the `items` argument
in the fixed formula redundant and therefore must be omitted from the
specification. If there are a mixture of dichotomous and polytomous
items the intercepts for the dichotomous models are also estimated for
consistency.

The decomposition of the \\\theta\\ parameters is also possible to form
latent regression and multilevel IRT models by using the `lr.fixed` and
`lr.random` inputs. These effects decompose \\\theta\\ such that

\$\$\theta = V \Gamma + W \zeta + \epsilon\$\$

where V and W are fixed and random effects design matrices for the
associated coefficients.

To simulate expected a posteriori predictions for the random effect
terms use the
[`randef`](https://philchalmers.github.io/mirt/reference/randef.md)
function.

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Chalmers, R. P. (2015). Extended Mixed-Effects Item Response Models with
the MH-RM Algorithm. *Journal of Educational Measurement, 52*, 200-222.
[doi:10.1111/jedm.12072](https://doi.org/10.1111/jedm.12072)

## See also

[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md),
[`randef`](https://philchalmers.github.io/mirt/reference/randef.md),
[`fixef`](https://philchalmers.github.io/mirt/reference/fixef.md),
[`boot.mirt`](https://philchalmers.github.io/mirt/reference/boot.mirt.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

# make some data
set.seed(1234)
N <- 750
a <- matrix(rlnorm(10,.3,1),10,1)
d <- matrix(rnorm(10), 10)
Theta <- matrix(sort(rnorm(N)))
pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
pseudoIQ <- (pseudoIQ - mean(pseudoIQ))/10  #rescale variable for numerical stability
group <- factor(rep(c('G1','G2','G3'), each = N/3))
data <- simdata(a,d,N, itemtype = rep('2PL',10), Theta=Theta)
covdata <- data.frame(group, pseudoIQ)

itemstats(data)
#> $overall
#>    N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  750            4.655          2.346 0.166 0.133 0.671     1.345
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  750 2 0.363 0.481   0.368         0.172       0.678
#> Item_2  750 2 0.335 0.472   0.631         0.485       0.617
#> Item_3  750 2 0.428 0.495   0.711         0.580       0.594
#> Item_4  750 2 0.512 0.500   0.234         0.021       0.708
#> Item_5  750 2 0.628 0.484   0.638         0.490       0.615
#> Item_6  750 2 0.472 0.500   0.669         0.523       0.607
#> Item_7  750 2 0.385 0.487   0.471         0.286       0.657
#> Item_8  750 2 0.301 0.459   0.483         0.312       0.651
#> Item_9  750 2 0.319 0.466   0.481         0.307       0.652
#> Item_10 750 2 0.912 0.283   0.295         0.180       0.670
#> 
#> $proportions
#>             0     1
#> Item_1  0.637 0.363
#> Item_2  0.665 0.335
#> Item_3  0.572 0.428
#> Item_4  0.488 0.512
#> Item_5  0.372 0.628
#> Item_6  0.528 0.472
#> Item_7  0.615 0.385
#> Item_8  0.699 0.301
#> Item_9  0.681 0.319
#> Item_10 0.088 0.912
#> 

# use parallel computing
if(interactive()) mirtCluster()

# specify IRT model
model <- 'Theta = 1-10'

# model with no person predictors
mod0 <- mirt(data, model, itemtype = 'Rasch')

# group as a fixed effect predictor (aka, uniform dif)
mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1804, Max-Change = 0.1582, Max-Change = 0.1359, Max-Change = 0.1172, Max-Change = 0.1000, Max-Change = 0.0848, Max-Change = 0.0698, Max-Change = 0.0570, Max-Change = 0.0459, Max-Change = 0.0370, Max-Change = 0.0293, Max-Change = 0.0246, Max-Change = 0.0201, Max-Change = 0.0183, Max-Change = 0.0156, Max-Change = 0.0103, Max-Change = 0.0106, Max-Change = 0.0098, Max-Change = 0.0093, Max-Change = 0.0055, Max-Change = 0.0062, Max-Change = 0.0108, Max-Change = 0.0065, Max-Change = 0.0076, Max-Change = 0.0137, Max-Change = 0.0067, Max-Change = 0.0036, Max-Change = 0.0023, Max-Change = 0.0095, Max-Change = 0.0108, Max-Change = 0.0141, Max-Change = 0.0083, Max-Change = 0.0064, Max-Change = 0.0056, Max-Change = 0.0151, Max-Change = 0.0066, Max-Change = 0.0049, Max-Change = 0.0097, Max-Change = 0.0052, Max-Change = 0.0064, Max-Change = 0.0080, Max-Change = 0.0028, Max-Change = 0.0076, Max-Change = 0.0055, Max-Change = 0.0063, Max-Change = 0.0047, Max-Change = 0.0077, Max-Change = 0.0048, Max-Change = 0.0053, Max-Change = 0.0095, Max-Change = 0.0085, Max-Change = 0.0055, Max-Change = 0.0095, Max-Change = 0.0115, Max-Change = 0.0062, Max-Change = 0.0076, Max-Change = 0.0043, Max-Change = 0.0053, Max-Change = 0.0035, Max-Change = 0.0064, Max-Change = 0.0078, Max-Change = 0.0051, Max-Change = 0.0046, Max-Change = 0.0057, Max-Change = 0.0020, Max-Change = 0.0050, Max-Change = 0.0058, Max-Change = 0.0027, Max-Change = 0.0093, Max-Change = 0.0032, Max-Change = 0.0042, Max-Change = 0.0052, Max-Change = 0.0036, Max-Change = 0.0018, Max-Change = 0.0045, Max-Change = 0.0024, Max-Change = 0.0042, Max-Change = 0.0082, Max-Change = 0.0086, Max-Change = 0.0041, Max-Change = 0.0053, Max-Change = 0.0045, Max-Change = 0.0170, Max-Change = 0.0099, Max-Change = 0.0101, Max-Change = 0.0036, Max-Change = 0.0036, Max-Change = 0.0038, Max-Change = 0.0076, Max-Change = 0.0062, Max-Change = 0.0065, Max-Change = 0.0065, Max-Change = 0.0056, Max-Change = 0.0076, Max-Change = 0.0084, Max-Change = 0.0055, Max-Change = 0.0083, Max-Change = 0.0023, Max-Change = 0.0065, Max-Change = 0.0067, Max-Change = 0.0066, Max-Change = 0.0082, Max-Change = 0.0057, Max-Change = 0.0032, Max-Change = 0.0022, Max-Change = 0.0035, Max-Change = 0.0064, Max-Change = 0.0047, Max-Change = 0.0093, Max-Change = 0.0073, Max-Change = 0.0025, Max-Change = 0.0040, Max-Change = 0.0090, Max-Change = 0.0089, Max-Change = 0.0055, Max-Change = 0.0073, Max-Change = 0.0068, Max-Change = 0.0070, Max-Change = 0.0012, Max-Change = 0.0116, Max-Change = 0.0067, Max-Change = 0.0072, Max-Change = 0.0057, Max-Change = 0.0018, Max-Change = 0.0095, Max-Change = 0.0058, Max-Change = 0.0035, Max-Change = 0.0024, Max-Change = 0.0059, Max-Change = 0.0059, Max-Change = 0.0058, Max-Change = 0.0080, Max-Change = 0.0033, Max-Change = 0.0074, Max-Change = 0.0073, Max-Change = 0.0078, Max-Change = 0.0066, Max-Change = 0.0064, Max-Change = 0.0035, Max-Change = 0.0048, Max-Change = 0.0064, Max-Change = 0.0050, Max-Change = 0.0039, Max-Change = 0.0038, Max-Change = 0.0065, Max-Change = 0.0068, Max-Change = 0.0092, Max-Change = 0.0043, Max-Change = 0.0068, Max-Change = 0.0037, Max-Change = 0.0050, Max-Change = 0.0054, Max-Change = 0.0032, Max-Change = 0.0052, Max-Change = 0.0028, Max-Change = 0.0067, Max-Change = 0.0052, Max-Change = 0.0035, Max-Change = 0.0057, Max-Change = 0.0049, Max-Change = 0.0056, Max-Change = 0.0030, Max-Change = 0.0007, Max-Change = 0.0016, Max-Change = 0.0075, Max-Change = 0.0040, Max-Change = 0.0058, Max-Change = 0.0023, Max-Change = 0.0035, Max-Change = 0.0011, Max-Change = 0.0046, Max-Change = 0.0061, Max-Change = 0.0018, Max-Change = 0.0095, Max-Change = 0.0076, Max-Change = 0.0073, Max-Change = 0.0022, Max-Change = 0.0035, Max-Change = 0.0048, Max-Change = 0.0053, Max-Change = 0.0102, Max-Change = 0.0054, Max-Change = 0.0048, Max-Change = 0.0028, Max-Change = 0.0055, Max-Change = 0.0081, Max-Change = 0.0010, Max-Change = 0.0050, Max-Change = 0.0072, Max-Change = 0.0047, Max-Change = 0.0100, Max-Change = 0.0074, Max-Change = 0.0054, Max-Change = 0.0040, Max-Change = 0.0051, Max-Change = 0.0031, Max-Change = 0.0062, Max-Change = 0.0077, Max-Change = 0.0071, Max-Change = 0.0033, Max-Change = 0.0062, Max-Change = 0.0026, Max-Change = 0.0040, Max-Change = 0.0086, Max-Change = 0.0046, Max-Change = 0.0067, Max-Change = 0.0048, Max-Change = 0.0070, Max-Change = 0.0062, Max-Change = 0.0078, Max-Change = 0.0077, Max-Change = 0.0041, Max-Change = 0.0018, Max-Change = 0.0074, Max-Change = 0.0023, Max-Change = 0.0036, Max-Change = 0.0037, Max-Change = 0.0045, Max-Change = 0.0059, Max-Change = 0.0057, Max-Change = 0.0042, Max-Change = 0.0059, Max-Change = 0.0054, Max-Change = 0.0096, Max-Change = 0.0061, Max-Change = 0.0041, Max-Change = 0.0048, Max-Change = 0.0029, Max-Change = 0.0087, Max-Change = 0.0017, Max-Change = 0.0022, Max-Change = 0.0036, Max-Change = 0.0031, Max-Change = 0.0031, Max-Change = 0.0025, Max-Change = 0.0085, Max-Change = 0.0045, Max-Change = 0.0064, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0048, gam = 0.1057, Max-Change = 0.0032, gam = 0.0780, Max-Change = 0.0010, gam = 0.0629, Max-Change = 0.0022, gam = 0.0532, Max-Change = 0.0012, gam = 0.0464, Max-Change = 0.0022, gam = 0.0413, Max-Change = 0.0008, gam = 0.0374, Max-Change = 0.0012, gam = 0.0342, Max-Change = 0.0016, gam = 0.0316, Max-Change = 0.0006, gam = 0.0294, Max-Change = 0.0006, gam = 0.0276, Max-Change = 0.0006
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
anova(mod0, mod1)
#>           AIC    SABIC       HQ      BIC    logLik      X2 df p
#> mod0 8799.543 8815.435 8819.126 8850.364 -4388.772             
#> mod1 8111.326 8130.107 8134.469 8171.387 -4042.663 692.217  2 0
summary(mod1)
#> 
#> Call:
#> mixedmirt(data = data, covdata = covdata, model = model, fixed = ~0 + 
#>     group + items)
#> 
#> --------------
#> FIXED EFFECTS:
#>         Estimate Std.Error z.value
#> groupG1   -1.858     0.100 -18.495
#> groupG2   -0.748     0.094  -7.979
#> groupG3    0.515     0.094   5.493
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>       Theta
#> Theta 0.118
#> 
coef(mod1)
#> $Item_1
#>         groupG1 groupG2 groupG3 a1  d  g  u
#> par      -1.858  -0.748   0.515  1  0  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA NA NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA NA NA NA
#> 
#> $Item_2
#>         groupG1 groupG2 groupG3 a1      d  g  u
#> par      -1.858  -0.748   0.515  1 -0.152  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA -0.388 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA  0.084 NA NA
#> 
#> $Item_3
#>         groupG1 groupG2 groupG3 a1     d  g  u
#> par      -1.858  -0.748   0.515  1 0.340  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA 0.109 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA 0.572 NA NA
#> 
#> $Item_4
#>         groupG1 groupG2 groupG3 a1     d  g  u
#> par      -1.858  -0.748   0.515  1 0.763  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA 0.532 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA 0.994 NA NA
#> 
#> $Item_5
#>         groupG1 groupG2 groupG3 a1     d  g  u
#> par      -1.858  -0.748   0.515  1 1.353  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA 1.119 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA 1.588 NA NA
#> 
#> $Item_6
#>         groupG1 groupG2 groupG3 a1     d  g  u
#> par      -1.858  -0.748   0.515  1 0.563  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA 0.332 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA 0.793 NA NA
#> 
#> $Item_7
#>         groupG1 groupG2 groupG3 a1      d  g  u
#> par      -1.858  -0.748   0.515  1  0.120  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA -0.113 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA  0.353 NA NA
#> 
#> $Item_8
#>         groupG1 groupG2 groupG3 a1      d  g  u
#> par      -1.858  -0.748   0.515  1 -0.339  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA -0.578 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA -0.101 NA NA
#> 
#> $Item_9
#>         groupG1 groupG2 groupG3 a1      d  g  u
#> par      -1.858  -0.748   0.515  1 -0.241  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA -0.478 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA -0.004 NA NA
#> 
#> $Item_10
#>         groupG1 groupG2 groupG3 a1     d  g  u
#> par      -1.858  -0.748   0.515  1 3.432  0  1
#> CI_2.5   -2.055  -0.931   0.331 NA 3.117 NA NA
#> CI_97.5  -1.661  -0.564   0.699 NA 3.747 NA NA
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0  0.118
#> CI_2.5      NA  0.072
#> CI_97.5     NA  0.164
#> 

if (FALSE) { # \dontrun{
# same model as above in lme4
wide <- data.frame(id=1:nrow(data),data,covdata)
long <- reshape2::melt(wide, id.vars = c('id', 'group', 'pseudoIQ'))
library(lme4)
lmod0 <- glmer(value ~ 0 + variable + (1|id), long, family = binomial)
lmod1 <- glmer(value ~ 0 + group + variable + (1|id), long, family = binomial)
anova(lmod0, lmod1)
} # }

# model using 2PL items instead of Rasch
mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items, itemtype = '2PL')
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1860, Max-Change = 0.1696, Max-Change = 0.1434, Max-Change = 0.1281, Max-Change = 0.1024, Max-Change = 0.0980, Max-Change = 0.0816, Max-Change = 0.0806, Max-Change = 0.0636, Max-Change = 0.0644, Max-Change = 0.0457, Max-Change = 0.0484, Max-Change = 0.0573, Max-Change = 0.0544, Max-Change = 0.0444, Max-Change = 0.0289, Max-Change = 0.0242, Max-Change = 0.0313, Max-Change = 0.0332, Max-Change = 0.0317, Max-Change = 0.0343, Max-Change = 0.0348, Max-Change = 0.0514, Max-Change = 0.0715, Max-Change = 0.0273, Max-Change = 0.0443, Max-Change = 0.0254, Max-Change = 0.0261, Max-Change = 0.0364, Max-Change = 0.0200, Max-Change = 0.0392, Max-Change = 0.0491, Max-Change = 0.0365, Max-Change = 0.0491, Max-Change = 0.0339, Max-Change = 0.0376, Max-Change = 0.0239, Max-Change = 0.0464, Max-Change = 0.0513, Max-Change = 0.0378, Max-Change = 0.0793, Max-Change = 0.0370, Max-Change = 0.0447, Max-Change = 0.0359, Max-Change = 0.0696, Max-Change = 0.0790, Max-Change = 0.0264, Max-Change = 0.0346, Max-Change = 0.0233, Max-Change = 0.0334, Max-Change = 0.0345, Max-Change = 0.0441, Max-Change = 0.0441, Max-Change = 0.0380, Max-Change = 0.0255, Max-Change = 0.0391, Max-Change = 0.0394, Max-Change = 0.0293, Max-Change = 0.0720, Max-Change = 0.0462, Max-Change = 0.0417, Max-Change = 0.0285, Max-Change = 0.0245, Max-Change = 0.0603, Max-Change = 0.0343, Max-Change = 0.0166, Max-Change = 0.0561, Max-Change = 0.0314, Max-Change = 0.0257, Max-Change = 0.0314, Max-Change = 0.0342, Max-Change = 0.0598, Max-Change = 0.0402, Max-Change = 0.0349, Max-Change = 0.0471, Max-Change = 0.0189, Max-Change = 0.0324, Max-Change = 0.0257, Max-Change = 0.0750, Max-Change = 0.0510, Max-Change = 0.0721, Max-Change = 0.0339, Max-Change = 0.0547, Max-Change = 0.0469, Max-Change = 0.0607, Max-Change = 0.0375, Max-Change = 0.0417, Max-Change = 0.0378, Max-Change = 0.0610, Max-Change = 0.0517, Max-Change = 0.0244, Max-Change = 0.0349, Max-Change = 0.0239, Max-Change = 0.0503, Max-Change = 0.0232, Max-Change = 0.0854, Max-Change = 0.0282, Max-Change = 0.0397, Max-Change = 0.0406, Max-Change = 0.0820, Max-Change = 0.0267, Max-Change = 0.0839, Max-Change = 0.0307, Max-Change = 0.0237, Max-Change = 0.0285, Max-Change = 0.1084, Max-Change = 0.1160, Max-Change = 0.0390, Max-Change = 0.0355, Max-Change = 0.0703, Max-Change = 0.0391, Max-Change = 0.0954, Max-Change = 0.0553, Max-Change = 0.0308, Max-Change = 0.0652, Max-Change = 0.0169, Max-Change = 0.0664, Max-Change = 0.1173, Max-Change = 0.0311, Max-Change = 0.0509, Max-Change = 0.1200, Max-Change = 0.0904, Max-Change = 0.1736, Max-Change = 0.0696, Max-Change = 0.0348, Max-Change = 0.0707, Max-Change = 0.0487, Max-Change = 0.0574, Max-Change = 0.0342, Max-Change = 0.0850, Max-Change = 0.0477, Max-Change = 0.0656, Max-Change = 0.0731, Max-Change = 0.0842, Max-Change = 0.1089, Max-Change = 0.0266, Max-Change = 0.1092, Max-Change = 0.0731, Max-Change = 0.0200, Max-Change = 0.0487, Max-Change = 0.0573, Max-Change = 0.1227, Max-Change = 0.0633, Max-Change = 0.0499, Max-Change = 0.0307, Max-Change = 0.1556, Max-Change = 0.0620, Max-Change = 0.0445, Max-Change = 0.1713, Max-Change = 0.0844, Max-Change = 0.0719, Max-Change = 0.0258, Max-Change = 0.0742, Max-Change = 0.0565, Max-Change = 0.0339, Max-Change = 0.0182, Max-Change = 0.0373, Max-Change = 0.1227, Max-Change = 0.0404, Max-Change = 0.0894, Max-Change = 0.0172, Max-Change = 0.1133, Max-Change = 0.0502, Max-Change = 0.0338, Max-Change = 0.0758, Max-Change = 0.0635, Max-Change = 0.0280, Max-Change = 0.0781, Max-Change = 0.0398, Max-Change = 0.0805, Max-Change = 0.1266, Max-Change = 0.1815, Max-Change = 0.0498, Max-Change = 0.0387, Max-Change = 0.0311, Max-Change = 0.0463, Max-Change = 0.1127, Max-Change = 0.0519, Max-Change = 0.1120, Max-Change = 0.2000, Max-Change = 0.0178, Max-Change = 0.0509, Max-Change = 0.1411, Max-Change = 0.1258, Max-Change = 0.0486, Max-Change = 0.0713, Max-Change = 0.0751, Max-Change = 0.1150, Max-Change = 0.0370, Max-Change = 0.1355, Max-Change = 0.2000, Max-Change = 0.1107, Max-Change = 0.0404, Max-Change = 0.1141, Max-Change = 0.1374, Max-Change = 0.1567, Max-Change = 0.0502, Max-Change = 0.0305, Max-Change = 0.0899, Max-Change = 0.0245, Max-Change = 0.0828, Max-Change = 0.1181, Max-Change = 0.0241, Max-Change = 0.0989, Max-Change = 0.1326, Max-Change = 0.0390, Max-Change = 0.1871, Max-Change = 0.0252, Max-Change = 0.0234, Max-Change = 0.0791, Max-Change = 0.1058, Max-Change = 0.0712, Max-Change = 0.1591, Max-Change = 0.0433, Max-Change = 0.0432, Max-Change = 0.0488, Max-Change = 0.1257, Max-Change = 0.2000, Max-Change = 0.1857, Max-Change = 0.0790, Max-Change = 0.0391, Max-Change = 0.0559, Max-Change = 0.2000, Max-Change = 0.1776, Max-Change = 0.0181, Max-Change = 0.0494, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1295, Max-Change = 0.1152, Max-Change = 0.2000, Max-Change = 0.1530, Max-Change = 0.0186, Max-Change = 0.1586, Max-Change = 0.0536, Max-Change = 0.0350, Max-Change = 0.0466, Max-Change = 0.0195, Max-Change = 0.0606, Max-Change = 0.1190, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.1087, gam = 0.1057, Max-Change = 0.0780, gam = 0.0780, Max-Change = 0.0283, gam = 0.0629, Max-Change = 0.0549, gam = 0.0532, Max-Change = 0.0181, gam = 0.0464, Max-Change = 0.0171, gam = 0.0413, Max-Change = 0.0038, gam = 0.0374, Max-Change = 0.0220, gam = 0.0342, Max-Change = 0.0112, gam = 0.0316, Max-Change = 0.0208, gam = 0.0294, Max-Change = 0.0122, gam = 0.0276, Max-Change = 0.0115, gam = 0.0260, Max-Change = 0.0075, gam = 0.0246, Max-Change = 0.0027, gam = 0.0233, Max-Change = 0.0056, gam = 0.0222, Max-Change = 0.0078, gam = 0.0212, Max-Change = 0.0072, gam = 0.0203, Max-Change = 0.0164, gam = 0.0195, Max-Change = 0.0044, gam = 0.0188, Max-Change = 0.0097, gam = 0.0181, Max-Change = 0.0156, gam = 0.0175, Max-Change = 0.0020, gam = 0.0169, Max-Change = 0.0043, gam = 0.0164, Max-Change = 0.0047, gam = 0.0159, Max-Change = 0.0035, gam = 0.0154, Max-Change = 0.0101, gam = 0.0150, Max-Change = 0.0149, gam = 0.0146, Max-Change = 0.0065, gam = 0.0142, Max-Change = 0.0078, gam = 0.0139, Max-Change = 0.0128, gam = 0.0135, Max-Change = 0.0101, gam = 0.0132, Max-Change = 0.0038, gam = 0.0129, Max-Change = 0.0029, gam = 0.0126, Max-Change = 0.0119, gam = 0.0124, Max-Change = 0.0033, gam = 0.0121, Max-Change = 0.0080, gam = 0.0119, Max-Change = 0.0169, gam = 0.0116, Max-Change = 0.0034, gam = 0.0114, Max-Change = 0.0056, gam = 0.0112, Max-Change = 0.0060, gam = 0.0110, Max-Change = 0.0023, gam = 0.0108, Max-Change = 0.0043, gam = 0.0106, Max-Change = 0.0020, gam = 0.0104, Max-Change = 0.0057, gam = 0.0102, Max-Change = 0.0074, gam = 0.0101, Max-Change = 0.0063, gam = 0.0099, Max-Change = 0.0016, gam = 0.0098, Max-Change = 0.0024, gam = 0.0096, Max-Change = 0.0020, gam = 0.0095, Max-Change = 0.0067, gam = 0.0093, Max-Change = 0.0017, gam = 0.0092, Max-Change = 0.0083, gam = 0.0091, Max-Change = 0.0015, gam = 0.0089, Max-Change = 0.0056, gam = 0.0088, Max-Change = 0.0027, gam = 0.0087, Max-Change = 0.0033, gam = 0.0086, Max-Change = 0.0011, gam = 0.0085, Max-Change = 0.0055, gam = 0.0084, Max-Change = 0.0023, gam = 0.0082, Max-Change = 0.0023, gam = 0.0081, Max-Change = 0.0026, gam = 0.0080, Max-Change = 0.0037, gam = 0.0080, Max-Change = 0.0017, gam = 0.0079, Max-Change = 0.0048, gam = 0.0078, Max-Change = 0.0043, gam = 0.0077, Max-Change = 0.0037, gam = 0.0076, Max-Change = 0.0094, gam = 0.0075, Max-Change = 0.0067, gam = 0.0074, Max-Change = 0.0068, gam = 0.0073, Max-Change = 0.0047, gam = 0.0073, Max-Change = 0.0012, gam = 0.0072, Max-Change = 0.0018, gam = 0.0071, Max-Change = 0.0062, gam = 0.0070, Max-Change = 0.0020, gam = 0.0070, Max-Change = 0.0033, gam = 0.0069, Max-Change = 0.0012, gam = 0.0068, Max-Change = 0.0033, gam = 0.0068, Max-Change = 0.0062, gam = 0.0067, Max-Change = 0.0011, gam = 0.0066, Max-Change = 0.0029, gam = 0.0066, Max-Change = 0.0041, gam = 0.0065, Max-Change = 0.0012, gam = 0.0065, Max-Change = 0.0057, gam = 0.0064, Max-Change = 0.0021, gam = 0.0064, Max-Change = 0.0071, gam = 0.0063, Max-Change = 0.0066, gam = 0.0062, Max-Change = 0.0013, gam = 0.0062, Max-Change = 0.0006, gam = 0.0061, Max-Change = 0.0019, gam = 0.0061, Max-Change = 0.0009, gam = 0.0060, Max-Change = 0.0007, gam = 0.0060, Max-Change = 0.0083, gam = 0.0059, Max-Change = 0.0034, gam = 0.0059, Max-Change = 0.0028, gam = 0.0058, Max-Change = 0.0063, gam = 0.0058, Max-Change = 0.0073, gam = 0.0058, Max-Change = 0.0037, gam = 0.0057, Max-Change = 0.0018, gam = 0.0057, Max-Change = 0.0017, gam = 0.0056, Max-Change = 0.0037, gam = 0.0056, Max-Change = 0.0027, gam = 0.0055, Max-Change = 0.0012, gam = 0.0055, Max-Change = 0.0006, gam = 0.0055, Max-Change = 0.0031, gam = 0.0054, Max-Change = 0.0048, gam = 0.0054, Max-Change = 0.0019, gam = 0.0053, Max-Change = 0.0032, gam = 0.0053, Max-Change = 0.0005, gam = 0.0053, Max-Change = 0.0049, gam = 0.0052, Max-Change = 0.0029, gam = 0.0052, Max-Change = 0.0023, gam = 0.0052, Max-Change = 0.0056, gam = 0.0051, Max-Change = 0.0060, gam = 0.0051, Max-Change = 0.0031, gam = 0.0051, Max-Change = 0.0050, gam = 0.0050, Max-Change = 0.0049, gam = 0.0050, Max-Change = 0.0009, gam = 0.0050, Max-Change = 0.0036, gam = 0.0049, Max-Change = 0.0055, gam = 0.0049, Max-Change = 0.0032, gam = 0.0049, Max-Change = 0.0043, gam = 0.0048, Max-Change = 0.0007, gam = 0.0048, Max-Change = 0.0024, gam = 0.0048, Max-Change = 0.0019, gam = 0.0048, Max-Change = 0.0025, gam = 0.0047, Max-Change = 0.0025, gam = 0.0047, Max-Change = 0.0007, gam = 0.0047, Max-Change = 0.0005, gam = 0.0046, Max-Change = 0.0008
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
anova(mod1, mod1b) #better with 2PL models using all criteria (as expected, given simdata pars)
#>            AIC    SABIC       HQ      BIC    logLik     X2 df p
#> mod1  8111.326 8130.107 8134.469 8171.387 -4042.663            
#> mod1b 7975.487 8007.270 8014.651 8077.128 -3965.743 153.84  9 0

# continuous predictor with group
mod2 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items + pseudoIQ)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1807, Max-Change = 0.1593, Max-Change = 0.1350, Max-Change = 0.1155, Max-Change = 0.1001, Max-Change = 0.0861, Max-Change = 0.0702, Max-Change = 0.0585, Max-Change = 0.0465, Max-Change = 0.0362, Max-Change = 0.0292, Max-Change = 0.0245, Max-Change = 0.0166, Max-Change = 0.0178, Max-Change = 0.0144, Max-Change = 0.0085, Max-Change = 0.0130, Max-Change = 0.0096, Max-Change = 0.0103, Max-Change = 0.0061, Max-Change = 0.0049, Max-Change = 0.0084, Max-Change = 0.0113, Max-Change = 0.0152, Max-Change = 0.0156, Max-Change = 0.0068, Max-Change = 0.0083, Max-Change = 0.0068, Max-Change = 0.0079, Max-Change = 0.0075, Max-Change = 0.0149, Max-Change = 0.0082, Max-Change = 0.0035, Max-Change = 0.0157, Max-Change = 0.0115, Max-Change = 0.0088, Max-Change = 0.0055, Max-Change = 0.0128, Max-Change = 0.0083, Max-Change = 0.0046, Max-Change = 0.0066, Max-Change = 0.0040, Max-Change = 0.0086, Max-Change = 0.0020, Max-Change = 0.0069, Max-Change = 0.0051, Max-Change = 0.0081, Max-Change = 0.0034, Max-Change = 0.0156, Max-Change = 0.0074, Max-Change = 0.0066, Max-Change = 0.0084, Max-Change = 0.0088, Max-Change = 0.0151, Max-Change = 0.0069, Max-Change = 0.0178, Max-Change = 0.0042, Max-Change = 0.0129, Max-Change = 0.0094, Max-Change = 0.0068, Max-Change = 0.0076, Max-Change = 0.0032, Max-Change = 0.0093, Max-Change = 0.0067, Max-Change = 0.0026, Max-Change = 0.0060, Max-Change = 0.0056, Max-Change = 0.0064, Max-Change = 0.0079, Max-Change = 0.0038, Max-Change = 0.0052, Max-Change = 0.0105, Max-Change = 0.0046, Max-Change = 0.0047, Max-Change = 0.0059, Max-Change = 0.0055, Max-Change = 0.0067, Max-Change = 0.0072, Max-Change = 0.0068, Max-Change = 0.0074, Max-Change = 0.0033, Max-Change = 0.0105, Max-Change = 0.0133, Max-Change = 0.0075, Max-Change = 0.0087, Max-Change = 0.0052, Max-Change = 0.0040, Max-Change = 0.0056, Max-Change = 0.0094, Max-Change = 0.0095, Max-Change = 0.0116, Max-Change = 0.0043, Max-Change = 0.0044, Max-Change = 0.0038, Max-Change = 0.0096, Max-Change = 0.0055, Max-Change = 0.0047, Max-Change = 0.0035, Max-Change = 0.0089, Max-Change = 0.0055, Max-Change = 0.0059, Max-Change = 0.0088, Max-Change = 0.0103, Max-Change = 0.0040, Max-Change = 0.0033, Max-Change = 0.0102, Max-Change = 0.0063, Max-Change = 0.0033, Max-Change = 0.0072, Max-Change = 0.0095, Max-Change = 0.0037, Max-Change = 0.0092, Max-Change = 0.0103, Max-Change = 0.0058, Max-Change = 0.0068, Max-Change = 0.0071, Max-Change = 0.0047, Max-Change = 0.0138, Max-Change = 0.0077, Max-Change = 0.0070, Max-Change = 0.0086, Max-Change = 0.0091, Max-Change = 0.0099, Max-Change = 0.0031, Max-Change = 0.0124, Max-Change = 0.0060, Max-Change = 0.0074, Max-Change = 0.0038, Max-Change = 0.0101, Max-Change = 0.0044, Max-Change = 0.0073, Max-Change = 0.0048, Max-Change = 0.0057, Max-Change = 0.0058, Max-Change = 0.0079, Max-Change = 0.0054, Max-Change = 0.0106, Max-Change = 0.0074, Max-Change = 0.0023, Max-Change = 0.0120, Max-Change = 0.0061, Max-Change = 0.0063, Max-Change = 0.0069, Max-Change = 0.0075, Max-Change = 0.0097, Max-Change = 0.0032, Max-Change = 0.0077, Max-Change = 0.0058, Max-Change = 0.0041, Max-Change = 0.0022, Max-Change = 0.0034, Max-Change = 0.0067, Max-Change = 0.0071, Max-Change = 0.0055, Max-Change = 0.0043, Max-Change = 0.0030, Max-Change = 0.0056, Max-Change = 0.0043, Max-Change = 0.0119, Max-Change = 0.0053, Max-Change = 0.0059, Max-Change = 0.0087, Max-Change = 0.0026, Max-Change = 0.0041, Max-Change = 0.0052, Max-Change = 0.0077, Max-Change = 0.0074, Max-Change = 0.0031, Max-Change = 0.0066, Max-Change = 0.0108, Max-Change = 0.0058, Max-Change = 0.0067, Max-Change = 0.0018, Max-Change = 0.0069, Max-Change = 0.0096, Max-Change = 0.0059, Max-Change = 0.0066, Max-Change = 0.0018, Max-Change = 0.0049, Max-Change = 0.0040, Max-Change = 0.0106, Max-Change = 0.0046, Max-Change = 0.0021, Max-Change = 0.0047, Max-Change = 0.0034, Max-Change = 0.0067, Max-Change = 0.0070, Max-Change = 0.0061, Max-Change = 0.0059, Max-Change = 0.0024, Max-Change = 0.0069, Max-Change = 0.0041, Max-Change = 0.0036, Max-Change = 0.0061, Max-Change = 0.0058, Max-Change = 0.0057, Max-Change = 0.0072, Max-Change = 0.0068, Max-Change = 0.0071, Max-Change = 0.0056, Max-Change = 0.0109, Max-Change = 0.0047, Max-Change = 0.0063, Max-Change = 0.0082, Max-Change = 0.0020, Max-Change = 0.0107, Max-Change = 0.0032, Max-Change = 0.0050, Max-Change = 0.0109, Max-Change = 0.0056, Max-Change = 0.0076, Max-Change = 0.0024, Max-Change = 0.0072, Max-Change = 0.0038, Max-Change = 0.0062, Max-Change = 0.0041, Max-Change = 0.0049, Max-Change = 0.0037, Max-Change = 0.0042, Max-Change = 0.0044, Max-Change = 0.0017, Max-Change = 0.0054, Max-Change = 0.0066, Max-Change = 0.0098, Max-Change = 0.0067, Max-Change = 0.0066, Max-Change = 0.0087, Max-Change = 0.0095, Max-Change = 0.0073, Max-Change = 0.0092, Max-Change = 0.0051, Max-Change = 0.0043, Max-Change = 0.0059, Max-Change = 0.0057, Max-Change = 0.0041, Max-Change = 0.0052, Max-Change = 0.0076, Max-Change = 0.0058, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0043, gam = 0.1057, Max-Change = 0.0022, gam = 0.0780, Max-Change = 0.0015, gam = 0.0629, Max-Change = 0.0023, gam = 0.0532, Max-Change = 0.0019, gam = 0.0464, Max-Change = 0.0010, gam = 0.0413, Max-Change = 0.0004, gam = 0.0374, Max-Change = 0.0012, gam = 0.0342, Max-Change = 0.0014, gam = 0.0316, Max-Change = 0.0007, gam = 0.0294, Max-Change = 0.0007, gam = 0.0276, Max-Change = 0.0006
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod2)
#> 
#> Call:
#> mixedmirt(data = data, covdata = covdata, model = model, fixed = ~0 + 
#>     group + items + pseudoIQ)
#> 
#> --------------
#> FIXED EFFECTS:
#>          Estimate Std.Error z.value
#> groupG1    -1.711     0.105 -16.317
#> groupG2    -0.750     0.093  -8.029
#> groupG3     0.366     0.099   3.714
#> pseudoIQ    0.268     0.058   4.594
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>       Theta
#> Theta 0.105
#> 
anova(mod1b, mod2)
#>            AIC    SABIC       HQ      BIC    logLik       X2 df   p
#> mod1b 7975.487 8007.270 8014.651 8077.128 -3965.743                
#> mod2  8090.269 8110.494 8115.192 8154.950 -4031.134 -130.782 -8 NaN

# view fixed design matrix with and without unique item level intercepts
withint <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, return.design = TRUE)
withoutint <- mixedmirt(data, covdata, model, fixed = ~ 0 + group, return.design = TRUE)

# notice that in result above, the intercepts 'items1 to items 10' were reassigned to 'd'
head(withint$X)
#>     items1 items2 items3 items4 items5 items6 items7 items8 items9 items10
#> 1.1      1      0      0      0      0      0      0      0      0       0
#> 2.1      1      0      0      0      0      0      0      0      0       0
#> 3.1      1      0      0      0      0      0      0      0      0       0
#> 4.1      1      0      0      0      0      0      0      0      0       0
#> 5.1      1      0      0      0      0      0      0      0      0       0
#> 6.1      1      0      0      0      0      0      0      0      0       0
#>     groupG2 groupG3
#> 1.1       0       0
#> 2.1       0       0
#> 3.1       0       0
#> 4.1       0       0
#> 5.1       0       0
#> 6.1       0       0
tail(withint$X)
#>        items1 items2 items3 items4 items5 items6 items7 items8 items9 items10
#> 745.10      0      0      0      0      0      0      0      0      0       1
#> 746.10      0      0      0      0      0      0      0      0      0       1
#> 747.10      0      0      0      0      0      0      0      0      0       1
#> 748.10      0      0      0      0      0      0      0      0      0       1
#> 749.10      0      0      0      0      0      0      0      0      0       1
#> 750.10      0      0      0      0      0      0      0      0      0       1
#>        groupG2 groupG3
#> 745.10       0       1
#> 746.10       0       1
#> 747.10       0       1
#> 748.10       0       1
#> 749.10       0       1
#> 750.10       0       1
head(withoutint$X) # no intercepts design here to be reassigned into item intercepts
#>     groupG1 groupG2 groupG3
#> 1.1       1       0       0
#> 2.1       1       0       0
#> 3.1       1       0       0
#> 4.1       1       0       0
#> 5.1       1       0       0
#> 6.1       1       0       0
tail(withoutint$X)
#>        groupG1 groupG2 groupG3
#> 745.10       0       0       1
#> 746.10       0       0       1
#> 747.10       0       0       1
#> 748.10       0       0       1
#> 749.10       0       0       1
#> 750.10       0       0       1

###################################################
### random effects
# make the number of groups much larger
covdata$group <- factor(rep(paste0('G',1:50), each = N/50))

# random groups
rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1864, Max-Change = 0.1623, Max-Change = 0.1436, Max-Change = 0.1243, Max-Change = 0.1035, Max-Change = 0.0953, Max-Change = 0.0862, Max-Change = 0.0730, Max-Change = 0.0612, Max-Change = 0.0490, Max-Change = 0.0426, Max-Change = 0.0342, Max-Change = 0.0359, Max-Change = 0.0255, Max-Change = 0.0251, Max-Change = 0.0175, Max-Change = 0.0105, Max-Change = 0.0121, Max-Change = 0.0184, Max-Change = 0.0194, Max-Change = 0.0092, Max-Change = 0.0064, Max-Change = 0.0083, Max-Change = 0.0065, Max-Change = 0.0026, Max-Change = 0.0066, Max-Change = 0.0110, Max-Change = 0.0004, Max-Change = 0.0044, Max-Change = 0.0104, Max-Change = 0.0118, Max-Change = 0.0071, Max-Change = 0.0056, Max-Change = 0.0073, Max-Change = 0.0134, Max-Change = 0.0043, Max-Change = 0.0143, Max-Change = 0.0042, Max-Change = 0.0109, Max-Change = 0.0061, Max-Change = 0.0096, Max-Change = 0.0130, Max-Change = 0.0158, Max-Change = 0.0124, Max-Change = 0.0075, Max-Change = 0.0080, Max-Change = 0.0056, Max-Change = 0.0025, Max-Change = 0.0109, Max-Change = 0.0030, Max-Change = 0.0135, Max-Change = 0.0156, Max-Change = 0.0136, Max-Change = 0.0085, Max-Change = 0.0083, Max-Change = 0.0118, Max-Change = 0.0017, Max-Change = 0.0074, Max-Change = 0.0278, Max-Change = 0.0030, Max-Change = 0.0086, Max-Change = 0.0053, Max-Change = 0.0006, Max-Change = 0.0139, Max-Change = 0.0106, Max-Change = 0.0060, Max-Change = 0.0048, Max-Change = 0.0084, Max-Change = 0.0076, Max-Change = 0.0094, Max-Change = 0.0090, Max-Change = 0.0120, Max-Change = 0.0036, Max-Change = 0.0029, Max-Change = 0.0089, Max-Change = 0.0036, Max-Change = 0.0030, Max-Change = 0.0074, Max-Change = 0.0044, Max-Change = 0.0153, Max-Change = 0.0043, Max-Change = 0.0295, Max-Change = 0.0054, Max-Change = 0.0105, Max-Change = 0.0154, Max-Change = 0.0094, Max-Change = 0.0021, Max-Change = 0.0077, Max-Change = 0.0113, Max-Change = 0.0035, Max-Change = 0.0045, Max-Change = 0.0065, Max-Change = 0.2000, Max-Change = 0.1344, Max-Change = 0.0744, Max-Change = 0.0786, Max-Change = 0.0581, Max-Change = 0.0426, Max-Change = 0.0315, Max-Change = 0.0381, Max-Change = 0.0227, Max-Change = 0.0228, Max-Change = 0.0159, Max-Change = 0.0193, Max-Change = 0.0112, Max-Change = 0.0092, Max-Change = 0.0150, Max-Change = 0.0126, Max-Change = 0.0093, Max-Change = 0.0185, Max-Change = 0.0088, Max-Change = 0.0108, Max-Change = 0.0167, Max-Change = 0.0128, Max-Change = 0.0118, Max-Change = 0.0094, Max-Change = 0.0056, Max-Change = 0.0258, Max-Change = 0.0316, Max-Change = 0.0118, Max-Change = 0.0106, Max-Change = 0.0198, Max-Change = 0.0252, Max-Change = 0.0309, Max-Change = 0.0159, Max-Change = 0.0162, Max-Change = 0.0087, Max-Change = 0.0045, Max-Change = 0.0110, Max-Change = 0.0057, Max-Change = 0.0069, Max-Change = 0.0039, Max-Change = 0.0225, Max-Change = 0.0115, Max-Change = 0.0144, Max-Change = 0.0102, Max-Change = 0.0078, Max-Change = 0.0102, Max-Change = 0.0065, Max-Change = 0.0275, Max-Change = 0.0248, Max-Change = 0.0045, Max-Change = 0.0114, Max-Change = 0.0065, Max-Change = 0.0160, Max-Change = 0.0145, Max-Change = 0.0068, Max-Change = 0.0129, Max-Change = 0.0116, Max-Change = 0.0161, Max-Change = 0.0138, Max-Change = 0.0076, Max-Change = 0.0139, Max-Change = 0.0116, Max-Change = 0.0067, Max-Change = 0.0162, Max-Change = 0.0228, Max-Change = 0.0090, Max-Change = 0.0079, Max-Change = 0.0054, Max-Change = 0.0223, Max-Change = 0.0171, Max-Change = 0.0075, Max-Change = 0.0105, Max-Change = 0.0246, Max-Change = 0.0028, Max-Change = 0.0061, Max-Change = 0.0095, Max-Change = 0.0026, Max-Change = 0.0123, Max-Change = 0.0207, Max-Change = 0.0063, Max-Change = 0.0044, Max-Change = 0.0084, Max-Change = 0.0051, Max-Change = 0.0050, Max-Change = 0.0085, Max-Change = 0.0051, Max-Change = 0.0207, Max-Change = 0.0080, Max-Change = 0.0121, Max-Change = 0.0258, Max-Change = 0.0036, Max-Change = 0.0144, Max-Change = 0.0108, Max-Change = 0.0057, Max-Change = 0.0055, Max-Change = 0.0196, Max-Change = 0.0071, Max-Change = 0.0062, Max-Change = 0.0077, Max-Change = 0.0068, Max-Change = 0.0187, Max-Change = 0.0174, Max-Change = 0.0161, Max-Change = 0.0129, Max-Change = 0.0047, Max-Change = 0.0088, Max-Change = 0.0096, Max-Change = 0.0171, Max-Change = 0.0073, Max-Change = 0.0098, Max-Change = 0.0079, Max-Change = 0.0198, Max-Change = 0.0111, Max-Change = 0.0425, Max-Change = 0.0146, Max-Change = 0.0073, Max-Change = 0.0110, Max-Change = 0.0258, Max-Change = 0.0168, Max-Change = 0.0345, Max-Change = 0.0010, Max-Change = 0.0187, Max-Change = 0.0113, Max-Change = 0.0134, Max-Change = 0.0198, Max-Change = 0.0027, Max-Change = 0.0127, Max-Change = 0.0065, Max-Change = 0.0132, Max-Change = 0.0139, Max-Change = 0.0045, Max-Change = 0.0400, Max-Change = 0.0098, Max-Change = 0.0053, Max-Change = 0.0015, Max-Change = 0.0103, Max-Change = 0.0102, Max-Change = 0.0054, Max-Change = 0.0199, Max-Change = 0.0100, Max-Change = 0.0041, Max-Change = 0.0184, Max-Change = 0.0145, Max-Change = 0.0103, Max-Change = 0.0038, Max-Change = 0.0143, Max-Change = 0.0152, Max-Change = 0.0110, Max-Change = 0.0068, Max-Change = 0.0118, Max-Change = 0.0066, Max-Change = 0.0130, Max-Change = 0.0208, Max-Change = 0.0037, Max-Change = 0.0121, Max-Change = 0.0141, Max-Change = 0.0072, Max-Change = 0.0143, Max-Change = 0.0123, Max-Change = 0.0070, Max-Change = 0.0035, Max-Change = 0.0055, Max-Change = 0.0196, Max-Change = 0.0027, Max-Change = 0.0068, Max-Change = 0.0078, Max-Change = 0.0070, Max-Change = 0.0023, Max-Change = 0.0069, Max-Change = 0.0069, Max-Change = 0.0115, Max-Change = 0.0052, Max-Change = 0.0048, Max-Change = 0.0232, Max-Change = 0.0055, Max-Change = 0.0128, Max-Change = 0.0051, Max-Change = 0.0093, Max-Change = 0.0252, Max-Change = 0.0142, Max-Change = 0.0101, Max-Change = 0.0092, Max-Change = 0.0333, Max-Change = 0.0122, Max-Change = 0.0084, Max-Change = 0.0098, Max-Change = 0.0081, Max-Change = 0.0300, Max-Change = 0.0306, Max-Change = 0.0112, Max-Change = 0.0048, Max-Change = 0.0091, Max-Change = 0.0088, Max-Change = 0.0263, Max-Change = 0.0100, Max-Change = 0.0118, Max-Change = 0.0103, Max-Change = 0.0021, Max-Change = 0.0084, Max-Change = 0.0130, Max-Change = 0.0080, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0151, gam = 0.1057, Max-Change = 0.0072, gam = 0.0780, Max-Change = 0.0048, gam = 0.0629, Max-Change = 0.0014, gam = 0.0532, Max-Change = 0.0047, gam = 0.0464, Max-Change = 0.0016, gam = 0.0413, Max-Change = 0.0021, gam = 0.0374, Max-Change = 0.0007, gam = 0.0342, Max-Change = 0.0026, gam = 0.0316, Max-Change = 0.0035, gam = 0.0294, Max-Change = 0.0018, gam = 0.0276, Max-Change = 0.0009, gam = 0.0260, Max-Change = 0.0013, gam = 0.0246, Max-Change = 0.0007, gam = 0.0233, Max-Change = 0.0016, gam = 0.0222, Max-Change = 0.0012, gam = 0.0212, Max-Change = 0.0011, gam = 0.0203, Max-Change = 0.0024, gam = 0.0195, Max-Change = 0.0020, gam = 0.0188, Max-Change = 0.0010, gam = 0.0181, Max-Change = 0.0013, gam = 0.0175, Max-Change = 0.0003, gam = 0.0169, Max-Change = 0.0011, gam = 0.0164, Max-Change = 0.0008, gam = 0.0159, Max-Change = 0.0007, gam = 0.0154, Max-Change = 0.0008
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(rmod1)
#> 
#> Call:
#> mixedmirt(data = data, covdata = covdata, model = 1, fixed = ~0 + 
#>     items, random = ~1 | group)
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>       F1
#> F1 0.062
#> 
#> $group
#>           COV_group
#> COV_group      1.11
#> 
coef(rmod1)
#> $Item_1
#>         a1      d  g  u
#> par      1 -0.713  0  1
#> CI_2.5  NA -0.927 NA NA
#> CI_97.5 NA -0.500 NA NA
#> 
#> $Item_2
#>         a1      d  g  u
#> par      1 -0.867  0  1
#> CI_2.5  NA -1.081 NA NA
#> CI_97.5 NA -0.653 NA NA
#> 
#> $Item_3
#>         a1      d  g  u
#> par      1 -0.370  0  1
#> CI_2.5  NA -0.583 NA NA
#> CI_97.5 NA -0.156 NA NA
#> 
#> $Item_4
#>         a1      d  g  u
#> par      1  0.057  0  1
#> CI_2.5  NA -0.161 NA NA
#> CI_97.5 NA  0.276 NA NA
#> 
#> $Item_5
#>         a1     d  g  u
#> par      1 0.654  0  1
#> CI_2.5  NA 0.422 NA NA
#> CI_97.5 NA 0.886 NA NA
#> 
#> $Item_6
#>         a1      d  g  u
#> par      1 -0.145  0  1
#> CI_2.5  NA -0.361 NA NA
#> CI_97.5 NA  0.071 NA NA
#> 
#> $Item_7
#>         a1      d  g  u
#> par      1 -0.592  0  1
#> CI_2.5  NA -0.805 NA NA
#> CI_97.5 NA -0.379 NA NA
#> 
#> $Item_8
#>         a1      d  g  u
#> par      1 -1.057  0  1
#> CI_2.5  NA -1.273 NA NA
#> CI_97.5 NA -0.842 NA NA
#> 
#> $Item_9
#>         a1      d  g  u
#> par      1 -0.957  0  1
#> CI_2.5  NA -1.172 NA NA
#> CI_97.5 NA -0.743 NA NA
#> 
#> $Item_10
#>         a1     d  g  u
#> par      1 2.769  0  1
#> CI_2.5  NA 2.425 NA NA
#> CI_97.5 NA 3.113 NA NA
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0  0.062
#> CI_2.5      NA  0.004
#> CI_97.5     NA  0.120
#> 
#> $group
#>         COV_group_group
#> par               1.112
#> CI_2.5            0.631
#> CI_97.5           1.592
#> 

# random groups and random items
rmod2 <- mixedmirt(data, covdata, 1, random = list(~ 1|group, ~ 1|items))
#> , Max-Change = 0.0274, Max-Change = 0.0217, Max-Change = 0.0087, Max-Change = 0.0228, Max-Change = 0.0182, Max-Change = 0.0177, Max-Change = 0.0220, Max-Change = 0.0255, Max-Change = 0.0132, Max-Change = 0.0106, Max-Change = 0.0070, Max-Change = 0.0106, Max-Change = 0.0120, Max-Change = 0.0060, Max-Change = 0.0116, Max-Change = 0.0019, Max-Change = 0.0029, Max-Change = 0.0069, Max-Change = 0.0130, Max-Change = 0.0098, Max-Change = 0.0232, Max-Change = 0.0104, Max-Change = 0.0034, Max-Change = 0.0054, Max-Change = 0.0075, Max-Change = 0.0051, Max-Change = 0.0113, Max-Change = 0.0033, Max-Change = 0.0032, Max-Change = 0.0056, Max-Change = 0.0037, Max-Change = 0.0047, Max-Change = 0.0051, Max-Change = 0.0028, Max-Change = 0.0024, Max-Change = 0.0006, Max-Change = 0.0045, Max-Change = 0.0037, Max-Change = 0.0007, Max-Change = 0.0013, Max-Change = 0.0057, Max-Change = 0.0169, Max-Change = 0.0131, Max-Change = 0.0095, Max-Change = 0.0038, Max-Change = 0.0032, Max-Change = 0.0049, Max-Change = 0.0047, Max-Change = 0.0137, Max-Change = 0.0022, Max-Change = 0.0077, Max-Change = 0.0074, Max-Change = 0.0048, Max-Change = 0.0044, Max-Change = 0.0044, Max-Change = 0.0059, Max-Change = 0.0017, Max-Change = 0.0037, Max-Change = 0.0062, Max-Change = 0.0062, Max-Change = 0.0026, Max-Change = 0.0063, Max-Change = 0.0057, Max-Change = 0.0151, Max-Change = 0.0044, Max-Change = 0.0059, Max-Change = 0.0086, Max-Change = 0.0036, Max-Change = 0.0044, Max-Change = 0.0054, Max-Change = 0.0030, Max-Change = 0.0035, Max-Change = 0.0056, Max-Change = 0.0004, Max-Change = 0.0009, Max-Change = 0.0027, Max-Change = 0.0004, Max-Change = 0.0020, Max-Change = 0.0066, Max-Change = 0.0013, Max-Change = 0.0023, Max-Change = 0.0045, Max-Change = 0.0076, Max-Change = 0.0024, Max-Change = 0.0015, Max-Change = 0.0012, Max-Change = 0.0033, Max-Change = 0.0040, Max-Change = 0.0105, Max-Change = 0.0004, Max-Change = 0.0064, Max-Change = 0.0121, Max-Change = 0.0016, Max-Change = 0.0090, Max-Change = 0.0018, Max-Change = 0.0069, Max-Change = 0.0070, Max-Change = 0.0055, Max-Change = 0.0037, Max-Change = 0.0504, Max-Change = 0.1569, Max-Change = 0.0713, Max-Change = 0.0281, Max-Change = 0.0327, Max-Change = 0.0567, Max-Change = 0.0160, Max-Change = 0.0127, Max-Change = 0.0138, Max-Change = 0.0226, Max-Change = 0.0120, Max-Change = 0.0178, Max-Change = 0.0122, Max-Change = 0.0057, Max-Change = 0.0087, Max-Change = 0.0137, Max-Change = 0.0088, Max-Change = 0.0196, Max-Change = 0.0088, Max-Change = 0.0167, Max-Change = 0.0123, Max-Change = 0.0119, Max-Change = 0.0139, Max-Change = 0.0092, Max-Change = 0.0199, Max-Change = 0.0132, Max-Change = 0.0030, Max-Change = 0.0045, Max-Change = 0.0065, Max-Change = 0.0194, Max-Change = 0.0114, Max-Change = 0.0146, Max-Change = 0.0089, Max-Change = 0.0129, Max-Change = 0.0156, Max-Change = 0.0131, Max-Change = 0.0130, Max-Change = 0.0087, Max-Change = 0.0086, Max-Change = 0.0054, Max-Change = 0.0052, Max-Change = 0.0106, Max-Change = 0.0146, Max-Change = 0.0035, Max-Change = 0.0120, Max-Change = 0.0061, Max-Change = 0.0061, Max-Change = 0.0023, Max-Change = 0.0062, Max-Change = 0.0086, Max-Change = 0.0083, Max-Change = 0.0035, Max-Change = 0.0077, Max-Change = 0.0064, Max-Change = 0.0033, Max-Change = 0.0061, Max-Change = 0.0065, Max-Change = 0.0129, Max-Change = 0.0230, Max-Change = 0.0068, Max-Change = 0.0033, Max-Change = 0.0121, Max-Change = 0.0015, Max-Change = 0.0018, Max-Change = 0.0051, Max-Change = 0.0016, Max-Change = 0.0140, Max-Change = 0.0063, Max-Change = 0.0193, Max-Change = 0.0179, Max-Change = 0.0064, Max-Change = 0.0051, Max-Change = 0.0061, Max-Change = 0.0067, Max-Change = 0.0110, Max-Change = 0.0131, Max-Change = 0.0151, Max-Change = 0.0057, Max-Change = 0.0098, Max-Change = 0.0081, Max-Change = 0.0068, Max-Change = 0.0101, Max-Change = 0.0025, Max-Change = 0.0068, Max-Change = 0.0015, Max-Change = 0.0160, Max-Change = 0.0077, Max-Change = 0.0060, Max-Change = 0.0057, Max-Change = 0.0068, Max-Change = 0.0053, Max-Change = 0.0143, Max-Change = 0.0049, Max-Change = 0.0012, Max-Change = 0.0075, Max-Change = 0.0059, Max-Change = 0.0126, Max-Change = 0.0084, Max-Change = 0.0163, Max-Change = 0.0068, Max-Change = 0.0511, Max-Change = 0.0088, Max-Change = 0.0075, Max-Change = 0.0063, Max-Change = 0.0070, Max-Change = 0.0058, Max-Change = 0.0064, Max-Change = 0.0086, Max-Change = 0.0082, Max-Change = 0.0145, Max-Change = 0.0049, Max-Change = 0.0087, Max-Change = 0.0072, Max-Change = 0.0055, Max-Change = 0.0123, Max-Change = 0.0124, Max-Change = 0.0067, Max-Change = 0.0308, Max-Change = 0.0055, Max-Change = 0.0237, Max-Change = 0.0222, Max-Change = 0.0140, Max-Change = 0.0075, Max-Change = 0.0180, Max-Change = 0.0111, Max-Change = 0.0081, Max-Change = 0.0086, Max-Change = 0.0068, Max-Change = 0.0053, Max-Change = 0.0135, Max-Change = 0.0094, Max-Change = 0.0109, Max-Change = 0.0093, Max-Change = 0.0060, Max-Change = 0.0103, Max-Change = 0.0121, Max-Change = 0.0125, Max-Change = 0.0084, Max-Change = 0.0025, Max-Change = 0.0033, Max-Change = 0.0045, Max-Change = 0.0197, Max-Change = 0.0107, Max-Change = 0.0052, Max-Change = 0.0106, Max-Change = 0.0027, Max-Change = 0.0060, Max-Change = 0.0156, Max-Change = 0.0083, Max-Change = 0.0089, Max-Change = 0.0035, Max-Change = 0.0071, Max-Change = 0.0058, Max-Change = 0.0088, Max-Change = 0.0106, Max-Change = 0.0061, Max-Change = 0.0024, Max-Change = 0.0024, Max-Change = 0.0082, Max-Change = 0.0089, Max-Change = 0.0078, Max-Change = 0.0087, Max-Change = 0.0008, Max-Change = 0.0146, Max-Change = 0.0082, Max-Change = 0.0179, Max-Change = 0.0087, Max-Change = 0.0155, Max-Change = 0.0085, Max-Change = 0.0156, Max-Change = 0.0047, Max-Change = 0.0068, Max-Change = 0.0045, Max-Change = 0.0033, Max-Change = 0.0102, Max-Change = 0.0098, Max-Change = 0.0087, Max-Change = 0.0054, Max-Change = 0.0097, Max-Change = 0.0142, Max-Change = 0.0110, Max-Change = 0.0107, Max-Change = 0.0042, Max-Change = 0.0030, Max-Change = 0.0054, Max-Change = 0.0063, Max-Change = 0.0173, Max-Change = 0.0087, Max-Change = 0.0167, Max-Change = 0.0098, Max-Change = 0.0046, Max-Change = 0.0170, Max-Change = 0.0173, Max-Change = 0.0085, Max-Change = 0.0076, Max-Change = 0.0287, Max-Change = 0.0082, Max-Change = 0.0072, Max-Change = 0.0086, Max-Change = 0.0074, Max-Change = 0.0086, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0048, gam = 0.1057, Max-Change = 0.0035, gam = 0.0780, Max-Change = 0.0012, gam = 0.0629, Max-Change = 0.0029, gam = 0.0532, Max-Change = 0.0044, gam = 0.0464, Max-Change = 0.0044, gam = 0.0413, Max-Change = 0.0006, gam = 0.0374, Max-Change = 0.0018, gam = 0.0342, Max-Change = 0.0023, gam = 0.0316, Max-Change = 0.0016, gam = 0.0294, Max-Change = 0.0028, gam = 0.0276, Max-Change = 0.0004, gam = 0.0260, Max-Change = 0.0009, gam = 0.0246, Max-Change = 0.0004
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(rmod2)
#> 
#> Call:
#> mixedmirt(data = data, covdata = covdata, model = 1, random = list(~1 | 
#>     group, ~1 | items))
#> 
#> --------------
#> FIXED EFFECTS:
#>             Estimate Std.Error z.value
#> (Intercept)   -0.605     0.013 -47.455
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>        F1
#> F1 0.0593
#> 
#> $group
#>           COV_group
#> COV_group     0.969
#> 
#> $items
#>           COV_items
#> COV_items     0.694
#> 
eff <- randef(rmod2) #estimate random effects

# random slopes with fixed intercepts (suppressed correlation)
rmod3 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ -1 + pseudoIQ|group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1864, Max-Change = 0.1623, Max-Change = 0.1436, Max-Change = 0.1243, Max-Change = 0.1035, Max-Change = 0.0953, Max-Change = 0.0862, Max-Change = 0.0730, Max-Change = 0.0612, Max-Change = 0.0490, Max-Change = 0.0426, Max-Change = 0.0342, Max-Change = 0.0359, Max-Change = 0.0255, Max-Change = 0.0251, Max-Change = 0.0175, Max-Change = 0.0105, Max-Change = 0.0121, Max-Change = 0.0184, Max-Change = 0.0194, Max-Change = 0.0092, Max-Change = 0.0064, Max-Change = 0.0083, Max-Change = 0.0065, Max-Change = 0.0026, Max-Change = 0.0066, Max-Change = 0.0110, Max-Change = 0.0004, Max-Change = 0.0044, Max-Change = 0.0104, Max-Change = 0.0118, Max-Change = 0.0071, Max-Change = 0.0056, Max-Change = 0.0073, Max-Change = 0.0134, Max-Change = 0.0043, Max-Change = 0.0143, Max-Change = 0.0042, Max-Change = 0.0109, Max-Change = 0.0061, Max-Change = 0.0096, Max-Change = 0.0130, Max-Change = 0.0158, Max-Change = 0.0124, Max-Change = 0.0075, Max-Change = 0.0080, Max-Change = 0.0056, Max-Change = 0.0025, Max-Change = 0.0109, Max-Change = 0.0030, Max-Change = 0.0135, Max-Change = 0.0156, Max-Change = 0.0136, Max-Change = 0.0085, Max-Change = 0.0083, Max-Change = 0.0118, Max-Change = 0.0017, Max-Change = 0.0074, Max-Change = 0.0278, Max-Change = 0.0030, Max-Change = 0.0086, Max-Change = 0.0053, Max-Change = 0.0006, Max-Change = 0.0139, Max-Change = 0.0106, Max-Change = 0.0060, Max-Change = 0.0048, Max-Change = 0.0084, Max-Change = 0.0076, Max-Change = 0.0094, Max-Change = 0.0090, Max-Change = 0.0120, Max-Change = 0.0036, Max-Change = 0.0029, Max-Change = 0.0089, Max-Change = 0.0036, Max-Change = 0.0030, Max-Change = 0.0074, Max-Change = 0.0044, Max-Change = 0.0153, Max-Change = 0.0043, Max-Change = 0.0295, Max-Change = 0.0054, Max-Change = 0.0105, Max-Change = 0.0154, Max-Change = 0.0094, Max-Change = 0.0021, Max-Change = 0.0077, Max-Change = 0.0113, Max-Change = 0.0035, Max-Change = 0.0045, Max-Change = 0.0065, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.0620, Max-Change = 0.0394, Max-Change = 0.0281, Max-Change = 0.0205, Max-Change = 0.0199, Max-Change = 0.0517, Max-Change = 0.0185, Max-Change = 0.0067, Max-Change = 0.0141, Max-Change = 0.0170, Max-Change = 0.0190, Max-Change = 0.0213, Max-Change = 0.0206, Max-Change = 0.0102, Max-Change = 0.0144, Max-Change = 0.0149, Max-Change = 0.0186, Max-Change = 0.0133, Max-Change = 0.0153, Max-Change = 0.0083, Max-Change = 0.0201, Max-Change = 0.0119, Max-Change = 0.0101, Max-Change = 0.0102, Max-Change = 0.0116, Max-Change = 0.0103, Max-Change = 0.0138, Max-Change = 0.0114, Max-Change = 0.0096, Max-Change = 0.0143, Max-Change = 0.0112, Max-Change = 0.0052, Max-Change = 0.0129, Max-Change = 0.0050, Max-Change = 0.0069, Max-Change = 0.0115, Max-Change = 0.0161, Max-Change = 0.0211, Max-Change = 0.0159, Max-Change = 0.0214, Max-Change = 0.0115, Max-Change = 0.0161, Max-Change = 0.0133, Max-Change = 0.0057, Max-Change = 0.0090, Max-Change = 0.0073, Max-Change = 0.0074, Max-Change = 0.0091, Max-Change = 0.0059, Max-Change = 0.0164, Max-Change = 0.0378, Max-Change = 0.0162, Max-Change = 0.0071, Max-Change = 0.0116, Max-Change = 0.0195, Max-Change = 0.0175, Max-Change = 0.0311, Max-Change = 0.0127, Max-Change = 0.0213, Max-Change = 0.0165, Max-Change = 0.0228, Max-Change = 0.0124, Max-Change = 0.0231, Max-Change = 0.0031, Max-Change = 0.0051, Max-Change = 0.0071, Max-Change = 0.0128, Max-Change = 0.0067, Max-Change = 0.0103, Max-Change = 0.0062, Max-Change = 0.0078, Max-Change = 0.0059, Max-Change = 0.0112, Max-Change = 0.0211, Max-Change = 0.0079, Max-Change = 0.0113, Max-Change = 0.0155, Max-Change = 0.0214, Max-Change = 0.0076, Max-Change = 0.0141, Max-Change = 0.0077, Max-Change = 0.0132, Max-Change = 0.0065, Max-Change = 0.0096, Max-Change = 0.0197, Max-Change = 0.0292, Max-Change = 0.0085, Max-Change = 0.0121, Max-Change = 0.0045, Max-Change = 0.0164, Max-Change = 0.0094, Max-Change = 0.0174, Max-Change = 0.0066, Max-Change = 0.0051, Max-Change = 0.0148, Max-Change = 0.0141, Max-Change = 0.0092, Max-Change = 0.0097, Max-Change = 0.0092, Max-Change = 0.0047, Max-Change = 0.0197, Max-Change = 0.0115, Max-Change = 0.0174, Max-Change = 0.0189, Max-Change = 0.0048, Max-Change = 0.0090, Max-Change = 0.0168, Max-Change = 0.0043, Max-Change = 0.0088, Max-Change = 0.0175, Max-Change = 0.0223, Max-Change = 0.0080, Max-Change = 0.0114, Max-Change = 0.0086, Max-Change = 0.0130, Max-Change = 0.0102, Max-Change = 0.0142, Max-Change = 0.0136, Max-Change = 0.0089, Max-Change = 0.0116, Max-Change = 0.0048, Max-Change = 0.0155, Max-Change = 0.0401, Max-Change = 0.0167, Max-Change = 0.0126, Max-Change = 0.0013, Max-Change = 0.0050, Max-Change = 0.0060, Max-Change = 0.0117, Max-Change = 0.0052, Max-Change = 0.0086, Max-Change = 0.0101, Max-Change = 0.0164, Max-Change = 0.0097, Max-Change = 0.0056, Max-Change = 0.0179, Max-Change = 0.0056, Max-Change = 0.0072, Max-Change = 0.0118, Max-Change = 0.0031, Max-Change = 0.0076, Max-Change = 0.0171, Max-Change = 0.0125, Max-Change = 0.0159, Max-Change = 0.0089, Max-Change = 0.0146, Max-Change = 0.0109, Max-Change = 0.0075, Max-Change = 0.0102, Max-Change = 0.0110, Max-Change = 0.0174, Max-Change = 0.0163, Max-Change = 0.0087, Max-Change = 0.0054, Max-Change = 0.0102, Max-Change = 0.0127, Max-Change = 0.0176, Max-Change = 0.0243, Max-Change = 0.0061, Max-Change = 0.0035, Max-Change = 0.0121, Max-Change = 0.0134, Max-Change = 0.0163, Max-Change = 0.0043, Max-Change = 0.0094, Max-Change = 0.0096, Max-Change = 0.0183, Max-Change = 0.0114, Max-Change = 0.0099, Max-Change = 0.0313, Max-Change = 0.0248, Max-Change = 0.0163, Max-Change = 0.0125, Max-Change = 0.0056, Max-Change = 0.0097, Max-Change = 0.0210, Max-Change = 0.0163, Max-Change = 0.0101, Max-Change = 0.0173, Max-Change = 0.0129, Max-Change = 0.0087, Max-Change = 0.0150, Max-Change = 0.0198, Max-Change = 0.0080, Max-Change = 0.0213, Max-Change = 0.0164, Max-Change = 0.0118, Max-Change = 0.0174, Max-Change = 0.0076, Max-Change = 0.0076, Max-Change = 0.0159, Max-Change = 0.0115, Max-Change = 0.0215, Max-Change = 0.0067, Max-Change = 0.0148, Max-Change = 0.0073, Max-Change = 0.0075, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0050, gam = 0.1057, Max-Change = 0.0160, gam = 0.0780, Max-Change = 0.0089, gam = 0.0629, Max-Change = 0.0076, gam = 0.0532, Max-Change = 0.0028, gam = 0.0464, Max-Change = 0.0033, gam = 0.0413, Max-Change = 0.0039, gam = 0.0374, Max-Change = 0.0038, gam = 0.0342, Max-Change = 0.0013, gam = 0.0316, Max-Change = 0.0024, gam = 0.0294, Max-Change = 0.0012, gam = 0.0276, Max-Change = 0.0031, gam = 0.0260, Max-Change = 0.0046, gam = 0.0246, Max-Change = 0.0009, gam = 0.0233, Max-Change = 0.0039, gam = 0.0222, Max-Change = 0.0026, gam = 0.0212, Max-Change = 0.0005, gam = 0.0203, Max-Change = 0.0010, gam = 0.0195, Max-Change = 0.0004
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(rmod3)
#> 
#> Call:
#> mixedmirt(data = data, covdata = covdata, model = 1, fixed = ~0 + 
#>     items, random = ~-1 + pseudoIQ | group)
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>       F1
#> F1 0.073
#> 
#> $group
#>              COV_group COV_pseudoIQ
#> COV_group        0.903         0.00
#> COV_pseudoIQ     0.000         0.16
#> 
eff <- randef(rmod3)
str(eff)
#> List of 2
#>  $ Theta: num [1:750, 1] -0.0802 0.0483 0.0384 0.0142 -0.0596 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr "F1"
#>  $ group: num [1:50, 1:2] -1.75 -1.23 -1.46 -1.65 -1.54 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:50] "G1" "G2" "G3" "G4" ...
#>   .. ..$ : chr [1:2] "group" "pseudoIQ"

###################################################
## LLTM, and 2PL version of LLTM
data(SAT12)
data <- key2binary(SAT12,
                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
model <- 'Theta = 1-32'

# for unconditional intercept comparison
mod <- mirt(data, model, itemtype='Rasch')
coef(mod, simplify=TRUE)
#> $items
#>         a1      d g u
#> Item.1   1 -1.077 0 1
#> Item.2   1  0.331 0 1
#> Item.3   1 -1.096 0 1
#> Item.4   1 -0.574 0 1
#> Item.5   1  0.581 0 1
#> Item.6   1 -1.912 0 1
#> Item.7   1  1.342 0 1
#> Item.8   1 -1.593 0 1
#> Item.9   1  2.324 0 1
#> Item.10  1 -0.362 0 1
#> Item.11  1  4.449 0 1
#> Item.12  1 -0.394 0 1
#> Item.13  1  0.791 0 1
#> Item.14  1  1.124 0 1
#> Item.15  1  1.724 0 1
#> Item.16  1 -0.402 0 1
#> Item.17  1  3.620 0 1
#> Item.18  1 -0.709 0 1
#> Item.19  1  0.237 0 1
#> Item.20  1  2.204 0 1
#> Item.21  1  2.684 0 1
#> Item.22  1  2.991 0 1
#> Item.23  1 -0.910 0 1
#> Item.24  1  1.153 0 1
#> Item.25  1 -0.590 0 1
#> Item.26  1 -0.179 0 1
#> Item.27  1  2.094 0 1
#> Item.28  1  0.150 0 1
#> Item.29  1 -0.769 0 1
#> Item.30  1 -0.274 0 1
#> Item.31  1  1.852 0 1
#> Item.32  1 -1.899 0 1
#> 
#> $means
#> Theta 
#>     0 
#> 
#> $cov
#>       Theta
#> Theta  0.82
#> 

# Suppose that the first 16 items were suspected to be easier than the last 16 items,
#   and we wish to test this item structure hypothesis (more intercept designs are possible
#   by including more columns).
itemdesign <- data.frame(itemorder = factor(c(rep('easier', 16), rep('harder', 16))))
rownames(itemdesign) <- colnames(data)
itemdesign
#>         itemorder
#> Item.1     easier
#> Item.2     easier
#> Item.3     easier
#> Item.4     easier
#> Item.5     easier
#> Item.6     easier
#> Item.7     easier
#> Item.8     easier
#> Item.9     easier
#> Item.10    easier
#> Item.11    easier
#> Item.12    easier
#> Item.13    easier
#> Item.14    easier
#> Item.15    easier
#> Item.16    easier
#> Item.17    harder
#> Item.18    harder
#> Item.19    harder
#> Item.20    harder
#> Item.21    harder
#> Item.22    harder
#> Item.23    harder
#> Item.24    harder
#> Item.25    harder
#> Item.26    harder
#> Item.27    harder
#> Item.28    harder
#> Item.29    harder
#> Item.30    harder
#> Item.31    harder
#> Item.32    harder

# notice that the 'fixed = ~ ... + items' argument is omitted
LLTM <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemdesign = itemdesign,
   SE = TRUE) # SE argument ensures that the information matrix is computed accurately
#> , Max-Change = 0.2000, Max-Change = 0.1211, Max-Change = 0.0434, Max-Change = 0.0359, Max-Change = 0.0327, Max-Change = 0.0388, Max-Change = 0.0141, Max-Change = 0.0267, Max-Change = 0.0141, Max-Change = 0.0104, Max-Change = 0.0137, Max-Change = 0.0148, Max-Change = 0.0114, Max-Change = 0.0116, Max-Change = 0.0112, Max-Change = 0.0118, Max-Change = 0.0093, Max-Change = 0.0096, Max-Change = 0.0059, Max-Change = 0.0056, Max-Change = 0.0114, Max-Change = 0.0124, Max-Change = 0.0062, Max-Change = 0.0059, Max-Change = 0.0072, Max-Change = 0.0089, Max-Change = 0.0046, Max-Change = 0.0061, Max-Change = 0.0051, Max-Change = 0.0022, Max-Change = 0.0050, Max-Change = 0.0068, Max-Change = 0.0016, Max-Change = 0.0028, Max-Change = 0.0049, Max-Change = 0.0067, Max-Change = 0.0019, Max-Change = 0.0038, Max-Change = 0.0054, Max-Change = 0.0082, Max-Change = 0.0038, Max-Change = 0.0019, Max-Change = 0.0055, Max-Change = 0.0044, Max-Change = 0.0065, Max-Change = 0.0020, Max-Change = 0.0034, Max-Change = 0.0029, Max-Change = 0.0062, Max-Change = 0.0021, Max-Change = 0.0029, Max-Change = 0.0023, Max-Change = 0.0030, Max-Change = 0.0008, Max-Change = 0.0130, Max-Change = 0.0023, Max-Change = 0.0062, Max-Change = 0.0028, Max-Change = 0.0057, Max-Change = 0.0035, Max-Change = 0.0029, Max-Change = 0.0017, Max-Change = 0.0020, Max-Change = 0.0012, Max-Change = 0.0021, Max-Change = 0.0028, Max-Change = 0.0012, Max-Change = 0.0020, Max-Change = 0.0019, Max-Change = 0.0036, Max-Change = 0.0029, Max-Change = 0.0029, Max-Change = 0.0051, Max-Change = 0.0044, Max-Change = 0.0053, Max-Change = 0.0021, Max-Change = 0.0016, Max-Change = 0.0037, Max-Change = 0.0021, Max-Change = 0.0044, Max-Change = 0.0073, Max-Change = 0.0023, Max-Change = 0.0070, Max-Change = 0.0025, Max-Change = 0.0043, Max-Change = 0.0035, Max-Change = 0.0026, Max-Change = 0.0011, Max-Change = 0.0024, Max-Change = 0.0019, Max-Change = 0.0024, Max-Change = 0.0031, Max-Change = 0.0037, Max-Change = 0.0065, Max-Change = 0.0016, Max-Change = 0.0041, Max-Change = 0.0015, Max-Change = 0.0028, Max-Change = 0.0038, Max-Change = 0.0031, Max-Change = 0.0016, Max-Change = 0.0069, Max-Change = 0.0031, Max-Change = 0.0026, Max-Change = 0.0034, Max-Change = 0.0028, Max-Change = 0.0035, Max-Change = 0.0050, Max-Change = 0.0047, Max-Change = 0.0036, Max-Change = 0.0028, Max-Change = 0.0040, Max-Change = 0.0049, Max-Change = 0.0024, Max-Change = 0.0062, Max-Change = 0.0044, Max-Change = 0.0023, Max-Change = 0.0040, Max-Change = 0.0037, Max-Change = 0.0031, Max-Change = 0.0024, Max-Change = 0.0037, Max-Change = 0.0040, Max-Change = 0.0030, Max-Change = 0.0023, Max-Change = 0.0020, Max-Change = 0.0022, Max-Change = 0.0064, Max-Change = 0.0008, Max-Change = 0.0028, Max-Change = 0.0063, Max-Change = 0.0073, Max-Change = 0.0020, Max-Change = 0.0055, Max-Change = 0.0041, Max-Change = 0.0032, Max-Change = 0.0033, Max-Change = 0.0030, Max-Change = 0.0027, Max-Change = 0.0036, Max-Change = 0.0049, Max-Change = 0.0020, Max-Change = 0.0024, Max-Change = 0.0041, Max-Change = 0.0026, Max-Change = 0.0031, Max-Change = 0.0018, Max-Change = 0.0002, Max-Change = 0.0031, Max-Change = 0.0031, Max-Change = 0.0054, Max-Change = 0.0037, Max-Change = 0.0008, Max-Change = 0.0040, Max-Change = 0.0019, Max-Change = 0.0045, Max-Change = 0.0029, Max-Change = 0.0032, Max-Change = 0.0029, Max-Change = 0.0013, Max-Change = 0.0032, Max-Change = 0.0034, Max-Change = 0.0036, Max-Change = 0.0046, Max-Change = 0.0044, Max-Change = 0.0075, Max-Change = 0.0045, Max-Change = 0.0043, Max-Change = 0.0046, Max-Change = 0.0012, Max-Change = 0.0027, Max-Change = 0.0060, Max-Change = 0.0018, Max-Change = 0.0031, Max-Change = 0.0019, Max-Change = 0.0085, Max-Change = 0.0035, Max-Change = 0.0013, Max-Change = 0.0054, Max-Change = 0.0023, Max-Change = 0.0016, Max-Change = 0.0031, Max-Change = 0.0021, Max-Change = 0.0027, Max-Change = 0.0059, Max-Change = 0.0024, Max-Change = 0.0027, Max-Change = 0.0041, Max-Change = 0.0032, Max-Change = 0.0014, Max-Change = 0.0021, Max-Change = 0.0075, Max-Change = 0.0044, Max-Change = 0.0061, Max-Change = 0.0026, Max-Change = 0.0014, Max-Change = 0.0076, Max-Change = 0.0022, Max-Change = 0.0014, Max-Change = 0.0052, Max-Change = 0.0035, Max-Change = 0.0014, Max-Change = 0.0029, Max-Change = 0.0028, Max-Change = 0.0006, Max-Change = 0.0050, Max-Change = 0.0027, Max-Change = 0.0054, Max-Change = 0.0038, Max-Change = 0.0051, Max-Change = 0.0040, Max-Change = 0.0028, Max-Change = 0.0055, Max-Change = 0.0050, Max-Change = 0.0020, Max-Change = 0.0037, Max-Change = 0.0067, Max-Change = 0.0031, Max-Change = 0.0011, Max-Change = 0.0024, Max-Change = 0.0075, Max-Change = 0.0025, Max-Change = 0.0044, Max-Change = 0.0043, Max-Change = 0.0024, Max-Change = 0.0034, Max-Change = 0.0008, Max-Change = 0.0025, Max-Change = 0.0021, Max-Change = 0.0031, Max-Change = 0.0038, Max-Change = 0.0012, Max-Change = 0.0064, Max-Change = 0.0045, Max-Change = 0.0014, Max-Change = 0.0030, Max-Change = 0.0028, Max-Change = 0.0024, Max-Change = 0.0055, Max-Change = 0.0039, Max-Change = 0.0063, Max-Change = 0.0022, Max-Change = 0.0013, Max-Change = 0.0015, Max-Change = 0.0027, Max-Change = 0.0006, Max-Change = 0.0032, Max-Change = 0.0020, Max-Change = 0.0071, Max-Change = 0.0010, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0050, gam = 0.1057, Max-Change = 0.0015, gam = 0.0780, Max-Change = 0.0014, gam = 0.0629, Max-Change = 0.0005, gam = 0.0532, Max-Change = 0.0005, gam = 0.0464, Max-Change = 0.0008
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(LLTM)
#> 
#> Call:
#> mixedmirt(data = data, model = model, fixed = ~0 + itemorder, 
#>     itemdesign = itemdesign, SE = TRUE)
#> 
#> --------------
#> FIXED EFFECTS:
#>                 Estimate Std.Error z.value
#> itemordereasier    0.165     0.029   5.746
#> itemorderharder    0.456     0.029  15.757
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>       Theta
#> Theta 0.359
#> 
coef(LLTM)
#> $Item.1
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.2
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.3
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.4
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.5
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.6
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.7
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.8
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.9
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.10
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.11
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.12
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.13
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.14
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.15
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.16
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.17
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.18
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.19
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.20
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.21
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.22
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.23
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.24
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.25
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.26
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.27
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.28
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.29
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.30
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.31
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $Item.32
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.165           0.456  1  0  0  1
#> CI_2.5            0.109           0.400 NA NA NA NA
#> CI_97.5           0.221           0.513 NA NA NA NA
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0  0.359
#> CI_2.5      NA  0.300
#> CI_97.5     NA  0.417
#> 
wald(LLTM)
#> itemordereasier.1.7.13.19.25.31.37.43.49.55.61.67.73.79.85.91.97.103.109.115.121.127.133.139.145.151.157.163.169.175.181.187 
#>                                                                                                                        0.165 
#> itemorderharder.2.8.14.20.26.32.38.44.50.56.62.68.74.80.86.92.98.104.110.116.122.128.134.140.146.152.158.164.170.176.182.188 
#>                                                                                                                        0.456 
#>                                                                                                                   COV_11.194 
#>                                                                                                                        0.359 
L <- matrix(c(-1, 1, 0), 1)
wald(LLTM, L) #first half different from second
#>        W df p
#> 1 92.085  1 0

# compare to items with estimated slopes (2PL)
twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL',
                   itemdesign = itemdesign)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1442, Max-Change = 0.1334, Max-Change = 0.1011, Max-Change = 0.1513, Max-Change = 0.1542, Max-Change = 0.1931, Max-Change = 0.1260, Max-Change = 0.1033, Max-Change = 0.0642, Max-Change = 0.0382, Max-Change = 0.0362, Max-Change = 0.0938, Max-Change = 0.1657, Max-Change = 0.0519, Max-Change = 0.0473, Max-Change = 0.0912, Max-Change = 0.1021, Max-Change = 0.0542, Max-Change = 0.0314, Max-Change = 0.0865, Max-Change = 0.1497, Max-Change = 0.0749, Max-Change = 0.0938, Max-Change = 0.0806, Max-Change = 0.0620, Max-Change = 0.1095, Max-Change = 0.0995, Max-Change = 0.0615, Max-Change = 0.0376, Max-Change = 0.0565, Max-Change = 0.1273, Max-Change = 0.0580, Max-Change = 0.0892, Max-Change = 0.1018, Max-Change = 0.1197, Max-Change = 0.0277, Max-Change = 0.0370, Max-Change = 0.0888, Max-Change = 0.1287, Max-Change = 0.0706, Max-Change = 0.0717, Max-Change = 0.1049, Max-Change = 0.0358, Max-Change = 0.0426, Max-Change = 0.0958, Max-Change = 0.0737, Max-Change = 0.0686, Max-Change = 0.0378, Max-Change = 0.0280, Max-Change = 0.0456, Max-Change = 0.0404, Max-Change = 0.0484, Max-Change = 0.0489, Max-Change = 0.0659, Max-Change = 0.0907, Max-Change = 0.0830, Max-Change = 0.0461, Max-Change = 0.0758, Max-Change = 0.0964, Max-Change = 0.0503, Max-Change = 0.1839, Max-Change = 0.0628, Max-Change = 0.0445, Max-Change = 0.0308, Max-Change = 0.1108, Max-Change = 0.0750, Max-Change = 0.0491, Max-Change = 0.1908, Max-Change = 0.0673, Max-Change = 0.1160, Max-Change = 0.0179, Max-Change = 0.0481, Max-Change = 0.1369, Max-Change = 0.0348, Max-Change = 0.0430, Max-Change = 0.0454, Max-Change = 0.0662, Max-Change = 0.1296, Max-Change = 0.0452, Max-Change = 0.0481, Max-Change = 0.1019, Max-Change = 0.0639, Max-Change = 0.0883, Max-Change = 0.0408, Max-Change = 0.1144, Max-Change = 0.0848, Max-Change = 0.0220, Max-Change = 0.0204, Max-Change = 0.1042, Max-Change = 0.0727, Max-Change = 0.0235, Max-Change = 0.1152, Max-Change = 0.1283, Max-Change = 0.2000, Max-Change = 0.1024, Max-Change = 0.0790, Max-Change = 0.0530, Max-Change = 0.1599, Max-Change = 0.0882, Max-Change = 0.0701, Max-Change = 0.0432, Max-Change = 0.0908, Max-Change = 0.0740, Max-Change = 0.0306, Max-Change = 0.1226, Max-Change = 0.0799, Max-Change = 0.0216, Max-Change = 0.0709, Max-Change = 0.0772, Max-Change = 0.0427, Max-Change = 0.0283, Max-Change = 0.1534, Max-Change = 0.0620, Max-Change = 0.1452, Max-Change = 0.0580, Max-Change = 0.1788, Max-Change = 0.0190, Max-Change = 0.0533, Max-Change = 0.0625, Max-Change = 0.0330, Max-Change = 0.0542, Max-Change = 0.1262, Max-Change = 0.0511, Max-Change = 0.0660, Max-Change = 0.0970, Max-Change = 0.0928, Max-Change = 0.0803, Max-Change = 0.0857, Max-Change = 0.0280, Max-Change = 0.1594, Max-Change = 0.0863, Max-Change = 0.0378, Max-Change = 0.0311, Max-Change = 0.0886, Max-Change = 0.0715, Max-Change = 0.1207, Max-Change = 0.0359, Max-Change = 0.0906, Max-Change = 0.1337, Max-Change = 0.0909, Max-Change = 0.0324, Max-Change = 0.1398, Max-Change = 0.1323, Max-Change = 0.0913, Max-Change = 0.0333, Max-Change = 0.0595, Max-Change = 0.0449, Max-Change = 0.1344, Max-Change = 0.2000, Max-Change = 0.0278, Max-Change = 0.0459, Max-Change = 0.1154, Max-Change = 0.0319, Max-Change = 0.0714, Max-Change = 0.1099, Max-Change = 0.0942, Max-Change = 0.0636, Max-Change = 0.0819, Max-Change = 0.0860, Max-Change = 0.0412, Max-Change = 0.2000, Max-Change = 0.0422, Max-Change = 0.1284, Max-Change = 0.0594, Max-Change = 0.0843, Max-Change = 0.0775, Max-Change = 0.0348, Max-Change = 0.1056, Max-Change = 0.0353, Max-Change = 0.0363, Max-Change = 0.1896, Max-Change = 0.0535, Max-Change = 0.0978, Max-Change = 0.0959, Max-Change = 0.0827, Max-Change = 0.0291, Max-Change = 0.1093, Max-Change = 0.0862, Max-Change = 0.1114, Max-Change = 0.0440, Max-Change = 0.0922, Max-Change = 0.0546, Max-Change = 0.0655, Max-Change = 0.1098, Max-Change = 0.0620, Max-Change = 0.0964, Max-Change = 0.1176, Max-Change = 0.0687, Max-Change = 0.0395, Max-Change = 0.1172, Max-Change = 0.1350, Max-Change = 0.1815, Max-Change = 0.0760, Max-Change = 0.0625, Max-Change = 0.0623, Max-Change = 0.1998, Max-Change = 0.0979, Max-Change = 0.1779, Max-Change = 0.0804, Max-Change = 0.1349, Max-Change = 0.0441, Max-Change = 0.0772, Max-Change = 0.1331, Max-Change = 0.1961, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0468, gam = 0.1057, Max-Change = 0.0182, gam = 0.0780, Max-Change = 0.0097, gam = 0.0629, Max-Change = 0.0159, gam = 0.0532, Max-Change = 0.0087, gam = 0.0464, Max-Change = 0.0160, gam = 0.0413, Max-Change = 0.0233, gam = 0.0374, Max-Change = 0.0305, gam = 0.0342, Max-Change = 0.0305, gam = 0.0316, Max-Change = 0.0061, gam = 0.0294, Max-Change = 0.0088, gam = 0.0276, Max-Change = 0.0146, gam = 0.0260, Max-Change = 0.0119, gam = 0.0246, Max-Change = 0.0078, gam = 0.0233, Max-Change = 0.0162, gam = 0.0222, Max-Change = 0.0076, gam = 0.0212, Max-Change = 0.0183, gam = 0.0203, Max-Change = 0.0075, gam = 0.0195, Max-Change = 0.0028, gam = 0.0188, Max-Change = 0.0026, gam = 0.0181, Max-Change = 0.0044, gam = 0.0175, Max-Change = 0.0026, gam = 0.0169, Max-Change = 0.0029, gam = 0.0164, Max-Change = 0.0043, gam = 0.0159, Max-Change = 0.0071, gam = 0.0154, Max-Change = 0.0045, gam = 0.0150, Max-Change = 0.0026, gam = 0.0146, Max-Change = 0.0028, gam = 0.0142, Max-Change = 0.0063, gam = 0.0139, Max-Change = 0.0102, gam = 0.0135, Max-Change = 0.0182, gam = 0.0132, Max-Change = 0.0051, gam = 0.0129, Max-Change = 0.0043, gam = 0.0126, Max-Change = 0.0050, gam = 0.0124, Max-Change = 0.0074, gam = 0.0121, Max-Change = 0.0018, gam = 0.0119, Max-Change = 0.0029, gam = 0.0116, Max-Change = 0.0115, gam = 0.0114, Max-Change = 0.0030, gam = 0.0112, Max-Change = 0.0032, gam = 0.0110, Max-Change = 0.0031, gam = 0.0108, Max-Change = 0.0034, gam = 0.0106, Max-Change = 0.0090, gam = 0.0104, Max-Change = 0.0032, gam = 0.0102, Max-Change = 0.0048, gam = 0.0101, Max-Change = 0.0014, gam = 0.0099, Max-Change = 0.0040, gam = 0.0098, Max-Change = 0.0042, gam = 0.0096, Max-Change = 0.0083, gam = 0.0095, Max-Change = 0.0041, gam = 0.0093, Max-Change = 0.0031, gam = 0.0092, Max-Change = 0.0018, gam = 0.0091, Max-Change = 0.0060, gam = 0.0089, Max-Change = 0.0014, gam = 0.0088, Max-Change = 0.0051, gam = 0.0087, Max-Change = 0.0093, gam = 0.0086, Max-Change = 0.0104, gam = 0.0085, Max-Change = 0.0017, gam = 0.0084, Max-Change = 0.0053, gam = 0.0082, Max-Change = 0.0017, gam = 0.0081, Max-Change = 0.0012, gam = 0.0080, Max-Change = 0.0032, gam = 0.0080, Max-Change = 0.0047, gam = 0.0079, Max-Change = 0.0016, gam = 0.0078, Max-Change = 0.0063, gam = 0.0077, Max-Change = 0.0021, gam = 0.0076, Max-Change = 0.0013, gam = 0.0075, Max-Change = 0.0016, gam = 0.0074, Max-Change = 0.0012, gam = 0.0073, Max-Change = 0.0017, gam = 0.0073, Max-Change = 0.0027, gam = 0.0072, Max-Change = 0.0019, gam = 0.0071, Max-Change = 0.0018, gam = 0.0070, Max-Change = 0.0038, gam = 0.0070, Max-Change = 0.0019, gam = 0.0069, Max-Change = 0.0021, gam = 0.0068, Max-Change = 0.0023, gam = 0.0068, Max-Change = 0.0064, gam = 0.0067, Max-Change = 0.0011, gam = 0.0066, Max-Change = 0.0022, gam = 0.0066, Max-Change = 0.0022, gam = 0.0065, Max-Change = 0.0008, gam = 0.0065, Max-Change = 0.0018, gam = 0.0064, Max-Change = 0.0019, gam = 0.0064, Max-Change = 0.0021, gam = 0.0063, Max-Change = 0.0035, gam = 0.0062, Max-Change = 0.0023, gam = 0.0062, Max-Change = 0.0014, gam = 0.0061, Max-Change = 0.0013, gam = 0.0061, Max-Change = 0.0028, gam = 0.0060, Max-Change = 0.0055, gam = 0.0060, Max-Change = 0.0030, gam = 0.0059, Max-Change = 0.0029, gam = 0.0059, Max-Change = 0.0025, gam = 0.0058, Max-Change = 0.0026, gam = 0.0058, Max-Change = 0.0021, gam = 0.0058, Max-Change = 0.0017, gam = 0.0057, Max-Change = 0.0016, gam = 0.0057, Max-Change = 0.0056, gam = 0.0056, Max-Change = 0.0021, gam = 0.0056, Max-Change = 0.0011, gam = 0.0055, Max-Change = 0.0018, gam = 0.0055, Max-Change = 0.0034, gam = 0.0055, Max-Change = 0.0019, gam = 0.0054, Max-Change = 0.0007, gam = 0.0054, Max-Change = 0.0020, gam = 0.0053, Max-Change = 0.0019, gam = 0.0053, Max-Change = 0.0032, gam = 0.0053, Max-Change = 0.0030, gam = 0.0052, Max-Change = 0.0024, gam = 0.0052, Max-Change = 0.0029, gam = 0.0052, Max-Change = 0.0014, gam = 0.0051, Max-Change = 0.0017, gam = 0.0051, Max-Change = 0.0045, gam = 0.0051, Max-Change = 0.0023, gam = 0.0050, Max-Change = 0.0024, gam = 0.0050, Max-Change = 0.0018, gam = 0.0050, Max-Change = 0.0020, gam = 0.0049, Max-Change = 0.0011, gam = 0.0049, Max-Change = 0.0010, gam = 0.0049, Max-Change = 0.0007, gam = 0.0048, Max-Change = 0.0018, gam = 0.0048, Max-Change = 0.0034, gam = 0.0048, Max-Change = 0.0037, gam = 0.0048, Max-Change = 0.0018, gam = 0.0047, Max-Change = 0.0018, gam = 0.0047, Max-Change = 0.0005, gam = 0.0047, Max-Change = 0.0036, gam = 0.0046, Max-Change = 0.0017, gam = 0.0046, Max-Change = 0.0010, gam = 0.0046, Max-Change = 0.0020, gam = 0.0046, Max-Change = 0.0016, gam = 0.0045, Max-Change = 0.0022, gam = 0.0045, Max-Change = 0.0008, gam = 0.0045, Max-Change = 0.0031, gam = 0.0045, Max-Change = 0.0011, gam = 0.0044, Max-Change = 0.0013, gam = 0.0044, Max-Change = 0.0016, gam = 0.0044, Max-Change = 0.0031, gam = 0.0044, Max-Change = 0.0027, gam = 0.0043, Max-Change = 0.0013, gam = 0.0043, Max-Change = 0.0033, gam = 0.0043, Max-Change = 0.0022, gam = 0.0043, Max-Change = 0.0024, gam = 0.0043, Max-Change = 0.0007, gam = 0.0042, Max-Change = 0.0015, gam = 0.0042, Max-Change = 0.0030, gam = 0.0042, Max-Change = 0.0014, gam = 0.0042, Max-Change = 0.0015, gam = 0.0041, Max-Change = 0.0012, gam = 0.0041, Max-Change = 0.0013, gam = 0.0041, Max-Change = 0.0032, gam = 0.0041, Max-Change = 0.0025, gam = 0.0041, Max-Change = 0.0008, gam = 0.0040, Max-Change = 0.0015, gam = 0.0040, Max-Change = 0.0012, gam = 0.0040, Max-Change = 0.0010, gam = 0.0040, Max-Change = 0.0017, gam = 0.0040, Max-Change = 0.0020, gam = 0.0040, Max-Change = 0.0017, gam = 0.0039, Max-Change = 0.0031, gam = 0.0039, Max-Change = 0.0007, gam = 0.0039, Max-Change = 0.0017, gam = 0.0039, Max-Change = 0.0013, gam = 0.0039, Max-Change = 0.0014, gam = 0.0038, Max-Change = 0.0016, gam = 0.0038, Max-Change = 0.0010, gam = 0.0038, Max-Change = 0.0015, gam = 0.0038, Max-Change = 0.0006, gam = 0.0038, Max-Change = 0.0020, gam = 0.0038, Max-Change = 0.0009, gam = 0.0037, Max-Change = 0.0018, gam = 0.0037, Max-Change = 0.0020, gam = 0.0037, Max-Change = 0.0009, gam = 0.0037, Max-Change = 0.0017, gam = 0.0037, Max-Change = 0.0011, gam = 0.0037, Max-Change = 0.0004, gam = 0.0036, Max-Change = 0.0011, gam = 0.0036, Max-Change = 0.0011, gam = 0.0036, Max-Change = 0.0013, gam = 0.0036, Max-Change = 0.0039, gam = 0.0036, Max-Change = 0.0023, gam = 0.0036, Max-Change = 0.0013, gam = 0.0036, Max-Change = 0.0013, gam = 0.0035, Max-Change = 0.0012, gam = 0.0035, Max-Change = 0.0028, gam = 0.0035, Max-Change = 0.0011, gam = 0.0035, Max-Change = 0.0021, gam = 0.0035, Max-Change = 0.0020, gam = 0.0035, Max-Change = 0.0018, gam = 0.0035, Max-Change = 0.0012, gam = 0.0034, Max-Change = 0.0007, gam = 0.0034, Max-Change = 0.0011, gam = 0.0034, Max-Change = 0.0011, gam = 0.0034, Max-Change = 0.0007, gam = 0.0034, Max-Change = 0.0034, gam = 0.0034, Max-Change = 0.0012, gam = 0.0034, Max-Change = 0.0008, gam = 0.0034, Max-Change = 0.0009, gam = 0.0033, Max-Change = 0.0007
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
# twoPL not mixing too well (AR should be between .2 and .5), decrease MHcand
twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL',
                  itemdesign = itemdesign, technical = list(MHcand = 0.8))
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1580, Max-Change = 0.1493, Max-Change = 0.1590, Max-Change = 0.1453, Max-Change = 0.1497, Max-Change = 0.2000, Max-Change = 0.1382, Max-Change = 0.1473, Max-Change = 0.0963, Max-Change = 0.0840, Max-Change = 0.0516, Max-Change = 0.0416, Max-Change = 0.0567, Max-Change = 0.0582, Max-Change = 0.0715, Max-Change = 0.0502, Max-Change = 0.0950, Max-Change = 0.0412, Max-Change = 0.0287, Max-Change = 0.0941, Max-Change = 0.1107, Max-Change = 0.1413, Max-Change = 0.1043, Max-Change = 0.0567, Max-Change = 0.0544, Max-Change = 0.0555, Max-Change = 0.0547, Max-Change = 0.0563, Max-Change = 0.0415, Max-Change = 0.0593, Max-Change = 0.0620, Max-Change = 0.0563, Max-Change = 0.1559, Max-Change = 0.1268, Max-Change = 0.0883, Max-Change = 0.1939, Max-Change = 0.0395, Max-Change = 0.0775, Max-Change = 0.1430, Max-Change = 0.2000, Max-Change = 0.0591, Max-Change = 0.0400, Max-Change = 0.0300, Max-Change = 0.0472, Max-Change = 0.0765, Max-Change = 0.0602, Max-Change = 0.0876, Max-Change = 0.0224, Max-Change = 0.0339, Max-Change = 0.0400, Max-Change = 0.0394, Max-Change = 0.0533, Max-Change = 0.0396, Max-Change = 0.0462, Max-Change = 0.0959, Max-Change = 0.0746, Max-Change = 0.1065, Max-Change = 0.1064, Max-Change = 0.1416, Max-Change = 0.1050, Max-Change = 0.0815, Max-Change = 0.0823, Max-Change = 0.2000, Max-Change = 0.0448, Max-Change = 0.0489, Max-Change = 0.0875, Max-Change = 0.0658, Max-Change = 0.1561, Max-Change = 0.1379, Max-Change = 0.1195, Max-Change = 0.1394, Max-Change = 0.1107, Max-Change = 0.0919, Max-Change = 0.0363, Max-Change = 0.0454, Max-Change = 0.0601, Max-Change = 0.0507, Max-Change = 0.0918, Max-Change = 0.0841, Max-Change = 0.0292, Max-Change = 0.0758, Max-Change = 0.0761, Max-Change = 0.1493, Max-Change = 0.0472, Max-Change = 0.1030, Max-Change = 0.0278, Max-Change = 0.0412, Max-Change = 0.0247, Max-Change = 0.0416, Max-Change = 0.1399, Max-Change = 0.1462, Max-Change = 0.2000, Max-Change = 0.0908, Max-Change = 0.2000, Max-Change = 0.1294, Max-Change = 0.0560, Max-Change = 0.0362, Max-Change = 0.1113, Max-Change = 0.0892, Max-Change = 0.0927, Max-Change = 0.0413, Max-Change = 0.0558, Max-Change = 0.0763, Max-Change = 0.1423, Max-Change = 0.0551, Max-Change = 0.0409, Max-Change = 0.0732, Max-Change = 0.1227, Max-Change = 0.0615, Max-Change = 0.1213, Max-Change = 0.1253, Max-Change = 0.0463, Max-Change = 0.0600, Max-Change = 0.1780, Max-Change = 0.1414, Max-Change = 0.0771, Max-Change = 0.1539, Max-Change = 0.0594, Max-Change = 0.0397, Max-Change = 0.0837, Max-Change = 0.0690, Max-Change = 0.0433, Max-Change = 0.0624, Max-Change = 0.0244, Max-Change = 0.0493, Max-Change = 0.1388, Max-Change = 0.0777, Max-Change = 0.0834, Max-Change = 0.0300, Max-Change = 0.0289, Max-Change = 0.0614, Max-Change = 0.1249, Max-Change = 0.0689, Max-Change = 0.0364, Max-Change = 0.0754, Max-Change = 0.1376, Max-Change = 0.0851, Max-Change = 0.0522, Max-Change = 0.0723, Max-Change = 0.0945, Max-Change = 0.1353, Max-Change = 0.0725, Max-Change = 0.0754, Max-Change = 0.1172, Max-Change = 0.0907, Max-Change = 0.2000, Max-Change = 0.0900, Max-Change = 0.0957, Max-Change = 0.1489, Max-Change = 0.0403, Max-Change = 0.0805, Max-Change = 0.1271, Max-Change = 0.0559, Max-Change = 0.0699, Max-Change = 0.0948, Max-Change = 0.0515, Max-Change = 0.0648, Max-Change = 0.1273, Max-Change = 0.0964, Max-Change = 0.0342, Max-Change = 0.1370, Max-Change = 0.0237, Max-Change = 0.0425, Max-Change = 0.0494, Max-Change = 0.0508, Max-Change = 0.0394, Max-Change = 0.1035, Max-Change = 0.0739, Max-Change = 0.0707, Max-Change = 0.0908, Max-Change = 0.0850, Max-Change = 0.0558, Max-Change = 0.1103, Max-Change = 0.1808, Max-Change = 0.1738, Max-Change = 0.1432, Max-Change = 0.0501, Max-Change = 0.2000, Max-Change = 0.0585, Max-Change = 0.1046, Max-Change = 0.1443, Max-Change = 0.1918, Max-Change = 0.0686, Max-Change = 0.0578, Max-Change = 0.0658, Max-Change = 0.0790, Max-Change = 0.1219, Max-Change = 0.1074, Max-Change = 0.1143, Max-Change = 0.1944, Max-Change = 0.2000, Max-Change = 0.0279, Max-Change = 0.1737, Max-Change = 0.0589, Max-Change = 0.1507, Max-Change = 0.1349, Max-Change = 0.0647, Max-Change = 0.1399, Max-Change = 0.0462, Max-Change = 0.1044, Max-Change = 0.0908, Max-Change = 0.0564, Max-Change = 0.0726, Max-Change = 0.1225, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0746, gam = 0.1057, Max-Change = 0.0272, gam = 0.0780, Max-Change = 0.0412, gam = 0.0629, Max-Change = 0.0188, gam = 0.0532, Max-Change = 0.0088, gam = 0.0464, Max-Change = 0.0202, gam = 0.0413, Max-Change = 0.0106, gam = 0.0374, Max-Change = 0.0108, gam = 0.0342, Max-Change = 0.0181, gam = 0.0316, Max-Change = 0.0109, gam = 0.0294, Max-Change = 0.0220, gam = 0.0276, Max-Change = 0.0024, gam = 0.0260, Max-Change = 0.0170, gam = 0.0246, Max-Change = 0.0048, gam = 0.0233, Max-Change = 0.0220, gam = 0.0222, Max-Change = 0.0166, gam = 0.0212, Max-Change = 0.0132, gam = 0.0203, Max-Change = 0.0060, gam = 0.0195, Max-Change = 0.0056, gam = 0.0188, Max-Change = 0.0016, gam = 0.0181, Max-Change = 0.0061, gam = 0.0175, Max-Change = 0.0037, gam = 0.0169, Max-Change = 0.0085, gam = 0.0164, Max-Change = 0.0039, gam = 0.0159, Max-Change = 0.0031, gam = 0.0154, Max-Change = 0.0112, gam = 0.0150, Max-Change = 0.0045, gam = 0.0146, Max-Change = 0.0044, gam = 0.0142, Max-Change = 0.0059, gam = 0.0139, Max-Change = 0.0052, gam = 0.0135, Max-Change = 0.0117, gam = 0.0132, Max-Change = 0.0068, gam = 0.0129, Max-Change = 0.0021, gam = 0.0126, Max-Change = 0.0019, gam = 0.0124, Max-Change = 0.0023, gam = 0.0121, Max-Change = 0.0069, gam = 0.0119, Max-Change = 0.0021, gam = 0.0116, Max-Change = 0.0033, gam = 0.0114, Max-Change = 0.0018, gam = 0.0112, Max-Change = 0.0051, gam = 0.0110, Max-Change = 0.0019, gam = 0.0108, Max-Change = 0.0025, gam = 0.0106, Max-Change = 0.0068, gam = 0.0104, Max-Change = 0.0052, gam = 0.0102, Max-Change = 0.0087, gam = 0.0101, Max-Change = 0.0064, gam = 0.0099, Max-Change = 0.0088, gam = 0.0098, Max-Change = 0.0082, gam = 0.0096, Max-Change = 0.0075, gam = 0.0095, Max-Change = 0.0048, gam = 0.0093, Max-Change = 0.0023, gam = 0.0092, Max-Change = 0.0030, gam = 0.0091, Max-Change = 0.0031, gam = 0.0089, Max-Change = 0.0021, gam = 0.0088, Max-Change = 0.0045, gam = 0.0087, Max-Change = 0.0025, gam = 0.0086, Max-Change = 0.0039, gam = 0.0085, Max-Change = 0.0046, gam = 0.0084, Max-Change = 0.0100, gam = 0.0082, Max-Change = 0.0049, gam = 0.0081, Max-Change = 0.0032, gam = 0.0080, Max-Change = 0.0046, gam = 0.0080, Max-Change = 0.0049, gam = 0.0079, Max-Change = 0.0015, gam = 0.0078, Max-Change = 0.0015, gam = 0.0077, Max-Change = 0.0033, gam = 0.0076, Max-Change = 0.0011, gam = 0.0075, Max-Change = 0.0021, gam = 0.0074, Max-Change = 0.0011, gam = 0.0073, Max-Change = 0.0030, gam = 0.0073, Max-Change = 0.0022, gam = 0.0072, Max-Change = 0.0013, gam = 0.0071, Max-Change = 0.0025, gam = 0.0070, Max-Change = 0.0018, gam = 0.0070, Max-Change = 0.0041, gam = 0.0069, Max-Change = 0.0017, gam = 0.0068, Max-Change = 0.0009, gam = 0.0068, Max-Change = 0.0009, gam = 0.0067, Max-Change = 0.0033, gam = 0.0066, Max-Change = 0.0023, gam = 0.0066, Max-Change = 0.0013, gam = 0.0065, Max-Change = 0.0049, gam = 0.0065, Max-Change = 0.0013, gam = 0.0064, Max-Change = 0.0012, gam = 0.0064, Max-Change = 0.0048, gam = 0.0063, Max-Change = 0.0051, gam = 0.0062, Max-Change = 0.0041, gam = 0.0062, Max-Change = 0.0081, gam = 0.0061, Max-Change = 0.0026, gam = 0.0061, Max-Change = 0.0031, gam = 0.0060, Max-Change = 0.0029, gam = 0.0060, Max-Change = 0.0040, gam = 0.0059, Max-Change = 0.0034, gam = 0.0059, Max-Change = 0.0018, gam = 0.0058, Max-Change = 0.0015, gam = 0.0058, Max-Change = 0.0046, gam = 0.0058, Max-Change = 0.0012, gam = 0.0057, Max-Change = 0.0017, gam = 0.0057, Max-Change = 0.0007, gam = 0.0056, Max-Change = 0.0009, gam = 0.0056, Max-Change = 0.0025, gam = 0.0055, Max-Change = 0.0017, gam = 0.0055, Max-Change = 0.0020, gam = 0.0055, Max-Change = 0.0016, gam = 0.0054, Max-Change = 0.0013, gam = 0.0054, Max-Change = 0.0015, gam = 0.0053, Max-Change = 0.0018, gam = 0.0053, Max-Change = 0.0012, gam = 0.0053, Max-Change = 0.0023, gam = 0.0052, Max-Change = 0.0026, gam = 0.0052, Max-Change = 0.0011, gam = 0.0052, Max-Change = 0.0024, gam = 0.0051, Max-Change = 0.0017, gam = 0.0051, Max-Change = 0.0011, gam = 0.0051, Max-Change = 0.0007, gam = 0.0050, Max-Change = 0.0030, gam = 0.0050, Max-Change = 0.0017, gam = 0.0050, Max-Change = 0.0009, gam = 0.0049, Max-Change = 0.0018, gam = 0.0049, Max-Change = 0.0022, gam = 0.0049, Max-Change = 0.0015, gam = 0.0048, Max-Change = 0.0026, gam = 0.0048, Max-Change = 0.0009, gam = 0.0048, Max-Change = 0.0028, gam = 0.0048, Max-Change = 0.0064, gam = 0.0047, Max-Change = 0.0041, gam = 0.0047, Max-Change = 0.0020, gam = 0.0047, Max-Change = 0.0019, gam = 0.0046, Max-Change = 0.0017, gam = 0.0046, Max-Change = 0.0010, gam = 0.0046, Max-Change = 0.0030, gam = 0.0046, Max-Change = 0.0026, gam = 0.0045, Max-Change = 0.0017, gam = 0.0045, Max-Change = 0.0027, gam = 0.0045, Max-Change = 0.0016, gam = 0.0045, Max-Change = 0.0031, gam = 0.0044, Max-Change = 0.0011, gam = 0.0044, Max-Change = 0.0027, gam = 0.0044, Max-Change = 0.0007, gam = 0.0044, Max-Change = 0.0012, gam = 0.0043, Max-Change = 0.0018, gam = 0.0043, Max-Change = 0.0026, gam = 0.0043, Max-Change = 0.0032, gam = 0.0043, Max-Change = 0.0051, gam = 0.0043, Max-Change = 0.0007, gam = 0.0042, Max-Change = 0.0020, gam = 0.0042, Max-Change = 0.0013, gam = 0.0042, Max-Change = 0.0013, gam = 0.0042, Max-Change = 0.0018, gam = 0.0041, Max-Change = 0.0016, gam = 0.0041, Max-Change = 0.0010, gam = 0.0041, Max-Change = 0.0010, gam = 0.0041, Max-Change = 0.0012, gam = 0.0041, Max-Change = 0.0021, gam = 0.0040, Max-Change = 0.0011, gam = 0.0040, Max-Change = 0.0021, gam = 0.0040, Max-Change = 0.0007, gam = 0.0040, Max-Change = 0.0016, gam = 0.0040, Max-Change = 0.0016, gam = 0.0040, Max-Change = 0.0009, gam = 0.0039, Max-Change = 0.0005, gam = 0.0039, Max-Change = 0.0017, gam = 0.0039, Max-Change = 0.0007, gam = 0.0039, Max-Change = 0.0017, gam = 0.0039, Max-Change = 0.0009, gam = 0.0038, Max-Change = 0.0022, gam = 0.0038, Max-Change = 0.0009, gam = 0.0038, Max-Change = 0.0010, gam = 0.0038, Max-Change = 0.0026, gam = 0.0038, Max-Change = 0.0039, gam = 0.0038, Max-Change = 0.0021, gam = 0.0037, Max-Change = 0.0011, gam = 0.0037, Max-Change = 0.0017, gam = 0.0037, Max-Change = 0.0013, gam = 0.0037, Max-Change = 0.0016, gam = 0.0037, Max-Change = 0.0011, gam = 0.0037, Max-Change = 0.0024, gam = 0.0036, Max-Change = 0.0034, gam = 0.0036, Max-Change = 0.0008, gam = 0.0036, Max-Change = 0.0013, gam = 0.0036, Max-Change = 0.0030, gam = 0.0036, Max-Change = 0.0009, gam = 0.0036, Max-Change = 0.0008, gam = 0.0036, Max-Change = 0.0018, gam = 0.0035, Max-Change = 0.0021, gam = 0.0035, Max-Change = 0.0004, gam = 0.0035, Max-Change = 0.0009, gam = 0.0035, Max-Change = 0.0005
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
anova(twoPL, LLTM) #much better fit
#>            AIC    SABIC       HQ      BIC    logLik        X2  df   p
#> twoPL 20484.33 20525.88 20542.52 20633.82 -10208.16                  
#> LLTM  25464.27 25467.94 25469.41 25477.47 -12729.14 -5041.948 -31 NaN
summary(twoPL)
#> 
#> Call:
#> mixedmirt(data = data, model = model, fixed = ~0 + itemorder, 
#>     itemtype = "2PL", itemdesign = itemdesign, technical = list(MHcand = 0.8))
#> 
#> --------------
#> FIXED EFFECTS:
#>                 Estimate Std.Error z.value
#> itemordereasier   -1.669     0.087 -19.127
#> itemorderharder   -1.644     0.095 -17.316
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>       Theta
#> Theta     1
#> 
coef(twoPL)
#> $Item.1
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 0.930  0  0  1
#> CI_2.5           -1.840          -1.831 0.696 NA NA NA
#> CI_97.5          -1.498          -1.458 1.165 NA NA NA
#> 
#> $Item.2
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 2.526  0  0  1
#> CI_2.5           -1.840          -1.831 2.201 NA NA NA
#> CI_97.5          -1.498          -1.458 2.851 NA NA NA
#> 
#> $Item.3
#>         itemordereasier itemorderharder   a1  d  g  u
#> par              -1.669          -1.644 1.00  0  0  1
#> CI_2.5           -1.840          -1.831 0.77 NA NA NA
#> CI_97.5          -1.498          -1.458 1.23 NA NA NA
#> 
#> $Item.4
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.288  0  0  1
#> CI_2.5           -1.840          -1.831 1.035 NA NA NA
#> CI_97.5          -1.498          -1.458 1.540 NA NA NA
#> 
#> $Item.5
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 2.704  0  0  1
#> CI_2.5           -1.840          -1.831 2.359 NA NA NA
#> CI_97.5          -1.498          -1.458 3.049 NA NA NA
#> 
#> $Item.6
#>         itemordereasier itemorderharder   a1  d  g  u
#> par              -1.669          -1.644 0.42  0  0  1
#> CI_2.5           -1.840          -1.831 0.19 NA NA NA
#> CI_97.5          -1.498          -1.458 0.65 NA NA NA
#> 
#> $Item.7
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 3.913  0  0  1
#> CI_2.5           -1.840          -1.831 3.461 NA NA NA
#> CI_97.5          -1.498          -1.458 4.365 NA NA NA
#> 
#> $Item.8
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 0.513  0  0  1
#> CI_2.5           -1.840          -1.831 0.282 NA NA NA
#> CI_97.5          -1.498          -1.458 0.744 NA NA NA
#> 
#> $Item.9
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 5.455  0  0  1
#> CI_2.5           -1.840          -1.831 4.840 NA NA NA
#> CI_97.5          -1.498          -1.458 6.071 NA NA NA
#> 
#> $Item.10
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.628  0  0  1
#> CI_2.5           -1.840          -1.831 1.352 NA NA NA
#> CI_97.5          -1.498          -1.458 1.905 NA NA NA
#> 
#> $Item.11
#>         itemordereasier itemorderharder     a1  d  g  u
#> par              -1.669          -1.644 13.512  0  0  1
#> CI_2.5           -1.840          -1.831 10.680 NA NA NA
#> CI_97.5          -1.498          -1.458 16.344 NA NA NA
#> 
#> $Item.12
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.245  0  0  1
#> CI_2.5           -1.840          -1.831 0.974 NA NA NA
#> CI_97.5          -1.498          -1.458 1.515 NA NA NA
#> 
#> $Item.13
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 3.042  0  0  1
#> CI_2.5           -1.840          -1.831 2.668 NA NA NA
#> CI_97.5          -1.498          -1.458 3.416 NA NA NA
#> 
#> $Item.14
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 3.488  0  0  1
#> CI_2.5           -1.840          -1.831 3.076 NA NA NA
#> CI_97.5          -1.498          -1.458 3.899 NA NA NA
#> 
#> $Item.15
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 4.749  0  0  1
#> CI_2.5           -1.840          -1.831 4.225 NA NA NA
#> CI_97.5          -1.498          -1.458 5.272 NA NA NA
#> 
#> $Item.16
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.469  0  0  1
#> CI_2.5           -1.840          -1.831 1.205 NA NA NA
#> CI_97.5          -1.498          -1.458 1.733 NA NA NA
#> 
#> $Item.17
#>         itemordereasier itemorderharder     a1  d  g  u
#> par              -1.669          -1.644  9.952  0  0  1
#> CI_2.5           -1.840          -1.831  8.451 NA NA NA
#> CI_97.5          -1.498          -1.458 11.454 NA NA NA
#> 
#> $Item.18
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.408  0  0  1
#> CI_2.5           -1.840          -1.831 1.145 NA NA NA
#> CI_97.5          -1.498          -1.458 1.670 NA NA NA
#> 
#> $Item.19
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 2.155  0  0  1
#> CI_2.5           -1.840          -1.831 1.825 NA NA NA
#> CI_97.5          -1.498          -1.458 2.486 NA NA NA
#> 
#> $Item.20
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 5.787  0  0  1
#> CI_2.5           -1.840          -1.831 5.118 NA NA NA
#> CI_97.5          -1.498          -1.458 6.456 NA NA NA
#> 
#> $Item.21
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 6.345  0  0  1
#> CI_2.5           -1.840          -1.831 5.595 NA NA NA
#> CI_97.5          -1.498          -1.458 7.094 NA NA NA
#> 
#> $Item.22
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 7.826  0  0  1
#> CI_2.5           -1.840          -1.831 6.827 NA NA NA
#> CI_97.5          -1.498          -1.458 8.824 NA NA NA
#> 
#> $Item.23
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 0.982  0  0  1
#> CI_2.5           -1.840          -1.831 0.732 NA NA NA
#> CI_97.5          -1.498          -1.458 1.231 NA NA NA
#> 
#> $Item.24
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 3.622  0  0  1
#> CI_2.5           -1.840          -1.831 3.186 NA NA NA
#> CI_97.5          -1.498          -1.458 4.057 NA NA NA
#> 
#> $Item.25
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.293  0  0  1
#> CI_2.5           -1.840          -1.831 1.021 NA NA NA
#> CI_97.5          -1.498          -1.458 1.565 NA NA NA
#> 
#> $Item.26
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.907  0  0  1
#> CI_2.5           -1.840          -1.831 1.608 NA NA NA
#> CI_97.5          -1.498          -1.458 2.205 NA NA NA
#> 
#> $Item.27
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 5.656  0  0  1
#> CI_2.5           -1.840          -1.831 5.018 NA NA NA
#> CI_97.5          -1.498          -1.458 6.293 NA NA NA
#> 
#> $Item.28
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 2.120  0  0  1
#> CI_2.5           -1.840          -1.831 1.807 NA NA NA
#> CI_97.5          -1.498          -1.458 2.433 NA NA NA
#> 
#> $Item.29
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.165  0  0  1
#> CI_2.5           -1.840          -1.831 0.910 NA NA NA
#> CI_97.5          -1.498          -1.458 1.420 NA NA NA
#> 
#> $Item.30
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 1.442  0  0  1
#> CI_2.5           -1.840          -1.831 1.156 NA NA NA
#> CI_97.5          -1.498          -1.458 1.727 NA NA NA
#> 
#> $Item.31
#>         itemordereasier itemorderharder    a1  d  g  u
#> par              -1.669          -1.644 5.235  0  0  1
#> CI_2.5           -1.840          -1.831 4.637 NA NA NA
#> CI_97.5          -1.498          -1.458 5.832 NA NA NA
#> 
#> $Item.32
#>         itemordereasier itemorderharder     a1  d  g  u
#> par              -1.669          -1.644  0.097  0  0  1
#> CI_2.5           -1.840          -1.831 -0.163 NA NA NA
#> CI_97.5          -1.498          -1.458  0.357 NA NA NA
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0      1
#> CI_2.5      NA     NA
#> CI_97.5     NA     NA
#> 

wald(twoPL)
#> itemordereasier.1.7.13.19.25.31.37.43.49.55.61.67.73.79.85.91.97.103.109.115.121.127.133.139.145.151.157.163.169.175.181.187 
#>                                                                                                                       -1.669 
#> itemorderharder.2.8.14.20.26.32.38.44.50.56.62.68.74.80.86.92.98.104.110.116.122.128.134.140.146.152.158.164.170.176.182.188 
#>                                                                                                                       -1.644 
#>                                                                                                                         a1.3 
#>                                                                                                                        0.930 
#>                                                                                                                         a1.9 
#>                                                                                                                        2.526 
#>                                                                                                                        a1.15 
#>                                                                                                                        1.000 
#>                                                                                                                        a1.21 
#>                                                                                                                        1.288 
#>                                                                                                                        a1.27 
#>                                                                                                                        2.704 
#>                                                                                                                        a1.33 
#>                                                                                                                        0.420 
#>                                                                                                                        a1.39 
#>                                                                                                                        3.913 
#>                                                                                                                        a1.45 
#>                                                                                                                        0.513 
#>                                                                                                                        a1.51 
#>                                                                                                                        5.455 
#>                                                                                                                        a1.57 
#>                                                                                                                        1.628 
#>                                                                                                                        a1.63 
#>                                                                                                                       13.512 
#>                                                                                                                        a1.69 
#>                                                                                                                        1.245 
#>                                                                                                                        a1.75 
#>                                                                                                                        3.042 
#>                                                                                                                        a1.81 
#>                                                                                                                        3.488 
#>                                                                                                                        a1.87 
#>                                                                                                                        4.749 
#>                                                                                                                        a1.93 
#>                                                                                                                        1.469 
#>                                                                                                                        a1.99 
#>                                                                                                                        9.952 
#>                                                                                                                       a1.105 
#>                                                                                                                        1.408 
#>                                                                                                                       a1.111 
#>                                                                                                                        2.155 
#>                                                                                                                       a1.117 
#>                                                                                                                        5.787 
#>                                                                                                                       a1.123 
#>                                                                                                                        6.345 
#>                                                                                                                       a1.129 
#>                                                                                                                        7.826 
#>                                                                                                                       a1.135 
#>                                                                                                                        0.982 
#>                                                                                                                       a1.141 
#>                                                                                                                        3.622 
#>                                                                                                                       a1.147 
#>                                                                                                                        1.293 
#>                                                                                                                       a1.153 
#>                                                                                                                        1.907 
#>                                                                                                                       a1.159 
#>                                                                                                                        5.656 
#>                                                                                                                       a1.165 
#>                                                                                                                        2.120 
#>                                                                                                                       a1.171 
#>                                                                                                                        1.165 
#>                                                                                                                       a1.177 
#>                                                                                                                        1.442 
#>                                                                                                                       a1.183 
#>                                                                                                                        5.235 
#>                                                                                                                       a1.189 
#>                                                                                                                        0.097 
L <- matrix(0, 1, 34)
L[1, 1] <- 1
L[1, 2] <- -1
wald(twoPL, L) # n.s., which is the correct conclusion. Rasch approach gave wrong inference
#>       W df     p
#> 1 0.075  1 0.785

## LLTM with item error term
LLTMwithError <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, random = ~ 1|items,
    itemdesign = itemdesign)
#> , Max-Change = 0.2000, Max-Change = 0.1211, Max-Change = 0.0434, Max-Change = 0.0359, Max-Change = 0.0327, Max-Change = 0.0388, Max-Change = 0.0141, Max-Change = 0.0267, Max-Change = 0.0141, Max-Change = 0.0104, Max-Change = 0.0137, Max-Change = 0.0148, Max-Change = 0.0114, Max-Change = 0.0116, Max-Change = 0.0112, Max-Change = 0.0118, Max-Change = 0.0093, Max-Change = 0.0096, Max-Change = 0.0059, Max-Change = 0.0056, Max-Change = 0.0114, Max-Change = 0.0124, Max-Change = 0.0062, Max-Change = 0.0059, Max-Change = 0.0072, Max-Change = 0.0089, Max-Change = 0.0046, Max-Change = 0.0061, Max-Change = 0.0051, Max-Change = 0.0022, Max-Change = 0.0050, Max-Change = 0.0068, Max-Change = 0.0016, Max-Change = 0.0028, Max-Change = 0.0049, Max-Change = 0.0067, Max-Change = 0.0019, Max-Change = 0.0038, Max-Change = 0.0054, Max-Change = 0.0082, Max-Change = 0.0038, Max-Change = 0.0019, Max-Change = 0.0055, Max-Change = 0.0044, Max-Change = 0.0065, Max-Change = 0.0020, Max-Change = 0.0034, Max-Change = 0.0029, Max-Change = 0.0062, Max-Change = 0.0021, Max-Change = 0.0029, Max-Change = 0.0023, Max-Change = 0.0030, Max-Change = 0.0008, Max-Change = 0.0130, Max-Change = 0.0023, Max-Change = 0.0062, Max-Change = 0.0028, Max-Change = 0.0057, Max-Change = 0.0035, Max-Change = 0.0029, Max-Change = 0.0017, Max-Change = 0.0020, Max-Change = 0.0012, Max-Change = 0.0021, Max-Change = 0.0028, Max-Change = 0.0012, Max-Change = 0.0020, Max-Change = 0.0019, Max-Change = 0.0036, Max-Change = 0.0029, Max-Change = 0.0029, Max-Change = 0.0051, Max-Change = 0.0044, Max-Change = 0.0053, Max-Change = 0.0021, Max-Change = 0.0016, Max-Change = 0.0037, Max-Change = 0.0021, Max-Change = 0.0044, Max-Change = 0.0073, Max-Change = 0.0023, Max-Change = 0.0070, Max-Change = 0.0025, Max-Change = 0.0043, Max-Change = 0.0035, Max-Change = 0.0026, Max-Change = 0.0011, Max-Change = 0.0024, Max-Change = 0.0019, Max-Change = 0.0024, Max-Change = 0.0031, Max-Change = 0.0037, Max-Change = 0.0065, Max-Change = 0.0016, Max-Change = 0.0041, Max-Change = 0.0015, Max-Change = 0.0028, Max-Change = 0.0038, Max-Change = 0.0162, Max-Change = 0.0703, Max-Change = 0.0721, Max-Change = 0.0756, Max-Change = 0.0767, Max-Change = 0.0753, Max-Change = 0.0733, Max-Change = 0.0733, Max-Change = 0.0746, Max-Change = 0.0813, Max-Change = 0.0766, Max-Change = 0.0805, Max-Change = 0.0598, Max-Change = 0.0810, Max-Change = 0.0682, Max-Change = 0.0495, Max-Change = 0.0376, Max-Change = 0.0505, Max-Change = 0.0389, Max-Change = 0.0431, Max-Change = 0.0425, Max-Change = 0.0365, Max-Change = 0.0017, Max-Change = 0.0223, Max-Change = 0.0212, Max-Change = 0.0175, Max-Change = 0.0118, Max-Change = 0.0094, Max-Change = 0.0110, Max-Change = 0.0050, Max-Change = 0.0071, Max-Change = 0.0469, Max-Change = 0.0090, Max-Change = 0.0240, Max-Change = 0.0292, Max-Change = 0.0250, Max-Change = 0.0413, Max-Change = 0.0280, Max-Change = 0.0109, Max-Change = 0.0112, Max-Change = 0.0110, Max-Change = 0.0070, Max-Change = 0.0045, Max-Change = 0.0100, Max-Change = 0.0423, Max-Change = 0.0076, Max-Change = 0.0067, Max-Change = 0.0157, Max-Change = 0.0166, Max-Change = 0.0028, Max-Change = 0.0170, Max-Change = 0.0066, Max-Change = 0.0326, Max-Change = 0.0109, Max-Change = 0.0608, Max-Change = 0.0353, Max-Change = 0.0181, Max-Change = 0.0166, Max-Change = 0.0089, Max-Change = 0.0236, Max-Change = 0.0192, Max-Change = 0.0196, Max-Change = 0.0130, Max-Change = 0.0084, Max-Change = 0.0269, Max-Change = 0.0053, Max-Change = 0.0115, Max-Change = 0.0035, Max-Change = 0.0131, Max-Change = 0.0098, Max-Change = 0.0172, Max-Change = 0.0192, Max-Change = 0.0064, Max-Change = 0.0090, Max-Change = 0.0053, Max-Change = 0.0098, Max-Change = 0.0136, Max-Change = 0.0036, Max-Change = 0.0086, Max-Change = 0.0073, Max-Change = 0.0156, Max-Change = 0.0275, Max-Change = 0.0075, Max-Change = 0.0196, Max-Change = 0.0265, Max-Change = 0.0065, Max-Change = 0.0024, Max-Change = 0.0069, Max-Change = 0.0211, Max-Change = 0.0112, Max-Change = 0.0036, Max-Change = 0.0203, Max-Change = 0.0169, Max-Change = 0.0201, Max-Change = 0.0116, Max-Change = 0.0366, Max-Change = 0.0321, Max-Change = 0.0214, Max-Change = 0.0076, Max-Change = 0.0036, Max-Change = 0.0295, Max-Change = 0.0102, Max-Change = 0.0080, Max-Change = 0.0075, Max-Change = 0.0052, Max-Change = 0.0053, Max-Change = 0.0056, Max-Change = 0.0135, Max-Change = 0.0038, Max-Change = 0.0131, Max-Change = 0.0118, Max-Change = 0.0025, Max-Change = 0.0028, Max-Change = 0.0016, Max-Change = 0.0049, Max-Change = 0.0037, Max-Change = 0.0056, Max-Change = 0.0049, Max-Change = 0.0053, Max-Change = 0.0093, Max-Change = 0.0041, Max-Change = 0.0051, Max-Change = 0.0041, Max-Change = 0.0045, Max-Change = 0.0218, Max-Change = 0.0083, Max-Change = 0.0033, Max-Change = 0.0025, Max-Change = 0.0034, Max-Change = 0.0044, Max-Change = 0.0031, Max-Change = 0.0012, Max-Change = 0.0034, Max-Change = 0.0059, Max-Change = 0.0066, Max-Change = 0.0054, Max-Change = 0.0035, Max-Change = 0.0013, Max-Change = 0.0056, Max-Change = 0.0026, Max-Change = 0.0028, Max-Change = 0.0015, Max-Change = 0.0034, Max-Change = 0.0090, Max-Change = 0.0061, Max-Change = 0.0072, Max-Change = 0.0051, Max-Change = 0.0036, Max-Change = 0.0027, Max-Change = 0.0047, Max-Change = 0.0005, Max-Change = 0.0095, Max-Change = 0.0015, Max-Change = 0.0021, Max-Change = 0.0031, Max-Change = 0.0042, Max-Change = 0.0020, Max-Change = 0.0016, Max-Change = 0.0061, Max-Change = 0.0032, Max-Change = 0.0063, Max-Change = 0.0035, Max-Change = 0.0016, Max-Change = 0.0042, Max-Change = 0.0032, Max-Change = 0.0039, Max-Change = 0.0085, Max-Change = 0.0131, Max-Change = 0.0041, Max-Change = 0.0046, Max-Change = 0.0041, Max-Change = 0.0125, Max-Change = 0.0036, Max-Change = 0.0048, Max-Change = 0.0058, Max-Change = 0.0059, Max-Change = 0.0079, Max-Change = 0.0025, Max-Change = 0.0065, Max-Change = 0.0044, Max-Change = 0.0013, Max-Change = 0.0040, Max-Change = 0.0020, Max-Change = 0.0077, Max-Change = 0.0023, Max-Change = 0.0028, Max-Change = 0.0068, Max-Change = 0.0068, Max-Change = 0.0090, Max-Change = 0.0024, Max-Change = 0.0017, Max-Change = 0.0051, Max-Change = 0.0049, Max-Change = 0.0061, Max-Change = 0.0052, Max-Change = 0.0062, Max-Change = 0.0030, Max-Change = 0.0112, Max-Change = 0.0029, Max-Change = 0.0061, Max-Change = 0.0036, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0037, gam = 0.1057, Max-Change = 0.0038, gam = 0.0780, Max-Change = 0.0033, gam = 0.0629, Max-Change = 0.0032, gam = 0.0532, Max-Change = 0.0016, gam = 0.0464, Max-Change = 0.0012, gam = 0.0413, Max-Change = 0.0016, gam = 0.0374, Max-Change = 0.0008, gam = 0.0342, Max-Change = 0.0003, gam = 0.0316, Max-Change = 0.0010
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(LLTMwithError)
#> 
#> Call:
#> mixedmirt(data = data, model = model, fixed = ~0 + itemorder, 
#>     random = ~1 | items, itemdesign = itemdesign)
#> 
#> --------------
#> FIXED EFFECTS:
#>                 Estimate Std.Error z.value
#> itemordereasier    0.147     0.153   0.961
#> itemorderharder    0.567     0.138   4.107
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>       Theta
#> Theta  0.77
#> 
#> $items
#>           COV_items
#> COV_items      2.36
#> 
# large item level variance after itemorder is regressed; not a great predictor of item difficulty
coef(LLTMwithError)
#> $Item.1
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.2
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.3
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.4
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.5
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.6
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.7
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.8
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.9
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.10
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.11
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.12
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.13
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.14
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.15
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.16
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.17
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.18
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.19
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.20
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.21
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.22
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.23
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.24
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.25
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.26
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.27
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.28
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.29
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.30
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.31
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $Item.32
#>         itemordereasier itemorderharder a1  d  g  u
#> par               0.147           0.567  1  0  0  1
#> CI_2.5           -0.153           0.296 NA NA NA NA
#> CI_97.5           0.446           0.838 NA NA NA NA
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0  0.770
#> CI_2.5      NA  0.656
#> CI_97.5     NA  0.884
#> 
#> $items
#>         COV_items_items
#> par               2.363
#> CI_2.5            1.229
#> CI_97.5           3.497
#> 

###################################################
### Polytomous example

# make an arbitrary group difference
covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))

# partial credit model
mod <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1962, Max-Change = 0.1682, Max-Change = 0.1133, Max-Change = 0.0674, Max-Change = 0.0992, Max-Change = 0.0598, Max-Change = 0.0391, Max-Change = 0.0767, Max-Change = 0.0756, Max-Change = 0.0498, Max-Change = 0.0341, Max-Change = 0.0267, Max-Change = 0.0180, Max-Change = 0.0456, Max-Change = 0.0402, Max-Change = 0.0274, Max-Change = 0.0393, Max-Change = 0.0554, Max-Change = 0.0284, Max-Change = 0.0071, Max-Change = 0.0122, Max-Change = 0.0220, Max-Change = 0.0069, Max-Change = 0.0146, Max-Change = 0.0283, Max-Change = 0.0373, Max-Change = 0.0308, Max-Change = 0.0094, Max-Change = 0.0156, Max-Change = 0.0788, Max-Change = 0.0158, Max-Change = 0.0216, Max-Change = 0.0166, Max-Change = 0.0117, Max-Change = 0.0380, Max-Change = 0.0473, Max-Change = 0.0301, Max-Change = 0.0510, Max-Change = 0.0303, Max-Change = 0.0534, Max-Change = 0.0369, Max-Change = 0.0366, Max-Change = 0.0164, Max-Change = 0.0112, Max-Change = 0.0184, Max-Change = 0.0241, Max-Change = 0.0132, Max-Change = 0.0130, Max-Change = 0.0322, Max-Change = 0.0117, Max-Change = 0.0257, Max-Change = 0.0167, Max-Change = 0.0103, Max-Change = 0.0230, Max-Change = 0.0250, Max-Change = 0.0330, Max-Change = 0.0080, Max-Change = 0.0185, Max-Change = 0.0132, Max-Change = 0.0584, Max-Change = 0.0139, Max-Change = 0.0264, Max-Change = 0.0230, Max-Change = 0.0302, Max-Change = 0.0573, Max-Change = 0.0190, Max-Change = 0.0394, Max-Change = 0.0492, Max-Change = 0.1002, Max-Change = 0.0201, Max-Change = 0.0600, Max-Change = 0.0258, Max-Change = 0.0563, Max-Change = 0.0101, Max-Change = 0.0197, Max-Change = 0.0245, Max-Change = 0.0342, Max-Change = 0.0495, Max-Change = 0.0175, Max-Change = 0.0301, Max-Change = 0.0319, Max-Change = 0.0743, Max-Change = 0.0160, Max-Change = 0.0126, Max-Change = 0.0248, Max-Change = 0.0238, Max-Change = 0.0409, Max-Change = 0.0393, Max-Change = 0.0248, Max-Change = 0.0402, Max-Change = 0.0136, Max-Change = 0.0208, Max-Change = 0.0514, Max-Change = 0.0607, Max-Change = 0.0259, Max-Change = 0.0346, Max-Change = 0.0313, Max-Change = 0.0295, Max-Change = 0.0630, Max-Change = 0.0644, Max-Change = 0.0463, Max-Change = 0.0292, Max-Change = 0.0367, Max-Change = 0.0403, Max-Change = 0.0443, Max-Change = 0.0677, Max-Change = 0.0219, Max-Change = 0.0139, Max-Change = 0.0540, Max-Change = 0.0469, Max-Change = 0.0200, Max-Change = 0.0142, Max-Change = 0.0141, Max-Change = 0.0355, Max-Change = 0.0165, Max-Change = 0.0266, Max-Change = 0.0173, Max-Change = 0.0206, Max-Change = 0.0325, Max-Change = 0.0488, Max-Change = 0.0090, Max-Change = 0.0374, Max-Change = 0.0068, Max-Change = 0.0110, Max-Change = 0.0160, Max-Change = 0.0323, Max-Change = 0.0151, Max-Change = 0.0225, Max-Change = 0.0657, Max-Change = 0.0597, Max-Change = 0.0197, Max-Change = 0.0156, Max-Change = 0.0257, Max-Change = 0.0453, Max-Change = 0.0137, Max-Change = 0.0347, Max-Change = 0.0347, Max-Change = 0.0123, Max-Change = 0.0305, Max-Change = 0.0093, Max-Change = 0.0487, Max-Change = 0.0134, Max-Change = 0.0418, Max-Change = 0.0272, Max-Change = 0.0292, Max-Change = 0.0300, Max-Change = 0.0573, Max-Change = 0.0277, Max-Change = 0.0086, Max-Change = 0.0474, Max-Change = 0.0222, Max-Change = 0.0189, Max-Change = 0.0280, Max-Change = 0.0296, Max-Change = 0.0060, Max-Change = 0.0725, Max-Change = 0.0465, Max-Change = 0.0170, Max-Change = 0.0170, Max-Change = 0.0480, Max-Change = 0.0664, Max-Change = 0.0191, Max-Change = 0.0223, Max-Change = 0.0139, Max-Change = 0.0465, Max-Change = 0.0288, Max-Change = 0.0670, Max-Change = 0.0177, Max-Change = 0.0684, Max-Change = 0.0415, Max-Change = 0.0246, Max-Change = 0.0121, Max-Change = 0.0224, Max-Change = 0.0199, Max-Change = 0.0073, Max-Change = 0.0258, Max-Change = 0.0154, Max-Change = 0.0388, Max-Change = 0.0320, Max-Change = 0.0259, Max-Change = 0.0501, Max-Change = 0.1062, Max-Change = 0.0204, Max-Change = 0.0840, Max-Change = 0.0111, Max-Change = 0.0397, Max-Change = 0.0191, Max-Change = 0.0168, Max-Change = 0.0256, Max-Change = 0.0382, Max-Change = 0.0726, Max-Change = 0.0670, Max-Change = 0.0159, Max-Change = 0.0099, Max-Change = 0.0558, Max-Change = 0.0305, Max-Change = 0.0040, Max-Change = 0.0227, Max-Change = 0.0270, Max-Change = 0.0105, Max-Change = 0.0755, Max-Change = 0.0284, Max-Change = 0.0125, Max-Change = 0.0220, Max-Change = 0.0491, Max-Change = 0.0176, Max-Change = 0.0107, Max-Change = 0.0212, Max-Change = 0.0090, Max-Change = 0.0251, Max-Change = 0.0252, Max-Change = 0.0374, Max-Change = 0.0572, Max-Change = 0.0189, Max-Change = 0.0261, Max-Change = 0.0190, Max-Change = 0.0336, Max-Change = 0.0352, Max-Change = 0.0230, Max-Change = 0.0160, Max-Change = 0.0347, Max-Change = 0.0476, Max-Change = 0.0099, Max-Change = 0.0521, Max-Change = 0.0576, Max-Change = 0.0144, Max-Change = 0.0152, Max-Change = 0.0292, Max-Change = 0.0164, Max-Change = 0.0307, Max-Change = 0.0539, Max-Change = 0.0581, Max-Change = 0.0183, Max-Change = 0.0614, Max-Change = 0.0283, Max-Change = 0.0063, Max-Change = 0.0422, Max-Change = 0.0528, Max-Change = 0.0289, Max-Change = 0.0444, Max-Change = 0.0230, Max-Change = 0.0127, Max-Change = 0.0158, Max-Change = 0.0254, Max-Change = 0.0421, Max-Change = 0.0088, Max-Change = 0.0437, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0321, gam = 0.1057, Max-Change = 0.0075, gam = 0.0780, Max-Change = 0.0111, gam = 0.0629, Max-Change = 0.0052, gam = 0.0532, Max-Change = 0.0125, gam = 0.0464, Max-Change = 0.0073, gam = 0.0413, Max-Change = 0.0042, gam = 0.0374, Max-Change = 0.0100, gam = 0.0342, Max-Change = 0.0023, gam = 0.0316, Max-Change = 0.0089, gam = 0.0294, Max-Change = 0.0029, gam = 0.0276, Max-Change = 0.0057, gam = 0.0260, Max-Change = 0.0063, gam = 0.0246, Max-Change = 0.0033, gam = 0.0233, Max-Change = 0.0053, gam = 0.0222, Max-Change = 0.0039, gam = 0.0212, Max-Change = 0.0055, gam = 0.0203, Max-Change = 0.0013, gam = 0.0195, Max-Change = 0.0041, gam = 0.0188, Max-Change = 0.0022, gam = 0.0181, Max-Change = 0.0027, gam = 0.0175, Max-Change = 0.0050, gam = 0.0169, Max-Change = 0.0055, gam = 0.0164, Max-Change = 0.0043, gam = 0.0159, Max-Change = 0.0011, gam = 0.0154, Max-Change = 0.0016, gam = 0.0150, Max-Change = 0.0014, gam = 0.0146, Max-Change = 0.0014, gam = 0.0142, Max-Change = 0.0016, gam = 0.0139, Max-Change = 0.0026, gam = 0.0135, Max-Change = 0.0041, gam = 0.0132, Max-Change = 0.0024, gam = 0.0129, Max-Change = 0.0024, gam = 0.0126, Max-Change = 0.0012, gam = 0.0124, Max-Change = 0.0010, gam = 0.0121, Max-Change = 0.0047, gam = 0.0119, Max-Change = 0.0034, gam = 0.0116, Max-Change = 0.0035, gam = 0.0114, Max-Change = 0.0021, gam = 0.0112, Max-Change = 0.0016, gam = 0.0110, Max-Change = 0.0016, gam = 0.0108, Max-Change = 0.0022, gam = 0.0106, Max-Change = 0.0027, gam = 0.0104, Max-Change = 0.0019, gam = 0.0102, Max-Change = 0.0012, gam = 0.0101, Max-Change = 0.0024, gam = 0.0099, Max-Change = 0.0007, gam = 0.0098, Max-Change = 0.0020, gam = 0.0096, Max-Change = 0.0020, gam = 0.0095, Max-Change = 0.0023, gam = 0.0093, Max-Change = 0.0021, gam = 0.0092, Max-Change = 0.0010, gam = 0.0091, Max-Change = 0.0006, gam = 0.0089, Max-Change = 0.0012, gam = 0.0088, Max-Change = 0.0012, gam = 0.0087, Max-Change = 0.0010, gam = 0.0086, Max-Change = 0.0014, gam = 0.0085, Max-Change = 0.0008, gam = 0.0084, Max-Change = 0.0007, gam = 0.0082, Max-Change = 0.0017, gam = 0.0081, Max-Change = 0.0018, gam = 0.0080, Max-Change = 0.0010, gam = 0.0080, Max-Change = 0.0007, gam = 0.0079, Max-Change = 0.0024, gam = 0.0078, Max-Change = 0.0003, gam = 0.0077, Max-Change = 0.0009, gam = 0.0076, Max-Change = 0.0014, gam = 0.0075, Max-Change = 0.0013, gam = 0.0074, Max-Change = 0.0022, gam = 0.0073, Max-Change = 0.0009, gam = 0.0073, Max-Change = 0.0011, gam = 0.0072, Max-Change = 0.0027, gam = 0.0071, Max-Change = 0.0010, gam = 0.0070, Max-Change = 0.0007, gam = 0.0070, Max-Change = 0.0008, gam = 0.0069, Max-Change = 0.0006
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
coef(mod)
#> $Comfort
#>         groupm a1 ak0 ak1 ak2 ak3 d0    d1    d2    d3
#> par     -0.084  1   0   1   2   3  0 3.098 5.718 4.364
#> CI_2.5  -0.336 NA  NA  NA  NA  NA NA 2.086 4.677 3.257
#> CI_97.5  0.167 NA  NA  NA  NA  NA NA 4.110 6.759 5.472
#> 
#> $Work
#>         groupm a1 ak0 ak1 ak2 ak3 d0    d1    d2    d3
#> par     -0.084  1   0   1   2   3  0 1.919 2.859 1.038
#> CI_2.5  -0.336 NA  NA  NA  NA  NA NA 1.455 2.311 0.352
#> CI_97.5  0.167 NA  NA  NA  NA  NA NA 2.382 3.406 1.723
#> 
#> $Future
#>         groupm a1 ak0 ak1 ak2 ak3 d0    d1    d2    d3
#> par     -0.084  1   0   1   2   3  0 2.665 4.114 3.013
#> CI_2.5  -0.336 NA  NA  NA  NA  NA NA 2.026 3.406 2.214
#> CI_97.5  0.167 NA  NA  NA  NA  NA NA 3.304 4.822 3.813
#> 
#> $Benefit
#>         groupm a1 ak0 ak1 ak2 ak3 d0    d1    d2    d3
#> par     -0.084  1   0   1   2   3  0 2.469 3.398 2.076
#> CI_2.5  -0.336 NA  NA  NA  NA  NA NA 1.932 2.779 1.349
#> CI_97.5  0.167 NA  NA  NA  NA  NA NA 3.006 4.016 2.803
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0  0.985
#> CI_2.5      NA  0.707
#> CI_97.5     NA  1.264
#> 

# gpcm to estimate slopes
mod2 <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group,
                 itemtype = 'gpcm')
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1776, Max-Change = 0.1068, Max-Change = 0.1180, Max-Change = 0.1121, Max-Change = 0.0996, Max-Change = 0.1305, Max-Change = 0.0338, Max-Change = 0.1041, Max-Change = 0.0588, Max-Change = 0.0460, Max-Change = 0.0885, Max-Change = 0.0954, Max-Change = 0.0450, Max-Change = 0.1303, Max-Change = 0.0976, Max-Change = 0.1155, Max-Change = 0.0919, Max-Change = 0.1057, Max-Change = 0.0780, Max-Change = 0.0979, Max-Change = 0.0675, Max-Change = 0.1021, Max-Change = 0.0427, Max-Change = 0.0814, Max-Change = 0.0625, Max-Change = 0.0662, Max-Change = 0.0919, Max-Change = 0.1059, Max-Change = 0.1302, Max-Change = 0.0629, Max-Change = 0.0813, Max-Change = 0.0782, Max-Change = 0.1577, Max-Change = 0.0302, Max-Change = 0.0778, Max-Change = 0.1032, Max-Change = 0.0458, Max-Change = 0.0519, Max-Change = 0.0532, Max-Change = 0.1538, Max-Change = 0.1558, Max-Change = 0.1330, Max-Change = 0.0485, Max-Change = 0.0976, Max-Change = 0.0576, Max-Change = 0.0236, Max-Change = 0.1744, Max-Change = 0.0174, Max-Change = 0.0716, Max-Change = 0.1091, Max-Change = 0.0652, Max-Change = 0.0722, Max-Change = 0.1248, Max-Change = 0.0696, Max-Change = 0.0928, Max-Change = 0.0199, Max-Change = 0.1060, Max-Change = 0.0868, Max-Change = 0.0901, Max-Change = 0.1216, Max-Change = 0.0595, Max-Change = 0.0360, Max-Change = 0.0834, Max-Change = 0.1512, Max-Change = 0.0561, Max-Change = 0.1249, Max-Change = 0.0744, Max-Change = 0.1373, Max-Change = 0.1111, Max-Change = 0.0304, Max-Change = 0.1097, Max-Change = 0.0966, Max-Change = 0.0830, Max-Change = 0.0232, Max-Change = 0.0628, Max-Change = 0.0403, Max-Change = 0.1226, Max-Change = 0.0310, Max-Change = 0.0457, Max-Change = 0.1221, Max-Change = 0.1108, Max-Change = 0.1636, Max-Change = 0.1611, Max-Change = 0.1112, Max-Change = 0.0778, Max-Change = 0.1025, Max-Change = 0.0633, Max-Change = 0.0344, Max-Change = 0.0296, Max-Change = 0.0352, Max-Change = 0.0954, Max-Change = 0.1280, Max-Change = 0.0501, Max-Change = 0.1303, Max-Change = 0.0549, Max-Change = 0.1357, Max-Change = 0.1346, Max-Change = 0.0546, Max-Change = 0.0694, Max-Change = 0.0904, Max-Change = 0.1042, Max-Change = 0.0233, Max-Change = 0.0545, Max-Change = 0.1438, Max-Change = 0.0567, Max-Change = 0.0701, Max-Change = 0.0924, Max-Change = 0.0772, Max-Change = 0.1517, Max-Change = 0.0887, Max-Change = 0.0719, Max-Change = 0.0368, Max-Change = 0.1098, Max-Change = 0.1057, Max-Change = 0.0267, Max-Change = 0.2000, Max-Change = 0.0747, Max-Change = 0.0943, Max-Change = 0.0283, Max-Change = 0.0432, Max-Change = 0.0995, Max-Change = 0.0564, Max-Change = 0.0250, Max-Change = 0.1342, Max-Change = 0.0971, Max-Change = 0.0706, Max-Change = 0.1048, Max-Change = 0.0733, Max-Change = 0.0640, Max-Change = 0.1025, Max-Change = 0.0336, Max-Change = 0.0938, Max-Change = 0.0412, Max-Change = 0.0341, Max-Change = 0.0649, Max-Change = 0.1457, Max-Change = 0.0507, Max-Change = 0.0572, Max-Change = 0.1275, Max-Change = 0.0559, Max-Change = 0.0260, Max-Change = 0.1585, Max-Change = 0.0572, Max-Change = 0.1437, Max-Change = 0.0925, Max-Change = 0.0352, Max-Change = 0.0378, Max-Change = 0.1587, Max-Change = 0.0881, Max-Change = 0.1472, Max-Change = 0.0677, Max-Change = 0.0775, Max-Change = 0.1794, Max-Change = 0.0888, Max-Change = 0.1212, Max-Change = 0.0466, Max-Change = 0.0483, Max-Change = 0.0559, Max-Change = 0.0440, Max-Change = 0.1637, Max-Change = 0.0334, Max-Change = 0.0264, Max-Change = 0.0531, Max-Change = 0.0262, Max-Change = 0.1146, Max-Change = 0.0283, Max-Change = 0.0489, Max-Change = 0.1149, Max-Change = 0.1332, Max-Change = 0.1170, Max-Change = 0.0390, Max-Change = 0.0484, Max-Change = 0.0567, Max-Change = 0.0678, Max-Change = 0.0505, Max-Change = 0.0526, Max-Change = 0.2000, Max-Change = 0.0915, Max-Change = 0.0471, Max-Change = 0.1616, Max-Change = 0.0679, Max-Change = 0.1139, Max-Change = 0.0587, Max-Change = 0.0723, Max-Change = 0.1630, Max-Change = 0.1455, Max-Change = 0.0390, Max-Change = 0.1358, Max-Change = 0.1216, Max-Change = 0.1024, Max-Change = 0.0847, Max-Change = 0.1856, Max-Change = 0.1123, Max-Change = 0.0340, Max-Change = 0.1251, Max-Change = 0.0451, Max-Change = 0.0858, Max-Change = 0.1119, Max-Change = 0.0819, Max-Change = 0.1402, Max-Change = 0.0619, Max-Change = 0.0700, Max-Change = 0.0546, Max-Change = 0.0847, Max-Change = 0.0283, Max-Change = 0.0607, Max-Change = 0.1731, Max-Change = 0.0976, Max-Change = 0.0916, Max-Change = 0.1646, Max-Change = 0.0688, Max-Change = 0.0275, Max-Change = 0.2000, Max-Change = 0.0603, Max-Change = 0.0934, Max-Change = 0.0944, Max-Change = 0.1044, Max-Change = 0.1739, Max-Change = 0.0828, Max-Change = 0.0794, Max-Change = 0.0359, Max-Change = 0.0265, Max-Change = 0.1483, Max-Change = 0.0836, Max-Change = 0.2000, Max-Change = 0.0339, Max-Change = 0.0317, Max-Change = 0.0652, Max-Change = 0.0273, Max-Change = 0.1724, Max-Change = 0.0454, Max-Change = 0.0369, Max-Change = 0.0521, Max-Change = 0.0760, Max-Change = 0.0256, Max-Change = 0.1141, Max-Change = 0.0382, Max-Change = 0.0813, Max-Change = 0.1710, Max-Change = 0.0684, Max-Change = 0.0900, Max-Change = 0.2000, Max-Change = 0.0473, Max-Change = 0.0458, Max-Change = 0.1378, Max-Change = 0.0302, Max-Change = 0.0392, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0665, gam = 0.1057, Max-Change = 0.0327, gam = 0.0780, Max-Change = 0.0223, gam = 0.0629, Max-Change = 0.0178, gam = 0.0532, Max-Change = 0.0253, gam = 0.0464, Max-Change = 0.0294, gam = 0.0413, Max-Change = 0.0359, gam = 0.0374, Max-Change = 0.0188, gam = 0.0342, Max-Change = 0.0164, gam = 0.0316, Max-Change = 0.0080, gam = 0.0294, Max-Change = 0.0069, gam = 0.0276, Max-Change = 0.0225, gam = 0.0260, Max-Change = 0.0079, gam = 0.0246, Max-Change = 0.0040, gam = 0.0233, Max-Change = 0.0105, gam = 0.0222, Max-Change = 0.0049, gam = 0.0212, Max-Change = 0.0051, gam = 0.0203, Max-Change = 0.0124, gam = 0.0195, Max-Change = 0.0096, gam = 0.0188, Max-Change = 0.0145, gam = 0.0181, Max-Change = 0.0066, gam = 0.0175, Max-Change = 0.0234, gam = 0.0169, Max-Change = 0.0129, gam = 0.0164, Max-Change = 0.0096, gam = 0.0159, Max-Change = 0.0039, gam = 0.0154, Max-Change = 0.0066, gam = 0.0150, Max-Change = 0.0084, gam = 0.0146, Max-Change = 0.0173, gam = 0.0142, Max-Change = 0.0095, gam = 0.0139, Max-Change = 0.0080, gam = 0.0135, Max-Change = 0.0065, gam = 0.0132, Max-Change = 0.0055, gam = 0.0129, Max-Change = 0.0156, gam = 0.0126, Max-Change = 0.0010, gam = 0.0124, Max-Change = 0.0063, gam = 0.0121, Max-Change = 0.0083, gam = 0.0119, Max-Change = 0.0064, gam = 0.0116, Max-Change = 0.0040, gam = 0.0114, Max-Change = 0.0010, gam = 0.0112, Max-Change = 0.0047, gam = 0.0110, Max-Change = 0.0012, gam = 0.0108, Max-Change = 0.0062, gam = 0.0106, Max-Change = 0.0039, gam = 0.0104, Max-Change = 0.0026, gam = 0.0102, Max-Change = 0.0022, gam = 0.0101, Max-Change = 0.0069, gam = 0.0099, Max-Change = 0.0022, gam = 0.0098, Max-Change = 0.0016, gam = 0.0096, Max-Change = 0.0025, gam = 0.0095, Max-Change = 0.0068, gam = 0.0093, Max-Change = 0.0016, gam = 0.0092, Max-Change = 0.0033, gam = 0.0091, Max-Change = 0.0024, gam = 0.0089, Max-Change = 0.0008, gam = 0.0088, Max-Change = 0.0018, gam = 0.0087, Max-Change = 0.0037, gam = 0.0086, Max-Change = 0.0033, gam = 0.0085, Max-Change = 0.0043, gam = 0.0084, Max-Change = 0.0015, gam = 0.0082, Max-Change = 0.0033, gam = 0.0081, Max-Change = 0.0017, gam = 0.0080, Max-Change = 0.0015, gam = 0.0080, Max-Change = 0.0023, gam = 0.0079, Max-Change = 0.0045, gam = 0.0078, Max-Change = 0.0084, gam = 0.0077, Max-Change = 0.0039, gam = 0.0076, Max-Change = 0.0057, gam = 0.0075, Max-Change = 0.0074, gam = 0.0074, Max-Change = 0.0028, gam = 0.0073, Max-Change = 0.0050, gam = 0.0073, Max-Change = 0.0041, gam = 0.0072, Max-Change = 0.0015, gam = 0.0071, Max-Change = 0.0017, gam = 0.0070, Max-Change = 0.0032, gam = 0.0070, Max-Change = 0.0028, gam = 0.0069, Max-Change = 0.0059, gam = 0.0068, Max-Change = 0.0013, gam = 0.0068, Max-Change = 0.0025, gam = 0.0067, Max-Change = 0.0022, gam = 0.0066, Max-Change = 0.0029, gam = 0.0066, Max-Change = 0.0010, gam = 0.0065, Max-Change = 0.0021, gam = 0.0065, Max-Change = 0.0022, gam = 0.0064, Max-Change = 0.0038, gam = 0.0064, Max-Change = 0.0062, gam = 0.0063, Max-Change = 0.0027, gam = 0.0062, Max-Change = 0.0034, gam = 0.0062, Max-Change = 0.0016, gam = 0.0061, Max-Change = 0.0009, gam = 0.0061, Max-Change = 0.0035, gam = 0.0060, Max-Change = 0.0053, gam = 0.0060, Max-Change = 0.0045, gam = 0.0059, Max-Change = 0.0027, gam = 0.0059, Max-Change = 0.0033, gam = 0.0058, Max-Change = 0.0043, gam = 0.0058, Max-Change = 0.0009, gam = 0.0058, Max-Change = 0.0033, gam = 0.0057, Max-Change = 0.0026, gam = 0.0057, Max-Change = 0.0025, gam = 0.0056, Max-Change = 0.0036, gam = 0.0056, Max-Change = 0.0012, gam = 0.0055, Max-Change = 0.0011, gam = 0.0055, Max-Change = 0.0044, gam = 0.0055, Max-Change = 0.0023, gam = 0.0054, Max-Change = 0.0014, gam = 0.0054, Max-Change = 0.0018, gam = 0.0053, Max-Change = 0.0028, gam = 0.0053, Max-Change = 0.0013, gam = 0.0053, Max-Change = 0.0060, gam = 0.0052, Max-Change = 0.0042, gam = 0.0052, Max-Change = 0.0029, gam = 0.0052, Max-Change = 0.0013, gam = 0.0051, Max-Change = 0.0012, gam = 0.0051, Max-Change = 0.0016, gam = 0.0051, Max-Change = 0.0022, gam = 0.0050, Max-Change = 0.0018, gam = 0.0050, Max-Change = 0.0008, gam = 0.0050, Max-Change = 0.0027, gam = 0.0049, Max-Change = 0.0021, gam = 0.0049, Max-Change = 0.0042, gam = 0.0049, Max-Change = 0.0008, gam = 0.0048, Max-Change = 0.0034, gam = 0.0048, Max-Change = 0.0020, gam = 0.0048, Max-Change = 0.0019, gam = 0.0048, Max-Change = 0.0020, gam = 0.0047, Max-Change = 0.0009, gam = 0.0047, Max-Change = 0.0027, gam = 0.0047, Max-Change = 0.0064, gam = 0.0046, Max-Change = 0.0032, gam = 0.0046, Max-Change = 0.0016, gam = 0.0046, Max-Change = 0.0038, gam = 0.0046, Max-Change = 0.0044, gam = 0.0045, Max-Change = 0.0015, gam = 0.0045, Max-Change = 0.0018, gam = 0.0045, Max-Change = 0.0019, gam = 0.0045, Max-Change = 0.0023, gam = 0.0044, Max-Change = 0.0010, gam = 0.0044, Max-Change = 0.0006, gam = 0.0044, Max-Change = 0.0012, gam = 0.0044, Max-Change = 0.0014, gam = 0.0043, Max-Change = 0.0025, gam = 0.0043, Max-Change = 0.0044, gam = 0.0043, Max-Change = 0.0040, gam = 0.0043, Max-Change = 0.0023, gam = 0.0043, Max-Change = 0.0011, gam = 0.0042, Max-Change = 0.0013, gam = 0.0042, Max-Change = 0.0020, gam = 0.0042, Max-Change = 0.0012, gam = 0.0042, Max-Change = 0.0027, gam = 0.0041, Max-Change = 0.0016, gam = 0.0041, Max-Change = 0.0035, gam = 0.0041, Max-Change = 0.0036, gam = 0.0041, Max-Change = 0.0044, gam = 0.0041, Max-Change = 0.0009, gam = 0.0040, Max-Change = 0.0007, gam = 0.0040, Max-Change = 0.0028, gam = 0.0040, Max-Change = 0.0029, gam = 0.0040, Max-Change = 0.0011, gam = 0.0040, Max-Change = 0.0026, gam = 0.0040, Max-Change = 0.0008, gam = 0.0039, Max-Change = 0.0011, gam = 0.0039, Max-Change = 0.0011, gam = 0.0039, Max-Change = 0.0017, gam = 0.0039, Max-Change = 0.0028, gam = 0.0039, Max-Change = 0.0012, gam = 0.0038, Max-Change = 0.0014, gam = 0.0038, Max-Change = 0.0020, gam = 0.0038, Max-Change = 0.0017, gam = 0.0038, Max-Change = 0.0020, gam = 0.0038, Max-Change = 0.0010, gam = 0.0038, Max-Change = 0.0028, gam = 0.0037, Max-Change = 0.0010, gam = 0.0037, Max-Change = 0.0015, gam = 0.0037, Max-Change = 0.0007, gam = 0.0037, Max-Change = 0.0008, gam = 0.0037, Max-Change = 0.0021, gam = 0.0037, Max-Change = 0.0017, gam = 0.0036, Max-Change = 0.0007, gam = 0.0036, Max-Change = 0.0020, gam = 0.0036, Max-Change = 0.0024, gam = 0.0036, Max-Change = 0.0015, gam = 0.0036, Max-Change = 0.0030, gam = 0.0036, Max-Change = 0.0018, gam = 0.0036, Max-Change = 0.0007, gam = 0.0035, Max-Change = 0.0012, gam = 0.0035, Max-Change = 0.0009, gam = 0.0035, Max-Change = 0.0007, gam = 0.0035, Max-Change = 0.0021, gam = 0.0035, Max-Change = 0.0040, gam = 0.0035, Max-Change = 0.0013, gam = 0.0035, Max-Change = 0.0032, gam = 0.0034, Max-Change = 0.0028, gam = 0.0034, Max-Change = 0.0010, gam = 0.0034, Max-Change = 0.0015, gam = 0.0034, Max-Change = 0.0009, gam = 0.0034, Max-Change = 0.0017, gam = 0.0034, Max-Change = 0.0024, gam = 0.0034, Max-Change = 0.0009, gam = 0.0034, Max-Change = 0.0003, gam = 0.0033, Max-Change = 0.0013, gam = 0.0033, Max-Change = 0.0014, gam = 0.0033, Max-Change = 0.0006, gam = 0.0033, Max-Change = 0.0013, gam = 0.0033, Max-Change = 0.0013, gam = 0.0033, Max-Change = 0.0011, gam = 0.0033, Max-Change = 0.0005, gam = 0.0033, Max-Change = 0.0046, gam = 0.0032, Max-Change = 0.0016, gam = 0.0032, Max-Change = 0.0018, gam = 0.0032, Max-Change = 0.0004, gam = 0.0032, Max-Change = 0.0020, gam = 0.0032, Max-Change = 0.0020, gam = 0.0032, Max-Change = 0.0007, gam = 0.0032, Max-Change = 0.0012, gam = 0.0032, Max-Change = 0.0015, gam = 0.0032, Max-Change = 0.0020, gam = 0.0031, Max-Change = 0.0005, gam = 0.0031, Max-Change = 0.0015, gam = 0.0031, Max-Change = 0.0017, gam = 0.0031, Max-Change = 0.0010, gam = 0.0031, Max-Change = 0.0008, gam = 0.0031, Max-Change = 0.0018, gam = 0.0031, Max-Change = 0.0024, gam = 0.0031, Max-Change = 0.0031, gam = 0.0031, Max-Change = 0.0023, gam = 0.0031, Max-Change = 0.0031, gam = 0.0030, Max-Change = 0.0004, gam = 0.0030, Max-Change = 0.0004, gam = 0.0030, Max-Change = 0.0009
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod2)
#> 
#> Call:
#> mixedmirt(data = Science, covdata = covdat, model = 1, fixed = ~0 + 
#>     group, itemtype = "gpcm")
#> 
#> --------------
#> FIXED EFFECTS:
#>        Estimate Std.Error z.value
#> groupm   -0.176     0.118  -1.494
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>    F1
#> F1  1
#> 
anova(mod, mod2)
#>           AIC    SABIC       HQ      BIC    logLik   X2 df     p
#> mod  3264.979 3276.155 3287.014 3320.577 -1618.489              
#> mod2 3258.379 3271.950 3285.136 3325.891 -1612.190 12.6  3 0.006

# graded model
mod3 <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group,
                 itemtype = 'graded')
#> , Max-Change = 0.0936, Max-Change = 0.0628, Max-Change = 0.0682, Max-Change = 0.0779, Max-Change = 0.0530, Max-Change = 0.0384, Max-Change = 0.0643, Max-Change = 0.0358, Max-Change = 0.1088, Max-Change = 0.0556, Max-Change = 0.0754, Max-Change = 0.0696, Max-Change = 0.0500, Max-Change = 0.0316, Max-Change = 0.0371, Max-Change = 0.0448, Max-Change = 0.0569, Max-Change = 0.0753, Max-Change = 0.0459, Max-Change = 0.0565, Max-Change = 0.0385, Max-Change = 0.0177, Max-Change = 0.0258, Max-Change = 0.0425, Max-Change = 0.0579, Max-Change = 0.0233, Max-Change = 0.0269, Max-Change = 0.0206, Max-Change = 0.0494, Max-Change = 0.0379, Max-Change = 0.0439, Max-Change = 0.0775, Max-Change = 0.0386, Max-Change = 0.0435, Max-Change = 0.0212, Max-Change = 0.0944, Max-Change = 0.0190, Max-Change = 0.0412, Max-Change = 0.0670, Max-Change = 0.0318, Max-Change = 0.0405, Max-Change = 0.0383, Max-Change = 0.0675, Max-Change = 0.0429, Max-Change = 0.0710, Max-Change = 0.0385, Max-Change = 0.0516, Max-Change = 0.0490, Max-Change = 0.0381, Max-Change = 0.0909, Max-Change = 0.0201, Max-Change = 0.0393, Max-Change = 0.0424, Max-Change = 0.0310, Max-Change = 0.0379, Max-Change = 0.0371, Max-Change = 0.1221, Max-Change = 0.0474, Max-Change = 0.0486, Max-Change = 0.0651, Max-Change = 0.0363, Max-Change = 0.0310, Max-Change = 0.0917, Max-Change = 0.0539, Max-Change = 0.0282, Max-Change = 0.0437, Max-Change = 0.0413, Max-Change = 0.0446, Max-Change = 0.0343, Max-Change = 0.0360, Max-Change = 0.0701, Max-Change = 0.0552, Max-Change = 0.0295, Max-Change = 0.0505, Max-Change = 0.0599, Max-Change = 0.0442, Max-Change = 0.0321, Max-Change = 0.0408, Max-Change = 0.0260, Max-Change = 0.0358, Max-Change = 0.0379, Max-Change = 0.0366, Max-Change = 0.0425, Max-Change = 0.0911, Max-Change = 0.0482, Max-Change = 0.0729, Max-Change = 0.0495, Max-Change = 0.0480, Max-Change = 0.0313, Max-Change = 0.0560, Max-Change = 0.0323, Max-Change = 0.0339, Max-Change = 0.0423, Max-Change = 0.0947, Max-Change = 0.0268, Max-Change = 0.0233, Max-Change = 0.0482, Max-Change = 0.0317, Max-Change = 0.0484, Max-Change = 0.0427, Max-Change = 0.0749, Max-Change = 0.0558, Max-Change = 0.0468, Max-Change = 0.0896, Max-Change = 0.0333, Max-Change = 0.0351, Max-Change = 0.0371, Max-Change = 0.0444, Max-Change = 0.0392, Max-Change = 0.0493, Max-Change = 0.0454, Max-Change = 0.0264, Max-Change = 0.0272, Max-Change = 0.0301, Max-Change = 0.0299, Max-Change = 0.0344, Max-Change = 0.0460, Max-Change = 0.0483, Max-Change = 0.0884, Max-Change = 0.0503, Max-Change = 0.0337, Max-Change = 0.0218, Max-Change = 0.0212, Max-Change = 0.0804, Max-Change = 0.0338, Max-Change = 0.0272, Max-Change = 0.0332, Max-Change = 0.0395, Max-Change = 0.0609, Max-Change = 0.0310, Max-Change = 0.0271, Max-Change = 0.0554, Max-Change = 0.0289, Max-Change = 0.0460, Max-Change = 0.0677, Max-Change = 0.0711, Max-Change = 0.0312, Max-Change = 0.0287, Max-Change = 0.0526, Max-Change = 0.0342, Max-Change = 0.0345, Max-Change = 0.0297, Max-Change = 0.0384, Max-Change = 0.0285, Max-Change = 0.0479, Max-Change = 0.0354, Max-Change = 0.1089, Max-Change = 0.0675, Max-Change = 0.0314, Max-Change = 0.0452, Max-Change = 0.1232, Max-Change = 0.0390, Max-Change = 0.0892, Max-Change = 0.0303, Max-Change = 0.0327, Max-Change = 0.1439, Max-Change = 0.0321, Max-Change = 0.0940, Max-Change = 0.0379, Max-Change = 0.0813, Max-Change = 0.0346, Max-Change = 0.0404, Max-Change = 0.0554, Max-Change = 0.0336, Max-Change = 0.0281, Max-Change = 0.0589, Max-Change = 0.0316, Max-Change = 0.0912, Max-Change = 0.0500, Max-Change = 0.0364, Max-Change = 0.1087, Max-Change = 0.0267, Max-Change = 0.0234, Max-Change = 0.0413, Max-Change = 0.0177, Max-Change = 0.0276, Max-Change = 0.0715, Max-Change = 0.0415, Max-Change = 0.0250, Max-Change = 0.1006, Max-Change = 0.0210, Max-Change = 0.0568, Max-Change = 0.1370, Max-Change = 0.0315, Max-Change = 0.0599, Max-Change = 0.0186, Max-Change = 0.0953, Max-Change = 0.0973, Max-Change = 0.0536, Max-Change = 0.0967, Max-Change = 0.0323, Max-Change = 0.0438, Max-Change = 0.0825, Max-Change = 0.0453, Max-Change = 0.0471, Max-Change = 0.0624, Max-Change = 0.0290, Max-Change = 0.0384, Max-Change = 0.0372, Max-Change = 0.0810, Max-Change = 0.0457, Max-Change = 0.0415, Max-Change = 0.0329, Max-Change = 0.0675, Max-Change = 0.0283, Max-Change = 0.0571, Max-Change = 0.0277, Max-Change = 0.0567, Max-Change = 0.0597, Max-Change = 0.1081, Max-Change = 0.0365, Max-Change = 0.0409, Max-Change = 0.0624, Max-Change = 0.0567, Max-Change = 0.0347, Max-Change = 0.1347, Max-Change = 0.0570, Max-Change = 0.0683, Max-Change = 0.0604, Max-Change = 0.0916, Max-Change = 0.0895, Max-Change = 0.0384, Max-Change = 0.0310, Max-Change = 0.0302, Max-Change = 0.0528, Max-Change = 0.0244, Max-Change = 0.0991, Max-Change = 0.0630, Max-Change = 0.0419, Max-Change = 0.0702, Max-Change = 0.0346, Max-Change = 0.0397, Max-Change = 0.0542, Max-Change = 0.0295, Max-Change = 0.0251, Max-Change = 0.0376, Max-Change = 0.0436, Max-Change = 0.0139, Max-Change = 0.0323, Max-Change = 0.0293, Max-Change = 0.0452, Max-Change = 0.1383, Max-Change = 0.0361, Max-Change = 0.0899, Max-Change = 0.0155, Max-Change = 0.0330, Max-Change = 0.0329, Max-Change = 0.0665, Max-Change = 0.0744, Max-Change = 0.0201, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0371, gam = 0.1057, Max-Change = 0.0347, gam = 0.0780, Max-Change = 0.0089, gam = 0.0629, Max-Change = 0.0091, gam = 0.0532, Max-Change = 0.0179, gam = 0.0464, Max-Change = 0.0041, gam = 0.0413, Max-Change = 0.0156, gam = 0.0374, Max-Change = 0.0204, gam = 0.0342, Max-Change = 0.0039, gam = 0.0316, Max-Change = 0.0069, gam = 0.0294, Max-Change = 0.0038, gam = 0.0276, Max-Change = 0.0136, gam = 0.0260, Max-Change = 0.0049, gam = 0.0246, Max-Change = 0.0066, gam = 0.0233, Max-Change = 0.0074, gam = 0.0222, Max-Change = 0.0075, gam = 0.0212, Max-Change = 0.0039, gam = 0.0203, Max-Change = 0.0053, gam = 0.0195, Max-Change = 0.0072, gam = 0.0188, Max-Change = 0.0108, gam = 0.0181, Max-Change = 0.0049, gam = 0.0175, Max-Change = 0.0108, gam = 0.0169, Max-Change = 0.0109, gam = 0.0164, Max-Change = 0.0079, gam = 0.0159, Max-Change = 0.0041, gam = 0.0154, Max-Change = 0.0025, gam = 0.0150, Max-Change = 0.0063, gam = 0.0146, Max-Change = 0.0015, gam = 0.0142, Max-Change = 0.0017, gam = 0.0139, Max-Change = 0.0084, gam = 0.0135, Max-Change = 0.0026, gam = 0.0132, Max-Change = 0.0021, gam = 0.0129, Max-Change = 0.0017, gam = 0.0126, Max-Change = 0.0022, gam = 0.0124, Max-Change = 0.0036, gam = 0.0121, Max-Change = 0.0035, gam = 0.0119, Max-Change = 0.0030, gam = 0.0116, Max-Change = 0.0023, gam = 0.0114, Max-Change = 0.0019, gam = 0.0112, Max-Change = 0.0028, gam = 0.0110, Max-Change = 0.0014, gam = 0.0108, Max-Change = 0.0045, gam = 0.0106, Max-Change = 0.0017, gam = 0.0104, Max-Change = 0.0016, gam = 0.0102, Max-Change = 0.0021, gam = 0.0101, Max-Change = 0.0019, gam = 0.0099, Max-Change = 0.0030, gam = 0.0098, Max-Change = 0.0021, gam = 0.0096, Max-Change = 0.0016, gam = 0.0095, Max-Change = 0.0009, gam = 0.0093, Max-Change = 0.0015, gam = 0.0092, Max-Change = 0.0024, gam = 0.0091, Max-Change = 0.0024, gam = 0.0089, Max-Change = 0.0015, gam = 0.0088, Max-Change = 0.0032, gam = 0.0087, Max-Change = 0.0023, gam = 0.0086, Max-Change = 0.0046, gam = 0.0085, Max-Change = 0.0013, gam = 0.0084, Max-Change = 0.0015, gam = 0.0082, Max-Change = 0.0019, gam = 0.0081, Max-Change = 0.0014, gam = 0.0080, Max-Change = 0.0016, gam = 0.0080, Max-Change = 0.0030, gam = 0.0079, Max-Change = 0.0009, gam = 0.0078, Max-Change = 0.0022, gam = 0.0077, Max-Change = 0.0009, gam = 0.0076, Max-Change = 0.0027, gam = 0.0075, Max-Change = 0.0026, gam = 0.0074, Max-Change = 0.0018, gam = 0.0073, Max-Change = 0.0044, gam = 0.0073, Max-Change = 0.0013, gam = 0.0072, Max-Change = 0.0018, gam = 0.0071, Max-Change = 0.0013, gam = 0.0070, Max-Change = 0.0015, gam = 0.0070, Max-Change = 0.0030, gam = 0.0069, Max-Change = 0.0032, gam = 0.0068, Max-Change = 0.0006, gam = 0.0068, Max-Change = 0.0016, gam = 0.0067, Max-Change = 0.0018, gam = 0.0066, Max-Change = 0.0017, gam = 0.0066, Max-Change = 0.0020, gam = 0.0065, Max-Change = 0.0007, gam = 0.0065, Max-Change = 0.0008, gam = 0.0064, Max-Change = 0.0013, gam = 0.0064, Max-Change = 0.0023, gam = 0.0063, Max-Change = 0.0024, gam = 0.0062, Max-Change = 0.0020, gam = 0.0062, Max-Change = 0.0009, gam = 0.0061, Max-Change = 0.0008, gam = 0.0061, Max-Change = 0.0016, gam = 0.0060, Max-Change = 0.0011, gam = 0.0060, Max-Change = 0.0005, gam = 0.0059, Max-Change = 0.0013, gam = 0.0059, Max-Change = 0.0007, gam = 0.0058, Max-Change = 0.0009, gam = 0.0058, Max-Change = 0.0007
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
coef(mod3)
#> $Comfort
#>         groupm    a1    d1    d2     d3
#> par     -0.248 1.000 4.951 2.735 -1.327
#> CI_2.5  -0.578 0.636 3.989 2.289 -1.682
#> CI_97.5  0.082 1.363 5.913 3.181 -0.971
#> 
#> $Work
#>         groupm    a1    d1    d2     d3
#> par     -0.248 1.212 3.034 1.020 -2.133
#> CI_2.5  -0.578 0.873 2.545 0.697 -2.554
#> CI_97.5  0.082 1.551 3.522 1.342 -1.712
#> 
#> $Future
#>         groupm    a1    d1    d2     d3
#> par     -0.248 2.581 5.761 2.524 -1.995
#> CI_2.5  -0.578 1.251 3.733 1.497 -2.807
#> CI_97.5  0.082 3.911 7.789 3.551 -1.182
#> 
#> $Benefit
#>         groupm    a1    d1    d2     d3
#> par     -0.248 1.075 3.456 1.110 -1.554
#> CI_2.5  -0.578 0.728 2.905 0.796 -1.925
#> CI_97.5  0.082 1.423 4.007 1.423 -1.183
#> 
#> $GroupPars
#>         MEAN_1 COV_11
#> par          0      1
#> CI_2.5      NA     NA
#> CI_97.5     NA     NA
#> 


###################################################
# latent regression with Rasch and 2PL models

set.seed(1)
n <- 300
a <- matrix(1, 10)
d <- matrix(rnorm(10))
Theta <- matrix(c(rnorm(n, 0), rnorm(n, 1), rnorm(n, 2)))
covdata <- data.frame(group=rep(c('g1','g2','g3'), each=n))
dat <- simdata(a, d, N=n*3, Theta=Theta, itemtype = '2PL')
itemstats(dat)
#> $overall
#>    N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  900            6.932          2.347 0.197 0.035 0.709     1.266
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  900 2 0.570 0.495   0.604         0.443       0.673
#> Item_2  900 2 0.689 0.463   0.518         0.351       0.690
#> Item_3  900 2 0.526 0.500   0.531         0.352       0.691
#> Item_4  900 2 0.892 0.310   0.422         0.305       0.698
#> Item_5  900 2 0.748 0.435   0.522         0.367       0.687
#> Item_6  900 2 0.543 0.498   0.569         0.398       0.682
#> Item_7  900 2 0.759 0.428   0.530         0.379       0.685
#> Item_8  900 2 0.803 0.398   0.516         0.375       0.686
#> Item_9  900 2 0.777 0.417   0.489         0.337       0.692
#> Item_10 900 2 0.626 0.484   0.548         0.378       0.685
#> 
#> $proportions
#>             0     1
#> Item_1  0.430 0.570
#> Item_2  0.311 0.689
#> Item_3  0.474 0.526
#> Item_4  0.108 0.892
#> Item_5  0.252 0.748
#> Item_6  0.457 0.543
#> Item_7  0.241 0.759
#> Item_8  0.197 0.803
#> Item_9  0.223 0.777
#> Item_10 0.374 0.626
#> 

# had we known the latent abilities, we could have computed the regression coefs
summary(lm(Theta ~ covdata$group))
#> 
#> Call:
#> lm(formula = Theta ~ covdata$group)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.9871 -0.6851 -0.0427  0.7170  3.8313 
#> 
#> Coefficients:
#>                 Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)      0.04434    0.05970   0.743    0.458    
#> covdata$groupg2  0.93468    0.08443  11.071   <2e-16 ***
#> covdata$groupg3  1.88134    0.08443  22.284   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 1.034 on 897 degrees of freedom
#> Multiple R-squared:  0.3563, Adjusted R-squared:  0.3549 
#> F-statistic: 248.3 on 2 and 897 DF,  p-value: < 2.2e-16
#> 

# but all we have is observed test data. Latent regression helps to recover these coefs
# Rasch model approach (and mirt equivalent)
rmod0 <- mirt(dat, 1, 'Rasch') # unconditional

# these two models are equivalent
rmod1a <- mirt(dat, 1, 'Rasch', covdata = covdata, formula = ~ group)
rmod1b <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items + group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1660, Max-Change = 0.1608, Max-Change = 0.1153, Max-Change = 0.1284, Max-Change = 0.1084, Max-Change = 0.0957, Max-Change = 0.0911, Max-Change = 0.0686, Max-Change = 0.0641, Max-Change = 0.0574, Max-Change = 0.0470, Max-Change = 0.0606, Max-Change = 0.0515, Max-Change = 0.0547, Max-Change = 0.0557, Max-Change = 0.0498, Max-Change = 0.0395, Max-Change = 0.0418, Max-Change = 0.0242, Max-Change = 0.0275, Max-Change = 0.0258, Max-Change = 0.0452, Max-Change = 0.0276, Max-Change = 0.0125, Max-Change = 0.0253, Max-Change = 0.0108, Max-Change = 0.0333, Max-Change = 0.0210, Max-Change = 0.0134, Max-Change = 0.0152, Max-Change = 0.0231, Max-Change = 0.0217, Max-Change = 0.0079, Max-Change = 0.0055, Max-Change = 0.0177, Max-Change = 0.0161, Max-Change = 0.0112, Max-Change = 0.0215, Max-Change = 0.0093, Max-Change = 0.0148, Max-Change = 0.0046, Max-Change = 0.0062, Max-Change = 0.0209, Max-Change = 0.0151, Max-Change = 0.0218, Max-Change = 0.0196, Max-Change = 0.0184, Max-Change = 0.0086, Max-Change = 0.0131, Max-Change = 0.0063, Max-Change = 0.0094, Max-Change = 0.0153, Max-Change = 0.0218, Max-Change = 0.0101, Max-Change = 0.0083, Max-Change = 0.0093, Max-Change = 0.0063, Max-Change = 0.0164, Max-Change = 0.0056, Max-Change = 0.0080, Max-Change = 0.0213, Max-Change = 0.0182, Max-Change = 0.0181, Max-Change = 0.0143, Max-Change = 0.0046, Max-Change = 0.0158, Max-Change = 0.0097, Max-Change = 0.0053, Max-Change = 0.0100, Max-Change = 0.0160, Max-Change = 0.0180, Max-Change = 0.0125, Max-Change = 0.0130, Max-Change = 0.0052, Max-Change = 0.0130, Max-Change = 0.0063, Max-Change = 0.0151, Max-Change = 0.0077, Max-Change = 0.0175, Max-Change = 0.0040, Max-Change = 0.0118, Max-Change = 0.0074, Max-Change = 0.0119, Max-Change = 0.0034, Max-Change = 0.0066, Max-Change = 0.0120, Max-Change = 0.0091, Max-Change = 0.0088, Max-Change = 0.0247, Max-Change = 0.0402, Max-Change = 0.0111, Max-Change = 0.0130, Max-Change = 0.0303, Max-Change = 0.0130, Max-Change = 0.0019, Max-Change = 0.0164, Max-Change = 0.0046, Max-Change = 0.0144, Max-Change = 0.0129, Max-Change = 0.0112, Max-Change = 0.0158, Max-Change = 0.0095, Max-Change = 0.0094, Max-Change = 0.0093, Max-Change = 0.0221, Max-Change = 0.0065, Max-Change = 0.0094, Max-Change = 0.0084, Max-Change = 0.0103, Max-Change = 0.0149, Max-Change = 0.0050, Max-Change = 0.0116, Max-Change = 0.0152, Max-Change = 0.0187, Max-Change = 0.0131, Max-Change = 0.0103, Max-Change = 0.0097, Max-Change = 0.0067, Max-Change = 0.0101, Max-Change = 0.0154, Max-Change = 0.0066, Max-Change = 0.0054, Max-Change = 0.0095, Max-Change = 0.0134, Max-Change = 0.0274, Max-Change = 0.0198, Max-Change = 0.0215, Max-Change = 0.0062, Max-Change = 0.0192, Max-Change = 0.0136, Max-Change = 0.0070, Max-Change = 0.0231, Max-Change = 0.0084, Max-Change = 0.0036, Max-Change = 0.0111, Max-Change = 0.0212, Max-Change = 0.0077, Max-Change = 0.0103, Max-Change = 0.0066, Max-Change = 0.0145, Max-Change = 0.0148, Max-Change = 0.0220, Max-Change = 0.0106, Max-Change = 0.0029, Max-Change = 0.0118, Max-Change = 0.0063, Max-Change = 0.0111, Max-Change = 0.0124, Max-Change = 0.0123, Max-Change = 0.0334, Max-Change = 0.0058, Max-Change = 0.0069, Max-Change = 0.0107, Max-Change = 0.0096, Max-Change = 0.0085, Max-Change = 0.0131, Max-Change = 0.0087, Max-Change = 0.0044, Max-Change = 0.0128, Max-Change = 0.0061, Max-Change = 0.0071, Max-Change = 0.0255, Max-Change = 0.0085, Max-Change = 0.0039, Max-Change = 0.0073, Max-Change = 0.0116, Max-Change = 0.0178, Max-Change = 0.0079, Max-Change = 0.0110, Max-Change = 0.0258, Max-Change = 0.0104, Max-Change = 0.0107, Max-Change = 0.0188, Max-Change = 0.0064, Max-Change = 0.0130, Max-Change = 0.0245, Max-Change = 0.0130, Max-Change = 0.0056, Max-Change = 0.0093, Max-Change = 0.0055, Max-Change = 0.0215, Max-Change = 0.0085, Max-Change = 0.0101, Max-Change = 0.0114, Max-Change = 0.0121, Max-Change = 0.0113, Max-Change = 0.0278, Max-Change = 0.0092, Max-Change = 0.0157, Max-Change = 0.0057, Max-Change = 0.0107, Max-Change = 0.0197, Max-Change = 0.0160, Max-Change = 0.0130, Max-Change = 0.0111, Max-Change = 0.0104, Max-Change = 0.0099, Max-Change = 0.0070, Max-Change = 0.0122, Max-Change = 0.0056, Max-Change = 0.0191, Max-Change = 0.0057, Max-Change = 0.0093, Max-Change = 0.0135, Max-Change = 0.0127, Max-Change = 0.0101, Max-Change = 0.0207, Max-Change = 0.0223, Max-Change = 0.0166, Max-Change = 0.0193, Max-Change = 0.0147, Max-Change = 0.0132, Max-Change = 0.0091, Max-Change = 0.0186, Max-Change = 0.0149, Max-Change = 0.0111, Max-Change = 0.0074, Max-Change = 0.0106, Max-Change = 0.0046, Max-Change = 0.0145, Max-Change = 0.0145, Max-Change = 0.0062, Max-Change = 0.0128, Max-Change = 0.0127, Max-Change = 0.0167, Max-Change = 0.0117, Max-Change = 0.0079, Max-Change = 0.0170, Max-Change = 0.0163, Max-Change = 0.0093, Max-Change = 0.0094, Max-Change = 0.0147, Max-Change = 0.0057, Max-Change = 0.0175, Max-Change = 0.0084, Max-Change = 0.0149, Max-Change = 0.0051, Max-Change = 0.0126, Max-Change = 0.0149, Max-Change = 0.0224, Max-Change = 0.0094, Max-Change = 0.0108, Max-Change = 0.0175, Max-Change = 0.0080, Max-Change = 0.0122, Max-Change = 0.0100, Max-Change = 0.0113, Max-Change = 0.0164, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0118, gam = 0.1057, Max-Change = 0.0087, gam = 0.0780, Max-Change = 0.0021, gam = 0.0629, Max-Change = 0.0056, gam = 0.0532, Max-Change = 0.0039, gam = 0.0464, Max-Change = 0.0022, gam = 0.0413, Max-Change = 0.0034, gam = 0.0374, Max-Change = 0.0035, gam = 0.0342, Max-Change = 0.0008, gam = 0.0316, Max-Change = 0.0007, gam = 0.0294, Max-Change = 0.0039, gam = 0.0276, Max-Change = 0.0024, gam = 0.0260, Max-Change = 0.0019, gam = 0.0246, Max-Change = 0.0009, gam = 0.0233, Max-Change = 0.0017, gam = 0.0222, Max-Change = 0.0011, gam = 0.0212, Max-Change = 0.0014, gam = 0.0203, Max-Change = 0.0008, gam = 0.0195, Max-Change = 0.0002, gam = 0.0188, Max-Change = 0.0021, gam = 0.0181, Max-Change = 0.0009, gam = 0.0175, Max-Change = 0.0015, gam = 0.0169, Max-Change = 0.0011, gam = 0.0164, Max-Change = 0.0010, gam = 0.0159, Max-Change = 0.0004, gam = 0.0154, Max-Change = 0.0009, gam = 0.0150, Max-Change = 0.0003
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
anova(rmod0, rmod1b)
#>             AIC    SABIC       HQ      BIC    logLik      X2 df p
#> rmod0  9698.483 9716.375 9718.663 9751.310 -4838.242             
#> rmod1b 9472.610 9493.755 9496.459 9535.041 -4723.305 229.873  2 0
coef(rmod1a, simplify=TRUE)
#> $items
#>         a1      d g u
#> Item_1   1 -0.466 0 1
#> Item_2   1  0.193 0 1
#> Item_3   1 -0.699 0 1
#> Item_4   1  1.795 0 1
#> Item_5   1  0.561 0 1
#> Item_6   1 -0.607 0 1
#> Item_7   1  0.636 0 1
#> Item_8   1  0.957 0 1
#> Item_9   1  0.759 0 1
#> Item_10  1 -0.167 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>       F1
#> F1 1.005
#> 
#> $lr.betas
#>                F1
#> (Intercept) 0.000
#> groupg2     0.797
#> groupg3     1.707
#> 
summary(rmod1b)
#> 
#> Call:
#> mixedmirt(data = dat, covdata = covdata, model = 1, fixed = ~0 + 
#>     items + group)
#> 
#> --------------
#> FIXED EFFECTS:
#>         Estimate Std.Error z.value
#> groupg2    0.821     0.118   6.972
#> groupg3    1.683     0.111  15.140
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>    F1
#> F1  1
#> 

# 2PL, requires different input to allow Theta variance to remain fixed
mod0 <- mirt(dat, 1) # unconditional
mod1a <- mirt(dat, 1, covdata = covdata, formula = ~ group, itemtype = '2PL')
mod1b <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items, lr.fixed = ~group, itemtype = '2PL')
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1937, Max-Change = 0.1590, Max-Change = 0.1455, Max-Change = 0.1217, Max-Change = 0.1043, Max-Change = 0.0948, Max-Change = 0.0799, Max-Change = 0.0691, Max-Change = 0.0633, Max-Change = 0.0481, Max-Change = 0.0417, Max-Change = 0.0393, Max-Change = 0.0264, Max-Change = 0.0549, Max-Change = 0.0198, Max-Change = 0.0260, Max-Change = 0.0205, Max-Change = 0.0201, Max-Change = 0.0274, Max-Change = 0.0250, Max-Change = 0.0381, Max-Change = 0.0146, Max-Change = 0.0227, Max-Change = 0.0456, Max-Change = 0.0278, Max-Change = 0.0151, Max-Change = 0.0267, Max-Change = 0.0274, Max-Change = 0.0234, Max-Change = 0.0229, Max-Change = 0.0192, Max-Change = 0.0254, Max-Change = 0.0336, Max-Change = 0.0151, Max-Change = 0.0182, Max-Change = 0.0177, Max-Change = 0.0190, Max-Change = 0.0151, Max-Change = 0.0118, Max-Change = 0.0229, Max-Change = 0.0186, Max-Change = 0.0140, Max-Change = 0.0197, Max-Change = 0.0176, Max-Change = 0.0128, Max-Change = 0.0154, Max-Change = 0.0254, Max-Change = 0.0187, Max-Change = 0.0143, Max-Change = 0.0170, Max-Change = 0.0234, Max-Change = 0.0362, Max-Change = 0.0341, Max-Change = 0.0320, Max-Change = 0.0249, Max-Change = 0.0158, Max-Change = 0.0274, Max-Change = 0.0170, Max-Change = 0.0274, Max-Change = 0.0203, Max-Change = 0.0240, Max-Change = 0.0293, Max-Change = 0.0185, Max-Change = 0.0166, Max-Change = 0.0328, Max-Change = 0.0220, Max-Change = 0.0272, Max-Change = 0.0278, Max-Change = 0.0197, Max-Change = 0.0209, Max-Change = 0.0114, Max-Change = 0.0254, Max-Change = 0.0274, Max-Change = 0.0209, Max-Change = 0.0179, Max-Change = 0.0176, Max-Change = 0.0122, Max-Change = 0.0209, Max-Change = 0.0207, Max-Change = 0.0190, Max-Change = 0.0240, Max-Change = 0.0207, Max-Change = 0.0182, Max-Change = 0.0190, Max-Change = 0.0265, Max-Change = 0.0289, Max-Change = 0.0183, Max-Change = 0.0164, Max-Change = 0.0246, Max-Change = 0.0345, Max-Change = 0.0196, Max-Change = 0.0266, Max-Change = 0.0163, Max-Change = 0.0237, Max-Change = 0.0174, Max-Change = 0.0131, Max-Change = 0.0152, Max-Change = 0.0165, Max-Change = 0.0232, Max-Change = 0.0191, Max-Change = 0.0154, Max-Change = 0.0209, Max-Change = 0.0144, Max-Change = 0.0282, Max-Change = 0.0247, Max-Change = 0.0117, Max-Change = 0.0103, Max-Change = 0.0224, Max-Change = 0.0156, Max-Change = 0.0179, Max-Change = 0.0159, Max-Change = 0.0229, Max-Change = 0.0217, Max-Change = 0.0193, Max-Change = 0.0191, Max-Change = 0.0189, Max-Change = 0.0414, Max-Change = 0.0325, Max-Change = 0.0230, Max-Change = 0.0280, Max-Change = 0.0241, Max-Change = 0.0183, Max-Change = 0.0318, Max-Change = 0.0178, Max-Change = 0.0209, Max-Change = 0.0166, Max-Change = 0.0181, Max-Change = 0.0169, Max-Change = 0.0252, Max-Change = 0.0146, Max-Change = 0.0267, Max-Change = 0.0209, Max-Change = 0.0270, Max-Change = 0.0175, Max-Change = 0.0263, Max-Change = 0.0149, Max-Change = 0.0140, Max-Change = 0.0228, Max-Change = 0.0180, Max-Change = 0.0383, Max-Change = 0.0207, Max-Change = 0.0251, Max-Change = 0.0118, Max-Change = 0.0148, Max-Change = 0.0181, Max-Change = 0.0310, Max-Change = 0.0256, Max-Change = 0.0140, Max-Change = 0.0130, Max-Change = 0.0195, Max-Change = 0.0199, Max-Change = 0.0365, Max-Change = 0.0221, Max-Change = 0.0294, Max-Change = 0.0219, Max-Change = 0.0207, Max-Change = 0.0159, Max-Change = 0.0211, Max-Change = 0.0287, Max-Change = 0.0139, Max-Change = 0.0225, Max-Change = 0.0134, Max-Change = 0.0103, Max-Change = 0.0168, Max-Change = 0.0301, Max-Change = 0.0250, Max-Change = 0.0232, Max-Change = 0.0402, Max-Change = 0.0227, Max-Change = 0.0194, Max-Change = 0.0263, Max-Change = 0.0164, Max-Change = 0.0164, Max-Change = 0.0213, Max-Change = 0.0324, Max-Change = 0.0187, Max-Change = 0.0201, Max-Change = 0.0249, Max-Change = 0.0327, Max-Change = 0.0269, Max-Change = 0.0392, Max-Change = 0.0311, Max-Change = 0.0137, Max-Change = 0.0276, Max-Change = 0.0189, Max-Change = 0.0218, Max-Change = 0.0157, Max-Change = 0.0336, Max-Change = 0.0302, Max-Change = 0.0216, Max-Change = 0.0214, Max-Change = 0.0190, Max-Change = 0.0213, Max-Change = 0.0266, Max-Change = 0.0149, Max-Change = 0.0213, Max-Change = 0.0228, Max-Change = 0.0215, Max-Change = 0.0333, Max-Change = 0.0090, Max-Change = 0.0261, Max-Change = 0.0156, Max-Change = 0.0148, Max-Change = 0.0208, Max-Change = 0.0142, Max-Change = 0.0172, Max-Change = 0.0320, Max-Change = 0.0210, Max-Change = 0.0249, Max-Change = 0.0289, Max-Change = 0.0183, Max-Change = 0.0131, Max-Change = 0.0257, Max-Change = 0.0169, Max-Change = 0.0257, Max-Change = 0.0114, Max-Change = 0.0163, Max-Change = 0.0158, Max-Change = 0.0116, Max-Change = 0.0185, Max-Change = 0.0238, Max-Change = 0.0268, Max-Change = 0.0185, Max-Change = 0.0219, Max-Change = 0.0204, Max-Change = 0.0216, Max-Change = 0.0216, Max-Change = 0.0324, Max-Change = 0.0184, Max-Change = 0.0172, Max-Change = 0.0248, Max-Change = 0.0177, Max-Change = 0.0218, Max-Change = 0.0126, Max-Change = 0.0175, Max-Change = 0.0239, Max-Change = 0.0197, Max-Change = 0.0267, Max-Change = 0.0247, Max-Change = 0.0338, Max-Change = 0.0156, Max-Change = 0.0186, Max-Change = 0.0153, Max-Change = 0.0425, Max-Change = 0.0265, Max-Change = 0.0145, Max-Change = 0.0254, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0225, gam = 0.1057, Max-Change = 0.0134, gam = 0.0780, Max-Change = 0.0108, gam = 0.0629, Max-Change = 0.0051, gam = 0.0532, Max-Change = 0.0074, gam = 0.0464, Max-Change = 0.0045, gam = 0.0413, Max-Change = 0.0046, gam = 0.0374, Max-Change = 0.0044, gam = 0.0342, Max-Change = 0.0056, gam = 0.0316, Max-Change = 0.0034, gam = 0.0294, Max-Change = 0.0035, gam = 0.0276, Max-Change = 0.0032, gam = 0.0260, Max-Change = 0.0024, gam = 0.0246, Max-Change = 0.0036, gam = 0.0233, Max-Change = 0.0029, gam = 0.0222, Max-Change = 0.0030, gam = 0.0212, Max-Change = 0.0034, gam = 0.0203, Max-Change = 0.0021, gam = 0.0195, Max-Change = 0.0019, gam = 0.0188, Max-Change = 0.0016, gam = 0.0181, Max-Change = 0.0019, gam = 0.0175, Max-Change = 0.0018, gam = 0.0169, Max-Change = 0.0011, gam = 0.0164, Max-Change = 0.0027, gam = 0.0159, Max-Change = 0.0013, gam = 0.0154, Max-Change = 0.0008, gam = 0.0150, Max-Change = 0.0019, gam = 0.0146, Max-Change = 0.0022, gam = 0.0142, Max-Change = 0.0011, gam = 0.0139, Max-Change = 0.0012, gam = 0.0135, Max-Change = 0.0015, gam = 0.0132, Max-Change = 0.0016, gam = 0.0129, Max-Change = 0.0012, gam = 0.0126, Max-Change = 0.0018, gam = 0.0124, Max-Change = 0.0013, gam = 0.0121, Max-Change = 0.0013, gam = 0.0119, Max-Change = 0.0008, gam = 0.0116, Max-Change = 0.0008, gam = 0.0114, Max-Change = 0.0011, gam = 0.0112, Max-Change = 0.0009, gam = 0.0110, Max-Change = 0.0022, gam = 0.0108, Max-Change = 0.0014, gam = 0.0106, Max-Change = 0.0014, gam = 0.0104, Max-Change = 0.0016, gam = 0.0102, Max-Change = 0.0011, gam = 0.0101, Max-Change = 0.0005, gam = 0.0099, Max-Change = 0.0009, gam = 0.0098, Max-Change = 0.0010
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
anova(mod0, mod1b)
#>            AIC    SABIC       HQ      BIC    logLik      X2 df p
#> mod0  9706.089 9738.620 9742.780 9802.136 -4833.044             
#> mod1b 9481.042 9516.826 9521.402 9586.694 -4718.521 229.047  2 0
coef(mod1a)$lr.betas
#>                    F1
#> (Intercept) 0.0000000
#> groupg2     0.7910307
#> groupg3     1.7064184
summary(mod1b)
#> 
#> Call:
#> mixedmirt(data = dat, covdata = covdata, model = 1, fixed = ~0 + 
#>     items, itemtype = "2PL", lr.fixed = ~group)
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>    F1
#> F1  1
#> 
#> --------------
#> LATENT REGRESSION FIXED EFFECTS:
#> 
#>                F1
#> (Intercept) 0.000
#> groupg2     0.767
#> groupg3     1.711
#> 
#>             Std.Error_F1 z_F1
#> (Intercept)           NA   NA
#> groupg2              NaN  NaN
#> groupg3              NaN  NaN

# specifying specific regression effects is accomplished by passing a list of formula
model <- 'F1 = 1-5
         F2 = 6-10'
covdata$contvar <- rnorm(nrow(covdata))
mod2 <- mirt(dat, model, itemtype = 'Rasch', covdata=covdata,
        formula = list(F1 = ~ group + contvar, F2 = ~ group))
coef(mod2)[11:12]
#> $GroupPars
#>     MEAN_1 MEAN_2   COV_11 COV_21  COV_22
#> par      0      0 1.041367      0 1.10377
#> 
#> $lr.betas
#>                     F1        F2
#> (Intercept)  0.0000000 0.0000000
#> groupg2      0.7211501 0.9034825
#> groupg3      1.7458886 1.6959637
#> contvar     -0.0220561 0.0000000
#> 
mod2b <- mixedmirt(dat, covdata, model, fixed = ~ 0 + items,
        lr.fixed = list(F1 = ~ group + contvar, F2 = ~ group))
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1871, Max-Change = 0.1644, Max-Change = 0.1436, Max-Change = 0.1201, Max-Change = 0.0966, Max-Change = 0.0836, Max-Change = 0.0668, Max-Change = 0.0564, Max-Change = 0.0473, Max-Change = 0.0478, Max-Change = 0.0398, Max-Change = 0.0320, Max-Change = 0.0277, Max-Change = 0.0222, Max-Change = 0.0242, Max-Change = 0.0219, Max-Change = 0.0384, Max-Change = 0.0196, Max-Change = 0.0110, Max-Change = 0.0102, Max-Change = 0.0227, Max-Change = 0.0374, Max-Change = 0.0230, Max-Change = 0.0220, Max-Change = 0.0179, Max-Change = 0.0126, Max-Change = 0.0147, Max-Change = 0.0112, Max-Change = 0.0132, Max-Change = 0.0187, Max-Change = 0.0215, Max-Change = 0.0175, Max-Change = 0.0136, Max-Change = 0.0183, Max-Change = 0.0191, Max-Change = 0.0203, Max-Change = 0.0139, Max-Change = 0.0104, Max-Change = 0.0128, Max-Change = 0.0171, Max-Change = 0.0076, Max-Change = 0.0106, Max-Change = 0.0171, Max-Change = 0.0120, Max-Change = 0.0079, Max-Change = 0.0159, Max-Change = 0.0150, Max-Change = 0.0168, Max-Change = 0.0141, Max-Change = 0.0118, Max-Change = 0.0125, Max-Change = 0.0125, Max-Change = 0.0175, Max-Change = 0.0294, Max-Change = 0.0151, Max-Change = 0.0190, Max-Change = 0.0217, Max-Change = 0.0161, Max-Change = 0.0126, Max-Change = 0.0155, Max-Change = 0.0153, Max-Change = 0.0110, Max-Change = 0.0178, Max-Change = 0.0240, Max-Change = 0.0150, Max-Change = 0.0126, Max-Change = 0.0265, Max-Change = 0.0196, Max-Change = 0.0156, Max-Change = 0.0130, Max-Change = 0.0129, Max-Change = 0.0152, Max-Change = 0.0131, Max-Change = 0.0238, Max-Change = 0.0202, Max-Change = 0.0158, Max-Change = 0.0153, Max-Change = 0.0320, Max-Change = 0.0169, Max-Change = 0.0065, Max-Change = 0.0112, Max-Change = 0.0151, Max-Change = 0.0143, Max-Change = 0.0103, Max-Change = 0.0204, Max-Change = 0.0275, Max-Change = 0.0153, Max-Change = 0.0142, Max-Change = 0.0090, Max-Change = 0.0178, Max-Change = 0.0169, Max-Change = 0.0168, Max-Change = 0.0179, Max-Change = 0.0116, Max-Change = 0.0154, Max-Change = 0.0159, Max-Change = 0.0073, Max-Change = 0.0198, Max-Change = 0.0088, Max-Change = 0.0131, Max-Change = 0.0153, Max-Change = 0.0103, Max-Change = 0.0208, Max-Change = 0.0164, Max-Change = 0.0185, Max-Change = 0.0169, Max-Change = 0.0157, Max-Change = 0.0116, Max-Change = 0.0152, Max-Change = 0.0185, Max-Change = 0.0235, Max-Change = 0.0216, Max-Change = 0.0143, Max-Change = 0.0147, Max-Change = 0.0131, Max-Change = 0.0133, Max-Change = 0.0147, Max-Change = 0.0245, Max-Change = 0.0225, Max-Change = 0.0147, Max-Change = 0.0129, Max-Change = 0.0187, Max-Change = 0.0161, Max-Change = 0.0104, Max-Change = 0.0095, Max-Change = 0.0215, Max-Change = 0.0172, Max-Change = 0.0230, Max-Change = 0.0131, Max-Change = 0.0124, Max-Change = 0.0204, Max-Change = 0.0158, Max-Change = 0.0079, Max-Change = 0.0159, Max-Change = 0.0189, Max-Change = 0.0205, Max-Change = 0.0081, Max-Change = 0.0107, Max-Change = 0.0105, Max-Change = 0.0112, Max-Change = 0.0151, Max-Change = 0.0155, Max-Change = 0.0216, Max-Change = 0.0065, Max-Change = 0.0172, Max-Change = 0.0146, Max-Change = 0.0218, Max-Change = 0.0178, Max-Change = 0.0161, Max-Change = 0.0150, Max-Change = 0.0168, Max-Change = 0.0229, Max-Change = 0.0081, Max-Change = 0.0251, Max-Change = 0.0164, Max-Change = 0.0177, Max-Change = 0.0141, Max-Change = 0.0165, Max-Change = 0.0094, Max-Change = 0.0180, Max-Change = 0.0202, Max-Change = 0.0258, Max-Change = 0.0112, Max-Change = 0.0084, Max-Change = 0.0145, Max-Change = 0.0118, Max-Change = 0.0186, Max-Change = 0.0136, Max-Change = 0.0161, Max-Change = 0.0147, Max-Change = 0.0220, Max-Change = 0.0173, Max-Change = 0.0165, Max-Change = 0.0153, Max-Change = 0.0142, Max-Change = 0.0109, Max-Change = 0.0173, Max-Change = 0.0212, Max-Change = 0.0155, Max-Change = 0.0207, Max-Change = 0.0158, Max-Change = 0.0218, Max-Change = 0.0238, Max-Change = 0.0150, Max-Change = 0.0181, Max-Change = 0.0211, Max-Change = 0.0203, Max-Change = 0.0145, Max-Change = 0.0190, Max-Change = 0.0110, Max-Change = 0.0135, Max-Change = 0.0153, Max-Change = 0.0077, Max-Change = 0.0108, Max-Change = 0.0227, Max-Change = 0.0123, Max-Change = 0.0181, Max-Change = 0.0098, Max-Change = 0.0138, Max-Change = 0.0080, Max-Change = 0.0181, Max-Change = 0.0157, Max-Change = 0.0196, Max-Change = 0.0087, Max-Change = 0.0069, Max-Change = 0.0177, Max-Change = 0.0112, Max-Change = 0.0124, Max-Change = 0.0177, Max-Change = 0.0126, Max-Change = 0.0113, Max-Change = 0.0161, Max-Change = 0.0148, Max-Change = 0.0178, Max-Change = 0.0090, Max-Change = 0.0252, Max-Change = 0.0165, Max-Change = 0.0113, Max-Change = 0.0218, Max-Change = 0.0137, Max-Change = 0.0158, Max-Change = 0.0125, Max-Change = 0.0205, Max-Change = 0.0102, Max-Change = 0.0209, Max-Change = 0.0157, Max-Change = 0.0233, Max-Change = 0.0130, Max-Change = 0.0194, Max-Change = 0.0069, Max-Change = 0.0165, Max-Change = 0.0131, Max-Change = 0.0244, Max-Change = 0.0131, Max-Change = 0.0201, Max-Change = 0.0140, Max-Change = 0.0105, Max-Change = 0.0267, Max-Change = 0.0181, Max-Change = 0.0254, Max-Change = 0.0143, Max-Change = 0.0081, Max-Change = 0.0130, Max-Change = 0.0250, Max-Change = 0.0111, Max-Change = 0.0075, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0137, gam = 0.1057, Max-Change = 0.0092, gam = 0.0780, Max-Change = 0.0048, gam = 0.0629, Max-Change = 0.0066, gam = 0.0532, Max-Change = 0.0055, gam = 0.0464, Max-Change = 0.0020, gam = 0.0413, Max-Change = 0.0031, gam = 0.0374, Max-Change = 0.0024, gam = 0.0342, Max-Change = 0.0043, gam = 0.0316, Max-Change = 0.0020, gam = 0.0294, Max-Change = 0.0018, gam = 0.0276, Max-Change = 0.0026, gam = 0.0260, Max-Change = 0.0016, gam = 0.0246, Max-Change = 0.0026, gam = 0.0233, Max-Change = 0.0024, gam = 0.0222, Max-Change = 0.0025, gam = 0.0212, Max-Change = 0.0025, gam = 0.0203, Max-Change = 0.0010, gam = 0.0195, Max-Change = 0.0016, gam = 0.0188, Max-Change = 0.0020, gam = 0.0181, Max-Change = 0.0009, gam = 0.0175, Max-Change = 0.0008, gam = 0.0169, Max-Change = 0.0004
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod2b)
#> 
#> Call:
#> mixedmirt(data = dat, covdata = covdata, model = model, fixed = ~0 + 
#>     items, lr.fixed = list(F1 = ~group + contvar, F2 = ~group))
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>      F1   F2
#> F1 1.01 0.00
#> F2 0.00 1.12
#> 
#> --------------
#> LATENT REGRESSION FIXED EFFECTS:
#> 
#>                 F1    F2
#> (Intercept)  0.000 0.000
#> groupg2      0.683 0.816
#> groupg3      1.727 1.635
#> contvar     -0.003 0.000
#> 
#>             Std.Error_F1 Std.Error_F2   z_F1  z_F2
#> (Intercept)           NA           NA     NA    NA
#> groupg2            0.151        0.137  4.540 5.963
#> groupg3            0.138        0.178 12.476 9.193
#> contvar            0.056           NA -0.051    NA

####################################################
## Simulated Multilevel Rasch Model

set.seed(1)
N <- 2000
a <- matrix(rep(1,10),10,1)
d <- matrix(rnorm(10))
cluster = 100
random_intercept = rnorm(cluster,0,1)
Theta = numeric()
for (i in 1:cluster)
    Theta <- c(Theta, rnorm(N/cluster,0,1) + random_intercept[i])

group = factor(rep(paste0('G',1:cluster), each = N/cluster))
covdata <- data.frame(group)
dat <- simdata(a,d,N, itemtype = rep('2PL',10), Theta=matrix(Theta))
itemstats(dat)
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  2000            5.414          2.749  0.25 0.019 0.769     1.321
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  2000 2 0.424 0.494   0.573         0.433       0.750
#> Item_2  2000 2 0.560 0.497   0.602         0.467       0.745
#> Item_3  2000 2 0.373 0.484   0.547         0.405       0.754
#> Item_4  2000 2 0.781 0.414   0.524         0.402       0.754
#> Item_5  2000 2 0.577 0.494   0.585         0.447       0.748
#> Item_6  2000 2 0.378 0.485   0.571         0.434       0.750
#> Item_7  2000 2 0.598 0.491   0.568         0.428       0.751
#> Item_8  2000 2 0.652 0.476   0.568         0.432       0.750
#> Item_9  2000 2 0.620 0.486   0.577         0.440       0.749
#> Item_10 2000 2 0.453 0.498   0.584         0.445       0.748
#> 
#> $proportions
#>             0     1
#> Item_1  0.577 0.424
#> Item_2  0.440 0.560
#> Item_3  0.627 0.373
#> Item_4  0.219 0.781
#> Item_5  0.424 0.577
#> Item_6  0.622 0.378
#> Item_7  0.402 0.598
#> Item_8  0.348 0.652
#> Item_9  0.380 0.620
#> Item_10 0.547 0.453
#> 

# null model
mod1 <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1825, Max-Change = 0.1571, Max-Change = 0.1285, Max-Change = 0.1118, Max-Change = 0.0969, Max-Change = 0.0772, Max-Change = 0.0708, Max-Change = 0.0561, Max-Change = 0.0421, Max-Change = 0.0388, Max-Change = 0.0371, Max-Change = 0.0292, Max-Change = 0.0279, Max-Change = 0.0280, Max-Change = 0.0298, Max-Change = 0.0281, Max-Change = 0.0155, Max-Change = 0.0170, Max-Change = 0.0196, Max-Change = 0.0178, Max-Change = 0.0237, Max-Change = 0.0187, Max-Change = 0.0213, Max-Change = 0.0099, Max-Change = 0.0244, Max-Change = 0.0166, Max-Change = 0.0162, Max-Change = 0.0050, Max-Change = 0.0021, Max-Change = 0.0102, Max-Change = 0.0067, Max-Change = 0.0116, Max-Change = 0.0072, Max-Change = 0.0132, Max-Change = 0.0120, Max-Change = 0.0046, Max-Change = 0.0212, Max-Change = 0.0191, Max-Change = 0.0018, Max-Change = 0.0017, Max-Change = 0.0211, Max-Change = 0.0152, Max-Change = 0.0106, Max-Change = 0.0155, Max-Change = 0.0061, Max-Change = 0.0041, Max-Change = 0.0070, Max-Change = 0.0047, Max-Change = 0.0054, Max-Change = 0.0031, Max-Change = 0.0151, Max-Change = 0.0071, Max-Change = 0.0130, Max-Change = 0.0058, Max-Change = 0.0035, Max-Change = 0.0044, Max-Change = 0.0068, Max-Change = 0.0106, Max-Change = 0.0058, Max-Change = 0.0065, Max-Change = 0.0030, Max-Change = 0.0078, Max-Change = 0.0160, Max-Change = 0.0089, Max-Change = 0.0077, Max-Change = 0.0106, Max-Change = 0.0092, Max-Change = 0.0155, Max-Change = 0.0181, Max-Change = 0.0104, Max-Change = 0.0151, Max-Change = 0.0145, Max-Change = 0.0030, Max-Change = 0.0058, Max-Change = 0.0041, Max-Change = 0.0064, Max-Change = 0.0187, Max-Change = 0.0024, Max-Change = 0.0019, Max-Change = 0.0066, Max-Change = 0.0035, Max-Change = 0.0057, Max-Change = 0.0059, Max-Change = 0.0063, Max-Change = 0.0050, Max-Change = 0.0150, Max-Change = 0.0118, Max-Change = 0.0072, Max-Change = 0.0049, Max-Change = 0.0037, Max-Change = 0.0012, Max-Change = 0.0023, Max-Change = 0.0076, Max-Change = 0.0033, Max-Change = 0.0071, Max-Change = 0.0045, Max-Change = 0.0169, Max-Change = 0.1604, Max-Change = 0.0823, Max-Change = 0.0655, Max-Change = 0.0506, Max-Change = 0.0506, Max-Change = 0.0465, Max-Change = 0.0366, Max-Change = 0.0268, Max-Change = 0.0292, Max-Change = 0.0274, Max-Change = 0.0089, Max-Change = 0.0328, Max-Change = 0.0200, Max-Change = 0.0162, Max-Change = 0.0139, Max-Change = 0.0238, Max-Change = 0.0122, Max-Change = 0.0054, Max-Change = 0.0121, Max-Change = 0.0102, Max-Change = 0.0059, Max-Change = 0.0210, Max-Change = 0.0241, Max-Change = 0.0091, Max-Change = 0.0027, Max-Change = 0.0034, Max-Change = 0.0108, Max-Change = 0.0070, Max-Change = 0.0126, Max-Change = 0.0176, Max-Change = 0.0136, Max-Change = 0.0162, Max-Change = 0.0077, Max-Change = 0.0107, Max-Change = 0.0056, Max-Change = 0.0137, Max-Change = 0.0097, Max-Change = 0.0057, Max-Change = 0.0045, Max-Change = 0.0058, Max-Change = 0.0074, Max-Change = 0.0012, Max-Change = 0.0096, Max-Change = 0.0073, Max-Change = 0.0058, Max-Change = 0.0048, Max-Change = 0.0091, Max-Change = 0.0057, Max-Change = 0.0058, Max-Change = 0.0106, Max-Change = 0.0070, Max-Change = 0.0027, Max-Change = 0.0044, Max-Change = 0.0036, Max-Change = 0.0032, Max-Change = 0.0063, Max-Change = 0.0053, Max-Change = 0.0049, Max-Change = 0.0114, Max-Change = 0.0080, Max-Change = 0.0054, Max-Change = 0.0143, Max-Change = 0.0183, Max-Change = 0.0110, Max-Change = 0.0025, Max-Change = 0.0047, Max-Change = 0.0071, Max-Change = 0.0039, Max-Change = 0.0101, Max-Change = 0.0093, Max-Change = 0.0067, Max-Change = 0.0014, Max-Change = 0.0041, Max-Change = 0.0047, Max-Change = 0.0073, Max-Change = 0.0108, Max-Change = 0.0033, Max-Change = 0.0120, Max-Change = 0.0109, Max-Change = 0.0149, Max-Change = 0.0138, Max-Change = 0.0082, Max-Change = 0.0119, Max-Change = 0.0063, Max-Change = 0.0074, Max-Change = 0.0123, Max-Change = 0.0063, Max-Change = 0.0044, Max-Change = 0.0047, Max-Change = 0.0098, Max-Change = 0.0102, Max-Change = 0.0207, Max-Change = 0.0144, Max-Change = 0.0101, Max-Change = 0.0039, Max-Change = 0.0079, Max-Change = 0.0092, Max-Change = 0.0020, Max-Change = 0.0040, Max-Change = 0.0080, Max-Change = 0.0104, Max-Change = 0.0072, Max-Change = 0.0118, Max-Change = 0.0109, Max-Change = 0.0041, Max-Change = 0.0042, Max-Change = 0.0096, Max-Change = 0.0043, Max-Change = 0.0046, Max-Change = 0.0031, Max-Change = 0.0090, Max-Change = 0.0019, Max-Change = 0.0042, Max-Change = 0.0068, Max-Change = 0.0022, Max-Change = 0.0023, Max-Change = 0.0132, Max-Change = 0.0032, Max-Change = 0.0025, Max-Change = 0.0070, Max-Change = 0.0065, Max-Change = 0.0040, Max-Change = 0.0055, Max-Change = 0.0056, Max-Change = 0.0031, Max-Change = 0.0085, Max-Change = 0.0067, Max-Change = 0.0029, Max-Change = 0.0080, Max-Change = 0.0046, Max-Change = 0.0137, Max-Change = 0.0070, Max-Change = 0.0052, Max-Change = 0.0119, Max-Change = 0.0087, Max-Change = 0.0039, Max-Change = 0.0069, Max-Change = 0.0050, Max-Change = 0.0104, Max-Change = 0.0009, Max-Change = 0.0057, Max-Change = 0.0095, Max-Change = 0.0100, Max-Change = 0.0024, Max-Change = 0.0044, Max-Change = 0.0044, Max-Change = 0.0057, Max-Change = 0.0026, Max-Change = 0.0058, Max-Change = 0.0074, Max-Change = 0.0071, Max-Change = 0.0028, Max-Change = 0.0094, Max-Change = 0.0088, Max-Change = 0.0034, Max-Change = 0.0017, Max-Change = 0.0107, Max-Change = 0.0048, Max-Change = 0.0058, Max-Change = 0.0038, Max-Change = 0.0046, Max-Change = 0.0038, Max-Change = 0.0072, Max-Change = 0.0089, Max-Change = 0.0011, Max-Change = 0.0016, Max-Change = 0.0039, Max-Change = 0.0060, Max-Change = 0.0053, Max-Change = 0.0025, Max-Change = 0.0028, Max-Change = 0.0076, Max-Change = 0.0029, Max-Change = 0.0118, Max-Change = 0.0030, Max-Change = 0.0035, Max-Change = 0.0043, Max-Change = 0.0055, Max-Change = 0.0024, Max-Change = 0.0072, Max-Change = 0.0038, Max-Change = 0.0030, Max-Change = 0.0019, Max-Change = 0.0023, Max-Change = 0.0024, Max-Change = 0.0070, Max-Change = 0.0108, Max-Change = 0.0044, Max-Change = 0.0071, Max-Change = 0.0100, Max-Change = 0.0029, Max-Change = 0.0073, Max-Change = 0.0056, Max-Change = 0.0049, Max-Change = 0.0081, Max-Change = 0.0053, Max-Change = 0.0076, Max-Change = 0.0016, Max-Change = 0.0056, Max-Change = 0.0062, Max-Change = 0.0020, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0020, gam = 0.1057, Max-Change = 0.0083, gam = 0.0780, Max-Change = 0.0013, gam = 0.0629, Max-Change = 0.0017, gam = 0.0532, Max-Change = 0.0011, gam = 0.0464, Max-Change = 0.0002, gam = 0.0413, Max-Change = 0.0012, gam = 0.0374, Max-Change = 0.0004, gam = 0.0342, Max-Change = 0.0012, gam = 0.0316, Max-Change = 0.0003, gam = 0.0294, Max-Change = 0.0004, gam = 0.0276, Max-Change = 0.0010, gam = 0.0260, Max-Change = 0.0006, gam = 0.0246, Max-Change = 0.0004, gam = 0.0233, Max-Change = 0.0003
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod1)
#> 
#> Call:
#> mixedmirt(data = dat, covdata = covdata, model = 1, fixed = ~0 + 
#>     items, random = ~1 | group)
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>      F1
#> F1 1.02
#> 
#> $group
#>           COV_group
#> COV_group     0.853
#> 

# include level 2 predictor for 'group' variance
covdata$group_pred <- rep(random_intercept, each = N/cluster)
mod2 <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items + group_pred, random = ~ 1|group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1797, Max-Change = 0.1540, Max-Change = 0.1292, Max-Change = 0.1114, Max-Change = 0.0981, Max-Change = 0.0797, Max-Change = 0.0726, Max-Change = 0.0565, Max-Change = 0.0436, Max-Change = 0.0397, Max-Change = 0.0353, Max-Change = 0.0313, Max-Change = 0.0278, Max-Change = 0.0326, Max-Change = 0.0265, Max-Change = 0.0204, Max-Change = 0.0212, Max-Change = 0.0173, Max-Change = 0.0144, Max-Change = 0.0184, Max-Change = 0.0145, Max-Change = 0.0115, Max-Change = 0.0066, Max-Change = 0.0137, Max-Change = 0.0068, Max-Change = 0.0131, Max-Change = 0.0088, Max-Change = 0.0131, Max-Change = 0.0121, Max-Change = 0.0098, Max-Change = 0.0085, Max-Change = 0.0064, Max-Change = 0.0127, Max-Change = 0.0081, Max-Change = 0.0079, Max-Change = 0.0073, Max-Change = 0.0062, Max-Change = 0.0087, Max-Change = 0.0085, Max-Change = 0.0046, Max-Change = 0.0091, Max-Change = 0.0039, Max-Change = 0.0077, Max-Change = 0.0035, Max-Change = 0.0032, Max-Change = 0.0047, Max-Change = 0.0046, Max-Change = 0.0074, Max-Change = 0.0039, Max-Change = 0.0037, Max-Change = 0.0104, Max-Change = 0.0064, Max-Change = 0.0067, Max-Change = 0.0035, Max-Change = 0.0043, Max-Change = 0.0088, Max-Change = 0.0028, Max-Change = 0.0047, Max-Change = 0.0040, Max-Change = 0.0024, Max-Change = 0.0078, Max-Change = 0.0162, Max-Change = 0.0064, Max-Change = 0.0054, Max-Change = 0.0070, Max-Change = 0.0016, Max-Change = 0.0066, Max-Change = 0.0069, Max-Change = 0.0076, Max-Change = 0.0078, Max-Change = 0.0057, Max-Change = 0.0103, Max-Change = 0.0023, Max-Change = 0.0051, Max-Change = 0.0057, Max-Change = 0.0055, Max-Change = 0.0184, Max-Change = 0.0044, Max-Change = 0.0037, Max-Change = 0.0050, Max-Change = 0.0042, Max-Change = 0.0070, Max-Change = 0.0010, Max-Change = 0.0078, Max-Change = 0.0022, Max-Change = 0.0104, Max-Change = 0.0058, Max-Change = 0.0075, Max-Change = 0.0049, Max-Change = 0.0010, Max-Change = 0.0014, Max-Change = 0.0033, Max-Change = 0.0029, Max-Change = 0.0033, Max-Change = 0.0053, Max-Change = 0.0015, Max-Change = 0.0026, Max-Change = 0.0172, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1757, Max-Change = 0.1798, Max-Change = 0.1665, Max-Change = 0.0043, Max-Change = 0.0036, Max-Change = 0.0047, Max-Change = 0.0044, Max-Change = 0.0095, Max-Change = 0.0068, Max-Change = 0.0034, Max-Change = 0.0076, Max-Change = 0.0038, Max-Change = 0.0058, Max-Change = 0.0064, Max-Change = 0.0070, Max-Change = 0.0053, Max-Change = 0.0115, Max-Change = 0.0084, Max-Change = 0.0088, Max-Change = 0.0042, Max-Change = 0.0051, Max-Change = 0.0060, Max-Change = 0.0075, Max-Change = 0.0062, Max-Change = 0.0056, Max-Change = 0.0088, Max-Change = 0.0070, Max-Change = 0.0040, Max-Change = 0.0161, Max-Change = 0.0061, Max-Change = 0.0044, Max-Change = 0.0038, Max-Change = 0.0070, Max-Change = 0.0018, Max-Change = 0.0071, Max-Change = 0.0095, Max-Change = 0.0081, Max-Change = 0.0092, Max-Change = 0.0070, Max-Change = 0.0144, Max-Change = 0.0061, Max-Change = 0.0084, Max-Change = 0.0076, Max-Change = 0.0065, Max-Change = 0.0046, Max-Change = 0.0049, Max-Change = 0.0092, Max-Change = 0.0063, Max-Change = 0.0027, Max-Change = 0.0061, Max-Change = 0.0085, Max-Change = 0.0043, Max-Change = 0.0063, Max-Change = 0.0052, Max-Change = 0.0028, Max-Change = 0.0054, Max-Change = 0.0076, Max-Change = 0.0027, Max-Change = 0.0039, Max-Change = 0.0033, Max-Change = 0.0041, Max-Change = 0.0082, Max-Change = 0.0050, Max-Change = 0.0104, Max-Change = 0.0078, Max-Change = 0.0037, Max-Change = 0.0060, Max-Change = 0.0080, Max-Change = 0.0045, Max-Change = 0.0068, Max-Change = 0.0039, Max-Change = 0.0078, Max-Change = 0.0041, Max-Change = 0.0038, Max-Change = 0.0019, Max-Change = 0.0061, Max-Change = 0.0112, Max-Change = 0.0074, Max-Change = 0.0025, Max-Change = 0.0039, Max-Change = 0.0087, Max-Change = 0.0016, Max-Change = 0.0032, Max-Change = 0.0020, Max-Change = 0.0066, Max-Change = 0.0022, Max-Change = 0.0048, Max-Change = 0.0050, Max-Change = 0.0066, Max-Change = 0.0036, Max-Change = 0.0040, Max-Change = 0.0049, Max-Change = 0.0049, Max-Change = 0.0078, Max-Change = 0.0041, Max-Change = 0.0093, Max-Change = 0.0046, Max-Change = 0.0109, Max-Change = 0.0053, Max-Change = 0.0105, Max-Change = 0.0050, Max-Change = 0.0022, Max-Change = 0.0057, Max-Change = 0.0064, Max-Change = 0.0040, Max-Change = 0.0044, Max-Change = 0.0053, Max-Change = 0.0095, Max-Change = 0.0033, Max-Change = 0.0030, Max-Change = 0.0066, Max-Change = 0.0019, Max-Change = 0.0035, Max-Change = 0.0124, Max-Change = 0.0072, Max-Change = 0.0072, Max-Change = 0.0075, Max-Change = 0.0074, Max-Change = 0.0034, Max-Change = 0.0055, Max-Change = 0.0036, Max-Change = 0.0020, Max-Change = 0.0090, Max-Change = 0.0070, Max-Change = 0.0048, Max-Change = 0.0064, Max-Change = 0.0046, Max-Change = 0.0049, Max-Change = 0.0049, Max-Change = 0.0042, Max-Change = 0.0090, Max-Change = 0.0049, Max-Change = 0.0034, Max-Change = 0.0071, Max-Change = 0.0040, Max-Change = 0.0124, Max-Change = 0.0072, Max-Change = 0.0065, Max-Change = 0.0072, Max-Change = 0.0140, Max-Change = 0.0046, Max-Change = 0.0095, Max-Change = 0.0052, Max-Change = 0.0078, Max-Change = 0.0026, Max-Change = 0.0081, Max-Change = 0.0097, Max-Change = 0.0031, Max-Change = 0.0037, Max-Change = 0.0112, Max-Change = 0.0082, Max-Change = 0.0051, Max-Change = 0.0026, Max-Change = 0.0031, Max-Change = 0.0055, Max-Change = 0.0118, Max-Change = 0.0060, Max-Change = 0.0030, Max-Change = 0.0054, Max-Change = 0.0032, Max-Change = 0.0106, Max-Change = 0.0020, Max-Change = 0.0077, Max-Change = 0.0047, Max-Change = 0.0053, Max-Change = 0.0060, Max-Change = 0.0062, Max-Change = 0.0019, Max-Change = 0.0039, Max-Change = 0.0009, Max-Change = 0.0143, Max-Change = 0.0031, Max-Change = 0.0069, Max-Change = 0.0050, Max-Change = 0.0027, Max-Change = 0.0020, Max-Change = 0.0085, Max-Change = 0.0023, Max-Change = 0.0040, Max-Change = 0.0049, Max-Change = 0.0040, Max-Change = 0.0048, Max-Change = 0.0056, Max-Change = 0.0050, Max-Change = 0.0046, Max-Change = 0.0065, Max-Change = 0.0133, Max-Change = 0.0035, Max-Change = 0.0059, Max-Change = 0.0067, Max-Change = 0.0044, Max-Change = 0.0060, Max-Change = 0.0030, Max-Change = 0.0072, Max-Change = 0.0066, Max-Change = 0.0067, Max-Change = 0.0045, Max-Change = 0.0033, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0026, gam = 0.1057, Max-Change = 0.0077, gam = 0.0780, Max-Change = 0.0018, gam = 0.0629, Max-Change = 0.0029, gam = 0.0532, Max-Change = 0.0016, gam = 0.0464, Max-Change = 0.0007, gam = 0.0413, Max-Change = 0.0006, gam = 0.0374, Max-Change = 0.0004
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...

# including group means predicts nearly all variability in 'group'
summary(mod2)
#> 
#> Call:
#> mixedmirt(data = dat, covdata = covdata, model = 1, fixed = ~0 + 
#>     items + group_pred, random = ~1 | group)
#> 
#> --------------
#> FIXED EFFECTS:
#>            Estimate Std.Error z.value
#> group_pred    1.017      0.05  20.282
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>      F1
#> F1 1.04
#> 
#> $group
#>           COV_group
#> COV_group    0.0335
#> 
anova(mod1, mod2)
#>           AIC    SABIC       HQ      BIC    logLik     X2 df p
#> mod1 23548.80 23577.89 23573.48 23616.01 -11762.40            
#> mod2 22744.85 22776.36 22771.59 22817.66 -11359.43 805.95  1 0

# can also be fit for Rasch/non-Rasch models with the lr.random input
mod1b <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items, lr.random = ~ 1|group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1825, Max-Change = 0.1571, Max-Change = 0.1285, Max-Change = 0.1118, Max-Change = 0.0969, Max-Change = 0.0772, Max-Change = 0.0708, Max-Change = 0.0561, Max-Change = 0.0421, Max-Change = 0.0388, Max-Change = 0.0371, Max-Change = 0.0292, Max-Change = 0.0279, Max-Change = 0.0280, Max-Change = 0.0298, Max-Change = 0.0281, Max-Change = 0.0155, Max-Change = 0.0170, Max-Change = 0.0196, Max-Change = 0.0178, Max-Change = 0.0237, Max-Change = 0.0187, Max-Change = 0.0213, Max-Change = 0.0099, Max-Change = 0.0244, Max-Change = 0.0166, Max-Change = 0.0162, Max-Change = 0.0050, Max-Change = 0.0021, Max-Change = 0.0102, Max-Change = 0.0067, Max-Change = 0.0116, Max-Change = 0.0072, Max-Change = 0.0132, Max-Change = 0.0120, Max-Change = 0.0046, Max-Change = 0.0212, Max-Change = 0.0191, Max-Change = 0.0018, Max-Change = 0.0017, Max-Change = 0.0211, Max-Change = 0.0152, Max-Change = 0.0106, Max-Change = 0.0155, Max-Change = 0.0061, Max-Change = 0.0041, Max-Change = 0.0070, Max-Change = 0.0047, Max-Change = 0.0054, Max-Change = 0.0031, Max-Change = 0.0151, Max-Change = 0.0071, Max-Change = 0.0130, Max-Change = 0.0058, Max-Change = 0.0035, Max-Change = 0.0044, Max-Change = 0.0068, Max-Change = 0.0106, Max-Change = 0.0058, Max-Change = 0.0065, Max-Change = 0.0030, Max-Change = 0.0078, Max-Change = 0.0160, Max-Change = 0.0089, Max-Change = 0.0077, Max-Change = 0.0106, Max-Change = 0.0092, Max-Change = 0.0155, Max-Change = 0.0181, Max-Change = 0.0104, Max-Change = 0.0151, Max-Change = 0.0145, Max-Change = 0.0030, Max-Change = 0.0058, Max-Change = 0.0041, Max-Change = 0.0064, Max-Change = 0.0187, Max-Change = 0.0024, Max-Change = 0.0019, Max-Change = 0.0066, Max-Change = 0.0035, Max-Change = 0.0057, Max-Change = 0.0059, Max-Change = 0.0063, Max-Change = 0.0050, Max-Change = 0.0150, Max-Change = 0.0118, Max-Change = 0.0072, Max-Change = 0.0049, Max-Change = 0.0037, Max-Change = 0.0012, Max-Change = 0.0023, Max-Change = 0.0076, Max-Change = 0.0033, Max-Change = 0.0071, Max-Change = 0.0045, Max-Change = 0.0169, Max-Change = 0.0704, Max-Change = 0.0689, Max-Change = 0.0489, Max-Change = 0.0433, Max-Change = 0.0410, Max-Change = 0.0380, Max-Change = 0.0347, Max-Change = 0.0321, Max-Change = 0.0280, Max-Change = 0.0249, Max-Change = 0.0214, Max-Change = 0.0182, Max-Change = 0.0159, Max-Change = 0.0134, Max-Change = 0.0112, Max-Change = 0.0144, Max-Change = 0.0072, Max-Change = 0.0063, Max-Change = 0.0066, Max-Change = 0.0080, Max-Change = 0.0093, Max-Change = 0.0062, Max-Change = 0.0128, Max-Change = 0.0042, Max-Change = 0.0035, Max-Change = 0.0051, Max-Change = 0.0045, Max-Change = 0.0125, Max-Change = 0.0124, Max-Change = 0.0186, Max-Change = 0.0070, Max-Change = 0.0146, Max-Change = 0.0084, Max-Change = 0.0090, Max-Change = 0.0087, Max-Change = 0.0026, Max-Change = 0.0069, Max-Change = 0.0104, Max-Change = 0.0086, Max-Change = 0.0023, Max-Change = 0.0116, Max-Change = 0.0048, Max-Change = 0.0024, Max-Change = 0.0099, Max-Change = 0.0040, Max-Change = 0.0122, Max-Change = 0.0177, Max-Change = 0.0101, Max-Change = 0.0027, Max-Change = 0.0057, Max-Change = 0.0282, Max-Change = 0.0083, Max-Change = 0.0079, Max-Change = 0.0087, Max-Change = 0.0150, Max-Change = 0.0057, Max-Change = 0.0038, Max-Change = 0.0046, Max-Change = 0.0086, Max-Change = 0.0077, Max-Change = 0.0059, Max-Change = 0.0037, Max-Change = 0.0041, Max-Change = 0.0053, Max-Change = 0.0033, Max-Change = 0.0042, Max-Change = 0.0113, Max-Change = 0.0037, Max-Change = 0.0240, Max-Change = 0.0113, Max-Change = 0.0101, Max-Change = 0.0107, Max-Change = 0.0018, Max-Change = 0.0060, Max-Change = 0.0068, Max-Change = 0.0036, Max-Change = 0.0208, Max-Change = 0.0056, Max-Change = 0.0050, Max-Change = 0.0046, Max-Change = 0.0040, Max-Change = 0.0077, Max-Change = 0.0057, Max-Change = 0.0027, Max-Change = 0.0062, Max-Change = 0.0106, Max-Change = 0.0084, Max-Change = 0.0156, Max-Change = 0.0045, Max-Change = 0.0033, Max-Change = 0.0047, Max-Change = 0.0060, Max-Change = 0.0102, Max-Change = 0.0108, Max-Change = 0.0071, Max-Change = 0.0064, Max-Change = 0.0054, Max-Change = 0.0063, Max-Change = 0.0056, Max-Change = 0.0072, Max-Change = 0.0031, Max-Change = 0.0037, Max-Change = 0.0040, Max-Change = 0.0152, Max-Change = 0.0039, Max-Change = 0.0064, Max-Change = 0.0041, Max-Change = 0.0056, Max-Change = 0.0082, Max-Change = 0.0027, Max-Change = 0.0019, Max-Change = 0.0034, Max-Change = 0.0036, Max-Change = 0.0074, Max-Change = 0.0074, Max-Change = 0.0067, Max-Change = 0.0098, Max-Change = 0.0066, Max-Change = 0.0039, Max-Change = 0.0180, Max-Change = 0.0060, Max-Change = 0.0005, Max-Change = 0.0055, Max-Change = 0.0100, Max-Change = 0.0078, Max-Change = 0.0070, Max-Change = 0.0094, Max-Change = 0.0027, Max-Change = 0.0117, Max-Change = 0.0056, Max-Change = 0.0116, Max-Change = 0.0065, Max-Change = 0.0048, Max-Change = 0.0084, Max-Change = 0.0158, Max-Change = 0.0072, Max-Change = 0.0078, Max-Change = 0.0027, Max-Change = 0.0103, Max-Change = 0.0035, Max-Change = 0.0060, Max-Change = 0.0029, Max-Change = 0.0074, Max-Change = 0.0067, Max-Change = 0.0049, Max-Change = 0.0053, Max-Change = 0.0152, Max-Change = 0.0037, Max-Change = 0.0036, Max-Change = 0.0120, Max-Change = 0.0011, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0146, gam = 0.1057, Max-Change = 0.0034, gam = 0.0780, Max-Change = 0.0012, gam = 0.0629, Max-Change = 0.0018, gam = 0.0532, Max-Change = 0.0005, gam = 0.0464, Max-Change = 0.0018, gam = 0.0413, Max-Change = 0.0019, gam = 0.0374, Max-Change = 0.0007, gam = 0.0342, Max-Change = 0.0006, gam = 0.0316, Max-Change = 0.0006
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod1b)
#> 
#> Call:
#> mixedmirt(data = dat, covdata = covdata, model = 1, fixed = ~0 + 
#>     items, lr.random = ~1 | group)
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>      F1
#> F1 1.54
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $group
#>           COV_group
#> COV_group      1.45
#> 

mod2b <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items + group_pred, lr.random = ~ 1|group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1797, Max-Change = 0.1540, Max-Change = 0.1292, Max-Change = 0.1114, Max-Change = 0.0981, Max-Change = 0.0797, Max-Change = 0.0726, Max-Change = 0.0565, Max-Change = 0.0436, Max-Change = 0.0397, Max-Change = 0.0353, Max-Change = 0.0313, Max-Change = 0.0278, Max-Change = 0.0326, Max-Change = 0.0265, Max-Change = 0.0204, Max-Change = 0.0212, Max-Change = 0.0173, Max-Change = 0.0144, Max-Change = 0.0184, Max-Change = 0.0145, Max-Change = 0.0115, Max-Change = 0.0066, Max-Change = 0.0137, Max-Change = 0.0068, Max-Change = 0.0131, Max-Change = 0.0088, Max-Change = 0.0131, Max-Change = 0.0121, Max-Change = 0.0098, Max-Change = 0.0085, Max-Change = 0.0064, Max-Change = 0.0127, Max-Change = 0.0081, Max-Change = 0.0079, Max-Change = 0.0073, Max-Change = 0.0062, Max-Change = 0.0087, Max-Change = 0.0085, Max-Change = 0.0046, Max-Change = 0.0091, Max-Change = 0.0039, Max-Change = 0.0077, Max-Change = 0.0035, Max-Change = 0.0032, Max-Change = 0.0047, Max-Change = 0.0046, Max-Change = 0.0074, Max-Change = 0.0039, Max-Change = 0.0037, Max-Change = 0.0104, Max-Change = 0.0064, Max-Change = 0.0067, Max-Change = 0.0035, Max-Change = 0.0043, Max-Change = 0.0088, Max-Change = 0.0028, Max-Change = 0.0047, Max-Change = 0.0040, Max-Change = 0.0024, Max-Change = 0.0078, Max-Change = 0.0162, Max-Change = 0.0064, Max-Change = 0.0054, Max-Change = 0.0070, Max-Change = 0.0016, Max-Change = 0.0066, Max-Change = 0.0069, Max-Change = 0.0076, Max-Change = 0.0078, Max-Change = 0.0057, Max-Change = 0.0103, Max-Change = 0.0023, Max-Change = 0.0051, Max-Change = 0.0057, Max-Change = 0.0055, Max-Change = 0.0184, Max-Change = 0.0044, Max-Change = 0.0037, Max-Change = 0.0050, Max-Change = 0.0042, Max-Change = 0.0070, Max-Change = 0.0010, Max-Change = 0.0078, Max-Change = 0.0022, Max-Change = 0.0104, Max-Change = 0.0058, Max-Change = 0.0075, Max-Change = 0.0049, Max-Change = 0.0010, Max-Change = 0.0014, Max-Change = 0.0033, Max-Change = 0.0029, Max-Change = 0.0033, Max-Change = 0.0053, Max-Change = 0.0015, Max-Change = 0.0026, Max-Change = 0.0504, Max-Change = 0.0563, Max-Change = 0.0438, Max-Change = 0.0414, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.0674, Max-Change = 0.0488, Max-Change = 0.1112, Max-Change = 0.0397, Max-Change = 0.0278, Max-Change = 0.0117, Max-Change = 0.0143, Max-Change = 0.0157, Max-Change = 0.0098, Max-Change = 0.0107, Max-Change = 0.0456, Max-Change = 0.0102, Max-Change = 0.0108, Max-Change = 0.0286, Max-Change = 0.0242, Max-Change = 0.0158, Max-Change = 0.0098, Max-Change = 0.0145, Max-Change = 0.0062, Max-Change = 0.0076, Max-Change = 0.0117, Max-Change = 0.0107, Max-Change = 0.0120, Max-Change = 0.0100, Max-Change = 0.0085, Max-Change = 0.0268, Max-Change = 0.0089, Max-Change = 0.0072, Max-Change = 0.0066, Max-Change = 0.0051, Max-Change = 0.0076, Max-Change = 0.0049, Max-Change = 0.0051, Max-Change = 0.0058, Max-Change = 0.0080, Max-Change = 0.0052, Max-Change = 0.0136, Max-Change = 0.0078, Max-Change = 0.0079, Max-Change = 0.0079, Max-Change = 0.0061, Max-Change = 0.0028, Max-Change = 0.0061, Max-Change = 0.0110, Max-Change = 0.0032, Max-Change = 0.0023, Max-Change = 0.0071, Max-Change = 0.0048, Max-Change = 0.0087, Max-Change = 0.0048, Max-Change = 0.0096, Max-Change = 0.0060, Max-Change = 0.0052, Max-Change = 0.0035, Max-Change = 0.0057, Max-Change = 0.0061, Max-Change = 0.0036, Max-Change = 0.0041, Max-Change = 0.0070, Max-Change = 0.0063, Max-Change = 0.0112, Max-Change = 0.0021, Max-Change = 0.0141, Max-Change = 0.0145, Max-Change = 0.0035, Max-Change = 0.0038, Max-Change = 0.0023, Max-Change = 0.0084, Max-Change = 0.0054, Max-Change = 0.0073, Max-Change = 0.0036, Max-Change = 0.0117, Max-Change = 0.0045, Max-Change = 0.0050, Max-Change = 0.0066, Max-Change = 0.0070, Max-Change = 0.0050, Max-Change = 0.0008, Max-Change = 0.0061, Max-Change = 0.0089, Max-Change = 0.0036, Max-Change = 0.0127, Max-Change = 0.0072, Max-Change = 0.0045, Max-Change = 0.0094, Max-Change = 0.0102, Max-Change = 0.0087, Max-Change = 0.0114, Max-Change = 0.0088, Max-Change = 0.0058, Max-Change = 0.0016, Max-Change = 0.0072, Max-Change = 0.0096, Max-Change = 0.0061, Max-Change = 0.0060, Max-Change = 0.0020, Max-Change = 0.0113, Max-Change = 0.0087, Max-Change = 0.0036, Max-Change = 0.0066, Max-Change = 0.0037, Max-Change = 0.0041, Max-Change = 0.0044, Max-Change = 0.0051, Max-Change = 0.0075, Max-Change = 0.0020, Max-Change = 0.0050, Max-Change = 0.0055, Max-Change = 0.0047, Max-Change = 0.0049, Max-Change = 0.0113, Max-Change = 0.0050, Max-Change = 0.0044, Max-Change = 0.0080, Max-Change = 0.0076, Max-Change = 0.0047, Max-Change = 0.0054, Max-Change = 0.0075, Max-Change = 0.0026, Max-Change = 0.0035, Max-Change = 0.0099, Max-Change = 0.0038, Max-Change = 0.0108, Max-Change = 0.0072, Max-Change = 0.0086, Max-Change = 0.0042, Max-Change = 0.0025, Max-Change = 0.0094, Max-Change = 0.0063, Max-Change = 0.0023, Max-Change = 0.0100, Max-Change = 0.0037, Max-Change = 0.0064, Max-Change = 0.0031, Max-Change = 0.0068, Max-Change = 0.0047, Max-Change = 0.0077, Max-Change = 0.0042, Max-Change = 0.0045, Max-Change = 0.0079, Max-Change = 0.0060, Max-Change = 0.0025, Max-Change = 0.0040, Max-Change = 0.0071, Max-Change = 0.0033, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0065, gam = 0.1057, Max-Change = 0.0048, gam = 0.0780, Max-Change = 0.0016, gam = 0.0629, Max-Change = 0.0011, gam = 0.0532, Max-Change = 0.0018, gam = 0.0464, Max-Change = 0.0011, gam = 0.0413, Max-Change = 0.0015, gam = 0.0374, Max-Change = 0.0010, gam = 0.0342, Max-Change = 0.0007, gam = 0.0316, Max-Change = 0.0013, gam = 0.0294, Max-Change = 0.0007, gam = 0.0276, Max-Change = 0.0015, gam = 0.0260, Max-Change = 0.0004, gam = 0.0246, Max-Change = 0.0005, gam = 0.0233, Max-Change = 0.0006
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod2b)
#> 
#> Call:
#> mixedmirt(data = dat, covdata = covdata, model = 1, fixed = ~0 + 
#>     items + group_pred, lr.random = ~1 | group)
#> 
#> --------------
#> FIXED EFFECTS:
#>            Estimate Std.Error z.value
#> group_pred    1.036     0.044  23.766
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>      F1
#> F1 1.12
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $group
#>           COV_group
#> COV_group       0.1
#> 
anova(mod1b, mod2b)
#>            AIC    SABIC       HQ      BIC    logLik     X2 df p
#> mod1b 23637.58 23664.24 23660.20 23699.19 -11807.79            
#> mod2b 22749.34 22778.43 22774.02 22816.55 -11362.67 890.24  1 0

mod3 <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items, lr.random = ~ 1|group, itemtype = '2PL')
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1738, Max-Change = 0.1500, Max-Change = 0.1263, Max-Change = 0.1114, Max-Change = 0.1014, Max-Change = 0.0822, Max-Change = 0.0827, Max-Change = 0.0623, Max-Change = 0.0533, Max-Change = 0.0476, Max-Change = 0.0405, Max-Change = 0.0350, Max-Change = 0.0288, Max-Change = 0.0341, Max-Change = 0.0200, Max-Change = 0.0251, Max-Change = 0.0278, Max-Change = 0.0296, Max-Change = 0.0227, Max-Change = 0.0251, Max-Change = 0.0128, Max-Change = 0.0190, Max-Change = 0.0346, Max-Change = 0.0151, Max-Change = 0.0156, Max-Change = 0.0152, Max-Change = 0.0234, Max-Change = 0.0312, Max-Change = 0.0126, Max-Change = 0.0161, Max-Change = 0.0104, Max-Change = 0.0176, Max-Change = 0.0121, Max-Change = 0.0184, Max-Change = 0.0109, Max-Change = 0.0192, Max-Change = 0.0312, Max-Change = 0.0234, Max-Change = 0.0202, Max-Change = 0.0160, Max-Change = 0.0187, Max-Change = 0.0125, Max-Change = 0.0146, Max-Change = 0.0210, Max-Change = 0.0128, Max-Change = 0.0280, Max-Change = 0.0095, Max-Change = 0.0179, Max-Change = 0.0180, Max-Change = 0.0198, Max-Change = 0.0162, Max-Change = 0.0227, Max-Change = 0.0190, Max-Change = 0.0175, Max-Change = 0.0161, Max-Change = 0.0179, Max-Change = 0.0180, Max-Change = 0.0196, Max-Change = 0.0331, Max-Change = 0.0150, Max-Change = 0.0124, Max-Change = 0.0152, Max-Change = 0.0151, Max-Change = 0.0181, Max-Change = 0.0116, Max-Change = 0.0185, Max-Change = 0.0169, Max-Change = 0.0155, Max-Change = 0.0132, Max-Change = 0.0220, Max-Change = 0.0219, Max-Change = 0.0216, Max-Change = 0.0183, Max-Change = 0.0150, Max-Change = 0.0161, Max-Change = 0.0175, Max-Change = 0.0191, Max-Change = 0.0234, Max-Change = 0.0195, Max-Change = 0.0121, Max-Change = 0.0264, Max-Change = 0.0156, Max-Change = 0.0221, Max-Change = 0.0141, Max-Change = 0.0168, Max-Change = 0.0245, Max-Change = 0.0215, Max-Change = 0.0192, Max-Change = 0.0181, Max-Change = 0.0134, Max-Change = 0.0205, Max-Change = 0.0155, Max-Change = 0.0164, Max-Change = 0.0185, Max-Change = 0.0150, Max-Change = 0.0166, Max-Change = 0.0210, Max-Change = 0.0473, Max-Change = 0.0461, Max-Change = 0.0412, Max-Change = 0.0413, Max-Change = 0.0340, Max-Change = 0.0316, Max-Change = 0.0291, Max-Change = 0.0308, Max-Change = 0.0263, Max-Change = 0.0238, Max-Change = 0.0303, Max-Change = 0.0196, Max-Change = 0.0270, Max-Change = 0.0153, Max-Change = 0.0199, Max-Change = 0.0317, Max-Change = 0.0225, Max-Change = 0.0130, Max-Change = 0.0198, Max-Change = 0.0186, Max-Change = 0.0191, Max-Change = 0.0128, Max-Change = 0.0198, Max-Change = 0.0147, Max-Change = 0.0193, Max-Change = 0.0131, Max-Change = 0.0147, Max-Change = 0.0083, Max-Change = 0.0112, Max-Change = 0.0137, Max-Change = 0.0113, Max-Change = 0.0099, Max-Change = 0.0167, Max-Change = 0.0127, Max-Change = 0.0201, Max-Change = 0.0106, Max-Change = 0.0136, Max-Change = 0.0101, Max-Change = 0.0236, Max-Change = 0.0105, Max-Change = 0.0137, Max-Change = 0.0112, Max-Change = 0.0255, Max-Change = 0.0077, Max-Change = 0.0084, Max-Change = 0.0166, Max-Change = 0.0123, Max-Change = 0.0106, Max-Change = 0.0126, Max-Change = 0.0116, Max-Change = 0.0607, Max-Change = 0.0115, Max-Change = 0.0106, Max-Change = 0.0147, Max-Change = 0.0134, Max-Change = 0.0115, Max-Change = 0.0167, Max-Change = 0.0107, Max-Change = 0.0197, Max-Change = 0.0114, Max-Change = 0.0129, Max-Change = 0.0109, Max-Change = 0.0091, Max-Change = 0.0131, Max-Change = 0.0050, Max-Change = 0.0086, Max-Change = 0.0109, Max-Change = 0.0122, Max-Change = 0.0118, Max-Change = 0.0093, Max-Change = 0.0092, Max-Change = 0.0115, Max-Change = 0.0121, Max-Change = 0.0123, Max-Change = 0.0128, Max-Change = 0.0101, Max-Change = 0.0093, Max-Change = 0.0131, Max-Change = 0.0152, Max-Change = 0.0073, Max-Change = 0.0066, Max-Change = 0.0100, Max-Change = 0.0141, Max-Change = 0.0091, Max-Change = 0.0114, Max-Change = 0.0183, Max-Change = 0.0095, Max-Change = 0.0086, Max-Change = 0.0208, Max-Change = 0.0086, Max-Change = 0.0149, Max-Change = 0.0069, Max-Change = 0.0095, Max-Change = 0.0091, Max-Change = 0.0139, Max-Change = 0.0117, Max-Change = 0.0094, Max-Change = 0.0104, Max-Change = 0.0130, Max-Change = 0.0095, Max-Change = 0.0111, Max-Change = 0.0144, Max-Change = 0.0103, Max-Change = 0.0111, Max-Change = 0.0168, Max-Change = 0.0082, Max-Change = 0.0091, Max-Change = 0.0100, Max-Change = 0.0122, Max-Change = 0.0162, Max-Change = 0.0082, Max-Change = 0.0097, Max-Change = 0.0077, Max-Change = 0.0198, Max-Change = 0.0207, Max-Change = 0.0113, Max-Change = 0.0158, Max-Change = 0.0121, Max-Change = 0.0188, Max-Change = 0.0107, Max-Change = 0.0114, Max-Change = 0.0156, Max-Change = 0.0122, Max-Change = 0.0133, Max-Change = 0.0126, Max-Change = 0.0077, Max-Change = 0.0111, Max-Change = 0.0134, Max-Change = 0.0164, Max-Change = 0.0138, Max-Change = 0.0115, Max-Change = 0.0137, Max-Change = 0.0124, Max-Change = 0.0099, Max-Change = 0.0106, Max-Change = 0.0096, Max-Change = 0.0151, Max-Change = 0.0135, Max-Change = 0.0120, Max-Change = 0.0106, Max-Change = 0.0117, Max-Change = 0.0133, Max-Change = 0.0083, Max-Change = 0.0150, Max-Change = 0.0153, Max-Change = 0.0076, Max-Change = 0.0113, Max-Change = 0.0138, Max-Change = 0.0088, Max-Change = 0.0066, Max-Change = 0.0110, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0088, gam = 0.1057, Max-Change = 0.0072, gam = 0.0780, Max-Change = 0.0034, gam = 0.0629, Max-Change = 0.0035, gam = 0.0532, Max-Change = 0.0031, gam = 0.0464, Max-Change = 0.0023, gam = 0.0413, Max-Change = 0.0020, gam = 0.0374, Max-Change = 0.0019, gam = 0.0342, Max-Change = 0.0024, gam = 0.0316, Max-Change = 0.0020, gam = 0.0294, Max-Change = 0.0011, gam = 0.0276, Max-Change = 0.0026, gam = 0.0260, Max-Change = 0.0019, gam = 0.0246, Max-Change = 0.0019, gam = 0.0233, Max-Change = 0.0010, gam = 0.0222, Max-Change = 0.0012, gam = 0.0212, Max-Change = 0.0015, gam = 0.0203, Max-Change = 0.0012, gam = 0.0195, Max-Change = 0.0008, gam = 0.0188, Max-Change = 0.0008, gam = 0.0181, Max-Change = 0.0008
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod3)
#> 
#> Call:
#> mixedmirt(data = dat, covdata = covdata, model = 1, fixed = ~0 + 
#>     items, itemtype = "2PL", lr.random = ~1 | group)
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>    F1
#> F1  1
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $group
#>           COV_group
#> COV_group         1
#> 
anova(mod1b, mod3)
#>            AIC    SABIC       HQ      BIC    logLik     X2 df p
#> mod1b 23637.58 23664.24 23660.20 23699.19 -11807.79            
#> mod3  23571.39 23619.87 23612.52 23683.41 -11765.69 84.189  9 0

head(cbind(randef(mod3)$group, random_intercept))
#>         group random_intercept
#> G1  1.1060541       1.51178117
#> G2 -0.5373113       0.38984324
#> G3 -0.4789119      -0.62124058
#> G4 -2.3803749      -2.21469989
#> G5  0.6505437       1.12493092
#> G6 -0.6378008      -0.04493361

# }
```
