# Simulate response patterns

Simulates response patterns for compensatory and noncompensatory MIRT
models from multivariate normally distributed factor (\\\theta\\)
scores, or from a user input matrix of \\\theta\\'s.

## Usage

``` r
simdata(
  a,
  d,
  N,
  itemtype,
  sigma = NULL,
  mu = NULL,
  guess = 0,
  upper = 1,
  rho = NULL,
  nominal = NULL,
  t = NULL,
  Theta = NULL,
  gpcm_mats = list(),
  returnList = FALSE,
  model = NULL,
  equal.K = TRUE,
  which.items = NULL,
  mins = 0,
  lca_cats = NULL,
  prob.list = NULL
)
```

## Arguments

- a:

  a matrix/vector of slope parameters. If slopes are to be constrained
  to zero then use `NA` or simply set them equal to 0

- d:

  a matrix/vector of intercepts. The matrix should have as many columns
  as the item with the largest number of categories, and filled empty
  locations with `NA`. When a vector is used the test is assumed to
  consist only of dichotomous items (because only one intercept per item
  is provided). When `itemtype = 'lca'` intercepts will not be used

- N:

  sample size

- itemtype:

  a character vector of length `nrow(a)` (or 1, if all the item types
  are the same) specifying the type of items to simulate. Inputs can
  either be the same as the inputs found in the `itemtype` argument in
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) or the
  internal classes defined by the package. Typical `itemtype` inputs
  that are passed to
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) are
  used then these will be converted into the respective internal classes
  automatically.

  If the internal class of the object is specified instead, the inputs
  can be
  `'dich', 'graded', 'gpcm', 'sequential', 'nominal', 'nestlogit', 'partcomp', 'gumm'`,
  `'lca'`, and all the models under the Luo (2001) family (see
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)), for
  dichotomous, graded, generalized partial credit, sequential, nominal,
  nested logit, partially compensatory, generalized graded unfolding
  model, latent class analysis model, and ordered unfolding models. Note
  that for the gpcm, nominal, and nested logit models there should be as
  many parameters as desired categories, however to parametrized them
  for meaningful interpretation the first category intercept should
  equal 0 for these models (second column for `'nestlogit'`, since first
  column is for the correct item traceline). For nested logit models the
  'correct' category is always the lowest category (i.e., == 1). It may
  be helpful to use
  [`mod2values`](https://philchalmers.github.io/mirt/reference/mod2values.md)
  on data-sets that have already been estimated to understand the
  itemtypes more intimately

- sigma:

  a covariance matrix of the underlying distribution. Default is the
  identity matrix. Used when `Theta` is not supplied

- mu:

  a mean vector of the underlying distribution. Default is a vector of
  zeros. Used when `Theta` is not supplied

- guess:

  a vector of guessing parameters for each item; only applicable for
  dichotomous items. Must be either a scalar value that will affect all
  of the dichotomous items, or a vector with as many values as to be
  simulated items

- upper:

  same as `guess`, but for upper bound parameters

- rho:

  a matrix of `rho` values to be used for the Lui (2001) family of
  ordered unfolding models (see
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)) to
  control the latitude of acceptance. All values must be positive

- nominal:

  a matrix of specific item category slopes for nominal models. Should
  be the dimensions as the intercept specification with one less column,
  with `NA` in locations where not applicable. Note that during
  estimation the first slope will be constrained to 0 and the last will
  be constrained to the number of categories minus 1, so it is best to
  set these as the values for the first and last categories as well

- t:

  matrix of t-values for the 'ggum' itemtype, where each row corresponds
  to a given item. Also determines the number of categories, where `NA`
  can be used for non-applicable categories

- Theta:

  a user specified matrix of the underlying ability parameters, where
  `nrow(Theta) == N` and `ncol(Theta) == ncol(a)`. When this is supplied
  the `N` input is not required

- gpcm_mats:

  a list of matrices specifying the scoring scheme for generalized
  partial credit models (see
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  details)

- returnList:

  logical; return a list containing the data, item objects defined by
  `mirt` containing the population parameters and item structure, and
  the latent trait matrix `Theta`? Default is FALSE

- model:

  a single group object, typically returned by functions such as
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) or
  [`bfactor`](https://philchalmers.github.io/mirt/reference/bfactor.md).
  Supplying this will render all other parameter elements (excluding the
  `Theta`, `N`, `mu`, and `sigma` inputs) redundant (unless explicitly
  provided). This input can therefore be used to create parametric
  bootstrap data whereby plausible data implied by the estimated model
  can be generated and evaluated

- equal.K:

  logical; when a `model` input is supplied, should the generated data
  contain the same number of categories as the original data indicated
  by `extract.mirt(model, 'K')`? Default is TRUE, which will redrawn
  data until this condition is satisfied

- which.items:

  an integer vector used to indicate which items to simulate when a
  `model` input is included. Default simulates all items

- mins:

  an integer vector (or single value to be used for each item)
  indicating what the lowest category should be. If `model` is supplied
  then this will be extracted from `slot(mod, 'Data')$mins`, otherwise
  the default is 0

- lca_cats:

  a vector indicating how many categories each lca item should have. If
  not supplied then it is assumed that 2 categories should be generated
  for each item

- prob.list:

  an optional list containing matrix/data.frames of probabilities values
  for each category to be simulated. This is useful when creating
  customized probability functions to be sampled from

## Details

Returns a data matrix simulated from the parameters, or a list
containing the data, item objects, and Theta matrix.

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Reckase, M. D. (2009). *Multidimensional Item Response Theory*. New
York: Springer.

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
### Parameters from Reckase (2009), p. 153

set.seed(1234)

a <- matrix(c(
 .7471, .0250, .1428,
 .4595, .0097, .0692,
 .8613, .0067, .4040,
1.0141, .0080, .0470,
 .5521, .0204, .1482,
1.3547, .0064, .5362,
1.3761, .0861, .4676,
 .8525, .0383, .2574,
1.0113, .0055, .2024,
 .9212, .0119, .3044,
 .0026, .0119, .8036,
 .0008, .1905,1.1945,
 .0575, .0853, .7077,
 .0182, .3307,2.1414,
 .0256, .0478, .8551,
 .0246, .1496, .9348,
 .0262, .2872,1.3561,
 .0038, .2229, .8993,
 .0039, .4720, .7318,
 .0068, .0949, .6416,
 .3073, .9704, .0031,
 .1819, .4980, .0020,
 .4115,1.1136, .2008,
 .1536,1.7251, .0345,
 .1530, .6688, .0020,
 .2890,1.2419, .0220,
 .1341,1.4882, .0050,
 .0524, .4754, .0012,
 .2139, .4612, .0063,
 .1761,1.1200, .0870),30,3,byrow=TRUE)*1.702

d <- matrix(c(.1826,-.1924,-.4656,-.4336,-.4428,-.5845,-1.0403,
  .6431,.0122,.0912,.8082,-.1867,.4533,-1.8398,.4139,
  -.3004,-.1824,.5125,1.1342,.0230,.6172,-.1955,-.3668,
  -1.7590,-.2434,.4925,-.3410,.2896,.006,.0329),ncol=1)*1.702

mu <- c(-.4, -.7, .1)
sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)

dataset1 <- simdata(a, d, 2000, itemtype = '2PL')
dataset2 <- simdata(a, d, 2000, itemtype = '2PL', mu = mu, sigma = sigma)

#mod <- mirt(dataset1, 3, method = 'MHRM')
#coef(mod)

# \donttest{

### Unidimensional graded response model with 5 categories each

a <- matrix(rlnorm(20,.2,.3))

# for the graded model, ensure that there is enough space between the intercepts,
# otherwise closer categories will not be selected often (minimum distance of 0.3 here)
diffs <- t(apply(matrix(runif(20*4, .3, 1), 20), 1, cumsum))
diffs <- -(diffs - rowMeans(diffs))
d <- diffs + rnorm(20)

dat <- simdata(a, d, 500, itemtype = 'graded')
# mod <- mirt(dat, 1)

### An example of a mixed item, bifactor loadings pattern with correlated specific factors

a <- matrix(c(
.8,.4,NA,
.4,.4,NA,
.7,.4,NA,
.8,NA,.4,
.4,NA,.4,
.7,NA,.4),ncol=3,byrow=TRUE)

d <- matrix(c(
-1.0,NA,NA,
 1.5,NA,NA,
 0.0,NA,NA,
0.0,-1.0,1.5,  #the first 0 here is the recommended constraint for nominal
0.0,1.0,-1, #the first 0 here is the recommended constraint for gpcm
2.0,0.0,NA),ncol=3,byrow=TRUE)

nominal <- matrix(NA, nrow(d), ncol(d))
# the first 0 and last (ncat - 1) = 2 values are the recommended constraints
nominal[4, ] <- c(0,1.2,2)

sigma <- diag(3)
sigma[2,3] <- sigma[3,2] <- .25
items <- c('2PL','2PL','2PL','nominal','gpcm','graded')

dataset <- simdata(a,d,2000,items,sigma=sigma,nominal=nominal)

#mod <- bfactor(dataset, c(1,1,1,2,2,2), itemtype=c(rep('2PL', 3), 'nominal', 'gpcm','graded'))
#coef(mod)

#### Convert standardized factor loadings to slopes

F2a <- function(F, D=1.702){
    h2 <- rowSums(F^2)
    a <- (F / sqrt(1 - h2)) * D
    a
}

(F <- matrix(c(rep(.7, 5), rep(.5,5))))
#>       [,1]
#>  [1,]  0.7
#>  [2,]  0.7
#>  [3,]  0.7
#>  [4,]  0.7
#>  [5,]  0.7
#>  [6,]  0.5
#>  [7,]  0.5
#>  [8,]  0.5
#>  [9,]  0.5
#> [10,]  0.5
(a <- F2a(F))
#>            [,1]
#>  [1,] 1.6682937
#>  [2,] 1.6682937
#>  [3,] 1.6682937
#>  [4,] 1.6682937
#>  [5,] 1.6682937
#>  [6,] 0.9826502
#>  [7,] 0.9826502
#>  [8,] 0.9826502
#>  [9,] 0.9826502
#> [10,] 0.9826502
d <- rnorm(10)

dat <- simdata(a, d, 5000, itemtype = '2PL')
mod <- mirt(dat, 1)
coef(mod, simplify=TRUE)$items
#>                a1           d g u
#> Item_1  1.5912124  0.09329250 0 1
#> Item_2  1.6611479  0.77704003 0 1
#> Item_3  1.6100289  0.04348744 0 1
#> Item_4  1.5251147 -0.39866322 0 1
#> Item_5  1.7126206 -0.69888431 0 1
#> Item_6  0.9352647 -1.71863199 0 1
#> Item_7  0.9211469  0.73522267 0 1
#> Item_8  0.9893392 -1.25501247 0 1
#> Item_9  0.9999014  0.41687594 0 1
#> Item_10 0.9416579 -1.00423353 0 1
summary(mod)
#>            F1    h2
#> Item_1  0.683 0.466
#> Item_2  0.698 0.488
#> Item_3  0.687 0.472
#> Item_4  0.667 0.445
#> Item_5  0.709 0.503
#> Item_6  0.482 0.232
#> Item_7  0.476 0.227
#> Item_8  0.503 0.253
#> Item_9  0.507 0.257
#> Item_10 0.484 0.234
#> 
#> SS loadings:  3.577 
#> Proportion Var:  0.358 
#> 
#> Factor correlations: 
#> 
#>    F1
#> F1  1

mod2 <- mirt(dat, 'F1 = 1-10
                   CONSTRAIN = (1-5, a1), (6-10, a1)')
summary(mod2)
#>            F1    h2
#> Item_1  0.689 0.475
#> Item_2  0.689 0.475
#> Item_3  0.689 0.475
#> Item_4  0.689 0.475
#> Item_5  0.689 0.475
#> Item_6  0.491 0.241
#> Item_7  0.491 0.241
#> Item_8  0.491 0.241
#> Item_9  0.491 0.241
#> Item_10 0.491 0.241
#> 
#> SS loadings:  3.576 
#> Proportion Var:  0.358 
#> 
#> Factor correlations: 
#> 
#>    F1
#> F1  1
anova(mod2, mod)
#>           AIC    SABIC       HQ      BIC    logLik    X2 df   p
#> mod2 58058.43 58098.51 58085.85 58136.64 -29017.22             
#> mod  58068.01 58134.80 58113.69 58198.35 -29014.01 6.425  8 0.6

#### Convert classical 3PL paramerization into slope-intercept form
nitems <- 50
as <- rlnorm(nitems, .2, .2)
bs <- rnorm(nitems, 0, 1)
gs <- rbeta(nitems, 5, 17)

# convert first item (only intercepts differ in resulting transformation)
traditional2mirt(c('a'=as[1], 'b'=bs[1], 'g'=gs[1], 'u'=1), cls='3PL')
#>         a1          d          g          u 
#>  1.2795115 -0.1107008  0.2525144  1.0000000 

# convert all difficulties to intercepts
ds <- numeric(nitems)
for(i in 1:nitems)
   ds[i] <- traditional2mirt(c('a'=as[i], 'b'=bs[i], 'g'=gs[i], 'u'=1),
                             cls='3PL')[2]

dat <- simdata(as, ds, N=5000, guess=gs, itemtype = '3PL')

# estimate with beta prior for guessing parameters
# mod <- mirt(dat, model="Theta = 1-50
#                         PRIOR = (1-50, g, expbeta, 5, 17)", itemtype = '3PL')
# coef(mod, simplify=TRUE, IRTpars=TRUE)$items
# data.frame(as, bs, gs, us=1)


#### Unidimensional nonlinear factor pattern

theta <- rnorm(2000)
Theta <- cbind(theta,theta^2)

a <- matrix(c(
.8,.4,
.4,.4,
.7,.4,
.8,NA,
.4,NA,
.7,NA),ncol=2,byrow=TRUE)
d <- matrix(rnorm(6))
itemtype <- rep('2PL',6)

nonlindata <- simdata(a=a, d=d, itemtype=itemtype, Theta=Theta)

#model <- '
#F1 = 1-6
#(F1 * F1) = 1-3'
#mod <- mirt(nonlindata, model)
#coef(mod)

#### 2PLNRM model for item 4 (with 4 categories), 2PL otherwise

a <- matrix(rlnorm(4,0,.2))

# first column of item 4 is the intercept for the correct category of 2PL model,
#    otherwise nominal model configuration
d <- matrix(c(
-1.0,NA,NA,NA,
 1.5,NA,NA,NA,
 0.0,NA,NA,NA,
 1, 0.0,-0.5,0.5),ncol=4,byrow=TRUE)

nominal <- matrix(NA, nrow(d), ncol(d))
nominal[4, ] <- c(NA,0,.5,.6)

items <- c(rep('2PL',3),'nestlogit')

dataset <- simdata(a,d,2000,items,nominal=nominal)

#mod <- mirt(dataset, 1, itemtype = c('2PL', '2PL', '2PL', '2PLNRM'), key=c(NA,NA,NA,0))
#coef(mod)
#itemplot(mod,4)

# return list of simulation parameters
listobj <- simdata(a,d,2000,items,nominal=nominal, returnList=TRUE)
str(listobj)
#> List of 3
#>  $ itemobjects:List of 4
#>   ..$ :Formal class 'dich' [package "mirt"] with 23 slots
#>   .. .. ..@ par          : num [1:4] 0.863 -1 -999 999
#>   .. .. .. ..- attr(*, "na.action")= 'omit' int [1:3] 3 4 5
#>   .. .. ..@ SEpar        : num(0) 
#>   .. .. ..@ parnames     : chr(0) 
#>   .. .. ..@ est          : logi(0) 
#>   .. .. ..@ dps          :function ()  
#>   .. .. ..@ dps2         :function ()  
#>   .. .. ..@ constr       : logi(0) 
#>   .. .. ..@ itemclass    : int(0) 
#>   .. .. ..@ parnum       : num(0) 
#>   .. .. ..@ nfact        : int 1
#>   .. .. ..@ nfixedeffects: num(0) 
#>   .. .. ..@ fixed.design : num[0 , 0 ] 
#>   .. .. ..@ dat          : num[0 , 0 ] 
#>   .. .. ..@ ncat         : int 2
#>   .. .. ..@ gradient     : num(0) 
#>   .. .. ..@ hessian      : num[0 , 0 ] 
#>   .. .. ..@ itemtrace    : num[0 , 0 ] 
#>   .. .. ..@ lbound       : num(0) 
#>   .. .. ..@ ubound       : num(0) 
#>   .. .. ..@ any.prior    : logi(0) 
#>   .. .. ..@ prior.type   : int(0) 
#>   .. .. ..@ prior_1      : num(0) 
#>   .. .. ..@ prior_2      : num(0) 
#>   ..$ :Formal class 'dich' [package "mirt"] with 23 slots
#>   .. .. ..@ par          : num [1:4] 1.25 1.5 -999 999
#>   .. .. .. ..- attr(*, "na.action")= 'omit' int [1:3] 3 4 5
#>   .. .. ..@ SEpar        : num(0) 
#>   .. .. ..@ parnames     : chr(0) 
#>   .. .. ..@ est          : logi(0) 
#>   .. .. ..@ dps          :function ()  
#>   .. .. ..@ dps2         :function ()  
#>   .. .. ..@ constr       : logi(0) 
#>   .. .. ..@ itemclass    : int(0) 
#>   .. .. ..@ parnum       : num(0) 
#>   .. .. ..@ nfact        : int 1
#>   .. .. ..@ nfixedeffects: num(0) 
#>   .. .. ..@ fixed.design : num[0 , 0 ] 
#>   .. .. ..@ dat          : num[0 , 0 ] 
#>   .. .. ..@ ncat         : int 2
#>   .. .. ..@ gradient     : num(0) 
#>   .. .. ..@ hessian      : num[0 , 0 ] 
#>   .. .. ..@ itemtrace    : num[0 , 0 ] 
#>   .. .. ..@ lbound       : num(0) 
#>   .. .. ..@ ubound       : num(0) 
#>   .. .. ..@ any.prior    : logi(0) 
#>   .. .. ..@ prior.type   : int(0) 
#>   .. .. ..@ prior_1      : num(0) 
#>   .. .. ..@ prior_2      : num(0) 
#>   ..$ :Formal class 'dich' [package "mirt"] with 23 slots
#>   .. .. ..@ par          : num [1:4] 1.02 0 -999 999
#>   .. .. .. ..- attr(*, "na.action")= 'omit' int [1:3] 3 4 5
#>   .. .. ..@ SEpar        : num(0) 
#>   .. .. ..@ parnames     : chr(0) 
#>   .. .. ..@ est          : logi(0) 
#>   .. .. ..@ dps          :function ()  
#>   .. .. ..@ dps2         :function ()  
#>   .. .. ..@ constr       : logi(0) 
#>   .. .. ..@ itemclass    : int(0) 
#>   .. .. ..@ parnum       : num(0) 
#>   .. .. ..@ nfact        : int 1
#>   .. .. ..@ nfixedeffects: num(0) 
#>   .. .. ..@ fixed.design : num[0 , 0 ] 
#>   .. .. ..@ dat          : num[0 , 0 ] 
#>   .. .. ..@ ncat         : int 2
#>   .. .. ..@ gradient     : num(0) 
#>   .. .. ..@ hessian      : num[0 , 0 ] 
#>   .. .. ..@ itemtrace    : num[0 , 0 ] 
#>   .. .. ..@ lbound       : num(0) 
#>   .. .. ..@ ubound       : num(0) 
#>   .. .. ..@ any.prior    : logi(0) 
#>   .. .. ..@ prior.type   : int(0) 
#>   .. .. ..@ prior_1      : num(0) 
#>   .. .. ..@ prior_2      : num(0) 
#>   ..$ :Formal class 'nestlogit' [package "mirt"] with 24 slots
#>   .. .. ..@ correctcat   : int 1
#>   .. .. ..@ par          : num [1:10] 0.771 1 -999 999 0 ...
#>   .. .. ..@ SEpar        : num(0) 
#>   .. .. ..@ parnames     : chr(0) 
#>   .. .. ..@ est          : logi(0) 
#>   .. .. ..@ dps          :function ()  
#>   .. .. ..@ dps2         :function ()  
#>   .. .. ..@ constr       : logi(0) 
#>   .. .. ..@ itemclass    : int(0) 
#>   .. .. ..@ parnum       : num(0) 
#>   .. .. ..@ nfact        : int 1
#>   .. .. ..@ nfixedeffects: num(0) 
#>   .. .. ..@ fixed.design : num[0 , 0 ] 
#>   .. .. ..@ dat          : num[0 , 0 ] 
#>   .. .. ..@ ncat         : int 4
#>   .. .. ..@ gradient     : num(0) 
#>   .. .. ..@ hessian      : num[0 , 0 ] 
#>   .. .. ..@ itemtrace    : num[0 , 0 ] 
#>   .. .. ..@ lbound       : num(0) 
#>   .. .. ..@ ubound       : num(0) 
#>   .. .. ..@ any.prior    : logi(0) 
#>   .. .. ..@ prior.type   : int(0) 
#>   .. .. ..@ prior_1      : num(0) 
#>   .. .. ..@ prior_2      : num(0) 
#>  $ data       : num [1:2000, 1:4] 1 1 1 0 0 1 0 1 0 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:4] "Item_1" "Item_2" "Item_3" "Item_4"
#>  $ Theta      : num [1:2000, 1] -0.399 1.931 1.621 -0.424 1.343 ...

# generate dataset from converged model
mod <- mirt(Science, 1, itemtype = c(rep('gpcm', 3), 'nominal'))
sim <- simdata(model=mod, N=1000)
head(sim)
#>      Comfort Work Future Benefit
#> [1,]       4    3      3       4
#> [2,]       3    2      3       3
#> [3,]       3    2      4       2
#> [4,]       4    3      4       3
#> [5,]       3    3      2       3
#> [6,]       2    3      3       2

Theta <- matrix(rnorm(100))
sim <- simdata(model=mod, Theta=Theta)
head(sim)
#>      Comfort Work Future Benefit
#> [1,]       3    3      4       4
#> [2,]       4    2      3       2
#> [3,]       3    2      3       2
#> [4,]       3    3      3       3
#> [5,]       3    3      3       3
#> [6,]       3    3      4       3

# alternatively, define a suitable object with functions from the mirtCAT package
# help(generate.mirt_object)
library(mirtCAT)
#> Loading required package: shiny

nitems <- 50
a1 <- rlnorm(nitems, .2,.2)
d <- rnorm(nitems)
g <- rbeta(nitems, 20, 80)
pars <- data.frame(a1=a1, d=d, g=g)
head(pars)
#>          a1          d         g
#> 1 1.1271598 -0.6479570 0.2637693
#> 2 1.3539630 -1.4832170 0.1903528
#> 3 1.1640507  0.3390045 0.1782070
#> 4 0.9568311 -1.3578148 0.1389409
#> 5 1.0394247  0.3835373 0.2361323
#> 6 1.1479152  1.3522255 0.1990177

obj <- generate.mirt_object(pars, '3PL')
dat <- simdata(N=200, model=obj)

#### 10 item GGUMs test with 4 categories each
a <- rlnorm(10, .2, .2)
b <- rnorm(10) #passed to d= input, but used as the b parameters
diffs <- t(apply(matrix(runif(10*3, .3, 1), 10), 1, cumsum))
t <- -(diffs - rowMeans(diffs))

dat <- simdata(a, b, 1000, 'ggum', t=t)
apply(dat, 2, table)
#>   Item_1 Item_2 Item_3 Item_4 Item_5 Item_6 Item_7 Item_8 Item_9 Item_10
#> 0    477    496    674    332    438    417    558    458    475     513
#> 1    299    309    218    388    364    365    288    320    314     305
#> 2    153    155     79    232    158    175    109    167    162     138
#> 3     71     40     29     48     40     43     45     55     49      44
# mod <- mirt(dat, 1, 'ggum')
# coef(mod)

### 10 items with the hyperbolic cosine model with differing category lengths
a <- matrix(1, 10)
d <- rnorm(10)
rho <- matrix(1:2, nrow=10, ncol=2, byrow=TRUE)
rho[1:2,2] <- NA   # first two items have K=2 categories

dat <- simdata(a, d, 1000, 'hcm', rho=rho)
itemstats(dat)
#> $overall
#>     N mean_total.score sd_total.score ave.r  sd.r alpha SEM.alpha
#>  1000            8.522          2.709 0.015 0.212 0.119     2.543
#> 
#> $itemstats
#>            N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1  1000 2 0.475 0.500   0.128        -0.057       0.146
#> Item_2  1000 2 0.343 0.475   0.352         0.186       0.056
#> Item_3  1000 3 0.574 0.796   0.367         0.079       0.082
#> Item_4  1000 3 0.989 0.895   0.531         0.231      -0.038
#> Item_5  1000 3 1.159 0.867   0.163        -0.158       0.234
#> Item_6  1000 3 0.911 0.874   0.511         0.214      -0.021
#> Item_7  1000 3 1.001 0.884   0.525         0.228      -0.034
#> Item_8  1000 3 0.965 0.896   0.016        -0.300       0.323
#> Item_9  1000 3 0.886 0.895   0.517         0.213      -0.023
#> Item_10 1000 3 1.219 0.861   0.210        -0.110       0.204
#> 
#> $proportions
#>             0     1     2
#> Item_1  0.525 0.475    NA
#> Item_2  0.657 0.343    NA
#> Item_3  0.620 0.186 0.194
#> Item_4  0.406 0.199 0.395
#> Item_5  0.309 0.223 0.468
#> Item_6  0.430 0.229 0.341
#> Item_7  0.390 0.219 0.391
#> Item_8  0.419 0.197 0.384
#> Item_9  0.464 0.186 0.350
#> Item_10 0.285 0.211 0.504
#> 
# mod <- mirt(dat, 1, 'hcm')
# list(est=coef(mod, simplify=TRUE)$items, pop=cbind(a, d, log(rho)))


######
# prob.list example

# custom probability function that returns a matrix
fun <- function(a, b, theta){
    P <- 1 / (1 + exp(-a * (theta-b)))
    cbind(1-P, P)
}

set.seed(1)
theta <- matrix(rnorm(100))
prob.list <- list()
nitems <- 5
a <- rlnorm(nitems, .2, .2); b <- rnorm(nitems, 0, 1/2)
for(i in 1:nitems) prob.list[[i]] <- fun(a[i], b[i], theta)
str(prob.list)
#> List of 5
#>  $ : num [1:100, 1:2] 0.836 0.68 0.865 0.317 0.645 ...
#>  $ : num [1:100, 1:2] 0.771 0.554 0.813 0.179 0.509 ...
#>  $ : num [1:100, 1:2] 0.75 0.569 0.788 0.239 0.532 ...
#>  $ : num [1:100, 1:2] 0.737 0.503 0.785 0.146 0.457 ...
#>  $ : num [1:100, 1:2] 0.828 0.669 0.858 0.308 0.634 ...

dat <- simdata(prob.list=prob.list)
head(dat)
#>      Item_1 Item_2 Item_3 Item_4 Item_5
#> [1,]      0      0      0      1      1
#> [2,]      0      0      0      1      0
#> [3,]      0      0      0      0      0
#> [4,]      1      1      0      0      0
#> [5,]      1      1      0      0      0
#> [6,]      0      1      0      1      0

# prob.list input is useful when defining custom items as well
name <- 'old2PL'
par <- c(a = .5, b = -2)
est <- c(TRUE, TRUE)
P.old2PL <- function(par,Theta, ncat){
     a <- par[1]
     b <- par[2]
     P1 <- 1 / (1 + exp(-1*a*(Theta - b)))
     cbind(1-P1, P1)
}

x <- createItem(name, par=par, est=est, P=P.old2PL)

prob.list[[1]] <- x@P(x@par, theta)


# }
```
