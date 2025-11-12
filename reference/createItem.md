# Create a user defined item with correct generic functions

Initializes the proper S4 class and methods necessary for
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md)
functions to use in estimation. To use the defined objects pass to the
`mirt(..., customItems = list())` command, and ensure that the classes
are properly labelled and unique in the list. Additionally, the input
`mirt(..., customItemsData = list())` can also be included to specify
additional item-level information to better recycle custom-item
definitions (e.g., for supplying varying Q-matrices), where the `list`
input must have the same length as the number of items. For further
examples regarding how this function can be used for fitting
unfolding-type models see Liu and Chalmers (2018).

## Usage

``` r
createItem(
  name,
  par,
  est,
  P,
  gr = NULL,
  hss = NULL,
  gen = NULL,
  lbound = NULL,
  ubound = NULL,
  derivType = "Richardson",
  derivType.hss = "Richardson",
  bytecompile = TRUE
)
```

## Arguments

- name:

  a character indicating the item class name to be defined

- par:

  a named vector of the starting values for the parameters

- est:

  a logical vector indicating which parameters should be freely
  estimated by default

- P:

  the probability trace function for all categories (first column is
  category 1, second category two, etc). First input contains a vector
  of all the item parameters, the second input must be a matrix called
  `Theta`, the third input must be the number of categories called
  `ncat`, and (optionally) a fourth argument termed `itemdata` may be
  included containing further users specification information. The last
  optional input is to be utilized within the estimation functions such
  as [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) via
  the list input `customItemsData` to more naturally recycle custom-item
  definitions. Therefore, these inputs must be of the form

  `function(par, Theta, ncat){...}`

  or

  `function(par, Theta, ncat, itemdata){...}`

  to be valid; however, the names of the arguements is not relavent.

  Finally, this function must return a `matrix` object of category
  probabilities, where the columns represent each respective category

- gr:

  gradient function (vector of first derivatives) of the log-likelihood
  used in estimation. The function must be of the form `gr(x, Theta)`,
  where `x` is the object defined by `createItem()` and `Theta` is a
  matrix of latent trait parameters. Tabulated (EM) or raw (MHRM) data
  are located in the `x@dat` slot, and are used to form the complete
  data log-likelihood. If not specified a numeric approximation will be
  used

- hss:

  Hessian function (matrix of second derivatives) of the log-likelihood
  used in estimation. If not specified a numeric approximation will be
  used (required for the MH-RM algorithm only). The input is identical
  to the `gr` argument

- gen:

  a function used when `GenRandomPars = TRUE` is passed to the
  estimation function to generate random starting values. Function must
  be of the form `function(object) ...` and must return a vector with
  properties equivalent to the `par` object. If NULL, parameters will
  remain at the defined starting values by default

- lbound:

  optional vector indicating the lower bounds of the parameters. If not
  specified then the bounds will be set to -Inf

- ubound:

  optional vector indicating the lower bounds of the parameters. If not
  specified then the bounds will be set to Inf

- derivType:

  if the `gr` term is not specified this type will be used to obtain the
  gradient numerically or symbolically. Default is the 'Richardson'
  extrapolation method; see
  [`numerical_deriv`](https://philchalmers.github.io/mirt/reference/numerical_deriv.md)
  for details and other options. If `'symbolic'` is supplied then the
  gradient will be computed using a symbolical approach (potentially the
  most accurate method, though may fail depending on how the `P`
  function was defined)

- derivType.hss:

  if the `hss` term is not specified this type will be used to obtain
  the Hessian numerically. Default is the 'Richardson' extrapolation
  method; see
  [`numerical_deriv`](https://philchalmers.github.io/mirt/reference/numerical_deriv.md)
  for details and other options. If `'symbolic'` is supplied then the
  Hessian will be computed using a symbolical approach (potentially the
  most accurate method, though may fail depending on how the `P`
  function was defined)

- bytecompile:

  logical; where applicable, byte compile the functions provided?
  Default is `TRUE` to provide

## Details

The [`summary()`](https://rdrr.io/r/base/summary.html) function will not
return proper standardized loadings since the function is not sure how
to handle them (no slopes could be defined at all!). Instead loadings of
.001 are filled in as place-holders.

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Liu, C.-W. and Chalmers, R. P. (2018). Fitting item response unfolding
models to Likert-scale data using mirt in R. *PLoS ONE, 13*, 5.
[doi:10.1371/journal.pone.0196292](https://doi.org/10.1371/journal.pone.0196292)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

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

# So, let's estimate it!
dat <- expand.table(LSAT7)
sv <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), pars = 'values')
tail(sv) #looks good
#>    group   item     class   name parnum value lbound ubound   est const nconst
#> 15   all Item.4      dich      g     15   0.0      0      1 FALSE  none   none
#> 16   all Item.4      dich      u     16   1.0      0      1 FALSE  none   none
#> 17   all Item.5    custom      a     17   0.5   -Inf    Inf  TRUE  none   none
#> 18   all Item.5    custom      b     18  -2.0   -Inf    Inf  TRUE  none   none
#> 19   all  GROUP GroupPars MEAN_1     19   0.0   -Inf    Inf FALSE  none   none
#> 20   all  GROUP GroupPars COV_11     20   1.0      0    Inf FALSE  none   none
#>    prior.type prior_1 prior_2
#> 15       none     NaN     NaN
#> 16       none     NaN     NaN
#> 17       none     NaN     NaN
#> 18       none     NaN     NaN
#> 19       none     NaN     NaN
#> 20       none     NaN     NaN
mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x))
coef(mod)
#> $Item.1
#>        a1     d g u
#> par 0.989 1.856 0 1
#> 
#> $Item.2
#>        a1     d g u
#> par 1.081 0.808 0 1
#> 
#> $Item.3
#>        a1     d g u
#> par 1.703 1.803 0 1
#> 
#> $Item.4
#>        a1     d g u
#> par 0.766 0.486 0 1
#> 
#> $Item.5
#>         a      b
#> par 0.737 -2.518
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> 
mod2 <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), method = 'MHRM')
coef(mod2)
#> $Item.1
#>        a1     d g u
#> par 0.965 1.842 0 1
#> 
#> $Item.2
#>        a1     d g u
#> par 1.093 0.809 0 1
#> 
#> $Item.3
#>        a1     d g u
#> par 1.747 1.822 0 1
#> 
#> $Item.4
#>        a1     d g u
#> par 0.779 0.487 0 1
#> 
#> $Item.5
#>         a      b
#> par 0.752 -2.474
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> 

# same definition as above, but using symbolic derivative computations
# (can be more accurate/stable)
xs <- createItem(name, par=par, est=est, P=P.old2PL, derivType = 'symbolic')
mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=xs))
coef(mod, simplify=TRUE)
#> $items
#>           a1     d  g  u     a      b
#> Item.1 0.989 1.856  0  1    NA     NA
#> Item.2 1.081 0.808  0  1    NA     NA
#> Item.3 1.703 1.803  0  1    NA     NA
#> Item.4 0.766 0.486  0  1    NA     NA
#> Item.5    NA    NA NA NA 0.737 -2.518
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 

# several secondary functions supported
M2(mod, calcNull=FALSE)
#>           M2 df     p RMSEA RMSEA_5 RMSEA_95 SRMSR
#> stats 11.936  5 0.036 0.037   0.009    0.065 0.032
itemfit(mod)
#>     item   S_X2 df.S_X2 RMSEA.S_X2 p.S_X2
#> 1 Item.1  4.750       2      0.037  0.093
#> 2 Item.2 14.441       2      0.079  0.001
#> 3 Item.3  1.266       2      0.000  0.531
#> 4 Item.4  5.241       2      0.040  0.073
#> 5 Item.5  0.941       2      0.000  0.625
fscores(mod, full.scores=FALSE)
#>       Item.1 Item.2 Item.3 Item.4 Item.5     F1 SE_F1
#>  [1,]      0      0      0      0      0 -1.870 0.693
#>  [2,]      0      0      0      0      1 -1.527 0.674
#>  [3,]      0      0      0      1      0 -1.514 0.673
#>  [4,]      0      0      0      1      1 -1.185 0.665
#>  [5,]      0      0      1      0      0 -1.096 0.665
#>  [6,]      0      0      1      0      1 -0.767 0.672
#>  [7,]      0      0      1      1      0 -0.754 0.673
#>  [8,]      0      0      1      1      1 -0.412 0.692
#>  [9,]      0      1      0      0      0 -1.372 0.668
#> [10,]      0      1      0      0      1 -1.045 0.666
#> [11,]      0      1      0      1      0 -1.032 0.666
#> [12,]      0      1      0      1      1 -0.702 0.675
#> [13,]      0      1      1      0      0 -0.610 0.680
#> [14,]      0      1      1      0      1 -0.258 0.704
#> [15,]      0      1      1      1      0 -0.244 0.705
#> [16,]      0      1      1      1      1  0.141 0.741
#> [17,]      1      0      0      0      0 -1.413 0.670
#> [18,]      1      0      0      0      1 -1.086 0.665
#> [19,]      1      0      0      1      0 -1.073 0.665
#> [20,]      1      0      0      1      1 -0.744 0.673
#> [21,]      1      0      1      0      0 -0.653 0.678
#> [22,]      1      0      1      0      1 -0.304 0.701
#> [23,]      1      0      1      1      0 -0.290 0.702
#> [24,]      1      0      1      1      1  0.090 0.736
#> [25,]      1      1      0      0      0 -0.933 0.667
#> [26,]      1      1      0      0      1 -0.600 0.680
#> [27,]      1      1      0      1      0 -0.587 0.681
#> [28,]      1      1      0      1      1 -0.233 0.706
#> [29,]      1      1      1      0      0 -0.132 0.715
#> [30,]      1      1      1      0      1  0.265 0.754
#> [31,]      1      1      1      1      0  0.282 0.755
#> [32,]      1      1      1      1      1  0.727 0.801
plot(mod)


# fit the same model, but specify gradient function explicitly (use of a browser() may be helpful)
gr <- function(x, Theta){
     # browser()
     a <- x@par[1]
     b <- x@par[2]
     P <- probtrace(x, Theta)
     PQ <- apply(P, 1, prod)
     r_P <- x@dat / P
     grad <- numeric(2)
     grad[2] <- sum(-a * PQ * (r_P[,2] - r_P[,1]))
     grad[1] <- sum((Theta - b) * PQ * (r_P[,2] - r_P[,1]))

     ## check with internal numerical form to be safe
     # numerical_deriv(x@par[x@est], mirt:::EML, obj=x, Theta=Theta)
     grad
}

x <- createItem(name, par=par, est=est, P=P.old2PL, gr=gr)
mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x))
coef(mod, simplify=TRUE)
#> $items
#>           a1     d  g  u     a      b
#> Item.1 0.989 1.856  0  1    NA     NA
#> Item.2 1.081 0.808  0  1    NA     NA
#> Item.3 1.703 1.803  0  1    NA     NA
#> Item.4 0.766 0.486  0  1    NA     NA
#> Item.5    NA    NA NA NA 0.737 -2.518
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>    F1
#> F1  1
#> 

### non-linear
name <- 'nonlin'
par <- c(a1 = .5, a2 = .1, d = 0)
est <- c(TRUE, TRUE, TRUE)
P.nonlin <- function(par,Theta, ncat=2){
     a1 <- par[1]
     a2 <- par[2]
     d <- par[3]
     P1 <- 1 / (1 + exp(-1*(a1*Theta + a2*Theta^2 + d)))
     cbind(1-P1, P1)
}

x2 <- createItem(name, par=par, est=est, P=P.nonlin)

mod <- mirt(dat, 1, c(rep('2PL',4), 'nonlin'), customItems=list(nonlin=x2))
coef(mod)
#> $Item.1
#>        a1     d g u
#> par 0.984 1.854 0 1
#> 
#> $Item.2
#>        a1     d g u
#> par 1.087 0.809 0 1
#> 
#> $Item.3
#>        a1     d g u
#> par 1.704 1.803 0 1
#> 
#> $Item.4
#>        a1     d g u
#> par 0.762 0.486 0 1
#> 
#> $Item.5
#>        a1    a2     d
#> par 0.806 0.065 1.818
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> 

### nominal response model (Bock 1972 version)
Tnom.dev <- function(ncat) {
   T <- matrix(1/ncat, ncat, ncat - 1)
   diag(T[-1, ]) <-  diag(T[-1, ]) - 1
   return(T)
}

name <- 'nom'
par <- c(alp=c(3,0,-3),gam=rep(.4,3))
est <- rep(TRUE, length(par))
P.nom <- function(par, Theta, ncat){
   alp <- par[1:(ncat-1)]
   gam <- par[ncat:length(par)]
   a <- Tnom.dev(ncat) %*% alp
   c <- Tnom.dev(ncat) %*% gam
   z <- matrix(0, nrow(Theta), ncat)
   for(i in 1:ncat)
       z[,i] <- a[i] * Theta + c[i]
   P <- exp(z) / rowSums(exp(z))
   P
}

nom1 <- createItem(name, par=par, est=est, P=P.nom)
nommod <- mirt(Science, 1, 'nom1', customItems=list(nom1=nom1))
coef(nommod)
#> $Comfort
#>       alp1   alp2   alp3   gam1   gam2   gam3
#> par -1.552 -2.015 -3.024 -3.639 -5.905 -4.533
#> 
#> $Work
#>      alp1   alp2   alp3   gam1   gam2   gam3
#> par -0.58 -1.262 -2.523 -1.464 -2.327 -0.326
#> 
#> $Future
#>       alp1 alp2   alp3   gam1   gam2  gam3
#> par -1.559 -3.8 -6.118 -3.676 -5.875 -3.96
#> 
#> $Benefit
#>       alp1   alp2   alp3   gam1   gam2   gam3
#> par -0.808 -1.358 -2.338 -2.145 -2.912 -1.622
#> 
#> $GroupPars
#>     MEAN_1 COV_11
#> par      0      1
#> 
Tnom.dev(4) %*% coef(nommod)[[1]][1:3] #a
#>             [,1]
#> [1,] -1.64770841
#> [2,] -0.09534806
#> [3,]  0.36680247
#> [4,]  1.37625400
Tnom.dev(4) %*% coef(nommod)[[1]][4:6] #d
#>            [,1]
#> [1,] -3.5191097
#> [2,]  0.1195514
#> [3,]  2.3861166
#> [4,]  1.0134416

# }
```
