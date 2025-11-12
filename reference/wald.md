# Wald statistics for mirt models

Compute a Wald test given an `L` vector or matrix of numeric contrasts.
Requires that the model information matrix be computed (by passing
`SE = TRUE` when estimating the model). Use `wald(model)` to observe how
the information matrix columns are named, especially if the estimated
model contains constrained parameters (e.g., 1PL).

## Usage

``` r
wald(object, L, C = NULL)
```

## Arguments

- object:

  estimated object from `mirt`, `bfactor`, `multipleGroup`, `mixedmirt`,
  or `mdirt`

- L:

  a coefficient matrix with dimensions `nconstrasts x npars.estimated`,
  or a character vector giving the hypothesis in symbolic form (syntax
  format borrowed from the `car` package; see `Details` below). Omitting
  this value will return the column names of the information matrix used
  to identify the (potentially constrained) parameters

- C:

  a constant vector of population parameters to be compared along side
  L, where `length(C) == row(L)`. By default a vector of 0's is
  constructed. Note that when using the syntax input for `L` this
  argument is ignored

  The following description is borrowed from `car` package documentation
  pertaining to the character vector input to the argument `L`: "The
  hypothesis matrix can be supplied as a numeric matrix (or vector), the
  rows of which specify linear combinations of the model coefficients,
  which are tested equal to the corresponding entries in the
  right-hand-side vector, which defaults to a vector of zeroes.

  Alternatively, the hypothesis can be specified symbolically as a
  character vector with one or more elements, each of which gives either
  a linear combination of coefficients, or a linear equation in the
  coefficients (i.e., with both a left and right side separated by an
  equals sign). Components of a linear expression or linear equation can
  consist of numeric constants, or numeric constants multiplying
  coefficient names (in which case the number precedes the coefficient,
  and may be separated from it by spaces or an asterisk); constants of 1
  or -1 may be omitted. Spaces are always optional. Components are
  separated by plus or minus signs. Newlines or tabs in hypotheses will
  be treated as spaces. See the examples below."

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

# View parnumber index
data(LSAT7)
data <- expand.table(LSAT7)
mod <- mirt(data, 1, SE = TRUE)
coef(mod)
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

# see how the information matrix relates to estimated parameters, and how it lines up
#   with the parameter index
(infonames <- wald(mod))
#>  a1.1   d.2  a1.5   d.6  a1.9  d.10 a1.13  d.14 a1.17  d.18 
#> 0.988 1.856 1.081 0.808 1.706 1.804 0.765 0.486 0.736 1.855 
index <- mod2values(mod)
index[index$est, ]
#>    group   item class name parnum value lbound ubound  est const nconst
#> 1    all Item.1  dich   a1      1 0.988   -Inf    Inf TRUE  none   none
#> 2    all Item.1  dich    d      2 1.856   -Inf    Inf TRUE  none   none
#> 5    all Item.2  dich   a1      5 1.081   -Inf    Inf TRUE  none   none
#> 6    all Item.2  dich    d      6 0.808   -Inf    Inf TRUE  none   none
#> 9    all Item.3  dich   a1      9 1.706   -Inf    Inf TRUE  none   none
#> 10   all Item.3  dich    d     10 1.804   -Inf    Inf TRUE  none   none
#> 13   all Item.4  dich   a1     13 0.765   -Inf    Inf TRUE  none   none
#> 14   all Item.4  dich    d     14 0.486   -Inf    Inf TRUE  none   none
#> 17   all Item.5  dich   a1     17 0.736   -Inf    Inf TRUE  none   none
#> 18   all Item.5  dich    d     18 1.855   -Inf    Inf TRUE  none   none
#>    prior.type prior_1 prior_2
#> 1        none     NaN     NaN
#> 2        none     NaN     NaN
#> 5        none     NaN     NaN
#> 6        none     NaN     NaN
#> 9        none     NaN     NaN
#> 10       none     NaN     NaN
#> 13       none     NaN     NaN
#> 14       none     NaN     NaN
#> 17       none     NaN     NaN
#> 18       none     NaN     NaN

# second item slope equal to 0?
L <- matrix(0, 1, 10)
L[1,3] <- 1
wald(mod, L)
#>    W df p
#> 1 41  1 0

# same as above using character syntax input
infonames
#>  a1.1   d.2  a1.5   d.6  a1.9  d.10 a1.13  d.14 a1.17  d.18 
#> 0.988 1.856 1.081 0.808 1.706 1.804 0.765 0.486 0.736 1.855 
wald(mod, "a1.5 = 0")
#>    W df p
#> 1 41  1 0

# simultaneously test equal factor slopes for item 1 and 2, and 4 and 5
L <- matrix(0, 2, 10)
L[1,1] <- L[2, 7] <- 1
L[1,3] <- L[2, 9] <- -1
L
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    1    0   -1    0    0    0    0    0    0     0
#> [2,]    0    0    0    0    0    0    1    0   -1     0
wald(mod, L)
#>       W df     p
#> 1 0.153  2 0.926

# Again, using more efficient syntax
infonames
#>  a1.1   d.2  a1.5   d.6  a1.9  d.10 a1.13  d.14 a1.17  d.18 
#> 0.988 1.856 1.081 0.808 1.706 1.804 0.765 0.486 0.736 1.855 
wald(mod, c("a1.1 = a1.5", "a1.13 = a1.17"))
#>       W df     p
#> 1 0.153  2 0.926

# log-Liklihood tests (requires estimating a new model)
cmodel <- 'theta = 1-5
           CONSTRAIN = (1,2, a1), (4,5, a1)'
mod2 <- mirt(data, cmodel)
# or, equivalently
#mod2 <- mirt(data, 1, constrain = list(c(1,5), c(13,17)))
anova(mod2, mod)
#>           AIC    SABIC       HQ      BIC    logLik    X2 df     p
#> mod2 5333.763 5347.616 5348.685 5373.025 -2658.881               
#> mod  5337.610 5354.927 5356.263 5386.688 -2658.805 0.152  2 0.927

#####
# test equality of means in multi-group model:
#    H0: (mu1 - mu2) = (mu3 - mu4)

set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
itemtype <- rep('2PL', nrow(a))
N <- 500
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .5)
dataset3 <- simdata(a, d, N, itemtype, mu = -1)
dataset4 <- simdata(a, d, N, itemtype, mu = -.5)
dat <- rbind(dataset1, dataset2, dataset3, dataset4)
group <- factor(rep(paste0('D', 1:4), each=N))
levels(group)
#> [1] "D1" "D2" "D3" "D4"
models <- 'F1 = 1-15'

# 3 means estimated
mod_free <- multipleGroup(dat, models, group = group, SE=TRUE,
                          invariance=c('slopes', 'intercepts', 'free_var','free_means'))
wald(mod_free) # obtain parameter names
#>   a1.1.63.125.187    d.2.64.126.188   a1.5.67.129.191    d.6.68.130.192 
#>             1.234             0.643             1.189            -0.573 
#>   a1.9.71.133.195   d.10.72.134.196  a1.13.75.137.199   d.14.76.138.200 
#>             0.981            -0.141             1.058             0.864 
#>  a1.17.79.141.203   d.18.80.142.204  a1.21.83.145.207   d.22.84.146.208 
#>             1.146             0.248             0.499             0.537 
#>  a1.25.87.149.211   d.26.88.150.212  a1.29.91.153.215   d.30.92.154.216 
#>             1.274             1.204             0.997            -0.368 
#>  a1.33.95.157.219   d.34.96.158.220  a1.37.99.161.223  d.38.100.162.224 
#>             0.899            -1.032             0.735            -1.087 
#> a1.41.103.165.227  d.42.104.166.228 a1.45.107.169.231  d.46.108.170.232 
#>             1.031             1.345             1.603            -0.160 
#> a1.49.111.173.235  d.50.112.174.236 a1.53.115.177.239  d.54.116.178.240 
#>             1.272             0.578             1.228             0.501 
#> a1.57.119.181.243  d.58.120.182.244        MEAN_1.123        COV_11.124 
#>             0.850            -0.080             0.382             0.865 
#>        MEAN_1.185        COV_11.186        MEAN_1.247        COV_11.248 
#>            -1.062             0.864            -0.594             0.936 
# View(mod2values(mod_free))

# reference group mean = 0 by default
wald(mod_free, c("0 - MEAN_1.123 = MEAN_1.185 - MEAN_1.247"))
#>       W df     p
#> 1 0.714  1 0.398


# }
```
