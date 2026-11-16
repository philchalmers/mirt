# Differential test functioning statistics

Function performs various omnibus differential test functioning
procedures on an object estimated with
[`multipleGroup()`](https://philchalmers.github.io/mirt/reference/multipleGroup.md).
If the latent means/covariances are suspected to differ then the input
object should contain a set of 'anchor' items to ensure that only
differential test features are being detected rather than group
differences. Returns signed (average area above and below) and unsigned
(total area) statistics, with descriptives such as the percent average
bias between group total scores for each statistic. If a grid of Theta
values is passed, these can be evaluated as well to determine specific
DTF location effects. For best results, the baseline model should
contain a set of 'anchor' items and have freely estimated
hyper-parameters in the focal groups. See
[`DIF`](https://philchalmers.github.io/mirt/reference/DIF.md) for
details.

## Usage

``` r
DTF(
  mod,
  draws = NULL,
  CI = 0.95,
  npts = 1000,
  theta_lim = c(-6, 6),
  Theta_nodes = NULL,
  plot = "none",
  auto.key = list(space = "right", points = FALSE, lines = TRUE),
  ...
)
```

## Arguments

- mod:

  a multipleGroup object which estimated only 2 groups

- draws:

  a number indicating how many draws to take to form a suitable multiple
  imputation estimate of the expected test scores (usually 100 or more).
  Returns a list containing the imputation distribution and null
  hypothesis test for the sDTF statistic

- CI:

  range of confidence interval when using draws input

- npts:

  number of points to use in the integration. Default is 1000

- theta_lim:

  lower and upper limits of the latent trait (theta) to be evaluated,
  and is used in conjunction with `npts`

- Theta_nodes:

  an optional matrix of Theta values to be evaluated in the draws for
  the sDTF statistic. However, these values are not averaged across, and
  instead give the bootstrap confidence intervals at the respective
  Theta nodes. Useful when following up a large uDTF/sDTF statistic to
  determine where the difference between the test curves are large
  (while still accounting for sampling variability). Returns a matrix
  with observed variability

- plot:

  a character vector indicating which plot to draw. Possible values are
  'none', 'func' for the test score functions, and 'sDTF' for the
  evaluated sDTF values across the integration grid. Each plot is drawn
  with imputed confidence envelopes

- auto.key:

  logical; automatically generate key in lattice plot?

- ...:

  additional arguments to be passed to `lattice` and `boot`

## References

Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

Chalmers, R. P., Counsell, A., and Flora, D. B. (2016). It might not
make a big DIF: Improved Differential Test Functioning statistics that
account for sampling variability. *Educational and Psychological
Measurement, 76*, 114-140.
[doi:10.1177/0013164415584576](https://doi.org/10.1177/0013164415584576)

## See also

[`multipleGroup`](https://philchalmers.github.io/mirt/reference/multipleGroup.md),
[`DIF`](https://philchalmers.github.io/mirt/reference/DIF.md)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
set.seed(1234)
n <- 30
N <- 500

# only first 5 items as anchors
model <- 'F = 1-30
          CONSTRAINB = (1-5, a1), (1-5, d)'

a <- matrix(1, n)
d <- matrix(rnorm(n), n)
group <- c(rep('Group_1', N), rep('Group_2', N))

## -------------
# groups completely equal
dat1 <- simdata(a, d, N, itemtype = '2PL')
dat2 <- simdata(a, d, N, itemtype = '2PL')
dat <- rbind(dat1, dat2)
mod <- multipleGroup(dat, model, group=group, SE=TRUE,
                     invariance=c('free_means', 'free_var'))
plot(mod)


DTF(mod)
#>    sDTF.score sDTF(%).score    uDTF.score uDTF(%).score 
#>     0.1743929     0.5813098     0.1771135     0.5903783 
if(interactive()) mirtCluster()
DTF(mod, draws = 1000) #95% C.I. for sDTF containing 0. uDTF is very small
#> $observed
#>    sDTF.score sDTF(%).score    uDTF.score uDTF(%).score 
#>     0.1743929     0.5813098     0.1771135     0.5903783 
#> 
#> $CIs
#>         sDTF.score sDTF(%).score uDTF.score uDTF(%).score
#> CI_97.5  0.4730704     1.5769014 0.76493827     2.5497942
#> CI_2.5  -0.1232603    -0.4108676 0.08069275     0.2689758
#> 
#> $tests
#> P(sDTF.score = 0) 
#>         0.2719623 
#> 
DTF(mod, draws = 1000, plot='sDTF') #sDTF 95% C.I.'s across Theta always include 0


## -------------
## random slopes and intercepts for 15 items, and latent mean difference
##    (no systematic DTF should exist, but DIF will be present)
set.seed(1234)
dat1 <- simdata(a, d, N, itemtype = '2PL', mu=.50, sigma=matrix(1.5))
dat2 <- simdata(a + c(numeric(15), runif(n-15, -.2, .2)),
                d + c(numeric(15), runif(n-15, -.5, .5)), N, itemtype = '2PL')
dat <- rbind(dat1, dat2)
mod1 <- multipleGroup(dat, 1, group=group)
plot(mod1) #does not account for group differences! Need anchors


mod2 <- multipleGroup(dat, model, group=group, SE=TRUE,
                      invariance=c('free_means', 'free_var'))
plot(mod2)


# significant DIF in multiple items....
# DIF(mod2, which.par=c('a1', 'd'), items2test=16:30)
DTF(mod2)
#>    sDTF.score sDTF(%).score    uDTF.score uDTF(%).score 
#>     0.1937336     0.6457787     0.1943930     0.6479765 
DTF(mod2, draws=1000) #non-sig DTF due to item cancellation
#> $observed
#>    sDTF.score sDTF(%).score    uDTF.score uDTF(%).score 
#>     0.1937336     0.6457787     0.1943930     0.6479765 
#> 
#> $CIs
#>          sDTF.score sDTF(%).score uDTF.score uDTF(%).score
#> CI_97.5  0.49490356     1.6496785 0.66763070     2.2254357
#> CI_2.5  -0.06004841    -0.2001614 0.08144079     0.2714693
#> 
#> $tests
#> P(sDTF.score = 0) 
#>         0.1694727 
#> 

## -------------
## systematic differing slopes and intercepts (clear DTF)
dat1 <- simdata(a, d, N, itemtype = '2PL', mu=.50, sigma=matrix(1.5))
dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)), d + c(numeric(15), rnorm(n-15, 1, .5)),
                N, itemtype = '2PL')
dat <- rbind(dat1, dat2)
mod3 <- multipleGroup(dat, model, group=group, SE=TRUE,
                      invariance=c('free_means', 'free_var'))
plot(mod3) #visable DTF happening


# DIF(mod3, c('a1', 'd'), items2test=16:30)
DTF(mod3) #unsigned bias. Signed bias indicates group 2 scores generally higher on average
#>    sDTF.score sDTF(%).score    uDTF.score uDTF(%).score 
#>    -0.8316425    -2.7721417     1.0763924     3.5879746 
DTF(mod3, draws=1000)
#> $observed
#>    sDTF.score sDTF(%).score    uDTF.score uDTF(%).score 
#>    -0.8316425    -2.7721417     1.0763924     3.5879746 
#> 
#> $CIs
#>         sDTF.score sDTF(%).score uDTF.score uDTF(%).score
#> CI_97.5 -0.5502245     -1.834082  1.3842485      4.614162
#> CI_2.5  -1.1247665     -3.749222  0.7954355      2.651452
#> 
#> $tests
#> P(sDTF.score = 0) 
#>      2.746624e-08 
#> 
DTF(mod3, draws=1000, plot='func')

DTF(mod3, draws=1000, plot='sDTF') #multiple DTF areas along Theta


# evaluate specific values for sDTF
Theta_nodes <- matrix(seq(-6,6,length.out = 100))
sDTF <- DTF(mod3, Theta_nodes=Theta_nodes)
head(sDTF)
#>             Theta       sDTF
#> score.1 -6.000000 0.03153745
#> score.2 -5.878788 0.03604145
#> score.3 -5.757576 0.04118127
#> score.4 -5.636364 0.04704308
#> score.5 -5.515152 0.05372335
#> score.6 -5.393939 0.06132967
sDTF <- DTF(mod3, Theta_nodes=Theta_nodes, draws=100)
head(sDTF)
#>             Theta       sDTF   CI_97.5      CI_2.5
#> score.1 -6.000000 0.03153745 0.1630065 -0.01734414
#> score.2 -5.878788 0.03604145 0.1785908 -0.01768205
#> score.3 -5.757576 0.04118127 0.1955824 -0.01792158
#> score.4 -5.636364 0.04704308 0.2140928 -0.01804281
#> score.5 -5.515152 0.05372335 0.2342399 -0.01802425
#> score.6 -5.393939 0.06132967 0.2561474 -0.01784329

# }
```
