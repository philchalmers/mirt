% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DRF.R
\name{DRF}
\alias{DRF}
\title{Differential Response Functioning statistics}
\usage{
DRF(
  mod,
  draws = NULL,
  focal_items = 1L:extract.mirt(mod, "nitems"),
  param_set = NULL,
  den.type = "marginal",
  best_fitting = FALSE,
  CI = 0.95,
  npts = 1000,
  quadpts = NULL,
  theta_lim = c(-6, 6),
  Theta_nodes = NULL,
  plot = FALSE,
  DIF = FALSE,
  DIF.cats = FALSE,
  groups2test = "all",
  pairwise = FALSE,
  simplify = TRUE,
  p.adjust = "none",
  par.strip.text = list(cex = 0.7),
  par.settings = list(strip.background = list(col = "#9ECAE1"), strip.border = list(col =
    "black")),
  auto.key = list(space = "right", points = FALSE, lines = TRUE),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{mod}{a multipleGroup object which estimated only 2 groups}

\item{draws}{a number indicating how many draws to take to form a suitable multiple imputation
or bootstrap estimate of the expected test scores (100 or more). If \code{boot = FALSE},
requires an estimated parameter information matrix. Returns a list containing the
bootstrap/imputation distribution and null hypothesis test for the sDRF statistics}

\item{focal_items}{a character/numeric vector indicating which items to include in the DRF tests. The
default uses all of the items (note that including anchors in the focal items has no effect
because they are exactly equal across groups). Selecting fewer items will result in tests of
'differential bundle functioning'}

\item{param_set}{an N x p matrix of parameter values drawn from the posterior (e.g., using the
parametric sampling approach, bootstrap, of MCMC). If supplied, then these will be used to compute
the DRF measures. Can be much more efficient to pre-compute these values if DIF, DBF, or DTF are
being evaluated within the same model (especially when using the bootstrap method).
See \code{\link{draw_parameters}}}

\item{den.type}{character specifying how the density of the latent traits is computed.
Default is \code{'marginal'} to include the proportional information from both groups,
\code{'focal'} for just the focal group, and \code{'reference'} for the reference group}

\item{best_fitting}{logical; use the best fitting parametric distribution (Gaussian by default)
that was used at the time of model estimation? This will result in much fast computations, however
the results are more dependent upon the underlying modelling assumptions. Default is FALSE, which
uses the empirical histogram approach}

\item{CI}{range of confidence interval when using draws input}

\item{npts}{number of points to use for plotting. Default is 1000}

\item{quadpts}{number of quadrature nodes to use when constructing DRF statistics. Default is extracted from
the input model object}

\item{theta_lim}{lower and upper limits of the latent trait (theta) to be evaluated, and is
used in conjunction with \code{quadpts} and \code{npts}}

\item{Theta_nodes}{an optional matrix of Theta values to be evaluated in the draws for the
sDRF statistics. However, these values are not averaged across, and instead give the bootstrap
confidence intervals at the respective Theta nodes. Useful when following up a large
sDRF or uDRF statistic, for example, to determine where the difference between the test curves are large
(while still accounting for sampling variability). Returns a matrix with observed
variability}

\item{plot}{logical; plot the 'sDRF' functions for the evaluated sDBF or sDTF values across the
integration grid or, if \code{DIF = TRUE}, the selected items as a faceted plot of individual items?
If plausible parameter sets were obtained/supplied then imputed confidence intervals will be included}

\item{DIF}{logical; return a list of item-level imputation properties using the DRF statistics?
These can generally be used as a DIF detection method and as a graphical display for
understanding DIF within each item}

\item{DIF.cats}{logical; same as \code{DIF = TRUE}, however computations will
be performed on each item category probability functions rather than the score functions.
Only useful for understanding DIF in polytomous items}

\item{groups2test}{when more than 2 groups are being investigated which two groups
should be used in the effect size comparisons?}

\item{pairwise}{logical; perform pairwise computations when the applying to multi-group settings}

\item{simplify}{logical; attempt to simplify the output rather than returning larger lists?}

\item{p.adjust}{string to be passed to the \code{\link{p.adjust}} function to adjust p-values.
Adjustments are located in the \code{adj_pvals} element in the returned list. Only applicable when
\code{DIF = TRUE}}

\item{par.strip.text}{plotting argument passed to \code{\link[lattice]{lattice}}}

\item{par.settings}{plotting argument passed to \code{\link[lattice]{lattice}}}

\item{auto.key}{plotting argument passed to \code{\link[lattice]{lattice}}}

\item{verbose}{logical; include additional information in the console?}

\item{...}{additional arguments to be passed to \code{\link[lattice]{lattice}}}
}
\description{
Function performs various omnibus differential item (DIF), bundle (DBF), and test (DTF)
functioning procedures on an object
estimated with \code{multipleGroup()}. The compensatory and non-compensatory statistics provided
are described in Chalmers (2018), which generally can be interpreted as IRT generalizations
of the SIBTEST and CSIBTEST statistics. For hypothesis tests, these measures
require the ACOV matrix to be computed in the
fitted multiple-group model (otherwise, sets of plausible draws from the posterior are explicitly
required).
}
\details{
The effect sizes estimates by the DRF function are
\deqn{sDRF = \int [S(C|\bm{\Psi}^{(R)},\theta) S(C|\bm{\Psi}^{(F)},\theta)] f(\theta)d\theta,}
\deqn{uDRF = \int |S(C|\bm{\Psi}^{(R)},\theta) S(C|\bm{\Psi}^{(F)},\theta)| f(\theta)d\theta,}
and
\deqn{dDRF = \sqrt{\int [S(C|\bm{\Psi}^{(R)},\theta) S(C|\bm{\Psi}^{(F)},\theta)]^2 f(\theta)d\theta}}
where \eqn{S(.)} are the scoring equations used to evaluate the model-implied
difference between the focal and reference group. The \eqn{f(\theta)}
terms can either be estimated from the posterior via an empirical
histogram approach (default), or can use the best
fitting prior distribution that is obtain post-convergence (default is a Guassian
distribution). Note that, in comparison to Chalmers (2018), the focal group is
the leftmost scoring function while the reference group is the rightmost
scoring function. This is largely to keep consistent with similar effect
size statistics, such as SIBTEST, DFIT, Wainer's measures of impact, etc,
which in general can be seen as special-case estimators of this family.
}
\examples{
\donttest{

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
dat1 <- simdata(a, d, N, itemtype = 'dich')
dat2 <- simdata(a, d, N, itemtype = 'dich')
dat <- rbind(dat1, dat2)
mod <- multipleGroup(dat, model, group=group, SE=TRUE,
                     invariance=c('free_means', 'free_var'))
plot(mod)
plot(mod, which.items = 6:10) #DBF
plot(mod, type = 'itemscore')
plot(mod, type = 'itemscore', which.items = 10:15)

# empirical histogram approach
DRF(mod)
DRF(mod, focal_items = 6:10) #DBF
DRF(mod, DIF=TRUE)
DRF(mod, DIF=TRUE, focal_items = 10:15)

# Best-fitting Gaussian distributions
DRF(mod, best_fitting=TRUE)
DRF(mod, focal_items = 6:10, best_fitting=TRUE) #DBF
DRF(mod, DIF=TRUE, best_fitting=TRUE)
DRF(mod, DIF=TRUE, focal_items = 10:15, best_fitting=TRUE)

DRF(mod, plot = TRUE)
DRF(mod, focal_items = 6:10, plot = TRUE) #DBF
DRF(mod, DIF=TRUE, plot = TRUE)
DRF(mod, DIF=TRUE, focal_items = 10:15, plot = TRUE)

if(interactive()) mirtCluster()
DRF(mod, draws = 500)
DRF(mod, draws = 500, best_fitting=TRUE)
DRF(mod, draws = 500, plot=TRUE)

# pre-draw parameter set to save computations
#  (more useful when using non-parametric bootstrap)
param_set <- draw_parameters(mod, draws = 500)
DRF(mod, focal_items = 6, param_set=param_set) #DIF test
DRF(mod, DIF=TRUE, param_set=param_set) #DIF test
DRF(mod, focal_items = 6:10, param_set=param_set) #DBF test
DRF(mod, param_set=param_set) #DTF test

DRF(mod, focal_items = 6:10, draws=500) #DBF test
DRF(mod, focal_items = 10:15, draws=500) #DBF test

DIFs <- DRF(mod, draws = 500, DIF=TRUE)
print(DIFs)
DRF(mod, draws = 500, DIF=TRUE, plot=TRUE)

DIFs <- DRF(mod, draws = 500, DIF=TRUE, focal_items = 6:10)
print(DIFs)
DRF(mod, draws = 500, DIF=TRUE, focal_items = 6:10, plot = TRUE)

DRF(mod, DIF=TRUE, focal_items = 6)
DRF(mod, draws=500, DIF=TRUE, focal_items = 6)

# evaluate specific values for sDRF
Theta_nodes <- matrix(seq(-6,6,length.out = 100))

sDTF <- DRF(mod, Theta_nodes=Theta_nodes)
head(sDTF)
sDTF <- DRF(mod, Theta_nodes=Theta_nodes, draws=200)
head(sDTF)

# sDIF (isolate single item)
sDIF <- DRF(mod, Theta_nodes=Theta_nodes, focal_items=6)
head(sDIF)
sDIF <- DRF(mod, Theta_nodes=Theta_nodes, focal_items = 6, draws=200)
head(sDIF)

## -------------
## random slopes and intercepts for 15 items, and latent mean difference
##    (no systematic DTF should exist, but DIF will be present)
set.seed(1234)
dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 0, .25)),
                d + c(numeric(15), rnorm(n-15, 0, .5)), N, itemtype = 'dich')
dat <- rbind(dat1, dat2)
mod1 <- multipleGroup(dat, 1, group=group)
plot(mod1)
DRF(mod1) #does not account for group differences! Need anchors

mod2 <- multipleGroup(dat, model, group=group, SE=TRUE,
                      invariance=c('free_means', 'free_var'))
plot(mod2)

# significant DIF in multiple items....
# DIF(mod2, which.par=c('a1', 'd'), items2test=16:30)
DRF(mod2)
DRF(mod2, draws=500) #non-sig DTF due to item cancellation

## -------------
## systematic differing slopes and intercepts (clear DTF)
set.seed(1234)
dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)),
                d + c(numeric(15), rnorm(n-15, 1, .5)),
                N, itemtype = 'dich')
dat <- rbind(dat1, dat2)
mod3 <- multipleGroup(dat, model, group=group, SE=TRUE,
                      invariance=c('free_means', 'free_var'))
plot(mod3) #visable DTF happening

# DIF(mod3, c('a1', 'd'), items2test=16:30)
DRF(mod3) #unsigned bias. Signed bias (group 2 scores higher on average)
DRF(mod3, draws=500)
DRF(mod3, draws=500, plot=TRUE) #multiple DRF areas along Theta

# plot the DIF
DRF(mod3, draws=500, DIF=TRUE, plot=TRUE)

# evaluate specific values for sDRF
Theta_nodes <- matrix(seq(-6,6,length.out = 100))
sDTF <- DRF(mod3, Theta_nodes=Theta_nodes, draws=200)
head(sDTF)

# DIF
sDIF <- DRF(mod3, Theta_nodes=Theta_nodes, focal_items = 30, draws=200)
head(sDIF)

## ----------------------------------------------------------------
# polytomous example
# simulate data where group 2 has a different slopes/intercepts
set.seed(4321)
a1 <- a2 <- matrix(rlnorm(20,.2,.3))
a2[c(16:17, 19:20),] <- a1[c(16:17, 19:20),] + c(-.5, -.25, .25, .5)

# for the graded model, ensure that there is enough space between the intercepts,
# otherwise closer categories will not be selected often
diffs <- t(apply(matrix(runif(20*4, .3, 1), 20), 1, cumsum))
diffs <- -(diffs - rowMeans(diffs))
d1 <- d2 <- diffs + rnorm(20)
rownames(d1) <- rownames(d2) <- paste0('Item.', 1:20)
d2[16:20,] <- d1[16:20,] + matrix(c(-.5, -.5, -.5, -.5,
                                    1, 0, 0, -1,
                                    .5, .5, -.5, -.5,
                                    1, .5, 0, -1,
                                    .5, .5, .5, .5), byrow=TRUE, nrow=5)

tail(data.frame(a.group1 = a1, a.group2 = a2), 6)
list(d.group1 = d1[15:20,], d.group2 = d2[15:20,])

itemtype <- rep('graded', nrow(a1))
N <- 600
dataset1 <- simdata(a1, d1, N, itemtype)
dataset2 <- simdata(a2, d2, N, itemtype, mu = -.25, sigma = matrix(1.25))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))

# item 1-10 as anchors
mod <- multipleGroup(dat, group=group, SE=TRUE,
                     invariance=c(colnames(dat)[1:10], 'free_means', 'free_var'))
coef(mod, simplify=TRUE)
plot(mod)
plot(mod, type='itemscore')

# DIF tests vis Wald method
DIF(mod, items2test=11:20,
   which.par=c('a1', paste0('d', 1:4)),
   Wald=TRUE, p.adjust='holm')

DRF(mod)
DRF(mod, DIF=TRUE, focal_items=11:20)
DRF(mod, DIF.cats=TRUE, focal_items=11:20)

## ----------------------------------------------------------------
### multidimensional DTF

set.seed(1234)
n <- 50
N <- 1000

# only first 5 items as anchors within each dimension
model <- 'F1 = 1-25
          F2 = 26-50
          COV = F1*F2
          CONSTRAINB = (1-5, a1), (1-5, 26-30, d), (26-30, a2)'

a <- matrix(c(rep(1, 25), numeric(50), rep(1, 25)), n)
d <- matrix(rnorm(n), n)
group <- c(rep('Group_1', N), rep('Group_2', N))
Cov <- matrix(c(1, .5, .5, 1.5), 2)
Mean <- c(0, 0.5)

# groups completely equal
dat1 <- simdata(a, d, N, itemtype = 'dich', sigma = cov2cor(Cov))
dat2 <- simdata(a, d, N, itemtype = 'dich', sigma = Cov, mu = Mean)
dat <- rbind(dat1, dat2)
mod <- multipleGroup(dat, model, group=group, SE=TRUE,
                     invariance=c('free_means', 'free_var'))
coef(mod, simplify=TRUE)
plot(mod, degrees = c(45,45))
DRF(mod)

# some intercepts slightly higher in Group 2
d2 <- d
d2[c(10:15, 31:35)] <- d2[c(10:15, 31:35)] + 1
dat1 <- simdata(a, d, N, itemtype = 'dich', sigma = cov2cor(Cov))
dat2 <- simdata(a, d2, N, itemtype = 'dich', sigma = Cov, mu = Mean)
dat <- rbind(dat1, dat2)
mod <- multipleGroup(dat, model, group=group, SE=TRUE,
                     invariance=c('free_means', 'free_var'))
coef(mod, simplify=TRUE)
plot(mod, degrees = c(45,45))

DRF(mod)
DRF(mod, draws = 500)

}
}
\references{
Chalmers, R. P. (2018). Model-Based Measures for Detecting and Quantifying Response Bias.
  \emph{Psychometrika, 83}(3), 696-732. \doi{10.1007/s11336-018-9626-9}
}
\seealso{
\code{\link{multipleGroup}}, \code{\link{DIF}}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{differential}
\keyword{functioning}
\keyword{response}
