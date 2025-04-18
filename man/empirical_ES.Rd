% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/empirical_ES.R
\name{empirical_ES}
\alias{empirical_ES}
\title{Empirical effect sizes based on latent trait estimates}
\usage{
empirical_ES(
  mod,
  Theta.focal = NULL,
  focal_items = 1L:extract.mirt(mod, "nitems"),
  DIF = TRUE,
  npts = 61,
  theta_lim = c(-6, 6),
  plot = FALSE,
  type = "b",
  par.strip.text = list(cex = 0.7),
  par.settings = list(strip.background = list(col = "#9ECAE1"), strip.border = list(col =
    "black")),
  ...
)
}
\arguments{
\item{mod}{a multipleGroup object which estimated only 2 groups. The first group in this object
is assumed to be the reference group by default (i.e., \code{ref.group = 1}), which conforms to the
 \code{invariance} arguments in \code{\link{multipleGroup}}}

\item{Theta.focal}{an optional matrix of Theta values from the focal group to be evaluated. If not supplied
the default values to \code{\link{fscores}} will be used in conjunction with the \code{...}
arguments passed}

\item{focal_items}{a numeric vector indicating which items to include the tests. The
default uses all of the items. Selecting fewer items will result in tests of
'differential bundle functioning' when \code{DIF = FALSE}}

\item{DIF}{logical; return a data.frame of item-level imputation properties? If \code{FALSE},
only DBF and DTF statistics will be reported}

\item{npts}{number of points to use in the integration. Default is 61}

\item{theta_lim}{lower and upper limits of the latent trait (theta) to be evaluated, and is
used in conjunction with \code{npts}}

\item{plot}{logical; plot expected scores of items/test where expected scores are computed
using focal group thetas and both focal and reference group item parameters}

\item{type}{type of objects to draw in \code{lattice}; default plots both points and lines}

\item{par.strip.text}{plotting argument passed to \code{\link[lattice]{lattice}}}

\item{par.settings}{plotting argument passed to \code{\link[lattice]{lattice}}}

\item{...}{additional arguments to be passed to \code{\link{fscores}} and \code{\link[lattice]{xyplot}}}
}
\description{
Computes effect size measures of differential item functioning and differential
test/bundle functioning based on expected scores from Meade (2010).
Item parameters from both reference and focal group are used in conjunction with
focal group empirical theta estimates (and an assumed normally distributed theta)
to compute expected scores.
}
\section{DIF}{


The default \code{DIF = TRUE} produces several effect sizes indices at the item level.
Signed indices allow DIF favoring the focal group at one point on the theta
distribution to cancel DIF favoring the reference group at another point on the theta
distribution. Unsigned indices take the absolute value before summing or averaging,
thus not allowing cancellation of DIF across theta.

\describe{
  \item{SIDS}{Signed Item Difference in the Sample. The average difference in expected scores
across the focal sample using both focal and reference group item parameters.}
  \item{UIDS}{Unsigned Item Difference in the Sample. Same as SIDS except absolute value of
expected scores is taken prior to averaging across the sample.}
  \item{D-Max}{The maximum difference in expected scores in the sample.}
  \item{ESSD}{Expected Score Standardized Difference. Cohen's D for difference in expected scores.}
  \item{SIDN}{Signed Item Difference in a Normal distribution. Identical to SIDS but
averaged across a normal distribution rather than the sample.}
  \item{UIDN}{Unsigned Item Difference in a Normal distribution. Identical to UIDS but
averaged across a normal distribution rather than the sample.}
}
}

\section{DBF/DTF}{


\code{DIF = FALSE} produces a series of test/bundle-level indices that are based on item-level
indices.

\describe{
  \item{STDS}{Signed Test Differences in the Sample. The sum of the SIDS across items.}
  \item{UTDS}{Unsigned Test Differences in the Sample. The sum of the UIDS across items.}
  \item{Stark's DTFR}{Stark's version of STDS using a normal distribution rather than
sample estimated thetas.}
  \item{UDTFR}{Unsigned Expected Test Scores Differences in the Sample. The difference
in observed summed scale scores expected, on average, across a hypothetical focal
group with a normally distributed theta, had DF been uniform in nature for all items}
  \item{UETSDS}{Unsigned Expected Test Score Differences in the Sample.
The hypothetical difference expected scale scores that would have been present if
scale-level DF had been uniform across respondents (i.e., always favoring the
focal group).}
  \item{UETSDN}{Identical to UETSDS but computed using a normal distribution.}
  \item{Test D-Max}{Maximum expected test score differences in the sample.}
  \item{ETSSD}{Expected Test Score Standardized Difference. Cohen's D for expected
test scores.}
}
}

\examples{
\donttest{

# no DIF
set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
itemtype <- rep('2PL', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)

# ensure 'Ref' is the first group (and therefore reference group during estimation)
group <- factor(c(rep('Ref', N), rep('Focal', N)), levels = c('Ref', 'Focal'))

mod <- multipleGroup(dat, 1, group = group,
   invariance = c(colnames(dat)[1:5], 'free_means', 'free_var'))
coef(mod, simplify=TRUE)

empirical_ES(mod)
empirical_ES(mod, DIF=FALSE)
empirical_ES(mod, DIF=FALSE, focal_items = 10:15)

empirical_ES(mod, plot=TRUE)
empirical_ES(mod, plot=TRUE, DIF=FALSE)

###---------------------------------------------
# DIF
set.seed(12345)
a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
a2[10:15,] <- a2[10:15,] + rnorm(6, 0, .3)
d2[10:15,] <- d2[10:15,] + rnorm(6, 0, .3)
itemtype <- rep('dich', nrow(a1))
N <- 1000
dataset1 <- simdata(a1, d1, N, itemtype)
dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- factor(c(rep('Ref', N), rep('Focal', N)), levels = c('Ref', 'Focal'))

mod <- multipleGroup(dat, 1, group = group,
   invariance = c(colnames(dat)[1:5], 'free_means', 'free_var'))
coef(mod, simplify=TRUE)

empirical_ES(mod)
empirical_ES(mod, DIF = FALSE)
empirical_ES(mod, plot=TRUE)
empirical_ES(mod, plot=TRUE, DIF=FALSE)

}
}
\references{
Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
\doi{10.18637/jss.v048.i06}

Meade, A. W. (2010). A taxonomy of effect size measures for the differential functioning
of items and scales. \emph{Journal of Applied Psychology, 95}, 728-743.
}
\author{
Adam Meade, with contributions by Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
