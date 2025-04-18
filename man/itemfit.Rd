% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itemfit.R
\name{itemfit}
\alias{itemfit}
\title{Item fit statistics}
\usage{
itemfit(
  x,
  fit_stats = "S_X2",
  which.items = 1:extract.mirt(x, "nitems"),
  na.rm = FALSE,
  p.adjust = "none",
  group.bins = 10,
  group.size = NA,
  group.fun = mean,
  mincell = 1,
  mincell.X2 = 2,
  return.tables = FALSE,
  pv_draws = 30,
  boot = 1000,
  boot_dfapprox = 200,
  S_X2.plot = NULL,
  S_X2.plot_raw.score = TRUE,
  ETrange = c(-2, 2),
  ETpoints = 11,
  empirical.plot = NULL,
  empirical.CI = 0.95,
  empirical.poly.collapse = FALSE,
  method = "EAP",
  Theta = NULL,
  par.strip.text = list(cex = 0.7),
  par.settings = list(strip.background = list(col = "#9ECAE1"), strip.border = list(col =
    "black")),
  auto.key = list(space = "right", points = FALSE, lines = TRUE),
  ...
)
}
\arguments{
\item{x}{a computed model object of class \code{SingleGroupClass},
\code{MultipleGroupClass}, or \code{DiscreteClass}}

\item{fit_stats}{a character vector indicating which fit statistics should be computed.
  Supported inputs are:

\itemize{
  \item \code{'S_X2'} : Orlando and Thissen (2000, 2003) and
    Kang and Chen's (2007) signed chi-squared test (default)
  \item \code{'Zh'} : Drasgow, Levine, & Williams (1985) Zh
  \item \code{'X2'} : Bock's (1972) chi-squared method.
    The default inputs compute Yen's (1981) Q1 variant of the X2 statistic
    (i.e., uses a fixed \code{group.bins = 10}). However, Bock's group-size variable
    median-based method can be computed by passing \code{group.fun = median} and
    modifying the \code{group.size} input to the desired number of bins
  \item \code{'G2'} : McKinley & Mills (1985) G2 statistic (similar method to Q1,
    but with the likelihood-ratio test).
  \item \code{'PV_Q1'} : Chalmers and Ng's (2017) plausible-value variant
    of the Q1 statistic.
  \item \code{'PV_Q1*'} : Chalmers and Ng's (2017) plausible-value variant
    of the Q1 statistic that uses parametric bootstrapping to obtain a suitable empirical
    distribution.
  \item \code{'X2*'} : Stone's (2000) fit statistics that require parametric
    bootstrapping
  \item \code{'X2*_df'} : Stone's (2000) fit statistics that require parametric
    bootstrapping to obtain scaled versions of the X2* and degrees of freedom
  \item \code{'infit'} : Compute the infit and outfit statistics
}

Note that 'S_X2' and 'Zh' cannot be computed when there are missing response data
(i.e., will require multiple-imputation/row-removal techniques).}

\item{which.items}{an integer vector indicating which items to test for fit.
Default tests all possible items}

\item{na.rm}{logical; remove rows with any missing values? This is required for methods such
as S-X2 because they require the "EAPsum" method from \code{\link{fscores}}}

\item{p.adjust}{method to use for adjusting all p-values for each respective item fit
statistic (see \code{\link{p.adjust}} for available options). Default is \code{'none'}}

\item{group.bins}{the number of bins to use for X2 and G2. For example,
setting \code{group.bins = 10} will will compute Yen's (1981) Q1 statistic when \code{'X2'} is
requested}

\item{group.size}{approximate size of each group to be used in calculating the \eqn{\chi^2}
statistic. The default \code{NA}
disables this command and instead uses the \code{group.bins} input to try and construct
equally sized bins}

\item{group.fun}{function used when \code{'X2'} or \code{'G2'} are computed. Determines the central
tendency measure within each partitioned group. E.g., setting \code{group.fun = median} will
obtain the median of each respective ability estimate in each subgroup (this is what was used
by Bock, 1972)}

\item{mincell}{the minimum expected cell size to be used in the S-X2 computations. Tables will be
collapsed across items first if polytomous, and then across scores if necessary}

\item{mincell.X2}{the minimum expected cell size to be used in the X2 computations. Tables will be
collapsed if polytomous, however if this condition can not be met then the group block will
be omitted in the computations}

\item{return.tables}{logical; return tables when investigating \code{'X2'}, \code{'S_X2'},
and \code{'X2*'}?}

\item{pv_draws}{number of plausible-value draws to obtain for PV_Q1 and PV_Q1*}

\item{boot}{number of parametric bootstrap samples to create for PV_Q1* and X2*}

\item{boot_dfapprox}{number of parametric bootstrap samples to create for the X2*_df statistic
to approximate the scaling factor for X2* as well as the scaled degrees of freedom estimates}

\item{S_X2.plot}{argument input is the same as \code{empirical.plot}, however the resulting image
is constructed according to the S-X2 statistic's conditional sum-score information}

\item{S_X2.plot_raw.score}{logical; use the raw-score information in the plot in stead of the latent
trait scale score? Default is \code{FALSE}}

\item{ETrange}{range of integration nodes for Stone's X2* statistic}

\item{ETpoints}{number of integration nodes to use for Stone's X2* statistic}

\item{empirical.plot}{a single numeric value or character of the item name indicating which
item to plot (via \code{itemplot}) and overlay with the empirical \eqn{\theta} groupings (see
\code{empirical.CI}). Useful for plotting the expected bins based on the \code{'X2'} or
\code{'G2'} method}

\item{empirical.CI}{a numeric value indicating the width of the empirical confidence interval
ranging between 0 and 1 (default of 0 plots not interval). For example, a 95% confidence
interval would be plotted when \code{empirical.CI = .95}. Only applicable to dichotomous items}

\item{empirical.poly.collapse}{logical; collapse polytomous item categories to for expected scoring
functions for empirical plots? Default is \code{FALSE}}

\item{method}{type of factor score estimation method. See \code{\link{fscores}} for more detail}

\item{Theta}{a matrix of factor scores for each person used for statistics that require
empirical estimates. If supplied, arguments typically passed to \code{fscores()} will be
ignored and these values will be used instead. Also required when estimating statistics
with missing data via imputation}

\item{par.strip.text}{plotting argument passed to \code{\link[lattice]{lattice}}}

\item{par.settings}{plotting argument passed to \code{\link[lattice]{lattice}}}

\item{auto.key}{plotting argument passed to \code{\link[lattice]{lattice}}}

\item{...}{additional arguments to be passed to \code{fscores()} and \code{\link[lattice]{lattice}}}
}
\description{
Computes item-fit statistics for a variety of unidimensional and multidimensional models.
Poorly fitting items should be inspected with the empirical plots/tables
for unidimensional models, otherwise \code{\link{itemGAM}} can be used to diagnose
where the functional form of the IRT model was misspecified, or models can be refit using
more flexible semi-parametric response models (e.g., \code{itemtype = 'spline'}).
If the latent trait density was approximated (e.g., Davidian curves, Empirical histograms, etc)
then passing \code{use_dentype_estimate = TRUE} will use the internally saved quadrature and
density components (where applicable). Currently, only S-X2 statistic supported for
mixture IRT models. Finally, where applicable the root mean-square error of approximation (RMSEA)
is reported to help gauge the magnitude of item misfit.
}
\examples{

\donttest{

P <- function(Theta){exp(Theta^2 * 1.2 - 1) / (1 + exp(Theta^2 * 1.2 - 1))}

#make some data
set.seed(1234)
a <- matrix(rlnorm(20, meanlog=0, sdlog = .1),ncol=1)
d <- matrix(rnorm(20),ncol=1)
Theta <- matrix(rnorm(2000))
items <- rep('2PL', 20)
ps <- P(Theta)
baditem <- numeric(2000)
for(i in 1:2000)
   baditem[i] <- sample(c(0,1), 1, prob = c(1-ps[i], ps[i]))
data <- cbind(simdata(a,d, 2000, items, Theta=Theta), baditem=baditem)

x <- mirt(data, 1)
raschfit <- mirt(data, 1, itemtype='Rasch')
fit <- itemfit(x)
fit

# p-value adjustment
itemfit(x, p.adjust='fdr')

# two different fit stats (with/without p-value adjustment)
itemfit(x, c('S_X2' ,'X2'), p.adjust='fdr')
itemfit(x, c('S_X2' ,'X2'))

# Conditional sum-score plot from S-X2 information
itemfit(x, S_X2.plot = 1) # good fit
itemfit(x, S_X2.plot = 2) # good fit
itemfit(x, S_X2.plot = 21) # bad fit

itemfit(x, 'X2') # just X2
itemfit(x, 'X2', method = 'ML') # X2 with maximum-likelihood estimates for traits
itemfit(x, group.bins=15, empirical.plot = 1, method = 'ML') #empirical item plot with 15 points
itemfit(x, group.bins=15, empirical.plot = 21, method = 'ML')

# PV and X2* statistics (parametric bootstrap stats not run to save time)
itemfit(x, 'PV_Q1')

if(interactive()) mirtCluster() # improve speed of bootstrap samples by running in parallel
# itemfit(x, 'PV_Q1*')
# itemfit(x, 'X2*') # Stone's 1993 statistic
# itemfit(x, 'X2*_df') # Stone's 2000 scaled statistic with df estimate

# empirical tables for X2 statistic
tabs <- itemfit(x, 'X2', return.tables=TRUE, which.items = 1)
tabs

#infit/outfit statistics. method='ML' agrees better with eRm package
itemfit(raschfit, 'infit', method = 'ML') #infit and outfit stats

#same as above, but inputting ML estimates instead (saves time for re-use)
Theta <- fscores(raschfit, method = 'ML')
itemfit(raschfit, 'infit', Theta=Theta)
itemfit(raschfit, empirical.plot=1, Theta=Theta)
itemfit(raschfit, 'X2', return.tables=TRUE, Theta=Theta, which.items=1)

# fit a new more flexible model for the mis-fitting item
itemtype <- c(rep('2PL', 20), 'spline')
x2 <- mirt(data, 1, itemtype=itemtype)
itemfit(x2)
itemplot(x2, 21)
anova(x, x2)

#------------------------------------------------------------

#similar example to Kang and Chen 2007
a <- matrix(c(.8,.4,.7, .8, .4, .7, 1, 1, 1, 1))
d <- matrix(rep(c(2.0,0.0,-1,-1.5),10), ncol=4, byrow=TRUE)
dat <- simdata(a,d,2000, itemtype = rep('graded', 10))
head(dat)

mod <- mirt(dat, 1)
itemfit(mod)
itemfit(mod, 'X2') # less useful given inflated Type I error rates
itemfit(mod, empirical.plot = 1)
itemfit(mod, empirical.plot = 1, empirical.poly.collapse=TRUE)

# collapsed tables (see mincell.X2) for X2 and G2
itemfit(mod, 'X2', return.tables = TRUE, which.items = 1)

mod2 <- mirt(dat, 1, 'Rasch')
itemfit(mod2, 'infit', method = 'ML')

# massive list of tables for S-X2
tables <- itemfit(mod, return.tables = TRUE)

#observed and expected total score patterns for item 1 (post collapsing)
tables$O[[1]]
tables$E[[1]]

# can also select specific items
# itemfit(mod, return.tables = TRUE, which.items=1)

# fit stats with missing data (run in parallel using all cores)
dat[sample(1:prod(dim(dat)), 100)] <- NA
raschfit <- mirt(dat, 1, itemtype='Rasch')

# use only valid data by removing rows with missing terms
itemfit(raschfit, c('S_X2', 'infit'), na.rm = TRUE)

# note that X2, G2, PV-Q1, and X2* do not require complete datasets
thetas <- fscores(raschfit, method = 'ML') # save for faster computations
itemfit(raschfit, c('X2', 'G2'), Theta=thetas)
itemfit(raschfit, empirical.plot=1, Theta=thetas)
itemfit(raschfit, 'X2', return.tables=TRUE, which.items=1, Theta=thetas)

}

}
\references{
Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored
in two or more nominal categories. \emph{Psychometrika, 37}, 29-51.

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
\doi{10.18637/jss.v048.i06}

Chalmers, R. P. & Ng, V. (2017). Plausible-Value Imputation Statistics for Detecting
Item Misfit. \emph{Applied Psychological Measurement, 41}, 372-387.
\doi{10.1177/0146621617692079}

Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with
polychotomous item response models and standardized indices.
\emph{British Journal of Mathematical and Statistical Psychology, 38}, 67-86.

Kang, T. & Chen, Troy, T. (2007). An investigation of the performance of the generalized
S-X2 item-fit index for polytomous IRT models. ACT

McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit statistics.
Applied Psychological Measurement, 9, 49-57.

Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for dichotomous item
response theory models. \emph{Applied Psychological Measurement, 24}, 50-64.

Reise, S. P. (1990). A comparison of item- and person-fit methods of assessing model-data fit
in IRT. \emph{Applied Psychological Measurement, 14}, 127-137.

Stone, C. A. (2000). Monte Carlo Based Null Distribution for an Alternative Goodness-of-Fit
Test Statistics in IRT Models. \emph{Journal of Educational Measurement, 37}, 58-75.

Wright B. D. & Masters, G. N. (1982). \emph{Rating scale analysis}. MESA Press.

Yen, W. M. (1981). Using simulation results to choose a latent trait model.
\emph{Applied Psychological Measurement, 5}, 245-262.
}
\seealso{
\code{\link{personfit}}, \code{\link{itemGAM}}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{fit}
\keyword{item}
