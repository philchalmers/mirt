\name{itemfit}
\alias{itemfit}
\title{Item fit statistics}
\usage{
  itemfit(x, Zh = TRUE, X2 = FALSE, group.size = 150,
    mincell = 1, S_X2.tables = FALSE,
    empirical.plot = NULL, method = "EAP", ...)
}
\arguments{
  \item{x}{a computed model object of class
  \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
  \code{MultipleGroupClass}}

  \item{Zh}{logical; calculate Zh and associated statistics
  (infit/outfit)? Disable this is you are only interested
  in computing the S-X2 quickly}

  \item{X2}{logical; calculate the X2 statistic for
  unidimensional models?}

  \item{mincell}{the minimum expected cell size to be used
  in the S-X2 computations. Tables will be collapsed across
  items first if polytomous, and then across scores if
  necessary}

  \item{S_X2.tables}{logical; return the tables in a list
  format used to compute the S-X2 stats?}

  \item{group.size}{approximate size of each group to be
  used in calculating the \eqn{\chi^2} statistic}

  \item{empirical.plot}{a single numeric value or character
  of the item name indicating which item to plot (via
  \code{itemplot}) and overlay with the empirical
  \eqn{\theta} groupings. Only applicable when \code{type =
  'X2'}. The default is \code{NULL}, therefore no plots are
  drawn}

  \item{method}{type of factor score estimation method. See
  \code{\link{fscores}} for more detail}

  \item{...}{additional arguments to be passed to
  \code{fscores()}}
}
\description{
  \code{itemfit} calculates the Zh values from Drasgow,
  Levine and Williams (1985), \eqn{\chi^2} values for
  unidimensional models, and S-X2 statistics for
  unidimensional models (Kang & Chen, 2007; Orlando &
  Thissen, 2000). For Rasch, partial credit, and rating
  scale models infit and outfit statistics are also
  produced.
}
\examples{
\dontrun{

#make some data
set.seed(1234)
a <- matrix(rlnorm(20, meanlog=0, sdlog = .1),ncol=1)
d <- matrix(rnorm(20),ncol=1)
items <- rep('dich', 20)
data <- simdata(a,d, 2000, items)

x <- mirt(data, 1)
raschfit <- mirt(data, 1, itemtype='Rasch')
fit <- itemfit(x)
fit

itemfit(x, empirical.plot = 1) #empirical item plot
itemfit(raschfit, method = 'ML') #infit and outfit stats (method='ML' agrees better with eRm package)

#similar example to Kang and Chen 2007
a <- matrix(c(.8,.4,.7, .8, .4, .7, 1, 1, 1, 1))
d <- matrix(rep(c(2.0,0.0,-1,-1.5),10), ncol=4, byrow=TRUE)
dat <- simdata(a,d,2000, itemtype = rep('graded', 10)) - 1
head(dat)

mod <- mirt(dat, 1)
itemfit(mod)

mod2 <- mirt(dat, 1, 'Rasch')
itemfit(mod2)

#massive list of tables
tables <- itemfit(mod, S_X2.tables = TRUE)

#observed and expected total score patterns for item 1 (post collapsing)
tables$O[[1]]
tables$E[[1]]


  }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Drasgow, F., Levine, M. V., & Williams, E. A. (1985).
  Appropriateness measurement with polychotomous item
  response models and standardized indices. \emph{Journal
  of Mathematical and Statistical Psychology, 38}, 67-86.

  Kang, T. & Chen, Troy, T. (2007). An investigation of the
  performance of the generalized S-X2 item-fit index for
  polytomous IRT models. ACT

  Orlando, M. & Thissen, D. (2000). Likelihood-based item
  fit indices for dichotomous item response theory models.
  \emph{Applied Psychological Measurement, 24}, 50-64.

  Reise, S. P. (1990). A comparison of item- and person-fit
  methods of assessing model-data fit in IRT. \emph{Applied
  Psychological Measurement, 14}, 127-137.
}
\seealso{
  \code{\link{personfit}}
}
\keyword{fit}
\keyword{item}

