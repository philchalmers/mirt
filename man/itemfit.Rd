\name{itemfit}
\alias{itemfit}
\title{Item fit statistics}
\usage{
  itemfit(x, X2 = FALSE, group.size = 150,
    empirical.plot = NULL, degrees = NULL)
}
\arguments{
  \item{x}{a computed model object of class
  \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
  \code{MultipleGroupClass}}

  \item{X2}{logical; calculate the X2 statistic for
  unidimensional models?}

  \item{group.size}{approximate size of each group to be
  used in calculating the \eqn{\chi^2} statistic}

  \item{empirical.plot}{a single numeric value indicating
  which item to plot (via \code{itemplot}) and overlay with
  the empirical \eqn{\theta} groupings. Only applicable
  when \code{type = 'X2'}. The default is \code{NULL},
  therefore no plots are drawn}

  \item{degrees}{the degrees angle to be passed to the
  \code{\link{iteminfo}} function}
}
\description{
  \code{itemfit} calculates the Zh values from Drasgow,
  Levine and Williams (1985), infits, and outfits for
  unidimensional and multidimensional models, and
  \eqn{\chi^2} values for unidimensional models.
}
\examples{
\dontrun{

#make some data
set.seed(1234)
a <- matrix(rnorm(20),ncol=1)
d <- matrix(rnorm(20),ncol=1)
items <- rep('dich', 20)
data <- simdata(a,d, 2000, items)

x <- mirt(data, 1)
fit <- itemfit(x)
fit

itemfit(x, type = 'X2', empirical.plot = 1) #empirical item plot

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

  Reise, S. P. (1990). A comparison of item- and person-fit
  methods of assessing model-data fit in IRT. \emph{Applied
  Psychological Measurement, 14}, 127-137.
}
\seealso{
  \code{\link{personfit}}
}
\keyword{fit}
\keyword{item}

