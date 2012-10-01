\name{personfit}
\alias{personfit}
\title{Person fit statistics}
\usage{
  personfit(x, degrees = NULL)
}
\arguments{
  \item{x}{a computed model object of class
  \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
  \code{MultipleGroupClass}}

  \item{degrees}{the degrees angle to be passed to the
  \code{\link{iteminfo}} function}
}
\description{
  \code{personfit} calculates the Zh values from Drasgow,
  Levine and Williams (1985) for unidimensional and
  multidimensional models, as well as infit and outfit
  statistics. The returned object is a \code{data.frame}
  consisting either of the tabulated data or full data with
  the statistics appended to the last columns.
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
fit <- personfit(x)
head(fit)

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
}
\seealso{
  \code{\link{itemfit}}
}
\keyword{fit}
\keyword{person}

