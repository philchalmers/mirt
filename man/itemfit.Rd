\name{itemfit}
\alias{itemfit}
\title{Item fit statistics}
\usage{
  itemfit(x, type = "Zh", ngroups = 10)
}
\arguments{
  \item{x}{a computed model object of class
  \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
  \code{MultipleGroupClass}}

  \item{type}{a character specifying whether the Zh
  (\code{'Zh'}) or \eqn{\chi^2} (\code{'X2'}) statistic
  should be computed. Not that \code{'X2'} can only be used
  for unidimensional models}

  \item{ngroups}{the number of theta groupings to use when
  computing \code{'X2'}. Cells that have any expected
  values less than 5 are dropped and the degrees of freedom
  are adjusted accordingly}
}
\description{
  \code{itemfit} calculates the Zh values from Drasgow,
  Levine and Williams (1985) for unidimensional and
  multidimensional models, or \eqn{chi^2} values for
  unidimensional models.
}
\examples{
\dontrun{

#make some data
set.seed(1234)
a <- matrix(rnorm(20),ncol=1)
d <- matrix(rnorm(20),ncol=1)
items <- rep('dich', 20)
data <- simdata(a,d, 2000, items)

x <- mirt(data, 1, SE = FALSE)
fit <- itemfit(x)
fit

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
  \code{\link{personfit}}
}
\keyword{fit}
\keyword{item}

