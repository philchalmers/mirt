\docType{methods}
\name{calcLogLik}
\alias{calcLogLik}
\alias{calcLogLik,confmirtClass-method}
\alias{calcLogLik-method}
\title{Monte Carlo Log-Likelihood Calculation}
\arguments{
  \item{object}{a model of class \code{confmirtClass}}

  \item{draws}{the number of Monte Carlo draws}

  \item{G2}{logical; estimate the G2 model fit statistic?}

  \item{...}{parameters that are passed}
}
\value{
  Returns an object of class \code{confmirtClass} with the
  log-likelihood, standard error, and (possibly) the G^2
  model fit statistic if there is no missing data.
}
\description{
  Calculates a new object that contain the Monte Carlo
  estimated observed log-likelihood values for
  \code{polymirt} and \code{confmirt} objects
}
\section{Methods}{
  \describe{ \item{calcLogLik}{\code{signature(object =
  "confmirtClass")}} }
}
\examples{
\dontrun{

mod1withLogLik <- calcLogLik(mod1, draws = 5000)

  }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{confmirt}}
}
\keyword{calcLogLik}

