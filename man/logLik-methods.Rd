\docType{methods}
\name{logLik}
\alias{logLik}
\alias{logLik,confmirtClass-method}
\alias{logLik-method}
\alias{logLik,polymirtClass-method}
\title{Monte Carlo Log-Likelihood Calculation}
\arguments{
  \item{object}{a model of class \code{mirtClass} or
  \code{bfactorClass}}

  \item{draws}{the number of Monte Carlo draws}

  \item{G2}{logical; estimate the G2 model fit statistic?}

  \item{...}{parameters that are passed}
}
\value{
  Returns an object of class \code{polymirtClass} or
  \code{confmirtClass} with the log-likelihood, standard
  error, and (possibly) the G^2 model fit statistic.
}
\description{
  Calculates a new object that contain the Monte Carlo
  estimated observed log-likelihood values for
  \code{polymirt} and \code{confmirt} objects
}
\section{Methods}{
  \describe{ \item{logLik}{\code{signature(object =
  "polymirtClass")}} \item{logLik}{\code{signature(object =
  "confmirtClass")}} }
}
\examples{
\dontrun{

mod1withLogLik <- logLik(mod1, draws = 5000)

  }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{polymirt}}, \code{\link{confmirt}}
}
\keyword{logLik}

