\name{logLik}
\alias{logLik}
\alias{logLik,polymirt-method}
\alias{logLik,confmirt-method}

\title{Monte Carlo Log-Likelihood Calculation}
\description{
Calculates a new object that contain the Monte Carlo estimated observed log-likelihood values for \code{polymirt} and \code{confmirt} objects
}

\usage{
logLik(object, ...)

\S4method{logLik}{polymirt}(object, draws = 3000, G2 = TRUE)

\S4method{logLik}{confmirt}(object, draws = 3000, G2 = TRUE)

}
 \arguments{
   \item{object}{a model of class \code{mirtClass} or \code{bfactorClass}}
   \item{draws}{the number of Monte Carlo draws}
   \item{G2}{logical; estimate the G2 model fit statistic?}
   \item{...}{parameters that are passed}
 }

\value{
 returns an object of class \code{polymirtClass} or \code{confmirtClass} with the log-likelihood, standard error, and (possibly) the G^2 model fit statistic.
}

\examples{
\dontrun{

mod1withLogLik <- logLik(mod1, draws = 5000)

  }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{logLik}
