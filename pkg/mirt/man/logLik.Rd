\name{logLik}
\alias{logLik}
\alias{logLik,polymirt-method}
\alias{logLik,confmirt-method}

\title{Monte Carlo Log Likelihood Calculation}
\description{
Calculate new objects that contain the Monte Carlo estimated observed log-likelihood values for \code{polymirt} and \code{confmirt} objects
}

\usage{
logLik(object, ...)

\S4method{logLik}{polymirt}(object, draws = 2000)

\S4method{logLik}{confmirt}(object, draws = 2000)

}
 \arguments{
   \item{object}{a model of class \code{mirt} or \code{bfactor}}
   \item{draws}{the number of Monte Carlo draws}
 }

\value{
 returns an object of \code{polymirt} or \code{confmirt} with the log-likelihood and standard error. 
}

\examples{
\dontrun{

mod2 <- logLik(mod1, draws = 5000)

  }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{logLik}
