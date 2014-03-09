\docType{methods}
\name{calcLogLik}
\alias{calcLogLik}
\alias{calcLogLik,ConfirmatoryClass-method}
\alias{calcLogLik,ExploratoryClass-method}
\alias{calcLogLik,MixedClass-method}
\alias{calcLogLik-method}
\title{Monte Carlo Log-Likelihood Calculation}
\usage{
calcLogLik(object, ...)

\S4method{calcLogLik}{ExploratoryClass}(object,
   draws = 5000, G2 = TRUE)
\S4method{calcLogLik}{ConfirmatoryClass}(object,
   draws = 5000, G2 = TRUE)
\S4method{calcLogLik}{MixedClass}(object,
   draws = 5000)
}
\arguments{
  \item{object}{a model of class \code{ConfirmatoryClass}
  or \code{ExploratoryClass}}

  \item{draws}{the number of Monte Carlo draws}

  \item{G2}{logical; estimate the G2 model fit statistic?}

  \item{...}{parameters that are passed}
}
\value{
Returns an object with the log-likelihood and Monte Carlo
standard errors, and (possibly) the G^2 and other model fit
statistic if there is no missing data.
}
\description{
Calculates a new object that contain the Monte Carlo
estimated observed log-likelihood values for mirt objects
estimated with the MH-RM algorithm. Can be estimated in
parallel automatically by defining a parallel object with
\code{\link{mirtCluster}}.
}
\section{Methods}{
  \describe{ \item{calcLogLik}{\code{signature(object =
  "ConfirmatoryClass")}, \code{signature(object =
  "ExploratoryClass")}, \code{signature(object =
  "MixedClass")} } }
}
\examples{
\dontrun{

# no parallel
mod1withLogLik <- calcLogLik(mod1, draws=5000)

#with parallel using detected number of cores
library(parallel)
mirtCluster(detectCores())
mod1withLogLik <- calcLogLik(mod1, draws=5000)

  }
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
\code{\link{mirt}}, \code{\link{multipleGroup}},
\code{\link{mixedmirt}}
}
\keyword{calcLogLik}

