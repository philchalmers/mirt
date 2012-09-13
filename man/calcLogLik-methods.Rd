\docType{methods}
\name{calcLogLik}
\alias{calcLogLik}
\alias{calcLogLik,ConfirmatoryClass-method}
\alias{calcLogLik,ExploratoryClass-method}
\alias{calcLogLik-method}
\title{Monte Carlo Log-Likelihood Calculation}
\arguments{
  \item{object}{a model of class \code{ConfirmatoryClass}
  or \code{ExploratoryClass}}

  \item{draws}{the number of Monte Carlo draws}

  \item{G2}{logical; estimate the G2 model fit statistic?}

  \item{...}{parameters that are passed}
}
\value{
  Returns an object with the log-likelihood, standard
  errors, information matrix, and (possibly) the G^2 and
  other model fit statistic if there is no missing data.
}
\description{
  Calculates a new object that contain the Monte Carlo
  estimated observed log-likelihood values for mirt objects
  estimated with the MH-RM algorithm
}
\section{Methods}{
  \describe{ \item{calcLogLik}{\code{signature(object =
  "ConfirmatoryClass")}, \code{signature(object =
  "ExploratoryClass")} } }
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
  \code{\link{confmirt}}, \code{\link{multipleGroup}}
}
\keyword{calcLogLik}

