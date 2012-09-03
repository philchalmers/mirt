\docType{methods}
\name{fscores}
\alias{fscores}
\alias{fscores,ConfirmatoryClass-method}
\alias{fscores,ExploratoryClass-method}
\alias{fscores,MultipleGroupClass-method}
\alias{fscores-method}
\title{Methods for Function fscores}
\arguments{
  \item{object}{a model of class \code{ExploratoryClass} or
  \code{ConfirmatoryClass}}

  \item{full.scores}{if \code{FALSE} (default) then a
  summary table with factor scores for each unique pattern
  is displayed. Otherwise the original data matrix is
  returned with the computed factor scores appended to the
  rightmost column}

  \item{rotate}{rotation declaration to be used when
  estimating the factor scores. If \code{""} then the
  \code{object@rotate} default value is used}

  \item{method}{type of factor score estimation method. Can
  be expected a-posteriori (\code{"EAP"}), Bayes modal
  (\code{"MAP"}), or maximum likelihood (\code{"ML"})}

  \item{quadpts}{number of quadratures to use per
  dimension}

  \item{verbose}{logical; print verbose output messages?}
}
\description{
  Computes MAP, EAP, or ML factor scores with a
  multivariate normal prior distribution.
}
\section{Methods}{
  \describe{ \item{fscores}{\code{signature(object =
  "ExploratoryClass")}}
  \item{fscores}{\code{signature(object =
  "ConfirmatoryClass")}}
  \item{fscores}{\code{signature(object =
  "MultipleGroupClass")}}}
}
\examples{
\dontrun{

mod <- mirt(Science, 1)
tabscores <- fscores(mod)
fullscores <- fscores(mod, full.scores = TRUE)
fullscores <- fscores(mod, full.scores = TRUE, method='MAP')


  }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{factor.scores}
\keyword{methods}

