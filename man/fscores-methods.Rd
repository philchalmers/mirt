\docType{methods}
\name{fscores}
\alias{fscores}
\alias{fscores,ConfirmatoryClass-method}
\alias{fscores,ExploratoryClass-method}
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

  \item{ndraws}{number of MH samples to draw for each
  response pattern}

  \item{thin}{controls how much the chain should be thinned
  by, default collects every 5th draw (\code{thin = 5}).
  Note that \code{ndraws/thin} must be a whole number. for
  \code{confmirtClass} objects only}

  \item{verbose}{logical; print verbose output messages?}
}
\value{
  Returns either a summary table with the response patterns
  and expected factor scores, or a complete data matrix
  with factor scores appended to the last column.
}
\description{
  Computes MAP, EAP, or ML factor scores for \code{mirt}
  and \code{bfactor} models, or a stochastic approximation
  with a multivariate normal prior for \code{polymirt} and
  \code{confmirt}. Note that only the general factor scores
  are computed for bifactor models.
}
\section{Methods}{
  \describe{ \item{fscores}{\code{signature(object =
  "ExploratoryClass")}}
  \item{fscores}{\code{signature(object =
  "ConfirmatoryClass")}}}
}
\examples{
\dontrun{

mod <- mirt(Science, 1)
tabscores <- fscores(mod)
fullscores <- fscores(mod, full.scores = TRUE)
fullscores <- fscores(mod, full.scores = TRUE, method='MAP')

mod2 <- confmirt(Science, 1)
tabscores2 <- fscores(mod2, ndraws = 5000)

  }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{factor.scores}
\keyword{methods}

