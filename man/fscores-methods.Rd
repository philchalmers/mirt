\docType{methods}
\name{fscores}
\alias{fscores}
\alias{fscores-method}
\alias{fscores,bfactorClass-method}
\alias{fscores,confmirtClass-method}
\alias{fscores,mirtClass-method}
\alias{fscores,polymirtClass-method}
\title{Methods for Function fscores}
\arguments{
  \item{object}{a model of class \code{mirtClass},
  \code{bfactorClass}, \code{polymirtClass}, or
  \code{confmirtClass}}

  \item{full.scores}{if \code{FALSE} (default) then a
  summary table with factor scores for each unique pattern
  is displayed. Otherwise the original data matrix is
  returned with the computed factor scores appended to the
  rightmost column}

  \item{method}{type of factor score estimation method. Can
  be expected a-posteriori (\code{"EAP"}), Bayes modal
  (\code{"MAP"}), or maximum likelihood (\code{"ML"})}

  \item{ndraws}{number of MH samples to draw for each
  response pattern}

  \item{thin}{controls how much the chain should be thinned
  by, default collects every 5th draw (\code{thin = 5}).
  Note that \code{ndraws/thin} must be a whole number}

  \item{...}{additional arguments to be passed}
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
  "bfactorClass")}} \item{fscores}{\code{signature(object =
  "mirtClass")}} \item{fscores}{\code{signature(object =
  "polymirtClass")}} \item{fscores}{\code{signature(object
  = "confmirtClass")}} }
}
\examples{
\dontrun{

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

