\name{fscores}
\alias{fscores}
\alias{fscores,bfactor-method}
\alias{fscores,confmirt-method}
\alias{fscores,mirt-method}
\alias{fscores,polymirt-method}
\title{Compute factor scores}
\arguments{
  \item{object}{a model of class \code{mirtClass} or
  \code{bfactorClass}}

  \item{full.scores}{if \code{FALSE} (default) then a
  summary table with factor scores for each unique pattern
  is displayed. Otherwise the original data matrix is
  returned with the computed factor scores appended to the
  rightmost column}

  \item{method}{type of factor score estimation method. Can
  be expected a-posteriori (\code{"EAP"}) or Bayes modal
  (\code{"MAP"})}

  \item{ndraws}{number of MH samplers to draw for each
  response pattern}

  \item{thin}{controls how much the chain should be thinned
  by, default collects every 5th draw. Note that
  \code{ndraws/thin} must be a whole number}

  \item{...}{additional arguments to be passed}
}
\value{
  Returns either a summary table with the response patterns
  and expected factor scores, or a complete data matrix
  with factor scores appended to the last column.
}
\description{
  Computes MAP or EAP factor scores for \code{mirt} and
  \code{bfactor} models, or stocastic approximations for
  \code{polymirt} and \code{confmirt}. Note that only the
  general factor scores are computed for bifactor models.
}
\examples{
\dontrun{

tabscores <- fscores(mod)
fullscores <- fscores(mod, full.scores = TRUE)


  }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{factor.scores}

