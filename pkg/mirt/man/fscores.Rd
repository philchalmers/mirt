\name{fscores}
\alias{fscores}
\alias{fscores,mirt-method}
\alias{fscores,bfactor-method}

\title{Compute factor scores}
\description{
Computes MAP or EAP factor scores for \code{mirt} and \code{bfactor} models. Note that only the general factor scores are computed for bifactor models.
}

\usage{
fscores(object, ...)

\S4method{fscores}{mirt}(object, full.scores = FALSE,
  method = "EAP", ...)

\S4method{fscores}{bfactor}(object, full.scores = FALSE,
  method = "EAP", ...)

}
 \arguments{
   \item{object}{a model of class \code{mirtClass} or \code{bfactorClass}}
   \item{full.scores}{if \code{FALSE} (default) then a summary table with factor scores for each unique pattern is displayed. Otherwise the original data matrix is returned with the computed factor scores}
  \item{method}{type of factor score estimation method. Can be expected a-posteriori (\code{"EAP"}) or Bayes modal (\code{"MAP"})}
   \item{...}{additional arguments to be passed}
 }

\value{
 returns either a summary table with the response patterns and expected factor scores, or a complete data matrix with factor scores appended to the last column.
}

\examples{
\dontrun{

tabscores <- fscores(mod)
fullscores <- fscores(mod, full.scores = TRUE)


  }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{factor.scores}
