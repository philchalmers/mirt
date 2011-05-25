\name{f.scores}
\alias{f.scores}
\alias{f.scores.mirt}
\alias{f.scores.bfactor}

\title{Compute factor scores}
\description{
Computes MAP or EAP factor scores for \code{mirt} and \code{bfactor} models. Note that only the general factor scores are computed for bifactor models.
}

\usage{
f.scores(object, ...)

\method{f.scores}{mirt}(object, full.scores = FALSE,
  method = "EAP", rotate = 'varimax', ...)
\method{f.scores}{bfactor}(object, full.scores = FALSE,
  method = "EAP", ...)

}
 \arguments{
   \item{object}{a model of class \code{mirt} or \code{bfactor}}
   \item{full.scores}{if \code{FALSE} (default) then a summary table with factor scores for each unique pattern is displayed. Otherwise the original data matrix is returned with the computed factor scores}
  \item{method}{type of factor score estimation method. Can be expected a-posteriori (\code{"EAP"}) or Bayes modal (\code{"MAP"})}
   \item{rotate}{type of rotation to use for before computing factor scores. See \code{\link{mirt}} for a list of possible rotations}
   \item{...}{additional arguments to be passed}
 }

\value{
 returns either a summary table with the response patterns and expected factor scores, or a complete data matrix with factor scores appended to the last column.
}

\examples{
\dontrun{

scores <- f.scores(mod)
scores

  }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{factor.scores}
