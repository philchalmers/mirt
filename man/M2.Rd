\name{M2}
\alias{M2}
\title{Compute M2 statistic}
\usage{
M2(obj, calcNull = FALSE, quadpts = NULL)
}
\arguments{
  \item{obj}{an estimated model object from the mirt
  package}

  \item{quadpts}{number of quadrature points to use during
  estimation. If \code{NULL}, \code{quadpts} will be
  extracted from the \code{obj}; if not available, a
  suitable value will be chosen based on the rubric found
  in \code{\link{mirt}}}

  \item{calcNull}{logical; calculate statistics for the
  null model as well?  Allows for statistics such as the
  limited information TLI and CFI}
}
\description{
Computes the M2 (Maydeu-Olivares & Joe, 2006) statistic for
dichotomous data and the M2* statistic for polytomous data
(collapsing over response categories for better stability;
see Cai and Hansen, 2013), as well as associated fit
indicies that are based on fitting the null model.
}
\examples{
\dontrun{
dat <- expand.table(LSAT7)
(mod1 <- mirt(dat, 1))
M2(mod1)

}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
Cai, L. & Hansen, M. (2013). Limited-information
goodness-of-fit testing of hierarchical item factor models.
British Journal of Mathematical and Statistical Psychology,
66, 245-276.

Maydeu-Olivares, A. & Joe, H. (2006). Limited information
goodness-of-fit testing in multidimensional contingency
tables Psychometrika, 71, 713-732.
}
\keyword{fit}
\keyword{model}

