\name{M2}
\alias{M2}
\title{Compute M2 statistic}
\usage{
  M2(obj, calcNull = FALSE)
}
\arguments{
  \item{obj}{an estimated model object from the mirt
  package}

  \item{calcNull}{logical; calculate statistics for the
  null model as well?  Allows for statistics such as the
  limited information TLI and CFI}

  \item{collapse_poly}{logical; collapse across polytomous
  item categories to reduce sparceness? Will also helo to
  reduce the internal matrix sizes. THIS FEATURE IS
  CURRENTLY EXPERIMENTAL AND SHOULD NOT BE TRUSTED}

  \item{prompt}{logical; prompt user for input if the
  internal matrices are too large?}
}
\description{
  Computes the M2 (Maydeu-Olivares & Joe, 2006) statistic
  and associated fit indicies. For now, only dichotomous
  models are supported.
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
  Maydeu-Olivares, A. & Joe, H. (2006). Limited information
  goodness-of-fit testing in multidimensional contingency
  tables Psychometrika, 71, 713-732.
}
\keyword{fit}
\keyword{model}

