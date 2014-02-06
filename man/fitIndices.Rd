\name{fitIndices}
\alias{fitIndices}
\title{Compute Extra Model Fit Indices}
\usage{
  fitIndices(obj, calcNull = FALSE, prompt = TRUE)
}
\arguments{
  \item{obj}{an estimated model object from the mirt
  package}

  \item{calcNull}{logical; calculate statistics for the
  null model as well?  Allows for statistics such as the
  limited information TLI and CFI}

  \item{prompt}{logical; prompt user for input if the
  internal matrices are too large?}
}
\description{
  Compute additional model fit indices that do not come as
  direct results following parameter convergence. Will
  compute the M2 (Maydeu-Olivares & Joe, 2006) statistic by
  default, and returns a data.frame containing various
  model fit statistics.
}
\examples{
\dontrun{
#LSAT6 example
dat <- expand.table(LSAT6)
(mod1 <- mirt(dat, 1, itemtype = '2PL', constrain = list(c(1,5,9,13,17))))
fitIndices(mod1)

#Science data with computing the null model M2 stat
(mod2 <- mirt(Science, 1))
fitIndices(mod2, calcNull = TRUE)
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

