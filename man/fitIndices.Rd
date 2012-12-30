\name{fitIndices}
\alias{fitIndices}
\title{Compute Extra Model Fit Indices}
\usage{
  fitIndices(obj)
}
\arguments{
  \item{obj}{an estimated model object from the mirt
  package}
}
\description{
  Compute additional model fit indecies that do not come as
  direct results following parameter convergence. Will only
  compute the M2 (Maydeu-Olivares & Joe, 2006) statistic by
  default, and returns a list containing the requested
  statistics.
}
\examples{
\dontrun{
#LSAT6 example
dat <- expand.table(LSAT6)
(mod1 <- mirt(dat, 1, itemtype = '1PL'))
fitIndices(mod1)

#Science data, much more sparce so M2 would be more informative
(mod2 <- mirt(Science, 1))
fitIndices(mod2)
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

