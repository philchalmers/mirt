\name{polymirt}
\alias{polymirt}
\title{Full-Information Item Factor Analysis for Mixed Data Formats}
\usage{
  polymirt(...)
}
\arguments{
  \item{...}{arguments to be passed to the
  \code{\link{confmirt}} estimation engine}
}
\description{
  \code{polymirt} fits an unconditional (exploratory)
  full-information maximum-likelihood factor analysis model
  to dichotomous and polychotomous data under the item
  response theory paradigm using Cai's (2010)
  Metropolis-Hastings Robbins-Monro algorithm. If
  requested, lower and upper asymptote parameters are
  estimated with a beta priors included automatically.
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{expand.table}}, \code{\link{key2binary}},
  \code{\link{confmirt}}, \code{\link{itemplot}}
}

