\name{PLCI.mirt}
\alias{PLCI.mirt}
\title{Compute profiled-likelihood confidence intervals}
\usage{
  PLCI.mirt(mod, alpha = 0.05)
}
\arguments{
  \item{mod}{a converged mirt model}

  \item{alpha}{two-tailed alpha critical level}
}
\description{
  Computes profiled-likelihood based confidence intervals.
  Supports the inclusion of prior parameter distributions
  as well as equality constraints. Currently only supports
  unidimensional dichotomous response models.
}
\examples{
\dontrun{
mirtCluster() #use all available cores to estimate CI's in parallel
dat <- expand.table(LSAT7)
mod <- mirt(dat, 1)

result <- PLCI.mirt(mod)
result

}
}
\seealso{
  \code{\link{boot.mirt}}
}
\keyword{likelihood}
\keyword{profiled}

