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
  as well as equality constraints. Only supports
  unidimensional models.
}
\examples{
\dontrun{
mirtCluster() #use all available cores to estimate CI's in parallel
dat <- expand.table(LSAT7)
mod <- mirt(dat, 1)

result <- PLCI.mirt(mod)
result

mod2 <- mirt(Science, 1)
result2 <- PLCI.mirt(mod2)
result2

}
}
\seealso{
  \code{\link{boot.mirt}}
}
\keyword{likelihood}
\keyword{profiled}

