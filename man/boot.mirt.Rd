\name{boot.mirt}
\alias{boot.mirt}
\title{Calculate bootstrapped standard errors for estimated models}
\usage{
  boot.mirt(x, R = 100, ...)
}
\arguments{
  \item{x}{an estimated object from \code{mirt},
  \code{bfactor}, or \code{multipleGroup}}

  \item{R}{number of draws to use (passed to the
  \code{boot()} function)}

  \item{...}{additional arguments to be passed on to
  \code{boot(...)}}
}
\description{
  Given an internal mirt object estimate the bootstrapped
  standard errors. It may be beneficial to run the
  computations using multicore architecture (e.g., the
  \code{parallel} package).
}
\examples{
\dontrun{

#standard
mod <- mirt(Science, 1)
booted <- boot.mirt(mod)
plot(booted)
booted

#run in parallel using snow backend using all available cores
mod <- mirt(Science, 1)
booted <- boot.mirt(mod, parallel = 'snow', ncpus = parallel::detectCores())
booted


}
}
\seealso{
  \code{\link{PLCI.mirt}}
}
\keyword{bootstrapped}
\keyword{errors}
\keyword{standard}

