\name{boot.mirt}
\alias{boot.mirt}
\title{Calculate bootstrapped standard errors for estimated models}
\usage{
  boot.mirt(x, R = 1000, return.boot = TRUE, ...)
}
\arguments{
  \item{x}{an estimated object from \code{mirt},
  \code{bfactor}, or \code{multipleGroup}}

  \item{R}{number of draws to use (passed to the
  \code{boot()} function)}

  \item{return.boot}{logical; return the estimated object
  from the \code{boot} package? If \code{FALSE} the
  estimated model is returned with the bootstrapped
  standard errors}

  \item{...}{additional arguments to be passed on to
  \code{boot(...)}}
}
\description{
  Given an internal mirt object estimate the bootstrapped
  standard errors. If possible, it will be benifitial to
  run the computations using multicore architecher (e.g.,
  using the \code{parallel} package).
}
\examples{
\dontrun{
mod <- mirt(Science, 1)
booted <- boot.mirt(mod, R = 100)
booted
modwithboot <- boot.mirt(mod, R = 100, return.boot = FALSE)
coef(modwithboot)

}
}
\keyword{bootstrapped}
\keyword{errors}
\keyword{standard}

