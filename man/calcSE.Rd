\name{calcSE}
\alias{calcSE}
\title{Calculate standard errors for estimated model}
\usage{
  calcSE(x, SE.type = "boot", R = 1000,
    return.boot = FALSE, ...)
}
\arguments{
  \item{x}{an estimated object}

  \item{SE.type}{a character indicating which class of
  standard errors to compute. Can be \code{'boot'} for
  boostrapped standard errors (does not compute the
  information matrix), \code{'BL'} for Bock and Leiberman
  computed standard errors, or \code{'MC'} for Monte Carlo
  compution of the information matrix to obtain stanard
  errors. If \code{SE = TRUE} was used at estimation
  runtime (indicating that MH-RM standard errors were
  computed) then these will be overwritten in the returned
  object.}

  \item{R}{number of draws to use (also passed to the
  \code{boot()} function)}

  \item{return.boot}{logical; return the estimated object
  from the \code{boot} package?}

  \item{...}{additional arguments to be passed}
}
\description{
  Given an internal mirt object estimate the standard
  errors, and, if possible, the parameter information
  matrix.
}
\examples{
\dontrun{
mod <- mirt(Science, 1)
modwithSE <- calcSE(mod)
coef(modwithSE)
}
}
\keyword{errors}
\keyword{standard}

