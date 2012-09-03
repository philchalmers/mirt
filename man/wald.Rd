\name{wald}
\alias{print.wald}
\alias{wald}
\title{Wald test for mirt models}
\usage{
  wald(L, object, C = 0)

  \method{print}{wald} (x, ...)
}
\arguments{
  \item{L}{a coefficient matrix with dimensions nconstrasts
  x npars. Use \code{constrain = 'index'} on the initially
  estimated model to obtain the parameter indicators}

  \item{object}{estimated object from mirt, confmirt, or
  multipleGroup}

  \item{C}{a constant vector to be compared along side L.
  Default is 0}

  \item{x}{an object of class 'wald'}

  \item{...}{additional arguments to be passed}
}
\description{
  Compute a Wald test given an \code{L} vector or matrix of
  contrasts.
}
\examples{
\dontrun{
#View parnumber index
mirt(Science, 2, constrain = 'index')
mod <- mirt(Science, 2)

#all first factor slopes equal to 0
L <- rep(0, 25)
L[c(1,6,11,16)] <- 1
wald(L, mod)

#first two items same factor 1 slope, last two items same factor 1 slope
L <- matrix(0, 2, 25)
L[1,1] <- L[2, 11] <- 1
L[1,6] <- L[2, 16] <- -1
wald(L, mod)
}
}
\keyword{wald}

