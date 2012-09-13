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
data(LSAT7)
data <- expand.table(LSAT7)
cmodel <- confmirt.model()
   F1 = 1,4,5
   F2 = 2,3


mod <- mirt(data, cmodel)
coef(mod, allpars = TRUE)
index <- mirt(data, cmodel, constrain = 'index')
index

#second factor slopes equal to 0?
L <- rep(0, 30)
L[c(7, 12)] <- 1
wald(L, mod)

#last two items same factor 1 slope, and item 2 and 3 same factor slope
L <- matrix(0, 2, 30)
L[1,16] <- L[2, 7] <- 1
L[1,21] <- L[2, 12] <- -1
wald(L, mod)

#logLiklihood tests (requires estimating a new model)
mod2 <- mirt(data, cmodel, constrain = list(c(7,12), c(16,21)))
anova(mod2, mod)
}
}
\keyword{wald}

