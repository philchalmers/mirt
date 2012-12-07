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
  x npars, or a vector if only one set of contrasts is
  being tested. Use \code{pars = 'values'} on the initially
  estimated model to obtain the parameter indicators}

  \item{object}{estimated object from \code{mirt},
  \code{bfactor}, \code{confmirt}, or \code{multipleGroup}}

  \item{C}{a constant vector/matrix to be compared along
  side L}

  \item{x}{an object of class 'wald'}

  \item{...}{additional arguments to be passed}
}
\description{
  Compute a Wald test given an \code{L} vector or matrix of
  numeric contrasts. Requires that the model information
  matrix be computed (including \code{SE = TRUE} when using
  the EM method).
}
\examples{
\dontrun{
#View parnumber index
data(LSAT7)
data <- expand.table(LSAT7)
cmodel <- confmirt.model()
   F1 = 1,4,5
   F2 = 2,3


mod <- mirt(data, cmodel, SE = TRUE)
coef(mod, allpars = TRUE)
index <- mirt(data, cmodel, pars = 'values')
index


#second factor slopes equal to 0?
L <- rep(0, 30)
names(L) <- index$name
L[index$name == 'a2' & index$est] <- 1
wald(L, mod)

#simultaneously test equal factor slopes for item 2 and 3, and 4 and 5
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

