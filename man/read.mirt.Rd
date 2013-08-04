\name{read.mirt}
\alias{read.mirt}
\title{Translate mirt parameters for plink package}
\usage{
  read.mirt(x, as.irt.pars = TRUE, ...)
}
\arguments{
  \item{x}{an object returned from \code{mirt, bfactor}, or
  \code{multipleGroup}}

  \item{as.irt.pars}{if \code{TRUE}, the parameters will be
  output as an \code{irt.pars} object}

  \item{...}{additional arguments to be passed to
  \code{coef()}}
}
\description{
  This function exports item parameters from the
  \code{mirt} package to the \code{plink} package.
}
\examples{
\dontrun{
data <- expand.table(LSAT7)
(mod1 <- mirt(data, 1))
plinkpars <- read.mirt(mod1)

(mod2 <- mirt(data, 2))
plinkpars2 <- read.mirt(mod2, rotate = 'none')

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{plink}

