\name{read.mirt}
\alias{read.mirt}
\title{Translate mirt parameters for plink package}
\usage{
  read.mirt(x, as.irt.pars = TRUE)
}
\arguments{
  \item{x}{an object returned from \code{mirt, bfactor,
  confmirt}, or \code{multipleGroup}}

  \item{as.irt.pars}{if \code{TRUE}, the parameters will be
  output as an \code{irt.pars} object}
}
\description{
  This function exports item parameters from the
  \code{mirt} package to the \code{\link{plink}} package.
}
\examples{
\dontrun{
data(LSAT7)
data <- expand.table(LSAT7)
(mod1 <- mirt(data, 1))
plinkpars <- read.mirt(mod1)

(mod2 <- mirt(data, 2))
plinkpars2 <- read.mirt(mod2)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{plink}

