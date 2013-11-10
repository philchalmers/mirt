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
  Multidimensional extensions are currently not supported.
}
\examples{
\dontrun{

#unidimensional

data <- expand.table(LSAT7)
(mod1 <- mirt(data, 1))
plinkpars <- read.mirt(mod1)
plot(plinkpars)
itemplot(mod1, 1)

#graded
mod2 <- mirt(Science, 1)
plinkpars <- read.mirt(mod2)
plot(plinkpars)
itemplot(mod2, 1)

#gpcm
mod3 <- mirt(Science, 1, itemtype = 'gpcm')
plinkpars <- read.mirt(mod3)
plot(plinkpars)
itemplot(mod3, 1)

#nominal
mod4 <- mirt(Science, 1, itemtype = 'nominal')
plinkpars <- read.mirt(mod4)
plot(plinkpars)
itemplot(mod4, 1)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{plink}

