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
  A plotting function for displaying the individuals
  trajectories and their modelled functional form. Useful
  for detecting aberrant individual trajectories.
}
\examples{
\dontrun{
data(LSAT7)
data <- expand.table(LSAT7)
(mod1 <- mirt(data, 1))
plinkpars <- read.mirt(mod1)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{plink}

