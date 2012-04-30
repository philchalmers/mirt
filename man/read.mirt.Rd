\name{read.mirt}
\alias{read.mirt}
\title{Translate mirt parameters for plink package}
\usage{
  read.mirt(x, loc.out = FALSE, as.irt.pars = TRUE)
}
\arguments{
  \item{x}{an object returned from \code{mirt, bfactor,
  polymirt}, or \code{confmirt}}

  \item{loc.out}{if \code{TRUE}, the step/threshold
  parameters will be reformated to be deviations from a
  location parameter}

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

