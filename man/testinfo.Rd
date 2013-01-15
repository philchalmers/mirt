\name{testinfo}
\alias{testinfo}
\title{Function to calculate test information}
\usage{
  testinfo(x, Theta, degrees = NULL, group = NULL)
}
\arguments{
  \item{x}{an estimated mirt object}

  \item{Theta}{a matrix of latent trait values}

  \item{degrees}{a vector of angles in degrees that are
  between 0 and 90 that jointly sum to 90. Only applicable
  when the input object is multidimensional}

  \item{group}{a number signifying which group the item
  should be extracted from (applies to 'MultipleGroupClass'
  objects only)}
}
\description{
  Given an estimated model compute the test information.
}
\examples{
\dontrun{
dat <- expand.table(deAyala)
mod <- mirt(dat, 1, '1PL')

Theta <- matrix(seq(-4,4,.01))
tinfo <- testinfo(mod, Theta)
plot(Theta, tinfo, type = 'l')

}
}
\keyword{information}

