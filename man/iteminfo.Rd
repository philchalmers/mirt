\name{iteminfo}
\alias{iteminfo}
\title{Function to calculate item information}
\usage{
  iteminfo(x, Theta, degrees = NULL)
}
\arguments{
  \item{x}{an extracted internal mirt object containing
  item information}

  \item{Theta}{a matrix of latent trait values}

  \item{degrees}{a vector of angles in degrees that are
  between 0 and 90 that jointly sum to 90. Only applicable
  when the input object is multidimensional}
}
\description{
  Given an internal mirt object extracted from an estimated
  model compute the item information.
}
\examples{
\dontrun{
mod <- mirt(Science, 1, SE = FALSE)
extr.2 <- extract.item(mod, 2)
Theta <- matrix(seq(-4,4, by = .1))
info.2 <- iteminfo(extr.2, Theta)

#do something with the info?
plot(Theta, info.2, type = 'l', Main = 'Item information')
}
}
\keyword{information}

