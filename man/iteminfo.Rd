\name{iteminfo}
\alias{iteminfo}
\title{Function to calculate item information}
\usage{
  iteminfo(x, Theta, degrees = NULL)
}
\arguments{
  \item{x}{an extracted internal mirt object containing
  item information}

  \item{Theta}{a vector (unidimensional) or matrix
  (multidimensional) of latent trait values}

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
mod <- mirt(Science, 1)
extr.2 <- extract.item(mod, 2)
Theta <- matrix(seq(-4,4, by = .1))
info.2 <- iteminfo(extr.2, Theta)

#do something with the info?
plot(Theta, info.2, type = 'l', main = 'Item information')

## Customized test information plot
T1 <- T2 <- 0
dat <- expand.table(LSAT7)
mod1 <- mirt(dat, 1)
mod2 <- mirt(dat, 1, 'Rasch')
for(i in 1:5){
  T1 <- T1 + iteminfo(extract.item(mod1, i), Theta)
  T2 <- T2 + iteminfo(extract.item(mod2, i), Theta)
}
plot(Theta, T2/T1, type = 'l', ylab = 'Relative Test Information', las = 1)
lines(Theta, T1/T1, col = 'red')

}
}
\seealso{
  \code{\link{extract.item}}
}
\keyword{information}

