\name{probtrace}
\alias{probtrace}
\title{Function to calculate probability trace lines}
\usage{
  probtrace(x, Theta)
}
\arguments{
  \item{x}{an extracted internal mirt object containing
  item information}

  \item{Theta}{a matrix of latent trait values}
}
\description{
  Given an internal mirt object extracted from an estimated
  model compute the probility trace lines for all
  categories.
}
\examples{
\dontrun{
mod <- mirt(Science, 1)
extr.2 <- extract.item(mod, 2)
Theta <- matrix(seq(-4,4, by = .1))
traceline <- probtrace(extr.2, Theta)

head(data.frame(traceline, Theta=Theta))

}
}
\keyword{tracelines}

