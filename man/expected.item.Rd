\name{expected.item}
\alias{expected.item}
\title{Function to calculate expected value of item}
\usage{
expected.item(x, Theta, min = 0)
}
\arguments{
  \item{x}{an extracted internal mirt object containing
  item information}

  \item{Theta}{a vector (unidimensional) or matrix
  (multidimensional) of latent trait values}

  \item{min}{a constant value added to the expected values
  indicating the lowest theoretical category. Default is 0}
}
\description{
Given an internal mirt object extracted from an estimated
model compute the expected value for an item given the
ability parameter(s).
}
\examples{
\dontrun{
mod <- mirt(Science, 1)
extr.2 <- extract.item(mod, 2)
Theta <- fscores(mod, full.scores=TRUE, scores.only=TRUE)
expected <- expected.item(extr.2, Theta, min(Science[,1])) #min() of first item

head(data.frame(expected, Theta=Theta))

}
}
\seealso{
\code{\link{extract.item}}
}
\keyword{expected}
\keyword{value}

