\name{expected.test}
\alias{expected.test}
\title{Function to calculate expected test score}
\usage{
expected.test(x, Theta, group = NULL)
}
\arguments{
  \item{x}{an estimated mirt object}

  \item{Theta}{a matrix of latent trait values}

  \item{group}{a number signifying which group the item
  should be extracted from (applies to 'MultipleGroupClass'
  objects only)}
}
\description{
Given an estimated model compute the expected test score.
Returns the expected values in the same form as the data
used to estimate the model.
}
\examples{
\dontrun{
dat <- expand.table(deAyala)
model <- mirt.model('F = 1-5
                    CONSTRAIN = (1-5, a1)')
mod <- mirt(dat, model)

Theta <- matrix(seq(-6,6,.01))
tscore <- expected.test(mod, Theta)
tail(cbind(Theta, tscore))

}
}
\seealso{
\code{\link{extract.item}}, \code{\link{expected.item}}
}
\keyword{expected}
\keyword{score}

