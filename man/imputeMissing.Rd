\name{imputeMissing}
\alias{imputeMissing}
\title{Imputing plausible data for missing values}
\usage{
  imputeMissing(x, Theta)
}
\arguments{
  \item{x}{an estimated model x from the mirt package}

  \item{Theta}{a matrix containing the estimates of the
  latent trait scores (e.g., via \code{\link{fscores}})}
}
\description{
  Given an estimated model from any of mirt's model fitting
  functions and an estimate of the latent trait, impute
  plausible missing data values. Returns the original data
  in a \code{data.frame} without any NA values.
}
\examples{
\dontrun{
dat <- expand.table(LSAT7)
(original <- mirt(dat, 1))
NAperson <- sample(1:nrow(dat), 20, replace = TRUE)
NAitem <- sample(1:ncol(dat), 20, replace = TRUE)
for(i in 1:20)
    dat[NAperson[i], NAitem[i]] <- NA
(mod <- mirt(dat, 1))
scores <- fscores(mod, method = 'MAP', full.scores = TRUE, scores.only = TRUE)

#re-estimate imputed dataset (good to do this multiple times and average over)
fulldata <- imputeMissing(mod, scores)
(fullmod <- mirt(fulldata, 1))


}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{data}
\keyword{impute}

