\name{extract.item}
\alias{extract.item}
\title{Extract an item object from mirt objects}
\usage{
  extract.item(x, item, group = NULL, drop.zeros = FALSE)
}
\arguments{
  \item{x}{mirt model of class 'ExploratoryClass',
  'ConfirmatoryClass', or 'MultipleGroupClass'}

  \item{item}{a number or character signifying which item
  to extract}

  \item{group}{a number signifying which group the item
  should be extracted from (applies to 'MultipleGroupClass'
  only)}

  \item{drop.zeros}{logical; drop slope values that are
  numerically close to zero to reduce dimensionality?
  Useful in objects returned from \code{\link{bfactor}} or
  other confirmatory models that contain several zero
  slopes}
}
\description{
  Extract the internal mirt objects from any estimated
  model.
}
\examples{
\dontrun{
mod <- mirt(Science, 1)
extr.1 <- extract.item(mod, 1)
}
}
\keyword{extract}

