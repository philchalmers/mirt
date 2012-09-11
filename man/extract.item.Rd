\name{extract.item}
\alias{extract.item}
\title{Extract an item object from mirt objects}
\usage{
  extract.item(x, item, group = NULL)
}
\arguments{
  \item{x}{mirt model of class 'ExploratoryClass',
  'ConfirmatoryClass', or 'MultipleGroupClass'}

  \item{item}{a number signifying which item to extract}

  \item{group}{a number signifying which group the item
  should be extracted from (applies to 'MultipleGroupClass'
  only)}
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

