\name{itemplot}
\alias{itemplot}
\title{Displays item surface and information plots}
\usage{
  itemplot(object, item, type = "trace", degrees = 45, ...)
}
\arguments{
  \item{object}{a computed model object of class
  \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
  \code{MultipleGroupClass}. Input may also be a
  \code{list} for comparing similar item types (e.g., 1PL
  vs 2PL)}

  \item{item}{a single numeric value indicating which item
  to plot}

  \item{type}{plot type to use, information
  (\code{'info'}), item trace lines (\code{'trace'}),
  relative efficiency lines (\code{'RE'}), expected score
  \code{'score'}, or information contours
  \code{('infocontour')} (not for \code{MultipleGroupClass}
  objects)}

  \item{degrees}{the degrees argument to be used if there
  are exactly two factors. See \code{\link{iteminfo}} for
  more detail}

  \item{...}{additional arguments to be passed to lattice}
}
\description{
  \code{itemplot} displays various item based IRT plots.
}
\examples{
\dontrun{

data(LSAT7)
fulldata <- expand.table(LSAT7)
mod1 <- mirt(fulldata,1)
mod2 <- mirt(fulldata,1, itemtype = '1PL')
mod3 <- mirt(fulldata,2)

itemplot(mod1, 2)
itemplot(mod1, 2, type = 'info')

mods <- list(twoPL = mod1, onePL = mod2)
itemplot(mods, 1, type = 'RE')

itemplot(mod3, 3, type = 'info')


    }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{plot}

