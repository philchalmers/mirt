\name{itemplot}
\alias{itemplot}
\title{Displays item surface and information plots}
\usage{
  itemplot(object, item, type = "trace", degrees = 45,
    CE = FALSE, CEalpha = 0.05, CEdraws = 1000, ...)
}
\arguments{
  \item{object}{a computed model object of class
  \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
  \code{MultipleGroupClass}. Input may also be a
  \code{list} for comparing similar item types (e.g., 1PL
  vs 2PL)}

  \item{item}{a single numeric value, or the item name,
  indicating which item to plot}

  \item{type}{plot type to use, information
  (\code{'info'}), item trace lines (\code{'trace'}),
  relative efficiency lines (\code{'RE'}), expected score
  \code{'score'}, or information contours
  \code{('infocontour')} (not for \code{MultipleGroupClass}
  objects)}

  \item{degrees}{the degrees argument to be used if there
  are exactly two factors. See \code{\link{iteminfo}} for
  more detail}

  \item{CE}{logical; plot confidence envelope?}

  \item{CEalpha}{area remaining in the tail for confidence
  envelope. Default gives 95\% confidence region}

  \item{CEdraws}{draws number of draws to use for
  confidence envelope}

  \item{...}{additional arguments to be passed to lattice}
}
\description{
  \code{itemplot} displays various item based IRT plots.
}
\examples{
\dontrun{

data(LSAT7)
fulldata <- expand.table(LSAT7)
mod1 <- mirt(fulldata,1,SE=TRUE)
mod2 <- mirt(fulldata,1, itemtype = '1PL')
mod3 <- mirt(fulldata,2)

itemplot(mod1, 2)
itemplot(mod1, 2, CE = TRUE)
itemplot(mod1, 2, type = 'info')
itemplot(mod1, 2, type = 'info', CE = TRUE)

mods <- list(twoPL = mod1, onePL = mod2)
itemplot(mods, 1, type = 'RE')

#multidimensional info
itemplot(mod3, 3, type = 'info')

#polytomous items
pmod <- mirt(Science, 1, SE=TRUE, SE.type = 'MHRM')
itemplot(pmod, 3)
itemplot(pmod, 3, CE = TRUE)

    }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{plot}

