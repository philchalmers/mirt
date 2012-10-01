\name{itemplot}
\alias{itemplot}
\title{Displays item surface and information plots}
\usage{
  itemplot(object, item, type = "trace", degrees = 45, ...)
}
\arguments{
  \item{object}{a computed model object of class
  \code{ExploratoryClass} or \code{ConfirmatoryClass}}

  \item{item}{a single numeric value indicating which item
  to plot}

  \item{type}{plot type to use, information
  (\code{'info'}), information contours
  \code{('infocontour')}, or item trace lines
  (\code{'trace'})}

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

itemplot(mod1, 2)
itemplot(mod1, 2, type = 'info')
    }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{plot}

