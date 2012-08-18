\docType{methods}
\name{itemplot}
\alias{itemplot}
\alias{itemplot,bfactorClass-method}
\alias{itemplot,mirtClass-method}
\alias{itemplot-method}
\title{Displays item surface and information plots}
\arguments{
  \item{object}{a computed model of class
  \code{bfactorClass} or \code{mirtClass}}

  \item{type}{plot type to use, information
  (\code{'info'}), information contours
  \code{('infocontour')}, or item trace lines
  (\code{'trace'})}

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

itemplot(mod1)
    }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{plot}}, \code{\link{mirt}},
  \code{\link{bfactor}}, \code{\link{polymirt}}
}
\keyword{plot}

