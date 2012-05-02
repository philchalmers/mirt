\docType{methods}
\name{itemplot}
\alias{itemplot}
\alias{itemplot-method}
\alias{itemplot,bfactorClass-method}
\alias{itemplot,mirtClass-method}
\alias{itemplot,polymirtClass-method}
\title{Displays item surface and information plots}
\arguments{
  \item{object}{a computed model of class
  \code{bfactorClass}, \code{mirtClass}, or
  \code{polymirtClass}}

  \item{...}{additional arguments to be passed on to
  \code{\link[plink]{plink}} generic \code{plot()}. See the
  \code{\link[plink]{plink}} package for further details.}
}
\description{
  \code{itemplot} displays various item based IRT plots.
}
\section{Methods}{
  \describe{ \item{itemplot}{\code{signature(object =
  "bfactorClass")}} \item{itemplot}{\code{signature(object
  = "mirtClass")}} \item{itemplot}{\code{signature(object =
  "polymirtClass")}} }
}
\examples{
\dontrun{

data(LSAT7)
fulldata <- expand.table(LSAT7)
mod1 <- mirt(fulldata,1)
mod2 <- mirt(fulldata,2)

itemplot(mod1)
itemplot(mod1, combine = 5, auto.key=list(space="right"))

itemplot(mod2, drape = TRUE)
itemplot(mod2, type = "vectorplot1")

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

