\docType{methods}
\name{itemplot}
\alias{itemplot}
\alias{itemplot-method}
\alias{itemplot,bfactorClass,numeric-method}
\alias{itemplot,mirtClass,numeric-method}
\alias{itemplot,polymirtClass,numeric-method}
\title{Displays item surface and information plots}
\arguments{
  \item{object}{a computed model of class
  \code{bfactorClass}, \code{mirtClass}, or
  \code{polymirtClass}}

  \item{item}{a single numerical value indicating which
  item to plot}

  \item{type}{type of graphic to plot, may be either
  \code{'curve'} for the response surface curve,
  \code{'info'} for the item information, \code{'contour'}
  for item probability contours, or \code{'infocontour'}
  for information contour plot. Note that
  \code{polymirtClass} objects can only display information
  plots}

  \item{npts}{number of points to use per dimension. A
  higher value will make the graphs look smoother}

  \item{rot}{allows rotation of the 3D graphics. Default is
  \code{list(x=-70,y=30,z=10)}}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{itemplot} displays 3D surface plots if the number
  of factors is 2, or standard 2D plots if the number of
  factors is equal to one.
}
\section{Methods}{
  \describe{ \item{itemplot}{\code{signature(object =
  "bfactorClass", item = "numeric")}}
  \item{itemplot}{\code{signature(object = "mirtClass",
  item = "numeric")}}
  \item{itemplot}{\code{signature(object = "polymirtClass",
  item = "numeric")}} }
}
\examples{
\dontrun{

data(LSAT7)
fulldata <- expand.table(LSAT7)
mod1 <- mirt(fulldata,1)
mod2 <- mirt(fulldata,2)

#first item unidimensional
itemplot(mod1, 1)
itemplot(mod1, 1, type = 'info')

#first item multidimensional
itemplot(mod2, 1)
itemplot(mod2, 1, type = 'info')
#turn it a little
itemplot(mod2, 1, type = 'info', rot = list(x = -70, y = -15, z = -10))

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

