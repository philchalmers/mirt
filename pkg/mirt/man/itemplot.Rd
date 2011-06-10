\name{itemplot}
\alias{itemplot}
\alias{itemplot,mirtClass-method}
\alias{itemplot,bfactorClass-method}
\alias{itemplot,polymirtClass-method}

\title{Displays item surface and information plots}
\description{
\code{itemplot} displays 3D surface plots if the number of factors is 2, or standard 2D plots if the number of factors is equal to one.
}

\usage{
itemplot(object, ...)

\S4method{itemplot}{mirt}(object, item, type = 'info',
  npts = 50, rot = list(x = -70, y = 30, z = 10), ...)
\S4method{itemplot}{bfactor}(object, item, type = 'curve',
  npts = 50, rot = list(x = -70, y = 30, z = 10), ...)
\S4method{itemplot}{polymirt}(object, item,
  npts = 50, rot = list(x = -70, y = 30, z = 10), ...)
}
 \arguments{
   \item{object}{a computed model of class \code{bfactor} or \code{mirt}}
   \item{item}{a single numerical value indicating which item to plot}
   \item{type}{type of graphic to plot, may be either \code{'curve'} for the response surface curve or \code{'info'} for the item information}
   \item{npts}{number of points to use per dimension. A higher value will make the graphs look smoother}
   \item{rot}{allows rotation of the 3D graphics}
   \item{...}{additional arguments to be passed}
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
\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\seealso{\code{\link{plot.mirt}}}

\keyword{plot}
