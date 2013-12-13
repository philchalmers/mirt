\docType{methods}
\name{plot-method}
\alias{plot,ConfirmatoryClass-method}
\alias{plot,ExploratoryClass-method}
\alias{plot,MultipleGroupClass-method}
\alias{plot-method}
\title{Plot various test implied functions from models}
\arguments{
  \item{x}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass} or \code{MultipleGroupClass}}

  \item{type}{type of plot to view; can be \code{'info'} to
  show the test information function, \code{'infocontour'}
  for the test information contours, \code{'SE'} for the
  test standard error function, \code{'trace'} and
  \code{'infotrace'} for all item probability information
  or trace lines (only available when all items are
  dichotomous), \code{'infoSE'} for a combined test
  information and standard error plot, and \code{'score'}
  for the expected total score. If \code{empiricalhist =
  TRUE} was used then the type \code{'empiricalhist'} also
  will be available to generate the empirical histogram
  plot}

  \item{theta_angle}{numeric values ranging from 0 to 90
  used in \code{plot}.  If a vector is used then a bubble
  plot is created with the summed information across the
  angles specified (e.g., \code{theta_angle = seq(0, 90,
  by=10)})}

  \item{npts}{number of quadrature points to be used for
  plotting features.  Larger values make plots look
  smoother}

  \item{rot}{allows rotation of the 3D graphics}

  \item{which.items}{numeric vector indicating which items
  to be used when plotting. Default is to use all available
  items}

  \item{facet_items}{logical; apply grid of plots accross
  items or groups? Applicable only to
  \code{multipleGroupClass} objects. Default is
  \code{FALSE}}

  \item{auto.key}{logical parameter passed to the
  \code{lattice} package}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{plot(x, y, type = 'info', npts = 50, theta_angle =
  45, which.items = 1:ncol(x@data), rot = list(xaxis = -70,
  yaxis = 30, zaxis = 10), auto.key = TRUE, ...)}
}
\examples{
\dontrun{
x <- mirt(Science, 1)
plot(x)
plot(x, type = 'trace')
}
}

