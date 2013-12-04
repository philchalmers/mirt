\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,ConfirmatoryClass-method}
\alias{plot,ExploratoryClass-method}
\alias{plot,MultipleGroupClass-method}
\title{Plot various test implied functions from models}
\arguments{
  \item{x}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass} or \code{MultipleGroupClass}}

  \item{type}{default \code{type = 'info'}. Type of plot to
  view; can be \code{'info'} to show the test information
  function, \code{'infocontour'} for the test information
  contours, \code{'SE'} for the test standard error
  function, \code{'trace'} and \code{'infotrace'} for all
  item probability information or trace lines (only
  available when all items are dichotomous),
  \code{'infoSE'} for a combined test information and
  standard error plot, and \code{'score'} for the expected
  total score. If \code{empiricalhist = TRUE} was used then
  the type \code{'empiricalhist'} also will be available to
  generate the empirical histogram plot}

  \item{theta_angle}{default \code{theta_angle = 45}.
  Numeric values ranging from 0 to 90 used in \code{plot}.
  If a vector is used then a bubble plot is created with
  the summed information across the angles specified (e.g.,
  \code{theta_angle = seq(0, 90, by=10)})}

  \item{npts}{default \code{npts = 50}. Number of
  quadrature points to be used for plotting features.
  Larger values make plots look smoother}

  \item{rot}{default \code{rot = list(xaxis = -70, yaxis =
  30, zaxis = 10)}. Allows rotation of the 3D graphics}

  \item{which.items}{default \code{which.items =
  1:ncol(x@data)}.  A numeric vector indicating which items
  to be used when plotting. Default is to use all available
  items}

  \item{auto.key}{default \code{auto.key = TRUE}. Logical
  parameter passed to the \code{lattice} package}

  \item{...}{additional arguments to be passed}
}
\description{
  Plot various test implied functions from models
}
\examples{
\dontrun{
x <- mirt(Science, 1)
plot(x)
plot(x, type = 'trace')
}
}
\seealso{
  \code{\link{plot}}
}

