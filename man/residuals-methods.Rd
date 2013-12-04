\docType{methods}
\name{residuals}
\alias{residuals}
\alias{residuals,ConfirmatoryClass-method}
\alias{residuals,ExploratoryClass-method}
\alias{residuals,MultipleGroupClass-method}
\title{Compute model residuals}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass} or \code{MultipleGroupClass}}

  \item{restype}{default \code{restype = 'LD'}. Type of
  residuals to be displayed.  Can be either \code{'LD'} for
  a local dependence matrix (Chen & Thissen, 1997) or
  \code{'exp'} for the expected values for the frequencies
  of every response pattern}

  \item{digits}{default \code{digits = 3}. Number of
  significant digits to be rounded}

  \item{df.p}{default \code{df.p = FALSE}. logical; print
  the degrees of freedom and p-values?}

  \item{full.scores}{default \code{full.scores = FALSE}.
  logical; compute relavent statistics for each subject in
  the original data?}

  \item{printvalue}{default \code{printvalue = NULL}. A
  numeric value to be specified when using the
  \code{res='exp'} option. Only prints patterns that have
  standardized residuals greater than
  \code{abs(printvalue)}. The default (NULL) prints all
  response patterns}

  \item{verbose}{default \code{verbose = TRUE}. Logical;
  allow information to be printed to the console?}

  \item{...}{additional arguments to be passed}
}
\description{
  Compute model residuals
}
\examples{
\dontrun{
x <- mirt(Science, 1)
residuals(x)
residuals(x, restype = 'exp')
}
}
\seealso{
  \code{\link{residuals}}
}

