\docType{methods}
\name{residuals-method}
\alias{residuals,ConfirmatoryClass-method}
\alias{residuals,ExploratoryClass-method}
\alias{residuals,MultipleGroupClass-method}
\alias{residuals-method}
\title{Compute model residuals}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass} or \code{MultipleGroupClass}}

  \item{restype}{type of residuals to be displayed.  Can be
  either \code{'LD'} for a local dependence matrix (Chen &
  Thissen, 1997) or \code{'exp'} for the expected values
  for the frequencies of every response pattern}

  \item{tables}{logical; for LD restype, return the
  observed, expected, and standarized residual tables for
  each item combination?}

  \item{digits}{number of significant digits to be rounded}

  \item{df.p}{logical; print the degrees of freedom and
  p-values?}

  \item{full.scores}{logical; compute relavent statistics
  for each subject in the original data?}

  \item{printvalue}{a numeric value to be specified when
  using the \code{res='exp'} option. Only prints patterns
  that have standardized residuals greater than
  \code{abs(printvalue)}. The default (NULL) prints all
  response patterns}

  \item{verbose}{logical; allow information to be printed
  to the console?}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{residuals(object, restype = 'LD', digits = 3, df.p
  = FALSE, full.scores = FALSE, printvalue = NULL, tables =
  FALSE, verbose = TRUE, ...)}
}
\examples{
\dontrun{
x <- mirt(Science, 1)
residuals(x)
residuals(x, tables = TRUE)
residuals(x, restype = 'exp')
}
}

