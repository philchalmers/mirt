\docType{methods}
\name{summary}
\alias{summary}
\alias{summary,ConfirmatoryClass-method}
\alias{summary,ExploratoryClass-method}
\alias{summary,MixedClass-method}
\alias{summary,MultipleGroupClass-method}
\title{Summary of model object}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass}, \code{MultipleGroupClass}, or
  \code{MixedClass}}

  \item{rotate}{default code{rotate = ''}. See
  \code{\link{mirt}} for details}

  \item{Target}{default \code{Target = NULL}. A dummy
  variable matrix indicting a target rotation pattern}

  \item{suppress}{default \code{suppress = 0}. A numeric
  value indicating which (possibly rotated) factor loadings
  should be suppressed. Typical values are around .3 in
  most statistical software. Default is 0 for no
  suppression}

  \item{digits}{default \code{digits = 3}. Number of
  significant digits to be rounded}

  \item{verbose}{default \code{verbose = TRUE}. Logical;
  allow information to be printed to the console?}

  \item{...}{additional arguments to be passed}
}
\description{
  Tranforms coefficients into a standardized factor
  loading's metric. For \code{MixedClass} objects, the
  fixed and random coeffiicents are printed.
}
\examples{
\dontrun{
x <- mirt(Science, 2)
summary(x)
summary(x, rotate = 'varimax')
}
}
\seealso{
  \code{\link{coef}}
}

