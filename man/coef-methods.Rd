\docType{methods}
\name{coef}
\alias{coef}
\alias{coef,ConfirmatoryClass-method}
\alias{coef,ExploratoryClass-method}
\alias{coef,MixedClass-method}
\alias{coef,MultipleGroupClass-method}
\title{Extract raw coefs from model object}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass}, \code{MultipleGroupClass}, or
  \code{MixedClass}}

  \item{CI}{default \code{CI = .95}. The amount of
  converged used to compute confidence intervals; default
  is 95 percent confidence intervals}

  \item{IRTpars}{default \code{IRTpars = FALSE}. Logical;
  convert slope intercept parameters into traditional IRT
  parameters?  Only applicable to unidimensional models}

  \item{rotate}{default \code{rotate = ''}. See
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

  \item{rawug}{default \code{rawug = FALSE}. Logical;
  return the untranformed internal g and u parameters?  If
  \code{FALSE}, g and u's are converted with the original
  format along with delta standard errors}

  \item{...}{additional arguments to be passed}
}
\description{
  Extract raw coefs from model object
}
\examples{
\dontrun{
x <- mirt(Science, 2)
coef(x)
coef(x, rotate = 'varimax')
}
}
\seealso{
  \code{\link{coef}}
}

