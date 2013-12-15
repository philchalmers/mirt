\docType{methods}
\name{coef-method}
\alias{coef,ConfirmatoryClass-method}
\alias{coef,ExploratoryClass-method}
\alias{coef,MixedClass-method}
\alias{coef,MultipleGroupClass-method}
\alias{coef-method}
\title{Extract raw coefs from model object}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass}, \code{MultipleGroupClass}, or
  \code{MixedClass}}

  \item{CI}{the amount of converged used to compute
  confidence intervals; default is 95 percent confidence
  intervals}

  \item{IRTpars}{logical; convert slope intercept
  parameters into traditional IRT parameters?  Only
  applicable to unidimensional models}

  \item{rotate}{see \code{\link{mirt}} for details}

  \item{Target}{a dummy variable matrix indicting a target
  rotation pattern}

  \item{suppress}{a numeric value indicating which
  (possibly rotated) factor loadings should be suppressed.
  Typical values are around .3 in most statistical
  software. Default is 0 for no suppression}

  \item{digits}{number of significant digits to be rounded}

  \item{verbose}{logical; allow information to be printed
  to the console?}

  \item{rawug}{logical; return the untranformed internal g
  and u parameters?  If \code{FALSE}, g and u's are
  converted with the original format along with delta
  standard errors}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{coef(object, CI = .95, rotate = '', Target = NULL,
  digits = 3, IRTpars = FALSE, rawug = FALSE, verbose =
  TRUE, ...)}
}
\examples{
\dontrun{
x <- mirt(Science, 2)
coef(x)
coef(x, rotate = 'varimax')
}
}
\seealso{
  \code{\link{summary-method}}
}

