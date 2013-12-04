\docType{methods}
\name{fitted-method}
\alias{fitted,ConfirmatoryClass-method}
\alias{fitted,ExploratoryClass-method}
\alias{fitted,MultipleGroupClass-method}
\alias{fitted-method}
\title{Compute fitted values}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass}, or \code{MultipleGroupClass}}

  \item{digits}{number of significant digits to be rounded}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{fitted(object, digits = 3, ...)}
}
\examples{
\dontrun{
x <- mirt(Science, 1)
fitted(x)
}
}

