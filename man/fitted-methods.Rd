\docType{methods}
\name{fitted}
\alias{fitted}
\alias{fitted,ConfirmatoryClass-method}
\alias{fitted,ExploratoryClass-method}
\alias{fitted,MultipleGroupClass-method}
\title{Compute fitted values}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass}, or \code{MultipleGroupClass}}

  \item{digits}{default \code{digits = 3}. Number of
  significant digits to be rounded}

  \item{...}{additional arguments to be passed}
}
\description{
  Compute fitted values
}
\examples{
\dontrun{
x <- mirt(Science, 1)
fitted(x)
}
}
\seealso{
  \code{\link{fitted}}
}

