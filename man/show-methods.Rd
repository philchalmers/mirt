\docType{methods}
\name{show}
\alias{show}
\alias{show,ConfirmatoryClass-method}
\alias{show,ExploratoryClass-method}
\alias{show,MixedClass-method}
\alias{show,MultipleGroupClass-method}
\title{Show model object}
\arguments{
  \item{x}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass}, \code{MultipleGroupClass}, or
  \code{MixedClass}}
}
\description{
  Print model object summaries to the console.
}
\examples{
\dontrun{
x <- mirt(Science, 1)
show(x)
}
}
\seealso{
  \code{\link{print}}
}

