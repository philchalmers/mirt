\docType{methods}
\name{print}
\alias{print}
\alias{print,ConfirmatoryClass-method}
\alias{print,ExploratoryClass-method}
\alias{print,MixedClass-method}
\alias{print,MultipleGroupClass-method}
\title{Print the model objects}
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
print(x)
}
}
\seealso{
  \code{\link{show}}
}

