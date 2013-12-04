\docType{methods}
\name{anova}
\alias{anova}
\alias{anova,ConfirmatoryClass-method}
\alias{anova,ExploratoryClass-method}
\alias{anova,MixedClass-method}
\alias{anova,MultipleGroupClass-method}
\title{Compare nested models}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass}, \code{MultipleGroupClass}, or
  \code{MixedClass}}

  \item{object2}{a second model estimated from any of the
  mirt package estimation methods}
}
\description{
  Compare nested models
}
\examples{
\dontrun{
x <- mirt(Science, 1)
x2 <- mirt(Science, 2)
anova(x, x2)
}
}
\seealso{
  \code{\link{anova}}
}

