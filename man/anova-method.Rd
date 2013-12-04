\docType{methods}
\name{anova-method}
\alias{anova,ConfirmatoryClass-method}
\alias{anova,ExploratoryClass-method}
\alias{anova,MixedClass-method}
\alias{anova,MultipleGroupClass-method}
\alias{anova-method}
\title{Compare nested models}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass}, \code{MultipleGroupClass}, or
  \code{MixedClass}}

  \item{object2}{a second model estimated from any of the
  mirt package estimation methods}
}
\description{
  \code{anova(object, object2)}
}
\examples{
\dontrun{
x <- mirt(Science, 1)
x2 <- mirt(Science, 2)
anova(x, x2)
}
}

