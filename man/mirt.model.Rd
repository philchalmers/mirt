\name{mirt.model}
\alias{mirt.model}
\title{Specify model loadings}
\usage{
  mirt.model(input = NULL, file = "", COV = NULL, ...)
}
\arguments{
  \item{input}{input for writing out the model syntax. Can
  either be a string declaration of class character or the
  so-called Q-matrix or class \code{matrix} that specifies
  the model either with integer or logical values. If the
  Q-matrix method is chosen covariances terms can be
  sepcified with the \code{COV} input}

  \item{file}{a input specifying an external file that
  declares the input.}

  \item{COV}{a symmetric, logical matrix used to declare
  which covariance terms are estimated}

  \item{...}{additional arguments for \code{scan()}}
}
\value{
  Returns a model specification object to be used in
  \code{\link{confmirt}}.
}
\description{
  The \code{mirt.model} function scans/reads user input to
  specify the confirmatory model.
}
\details{
  Factors are first named and then specify which numerical
  items they affect (i.e., where the slope is not equal to
  0), separated either by commas or by - to indicate a
  range of items. Products between factors may be specified
  by enclosing the left hand term within brackets. To
  finish the declaration of a model simply enter a blank
  line with only a carriage return (i.e., the 'enter' or
  'return' key), or instead read in a input version of the
  model syntax.

  There is an optional keyword for specifying the
  correlation between relationships between factors called
  \code{COV}, and non-linear factor products can be
  included by enclosing the product combination on the left
  hand side of the declaration (e.g., \code{(F1*F1)} would
  create a quadratic factor for \code{F1}).

  \describe{ \item{COV}{Specify the relationship between
  the latent factors. Estimating a correlation between
  factors is declared by joining the two factors with an
  asterisk (e.g., F1*F2).} }
}
\examples{
\dontrun{

model <- mirt.model()
  F1 = 1,2,3,4-10
  F2 = 10-20
  (F1*F2) = 1,2,3,4-10
  COV = F1*F2


#Or alternatively
s <- 'F1 = 1,2,3,4-10
      F2 = 10-20
      (F1*F2) = 1,2,3,4-10
      COV = F1*F2'
model <- mirt.model(s)


#Q-matrix specification
Q <- matrix(c(1,1,1,0,0,0,0,0,0,1,1,1), ncol=2, dimnames = list(NULL, c('Factor1', 'Factor2')))
COV <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2)
model <- mirt.model(Q, COV=COV)

    }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}

