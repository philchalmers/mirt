\name{mirt.model}
\alias{mirt.model}
\title{Specify model loadings}
\usage{
  mirt.model(input = NULL, file = "", COV = NULL,
    quiet = TRUE, ...)
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

  \item{quiet}{logical argument passed to \code{scan()} to
  suppress console read message}

  \item{...}{additional arguments for \code{scan()}}
}
\value{
  Returns a model specification object to be used in
  \code{\link{mirt}}, \code{\link{multipleGroup}}, or
  \code{\link{mixedmirt}}.
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
  'return' key), or instead read in an input version of the
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
  asterisk (e.g., F1*F2).} \item{CONSTRAIN}{A bracketed,
  comma seperate list specifying equality constrains
  between items. The input format is \code{CONSTRAIN =
  (items, ..., parameterName, OptionalGroup), (items, ...,
  parameterName, OptionalGroup)}. If \code{OptionalGroup}
  is ommited then the constraints are applied within all
  groups. For example, in a single group 10-item dichotmous
  tests, using the default 2PL model, the first and last 5
  item slopes can be constrained to be equal by using
  \code{CONSTRAIN = (1-5, a1), (6-10, a1)}, or some
  cobmination such as \code{CONSTRAIN = (1-3,4,5,a1),
  (6,7,8-10,a1)}} \item{CONSTRAINB}{A bracketed, comma
  seperate list specifying equality constrains between
  groups. The input format is \code{CONSTRAINB = (items,
  ..., parameterName), (items, ..., parameterName)}. For
  example, in a two group 10-item dichotmous tests, using
  the default 2PL model, the first 5 item slopes can be
  constrained to be equal accross both groups by using
  \code{CONSTRAINB = (1-5, a1)}, or some cobmination such
  as \code{CONSTRAINB = (1-3,4,5,a1)}} }
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

## constrain various items slopes and all intercepts in single group model to be equal
s <- 'F = 1-10
      CONSTRAIN = (1-3, 5, 6, a1), (1-10, d)'
model <- mirt.model(s)


## constrain various items slopes and intercepts accross groups for use in multipleGroup(),
#  and constrain first two slopes within 'group1' to be equal
s <- 'F = 1-10
      CONSTRAIN = (1-2, a1)
      CONSTRAINB = (1-3, 5, 6, a1), (1-10, d)'
model <- mirt.model(s)

    }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}

