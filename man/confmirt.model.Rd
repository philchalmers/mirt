\name{confmirt.model}
\alias{confmirt.model}
\title{Specify confirmatory loadings}
\usage{
  confmirt.model(file = "", ...)
}
\arguments{
  \item{file}{a string specifying an external file that
  declares the input.}

  \item{...}{additional arguments for \code{scan()}}
}
\value{
  Returns a model specification object to be used in
  \code{\link{confmirt}}.
}
\description{
  The \code{confmirt.model} function scans user input to
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
  'return' key).

  There are six optional keyword for specifying the
  relationships between factors, item parameters, and prior
  distributions. Note that all independent terms must be
  separated with a comma. These are:

  \describe{ \item{COV}{Specify the relationship between
  the latent factors. Estimating a correlation between
  factors is declared by joining the two factors with an
  asterisk (e.g., F1*F2), while fixing this value to a
  constant is performed by using the keyword 'eq' followed
  by some numeric value (e.g., \code{F1*F2 eq .5} fixes the
  correlation between F1 and F2 to .5).}
  \item{MEAN}{Specify the constant values of the latent
  means. For instance, \code{F1 eq 1} gives a latent mean
  of 1 for factor F1. By default means are all assumed to
  be fixed at 0.} \item{SLOPE}{Specify fixed and equality
  constraints on and between item slope parameters. These
  constraints require the declaration of the factor,
  followed by @ and the respective item, and then the
  desired constraint. Parameters may be fixed to be equal
  during estimation using 'eq' or may be set to a constant
  value using the same keyword. For instance, \code{F1@2 eq
  .5} fixes the first factor slope on item 2 to be .5, and
  \code{F1@2 eq F2@4 eq F2@5}, fixes three slope parameters
  to be equal during estimation.} \item{INT}{Specify fixed
  and equality relationships between the item intercepts.
  These constraints require the declaration of the item,
  followed by @, then the respective intercept (from lowest
  to highest), and finally the constraint behaviour.
  Parameters may be fixed to be equal during estimation
  using 'eq', or equal to a constant value using the same
  method. For instance, \code{1@2 eq .5} fixes the first
  item's second intercept to be .5, and \code{1@1 eq 2@1 eq
  3@1} fixes the three intercept parameters to be equal
  during estimation.} \item{PARTCOMP}{Specify which items
  should be treated as having a partially compensatory
  factor relationship rather than the traditional
  compensatory relationship. This is specified in the same
  fashion as how items load on factors. Note that use of
  PARTCOMP should be accompanied by items that 'anchor' the
  metric of the factors, otherwise the stochastic estimates
  may have trouble converging or converge to poor values.}
  \item{PRIOR}{Give normal prior distributions on select
  parameters that may be misbehaving during estimation.
  These are specified in the form \code{PRIOR =
  N(mean,sd)@parnumber}. To obtain the appropriate
  parameter number use the \code{returnindex = TRUE} option
  in the \code{confmirt()} function and locate the desired
  parameter.} }
}
\examples{
\dontrun{

model <- confmirt.model()
  F1 = 1,2,3,4-10
  F2 = 10-20
  (F1*F1) = 1,2,3,4-10
  COV = F1*F2
  PARTCOMP = 10
  SLOPE = F1@1 eq F1@2
  MEAN = F1 eq 1, F2 eq -1
  INT = 1@1 eq -0.5
  PRIOR = N(0,4)@1

    }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}

