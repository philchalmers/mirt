\docType{methods}
\name{residuals-method}
\alias{residuals,ConfirmatoryClass-method}
\alias{residuals,ExploratoryClass-method}
\alias{residuals,MultipleGroupClass-method}
\alias{residuals-method}
\title{Compute model residuals}
\arguments{
  \item{object}{an object of class \code{ExploratoryClass},
  \code{ConfirmatoryClass} or \code{MultipleGroupClass}}

  \item{type}{type of residuals to be displayed.  Can be
  either \code{'LD'} or \code{'LDG2'} for a local
  dependence matrix based on the X2 or G2 statistics (Chen
  & Thissen, 1997), \code{'Q3'} for the statistic proposed
  by Yen (1984), or \code{'exp'} for the expected values
  for the frequencies of every response pattern}

  \item{tables}{logical; for LD type, return the observed,
  expected, and standardized residual tables for each item
  combination?}

  \item{digits}{number of significant digits to be rounded}

  \item{df.p}{logical; print the degrees of freedom and
  p-values?}

  \item{full.scores}{logical; compute relevant statistics
  for each subject in the original data?}

  \item{printvalue}{a numeric value to be specified when
  using the \code{res='exp'} option. Only prints patterns
  that have standardized residuals greater than
  \code{abs(printvalue)}. The default (NULL) prints all
  response patterns}

  \item{verbose}{logical; allow information to be printed
  to the console?}

  \item{Theta}{a matrix of factor scores used for
  statistics that require empirical estimates (i.e., Q3).
  If supplied, arguments typically passed to
  \code{fscores()} will be ignored and these values will be
  used instead}

  \item{...}{additional arguments to be passed to
  \code{fscores()}}
}
\description{
  \code{residuals(object, type = 'LD', digits = 3, df.p =
  FALSE, full.scores = FALSE, printvalue = NULL, tables =
  FALSE, verbose = TRUE, Theta = NULL, ...)}
}
\examples{
\dontrun{

x <- mirt(Science, 1)
residuals(x)
residuals(x, tables = TRUE)
residuals(x, type = 'exp')

# with and without supplied factor scores
Theta <- fscores(x, full.scores=TRUE, scores.only=TRUE)
residuals(x, type = 'Q3', Theta=Theta)
residuals(x, type = 'Q3', method = 'ML')

}
}
\references{
  Chen, W. H. & Thissen, D. (1997). Local dependence
  indices for item pairs using item response theory.
  \emph{Journal of Educational and Behavioral Statistics,
  22}, 265-289.

  Yen, W. (1984). Effects of local item dependence on the
  fit and equating performance of the three parameter
  logistic model. \emph{Applied Psychological Measurement,
  8}, 125-145.
}

