\name{tetrachor}
\alias{tetrachor}
\title{Computes a tetrachoric matrix with Carroll's correction for guessing}
\description{
The \code{tetrachor} function computes a matrix of tetrachoric correlations between the variables, and if applicable computes the tetrachoric correlations on chance adjusted tables specified by Carroll (1945). If the matrix in not positive definite (and \code{smooth = TRUE}) then it is smoothed by setting the negative eigenvalues to zero and recomputing the correlation matrix. 
}

\usage{
tetrachor(fulldata, guess = 0, nowarn = TRUE)
}
 \arguments{
   \item{fulldata}{dichotomous data matrix}
   \item{guess}{pseudo-guessing parameter; if a single value is specified then the value is repeated for each item}
   \item{smooth}{logical; smooth the matrix to be positive definite?}
 }

\value{
  returns a square tetrachoric correlation matrix. If smoothed the trace of the matrix may not sum to the number of variables. Any values that were initially \code{NA} are recoded to 0.
}

\examples{
    \dontrun{
#use the computed matrix in model estimation
data(LSAT7)
LSAT7full <- expand.table(LSAT7)
tmat <- tetrachor(LSAT7full)

mirt(LSAT7full, 2, prev.cor = tmat)
   }
}

\references{

Carroll, J. B. (1945). The effect of difficulty and chance success on correlations between items and between tests. \emph{Psychometrika, 26}, 347-372.

}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{cor}
