\name{key2binary}
\alias{key2binary}
\title{Convert response patterns to binary data}
\description{
The \code{key2binary} function will convert response pattern data to a dichotomous format, given a response key.
}

\usage{
key2binary(fulldata, key)
}
 \arguments{
   \item{fulldata}{an object of class \code{data.frame}, \code{matrix}, or \code{table} with the response patterns}
   \item{key}{a vector consisting of the 'correct' response to the items. Each value corresponds to each column in \code{fulldata}}
 }

\value{
  returns a numeric matrix with all the response patterns in dichotomous format.
}

\examples{
\dontrun{
data(SAT12)
key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)

dicho.SAT12 <- key2binary(SAT12,key)
    }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}


