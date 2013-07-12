\name{mirtCluster}
\alias{mirtCluster}
\title{Define a cluster object to be used in mirt functions}
\usage{
  mirtCluster(ncores)
}
\arguments{
  \item{ncores}{number of cores to be used in the returned
  object which is passed to \code{parallel::makeCluster()}}
}
\description{
  This function defines a object that is placed in the
  users workspace (i.e., the\code{.GlobalEnv}). Relavent
  internal functions such as \code{calcLogLik},
  \code{fscores}, etc, will utilize this object
  automatically to capitilze on parallel processing
  architecture. The object defined is a call from
  \code{parallel::makeCluster()}.
}
\examples{
\dontrun{

#make 4 cores available for parallel computing
mirtCluster(4)

#use all available cores
library(parallel)
mirtCluster(detectCores())

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{parallel}

