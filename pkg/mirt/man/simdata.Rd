\name{simdata}
\alias{simdata}
\title{Simulate response patterns from an underlying multivariate normal distribution}
\description{
Simulates response patterns with multivariate normal distributions.
}

\usage{
simdata(a, d, N, guess = 0, sigma = NULL, mu = NULL)

}
 \arguments{
   \item{a}{a matrix of slope parameters}
   \item{d}{a vector of intercepts}
   \item{N}{sample size}
   \item{guess}{a vector of guessing parameters}
   \item{sigma}{a covariance matrix of the underlying distribution}
   \item{mu}{a mean vector of the underlying distribution}
   \item{logconst}{logical; transform the parameters to use log constant 1.702?}
 }

\details{
Returns a data matrix simulated from the parameters.
}

\references{
Reckase, M. D. (2009). \emph{Multidimensional Item Response Theory}. New York: Springer.
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\examples{
\dontrun{
#parameters from Reckase (2009), p. 153
a <- matrix(c(
 .7471, .0250, .1428,
 .4595, .0097, .0692,
 .8613, .0067, .4040,
1.0141, .0080, .0470,
 .5521, .0204, .1482,
1.3547, .0064, .5362,
1.3761, .0861, .4676,
 .8525, .0383, .2574,
1.0113, .0055, .2024,
 .9212, .0119, .3044,
 .0026, .0119, .8036,
 .0008, .1905,1.1945,
 .0575, .0853, .7077,
 .0182, .3307,2.1414,
 .0256, .0478, .8551,
 .0246, .1496, .9348,
 .0262, .2872,1.3561,
 .0038, .2229, .8993,
 .0039, .4720, .7318,
 .0068, .0949, .6416,
 .3073, .9704, .0031,
 .1819, .4980, .0020,
 .4115,1.1136, .2008,
 .1536,1.7251, .0345,
 .1530, .6688, .0020,
 .2890,1.2419, .0220,
 .1341,1.4882, .0050,
 .0524, .4754, .0012,
 .2139, .4612, .0063,
 .1761,1.1200, .0870),30,3,byrow=TRUE)

d <- c(.1826,-.1924,-.4656,-.4336,-.4428,-.5845,-1.0403,
  .6431,.0122,.0912,.8082,-.1867,.4533,-1.8398,.4139,
  -.3004,-.1824,.5125,1.1342,.0230,.6172,-.1955,-.3668,
  -1.7590,-.2434,.4925,-.3410,.2896,.006,.0329)

mu <- c(-.4, -.7, .1)
sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)

dataset1 <- simdata(a, d, 2000)
dataset2 <- simdata(a, d, 2000, mu = mu, sigma = sigma)

   }
}

\keyword{data} 