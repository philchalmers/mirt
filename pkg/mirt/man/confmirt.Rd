\name{confmirt}

\alias{confmirt}
\alias{print.confmirt}
\alias{coef.confmirt}

\title{Confirmatory Full-Information Item Factor Analysis for Mixed Data Formats}
\description{
  \code{confmirt} fits a conditional maximum likelihood factor analysis model to dichotomous and polychotomous data under the item response theory paradigm using Cai's (2010) Metropolis-Hastings Robbins-Monro algorithm.
}

\usage{
polymirt <- function(data, nfact, guess = 0, prev.cor = NULL,
	ncycles = 2000, SEM.cycles = 50, kdraws = 1, tol = .0005, debug = FALSE, ...)

\method{print}{polymirt}(x, ...)

\method{summary}{polymirt}(object, rotate='varimax', digits = 3, ...)

\method{coef}{polymirt}(object, SE = TRUE, digits = 3, ...)

\method{plot}{polymirt}(x, npts = 50,
  rot = list(x = -70, y = 30, z = 10), ...)

}

\arguments{
  \item{data}{a \code{matrix} or \code{data.frame} that consists of numerically ordered data.}
  \item{nfact}{number of factors to be extracted.}
  \item{guess}{starting values for the pseudo-guessing parameter. Can be entered as a single value to assign a global guessing parameter or may be entered as a numeric vector for each item. When the value is greater than zero the parameter will be estimated.}
  \item{prev.cor}{uses a previously computed correlation matrix to be used to estimate starting values the estimation. The input could be any correlation matrix, but it is advised to use a matrix of polychoric correlations.}
  \item{rotate}{type of rotation to perform after the initial orthogonal parameters have been extracted. See below for list of possible rotations}
  \item{ncycles}{the maximum number of iterations to be performed}
  \item{SE}{logical; display the standard errors?}
  \item{x}{an object of class \code{polymirt} to be plotted or printed}
  \item{object}{a model estimated from \code{polymirt} of class \code{polymirt}}
  \item{digits}{the number of significant digits to be rounded}
  \item{npts}{number of quadrature points to be used for plotting features. Larger values make plots look smoother}
  \item{rot}{allows rotation of the 3D graphics}
  \item{debug}{logical; turn on debugging features?}
  \item{...}{additional arguments to be passed}
}

\details{

  \code{polymirt} follows the item factor analysis strategy by a stochastic version of maximum likelihood estimation described by Cai (2010). The general equation used for multidimensional item response theory in this package is in the logistic form with a scaling correction of 1.702. This correction is applied to allow comparison to mainstream programs such as TESTFACT (2003) and POLYFACT.

  \code{summary.mirt} allows for various rotations available from the \code{GPArotation} package. These are:

  \describe{
    \item{orthogonal: }{\code{"varimax", "quartimax", "tandemI", "tandemII", "entropy", "mccammon"}}
    \item{oblique: }{\code{"promax", "oblimin", "quartimin", "oblimax", "simplimax"}}
  }

  Using \code{plot} will plot the test information function for 1 and 2 dimensional solutions. To examine individuals item plots use \code{\link{itemplot}} (although the \code{\link[plink]{plink}} package is much more suitable for IRT graphics) which will also plot information and surface functions.
}

\value{
  \code{polymirt} returns an object of class \code{polymirt}, with the following elements:
  \item{pars}{estimated parameters of the model, the rightmost columns being the intercepts}
  \item{guess}{a vector of the estimated pseudo-guessing parameters}
  \item{SEpars}{a matrix of standard errors for the parameters estimated}
  \item{Theta}{last iteration of the Metropolis-Hastings imputation of the factor scores}
  \item{fulldata}{complete data complete with the original item names and values}
  \item{F}{unrotated factor loadings matrix}
  \item{h2}{a vector of communalities}
  \item{cycles}{number of MHRM cycles performed; will be less than \code{ncycles} if the tolerance is reached}
  \item{itemloc}{a vector indicating where the item location begin}
  \item{converge}{logical; did the model converge?}
  \item{Call}{the function call}
}

\references{

Cai, L. (2010). High-Dimensional exploratory item factor analysis by a Metropolis-Hastings Robbins-Monro algorithm. \emph{Psychometrika, 75}, 33-57.

Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., & Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics, and Full-information Item Factor Analysis [Computer software]. Lincolnwood, IL: Scientific Software International.

}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\seealso{\code{\link{expand.table}},\code{\link{key2binary}}}

\examples{
\dontrun{
#load LSAT section 7 data and compute 1 and 2 factor models
data(LSAT7)
fulldata <- expand.table(LSAT7)

(mod1 <- polymirt(fulldata, 1))
summary(mod1)
(mod2 <- polymirt(fulldata, 2))
summary(mod2)
coef(mod2, SE = TRUE)
plot(mod2)

###########
data(Science)
(mod1 <- polymirt(Science,1))
summary(mod1)
coef(mod1)

(mod2 <- polymirt(Science,2))
summary(mod2,'promax')
coef(mod2)
plot(mod2)

     }
}

\keyword{models}
