\name{bfactor}
\alias{bfactor}
\alias{print.bfactor}
\alias{summary.bfactor}
\alias{coef.bfactor}
\alias{anova.bfactor}
\alias{fitted.bfactor}
\alias{residuals.bfactor}

\title{Full-Information Item Bifactor Analysis}
\description{
  \code{bfactor} fits a confirmatory maximum likelihood bifactor model to dichotomous data under the item response theory paradigm. Pseudo-guessing parameters may be included but must be declared as constant, since the estimation of these parameters often leads to unacceptable solutions.
}

\usage{
bfactor(fulldata, specific, guess = 0, prev.cor=NULL, par.prior = NULL,
  startvalues = NULL,  quadpts = NULL, ncycles = 50, EMtol=.005, nowarn = TRUE, ...)

\method{print}{bfactor}(x, ...)

\method{summary}{bfactor}(object, digits = 3, ...)

\method{coef}{bfactor}(object, digits = 3, ...)

\method{anova}{bfactor}(object, object2, ...)

\method{fitted}{bfactor}(object, digits = 3, ...)

\method{residuals}{bfactor}(object, digits = 3, res.cor=FALSE, ...)

}

\arguments{
  \item{fulldata}{a \code{matrix} or \code{data.frame} that consists of only 0, 1, and \code{NA} values to be factor analyzed. If scores have been recorded by the response pattern then they can be recoded to dichotomous format using the \code{\link{key2binary}} function. }
  \item{specific}{a numeric vector specifying where the which factor loads on which item. For example if for a 4 item test with two specific factors, the first specific factor loading on the first two items, then the vector may be specified as \code{c(1,1,2,2)}.}
  \item{guess}{fixed pseudo-guessing parameter. Can be entered as a single value to assign a global guessing parameter or may be entered as a numeric vector for each item.}
  \item{prev.cor}{uses a previously computed correlation matrix to be used to estimate starting values for the EM estimation. The input could be any correlation matrix, but it is advised to use a matrix of tetrachoric correlations (with Carroll's (1945) correction for guessing, if necessary) which can be computed by \code{\link{tetrachor}}. Otherwise, the value may be set to a previously computed \code{bfactor} model to obtain the correlation matrix therein.}
  \item{par.prior}{a list declaring which items should have assumed priors distributions, and what these prior weights are. Elements are \code{slope} and \code{int} to specify the coefficients beta prior for the slopes and normal prior for the intercepts, and \code{slope.items} and \code{int.items} to specify which items to constrain. The value in  \code{slope} is the \emph{p} meta-parameter for the beta distribution (where \emph{p} > 1 constrains the slopes), and the two values in \code{int} are the normal distribution intercept and variance. Larger values of the variance have less impact on the solution. For example, if items 2 and 3 were Heywood cases with no extreme item facilities, and item 4 had a very large item facility (say, greater than .95) then a possible constraint might be \code{par.prior = list(int = c(0,2), slope = 1.2, int.items = 4, slope.items = c(2,3))}}
  \item{startvalues}{user declared start values for parameters. }
  \item{quadpts}{number of quadrature points per dimension. If \code{NULL} then the number of quadrature points is set to 9. }
  \item{ncycles}{the number of EM iterations to be performed.}
  \item{EMtol}{if the largest change in the EM cycle is less than this value then the EM iteration are stopped early.}
  \item{x}{a model of class \code{bfactor} to be printed}
  \item{object}{a model estimated from \code{bfactor} of class \code{bfactor}.}
  \item{object2}{a model estimated from \code{mirt} of class \code{mirt} with more estimated parameters than \code{object}.}
  \item{digits}{number of significant digits to be rounded.}
  \item{nowarn}{logical; suppress warnings from dependent packages?}
  \item{res.cor}{logical; print the residual correlation matrix instead of residual values?}
  \item{...}{additional arguments to be passed.}
}

\details{

  \code{bfactor} follows the item factor analysis strategy explicated by Gibbons and Hedeker (1992). Nested models may be compared via the approximate chi-squared difference test or by a reduction in AIC (accessible via \code{\link{anova.bfactor}}); note that this only makes sense when comparing class \code{bfactor} models to class \code{mirt}. The general equation used for item bifactor analysis in this package is in the logistic form with a scaling correction of 1.702. This correction is applied to allow comparison to mainstream programs such as TESTFACT 4 (2003).

  Unlike TESTFACT 4 (2003) initial start values are computed by using information from the matrix of tetrachoric correlations, potentially with Carroll's (1945) adjustment for chance responses. To begin, a MINRES factor analysis with one factor is extracted, and the transformed loadings and intercepts (see \link{mirt} for more details) are used as starting values for the general factor loadings and item intercepts. Values for the specific factor loadings are taken to be half the magnitude of the extracted general factor loadings. Note that while the sign of the loading may be incorrect for specific factors (and possibly for some of the general factor loadings) the intercepts and general factor loadings will be relatively close to the final solution. These initial values should be an improvement over the TESTFACT 4 initial starting values of 1.414 for all the general factor slopes, 1 for all the specific factor slopes, and 0 for all the intercepts.

  Factor scores are estimated assuming a normal prior distribution and can be appended to the input data matrix (\code{full.data = TRUE}) or displayed in a summary table for all the unique response patterns. Fitted and residual values may be observed via \code{fitted.bfactor} and \code{residual.bfactor}. To examine individuals item plots use \code{\link{itemplot.bfactor}} which will also plot information and surface functions (although the \code{\link[plink]{plink}} package is much more suitable for IRT graphics).
}

\value{
  \code{bfactor} returns an object of class \code{bfactor}, with the following elements:
  \item{pars}{estimated parameters of the model, the rightmost column being the intercept}
  \item{guess}{a vector of the constant user-supplied pseudo-guessing parameters}
  \item{facility}{a proportion vector of item endorsement}
  \item{cormat}{if not initially specified, a tetrachoric matrix with Carroll's correction for guessing (if applicable); else the user supplied correlation matrix }
  \item{log.lik}{log-likelihood of the model}
  \item{AIC}{Akaike Information Criteria}
  \item{X2}{chi-squared statistic}
  \item{df}{degrees of freedom associated with \code{X2}}
  \item{p}{probability associated with \code{df} and \code{X2}}
  \item{tabdata}{a summary table of unique response patterns with the frequency of occurrence in the rightmost column}
  \item{Pl}{calculated expected probability vector for each unique response pattern}
  \item{Theta}{combination of quadrature points for each factor; contains \eqn{quadpts^{nfact}} number of rows}
  \item{fulldata}{complete data complete with the original item names}
  \item{empdist}{calculated empirical distribution}
  \item{h2}{a vector of communalities}
  \item{EMiter}{number of EM cycles performed; will be less than \code{ncycles} if the tolerance is reached}
  \item{specific}{a numeric vector specifying where the specific item loadings are located}
  \item{logicalfact}{a logical matrix indicating where the factor loadings are located}
  \item{itemnames}{the names of the items}
  \item{Call}{the function call}

}

\references{

Gibbons, R. D., & Hedeker, D. R. (1992). Full-information Item Bi-Factor Analysis. \emph{Psychometrika, 57}, 423-436.

Carroll, J. B. (1945). The effect of difficulty and chance success on correlations between items and between tests. \emph{Psychometrika, 26}, 347-372.

Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., & Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics, and Full-information Item Factor Analysis [Computer software]. Lincolnwood, IL: Scientific Software International.

}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\seealso{\code{\link{expand.table}},\code{\link{key2binary}}}

\examples{
\dontrun{

###load SAT12 and compute bifactor model with 3 specific factors
data(SAT12)
fulldata <- key2binary(SAT12,
  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
tmat <- tetrachor(fulldata) #this can take a while
specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
mod1 <- bfactor(fulldata, specific,prev.cor = tmat)

###Try with guessing parameters added
guess <- rep(.1,32)
mod2 <- bfactor(fulldata, specific, guess = guess, prev.cor = tmat)
coef(mod2) #didn't converge, item 6 item facility too close to guessing parameter
guess[6]<-0
mod2 <- bfactor(fulldata, specific, guess = guess, prev.cor = tmat)
anova(mod1,mod2) #although no change in df X2 decreased

    }
}

\keyword{models}
