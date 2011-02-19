\name{mirt}

\alias{mirt}
\alias{print.mirt}
\alias{summary.mirt}
\alias{coef.mirt}
\alias{anova.mirt}
\alias{fitted.mirt}
\alias{residuals.mirt}
\alias{plot.mirt}

\title{Full-Information Item Factor Analysis (Multidimensional Item Response Theory)}
\description{
  \code{mirt} fits an unconditional maximum likelihood factor analysis model to dichotomous data under the item response theory paradigm. Pseudo-guessing parameters may be included but must be declared as constant, since the estimation of these parameters often leads to unacceptable solutions.
}

\usage{
mirt(fulldata, nfact, guess = 0, prev.cor = NULL, par.prior = NULL,
  startvalues = NULL, quadpts = NULL, ncycles = 50, EMtol = .005, nowarn = TRUE, ...)

\method{print}{mirt}(x, ...)

\method{summary}{mirt}(object, digits = 3, rotate='varimax', ...)

\method{coef}{mirt}(object, digits = 3, ...)

\method{anova}{mirt}(object, object2, ...)

\method{fitted}{mirt}(object, digits = 3, ...)

\method{plot}{mirt}(x, type = 'curve', npts = 30,
  rot = list(x = -70, y = 30, z = 10), ...)

\method{residuals}{mirt}(object, digits = 3, res.cor=FALSE, ...)

}

\arguments{
  \item{fulldata}{a \code{matrix} or \code{data.frame} that consists of only 0, 1, and \code{NA} values to be factor analyzed. If scores have been recorded by the response pattern then they can be recoded to dichotomous format using the \code{\link{key2binary}} function}
  \item{nfact}{number of factors to be extracted. Maximum of 10}
  \item{guess}{fixed pseudo-guessing parameter. Can be entered as a single value to assign a global guessing parameter or may be entered as a numeric vector for each item}
  \item{prev.cor}{uses a previously computed correlation matrix to be used to estimate starting values for the EM estimation. The input could be any correlation matrix, but it is advised to use a matrix of tetrachoric correlations (with Carroll's (1945) correction for guessing, if necessary) which can be computed by \code{\link{tetrachor}}. Otherwise, the value may be set to a previously computed \code{mirt} model to obtain the correlation matrix therein}
  \item{par.prior}{a list declaring which items should have assumed priors distributions, and what these prior weights are. Elements are \code{slope} and \code{int} to specify the coefficients beta prior for the slopes and normal prior for the intercepts, and \code{slope.items} and \code{int.items} to specify which items to constrain. The value in  \code{slope} is the \emph{p} meta-parameter for the beta distribution (where \emph{p} > 1 constrains the slopes), and the two values in \code{int} are the normal distribution intercept and variance. Larger values of the variance have less impact on the solution. For example, if items 2 and 3 were Heywood cases with no extreme item facilities, and item 4 had a very large item facility (say, greater than .95) then a possible constraint might be \code{par.prior = list(int = c(0,2), slope = 1.2, int.items = 4, slope.items = c(2,3))}}
  \item{rotate}{type of rotation to perform after the initial orthogonal parameters have been extracted. See below for list of possible rotations}
  \item{startvalues}{user declared start values for parameters}
  \item{quadpts}{number of quadrature points per dimension. If \code{NULL} then they are determined by \code{ceiling{15/nfact}} for \code{nfact} < 6, else they are set to 3}
  \item{ncycles}{the number of EM iterations to be performed}
  \item{EMtol}{if the largest change in the EM cycle is less than this value then the EM iteration are stopped early}
  \item{x}{an object of class \code{mirt} to be plotted or printed}
  \item{object}{a model estimated from \code{mirt} of class \code{mirt}}
  \item{object2}{a second model estimated from \code{mirt} of class \code{mirt} with more estimated parameters than \code{object}}
  \item{digits}{number of significant digits to be rounded}
  \item{type}{type of plot to view; can be \code{'curve'} for the total test score as a function of two dimensions, or \code{'info'} to show the test information function for two dimensions}
  \item{npts}{number of quadrature points to be used for plotting features. Larger values make plots look smoother}
  \item{rot}{allows rotation of the 3D graphics}
  \item{nowarn}{logical; suppress warnings from dependent packages?}
  \item{res.cor}{logical; print the residual correlation matrix instead of residual values?}
  \item{...}{additional arguments to be passed}
}

\details{

  \code{mirt} follows the item factor analysis strategy by marginal maximum likelihood estimation (MML) that is outlined in Bock and Aiken (1981) and Bock, Gibbons and Muraki (1988). Nested models may be compared via the approximate chi-squared difference test or by a reduction in AIC (accessible via \code{\link{anova.mirt}}). The general equation used for multidimensional item response theory in this package is in the logistic form with a scaling correction of 1.702. This correction is applied to allow comparison to mainstream programs such as TESTFACT (2003). The target equation is

\deqn{P(X | \theta; \bold{a}_i; d_i; g_i) = g_i + (1 - g_i) * 1.702(\bold{a}_i' \bold{a}_i + d_i)/ (1 + 1.702(\bold{a}_i' \bold{a}_i + d_i))}

  where \emph{i} is the item index, \eqn{\bold{a}_i} is the vector of discrimination parameters (i.e. slopes), \eqn{d_i} is the intercept, and \eqn{g_i} is the pseudo-guessing parameter. To avoid estimation difficulties the \eqn{g_i}'s must be specified by the user.

  Estimation of the model begins with computing a matrix of tetrachoric correlations, potentially with Carroll's (1945) adjustment for chance responds. A MINRES factor analysis with \code{nfact} is then extracted and item parameters are estimated by \eqn{a_{ij} = f_{ij}/u_i}, where  \eqn{f_{ij}} is the factor loading for the \emph{i}th item on the \emph{j}th factor, and \eqn{u_i} is the square root of the factor uniqueness, \eqn{\sqrt{1 - h_i^2}}. The initially intercept parameters are determined by calculating the inverse normal of the item facility (i.e., item easiness), \eqn{q_i}, to obtain \eqn{d_i = q_i / u_i}. Following these initial estimates the model is iterated using the EM estimation strategy with non-adaptive quadrature points. Implicit equation accelerations (described in Ramsey (1975)) are also added to facilitate parameter convergence speed, and this is adjusted every third cycle.

  Factor scores are estimated assuming a normal prior distribution and can be appended to the input data matrix (\code{full.data = TRUE}) or displayed in a summary table for all the unique response patterns. \code{summary.mirt} allows for various rotations available from the \code{GPArotation} package. These are:

  \describe{
    \item{orthogonal: }{\code{"varimax", "quartimax", "tandemI", "tandemII", "entropy", "mccammon"}}
    \item{oblique: }{\code{"promax", "oblimin", "quartimin", "oblimax", "simplimax"}}
  }

  Using \code{plot} will plot the either the test surface function or the test information function for 1 and 2 dimensional solutions. To examine individuals item plots use \code{\link{itemplot}} (although the \code{\link[plink]{plink}} package is much more suitable for IRT graphics) which will also plot information and surface functions.
}

\value{
  \code{mirt} returns an object of class \code{mirt}, with the following elements:
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
  \item{Theta}{combination of quadrature points for each factor; contains \eqn{quadpts^nfact} number of rows}
  \item{fulldata}{complete data complete with the original item names}
  \item{empdist}{calculated empirical distribution}
  \item{F}{unrotated factor loadings matrix}
  \item{h2}{a vector of communalities}
  \item{EMiter}{number of EM cycles performed; will be less than \code{ncycles} if the tolerance is reached}
  \item{par.prior}{list of the assumed prior parameters used to constrain the solution. If no priors are specified then the solution will be a matrix of \code{c(1,0,0)}, indicating no influence on the parameters. For instance, if we wish to contrain the slopes of item 1 and 3 to 1.2, and the intercepts of item 4 with a weak normal prior (mean = 0) then \code{par.prior = list(slope = 1.2, slope.items = c(1,3), int = c(0,3), int.items = 4)}}
  \item{Call}{the function call}
}

\section{Convergence}{

Full information factor analysis is known to have problems with convergence, and some items may need to be constrained or removed entirely to allow for and acceptable solution. Be mindful of the item facility values that are printed with \code{coef.mirt} since these will be helpful in determining whether a guessing parameter should be removed (item facility value is too close to the guessing parameter) or if an item should be constrained or removed entirely (values too close to 0 or 1). As a general rule, items with facilities greater than .95, or items that are only .05 greater than the guessing parameter, should be considered for being dropped from the analysis. It is also important to pay attention to the log-likelihood values printed after the model is estimated, since if these values do not increase then there likely is a problem with the estimation. Also, increasing the number of quadrature points per dimension should help to stabilize the estimation process.

}

\references{

Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm. \emph{Psychometrika, 46}(4), 443-459.

Bock, R. D., Gibbons, R., & Muraki, E. (1988). Full-Information Item Factor Analysis. \emph{Applied Psychological Measurement, 12}(3), 261-280.

Carroll, J. B. (1945). The effect of difficulty and chance success on correlations between items and between tests. \emph{Psychometrika, 26}, 347-372.

Ramsay, J. O. (1975). Solving implicit equations in psychometric data analysis. \emph{Psychometrika, 40}(3), 337-360.

Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., & Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics, and Full-information Item Factor Analysis [Computer software]. Lincolnwood, IL: Scientific Software International.

}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\seealso{\code{\link{expand.table}},\code{\link{key2binary}}}

\examples{
\dontrun{
#load LSAT section 7 data and compute 1 and 2 factor models
data(LSAT7)
fulldata <- expand.table(LSAT7)

mod1 <- mirt(fulldata, 1)
mod1
summary(mod1)
mod2 <- mirt(fulldata, 2, prev.cor = mod1) #use previous tetrachoric matrix to save time
mod2
summary(mod2)
anova(mod1, mod2) #compare the two models

coef(mod2)
scores <- fscores(mod2) #save factor score table

###########
data(SAT12)
fulldata <- key2binary(SAT12,
  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
tmat <- tetrachor(fulldata, .1) #this can save some time if reused
mod1 <- mirt(fulldata, 1, guess = .1, prev.cor = tmat)
coef(mod1)
#didn't converge, item 6 and 32 appear too difficult (facility = .16)
#fix guessing pars
guess <- rep(.1,32)
guess[c(6,32)] <- 0

tmat <- tetrachor(fulldata, guess)
mod1 <- mirt(fulldata, 1, guess = guess, prev.cor = tmat)
mod2 <- mirt(fulldata, 2, guess = guess, prev.cor = tmat, ncycles = 100)
anova(mod1,mod2)

     }
}

\keyword{models}
