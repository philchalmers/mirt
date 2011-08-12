\name{mirt}
\alias{mirt}
\alias{summary,mirt-method}
\alias{coef,mirt-method}
\alias{anova,mirt-method}
\alias{fitted,mirt-method}
\alias{plot,mirt-method}
\alias{residuals,mirt-method}

\title{Full-Information Item Factor Analysis (Multidimensional Item Response Theory)}
\description{
  \code{mirt} fits an unconditional maximum likelihood factor analysis model to dichotomous data under the item response theory paradigm. Pseudo-guessing parameters may be included but must be declared as constant, since the estimation of these parameters often leads to unacceptable solutions. Missing values are automatically assumed to represent the lowest response category.
}

\usage{
mirt(fulldata, nfact, guess = 0, prev.cor = NULL, par.prior = FALSE,
  startvalues = NULL, quadpts = NULL, ncycles = 300, tol = .001, nowarn = TRUE, 
  debug = FALSE, ...)

\S4method{summary}{mirt}(object, rotate='varimax', suppress = 0, digits = 3, ...)

\S4method{coef}{mirt}(object, digits = 3, ...)

\S4method{anova}{mirt}(object, object2, ...)

\S4method{fitted}{mirt}(object, digits = 3, ...)

\S4method{plot}{mirt}(x, type = 'info', npts = 50, rot = list(x = -70, y = 30, z = 10), ...)

\S4method{residuals}{mirt}(object, restype = 'LD', digits = 3, ...)

}

\arguments{
  \item{fulldata}{a \code{matrix} or \code{data.frame} that consists of only 0, 1, and \code{NA} values to be factor analyzed. If scores have been recorded by the response pattern then they can be recoded to dichotomous format using the \code{\link{key2binary}} function}
  \item{nfact}{number of factors to be extracted}
  \item{guess}{fixed pseudo-guessing parameters. Can be entered as a single value to assign a global guessing parameter or may be entered as a numeric vector for each item}
  \item{prev.cor}{use a previously computed correlation matrix to be used to estimate starting values for the EM estimation? Default in \code{NULL}}
  \item{par.prior}{a list declaring which items should have assumed priors distributions, and what these prior weights are. Elements are \code{slope} and \code{int} to specify the coefficients beta prior for the slopes and normal prior for the intercepts, and \code{slope.items} and \code{int.items} to specify which items to constrain. The value in  \code{slope} is the \emph{p} meta-parameter for the beta distribution (where \emph{p} > 1 constrains the slopes), and the two values in \code{int} are the normal distribution intercept and variance. Larger values of the variance have less impact on the solution. For example, if items 2 and 3 were Heywood cases with no extreme item facilities, and item 4 had a very large item facility (say, greater than .95) then a possible constraint might be \code{par.prior = list(int = c(0,2), slope = 1.2, int.items = 4, slope.items = c(2,3))}}
  \item{rotate}{type of rotation to perform after the initial orthogonal parameters have been extracted. See below for list of possible rotations}
  \item{startvalues}{user declared start values for parameters}
  \item{quadpts}{number of quadrature points per dimension}
  \item{ncycles}{the number of EM iterations to be performed}
  \item{tol}{if the largest change in the EM cycle is less than this value then the EM iteration are stopped early}
  \item{x}{an object of class \code{mirt} to be plotted or printed}
  \item{object}{a model estimated from \code{mirt} of class \code{mirt}}
  \item{object2}{a second model estimated from \code{mirt} of class \code{mirt} with more estimated parameters than \code{object}}
  \item{suppress}{a numeric value indicating which (possibly rotated) factor loadings should be suppressed. Typical values are around .3 in most statistical software}
  \item{digits}{number of significant digits to be rounded}
  \item{type}{type of plot to view; can be \code{'curve'} for the total test score as a function of two dimensions, or \code{'info'} to show the test information function for two dimensions}
  \item{npts}{number of quadrature points to be used for plotting features. Larger values make plots look smoother}
  \item{rot}{allows rotation of the 3D graphics}
  \item{restype}{type of residuals to be displayed. Can be either \code{'LD'} for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the expected values for the frequencies of every response pattern}
  \item{nowarn}{logical; suppress warnings from dependent packages?}
  \item{debug}{logical; turn on debugging features?}
  \item{...}{additional arguments to be passed}
}

\details{

  \code{mirt} follows the item factor analysis strategy by marginal maximum likelihood estimation (MML) outlined in Bock and Aiken (1981) and Bock, Gibbons and Muraki (1988). Nested models may be compared via the approximate chi-squared difference test or by a reduction in AIC (comparison via \code{\link{anova}}). The general equation used for multidimensional item response theory in this package is in the logistic form with a scaling correction of 1.702. This correction is applied to allow comparison to mainstream programs such as TESTFACT (2003). The target equation is

\deqn{P(X | \theta; \bold{a}_i; d_i; g_i) = g_j + (1 - g_j) * exp(1.702(\bold{a}_j' \bold{a}_j + d_j))/ (1 + exp(1.702(\bold{a}_j' \bold{a}_j + d_j)))}

  where \emph{j} is the item index, \eqn{\bold{a}_j} is the vector of discrimination parameters (i.e. slopes), \eqn{d_j} is the intercept, and \eqn{g_j} is the pseudo-guessing parameter. To avoid estimation difficulties the \eqn{g_j}'s must be specified by the user.

  Estimation begins by computing a matrix of quasi-tetrachoric correlations, potentially with Carroll's (1945) adjustment for chance responds. A MINRES factor analysis with \code{nfact} is then extracted and item parameters are estimated by \eqn{a_{ij} = f_{ij}/u_j}, where  \eqn{f_{ij}} is the factor loading for the \emph{j}th item on the \emph{i}th factor, and \eqn{u_j} is the square root of the factor uniqueness, \eqn{\sqrt{1 - h_j^2}}. The initially intercept parameters are determined by calculating the inverse normal of the item facility (i.e., item easiness), \eqn{q_j}, to obtain \eqn{d_j = q_j / u_j}. Following these initial estimates the model is iterated using the EM estimation strategy with fixed quadrature points. Implicit equation accelerations (described in Ramsey (1975)) are also added to facilitate parameter convergence speed, and this is adjusted every third cycle.

  Factor scores are estimated assuming a normal prior distribution and can be appended to the input data matrix (\code{full.data = TRUE}) or displayed in a summary table for all the unique response patterns. \code{summary} allows for various rotations available from the \code{GPArotation} package. These are:

  \describe{
    \item{orthogonal: }{\code{"varimax", "quartimax", "tandemI", "tandemII", "entropy", "mccammon"}}
    \item{oblique: }{\code{"promax", "oblimin", "quartimin", "oblimax", "simplimax"}}
  }

  Using \code{plot} will plot the either the test surface function or the test information function for 1 and 2 dimensional solutions. To examine individuals item plots use \code{\link{itemplot}} (although the \code{\link[plink]{plink}} package is much more suitable for IRT graphics) which will also plot information and surface functions. Residuals are computed using the LD statistic (Chen \& Thissen, 1997) in the lower diagonal of the matrix returned by \code{residuals}, and Cramer's V above the diagonal.
}

\section{Convergence}{

Unrestricted full-information factor analysis is known to have problems with convergence, and some items may need to be constrained or removed entirely to allow for an acceptable solution. Be mindful of the item facility values that are printed with \code{coef} since these will be helpful in determining whether a guessing parameter should be removed (item facility value is too close to the guessing parameter) or if an item should be constrained or removed entirely (values too close to 0 or 1). As a general rule, items with facilities greater than .95, or items that are only .05 greater than the guessing parameter, should be considered for removal from the analysis. Also, increasing the number of quadrature points per dimension may help to stabilize the estimation process.

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

(mod1 <- mirt(fulldata, 1))
summary(mod1)
residuals(mod1)
plot(mod1) #test information function

(mod2 <- mirt(fulldata, 2))
summary(mod2)
coef(mod2)
residuals(mod2)

anova(mod1, mod2) #compare the two models
scores <- fscores(mod2) #save factor score table

###########
data(SAT12)
fulldata <- key2binary(SAT12,
  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))

#without guessing scree(tmat) #looks like a 2 factor solution
mod1 <- mirt(fulldata, 1)
mod2 <- mirt(fulldata, 2)
mod3 <- mirt(fulldata, 3)
anova(mod1,mod2)
anova(mod2, mod3) #negative AIC, 2 factors probably best

#with guessing
mod1g <- mirt(fulldata, 1, guess = .1)
coef(mod1g)
mod2g <- mirt(fulldata, 2, guess = .1)
coef(mod2g)
anova(mod1g, mod2g)
summary(mod2g, rotate='promax')
     }
}

\keyword{models}
