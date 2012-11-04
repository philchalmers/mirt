#' Full-Information Item Factor Analysis (Multidimensional Item Response
#' Theory)
#' 
#' \code{mirt} fits an unconditional maximum likelihood factor analysis model
#' to dichotomous and polytomous data under the item response theory paradigm. 
#' Fits univariate and multivariate Rasch, 1-4PL, graded, (generalized) partial credit, 
#' nominal, multiple choice, and partially compenatory models using the EM algorithm.
#' 
#' \code{mirt} follows the item factor analysis strategy by marginal maximum
#' likelihood estimation (MML) outlined in Bock and Aiken (1981), Bock,
#' Gibbons and Muraki (1988), and Muraki and Carlson (1995). 
#' Nested models may be compared via the approximate
#' chi-squared difference test or by a reduction in AIC/BIC values (comparison
#' via \code{\link{anova}}). The general equation used for 
#' multidimensional item response theory is a logistic form with a scaling
#' correction of 1.702. This correction is applied to allow comparison to
#' mainstream programs such as TESTFACT (2003) and POLYFACT. 
#' 
#' Factor scores are estimated assuming a normal prior distribution and can be
#' appended to the input data matrix (\code{full.data = TRUE}) or displayed in
#' a summary table for all the unique response patterns. \code{summary} and \code{coef} allow
#' for all the rotations available from the \code{GPArotation} package (e.g., \code{rotate = 'oblimin'})
#' as well as a \code{'promax'} rotation. 
#' 
#' Using \code{plot} will plot the either the test surface function or the test
#' information function for 1 and 2 dimensional solutions. To examine
#' individual item plots use \code{\link{itemplot}}. Residuals are
#' computed using the LD statistic (Chen & Thissen, 1997) in the lower
#' diagonal of the matrix returned by \code{residuals}, and Cramer's V above
#' the diagonal.
#' 
#' @section Confirmatory IRT:
#' 
#' Specification of the confirmatory item factor analysis model follows many of
#' the rules in the SEM framework for confirmatory factor analysis. The
#' variances of the latent factors are automatically fixed to 1 to help
#' facilitate model identification. All parameters may be fixed to constant
#' values or set equal to other parameters using the appropriate declarations.
#' If the model is confirmatory then the returned class will be 'ConfirmatoryClass'.
#' 
#' @section Exploratory IRT:
#' 
#' Specifying a number as the second input to confmirt an exploratory IRT model is estimated and 
#' can be viewed as a stochastic analogue of \code{mirt}, with much of the same behaviour and 
#' specifications. Rotation and target matrix options will be used in this subroutine and will be
#' passed to the returned object for use in generic functions such as \code{summary()} and 
#' \code{fscores}. Again, factor means and variances are fixed to ensure proper identification. 
#' If the model is confirmatory then the returned class will be 'ExploratoryClass'.
#' 
#' Estimation often begins by computing a matrix of quasi-tetrachoric correlations,
#' potentially with Carroll's (1945) adjustment for chance responds. A MINRES
#' factor analysis with \code{nfact} is then extracted and item parameters are
#' estimated by \eqn{a_{ij} = f_{ij}/u_j}, where \eqn{f_{ij}} is the factor
#' loading for the \emph{j}th item on the \emph{i}th factor, and \eqn{u_j} is
#' the square root of the factor uniqueness, \eqn{\sqrt{1 - h_j^2}}. The
#' initial intercept parameters are determined by calculating the inverse
#' normal of the item facility (i.e., item easiness), \eqn{q_j}, to obtain
#' \eqn{d_j = q_j / u_j}. A similar implementation is also used for obtaining 
#' initial values for polychotomous items. Following these initial estimates the model is
#' iterated using the EM estimation strategy with fixed quadrature points.
#' Implicit equation accelerations described by Ramsey (1975) are also added to
#' facilitate parameter convergence speed, and these are adjusted every third
#' cycle.
#' 
#' @aliases mirt summary,ExploratoryClass-method coef,ExploratoryClass-method anova,ExploratoryClass-method 
#' fitted,ExploratoryClass-method plot,ExploratoryClass-method residuals,ExploratoryClass-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data, with missing data coded as \code{NA}
#' @param model an object returned from \code{confmirt.model()} declaring how
#' the factor model is to be estimated, or a single numeric value indicating the number 
#' of exploratory factors to estimate. See \code{\link{confmirt.model}} for
#' more details
#' @param itemtype type of items to be modeled, declared as a vector for each item or a single value
#' which will be repeated globally. The NULL default assumes that the items follow a graded or 2PL structure,
#' however they may be changed to the following: 'Rasch', '1PL', '2PL', '3PL', '3PLu', 
#' '4PL', 'graded', 'grsm', 'gpcm', 'nominal', 'mcm', 'PC2PL', and 'PC3PL', for the Rasch/partial credit, 1 and 2 parameter logistic, 
#' 3 parameter logistic (lower asymptote and upper), 4 parameter logistic, graded response model, 
#' rating scale graded response model, generalized partial credit model, nominal model, 
#' multiple choice model, and 2-3PL partially compensatory model, respectively 
#' @param grsm.block an optional numeric vector indicating where the blocking should occur when using 
#' the grsm, NA represents items that do not belong to the grsm block (other items that may be estimated
#' in the test data). For example, to specify two blocks of 3 with a 2PL item for the last item:
#' \code{grsm.block = c(rep(1,3), rep(2,3), NA)}. If NULL the all items are assumed to be within the same 
#' group and therefore have the same number of item categories
#' @param SE logical, estimate the standard errors? Calls the MHRM subroutine for a stochastic approximation
#' @param SEtol tollerance value used to stop the MHRM estimation when \code{SE = TRUE}. Lower values
#' will take longer but may be more stable for computing the information matrix
#' @param guess fixed pseudo-guessing parameters. Can be entered as a single
#' value to assign a global guessing parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param upper fixed upper bound parameters for 4-PL model. Can be entered as a single
#' value to assign a global guessing parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param prev.cor use a previously computed correlation matrix to be used to
#' estimate starting values for the EM estimation? Default in \code{NULL} 
#' @param rotate type of rotation to perform after the initial orthogonal
#' parameters have been extracted by using \code{summary}; default is \code{'varimax'}. 
#' See below for list of possible rotations. If \code{rotate != ''} in the \code{summary} 
#' input then the default from the object is ignored and the new rotation from the list 
#' is used instead
#' @param D a numeric value used to adjust the logistic metric to be more similar to a normal
#' cumulative density curve. Default is 1.702
#' @param Target a dummy variable matrix indicting a target rotation pattern
#' @param constrain a list of user declared equality constraints. To see how to define the
#' parameters correctly use \code{pars = 'values'} initially to see how the parameters are labeled.
#' To constrain parameters to be equal create a list with separate concatenated vectors signifying which
#' parameters to constrain. For example, to set parameters 1 and 5 equal, and also set parameters 2, 6, and 10 equal
#' use \code{constrain = list(c(1,5), c(2,6,10))}
#' @param parprior a list of user declared prior item probabilities. To see how to define the
#' parameters correctly use \code{pars = 'values'} initially to see how the parameters are labeled.
#' Can define either normal (normally for slopes and intercepts) or beta (for guessing and upper bounds) prior
#' probabilities. To specify a prior the form is c('priortype', ...), where normal priors 
#' are \code{parprior = list(c(parnumbers, 'norm', mean, sd))} and betas are 
#' \code{parprior = list(c(parnumbers, 'beta', alpha, beta))} 
#' @param pars a data.frame with the structure of how the starting values, parameter numbers, and estimation
#' logical values are defined. The user may observe how the model defines the values by using \code{pars = 
#' 'values'}, and this object can in turn be modified and input back into the estimation with \code{pars = 
#' mymodifiedpars}
#' @param quadpts number of quadrature points per dimension
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#' option. Only prints patterns that have standardized residuals greater than 
#' \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param x an object of class \code{mirt} to be plotted or printed
#' @param y an unused variable to be ignored
#' @param object a model estimated from \code{mirt} of class \code{ExploratoryClass} or 
#' \code{ConfirmatoryClass}
#' @param object2 a second model estimated from any of the mirt package estimation methods
#' \code{ExploratoryClass} with more estimated parameters than \code{object}
#' @param suppress a numeric value indicating which (possibly rotated) factor
#' loadings should be suppressed. Typical values are around .3 in most
#' statistical software. Default is 0 for no suppression
#' @param digits number of significant digits to be rounded
#' @param type type of plot to view; can be \code{'curve'} for the total test
#' score as a function of two dimensions, or \code{'info'} to show the test
#' information function for two dimensions
#' @param theta_angle numeric values ranging from 0 to 90 used in \code{plot}. If a vector is 
#' used then a bubble plot is created with the summed information across the angles specified 
#' (e.g., \code{theta_angle = seq(0, 90, by=10)})
#' @param npts number of quadrature points to be used for plotting features.
#' Larger values make plots look smoother
#' @param rot allows rotation of the 3D graphics
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param df.p logical; print the degrees of freedom and p-values?
#' @param verbose logical; print observed log-likelihood value at each iteration?
#' @param debug logical; turn on debugging features?
#' @param technical a list containing lower level technical parameters for estimation. May be:
#' \describe{ 
#' \item{MAXQUAD}{maximum number of quadratures; default 10000}
#' \item{MSTEPMAXIT}{number of M-step iterations; default 25}
#' \item{TOL}{EM convergence threshold; default .001}
#' \item{NCYCLES}{maximum number of EM cycles; default 300}
#' }
#' @param ... additional arguments to be passed
#' @section Convergence:
#' 
#' Unrestricted full-information factor analysis is known to have problems with
#' convergence, and some items may need to be constrained or removed entirely
#' to allow for an acceptable solution. As a general rule dichotomous items with
#' means greater than .95, or items that are only .05 greater than the
#' guessing parameter, should be considered for removal from the analysis or
#' treated with prior distributions. The same type of reasoning is applicable when including 
#' upper bound parameters as well. Also, increasing the number of quadrature
#' points per dimension may help to stabilize the estimation process.
#' 
#' @section IRT Models:
#' 
#' The parameter labels use follow the convention, here using two factors and \eqn{k} as the number 
#' of categories. Throughout all models D is a constant (default 1.702):
#' 
#' \describe{ 
#' \item{Rasch}{
#' Only one intercept estimated. \deqn{P(x = 1|\theta, d) = \frac{1}{1 + 
#' exp(-D*(\theta + d))}}
#' }
#' \item{1-4PL}{
#' Depending on the model \eqn{u} may be equal to 1 and \eqn{g} may be equal to 0. 
#' \deqn{P(x = 1|\theta, \psi) = g + \frac{(u - g)}{1 + exp(-D * 
#' (a_1 * \theta_1 + a_2 * \theta_2 + d))}} 
#' }
#' \item{graded}{
#' The graded model consists of sequential 2PL models, and here \eqn{k} is 
#' the predicted category. 
#' \deqn{P(x = k | \theta, \psi) = P(x \ge k | \theta, \phi) - P(x \ge k + 1 | \theta, \phi)}
#' }
#' \item{grsm}{
#' A more constrained version of the graded model where graded spacing is equal accross item blocks
#' and only adjusted by a single 'difficulty' parameter (c). Again,
#' \deqn{P(x = k | \theta, \psi) = P(x \ge k | \theta, \phi) - P(x \ge k + 1 | \theta, \phi)} 
#' but now 
#' \deqn{P = \frac{1}{1 + exp(-D * (a_1 * \theta_1 + a_2 * \theta_2 + d_k + c))}} 
#' } 
#' \item{gpcm/nominal}{For the gpcm the \eqn{d_k} values are treated as fixed and orderd values 
#' from 0:(k-1) (in the nominal model \eqn{d_0} is also set to 0). Additionally, for identification 
#' in the nominal model \eqn{ak_0 = 1}, \eqn{ak_k = (k - 1)}.
#' \deqn{P(x = k | \theta, \psi) = \frac{exp(-D * ak_k * (a_1 * \theta_1 + a_2 * \theta_2) + d_k)}
#' {\sum_i^k exp(-D * ak_k * (a_1 * \theta_1 + a_2 * \theta_2) + d_k)}}
#' }
#' \item{mcm}{For identification \eqn{ak_0 = d_0 = 0} and \eqn{\sum_0^k t_k = 1}.
#' \deqn{P(x = k | \theta, \psi) = C_0 (\theta) * t_k  + (1 - C_0 (\theta)) * 
#' \frac{exp(-D * ak_k * (a_1 * \theta_1 + a_2 * \theta_2) + d_k)}  
#' {\sum_i^k exp(-D * ak_k * (a_1 * \theta_1 + a_2 * \theta_2) + d_k)}}
#'
#' where \eqn{C_0 (\theta) = \frac{exp(-D * ak_0 * (a_1 * \theta_1 + a_2 * \theta_2) + d_0)}  
#' {\sum_i^k exp(-D * ak_k * (a_1 * \theta_1 + a_2 * \theta_2) + d_k)}}
#' }
#' \item{partcomp}{Partially compensatory models consist of the products of 2PL probability curves. 
#' \deqn{P(x = 1 | \theta, \psi) = g + (1 - g) (\frac{1}{1 + exp(-D * (a_1 * \theta_1 + d_1))} * 
#' \frac{1}{1 + exp(-D * (a_2 * \theta_2 + d_2))})}
#' }
#' }
#' 
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, 
#' \code{\link{confmirt}}, \code{\link{bfactor}}, \code{\link{multipleGroup}}, \code{\link{wald}}
#' \code{\link{itemplot}}, \code{\link{fscores}}
#' 
#' @references
#' 
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of
#' item parameters: Application of an EM algorithm. \emph{Psychometrika,
#' 46}(4), 443-459.
#' 
#' Bock, R. D., Gibbons, R., & Muraki, E. (1988). Full-Information Item Factor
#' Analysis. \emph{Applied Psychological Measurement, 12}(3), 261-280.
#' 
#' Carroll, J. B. (1945). The effect of difficulty and chance success on
#' correlations between items and between tests. \emph{Psychometrika, 26},
#' 347-372.
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6),
#' 1-29.
#'
#' Muraki, E. & Carlson, E. B. (1995). Full-information factor analysis for polytomous 
#' item responses. \emph{Applied Psychological Measurement, 19}, 73-90.
#' 
#' Ramsay, J. O. (1975). Solving implicit equations in psychometric data
#' analysis. \emph{Psychometrika, 40}(3), 337-360.
#' 
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., &
#' Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics,
#' and Full-information Item Factor Analysis [Computer software]. Lincolnwood,
#' IL: Scientific Software International.
#' @keywords models
#' @usage 
#' mirt(data, model, itemtype = NULL, guess = 0, upper = 1, SE = FALSE, SEtol = .001, pars = NULL, 
#' constrain = NULL, parprior = NULL, rotate = 'varimax', Target = NaN, 
#' prev.cor = NULL, quadpts = NULL, grsm.block = NULL, D = 1.702, verbose = FALSE, debug = FALSE, 
#' technical = list(), ...)
#' 
#' \S4method{summary}{ExploratoryClass}(object, rotate = '', Target = NULL, suppress = 0, digits = 3, 
#' verbose = TRUE, ...)
#' 
#' \S4method{coef}{ExploratoryClass}(object, rotate = '', Target = NULL, digits = 3,  ...)
#' 
#' \S4method{anova}{ExploratoryClass}(object, object2)
#' 
#' \S4method{fitted}{ExploratoryClass}(object, digits = 3, ...)
#' 
#' \S4method{plot}{ExploratoryClass}(x, y, type = 'info', npts = 50, theta_angle = 45, 
#' rot = list(xaxis = -70, yaxis = 30, zaxis = 10), ...)
#' 
#' \S4method{residuals}{ExploratoryClass}(object, restype = 'LD', digits = 3, df.p = FALSE, printvalue = NULL, 
#' verbose = TRUE, ...)
#' @export mirt
#' @examples
#' 
#' \dontrun{
#' #load LSAT section 7 data and compute 1 and 2 factor models
#' data(LSAT7)
#' data <- expand.table(LSAT7)
#' 
#' (mod1 <- mirt(data, 1))
#' summary(mod1)
#' residuals(mod1)
#' plot(mod1) #test information function
#' 
#' #estimated 3PL model for item 5 only
#' (mod1.3PL <- mirt(data, 1, itemtype = c('2PL', '2PL', '2PL', '2PL', '3PL')))
#' coef(mod1.3PL)
#' 
#' (mod2 <- mirt(data, 2, SE = TRUE))
#' summary(mod2, rotate = 'oblimin')
#' coef(mod2)
#' residuals(mod2)
#' plot(mod2)
#' 
#' anova(mod1, mod2) #compare the two models
#' scores <- fscores(mod2) #save factor score table
#' 
#' #confirmatory
#' cmodel <- confmirt.model()
#'    F1 = 1,4,5
#'    F2 = 2,3
#'  
#' 
#' cmod <- mirt(data, cmodel)
#' coef(cmod)
#' anova(cmod, mod2)
#' 
#' ###########
#' #data from the 'ltm' package in numeric format
#' pmod1 <- mirt(Science, 1)
#' plot(pmod1)
#' summary(pmod1)
#'
#' #Constrain all slopes to be equal
#' #first obtain parameter index
#' values <- mirt(Science,1, pars = 'values') 
#' values #note that slopes are numbered 1,5,9,13
#' (pmod1_equalslopes <- mirt(Science, 1, constrain = list(c(1,5,9,13))))
#' coef(pmod1_equalslopes)
#' 
#' pmod2 <- mirt(Science, 2)
#' summary(pmod2)
#' residuals(pmod2)
#' plot(pmod2, theta_angle = seq(0,90, by = 5)) #sum across angles of theta 1
#' itemplot(pmod2, 1)
#' anova(pmod1, pmod2)
#' 
#' 
#' ###########
#' data(SAT12)
#' data <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' 
#' mod1 <- mirt(data, 1)
#' mod2 <- mirt(data, 2, quadpts = 15)
#' mod3 <- mirt(data, 3, quadpts = 10)
#' anova(mod1,mod2)
#' anova(mod2, mod3) #negative AIC, 2 factors probably best
#' 
#' #with fixed guessing parameters
#' mod1g <- mirt(data, 1, guess = .1)
#' coef(mod1g)
#' 
#' #with estimated guessing and beta priors (for better stability)
#' itemtype <- rep('3PL', 32)
#' sv <- mirt(data, 1, itemtype, pars = 'values')
#' gindex <- sv$parnum[sv$name == 'g']
#' parprior <- list(c(gindex, 'beta', 10, 90)) 
#' mod1wg <- mirt(data, 1, itemtype, guess = .1, parprior=parprior, verbose=TRUE)
#' coef(mod1wg)
#' anova(mod1g, mod1wg)
#' 
#' ###########
#' #graded rating scale example
#' 
#' #make some data
#' a <- matrix(rep(1/1.702, 10))
#' d <- matrix(c(1,0.5,-.5,-1), 10, 4, byrow = TRUE)
#' c <- seq(-1, 1, length.out=10)
#' data <- simdata(a, d + c, 2000, itemtype = rep('graded',10))
#'
#' #use much better start values to save iterations
#' sv <- mirt(data, 1, itemtype = 'grsm', pars = 'values')
#' sv[,5] <- c(as.vector(t(cbind(a,d,c))),0,1) 
#'
#' mod1 <- mirt(data, 1)
#' mod2 <- mirt(data, 1, itemtype = 'grsm', verbose = TRUE, pars = sv)
#' coef(mod2)
#' anova(mod2, mod1) #not sig, mod2 should be prefered 
#' }
#' 
mirt <- function(data, model, itemtype = NULL, guess = 0, upper = 1, SE = FALSE, SEtol = .001,
                 pars = NULL, constrain = NULL, parprior = NULL, rotate = 'varimax', Target = NaN, 
                 prev.cor = NULL, quadpts = NULL, grsm.block = NULL, D = 1.702, verbose = FALSE, 
                 debug = FALSE, technical = list(), ...)
{   
    if(debug == 'Main') browser()
    Call <- match.call()    
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), 
                      itemtype=itemtype, guess=guess, upper=upper, grsm.block=grsm.block,
                      pars=pars, method = 'EM', constrain=constrain, SE=SE, SEtol=SEtol,
                      parprior=parprior, quadpts=quadpts, rotate=rotate, Target=Target, D=D,
                      technical = technical, debug = debug, verbose = verbose, ...)
    if(is(mod, 'ExploratoryClass') || is(mod, 'ConfirmatoryClass'))
        mod@Call <- Call
    return(mod)    
}
