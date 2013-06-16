#' Full-Information Item Factor Analysis (Multidimensional Item Response
#' Theory)
#'
#' \code{mirt} fits an unconditional maximum likelihood factor analysis model
#' to dichotomous and polytomous data under the item response theory paradigm.
#' Fits univariate and multivariate Rasch, 1-4PL, graded, (generalized) partial credit,
#' nominal, graded rating scale, Rasch rating scale, nested logistic,
#' and partially compensatory models using the EM algorithm. User defined item classes
#' can also be defined using the \code{\link{createItem}} function. Models may also contain 'explanatory'
#' person or item level predictors, though these can only be included by using the
#' \code{\link{mixedmirt}} function.
#'
#' \code{mirt} follows the item factor analysis strategy by marginal maximum
#' likelihood estimation (MML) outlined in Bock and Aiken (1981), Bock,
#' Gibbons and Muraki (1988), and Muraki and Carlson (1995).
#' Nested models may be compared via the approximate
#' chi-squared difference test or by a reduction in AIC/BIC values (comparison
#' via \code{\link{anova}}). 
#'
#' \code{summary} and \code{coef} allow
#' for all the rotations available from the \code{GPArotation} package (e.g., \code{rotate = 'oblimin'})
#' as well as a \code{'promax'} rotation. Using \code{plot} will plot the test information function
#' or the test standard errors
#' for 1 and 2 dimensional solutions, or all item trace lines if only 1 dimensional the test is only
#' dichotomous items. To examine
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
#' If the model is confirmatory then the returned class will be 'ConfirmatoryClass'. Confirmatory
#' models may also contain 'explanatory' person or item level predictors, though including predictors
#' is limited only to the \code{\link{mixedmirt}} function.
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
#' initial values for polytomous items.
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
#' '4PL', 'graded', 'grsm', 'gpcm', 'rsm', 'nominal', 'PC2PL', 'PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM',
#' and '4PLNRM', for the Rasch/partial credit, 1 and 2 parameter logistic,
#' 3 parameter logistic (lower asymptote and upper), 4 parameter logistic, graded response model,
#' rating scale graded response model, generalized partial credit model, Rasch rating scale model, nominal model,
#' 2-3PL partially compensatory model, and 2-4 parameter nested logistic
#' models, respectively. User defined item classes
#' can also be defined using the \code{\link{createItem}} function
#' @param grsm.block an optional numeric vector indicating where the blocking should occur when using
#' the grsm, NA represents items that do not belong to the grsm block (other items that may be estimated
#' in the test data). For example, to specify two blocks of 3 with a 2PL item for the last item:
#' \code{grsm.block = c(rep(1,3), rep(2,3), NA)}. If NULL the all items are assumed to be within the same
#' group and therefore have the same number of item categories
#' @param rsm.block same as \code{grsm.block}, but for \code{'rsm'} blocks
#' @param key a numeric vector of the response scoring key. Required when using nested logit item types, and
#' must be the same length as the number of items used. Items that are not nested logit will ignore this vector,
#' so use \code{NA} in item locations that are not applicable
#' @param SE logical; estimate the standard errors? Calculates the information matrix from MHRM subroutine for
#' stochastic approximation, Bock and Lieberman style information (use only with small number of items), or
#' supplemented EM (SEM) computations for Bock and Lieberman style information matrix
#' @param SE.type type of estimation method to use for calculating the parameter information matrix.
#' Can be \code{'MHRM'} for stochastic estimation, \code{'BL'} for the Bock and Lieberman approach 
#' (EM only), or \code{'SEM'} for the supplemented EM. Note that for the \code{'SEM'} option increasing 
#' the number of EM cycles (\code{NCYCLES}, see below) will help to improve the accuracy. 
#' Bootstrapped standard errors are also possible but must be run with the \code{\link{boot.mirt}} function
#' @param guess fixed pseudo-guessing parameters. Can be entered as a single
#' value to assign a global guessing parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param upper fixed upper bound parameters for 4-PL model. Can be entered as a single
#' value to assign a global guessing parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param rotate type of rotation to perform after the initial orthogonal
#' parameters have been extracted by using \code{summary}; default is \code{'oblimin'}.
#' See below for list of possible rotations. If \code{rotate != ''} in the \code{summary}
#' input then the default from the object is ignored and the new rotation from the list
#' is used instead
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
#' @param calcNull logical; calculate the Null model for fit statics (e.g., TLI)? Only applicable if the
#' data contains no NA's
#' @param cl a cluster object from the \code{parallel} package (set from using \code{makeCluster(ncores)})
#' @param quadpts number of quadrature points per dimension. By default the number of quadrature uses the 
#' following scheme: \code{switch(as.character(nfact), '1'=40, '2'=20, '3'=10, '4'=7, '5'=5, 3)}
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
#' @param type type of plot to view; can be \code{'info'} to show the test
#' information function, \code{'infocontour'} for the test information contours,
#' \code{'SE'} for the test standard error function, \code{'trace'} and \code{'infotrace'}
#' for all item probability information or trace lines (only available when all items are dichotomous),
#' \code{'infoSE'} for a combined test information and standard error plot, and \code{'score'} for
#' the expected total score
#' @param theta_angle numeric values ranging from 0 to 90 used in \code{plot}. If a vector is
#' used then a bubble plot is created with the summed information across the angles specified
#' (e.g., \code{theta_angle = seq(0, 90, by=10)})
#' @param npts number of quadrature points to be used for plotting features.
#' Larger values make plots look smoother
#' @param rot allows rotation of the 3D graphics
#' @param large a logical, indicating whether the internal collapsed data should be returned,
#' or list of internally computed mirt parameters containing the data. If \code{TRUE} a list containing
#' the organized data used prior to estimation is returned. This list object can then be passed back into
#' \code{large} to avoid reorganizing the data in every estimation (useful when the dataset are very large
#' and computing the tabularized data is computationally burdensome). Therefore, the best strategy for large data
#' is to always pass the internal data to the estimation function, shown below:
#' \describe{
#' \item{Compute organized data}{e.g., \code{internaldat <- mirt(Science, 1, large = TRUE)}}
#' \item{Pass the organized data to all estimation functions}{e.g.,
#' \code{mod <- mirt(Science, 1, large = internaldat)}}
#' }
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param df.p logical; print the degrees of freedom and p-values?
#' @param nominal.highlow optional matrix indicating the highest (row 1) and lowest (row 2) categories
#' to be used for the nominal response model. Using this input may result in better numerical stability.
#' The matrix input should be a 2 by nitems numeric matrix, where each number represents the \emph{reduced}
#' category representation (mirt omits categories that are missing, so if the unique values for an item
#' are c(1,2,5,6) they are treated as being the same as c(1,2,3,4). Viewing the starting values will help
#' to identify the categories)
#' @param verbose logical; print observed log-likelihood value at each iteration?
#' @param technical a list containing lower level technical parameters for estimation. May be:
#' \describe{
#' \item{MAXQUAD}{maximum number of quadratures; default 10000}
#' \item{TOL}{EM convergence threshold; default .0001}
#' \item{MSTEPTOL}{convergence threshold for Mstep; default is \code{TOL/1000}}
#' \item{SEtol}{tolerance value used to stop the MHRM estimation when \code{SE = TRUE}
#' and \code{SE.type = 'MHRM'}. Lower values will take longer but may be more
#' stable for computing the information matrix. Default is .001}
#' \item{NCYCLES}{maximum number of EM cycles; default 500}
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
#' The parameter labels use the follow convention, here using two factors and \eqn{k} as the number
#' of categories. 
#'
#' \describe{
#' \item{Rasch}{
#' Only one intercept estimated, and the latent variance of \eqn{\theta} is freely estimated. If
#' the data have more than two categories then a partial credit model is used instead (see 'gpcm' below).
#' \deqn{P(x = 1|\theta, d) = \frac{1}{1 + exp(-(\theta + d))}}
#' }
#' \item{1-4PL}{
#' Depending on the model \eqn{u} may be equal to 1 and \eqn{g} may be equal to 0.
#' \deqn{P(x = 1|\theta, \psi) = g + \frac{(u - g)}{1 + exp(-(a_1 * \theta_1 + a_2 * \theta_2 + d))}}
#' For the 1PL model the number of factors must equal 1, and if so all the \eqn{a_1} values 
#' are constrained to be equal accross all items and the latent variance of \eqn{\theta} is 
#' freely estimated.
#' }
#' \item{graded}{
#' The graded model consists of sequential 2PL models, and here \eqn{k} is
#' the predicted category.
#' \deqn{P(x = k | \theta, \psi) = P(x \ge k | \theta, \phi) - P(x \ge k + 1 | \theta, \phi)}
#' }
#' \item{grsm}{
#' A more constrained version of the graded model where graded spacing is equal accross item blocks
#' and only adjusted by a single 'difficulty' parameter (c) while the latent variance 
#' of \eqn{\theta} is freely estimated. Again,
#' \deqn{P(x = k | \theta, \psi) = P(x \ge k | \theta, \phi) - P(x \ge k + 1 | \theta, \phi)}
#' but now
#' \deqn{P = \frac{1}{1 + exp(-(a_1 * \theta_1 + a_2 * \theta_2 + d_k + c))}}
#' }
#' \item{gpcm/nominal}{For the gpcm the \eqn{d_k} values are treated as fixed and orderd values
#' from 0:(k-1) (in the nominal model \eqn{d_0} is also set to 0). Additionally, for identification
#' in the nominal model \eqn{ak_0 = 1}, \eqn{ak_k = (k - 1)}.
#' \deqn{P(x = k | \theta, \psi) = \frac{exp(-ak_k * (a_1 * \theta_1 + a_2 * \theta_2) + d_k)}
#' {\sum_i^k exp(-ak_k * (a_1 * \theta_1 + a_2 * \theta_2) + d_k)}}
#' 
#' For partial credit model (when \code{itemtype = 'Rasch'}; unidimenional only) the above model 
#' is further constrained so that \eqn{ak_k = (0,1,\ldots, k-1)}, \eqn{a_1 = 1}, and the latent 
#' variance of \eqn{\theta_1} is freely estimated. 
#' 
#' In the nominal model this parametrizations helps to identify the empirical ordering of the 
#' categories by inspecting the \eqn{ak} values. Larger values indicate that the item category is 
#' more positively related to the latent trait(s) being measured. For instance, if an item was 
#' truly ordinal (such as a Likert scale), and had 4 response categories, we would expect 
#' to see \eqn{ak_0 < ak_1 < ak_2 < ak_3} following estimation. If on the other hand \eqn{ak_0 > ak_1}
#' then it would appear that the second category is less related to to the trait than the first, and 
#' therefore the second category should be understood as the 'lowest score'.
#' 
#' }
#' \item{rsm}{
#' A more constrained version of the partial credit model where the spacing is equal
#' accross item blocks and only adjusted by a single 'difficulty' parameter (c). Note that this is 
#' analogous to the relationship between the graded model and the grsm (with an additional 
#' constraint regarding the fixed discrimination parameters; the discrimination constraint can, 
#' however, be relaxed by adjusting the starting values specifications manually and applying 
#' additional equality constraints).
#' }
#' \item{partcomp}{Partially compensatory models consist of the products of 2PL probability curves.
#' \deqn{P(x = 1 | \theta, \psi) = g + (1 - g) (\frac{1}{1 + exp(-(a_1 * \theta_1 + d_1))} *
#' \frac{1}{1 + exp(-(a_2 * \theta_2 + d_2))})}
#' }
#' \item{1-4PLNRM}{Nested logistic curves for modeling distractor items. Requires a scoring key.
#' The model is broken into two components for the probability of endorsement. For successful endorsement
#' the probability trace is the 1-4PL model, while for unsuccessful endorsement:
#' \deqn{P(x = 0 | \theta, \psi) = (1 - P_{1-4PL}(x = 1 | \theta, \psi)) * P_{nominal}(x = k | \theta, \psi)}
#' which is the product of the compliment of the dichotomous trace line with the nominal
#' response model. In the nominal model, the slope parameters defined above are constrained to be 1's,
#' while the last value of the \eqn{ak} is freely estimated.
#' }
#' }
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{confmirt.model}}, \code{\link{mirt}},
#' \code{\link{confmirt}}, \code{\link{bfactor}}, \code{\link{multipleGroup}}, \code{\link{mixedmirt}},
#' \code{\link{wald}}, \code{\link{itemplot}}, \code{\link{fscores}}, \code{\link{fitIndices}},
#' \code{\link{extract.item}}, \code{\link{iteminfo}}, \code{\link{testinfo}}, \code{\link{probtrace}},
#' \code{\link{boot.mirt}}, \code{\link{imputeMissing}}, \code{\link{itemfit}}, \code{\link{mod2values}},
#' \code{\link{read.mirt}}, \code{\link{simdata}}, \code{\link{createItem}}
#'
#' @references
#'
#' Andrich, D. (1978). A rating scale formulation for ordered response categories.
#' \emph{Psychometrika, 43}, 561-573.
#'
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of
#' item parameters: Application of an EM algorithm. \emph{Psychometrika,
#' 46}(4), 443-459.
#'
#' Bock, R. D., Gibbons, R., & Muraki, E. (1988). Full-Information Item Factor
#' Analysis. \emph{Applied Psychological Measurement, 12}(3), 261-280.
#'
#' Bock, R. D. & Lieberman, M. (1970). Fitting a response model for n dichotomously
#' scored items. \emph{Psychometrika, 35}, 179-197.
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#'
#' Lord, F. M. & Novick, M. R. (1968). Statistical theory of mental test scores. Addison-Wesley.
#'
#' Rasch, G. (1960). Probabilistic models for some intelligence and attainment tests.
#' \emph{Danish Institute for Educational Research}.
#'
#' Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm.
#' \emph{Applied Psychological Measurement, 16}, 159-176.
#'
#' Muraki, E. & Carlson, E. B. (1995). Full-information factor analysis for polytomous
#' item responses. \emph{Applied Psychological Measurement, 19}, 73-90.
#'
#' Samejima, F. (1969). Estimation of latent ability using a response pattern of
#' graded scores. \emph{Psychometrika Monographs}, 34.
#'
#' Suh, Y. & Bolt, D. (2010). Nested logit models for multiple-choice item response data.
#' \emph{Psychometrika, 75}, 454-473.
#'
#' Sympson, J. B. (1977). A model for testing with multidimensional items.
#' Proceedings of the 1977 Computerized Adaptive Testing Conference.
#'
#' Thissen, D. (1982). Marginal maximum likelihood estimation for the one-parameter logistic model.
#' \emph{Psychometrika, 47}, 175-186.
#'
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., &
#' Bock, R. D. (2003). \emph{TESTFACT 4 for Windows: Test Scoring, Item Statistics,
#' and Full-information Item Factor Analysis} [Computer software]. Lincolnwood,
#' IL: Scientific Software International.
#'
#' @keywords models
#' @usage
#' mirt(data, model, itemtype = NULL, guess = 0, upper = 1, SE = FALSE, SE.type = 'SEM', pars = NULL,
#' constrain = NULL, parprior = NULL, calcNull = TRUE, rotate = 'oblimin', Target = NaN,
#' quadpts = NULL, grsm.block = NULL, rsm.block = NULL, key=  NULL, nominal.highlow = NULL,
#' cl = NULL, large = FALSE, verbose = TRUE, technical = list(), ...)
#'
#' \S4method{summary}{ExploratoryClass}(object, rotate = '', Target = NULL, suppress = 0, digits = 3,
#' verbose = TRUE, ...)
#'
#' \S4method{coef}{ExploratoryClass}(object, rotate = '', Target = NULL, digits = 3, 
#' verbose = TRUE, ...)
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
#' data <- expand.table(LSAT7)
#'
#' (mod1 <- mirt(data, 1))
#' coef(mod1)
#' coef(mod2 <- mirt(data, 1, SE = TRUE)) #standard errors with SEM method
#' coef(mod3 <- mirt(data, 1, SE = TRUE, SE.type = 'BL')) #standard errors with BL method
#' residuals(mod1)
#' plot(mod1) #test information function
#' plot(mod1, type = 'trace') #trace lines
#'
#' #estimated 3PL model for item 5 only
#' (mod1.3PL <- mirt(data, 1, itemtype = c('2PL', '2PL', '2PL', '2PL', '3PL')))
#' coef(mod1.3PL)
#'
#' #two factors (exploratory)
#' mod2 <- mirt(data, 2) 
#' 
#' #too few iterations, try running more using current model as new starting 
#' #   values (could also increase NCYCLES and rerun) 
#' mod2 <- mirt(data, 2, pars = mod2values(mod2))
#' coef(mod2)
#' summary(mod2, rotate = 'oblimin') #oblimin rotation
#' residuals(mod2)
#' plot(mod2)
#'
#' anova(mod1, mod2) #compare the two models
#' scores <- fscores(mod2) #save factor score table
#' scoresfull <- fscores(mod2, full.scores = TRUE, scores.only = TRUE) #factor scores for original data
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
#' fitIndices(pmod1) #M2 limited information statistic
#'
#' #Constrain all slopes to be equal with the constrain = list() input
#' #first obtain parameter index
#' values <- mirt(Science,1, pars = 'values')
#' values #note that slopes are numbered 1,5,9,13, or index with values$parnum[values$name == 'a1']
#' (pmod1_equalslopes <- mirt(Science, 1, constrain = list(c(1,5,9,13))))
#' 
#' coef(pmod1_equalslopes)
#' anova(pmod1_equalslopes, pmod1) #significantly worse fit with almost all criteria
#'
#' pmod2 <- mirt(Science, 2, technical = list(NCYCLES = 2000))
#' summary(pmod2)
#' plot(pmod2)
#' itemplot(pmod2, 1)
#' anova(pmod1, pmod2)
#'
#' #unidimensional fit with a generalized partial credit and nominal model
#' (gpcmod <- mirt(Science, 1, 'gpcm'))
#' coef(gpcmod)
#' 
#' #for the nominal model the lowest and highest categories are assumed to be the theoretically lowest
#' #  and highest categories that related to the latetent trait(s), however a custom nominal.highlow matrix 
#' #  can be passed to declare which item category should be treated as the 'highest' and 'lowest' instead 
#' (nomod <- mirt(Science, 1, 'nominal'))
#' coef(nomod) #ordering of ak values suggest that the items are indeed ordinal 
#' anova(gpcmod, nomod)
#' itemplot(nomod, 3)
#'
#' ###########
#' #empirical dimensionality testing that includes 'guessing'
#'
#' data(SAT12)
#' data <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#'
#' mod1 <- mirt(data, 1)
#' mod2 <- mirt(data, 2)
#' mod3 <- mirt(data, 3, technical = list(TOL = 1e-3)) #difficulty converging to 1e-4 with reduced quadpts
#' anova(mod1,mod2)
#' anova(mod2, mod3) #negative AIC, 2 factors probably best
#'
#' #with fixed guessing parameters
#' mod1g <- mirt(data, 1, guess = .1)
#' coef(mod1g)
#'
#' ###########
#' #graded rating scale example
#'
#' #make some data
#' set.seed(1234)
#' a <- matrix(rep(1, 10))
#' d <- matrix(c(1,0.5,-.5,-1), 10, 4, byrow = TRUE)
#' c <- seq(-1, 1, length.out=10)
#' data <- simdata(a, d + c, 2000, itemtype = rep('graded',10))
#'
#' #use much better start values to save iterations
#' sv <- mirt(data, 1, itemtype = 'grsm', pars = 'values')
#' sv[,'value'] <- c(as.vector(t(cbind(a,d,c))),0,1)
#' 
#' #also possible to edit start values with a GUI approach with
#' #   sv <- edit(sv)
#'
#' mod1 <- mirt(data, 1)
#' mod2 <- mirt(data, 1, itemtype = 'grsm', pars = sv)
#' coef(mod2)
#' anova(mod2, mod1) #not sig, mod2 should be preferred
#'
#' ###########
#' # 2PL nominal response model example (Suh and Bolt, 2010)
#' data(SAT12)
#' SAT12[SAT12 == 8] <- NA
#' head(SAT12)
#'
#' #correct answer key
#' key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
#' scoredSAT12 <- key2binary(SAT12, key)
#' mod0 <- mirt(scoredSAT12, 1)
#'
#' #for first 5 items use 2PLNRM and nominal
#' scoredSAT12[,1:5] <- as.matrix(SAT12[,1:5])
#' mod1 <- mirt(scoredSAT12, 1, c(rep('nominal',5),rep('2PL', 27)))
#' mod2 <- mirt(scoredSAT12, 1, c(rep('2PLNRM',5),rep('2PL', 27)), key=key)
#' coef(mod0)$Item.1
#' coef(mod1)$Item.1
#' coef(mod2)$Item.1
#' itemplot(mod0, 1)
#' itemplot(mod1, 1)
#' itemplot(mod2, 1)
#'
#' #compare added information from distractors
#' Theta <- matrix(seq(-4,4,.01))
#' par(mfrow = c(2,3))
#' for(i in 1:5){
#'     info <- iteminfo(extract.item(mod0,i), Theta)
#'     info2 <- iteminfo(extract.item(mod2,i), Theta)
#'     plot(Theta, info2, type = 'l', main = paste('Information for item', i), ylab = 'Information')
#'     lines(Theta, info, col = 'red')
#' }
#'
#' #test information
#' par(mfrow = c(1,1))
#' plot(Theta, testinfo(mod2, Theta), type = 'l', main = 'Test information', ylab = 'Information')
#' lines(Theta, testinfo(mod0, Theta), col = 'red')
#'
#'
#' }
mirt <- function(data, model, itemtype = NULL, guess = 0, upper = 1, SE = FALSE, SE.type = 'SEM',
                 pars = NULL, constrain = NULL, parprior = NULL, calcNull = TRUE, rotate = 'oblimin',
                 Target = NaN, quadpts = NULL, grsm.block = NULL, rsm.block = NULL,
                 key = NULL, nominal.highlow = NULL, cl = NULL, 
                 large = FALSE, verbose = TRUE, technical = list(), ...)
{
    Call <- match.call()
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)),
                      itemtype=itemtype, guess=guess, upper=upper, grsm.block=grsm.block,
                      pars=pars, method = 'EM', constrain=constrain, SE=SE,
                      parprior=parprior, quadpts=quadpts, rotate=rotate, Target=Target, 
                      rsm.block=rsm.block, technical=technical, verbose=verbose,
                      calcNull=calcNull, SE.type=SE.type, cl=cl, large=large, key=key, 
                      nominal.highlow=nominal.highlow, ...)
    if(is(mod, 'ExploratoryClass') || is(mod, 'ConfirmatoryClass'))
        mod@Call <- Call
    return(mod)
}
