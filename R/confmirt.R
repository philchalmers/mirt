#' Confirmatory Full-Information Item Factor Analysis
#' 
#' \code{confmirt} fits a conditional (i.e., confirmatory) full-information
#' maximum-likelihood factor analysis model to dichotomous and polytomous
#' data under the item response theory paradigm using Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm. Fits univariate and multivariate Rasch, 
#' 1-4PL, graded, (generalized) partial credit, nominal, multiple choice, graded rating scale, Rasch rating scale,
#' and partially-compensatory models, potentially with polynomial and product constructed latent traits. 
#' Models may also contain 'explanatory' 
#' person or item level predictors, though these can only be included by using the 
#' \code{\link{mixedmirt}} function.
#'  
#' \code{confmirt} follows a confirmatory and exploratory item factor analysis strategy that
#' uses a stochastic version of maximum likelihood estimation described by Cai
#' (2010a, 2010b). The general equation used for multidimensional item response theory
#' in this function is in the logistic form with a scaling correction of 1.702.
#' This correction is applied to allow comparison to mainstream programs such
#' as TESTFACT (2003) and POLYFACT. Missing data are treated as 'missing at
#' random' so that each response vector is included in the estimation (i.e.,
#' full-information). Residuals are computed using the LD statistic (Chen &
#' Thissen, 1997) in the lower diagonal of the matrix returned by
#' \code{residuals}, and Cramer's V above the diagonal. For computing the
#' log-likelihood more accurately see \code{\link{calcLogLik}}.
#' 
#' \code{coef} displays the item parameters with their associated standard
#' errors, while use of \code{summary} transforms the slopes into a factor
#' loadings metric and if the model is exploratory allows for rotating the parameters. 
#' Also, nested models may be compared by using the
#' \code{anova} function, where a Chi-squared difference test and AIC/BIC
#' difference values are displayed.
#' 
#' @section Convergence:
#' 
#' The MHRM algorithm often is more stable than the EM counterpart in \code{mirt} but 
#' convergence of the algorithm should be interpreted with caution. When the number of iterations 
#' grows very high (e.g., greater than 1500) or when \code{Max Change = .2500} values are repeatedly printed
#' to the console too often (indicating that the parameters were being constrained since they are naturally 
#' moving in steps greater than 0.25) then the model may either be ill defined or have a 
#' very flat likelihood surface, and genuine maximum likelihood parameter estimates may be difficult to find. 
#' 
#' @section Confirmatory IRT:
#' 
#' Specification of the confirmatory item factor analysis model follows many of
#' the rules in the SEM framework for confirmatory factor analysis. The
#' variances of the latent factors are automatically fixed to 1 to help
#' facilitate model identification. All parameters may be fixed to constant
#' values or set equal to other parameters using the appropriate declarations.
#' 
#' @section Exploratory IRT:
#' 
#' Specifying a number as the second input to confmirt an exploratory IRT model is estimated and 
#' can be viewed as a stochastic analogue of \code{mirt}, with much of the same behaviour and 
#' specifications. Rotation and target matrix options will be used in this subroutine and will be
#' passed to the returned object for use in generic functions such as \code{summary()} and 
#' \code{fscores}. Again, factor means and variances are fixed to ensure proper identification. See
#' \code{\link{mirt}} for more details.
#' 
#' 
#' @aliases confmirt coef,ConfirmatoryClass-method summary,ConfirmatoryClass-method
#' residuals,ConfirmatoryClass-method anova,ConfirmatoryClass-method fitted,ConfirmatoryClass-method
#' plot,ConfirmatoryClass-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data, with missing data coded as \code{NA}
#' @param model an object returned from \code{confmirt.model()} declaring how
#' the factor model is to be estimated, or a single numeric value indicating the number 
#' of exploratory factors to estimate. See \code{\link{confmirt.model}} for
#' more details
#' @param guess initial (or fixed) values for the pseudo-guessing parameter. Can be 
#' entered as a single value to assign a global guessing parameter or may be entered as
#' a numeric vector for each item
#' @param upper initial (or fixed) upper bound parameters for 4-PL model. Can be 
#' entered as a single value to assign a global upper bound parameter or may be entered as a 
#' numeric vector corresponding to each item
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#' option. Only prints patterns that have standardized residuals greater than 
#' \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param verbose logical; display iteration history during estimation?
#' @param draws the number of Monte Carlo draws to estimate the log-likelihood
#' @param calcNull logical; calculate the Null model for fit statics (e.g., TLI)?
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param itemtype see \code{\link{mirt}} for details
#' @param grsm.block see \code{\link{mirt}} for details
#' @param rsm.block see \code{\link{mirt}} for details
#' @param constrain see \code{\link{mirt}} for details
#' @param parprior see \code{\link{mirt}} for details
#' @param pars see \code{\link{mirt}} for details
#' @param debug logical; turn on debugging features?
#' @param object an object of class \code{ConfirmatoryClass}
#' @param object2 an object of class \code{ConfirmatoryClass}
#' @param digits the number of significant digits to be rounded
#' @param rotate if \code{model} is numeric (indicating an exploratory item FA) then this 
#' rotation is used. Default is \code{'varimax'}
#' @param Target a dummy variable matrix indicting a target rotation pattern
#' @param suppress a numeric value indicating which factor
#' loadings should be suppressed. Typical values are around .3 in most
#' statistical software. Default is 0 for no suppression
#' @param D a numeric value used to adjust the logistic metric to be more similar to a normal
#' cumulative density curve. Default is 1.702
#' @param technical list specifying subtle parameters that can be adjusted. These 
#' values are 
#' @param df.p logical; print the degrees of freedom and p-values?
#' @param x an object of class \code{mirt} to be plotted or printed
#' @param y an unused variable to be ignored
#' @param type type of plot to view; can be \code{'info'} to show the test
#' information function, \code{'infocontour'} for the test information contours, 
#' or \code{'SE'} for the test standard error function
#' @param theta_angle numeric values ranging from 0 to 90 used in \code{plot}. If a vector is 
#' used then a bubble plot is created with the summed information across the angles specified 
#' (e.g., \code{theta_angle = seq(0, 90, by=10)})
#' @param npts number of quadrature points to be used for plotting features.
#' Larger values make plots look smoother
#' @param rot allows rotation of the 3D graphics
#' \describe{
#' \item{NCYCLES}{max number of MH-RM cycles; default 2000}
#' \item{BURNIN}{number of burn in cycles (stage 1); default 150}
#' \item{SEMCYCLES}{number of SEM cycles (stage 2); default 50}
#' \item{KDRAWS}{number of parallel MH sets to be drawn; default 1}
#' \item{TOL}{minimum threshold tolerance for convergence of MH-RM, must occur on three consecutive
#' occations; default .001} 
#'   \item{set.seed}{seed number used during estimation. Default is 12345} 	 
#'   \item{gain}{a vector of three values specifying the numerator, exponent, and subtracted
#'      values for the RM gain value. Default is \code{c(0.05,0.5,0.004)}}   	
#' }
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{simdata}},
#' \code{\link{fscores}}, \code{\link{confmirt.model}}, \code{\link{wald}}, 
#' \code{\link{multipleGroup}}, \code{\link{itemplot}}, \code{\link{fitIndices}}, 
#' \code{\link{mixedmirt}}, \code{\link{testinfo}}, \code{\link{iteminfo}},
#' @references
#' 
#' Cai, L. (2010a). High-Dimensional exploratory item factor analysis by a
#' Metropolis-Hastings Robbins-Monro algorithm. \emph{Psychometrika, 75},
#' 33-57.
#' 
#' Cai, L. (2010b). Metropolis-Hastings Robbins-Monro algorithm for confirmatory
#' item factor analysis. \emph{Journal of Educational and Behavioral
#' Statistics, 35}, 307-335.
#' 
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6),
#' 1-29.
#'
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., &
#' Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics,
#' and Full-information Item Factor Analysis [Computer software]. Lincolnwood,
#' IL: Scientific Software International.
#' @keywords models
#' @usage 
#' confmirt(data, model, itemtype = NULL, guess = 0, upper = 1, pars = NULL, 
#' constrain = NULL, parprior = NULL, calcNull = TRUE, grsm.block = NULL, rsm.block = NULL, verbose = TRUE, 
#' draws = 2000, debug = FALSE, rotate = 'varimax', Target = NULL, D = 1.702, 
#' technical = list(),  ...)
#' 
#' \S4method{summary}{ConfirmatoryClass}(object, suppress = 0, digits = 3, verbose = TRUE, ...)
#' 
#' \S4method{coef}{ConfirmatoryClass}(object, digits = 3, ...)
#' 
#' \S4method{anova}{ConfirmatoryClass}(object, object2)
#' 
#' \S4method{fitted}{ConfirmatoryClass}(object, digits = 3, ...)
#' 
#' \S4method{plot}{ConfirmatoryClass}(x, y, type = 'info', npts = 50, theta_angle = 45, 
#' rot = list(xaxis = -70, yaxis = 30, zaxis = 10), ...)
#' 
#' \S4method{residuals}{ConfirmatoryClass}(object, restype = 'LD', digits = 3, df.p = FALSE, printvalue = NULL, 
#' verbose = TRUE, ...)
#'
#' @export confmirt
#' @examples
#'  
#' \dontrun{
#' #Exploratory model estimation, similar to mirt()
#' data(LSAT7)
#' fulldata <- expand.table(LSAT7) 
#' (mod1 <- confmirt(fulldata, 1))
#' 
#' #Confirmatory models
#' 
#' #simulate data
#' a <- matrix(c(
#' 1.5,NA,
#' 0.5,NA,
#' 1.0,NA,
#' 1.0,0.5,
#'  NA,1.5,
#'  NA,0.5,
#'  NA,1.0,
#'  NA,1.0),ncol=2,byrow=TRUE)
#' 
#' d <- matrix(c(
#' -1.0,NA,NA,
#' -1.5,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA,
#' 3.0,2.0,-0.5,
#' 2.5,1.0,-1,
#' 2.0,0.0,NA,
#' 1.0,NA,NA),ncol=3,byrow=TRUE)
#' 
#' sigma <- diag(2)
#' sigma[1,2] <- sigma[2,1] <- .4
#' items <- c(rep('dich',4), rep('graded',3), 'dich')
#' dataset <- simdata(a,d,2000,items,sigma)
#' 
#' #analyses
#' #CIFA for 2 factor crossed structure
#' 
#' model.1 <- confmirt.model()
#'   F1 = 1-4
#'   F2 = 4-8
#'   COV = F1*F2
#' 
#' 
#' #compute model, and use parallel computation of the log-likelhood
#' library(parallel)
#' cl <- makeCluster(4)
#' mod1 <- confmirt(dataset, model.1, cl=cl)
#' coef(mod1)
#' summary(mod1)
#' residuals(mod1)
#' 
#' #####
#' #bifactor 	
#' model.3 <- confmirt.model()
#'   G = 1-8
#'   F1 = 1-4
#'   F2 = 5-8
#' 
#' 
#' mod3 <- confmirt(dataset,model.3)
#' coef(mod3)
#' summary(mod3)
#' residuals(mod3)
#' anova(mod1,mod3)
#'
#' #####
#' #polynomial/combinations
#' data(SAT12)
#' data <- key2binary(SAT12,
#'                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' 
#' model.quad <- confmirt.model()
#'        F1 = 1-32
#'   (F1*F1) = 1-32       
#'   
#' 
#' model.combo <- confmirt.model()
#'        F1 = 1-16
#'        F2 = 17-32
#'   (F1*F2) = 1-8
#'
#' 
#' (mod.quad <- confmirt(data, model.quad))
#' (mod.combo <- confmirt(data, model.combo))
#' anova(mod.quad, mod.combo)
#' 
#' }
#' 
confmirt <- function(data, model, itemtype = NULL, guess = 0, upper = 1, pars = NULL, 
                     constrain = NULL, parprior = NULL, calcNull = TRUE, grsm.block = NULL, rsm.block = NULL, 
                     verbose = TRUE, draws = 3000, debug = FALSE, rotate = 'varimax', Target = NULL, 
                     D = 1.702, technical = list(),  ...)
{   
    if(debug == 'Main') browser()
    Call <- match.call()    
    mod <- ESTIMATION(data=data, model=model, group = rep('all', nrow(data)), itemtype=itemtype, 
                      guess=guess, upper=upper, grsm.block=grsm.block, D=D, calcNull=calcNull,
                      pars=pars, constrain=constrain, parprior=parprior, verbose=verbose, 
                      rsm.block=rsm.block, draws=draws, debug=debug, technical = list(),  ...)
    if(is(mod, 'ExploratoryClass') || is(mod, 'ConfirmatoryClass'))
        mod@Call <- Call
    return(mod)    
}
