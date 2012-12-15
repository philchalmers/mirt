#' Full-Information Item Bi-factor Analysis
#' 
#' \code{bfactor} fits a confirmatory maximum likelihood bi-factor model to
#' dichotomous and polytomous data under the item response theory paradigm. 
#' Fits univariate and multivariate 1-4PL, graded, (generalized) partial credit, 
#' nominal, multiple choice, and partially compensatory models using 
#' a dimensional reduction EM algorithm so that regardless 
#' of the number of specific factors estimated the model only uses a two-dimensional quadrature grid
#' for integration. See \code{\link{confmirt}} for appropriate methods to be used
#' on the objects returned from the estimation.
#' 
#' 
#' \code{bfactor} follows the item factor analysis strategy explicated by
#' Gibbons and Hedeker (1992) and Gibbons et al. (2007). 
#' Nested models may be compared via an approximate
#' chi-squared difference test or by a reduction in AIC or BIC (accessible via
#' \code{\link{anova}}). See \code{\link{mirt}} for more details regarding the 
#' IRT estimation approach used in this package.
#' 
#' 
#' @aliases bfactor 
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data, with missing data coded as \code{NA}
#' @param model a numeric vector specifying which factor loads on which
#' item. For example, if for a 4 item test with two specific factors, the first
#' specific factor loads on the first two items and the second specific factor
#' on the last two, then the vector is \code{c(1,1,2,2)}.
#' @param itemtype see \code{\link{mirt}} for details
#' @param grsm.block see \code{\link{mirt}} for details
#' @param rsm.block see \code{\link{mirt}} for details
#' @param guess fixed pseudo-guessing parameter. Can be entered as a single
#' value to assign a global value or may be entered as a numeric vector for
#' each item of length \code{ncol(data)}.
#' @param upper fixed upper bound parameters for 4-PL model. Can be entered as a single
#' value to assign a global guessing parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param SE logical, estimate the standard errors? Calls the MHRM subroutine for a stochastic approximation
#' @param SEtol tollerance value used to stop the MHRM estimation when \code{SE = TRUE}. Lower values
#' will take longer but may be more stable for computing the information matrix
#' @param constrain see \code{\link{mirt}} for details
#' @param parprior see \code{\link{mirt}} for details
#' @param pars see \code{\link{mirt}} for details
#' @param D a numeric value used to adjust the logistic metric to be more similar to a normal
#' cumulative density curve. Default is 1.702
#' @param prev.cor uses a previously computed correlation matrix to be used to
#' estimate starting values for the EM estimation
#' @param quadpts number of quadrature points per dimension. 
#' @param verbose logical; print observed log-likelihood value at each iteration?
#' @param debug logical; turn on debugging features?
#' @param technical see \code{\link{mirt}} for details
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{confmirt}},
#' \code{\link{fscores}}, \code{\link{multipleGroup}}, \code{\link{wald}}
#' @references
#' 
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6),
#' 1-29.
#'
#' Gibbons, R. D., & Hedeker, D. R. (1992). Full-information Item Bi-Factor
#' Analysis. \emph{Psychometrika, 57}, 423-436.
#'
#' Gibbons, R. D., Darrell, R. B., Hedeker, D., Weiss, D. J., Segawa, E., Bhaumik, D. K., 
#' Kupfer, D. J., Frank, E., Grochocinski, V. J., & Stover, A. (2007).
#' Full-Information item bifactor analysis of graded response data. 
#' \emph{Applied Psychological Measurement, 31}, 4-19
#' 
#' @keywords models
#' @usage
#' bfactor(data, model, itemtype = NULL, guess = 0, upper = 1, SE = FALSE, SEtol = .001, pars = NULL,
#' constrain = NULL, parprior = NULL, prev.cor = NULL, quadpts = 20, grsm.block = NULL, 
#' rsm.block = NULL, D = 1.702, verbose = FALSE, debug = FALSE, technical = list(), ...)
#' 
#'
#' @export bfactor
#' @examples
#' 
#' \dontrun{
#' 
#' ###load SAT12 and compute bifactor model with 3 specific factors
#' data(SAT12)
#' data <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
#' mod1 <- bfactor(data, specific)
#' summary(mod1)
#' 
#' ###Try with fixed guessing parameters added
#' guess <- rep(.1,32)
#' mod2 <- bfactor(data, specific, guess = guess)
#' coef(mod2) 
#'
#' #########
#' #simulate data
#' a <- matrix(c(
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5),ncol=3,byrow=TRUE)
#' 
#' d <- matrix(c(
#' -1.0,NA,NA,
#' -1.5,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA,
#' 2.5,1.0,-1,
#' 3.0,2.0,-0.5,
#' 3.0,2.0,-0.5,
#' 3.0,2.0,-0.5,
#' 2.5,1.0,-1,
#' 2.0,0.0,NA,
#' -1.0,NA,NA,
#' -1.5,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA),ncol=3,byrow=TRUE)
#' items <- rep('dich', 14)
#' items[5:10] <- 'graded'
#' 
#' sigma <- diag(3)
#' dataset <- simdata(a,d,2000,itemtype=items,sigma=sigma)
#'
#' specific <- c(rep(1,7),rep(2,7))
#' simmod <- bfactor(dataset, specific)
#' coef(simmod)
#'
#'     }
#' 
bfactor <- function(data, model, itemtype = NULL, guess = 0, upper = 1, SE = FALSE, SEtol = .001,
                    pars = NULL, constrain = NULL, parprior = NULL, prev.cor = NULL, quadpts = 20, 
                    grsm.block = NULL, rsm.block = NULL, D = 1.702, verbose = FALSE, debug = FALSE, 
                    technical = list(), ...)
{         
    if(debug == 'Main') browser()
    Call <- match.call()		    
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), 
                      itemtype=itemtype, guess=guess, upper=upper, grsm.block=grsm.block,
                      pars=pars, method = 'EM', constrain=constrain, SE = SE, SEtol=SEtol,
                      parprior=parprior, quadpts=quadpts, D=D, rsm.block=rsm.block,
                      technical = technical, debug = debug, verbose = verbose, 
                      BFACTOR = TRUE, ...)
    if(is(mod, 'ConfirmatoryClass') || is(mod, 'MultipleGroupClass'))
        mod@Call <- Call
    return(mod)
} 
