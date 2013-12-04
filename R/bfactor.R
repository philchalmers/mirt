#' Full-Information Item Bi-factor and Two-Tier Analysis
#'
#' \code{bfactor} fits a confirmatory maximum likelihood two-tier/bifactor model to
#' dichotomous and polytomous data under the item response theory paradigm.
#' The IRT models are fit using a dimensional reduction EM algorithm so that regardless
#' of the number of specific factors estimated the model only uses the number of 
#' factors in the second-tier structure plus 1. For the bifactor model the maximum 
#' number of dimensions is only 2 since the second-tier only consists of a  
#' ubiquitous unidimensional factor. See \code{\link{mirt}} for appropriate methods to be used
#' on the objects returned from the estimation.
#'
#' \code{bfactor} follows the item factor analysis strategy explicated by
#' Gibbons and Hedeker (1992), Gibbons et al. (2007), and Cai (2010).
#' Nested models may be compared via an approximate
#' chi-squared difference test or by a reduction in AIC or BIC (accessible via
#' \code{\link{anova}}). See \code{\link{mirt}} for more details regarding the
#' IRT estimation approach used in this package. The default is to use 21 quadrature 
#' for each dimensions, but this can be over-written by passing a \code{quadpts = #}
#' argument. 
#' 
#' Note: for multiple group two-tier analyses only the second-tier means and variances 
#' should be freed since the specific factors are not treated independently due to the 
#' dimension reduction technique.  
#'
#' @aliases bfactor
#' @param data a \code{matrix} or \code{data.frame} that consists of
#'   numerically ordered data, with missing data coded as \code{NA}
#' @param model a numeric vector specifying which factor loads on which
#'   item. For example, if for a 4 item test with two specific factors, the first
#'   specific factor loads on the first two items and the second specific factor
#'   on the last two, then the vector is \code{c(1,1,2,2)}. For items that should only load 
#'   on the second-tier factors (have no specific component) \code{NA} values may 
#'   be used as place-holders. These numbers will be translated into a format suitable for
#'   \code{mirt.model()}, combined with the definition in \code{model2}, with the letter 'S' added to the 
#'   respective factor number
#' @param model2 a two-tier model specification object defined by \code{mirt.model()}. By default
#'   the model will fit a unidimensional model in the second-tier, and therefore be equivalent to 
#'   the bifactor model
#' @param group a factor variable indicating group membership used for multiple group analyses
#' @param SE logical; calculate information matrix and standard errors?
#' @param SE.type type of standard errors to calculate. See \code{\link{mirt}} for details
#' @param verbose logical; print observed log-likelihood value at each iteration?
#' @param ... additional arguments to be passed to the main estimation function. See \code{\link{mirt}}
#'   for more details
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{anova-method}}, \code{\link{coef-method}}, \code{\link{summary-method}},
#'   \code{\link{residuals-method}}, \code{\link{plot-method}}, \code{\link{fitted-method}},
#'   \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{mirt.model}}, \code{\link{mirt}},
#'   \code{\link{bfactor}}, \code{\link{multipleGroup}}, \code{\link{mixedmirt}},
#'   \code{\link{wald}}, \code{\link{itemplot}}, \code{\link{fscores}}, \code{\link{fitIndices}},
#'   \code{\link{extract.item}}, \code{\link{iteminfo}}, \code{\link{testinfo}}, \code{\link{probtrace}},
#'   \code{\link{boot.mirt}}, \code{\link{imputeMissing}}, \code{\link{itemfit}}, \code{\link{mod2values}},
#'   \code{\link{simdata}}, \code{\link{createItem}}
#' @references
#' 
#' Cai, L. (2010). A two-tier full-information item factor analysis model with applications.
#' \emph{Psychometrika, 75}, 581-612.
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
#' bfactor(data, model, model2 = mirt.model(paste0('G = 1-', ncol(data))), 
#' SE = FALSE, SE.type = 'SEM', group = NULL, verbose = TRUE, ...)
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
#' itemplot(mod1, 18, drop.zeros = TRUE) #drop the zero slopes to allow plotting
#'
#' ###Try with fixed guessing parameters added
#' guess <- rep(.1,32)
#' mod2 <- bfactor(data, specific, guess = guess)
#' coef(mod2)
#' anova(mod1, mod2)
#' 
#' ## don't estimate specific factor for item 32
#' specific[32] <- NA
#' mod3 <- bfactor(data, specific) 
#' anova(mod1, mod3)
#' 
#' # same, but decalred manually (not run)
#' #sv <- mod2values(mod1)
#' #sv$value[220] <- 0 #parameter 220 is the 32 items specific slope
#' #sv$est[220] <- FALSE
#' #mod3 <- bfactor(data, specific, pars = sv) #with excellent starting values 
#' 
#' 
#' #########
#' # mixed itemtype example 
#' 
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
#' #########
#' # Two-tier model
#' 
#' #simulate data
#' a <- matrix(c(
#' 0,1,0.5,NA,NA,
#' 0,1,0.5,NA,NA,
#' 0,1,0.5,NA,NA,
#' 0,1,0.5,NA,NA,
#' 0,1,0.5,NA,NA,
#' 0,1,NA,0.5,NA,
#' 0,1,NA,0.5,NA,
#' 0,1,NA,0.5,NA,
#' 1,0,NA,0.5,NA,
#' 1,0,NA,0.5,NA,
#' 1,0,NA,0.5,NA,
#' 1,0,NA,NA,0.5,
#' 1,0,NA,NA,0.5,
#' 1,0,NA,NA,0.5,
#' 1,0,NA,NA,0.5,
#' 1,0,NA,NA,0.5),ncol=5,byrow=TRUE)
#'
#' d <- matrix(rnorm(16))
#' items <- rep('dich', 16)
#'
#' sigma <- diag(5)
#' sigma[1,2] <- sigma[2,1] <- .4
#' set.seed(1234)
#' dataset <- simdata(a,d,2000,itemtype=items,sigma=sigma)
#'
#' specific <- c(rep(1,5),rep(2,6),rep(3,5))
#' model <- mirt.model('
#'     G1 = 1-8
#'     G2 = 9-16
#'     COV = G1*G2')
#'     
#' simmod <- bfactor(dataset, specific, model, quadpts = 9, technical = list(TOL = 1e-3))
#' coef(simmod)
#' summary(simmod)
#'
#'     }
#'
bfactor <- function(data, model, model2 = mirt.model(paste0('G = 1-', ncol(data))), 
                    SE = FALSE, SE.type = 'SEM', group = NULL, verbose = TRUE, ...)
{
    Call <- match.call()
    if(!is.numeric(model))
        stop('model must be a numeric vector')
    if(is.numeric(model))
        if(length(model) != ncol(data)) 
            stop('length of model must equal the number of items')
    nspec <- length(na.omit(unique(model)))
    specific <- model
    if(!is(model2, 'mirt.model')) 
        stop('model2 must be an appropriate second-tier model defined with mirt.model()')
    model <- bfactor2mod(model, ncol(data))
    model$x <- rbind(model2$x, model$x)
    attr(model, 'nspec') <- nspec
    attr(model, 'specific') <- specific
    if(is.null(group)) group <- rep('all', nrow(data))
    mod <- ESTIMATION(data=data, model=model, group=group,
                      method = 'EM', verbose=verbose,
                      BFACTOR = TRUE, SE=SE, ...)
    if(is(mod, 'ConfirmatoryClass') || is(mod, 'MultipleGroupClass'))
        mod@Call <- Call
    return(mod)
}
