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
#' on the last two, then the vector is \code{c(1,1,2,2)}. For items that should only load 
#' on the general factor (have no specific component) \code{NA} values may be used as placeholders
#' @param quadpts number of quadrature points per dimension (default 20).
#' @param SE logical; calculate information matrix and standard errors?
#' @param SE.type type of standard errors to calculate. See \code{\link{mirt}} for details
#' @param verbose logical; print observed log-likelihood value at each iteration?
#' @param ... additional arguments to be passed to the main estimation function. See \code{\link{mirt}}
#' for more details
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{confmirt.model}}, \code{\link{mirt}},
#' \code{\link{confmirt}}, \code{\link{bfactor}}, \code{\link{multipleGroup}}, \code{\link{mixedmirt}},
#' \code{\link{wald}}, \code{\link{itemplot}}, \code{\link{fscores}}, \code{\link{fitIndices}},
#' \code{\link{extract.item}}, \code{\link{iteminfo}}, \code{\link{testinfo}}, \code{\link{probtrace}},
#' \code{\link{boot.mirt}}, \code{\link{imputeMissing}}, \code{\link{itemfit}}, \code{\link{mod2values}},
#' \code{\link{read.mirt}}, \code{\link{simdata}}, \code{\link{createItem}}
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
#' bfactor(data, model, quadpts = 20, SE = FALSE, SE.type = 'SEM', verbose = TRUE, ...)
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
#' ## don't estimate specific factor for item 32
#' specific[32] <- NA
#' mod3 <- bfactor(data, specific)
#' anova(mod1, mod3)
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
bfactor <- function(data, model, quadpts = 20, SE = FALSE, SE.type = 'SEM', verbose = TRUE, ...)
{
    Call <- match.call()
    if(any(is.na(model))){
        tmpmodel <- model
        tmpmodel[is.na(tmpmodel)] <- min(model, na.rm = TRUE)
        if(!is.null(list(...)$pars))
            stop('\'pars\' argument cannot be used since model contains NA values.
                 Remove NAs from model specification.')
        tmp <- bfactor(data, tmpmodel, pars = 'values', ...)
        name <- 'a2'
        vals <- tmp$value[tmp$name == name]
        est <- tmp$est[tmp$name == name]
        vals[is.na(model)] <- 0
        est[is.na(model)] <- FALSE
        tmp$value[tmp$name == name] <- vals
        tmp$est[tmp$name == name] <- est        
        mod <- bfactor(data, tmpmodel, pars=tmp, quadpts=quadpts, SE=SE, verbose=verbose, ...)
        if(is(mod, 'ConfirmatoryClass') || is(mod, 'MultipleGroupClass'))
            mod@Call <- Call
        return(mod)
    }
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)),
                      method = 'EM', quadpts=quadpts, verbose=verbose,
                      BFACTOR = TRUE, SE=SE, ...)
    if(is(mod, 'ConfirmatoryClass') || is(mod, 'MultipleGroupClass'))
        mod@Call <- Call
    return(mod)
}
