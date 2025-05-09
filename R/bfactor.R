#' Full-Information Item Bi-factor and Two-Tier Analysis
#'
#' \code{bfactor} fits a confirmatory maximum likelihood two-tier/bifactor/testlet model to
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
#' IRT estimation approach used in this package.
#'
#' The two-tier model has a specific block diagonal covariance structure between the primary and
#' secondary latent traits. Namely, the secondary latent traits are assumed to be orthogonal to
#' all traits and have a fixed variance of 1, while the primary traits can be organized to vary
#' and covary with other primary traits in the model.
#'
#' \deqn{\Sigma_{two-tier} = \left(\begin{array}{cc} G & 0 \\ 0 & diag(S) \end{array} \right)}
#'
#' The bifactor model is a special case of the two-tier model when \eqn{G} above is a 1x1 matrix,
#' and therefore only 1 primary factor is being modeled. Evaluation of the numerical integrals
#' for the two-tier model requires only \eqn{ncol(G) + 1} dimensions for integration since the
#' \eqn{S} second order (or 'specific') factors require only 1 integration grid due to the
#' dimension reduction technique.
#'
#' Note: for multiple group two-tier analyses only the second-tier means and variances
#' should be freed since the specific factors are not treated independently due to the
#' dimension reduction technique.
#'
#' @return function returns an object of class \code{SingleGroupClass}
#'   (\link{SingleGroupClass-class}) or \code{MultipleGroupClass}(\link{MultipleGroupClass-class}).
#'
#' @aliases bfactor
#' @param data a \code{matrix} or \code{data.frame} that consists of
#'   numerically ordered data, organized in the form of integers,
#'   with missing data coded as \code{NA}
#' @param model a numeric vector specifying which factor loads on which
#'   item. For example, if for a 4 item test with two specific factors, the first
#'   specific factor loads on the first two items and the second specific factor
#'   on the last two, then the vector is \code{c(1,1,2,2)}. For items that should only load
#'   on the second-tier factors (have no specific component) \code{NA} values may
#'   be used as place-holders. These numbers will be translated into a format suitable for
#'   \code{mirt.model()}, combined with the definition in \code{model2}, with the letter 'S'
#'   added to the respective factor number
#'
#'   Alternatively, input can be specified using the \code{\link{mirt.model}} syntax with the
#'   restriction that each item must load on exactly one specific factor (or no specific factors,
#'   if it is only predicted by the general factor specified in \code{model2})
#' @param model2 a two-tier model specification object defined by \code{mirt.model()} or
#'   a string to be passed to \code{\link{mirt.model}}. By default
#'   the model will fit a unidimensional model in the second-tier, and therefore be equivalent to
#'   the bifactor model
#' @param group a factor variable indicating group membership used for multiple group analyses
#' @param quadpts number of quadrature nodes to use after accounting for the reduced number of dimensions.
#'   Scheme is the same as the one used in \code{\link{mirt}}, however it is in regards to the reduced
#'   dimensions (e.g., a bifactor model has 2 dimensions to be integrated)
#' @param invariance see \code{\link{multipleGroup}} for details, however, the specific factor variances
#'   and means will be constrained according to the dimensional reduction algorithm
#' @param ... additional arguments to be passed to the estimation engine. See \code{\link{mirt}}
#'   for more details and examples
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{mirt}}
#' @references
#'
#' Cai, L. (2010). A two-tier full-information item factor analysis model with applications.
#' \emph{Psychometrika, 75}, 581-612.
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Bradlow, E.T., Wainer, H., & Wang, X. (1999). A Bayesian random effects model for testlets.
#' \emph{Psychometrika, 64}, 153-168.
#'
#' Gibbons, R. D., & Hedeker, D. R. (1992). Full-information Item Bi-Factor
#' Analysis. \emph{Psychometrika, 57}, 423-436.
#'
#' Gibbons, R. D., Darrell, R. B., Hedeker, D., Weiss, D. J., Segawa, E., Bhaumik, D. K.,
#' Kupfer, D. J., Frank, E., Grochocinski, V. J., & Stover, A. (2007).
#' Full-Information item bifactor analysis of graded response data.
#' \emph{Applied Psychological Measurement, 31}, 4-19.
#'
#' Wainer, H., Bradlow, E.T., & Wang, X. (2007). Testlet response theory and its applications.
#' New York, NY: Cambridge University Press.
#'
#' @keywords models
#' @export bfactor
#' @examples
#'
#' \donttest{
#'
#' ### load SAT12 and compute bifactor model with 3 specific factors
#' data(SAT12)
#' data <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
#' mod1 <- bfactor(data, specific)
#' summary(mod1)
#' itemplot(mod1, 18, drop.zeros = TRUE) #drop the zero slopes to allow plotting
#'
#' # alternative model definition via ?mirt.model syntax
#' specific2 <- "S1 = 7,9,10,11,13,15,17,18,21,22,24,27,31
#'               S2 = 1,3,6,8,16,29,32
#'               S3 = 2,4,5,12,14,19,20,23,25,26,28,30"
#' mod2 <- bfactor(data, specific2)
#' anova(mod1, mod2) # same
#'
#' # also equivalent using item names instead (not run)
#' specific3 <- "S1 = Item.7, Item.9, Item.10, Item.11, Item.13, Item.15,
#'                 Item.17, Item.18, Item.21, Item.22, Item.24, Item.27, Item.31
#'               S2 = Item.1, Item.3, Item.6, Item.8, Item.16, Item.29, Item.32
#'               S3 = Item.2, Item.4, Item.5, Item.12, Item.14, Item.19,
#'                 Item.20, Item.23, Item.25, Item.26, Item.28, Item.30"
#' # mod3 <- bfactor(data, specific3)
#' # anova(mod1, mod2, mod3)  # all same
#'
#' ### Try with fixed guessing parameters added
#' guess <- rep(.1,32)
#' mod2 <- bfactor(data, specific, guess = guess)
#' coef(mod2)
#' anova(mod1, mod2)
#'
#' ## don't estimate specific factor for item 32
#' specific[32] <- NA
#' mod3 <- bfactor(data, specific)
#' anova(mod3, mod1)
#'
#' # same, but with syntax (not run)
#' specific3 <- "S1 = 7,9,10,11,13,15,17,18,21,22,24,27,31
#'               S2 = 1,3,6,8,16,29
#'               S3 = 2,4,5,12,14,19,20,23,25,26,28,30"
#' # mod3b <- bfactor(data, specific3)
#' # anova(mod3b)
#'
#'
#' #########
#' # mixed itemtype example
#'
#' # simulate data
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
#' items <- rep('2PL', 14)
#' items[5:10] <- 'graded'
#'
#' sigma <- diag(3)
#' dataset <- simdata(a,d,5000,itemtype=items,sigma=sigma)
#' itemstats(dataset)
#'
#' specific <- "S1 = 1-7
#'              S2 = 8-14"
#' simmod <- bfactor(dataset, specific)
#' coef(simmod, simplify=TRUE)
#'
#'
#' #########
#' # General testlet response model (Wainer, 2007)
#'
#' # simulate data
#' set.seed(1234)
#' a <- matrix(0, 12, 4)
#' a[,1] <- rlnorm(12, .2, .3)
#' ind <- 1
#' for(i in 1:3){
#'    a[ind:(ind+3),i+1] <- a[ind:(ind+3),1]
#'    ind <- ind+4
#' }
#' print(a)
#' d <- rnorm(12, 0, .5)
#' sigma <- diag(c(1, .5, 1, .5))
#' dataset <- simdata(a,d,2000,itemtype=rep('2PL', 12),sigma=sigma)
#' itemstats(dataset)
#'
#' # estimate by applying constraints and freeing the latent variances
#' specific <- "S1 = 1-4
#'              S2 = 5-8
#'              S3 = 9-12"
#' model <- "G = 1-12
#'           CONSTRAIN = (1, a1, a2), (2, a1, a2), (3, a1, a2), (4, a1, a2),
#'             (5, a1, a3), (6, a1, a3), (7, a1, a3), (8, a1, a3),
#'             (9, a1, a4), (10, a1, a4), (11, a1, a4), (12, a1, a4)
#'           COV = S1*S1, S2*S2, S3*S3"
#'
#' simmod <- bfactor(dataset, specific, model)
#' coef(simmod, simplify=TRUE)
#'
#' # Constrained testlet model (Bradlow, 1999)
#' model2 <- "G = 1-12
#'           CONSTRAIN = (1, a1, a2), (2, a1, a2), (3, a1, a2), (4, a1, a2),
#'             (5, a1, a3), (6, a1, a3), (7, a1, a3), (8, a1, a3),
#'             (9, a1, a4), (10, a1, a4), (11, a1, a4), (12, a1, a4),
#'             (GROUP, COV_22, COV_33, COV_44)
#'           COV = S1*S1, S2*S2, S3*S3"
#'
#' simmod2 <- bfactor(dataset, specific, model2)
#' coef(simmod2, simplify=TRUE)
#' anova(simmod2, simmod)
#'
#'
#' #########
#' # Two-tier model
#'
#' # simulate data
#' set.seed(1234)
#' a <- matrix(c(
#'   0,1,0.5,NA,NA,
#'   0,1,0.5,NA,NA,
#'   0,1,0.5,NA,NA,
#'   0,1,0.5,NA,NA,
#'   0,1,0.5,NA,NA,
#'   0,1,NA,0.5,NA,
#'   0,1,NA,0.5,NA,
#'   0,1,NA,0.5,NA,
#'   1,0,NA,0.5,NA,
#'   1,0,NA,0.5,NA,
#'   1,0,NA,0.5,NA,
#'   1,0,NA,NA,0.5,
#'   1,0,NA,NA,0.5,
#'   1,0,NA,NA,0.5,
#'   1,0,NA,NA,0.5,
#'   1,0,NA,NA,0.5),ncol=5,byrow=TRUE)
#'
#' d <- matrix(rnorm(16))
#' items <- rep('2PL', 16)
#'
#' sigma <- diag(5)
#' sigma[1,2] <- sigma[2,1] <- .4
#' dataset <- simdata(a,d,2000,itemtype=items,sigma=sigma)
#' itemstats(dataset)
#'
#' specific <- "S1 = 1-5
#'              S2 = 6-11
#'              S3 = 12-16"
#' model <- '
#'     G1 = 1-8
#'     G2 = 9-16
#'     COV = G1*G2'
#'
#' # quadpts dropped for faster estimation, but not as precise
#' simmod <- bfactor(dataset, specific, model, quadpts = 9, TOL = 1e-3)
#' coef(simmod, simplify=TRUE)
#' summary(simmod)
#' itemfit(simmod, QMC=TRUE)
#' M2(simmod, QMC=TRUE)
#' residuals(simmod, QMC=TRUE)
#'
#' }
#'
bfactor <- function(data, model, model2 = paste0('G = 1-', ncol(data)),
                    group = NULL, quadpts = NULL, invariance = '', ...)
{
    Call <- match.call()
    dots <- list(...)
    if(!is.null(dots$dentype))
        stop('bfactor does not currently support changing the dentype input', call.=FALSE)
    if(!is.null(dots$method))
        stop('method cannot be changed for bifactor models', call.=FALSE)
    if(!is.null(dots$formula))
        stop('bfactor does not currently support latent regression models', call.=FALSE) #TODO
    if(missing(model)) missingMsg('model')
    if(is.character(model))
        model <- mirt.model(model, itemnames=colnames(data))
    if(is(model, 'mirt.model')){
        tmp <- rep(NA, ncol(data))
        toparse <- model$x
        toparse <- toparse[!(toparse[,1] %in% mirt.model_keywords()), 2]
        toparse <- replace_dash(toparse)
        loads <- lapply(toparse, \(x) as.integer(strsplit(x, split=',')[[1L]]))
        stopifnot("bifactor model specific loadings not unique; please fix" =
                      all(table(do.call(c, loads)) == 1L))
        for(i in 1L:length(loads)) tmp[loads[[i]]] <- i
        model <- tmp
    }
    if(!is.numeric(model))
        stop('model must be a numeric vector', call.=FALSE)
    if(is.numeric(model))
        if(length(model) != ncol(data))
            stop('length of model must equal the number of items', call.=FALSE)
    uniq_vals <- sort(na.omit(unique(model)))
    specific <- model
    for(i in seq_len(length(uniq_vals)))
        specific[uniq_vals[i] == model & !is.na(model)] <- i
    nspec <- length(uniq_vals)
    if(is.character(model2)){
        tmp <- any(sapply(colnames(data), grepl, x=model2))
        model2 <- mirt.model(model2, itemnames = if(tmp) colnames(data) else NULL)
    }
    if(!is(model2, 'mirt.model'))
        stop('model2 must be an appropriate second-tier model', call.=FALSE)
    model <- bfactor2mod(specific, ncol(data))
    model$x <- rbind(model2$x, model$x)
    attr(model, 'nspec') <- nspec
    attr(model, 'specific') <- specific
    if(is.null(group)) group <- rep('all', nrow(data))
    mod <- ESTIMATION(data=data, model=model, group=group,
                      method='EM', quadpts=quadpts,
                      BFACTOR = TRUE, invariance=invariance, ...)
    if(is(mod, 'SingleGroupClass') || is(mod, 'MultipleGroupClass'))
        mod@Call <- Call
    return(mod)
}
