#' Multidimensional discrete item response theory
#'
#' \code{mdirt} fits a variety of item response models with discrete latent variables.
#' Posterior classification accuracy for each response pattern may be obtained
#' via the \code{\link{fscores}} function. The \code{summary()} function will display
#' the category probability values given the class membership, which can also
#' be displayed graphically with \code{plot()}, while \code{coef()}
#' displays the raw coefficient values (and their standard errors, if estimated). Finally,
#' \code{anova()} is used to compare nested models.
#'
#' @section 'lca' model definition:
#'
#' The latent class IRT model with two latent classes has the form
#'
#' \deqn{P(x = k|\theta_1, \theta_2, a1, a2) = \frac{exp(s_k (a1 \theta_1 + a2 \theta_2))}{
#'   \sum_j^K exp(s_j (a1 \theta_1 + a2 \theta_2))}}
#'
#' where the \eqn{\theta} values generally take on discrete points (such as 0 or 1), and
#' the \eqn{s_k}'s are the scoring values for each category. If the model is selected to be
#' \code{'lca'} then the \eqn{s_k} values are fixed to \code{s_k = 0:(ncat - 1)}, whereas if
#' the model is \code{'nlca'} the \eqn{s_k} are all fixed to 1. For proper identification, the
#' first category slope parameters (\eqn{a1} and \eqn{a2}) are never freely estimated.
#'
#' @param data a \code{matrix} or \code{data.frame} that consists of
#'   numerically ordered data, with missing data coded as \code{NA}
#' @param model number of classes to fit, or alternatively a \code{\link{mirt.model}} definition
#' @param itemtype item types to use. Can be the \code{'lca'} model for defining ordinal
#'   item response models (dichotomous items are a special case), \code{'nlca'} for the
#'   unordered latent class model
#' @param group a factor variable indicating group membership used for multiple group analyses
#' @param GenRandomPars logical; use random starting values
#' @param technical technical input list, most interesting for discrete latent models
#'   by building a \code{customTheta} input. The default builds the integration grid for the
#'   latent class model with \code{customTheta = diag(nclasses)}; see \code{\link{mirt}} for
#'   further details
#' @param nruns a numeric value indicating how many times the model should be fit to the data
#'   when using random starting values. If greater than 1, \code{GenRandomPars} is set to true
#'   by default
#' @param return_max logical; when \code{nruns > 1}, return the model that has the most optimal
#'   maximum likelihood criteria? If FALSE, returns a list of all the estimated objects
#' @param pars used for modifying starting values; see \code{\link{mirt}} for details
#' @param verbose logical; turn on messages to the R console
#' @param ... additional arguments to be passed to the estimation engine. See \code{\link{mirt}}
#'   for more details and examples
#'
#' @seealso \code{\link{fscores}}, \code{\link{mirt.model}}, \code{\link{M2}},
#'   \code{\link{itemfit}}, \code{\link{boot.mirt}}, \code{\link{mirtCluster}},
#'   \code{\link{wald}}, \code{\link{coef-method}}, \code{\link{summary-method}},
#'   \code{\link{anova-method}}, \code{\link{residuals-method}}
#' @keywords models
#' @export mdirt
#' @examples
#'
#' \dontrun{
#' #LSAT6 dataset
#' dat <- expand.table(LSAT6)
#'
#' # fit with 2-3 latent classes
#' (mod2 <- mdirt(dat, 2))
#' (mod3 <- mdirt(dat, 3))
#' summary(mod2)
#' residuals(mod2)
#' residuals(mod2, type = 'exp')
#' anova(mod2, mod3)
#' M2(mod2)
#' itemfit(mod2)
#'
#' #generate classification plots
#' plot(mod2)
#' plot(mod2, facet_items = FALSE)
#' plot(mod2, profile = TRUE)
#'
#' # available for polytomous data
#' mod <- mdirt(Science, 2)
#' summary(mod)
#' plot(mod)
#' plot(mod, profile=TRUE)
#'
#' # classification based on response patterns
#' fscores(mod2, full.scores = FALSE)
#'
#' # classify individuals either with the largest posterior probability.....
#' fs <- fscores(mod2)
#' head(fs)
#' classes <- matrix(1:2, nrow(fs), 2, byrow=TRUE)
#' class_max <- classes[t(apply(fs, 1, max) == fs)]
#' table(class_max)
#'
#' # ... or by probability sampling (closer to estimated class proportions)
#' class_prob <- apply(fs, 1, function(x) sample(1:2, 1, prob=x))
#' table(class_prob)
#'
#' # fit with random starting points (run in parallel to save time)
#' mirtCluster()
#' mod <- mdirt(dat, 2, nruns=10)
#'
#' #--------------------------
#' # Grade of measurement model
#'
#' # define a custom Theta grid for including a 'fuzzy' class membership
#' (Theta <- matrix(c(1, 0, .5, .5, 0, 1), nrow=3 , ncol=2, byrow=TRUE))
#' (mod_gom <- mdirt(dat, 2, technical = list(customTheta = Theta)))
#' summary(mod_gom)
#'
#' #-----------------
#' # Multidimensional discrete model
#'
#' dat <- key2binary(SAT12,
#'      key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#'
#' # define Theta grid for three latent classes
#' (Theta <- matrix(c(0,0,0, 1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,0,1, 0,1,1, 1,1,1),
#'    ncol=3, byrow=TRUE))
#' (mod_discrete <- mdirt(dat, 3, technical = list(customTheta = Theta)))
#' summary(mod_discrete)
#'
#' }
mdirt <- function(data, model, itemtype = 'lca', nruns = 1,
                  return_max = TRUE, group = NULL, GenRandomPars = FALSE,
                  verbose = TRUE, pars = NULL, technical = list(), ...)
{
    Call <- match.call()
    if(!all(itemtype %in% c('lca', 'nlca')))
        stop('Selected itemtype not supported. Please use itemtype \'lca\' or \'nlca\'', call.=FALSE)
    if(nruns > 1) GenRandomPars <- TRUE
    if(is.null(group)) group <- rep('all', nrow(data))
    mods <- myLapply(1:nruns, function(x, ...) return(ESTIMATION(...)),
                     data=data, model=model, group=group, itemtype=itemtype, method='EM',
                     technical=technical, calcNull=FALSE, GenRandomPars=GenRandomPars,
                     discrete=TRUE, verbose=ifelse(nruns > 1L, FALSE, verbose), pars=pars, ...)
    if(is(mods[[1L]], 'DiscreteClass')){
        for(i in 1:length(mods)) mods[[i]]@Call <- Call
    }
    if(!return_max){
        return(mods)
    } else {
        if(is(mods[[1L]], 'DiscreteClass')){
            LL <- sapply(mods, function(x) x@logLik)
            if(verbose && nruns > 1L){
                cat('Model log-likelihoods:\n')
                print(round(LL, 4))
            }
            mods <- mods[[which(max(LL) == LL)[1L]]]
        }
    }
    if(!is.null(pars) && pars == 'values') mods <- mods[[1L]]
    return(mods)
}
