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
#' \deqn{P(x = k|\theta_1, \theta_2, a1, a2) = \frac{exp(a1 \theta_1 + a2 \theta_2)}{
#'   \sum_j^K exp(a1 \theta_1 + a2 \theta_2)}}
#'
#' where the \eqn{\theta} values generally take on discrete points (such as 0 or 1).
#' For proper identification, the first category slope parameters
#' (\eqn{a1} and \eqn{a2}) are never freely estimated. Alternatively, supplying a different
#' grid of \eqn{\theta} values will allow the estimation of similar models (multidimensional
#' discrete models, grade of membership, etc.). See the examples below.
#'
#' @param data a \code{matrix} or \code{data.frame} that consists of
#'   numerically ordered data, with missing data coded as \code{NA}
#' @param model number of classes to fit, or alternatively a \code{\link{mirt.model}} definition
#' @param method estimation method. Can be 'EM' or 'BL' (see \code{\link{mirt}} for more details)
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
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
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
#' # generate classification plots
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
#' # plausible value imputations for stocastic classification in both classes
#' pvs <- fscores(mod2, plausible.draws=10)
#' tabs <- lapply(pvs, function(x) apply(x, 2, table))
#' tabs[[1]]
#'
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
#' # Multidimensional discrete latent class model
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
#' # Located latent class model
#' model <- mirt.model('C1 = 1-32
#'                      C2 = 1-32
#'                      C3 = 1-32
#'                      CONSTRAIN = (1-32, a1), (1-32, a2), (1-32, a3)')
#' (mod_located <- mdirt(dat, model,
#'                      technical = list(customTheta = diag(3))))
#' summary(mod_located)
#'
#' #-----------------
#' ### DINA model example
#' # generate some suitable data for a two dimensional DINA application
#' #     (first columns are intercepts)
#' set.seed(1)
#' Theta <- expand.table(matrix(c(1,0,0,0, 200,
#'                                1,1,0,0, 200,
#'                                1,0,1,0, 100,
#'                                1,1,1,1, 500), 4, 5, byrow=TRUE))
#' a <- matrix(c(rnorm(15, -1.5, .5), rlnorm(5, .2, .3), numeric(15), rlnorm(5, .2, .3),
#'               numeric(15), rlnorm(5, .2, .3)), 15, 4)
#'
#' guess <- plogis(a[11:15,1]) # population guess
#' slip <- 1 - plogis(rowSums(a[11:15,])) # population slip
#'
#' dat <- simdata(a, Theta=Theta, itemtype = 'lca')
#'
#' # first column is the intercept, 2nd and 3rd are attributes
#' theta <- matrix(c(1,0,0,
#'                   1,1,0,
#'                   1,0,1,
#'                   1,1,1), 4, 3, byrow=TRUE)
#' theta <- cbind(theta, theta[,2] * theta[,3]) #DINA interaction of main attributes
#' model <- mirt.model('Intercept = 1-15
#'                      A1 = 1-5
#'                      A2 = 6-10
#'                      A1A2 = 11-15')
#'
#' mod <- mdirt(dat, model, technical = list(customTheta = theta))
#' coef(mod)
#' summary(mod)
#' M2(mod) # fits well
#'
#' cfs <- coef(mod, simplify=TRUE)$items[11:15,]
#' cbind(guess, estguess = plogis(cfs[,1]))
#' cbind(slip, estslip = 1 - plogis(rowSums(cfs)))
#'
#'
#' ### DINO model example
#' theta <- matrix(c(1,0,0,
#'                   1,1,0,
#'                   1,0,1,
#'                   1,1,1), 4, 3, byrow=TRUE)
#' # define theta matrix with negative interaction term
#' theta <- cbind(theta, -theta[,2] * theta[,3])
#'
#' model <- mirt.model('Intercept = 1-15
#'                      A1 = 1-5, 11-15
#'                      A2 = 6-15
#'                      Yoshi = 11-15
#'                      CONSTRAIN = (11,a2,a3,a4), (12,a2,a3,a4), (13,a2,a3,a4),
#'                                  (14,a2,a3,a4), (15,a2,a3,a4)')
#'
#' mod <- mdirt(dat, model, technical = list(customTheta = theta))
#' coef(mod, simplify=TRUE)
#' summary(mod)
#' M2(mod) #doesn't fit as well, because not the generating model
#'
#'
#' #------------------
#' #multidimensional latent class model
#'
#' dat <- key2binary(SAT12,
#'      key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#'
#' # 5 latent classes within 2 different sets of items
#' model <- mirt.model('C1 = 1-16
#'                      C2 = 1-16
#'                      C3 = 1-16
#'                      C4 = 1-16
#'                      C5 = 1-16
#'                      C6 = 17-32
#'                      C7 = 17-32
#'                      C8 = 17-32
#'                      C9 = 17-32
#'                      C10 = 17-32
#'                      CONSTRAIN = (1-16, a1), (1-16, a2), (1-16, a3), (1-16, a4), (1-16, a5),
#'                        (17-32, a6), (17-32, a7), (17-32, a8), (17-32, a9), (17-32, a10)')
#'
#' theta <- diag(10)
#' mod <- mdirt(dat, model, technical = list(customTheta = theta))
#' coef(mod, simplify=TRUE)
#' summary(mod)
#'
#' #------------------
#' # multiple group with constrained group probabilities
#'  dat <- key2binary(SAT12,
#'    key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' group <- rep(c('G1', 'G2'), each = nrow(SAT12)/2)
#' Theta <- diag(2)
#'
#' # the latent class parameters are technically located in the (nitems + 1) location
#' model <- mirt.model('A1 = 1-32
#'                      A2 = 1-32
#'                      CONSTRAINB = (33, c1)')
#' mod <- mdirt(dat, model, group = group,
#'              technical = list(customTheta = Theta))
#' coef(mod, simplify=TRUE)
#' summary(mod)
#'
#'
#' }
mdirt <- function(data, model, nruns = 1, method = 'EM',
                  return_max = TRUE, group = NULL, GenRandomPars = FALSE,
                  verbose = TRUE, pars = NULL, technical = list(), ...)
{
    Call <- match.call()
    itemtype <- 'lca'
    stopifnot(method %in% c('EM', 'BL'))
    if(nruns > 1) GenRandomPars <- TRUE
    if(is.null(group)) group <- rep('all', nrow(data))
    if((is(model, 'mirt.model') || is.character(model)) && is.null(technical$customTheta))
        stop('customTheta input required when using a mirt.model type input')
    mods <- myLapply(1:nruns, function(x, ...) return(ESTIMATION(...)), method=method,
                     data=data, model=model, group=group, itemtype=itemtype,
                     technical=technical, calcNull=FALSE, GenRandomPars=GenRandomPars,
                     discrete=TRUE, verbose=ifelse(nruns > 1L, FALSE, verbose), pars=pars, ...)
    if(is(mods[[1L]], 'DiscreteClass')){
        for(i in 1:length(mods)) mods[[i]]@Call <- Call
    }
    if(!return_max){
        return(mods)
    } else {
        if(is(mods[[1L]], 'DiscreteClass')){
            LL <- sapply(mods, function(x) x@Fit$logLik)
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
