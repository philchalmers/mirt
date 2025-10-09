#' Multidimensional discrete item response theory
#'
#' \code{mdirt} fits a variety of item response models with discrete latent variables.
#' These include, but are not limited to, latent class analysis, multidimensional latent
#' class models, multidimensional discrete latent class models, DINA/DINO models,
#' grade of measurement models, C-RUM, and so on. If response models are not defined explicitly
#' then customized models can defined using the \code{\link{createItem}} function.
#'
#' Posterior classification accuracy for each response pattern may be obtained
#' via the \code{\link{fscores}} function. The \code{summary()} function will display
#' the category probability values given the class membership, which can also
#' be displayed graphically with \code{plot()}, while \code{coef()}
#' displays the raw coefficient values (and their standard errors, if estimated). Finally,
#' \code{anova()} is used to compare nested models, while
#' \code{\link{M2}} and \code{\link{itemfit}} may be used for model fitting purposes.
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
#' When the \code{item.Q} for is utilized, the above equation can be understood as
#'
#' \deqn{P(x = k|\theta_1, \theta_2, a1, a2) = \frac{exp(a1 \theta_1 Q_{j1} + a2 \theta_2 Q_{j2})}{
#'   \sum_j^K exp(a1 \theta_1 Q_{j1} + a2 \theta_2 Q_{j2})}}
#'
#' where by construction \code{Q} is a \eqn{K_i \times A} matrix indicating whether the category should
#' be modeled according to the latent class structure. For the standard latent class model, the Q-matrix
#' has as many rows as categories, as many columns as the number of classes/attributes modeled,
#' and consist of 0's in the first row and 1's elsewhere. This of course can be over-written by passing
#' an alternative \code{item.Q} definition for each respective item.
#'
#' @param data a \code{matrix} or \code{data.frame} that consists of
#'   numerically ordered data, organized in the form of integers,
#'    with missing data coded as \code{NA}
#' @param model number of mutually exclusive classes to fit, or alternatively a more specific
#'   \code{\link{mirt.model}} definition (which reflects the so-called Q-matrix).
#'   Note that when using a \code{\link{mirt.model}},
#'   the order with which the syntax factors/attributes are defined are associated with the
#'   columns in the \code{customTheta} input
#' @param method estimation method. Can be 'EM' or 'BL' (see \code{\link{mirt}} for more details)
#' @param optimizer optimizer used for the M-step, set to \code{'nlminb'} by default.
#'   See \code{\link{mirt}} for more details
#' @param group a factor variable indicating group membership used for multiple group analyses
#' @param itemtype a vector indicating the itemtype associated with each item.
#'   For discrete models this is limited to only 'lca' or items defined using a
#'   \code{\link{createItem}} definition
#' @param item.Q a list of item-level Q-matrices indicating how the respective categories should be
#'   modeled by the underlying attributes. Each matrix must represent a \eqn{K_i \times A} matrix,
#'   where \eqn{K_i} represents the number of categories for the ith item, and \eqn{A} is the number
#'   of attributes included in the \code{Theta} matrix; otherwise, a value of\code{NULL} will default
#'   to a matrix consisting of 1's for each \eqn{K_i \times A} element except for the first row, which
#'   contains only 0's for proper identification. Incidentally, the first row of each matrix \code{must}
#'   contain only 0's so that the first category represents the reference category for identification
#' @param GenRandomPars logical; use random starting values
#' @param customTheta input passed to \code{technical = list(customTheta = ...)}, but is included
#'   directly in this function for convenience. This input is most interesting for discrete latent models
#'   because it allows customized patterns of latent classes (i.e., defines the possible combinations
#'   of the latent attribute profile). The default builds the pattern \code{customTheta = diag(model)},
#'   which is the typical pattern for the traditional latent class analysis whereby class
#'   membership mutually distinct and exhaustive. See \code{\link{thetaComb}} for a quick method
#'   to generate a matrix with all possible combinations
#' @param nruns a numeric value indicating how many times the model should be fit to the data
#'   when using random starting values. If greater than 1, \code{GenRandomPars} is set to \code{TRUE}
#'   by default. Using this returns a list of fitted model objects, where the model
#'   with the highest log-likelihood should generally be selected as the model
#'   best associated with the MLE. Note that if a \code{\link{mirtCluster}} was
#'   defined earlier then the runs will be run in parallel
#' @param return_max logical; when \code{nruns > 1}, return the model that has the most optimal
#'   maximum likelihood criteria? If FALSE, returns a list of all the estimated objects
#' @param covdata a data.frame of data used for latent regression models
#' @param formula an R formula (or list of formulas) indicating how the latent traits
#'   can be regressed using external covariates in \code{covdata}. If a named list
#'   of formulas is supplied (where the names correspond to the latent trait/attribute names in \code{model})
#'   then specific regression effects can be estimated for each factor. Supplying a single formula
#'   will estimate the regression parameters for all latent variables by default
#' @param structure an R formula allowing the profile probability patterns (i.e., the structural component of
#'   the model) to be fitted according to a log-linear model. When \code{NULL}, all profile probabilities
#'   (except one) will be estimated. Use of this input requires that the \code{customTheta} input is supplied,
#'   and that the column names in this matrix match the names found within this formula
#' @param pars used for modifying starting values; see \code{\link{mirt}} for details
#' @param verbose logical; turn on messages to the R console
#' @param technical list of lower-level inputs. See \code{\link{mirt}} for details
#' @param ... additional arguments to be passed to the estimation engine. See \code{\link{mirt}}
#'   for more details and examples
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#'
#' Proctor, C. H. (1970). A probabilistic formulation and statistical analysis for Guttman scaling.
#'   \emph{Psychometrika, 35}, 73-78.
#' \doi{10.18637/jss.v048.i06}
#' @seealso \code{\link{thetaComb}}, \code{\link{fscores}}, \code{\link{mirt.model}}, \code{\link{M2}},
#'   \code{\link{itemfit}}, \code{\link{boot.mirt}}, \code{\link{mirtCluster}},
#'   \code{\link{wald}}, \code{\link{coef-method}}, \code{\link{summary-method}},
#'   \code{\link{anova-method}}, \code{\link{residuals-method}}
#' @keywords models
#' @export mdirt
#' @examples
#'
#' # LSAT6 dataset
#' dat <- expand.table(LSAT6)
#'
#' # fit with 2-3 latent classes
#' (mod2 <- mdirt(dat, 2))
#' \donttest{
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
#' classes <- 1:2
#' class_max <- classes[apply(apply(fs, 1, max) == fs, 1, which)]
#' table(class_max)
#'
#' # ... or by probability sampling (i.e., plausible value draws)
#' class_prob <- apply(fs, 1, function(x) sample(1:2, 1, prob=x))
#' table(class_prob)
#'
#' # plausible value imputations for stochastic classification in both classes
#' pvs <- fscores(mod2, plausible.draws=10)
#' tabs <- lapply(pvs, function(x) apply(x, 2, table))
#' tabs[[1]]
#'
#'
#' # fit with random starting points (run in parallel to save time)
#' if(interactive()) mirtCluster()
#' mod <- mdirt(dat, 2, nruns=10)
#'
#' #--------------------------
#' # Grade of measurement model
#'
#' # define a custom Theta grid for including a 'fuzzy' class membership
#' (Theta <- matrix(c(1, 0, .5, .5, 0, 1), nrow=3 , ncol=2, byrow=TRUE))
#' (mod_gom <- mdirt(dat, 2, customTheta = Theta))
#' summary(mod_gom)
#'
#' #-----------------
#' # Multidimensional discrete latent class model
#'
#' dat <- key2binary(SAT12,
#'      key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#'
#' # define Theta grid for three latent classes
#' (Theta <- thetaComb(0:1, 3))
#' (mod_discrete <- mdirt(dat, 3, customTheta = Theta))
#' summary(mod_discrete)
#'
#' # Located latent class model
#' model <- mirt.model('C1 = 1-32
#'                      C2 = 1-32
#'                      C3 = 1-32
#'                      CONSTRAIN = (1-32, a1), (1-32, a2), (1-32, a3)')
#' (mod_located <- mdirt(dat, model, customTheta = diag(3)))
#' summary(mod_located)
#'
#' #-----------------
#' ### DINA model example
#' # generate some suitable data for a two dimensional DINA application
#' #     (first columns are intercepts)
#' set.seed(1)
#' Theta <- expand.table(matrix(c(1,0,0,0,
#'                                1,1,0,0,
#'                                1,0,1,0,
#'                                1,1,1,1), 4, 4, byrow=TRUE),
#'                       freq = c(200,200,100,500))
#' a <- matrix(c(rnorm(15, -1.5, .5), rlnorm(5, .2, .3), numeric(15), rlnorm(5, .2, .3),
#'               numeric(15), rlnorm(5, .2, .3)), 15, 4)
#'
#' guess <- plogis(a[11:15,1]) # population guess
#' slip <- 1 - plogis(rowSums(a[11:15,])) # population slip
#'
#' dat <- simdata(a, Theta=Theta, itemtype = 'lca')
#'
#' # first column is the intercept, 2nd and 3rd are attributes
#' theta <- cbind(1, thetaComb(0:1, 2))
#' theta <- cbind(theta, theta[,2] * theta[,3]) #DINA interaction of main attributes
#' model <- mirt.model('Intercept = 1-15
#'                      A1 = 1-5
#'                      A2 = 6-10
#'                      A1A2 = 11-15')
#'
#' # last 5 items are DINA (first 10 are unidimensional C-RUMs)
#' DINA <- mdirt(dat, model, customTheta = theta)
#' coef(DINA, simplify=TRUE)
#' summary(DINA)
#' M2(DINA) # fits well (as it should)
#'
#' cfs <- coef(DINA, simplify=TRUE)$items[11:15,]
#' cbind(guess, estguess = plogis(cfs[,1]))
#' cbind(slip, estslip = 1 - plogis(rowSums(cfs)))
#'
#'
#' ### DINO model example
#' theta <- cbind(1, thetaComb(0:1, 2))
#' # define theta matrix with negative interaction term
#' (theta <- cbind(theta, -theta[,2] * theta[,3]))
#'
#' model <- mirt.model('Intercept = 1-15
#'                      A1 = 1-5, 11-15
#'                      A2 = 6-15
#'                      Yoshi = 11-15
#'                      CONSTRAIN = (11,a2,a3,a4), (12,a2,a3,a4), (13,a2,a3,a4),
#'                                  (14,a2,a3,a4), (15,a2,a3,a4)')
#'
#' # last five items are DINOs (first 10 are unidimensional C-RUMs)
#' DINO <- mdirt(dat, model, customTheta = theta)
#' coef(DINO, simplify=TRUE)
#' summary(DINO)
#' M2(DINO) #doesn't fit as well, because not the generating model
#'
#' ## C-RUM (analogous to MIRT model)
#' theta <- cbind(1, thetaComb(0:1, 2))
#' model <- mirt.model('Intercept = 1-15
#'                      A1 = 1-5, 11-15
#'                      A2 = 6-15')
#'
#' CRUM <- mdirt(dat, model, customTheta = theta)
#' coef(CRUM, simplify=TRUE)
#' summary(CRUM)
#'
#' # good fit, but over-saturated (main effects for items 11-15 can be set to 0)
#' M2(CRUM)
#'
#' #------------------
#' # multidimensional latent class model
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
#' theta <- diag(10) # defined explicitly. Otherwise, this profile is assumed
#' mod <- mdirt(dat, model, customTheta = theta)
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
#' mod <- mdirt(dat, model, group = group, customTheta = Theta)
#' coef(mod, simplify=TRUE)
#' summary(mod)
#'
#'
#' #------------------
#' # Probabilistic Guttman Model (Proctor, 1970)
#'
#' # example analysis can also be found in the sirt package (see ?prob.guttman)
#' data(data.read, package = 'sirt')
#' head(data.read)
#'
#' Theta <- matrix(c(1,0,0,0,
#'                   1,1,0,0,
#'                   1,1,1,0,
#'                   1,1,1,1), 4, byrow=TRUE)
#'
#' model <- mirt.model("INTERCEPT = 1-12
#'                      C1 = 1,7,9,11
#'                      C2 = 2,5,8,10,12
#'                      C3 = 3,4,6")
#'
#' mod <- mdirt(data.read, model, customTheta=Theta)
#' summary(mod)
#'
#' M2(mod)
#' itemfit(mod)
#'
#'
#' }
mdirt <- function(data, model, customTheta = NULL, structure = NULL, item.Q = NULL,
                  nruns = 1, method = 'EM', covdata = NULL, formula = NULL, itemtype = 'lca',
                  optimizer = 'nlminb', return_max = TRUE, group = NULL, GenRandomPars = FALSE,
                  verbose = TRUE, pars = NULL, technical = list(), ...)
{
    Call <- match.call()
    dots <- list(...)
    if(!is.null(dots$dentype))
        stop('mdirt does not support changing the dentype input', call.=FALSE)
    latent.regression <- latentRegression_obj(data=data, covdata=covdata, formula=formula,
                                              dentype = 'discrete', method=method)
    technical$customTheta <- customTheta
    valid_itemtype <- 'lca'
    if(!is.null(dots$customItems))
        valid_itemtype <- c(valid_itemtype, names(dots$customItems))
    stopifnot(all(itemtype %in% valid_itemtype))
    stopifnot(method %in% c('EM', 'BL'))
    if(nruns > 1) GenRandomPars <- TRUE
    if(is.null(group)) group <- rep('all', nrow(data))
    if((is(model, 'mirt.model') || is.character(model)) && is.null(technical$customTheta))
        stop('customTheta input required when using a mirt.model type input')
    technical$omp <- FALSE
    mods <- myLapply(1:nruns, function(x, ...) return(ESTIMATION(...)), progress=verbose && nruns > 1L,
                     method=method, latent.regression=latent.regression, structure=structure,
                     data=data, model=model, group=group, itemtype=itemtype, optimizer=optimizer,
                     technical=technical, calcNull=FALSE, GenRandomPars=GenRandomPars, item.Q=item.Q,
                     dentype = 'discrete', verbose=ifelse(nruns > 1L, FALSE, verbose), pars=pars, ...)
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
    if(!is.null(pars) && !is.data.frame(pars) && pars == 'values') mods <- mods[[1L]]
    return(mods)
}
