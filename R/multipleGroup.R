#' Multiple Group Estimation
#'
#' \code{multipleGroup} performs a full-information
#' maximum-likelihood multiple group analysis for dichotomous and polytomous
#' data under the item response theory paradigm using either Cai's (2010)
#' Metropolis-Hastings Robbins-Monro (MHRM) algorithm or with an EM algorithm approach. This function
#' may be used for detecting differential item functioning (DIF), thought the \code{\link{DIF}} function
#' may provide a more convenient approach.
#'
#' By default the estimation in \code{multipleGroup} assumes that the models are maximally
#' independent, and therefore could initially be performed by sub-setting the data and running identical
#' models with \code{mirt} and aggregating the results (e.g., log-likelihood).
#' However, constrains may be imposed across groups by invoking various \code{invariance} keywords and
#' \code{constrain = ...} arguments, by inputting user specified design matrix from
#' \code{mod2values} or from passing \code{pars = 'values'}, or by supplying a \code{constrain} list
#' for user defined equality constraints between parameters.
#'
#' @aliases multipleGroup
#' @param data a \code{matrix} or \code{data.frame} that consists of
#'   numerically ordered data, with missing data coded as \code{NA}
#' @param model a single model object returned from \code{mirt.model()} declaring how
#'   the factor model is to be estimated. See \code{\link{mirt.model}} for more details
#' @param group a character vector indicating group membership
#' @param invariance a character vector containing the following possible options:
#'   \describe{
#'     \item{\code{'free_means'}}{for freely estimating all latent means (reference group constrained to 0)}
#'     \item{\code{'free_var'}}{for freely estimating all latent variances (reference group constrained to 1's)}
#'     \item{\code{'free_cov'}}{for freely estimating all latent covariances (reference group constrained to an
#'     Identity matrix)}
#'     \item{\code{'free_varcov'}}{calls both \code{'free_var'} and \code{'free_cov'}}
#'     \item{\code{'slopes'}}{to constrain all the slopes to be equal across all groups}
#'     \item{\code{'intercepts'}}{to constrain all the intercepts to be equal across all groups, note for
#'     nominal models this also includes the category specific slope parameters}
#'    }
#'
#'   Additionally, specifying specific item name bundles (from \code{colnames(data)}) will
#'   constrain all freely estimated parameters in each item to be equal across groups. This is useful
#'   for selecting 'anchor' items for vertical and horizontal scaling, and for detecting differential item
#'   functioning (DIF) across groups
#' @param guess initial (or fixed) values for the pseudo-guessing parameter. Can be
#'   entered as a single value to assign a global guessing parameter or may be entered as
#'   a numeric vector for each item
#' @param upper initial (or fixed) upper bound parameters for 4-PL model. Can be
#'   entered as a single value to assign a global upper bound parameter or may be entered as a
#'   numeric vector corresponding to each item
#' @param accelerate see \code{\link{mirt}} for more details
#' @param SE logical; estimate the information matrix for standard errors?
#' @param SE.type see \code{\link{mirt}} for more details
#' @param verbose logical; display iteration history during estimation?
#' @param draws the number of Monte Carlo draws to estimate the log-likelihood
#' @param quadpts the number of quadratures to be used per dimensions when \code{method = 'EM'}
#' @param calcNull logical; calculate the Null model for fit statics (e.g., TLI)?
#' @param method a character indicating whether to use the EM (\code{'EM'}) or the MH-RM
#'   (\code{'MHRM'}) algorithm
#' @param type type of plot to view; can be \code{'info'} to show the test
#'   information function, \code{'infocontour'} for the test information contours,
#'   \code{'SE'} for the test standard error function, \code{'RE'} for the relative efficiency plot,
#'   and \code{'score'} for the expected total score plot
#' @param empiricalhist logical; estimate prior distribution using an empirical histogram approach.
#'   see \code{mirt} for details
#' @param GenRandomPars logical; generate random starting values prior to optimization instead of
#'          using the fixed internal starting values?
#' @param key see \code{\link{mirt}} for details
#' @param itemtype see \code{\link{mirt}} for details
#' @param constrain see \code{\link{mirt}} for details
#' @param grsm.block see \code{\link{mirt}} for details
# @param rsm.block see \code{\link{mirt}} for details
#' @param parprior see \code{\link{mirt}} for details
#' @param pars see \code{\link{mirt}} for details
#' @param TOL see \code{\link{mirt}} for details
#' @param ... additional arguments to be passed
#' @param technical list specifying subtle parameters that can be adjusted. See
#'   \code{\link{mirt}} for details
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{anova-method}}, \code{\link{coef-method}}, \code{\link{summary-method}},
#'   \code{\link{residuals-method}}, \code{\link{plot-method}},
#'   \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{mirt.model}}, \code{\link{mirt}},
#'   \code{\link{bfactor}}, \code{\link{multipleGroup}}, \code{\link{mixedmirt}},
#'   \code{\link{wald}}, \code{\link{itemplot}}, \code{\link{fscores}}, \code{\link{M2}},
#'   \code{\link{extract.item}}, \code{\link{iteminfo}}, \code{\link{testinfo}}, \code{\link{probtrace}},
#'   \code{\link{boot.mirt}}, \code{\link{imputeMissing}}, \code{\link{itemfit}}, \code{\link{mod2values}},
#'   \code{\link{simdata}}, \code{\link{createItem}}, \code{\link{mirtCluster}}, \code{\link{DIF}}
#   \code{\link{DTF}}
#' @keywords models
#' @export multipleGroup
#' @examples
#' \dontrun{
#' #single factor
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('dich', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#' models <- mirt.model('F1 = 1-15')
#'
#' mod_configural <- multipleGroup(dat, models, group = group) #completely separate analyses
#' #limited information fit statistics
#' M2(mod_configural)
#'
#' mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes')) #equal slopes
#' #equal intercepts, free variance and means
#' mod_scalar2 <- multipleGroup(dat, models, group = group,
#'                              invariance=c('slopes', 'intercepts', 'free_var','free_means'))
#' mod_scalar1 <- multipleGroup(dat, models, group = group,  #fixed means
#'                              invariance=c('slopes', 'intercepts', 'free_var'))
#' mod_fullconstrain <- multipleGroup(dat, models, group = group,
#'                              invariance=c('slopes', 'intercepts'))
#'
#' summary(mod_scalar2)
#' coef(mod_scalar2)
#' residuals(mod_scalar2)
#' plot(mod_configural)
#' plot(mod_configural, type = 'score')
#' plot(mod_configural, type = 'trace')
#' plot(mod_configural, type = 'trace', which.items = 1:4)
#' itemplot(mod_configural, 2)
#' itemplot(mod_configural, 2, type = 'RE')
#'
#' anova(mod_metric, mod_configural) #equal slopes only
#' anova(mod_scalar2, mod_metric) #equal intercepts, free variance and mean
#' anova(mod_scalar1, mod_scalar2) #fix mean
#' anova(mod_fullconstrain, mod_scalar1) #fix variance
#'
#'
#' #test whether first 6 slopes should be equal across groups
#' values <- multipleGroup(dat, models, group = group, pars = 'values')
#' values
#' constrain <- list(c(1, 63), c(5,67), c(9,71), c(13,75), c(17,79), c(21,83))
#' equalslopes <- multipleGroup(dat, models, group = group, constrain = constrain)
#' anova(equalslopes, mod_configural)
#'
#' #same as above, but using mirt.model syntax
#' newmodel <- mirt.model('
#'     F = 1-15
#'     CONSTRAINB = (1-6, a1)')
#' equalslopes <- multipleGroup(dat, newmodel, group = group)
#' coef(equalslopes)
#'
#' #############
#' #DIF test for each item (using all other items as anchors)
#' itemnames <- colnames(dat)
#' refmodel <- multipleGroup(dat, models, group = group, SE=TRUE,
#'                              invariance=c('free_means', 'free_varcov', itemnames))
#'
#' #loop over items (in practice, run in parallel to increase speed). May be better to use ?DIF
#' estmodels <- vector('list', ncol(dat))
#' for(i in 1:ncol(dat))
#'     estmodels[[i]] <- multipleGroup(dat, models, group = group, verbose = FALSE, calcNull=FALSE,
#'                              invariance=c('free_means', 'free_varcov', itemnames[-i]))
#'
#' (anovas <- lapply(estmodels, anova, object2=refmodel, verbose=FALSE))
#'
#' #family-wise error control
#' p <- do.call(rbind, lapply(anovas, function(x) x[2, 'p']))
#' p.adjust(p, method = 'BH')
#'
#' #same as above, except only test if slopes vary (1 df)
#' #constrain all intercepts
#' estmodels <- vector('list', ncol(dat))
#' for(i in 1:ncol(dat))
#'     estmodels[[i]] <- multipleGroup(dat, models, group = group, verbose = FALSE, calcNull=FALSE,
#'                              invariance=c('free_means', 'free_varcov', 'intercepts',
#'                              itemnames[-i]))
#'
#' (anovas <- lapply(estmodels, anova, object2=refmodel, verbose=FALSE))
#'
#' #quickly test with Wald test using DIF()
#' mod_configural2 <- multipleGroup(dat, models, group = group, SE=TRUE)
#' DIF(mod_configural2, which.par = c('a1', 'd'), Wald=TRUE, p.adjust = 'fdr')
#'
#' #############
#' #multiple factors
#'
#' a <- matrix(c(abs(rnorm(5,1,.3)), rep(0,15),abs(rnorm(5,1,.3)),
#' rep(0,15),abs(rnorm(5,1,.3))), 15, 3)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' mu <- c(-.4, -.7, .1)
#' sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)
#' itemtype <- rep('dich', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = mu, sigma = sigma)
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#'
#' #group models
#' model <- mirt.model('
#'    F1 = 1-5
#'    F2 = 6-10
#'    F3 = 11-15')
#'
#' #define mirt cluster to use parallel architecture
#' mirtCluster()
#'
#' #EM approach (not as accurate with 3 factors, but generally good for quick model comparisons)
#' mod_configural <- multipleGroup(dat, model, group = group) #completely separate analyses
#' mod_metric <- multipleGroup(dat, model, group = group, invariance=c('slopes')) #equal slopes
#' mod_fullconstrain <- multipleGroup(dat, model, group = group, #equal means, slopes, intercepts
#'                              invariance=c('slopes', 'intercepts'))
#'
#' anova(mod_metric, mod_configural)
#' anova(mod_fullconstrain, mod_metric)
#'
#' #same as above, but with MHRM (generally  more accurate with 3+ factors, but slower)
#' mod_configural <- multipleGroup(dat, model, group = group, method = 'MHRM')
#' mod_metric <- multipleGroup(dat, model, group = group, invariance=c('slopes'), method = 'MHRM')
#' mod_fullconstrain <- multipleGroup(dat, model, group = group, method = 'MHRM',
#'                              invariance=c('slopes', 'intercepts'))
#'
#' anova(mod_metric, mod_configural)
#' anova(mod_fullconstrain, mod_metric)
#'
#' ############
#' #polytomous item example
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' d <- cbind(d, d-1, d-2)
#' itemtype <- rep('graded', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#' model <- mirt.model('F1 = 1-15')
#'
#' mod_configural <- multipleGroup(dat, model, group = group)
#' plot(mod_configural)
#' plot(mod_configural, type = 'SE')
#' itemplot(mod_configural, 1)
#' itemplot(mod_configural, 1, type = 'info')
#' fs <- fscores(mod_configural)
#' head(fs[["D1"]])
#' fscores(mod_configural, method = 'EAPsum')
#'
#' # constrain slopes within each group to be equal (but not across groups)
#' model2 <- mirt.model('F1 = 1-15
#'                       CONSTRAIN = (1-15, a1)')
#' mod_configural2 <- multipleGroup(dat, model2, group = group)
#' plot(mod_configural2, type = 'SE')
#' plot(mod_configural2, type = 'RE')
#' itemplot(mod_configural2, 10)
#'
#' ############
#' ## empirical histogram example (normal and bimodal groups)
#' set.seed(1234)
#' a <- matrix(rlnorm(50, .2, .2))
#' d <- matrix(rnorm(50))
#' ThetaNormal <- matrix(rnorm(2000))
#' ThetaBimodal <- scale(matrix(c(rnorm(1000, -2), rnorm(1000,2)))) #bimodal
#' Theta <- rbind(ThetaNormal, ThetaBimodal)
#' dat <- simdata(a, d, 4000, itemtype = 'dich', Theta=Theta)
#' group <- rep(c('G1', 'G2'), each=2000)
#'
#' EH <- multipleGroup(dat, 1, group=group, empiricalhist = TRUE, invariance = colnames(dat))
#' coef(EH)
#' plot(EH, type = 'empiricalhist', npts = 60)
#'
#' #dif test for item 1
#' EH1 <- multipleGroup(dat, 1, group=group, empiricalhist = TRUE, invariance = colnames(dat)[-1])
#' anova(EH, EH1)
#'
#' }
multipleGroup <- function(data, model, group, itemtype = NULL, guess = 0, upper = 1,
                          SE = FALSE, SE.type = 'crossprod', invariance = '', pars = NULL,
                          method = 'EM', constrain = NULL, parprior = NULL, calcNull = TRUE,
                          draws = 5000, quadpts = NULL, TOL = NULL, grsm.block = NULL, 
                          key = NULL, technical = list(), accelerate = TRUE, empiricalhist = FALSE,
                          GenRandomPars = FALSE, verbose = TRUE, ...)
{
    Call <- match.call()
    if(length(model) > 1L)
        stop('multipleGroup only supports single group inputs')
    invariance.check <- invariance %in% c('free_means', 'free_var', 'free_varcov')
    if(empiricalhist && any(invariance.check))
        stop('freeing group parameters not meaningful when estimating empirical histograms')
    if(sum(invariance.check == 2L) && length(constrain) == 0){
        warn <- TRUE
        if(is(model, 'mirt.model')){
            if(any(model$x[,1L] == 'CONSTRAINB'))
                warn <- FALSE
        }
        if(warn)
            stop('Model is not identified without further constrains (may require additional anchoring items).')
    }
    mod <- ESTIMATION(data=data, model=model, group=group, invariance=invariance,
                      itemtype=itemtype, guess=guess, upper=upper, empiricalhist=empiricalhist,
                      pars=pars, constrain=constrain, SE=SE, grsm.block=grsm.block,
                      parprior=parprior, quadpts=quadpts, method=method, 
                      technical = technical, verbose = verbose, calcNull=calcNull,
                      SE.type = SE.type, key=key, accelerate=accelerate, draws=draws,
                      GenRandomPars=GenRandomPars, TOL=TOL, ...)
    if(is(mod, 'MultipleGroupClass'))
        mod@Call <- Call
    return(mod)
}
