#' Mixed effects modeling for MIRT models
#'
#' \code{mixedmirt} fits MIRT models using FIML estimation to dichotomous and polytomous
#' IRT models conditional on fixed and random effect of person and item level covariates.
#' This can also be understood as 'explanatory IRT' if only fixed effects are modeled, or
#' multilevel/mixed IRT if random and fixed effects are included. The method uses the MH-RM
#' algorithm exclusively. Additionally, computation of the log-likelihood can be sped up by
#' using parallel estimation via \code{\link{mirtCluster}}.
#'
#' For dichotomous response models, \code{mixedmirt} follows the general form
#'
#'  \deqn{P(x = 1|\theta, \psi) = g + \frac{(u - g)}{1 + exp(-1 * [\theta a +
#'  X \beta + Z \delta])}}
#'
#'  where X is a design matrix with associated \eqn{\beta} fixed effect intercept coefficients,
#'  and Z is a design matrix with associated \eqn{\delta} random effects for the intercepts.
#'  For simplicity and easier interpretation, the unique item intercept values typically found in
#'  \eqn{X \beta}
#'  are extracted and reassigned within mirt's 'intercept' parameters (e.g., \code{'d'}).
#'  To observe how the design matrices are structured prior to reassignment and estimation pass
#'  the argument \code{return.design = TRUE}.
#'
#'  Polytomous IRT models follow a similar format except the item intercepts are automatically
#'  estimated internally, rendering the \code{items} argument in the fixed formula redundant and
#'  therefore must be omitted from the specification. If there are a mixture of dichotomous and
#'  polytomous items the intercepts for the dichotomous models are also estimated for consistency.
#'
#'  The decomposition of the \eqn{\theta} parameters is also possible to form
#'  latent regression and multilevel IRT models by using the \code{lr.fixed} and \code{lr.random}
#'  inputs. These effects decompose \eqn{\theta} such that
#'
#'  \deqn{\theta = V \Gamma + W \zeta + \epsilon}
#'
#'  where V and W are fixed and random effects design matrices for the associated coefficients.
#'
#'  To simulate expected a posteriori predictions for the random effect terms
#'  use the \code{\link{randef}} function.
#'
#' @return function returns an object of class \code{MixedClass}
#'   (\link{MixedClass-class}).
#'
#' @aliases mixedmirt
#' @param data a \code{matrix} or \code{data.frame} that consists of
#'   numerically ordered data, with missing data coded as \code{NA}
#' @param covdata a \code{data.frame} that consists of the \code{nrow(data)} by \code{K}
#'   'person level' fixed and random predictors. If missing data are present in this object
#'   then the observations from \code{covdata} and \code{data} will be removed row-wise
#'   via the \code{rowSums(is.na(covdata)) > 0}
#' @param model an object returned from, or a string to be passed to, \code{mirt.model()}
#'   to declare how the IRT model is to be estimated. See \code{\link{mirt.model}} and
#'   \code{\link{mirt}} for more detail
#' @param fixed a right sided R formula for specifying the fixed effect (aka 'explanatory')
#'   predictors from \code{covdata} and \code{itemdesign}. To estimate the intercepts for
#'   each item the keyword \code{items} is reserved and automatically added to the \code{itemdesign}
#'   input. If any polytomous items are being model the \code{items} are argument is not valid
#'   since all intercept parameters are freely estimated and identified with the parameterizations
#'   found in \code{\link{mirt}}, and the first column in the fixed design matrix
#'   (commonly the intercept or a reference group) is omitted
#' @param random a right sided formula or list of formulas containing crossed random effects
#'   of the form \code{v1 + ... v_n | G}, where \code{G} is the grouping variable and \code{v_n} are
#'   random numeric predictors within each group. If no intercept value is specified then by default
#'   the correlations between the \code{v}'s and \code{G} are estimated, but can be suppressed by
#'   including the \code{~ -1 + ...} or 0 constant. \code{G} may contain interaction terms, such as
#'   \code{group:items} to include cross or person-level interactions effects
#' @param itemtype same as itemtype in \code{\link{mirt}}, except when the \code{fixed}
#'   or \code{random} inputs are used does not support the following
#'   item types: \code{c('PC2PL', 'PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')}
#' @param itemdesign a \code{data.frame} object used to create a design matrix for the items, where
#'   each \code{nrow(itemdesign) == nitems} and the number of columns is equal to the number of
#'   fixed effect predictors (i.e., item intercepts). By default an \code{items} variable is
#'   reserved for modeling the item intercept parameters
#' @param lr.fixed an R formula (or list of formulas) to specify regression
#'   effects in the latent variables from the variables in \code{covdata}. This is used to construct models such as the so-called
#'   'latent regression model' to explain person-level ability/trait differences. If a named list
#'   of formulas is supplied (where the names correspond to the latent trait names in \code{model})
#'   then specific regression effects can be estimated for each factor. Supplying a single formula
#'   will estimate the regression parameters for all latent traits by default.
#' @param lr.random a list of random effect terms for modeling variability in the
#'   latent trait scores, where the syntax uses the same style as in the \code{random} argument.
#'   Useful for building so-called 'multilevel IRT' models which are non-Rasch (multilevel Rasch
#'   models do not technically require these because they can be built using the
#'   \code{fixed} and \code{random} inputs alone)
#' @param constrain a list indicating parameter equality constrains. See \code{\link{mirt}} for
#'   more detail
#' @param pars used for parameter starting values. See \code{\link{mirt}} for more detail
#' @param return.design logical; return the design matrices before they have (potentially)
#'   been reassigned?
#' @param SE logical; compute the standard errors by approximating the information matrix
#'   using the MHRM algorithm? Default is TRUE
#' @param internal_constraints logical; use the internally defined constraints for constraining
#'   effects across persons and items? Default is TRUE. Setting this to FALSE runs the risk of
#'   under-identification
#' @param technical the technical list passed to the MH-RM estimation engine, with the
#'   SEtol default increased to .0001. Additionally, the argument \code{RANDSTART} is available
#'   to indicate at which iteration (during the burn-in stage) the additional random effect
#'   variables should begin to be approximated (i.e.,
#'   elements in \code{lr.random} and \code{random}). The default for \code{RANDSTART} is to start
#'   at iteration 100, and when random effects are included the default number of burn-in iterations is increased
#'   from 150 to 200. See \code{\link{mirt}} for further details
#' @param ... additional arguments to be passed to the MH-RM estimation engine. See
#'   \code{\link{mirt}} for more details and examples
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{mirt}}, \code{\link{randef}}, \code{\link{fixef}}, \code{\link{boot.mirt}}
#' @export mixedmirt
#'
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Chalmers, R. P. (2015). Extended Mixed-Effects Item Response Models with the MH-RM Algorithm.
#' \emph{Journal of Educational Measurement, 52}, 200-222. \doi{10.1111/jedm.12072}
#'
#' @examples
#'
#' \dontrun{
#'
#' # make some data
#' set.seed(1234)
#' N <- 750
#' a <- matrix(rlnorm(10,.3,1),10,1)
#' d <- matrix(rnorm(10), 10)
#' Theta <- matrix(sort(rnorm(N)))
#' pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
#' pseudoIQ <- (pseudoIQ - mean(pseudoIQ))/10  #rescale variable for numerical stability
#' group <- factor(rep(c('G1','G2','G3'), each = N/3))
#' data <- simdata(a,d,N, itemtype = rep('2PL',10), Theta=Theta)
#' covdata <- data.frame(group, pseudoIQ)
#'
#' itemstats(data)
#'
#' # use parallel computing
#' if(interactive()) mirtCluster()
#'
#' # specify IRT model
#' model <- 'Theta = 1-10'
#'
#' # model with no person predictors
#' mod0 <- mirt(data, model, itemtype = 'Rasch')
#'
#' # group as a fixed effect predictor (aka, uniform dif)
#' mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items)
#' anova(mod0, mod1)
#' summary(mod1)
#' coef(mod1)
#'
#' # same model as above in lme4
#' wide <- data.frame(id=1:nrow(data),data,covdata)
#' long <- reshape2::melt(wide, id.vars = c('id', 'group', 'pseudoIQ'))
#' library(lme4)
#' lmod0 <- glmer(value ~ 0 + variable + (1|id), long, family = binomial)
#' lmod1 <- glmer(value ~ 0 + group + variable + (1|id), long, family = binomial)
#' anova(lmod0, lmod1)
#'
#' # model using 2PL items instead of Rasch
#' mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items, itemtype = '2PL')
#' anova(mod1, mod1b) #better with 2PL models using all criteria (as expected, given simdata pars)
#'
#' # continuous predictor with group
#' mod2 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items + pseudoIQ)
#' summary(mod2)
#' anova(mod1b, mod2)
#'
#' # view fixed design matrix with and without unique item level intercepts
#' withint <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, return.design = TRUE)
#' withoutint <- mixedmirt(data, covdata, model, fixed = ~ 0 + group, return.design = TRUE)
#'
#' # notice that in result above, the intercepts 'items1 to items 10' were reassigned to 'd'
#' head(withint$X)
#' tail(withint$X)
#' head(withoutint$X) # no intercepts design here to be reassigned into item intercepts
#' tail(withoutint$X)
#'
#' ###################################################
#' ### random effects
#' # make the number of groups much larger
#' covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
#'
#' # random groups
#' rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group)
#' summary(rmod1)
#' coef(rmod1)
#'
#' # random groups and random items
#' rmod2 <- mixedmirt(data, covdata, 1, random = list(~ 1|group, ~ 1|items))
#' summary(rmod2)
#' eff <- randef(rmod2) #estimate random effects
#'
#' # random slopes with fixed intercepts (suppressed correlation)
#' rmod3 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ -1 + pseudoIQ|group)
#' summary(rmod3)
#' eff <- randef(rmod3)
#' str(eff)
#'
#' ###################################################
#' ## LLTM, and 2PL version of LLTM
#' data(SAT12)
#' data <- key2binary(SAT12,
#'                    key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' model <- 'Theta = 1-32'
#'
#' # Suppose that the first 16 items were suspected to be easier than the last 16 items,
#' #   and we wish to test this item structure hypothesis (more intercept designs are possible
#' #   by including more columns).
#' itemdesign <- data.frame(itemorder = factor(c(rep('easier', 16), rep('harder', 16))))
#'
#' # notice that the 'fixed = ~ ... + items' argument is omitted
#' LLTM <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemdesign = itemdesign,
#'    SE = TRUE) # SE argument ensures that the information matrix is computed accurately
#' summary(LLTM)
#' coef(LLTM)
#' wald(LLTM)
#' L <- matrix(c(-1, 1, 0), 1)
#' wald(LLTM, L) #first half different from second
#'
#' # compare to items with estimated slopes (2PL)
#' twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL',
#'                    itemdesign = itemdesign)
#' # twoPL not mixing too well (AR should be between .2 and .5), decrease MHcand
#' twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL',
#'                   itemdesign = itemdesign, technical = list(MHcand = 0.8))
#' anova(twoPL, LLTM) #much better fit
#' summary(twoPL)
#' coef(twoPL)
#'
#' wald(twoPL)
#' L <- matrix(0, 1, 34)
#' L[1, 1] <- 1
#' L[1, 2] <- -1
#' wald(twoPL, L) # n.s., which is the correct conclusion. Rasch approach gave wrong inference
#'
#' ## LLTM with item error term
#' LLTMwithError <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, random = ~ 1|items,
#'     itemdesign = itemdesign)
#' summary(LLTMwithError)
#' # large item level variance after itemorder is regressed; not a great predictor of item difficulty
#' coef(LLTMwithError)
#'
#' ###################################################
#' ### Polytomous example
#'
#' # make an arbitrary group difference
#' covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
#'
#' # partial credit model
#' mod <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group)
#' coef(mod)
#'
#' # gpcm to estimate slopes
#' mod2 <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group,
#'                  itemtype = 'gpcm')
#' summary(mod2)
#' anova(mod, mod2)
#'
#' # graded model
#' mod3 <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group,
#'                  itemtype = 'graded')
#' coef(mod3)
#'
#'
#' ###################################################
#' # latent regression with Rasch and 2PL models
#'
#' set.seed(1)
#' n <- 300
#' a <- matrix(1, 10)
#' d <- matrix(rnorm(10))
#' Theta <- matrix(c(rnorm(n, 0), rnorm(n, 1), rnorm(n, 2)))
#' covdata <- data.frame(group=rep(c('g1','g2','g3'), each=n))
#' dat <- simdata(a, d, N=n*3, Theta=Theta, itemtype = '2PL')
#' itemstats(dat)
#'
#' # had we known the latent abilities, we could have computed the regression coefs
#' summary(lm(Theta ~ covdata$group))
#'
#' # but all we have is observed test data. Latent regression helps to recover these coefs
#' # Rasch model approach (and mirt equivalent)
#' rmod0 <- mirt(dat, 1, 'Rasch') # unconditional
#'
#' # these two models are equivalent
#' rmod1a <- mirt(dat, 1, 'Rasch', covdata = covdata, formula = ~ group)
#' rmod1b <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items + group)
#' anova(rmod0, rmod1b)
#' coef(rmod1a, simplify=TRUE)
#' summary(rmod1b)
#'
#' # 2PL, requires different input to allow Theta variance to remain fixed
#' mod0 <- mirt(dat, 1) # unconditional
#' mod1a <- mirt(dat, 1, covdata = covdata, formula = ~ group, itemtype = '2PL')
#' mod1b <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items, lr.fixed = ~group, itemtype = '2PL')
#' anova(mod0, mod1b)
#' coef(mod1a)$lr.betas
#' summary(mod1b)
#'
#' # specifying specific regression effects is accomplished by passing a list of formula
#' model <- 'F1 = 1-5
#'          F2 = 6-10'
#' covdata$contvar <- rnorm(nrow(covdata))
#' mod2 <- mirt(dat, model, itemtype = 'Rasch', covdata=covdata,
#'         formula = list(F1 = ~ group + contvar, F2 = ~ group))
#' coef(mod2)[11:12]
#' mod2b <- mixedmirt(dat, covdata, model, fixed = ~ 0 + items,
#'         lr.fixed = list(F1 = ~ group + contvar, F2 = ~ group))
#' summary(mod2b)
#'
#' ####################################################
#' ## Simulated Multilevel Rasch Model
#'
#' set.seed(1)
#' N <- 2000
#' a <- matrix(rep(1,10),10,1)
#' d <- matrix(rnorm(10))
#' cluster = 100
#' random_intercept = rnorm(cluster,0,1)
#' Theta = numeric()
#' for (i in 1:cluster)
#'     Theta <- c(Theta, rnorm(N/cluster,0,1) + random_intercept[i])
#'
#' group = factor(rep(paste0('G',1:cluster), each = N/cluster))
#' covdata <- data.frame(group)
#' dat <- simdata(a,d,N, itemtype = rep('2PL',10), Theta=matrix(Theta))
#' itemstats(dat)
#'
#' # null model
#' mod1 <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group)
#' summary(mod1)
#'
#' # include level 2 predictor for 'group' variance
#' covdata$group_pred <- rep(random_intercept, each = N/cluster)
#' mod2 <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items + group_pred, random = ~ 1|group)
#'
#' # including group means predicts nearly all variability in 'group'
#' summary(mod2)
#' anova(mod1, mod2)
#'
#' # can also be fit for Rasch/non-Rasch models with the lr.random input
#' mod1b <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items, lr.random = ~ 1|group)
#' summary(mod1b)
#'
#' mod2b <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items + group_pred, lr.random = ~ 1|group)
#' summary(mod2b)
#' anova(mod1b, mod2b)
#'
#' mod3 <- mixedmirt(dat, covdata, 1, fixed = ~ 0 + items, lr.random = ~ 1|group, itemtype = '2PL')
#' summary(mod3)
#' anova(mod1b, mod3)
#'
#' head(cbind(randef(mod3)$group, random_intercept))
#'
#' }
mixedmirt <- function(data, covdata = NULL, model = 1, fixed = ~ 1, random = NULL, itemtype = 'Rasch',
                      lr.fixed = ~ 1, lr.random = NULL, itemdesign = NULL, constrain = NULL,
                      pars = NULL, return.design = FALSE, SE = TRUE, internal_constraints = TRUE,
                      technical = list(SEtol = 1e-4), ...)
{
    Call <- match.call()
    dots <- list(...)
    opars <- pars
    if(!is.null(dots$method))
        stop('method cannot be changed in mixedmirt() function', call.=FALSE)
    svinput <- pars
    iconstrain <- constrain
    covdataold <- covdata
    completely_missing <- which(rowSums(is.na(data)) == ncol(data))
    if(length(completely_missing)){
        data <- data[-completely_missing, , drop=FALSE]
        if(!is.null(covdata))
            covdata <- covdata[-completely_missing, , drop=FALSE]
    }
    itemdesignold <- if(is.null(itemdesign)) data.frame() else itemdesign
    if(length(itemtype) == 1L) itemtype <- rep(itemtype, ncol(data))
    if(any(itemtype %in% c('spline', 'ideal'))){
        if(fixed != ~ -1){
            warning(paste0('Unsupported itemtype detected for modeling intercepts.\n',
                    'fixed = ~ -1 used by default to remove intercept'))
            fixed <- ~ -1
        }
        if(!is.null(random)){
            warning('random set to NULL due to unsupported itemtypes')
            random <- NULL
        }
    }
    if(fixed != ~ 1 || !is.null(random))
        if(any(itemtype %in% c('PC2PL', 'PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')))
            stop('itemtype contains unsupported classes of items when fixed or random are used', call.=FALSE)
    if(is(random, 'formula')) {
        random <- list(random)
    } else if(is.null(random)) random <- list()
    RETVALUES <- ifelse(is.character(pars), TRUE, FALSE)
    if(!is.list(random)) stop('Incorrect input for random argument', call.=FALSE)
    if(is.null(covdata)) covdata <- data.frame(UsElEsSvAR = factor(rep(1L, nrow(data))))
    if(is.null(itemdesign)){
        itemdesign <- data.frame(items = factor(1L:ncol(data)))
    } else {
        if(!is.null(itemdesign$items))
            stop('itemdesign internally reserves the predictor name \'items\'. Please Change ',
                 call.=FALSE)
        itemdesign$items <- factor(1L:ncol(data))
    }
    if(!is.data.frame(covdata) || ! is.data.frame(itemdesign))
        stop('Predictor variable inputs must be data.frame objects', call.=FALSE)
    if(nrow(covdata) != nrow(data))
        stop('number of rows in covdata do not match number of rows in data', call.=FALSE)
    dropcases <- which(rowSums(is.na(covdata)) != 0)
    if(length(dropcases) > 0L){
        if(is.null(technical$warn) || technical$warn)
            warning("Missing values in covdata not permitted. Removing all observations row-wise
                     when rowSums(is.na(covdata)) > 0")
        data <- data[-dropcases, ]
        covdata <- covdata[-dropcases, ,drop=FALSE]
    }
    pick <- sapply(covdata, is.numeric)
    if(any(pick)){
        if(any(covdata[,pick] > 10 | covdata[,pick] < -10))
            warning('Continuous variables in covdata should be rescaled to fall
                     between -10 and 10 for better numerical stability', call.=FALSE)
    }
    longdata <- reshape(data.frame(ID=1L:nrow(data), data, covdata), idvar='ID',
                        varying=list(1L:ncol(data) + 1L), direction='long')
    colnames(longdata) <- c('ID', colnames(covdata), 'items', 'response')
    for(i in 1L:ncol(itemdesign))
        longdata[, colnames(itemdesign)[i]] <- rep(itemdesign[ ,i], each=nrow(data))
    mf <- model.frame(fixed, longdata)
    mm <- model.matrix(fixed, mf)
    K <- sapply(as.data.frame(data), function(x) length(na.omit(unique(x))))
    if(any(K > 2)){
        if(any(colnames(mm) %in% paste0('items', 1:ncol(data))))
            stop('fixed formulas do no support the \'items\' internal variable for
                 polytomous items. Please remove', call.=FALSE)
        mm <- mm[ , -1L, drop = FALSE]
    }
    if(return.design) return(list(X=mm, Z=NaN))
    itemindex <- colnames(mm) %in% paste0('items', 1L:ncol(data))
    mmitems <- mm[, itemindex]
    mm <- mm[ ,!itemindex, drop = FALSE]
    if(length(random) > 0L){
        mr <- make.randomdesign(random=random, longdata=longdata, covnames=colnames(covdata),
                                itemdesign=itemdesign, N=nrow(covdata))
    } else mr <- list()
    mixed.design <- list(fixed=mm, random=mr)
    if(is(lr.random, 'formula')) lr.random <- list(lr.random)
    if(length(lr.random) > 0L){
        lr.random <- make.randomdesign(random=lr.random, longdata=covdata,
                                       covnames=colnames(covdata), itemdesign=NULL,
                                       N=nrow(covdata), LR=TRUE)
    } else lr.random <-
    if(is.null(constrain)) constrain <- list()
    if(is.list(lr.fixed) || (lr.fixed != ~ 1) || length(lr.random) > 0L){
        latent.regression <- list(df=covdata, formula=lr.fixed,
                                  EM=FALSE, lr.random=lr.random)
    } else latent.regression <- NULL
    sv <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), itemtype=itemtype,
                     D=1, mixed.design=mixed.design, method='MIXED', constrain=NULL, pars='values',
                     latent.regression=latent.regression, ...)
    mmnames <- colnames(mm)
    if(ncol(mm) > 0L){
        for(i in 1L:ncol(mm)){
            mmparnum <- sv$parnum[sv$name == mmnames[i]]
            constrain[[length(constrain) + 1L]] <- mmparnum
        }
    }
    if(ncol(mmitems) > 0L){
        itemindex <- colnames(data)[which(paste0('items', 1L:ncol(data)) %in% colnames(mmitems))]
        for(i in itemindex){
            tmp <- sv[i == sv$item, ]
            tmp$est[tmp$name %in% paste0('d', 1L:50L)] <- TRUE
            tmp$est[tmp$name == 'd'] <- TRUE
            sv[i == sv$item, ] <- tmp
        }
        attr(sv, 'values') <- pars
        pars <- sv
    } else {
        if(sum(sv$name == 'd') == ncol(data)){
            sv$value[sv$name == 'd'] <- 0
            sv$est[sv$name == 'd'] <- FALSE
        }
        pars <- sv
    }
    if(RETVALUES){
        attr(pars, 'values') <- NULL
        return(pars)
    }
    if(is.data.frame(svinput)) pars <- svinput
    if(!internal_constraints) constrain <- iconstrain
    attr(mixed.design, 'covdata') <- if(!is.null(covdataold)) covdataold else data.frame()
    attr(mixed.design, 'itemdesign') <- itemdesignold
    attr(mixed.design, 'formula') <- list(fixed=fixed, random=random, lr.fixed=lr.fixed,
                                          lr.random=lr.random)
    if(!is.null(opars)) pars <- opars
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), itemtype=itemtype,
                      mixed.design=mixed.design, method='MIXED', constrain=constrain, pars=pars,
                      SE=SE, latent.regression=latent.regression, technical=technical, ...)
    if(is(mod, 'MixedClass'))
        mod@Call <- Call
    mod@Data$completely_missing <- completely_missing
    return(mod)
}
