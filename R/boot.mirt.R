#' Calculate bootstrapped standard errors for estimated models
#'
#' Given an internal mirt object estimate the bootstrapped standard errors. It may
#' be beneficial to run the computations using multi-core architecture (e.g., the \code{parallel}
#' package). Parameters are organized from the freely estimated values in \code{mod2values(x)}
#' (equality constraints will also be returned in the bootstrapped estimates).
#'
#' @aliases boot.mirt
#' @param x an estimated model object
#' @param R number of draws to use (passed to the \code{boot()} function)
#' @param technical technical arguments passed to estimation engine. See \code{\link{mirt}}
#'   for details
#' @param boot.fun a user-defined function used to extract the information from the bootstrap
#'   fitted models. Must be of the form \code{boot.fun(x)}, where \code{x} is the
#'   bootstrap fitted model under investigation, and the return must be a numeric vector. If
#'   omitted a default function will be defined internally that returns the estimated
#'   parameters from the \code{mod} object, resulting in bootstrapped parameter estimate
#'   results
#' @param ... additional arguments to be passed on to \code{boot(...)} and mirt's
#'   estimation engine
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords bootstrapped standard errors
#' @export boot.mirt
#' @examples
#'
#' \donttest{
#'
#' # standard
#' mod <- mirt(Science, 1)
#' booted <- boot.mirt(mod, R=20)
#' plot(booted)
#' booted
#'
#' #run in parallel using snow back-end using all available cores
#' mod <- mirt(Science, 1)
#' booted <- boot.mirt(mod, parallel = 'snow', ncpus = parallel::detectCores())
#' booted
#'
#' ####
#' # bootstrapped CIs for standardized factor loadings
#' boot.fun <- function(mod){
#'   so <- summary(mod, verbose=FALSE)
#'   as.vector(so$rotF)
#' }
#'
#' # test to see if it works before running
#' boot.fun(mod)
#'
#' # run
#' booted.loads <- boot.mirt(mod, boot.fun=boot.fun)
#' booted.loads
#'
#' }
boot.mirt <- function(x, R = 100, boot.fun = NULL, technical = NULL, ...){
    boot.draws <- function(orgdat, ind, boot.fun, invariance,
                           npars, constrain, parprior, model, itemtype, group,
                           class, LR, obj, DTF = NULL, technical, ...) {
        ngroup <- length(unique(group))
        dat <- orgdat[ind, ]
        rownames(dat) <- NULL
        g <- group[ind]
        if(length(unique(g)) != ngroup) return(rep(NA, npars))
        mod <- NULL
        if(class == 'MixedClass'){
            fm <- obj@Model$formulas
            itemdesign <- if(length(obj@Data$itemdesign)) obj@Data$itemdesign else NULL
            covdata <- if(length(obj@Data$covdata)) obj@Data$covdata[ind, , drop=FALSE] else NULL
            mod <- try(mixedmirt(data=dat, model=model, covdata=covdata, itemtype=itemtype,
                                 itemdesign=itemdesign, fixed=fm$fixed, random=fm$random,
                                 lr.fixed=fm$lr.fixed, lr.random=fm$lr.random,
                                 constrain=constrain, parprior=parprior, draws=1,
                                 verbose=FALSE, technical=technical,
                                 SE=FALSE, ...))
        } else if(class == 'DiscreteClass'){
            mod <- try(mdirt(data=dat, model=model, itemtype=itemtype, group=g,
                             constrain=constrain, parprior=parprior,
                             verbose=FALSE, technical=technical,
                             SE=FALSE, ...))
        } else {
            if(class == 'MultipleGroupClass'){
                mod <- try(multipleGroup(data=dat, model=model, itemtype=itemtype, group=g,
                                     constrain=constrain, parprior=parprior, invariance=invariance,
                                     calcNull=FALSE, verbose=FALSE, technical=technical,
                                     method=x@Options$method, draws=1, SE=FALSE, ...))
            } else {
                if(.hasSlot(LR, 'beta')){
                    formula <- LR@formula
                    if(length(formula) == 1L) formula <- formula[[1]]
                    df <- LR@df[ind, ]
                } else {
                    formula = ~ 1
                    df <- NULL
                }
                mod <- try(mirt(data=dat, model=model, itemtype=itemtype, constrain=constrain,
                            parprior=parprior, calcNull=FALSE, verbose=FALSE,
                            technical=technical, formula=formula,
                            covdata=df, method=x@Options$method, draws=1, SE=FALSE, ...))
            }
        }
        if(!is.null(DTF)){
            if(is(mod, 'try-error')) return(rep(NA, 4L))
            return(calc_DTFs(mod=mod, Theta = DTF$Theta, max_score = DTF$max_score, plot='none',
                             type=DTF$type))
        }
        if(is(mod, 'try-error') || !extract.mirt(mod, 'converged'))
            return(rep(NA, npars))
        ret <- boot.fun(mod)
        ret
    }

    if(missing(x)) missingMsg('x')
    if(x@Options$exploratory)
        warning('Note: bootstrapped standard errors for slope parameters in exploratory
                       models are not meaningful.', call.=FALSE)
    dat <- x@Data$data
    itemtype <- x@Model$itemtype
    class <- class(x)
    group <- if(class == 'MultipleGroupClass') x@Data$group else NULL
    model <- x@Model$model
    parprior <- x@Model$parprior
    constrain <- x@Model$constrain
    customGroup <- extract.mirt(x, 'customGroup')
    customItems <- extract.mirt(x, 'customItems')
    LR <- x@Model$lrPars
    if(length(parprior) == 0L) parprior <- NULL
    if(length(constrain) == 0L) constrain <- NULL
    structure <- mod2values(x)
    longpars <- structure$value
    npars <- sum(structure$est)
    invariance <- extract.mirt(x, 'invariance')
    if(is.null(boot.fun)){
        boot.fun <- function(x){
            structure <- mod2values(x)
            longpars <- structure$value[structure$est]
            if(length(longpars) != npars) return(rep(NA, npars)) #in case intercepts dropped
            longpars
        }
    } else npars <- length(boot.fun(x))
    if(is.null(technical)) technical <- list(parallel=FALSE)
    else technical$parallel <- FALSE
    if(requireNamespace("boot", quietly = TRUE)){
      boots <- boot::boot(dat, boot.draws, R=R, npars=npars,
                          constrain=constrain, class=class, invariance=invariance,
                          customGroup=customGroup, customItems=customItems,
                          parprior=parprior, model=model, itemtype=itemtype, group=group, LR=LR,
                          obj=x, technical=technical, boot.fun=boot.fun, ...)
      boots$call <- match.call()
    }
    if(!is.null(DTF)) return(boots)
    names(boots$t0) <- paste(paste(structure$item[structure$est],
                             structure$name[structure$est], sep='.'),
                             structure$parnum[structure$est], sep='_')
    return(boots)
}
