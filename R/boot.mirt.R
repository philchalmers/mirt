#' Calculate bootstrapped standard errors for estimated models
#'
#' Given an internal mirt object estimate the bootstrapped standard errors. It may
#' be beneficial to run the computations using multi-core architecture (e.g., the \code{parallel}
#' package).
#'
#' @aliases boot.mirt
#' @param x an estimated model object
#' @param R number of draws to use (passed to the \code{boot()} function)
#' @param ... additional arguments to be passed on to \code{boot(...)}
#' @keywords bootstrapped standard errors
#' @export boot.mirt
#' @seealso
#' \code{\link{PLCI.mirt}}
#' @examples
#'
#' \dontrun{
#'
#' #standard
#' mod <- mirt(Science, 1)
#' booted <- boot.mirt(mod)
#' plot(booted)
#' booted
#'
#' #run in parallel using snow back-end using all available cores
#' mod <- mirt(Science, 1)
#' booted <- boot.mirt(mod, parallel = 'snow', ncpus = parallel::detectCores())
#' booted
#'
#'
#' }
boot.mirt <- function(x, R = 100, ...){
    boot.draws <- function(orgdat, ind, npars, constrain, parprior, model, itemtype, group,
                           class, LR, obj, ...) {
        ngroup <- length(unique(group))
        dat <- orgdat[ind, ]
        rownames(dat) <- NULL
        g <- group[ind]
        if(length(unique(g)) != ngroup) return(rep(NA, npars))
        if(class == 'MixedClass'){
            fm <- obj@formulas
            itemdesign <- if(length(obj@itemdesign)) obj@itemdesign else NULL
            covdata <- if(length(obj@covdata)) obj@covdata[ind, , drop=FALSE] else NULL
            mod <- try(mixedmirt(data=dat, model=model, covdata=covdata, itemtype=itemtype,
                                 itemdesign=itemdesign, fixed=fm$fixed, random=fm$random,
                                 lr.fixed=fm$lr.fixed, lr.random=fm$lr.random,
                                 constrain=constrain, parprior=parprior,
                                 verbose=FALSE, technical=list(parallel=FALSE),
                                 SE=FALSE, ...))
        } else if(class == 'DiscreteClass'){
            mod <- try(mdirt(data=dat, model=model, itemtype=itemtype, group=g,
                             constrain=constrain, parprior=parprior,
                             verbose=FALSE, technical=list(parallel=FALSE),
                             SE=FALSE, ...))
        } else {
            if(class == 'MultipleGroupClass'){
                mod <- try(multipleGroup(data=dat, model=model, itemtype=itemtype, group=g,
                                     constrain=constrain, parprior=parprior,
                                     calcNull=FALSE, verbose=FALSE, technical=list(parallel=FALSE),
                                     method=x@method, draws=1, SE=FALSE, ...))
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
                            technical=list(parallel=FALSE), formula=formula,
                            covdata=df, method=x@method, draws=1, SE=FALSE, ...))
            }
        }
        if(is(mod, 'try-error')) return(rep(NA, npars))
        structure <- mod2values(mod)
        longpars <- structure$value[structure$est]
        if(length(longpars) != npars) return(rep(NA, npars)) #in case intercepts dropped
        return(longpars)
    }

    return.boot <- TRUE
    dat <- x@Data$data
    method <- x@method
    itemtype <- x@itemtype
    class <- class(x)
    group <- if(class == 'MultipleGroupClass') x@Data$group else NULL
    model <- x@model[[1L]]
    parprior <- x@parprior
    constrain <- x@constrain
    LR <- x@lrPars
    if(length(parprior) == 0L) parprior <- NULL
    if(length(constrain) == 0L) constrain <- NULL
    prodlist <- x@prodlist
    ret <- x
    if(!require(boot)) require('boot')
    structure <- mod2values(x)
    longpars <- structure$value
    npars <- length(longpars)
    boots <- boot::boot(dat, boot.draws, R=R, npars=npars, constrain=constrain, class=class,
                  parprior=parprior, model=model, itemtype=itemtype, group=group, LR=LR,
                  obj=x, ...)
    if(class == 'ExploratoryClass')
        message('Note: bootstrapped standard errors for slope parameters for exploratory
                       models are not meaningful.')
    return(boots)
}
