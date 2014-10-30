#' Calculate bootstrapped standard errors for estimated models
#'
#' Given an internal mirt object estimate the bootstrapped standard errors. It may
#' be beneficial to run the computations using multi-core architecture (e.g., the \code{parallel}
#' package).
#'
#' @aliases boot.mirt
#' @param x an estimated object from \code{mirt}, \code{bfactor}, \code{multipleGroup},
#'   or \code{mdirt}
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
                           discrete, LR, ...) {
        ngroup <- length(unique(group))
        dat <- orgdat[ind, ]
        g <- group[ind]
        if(length(unique(g)) != ngroup) return(rep(NA, npars))
        if(discrete){
            mod <- try(mdirt(data=dat, model=model, itemtype=itemtype, group=g,
                                     constrain=constrain, parprior=parprior,
                                     verbose=FALSE, technical=list(parallel=FALSE),
                                     ...))
        } else {
            if(!is.null(group)){
                mod <- try(multipleGroup(data=dat, model=model, itemtype=itemtype, group=g,
                                     constrain=constrain, parprior=parprior,
                                     calcNull=FALSE, verbose=FALSE, technical=list(parallel=FALSE),
                                     ...))
            } else {
                mod <- try(mirt(data=dat, model=model, itemtype=itemtype, constrain=constrain,
                            parprior=parprior, calcNull=FALSE, verbose=FALSE,
                            technical=list(parallel=FALSE), formula=LR$formula,
                            covdata=LR$X, ...))
            }
        }
        if(is(mod, 'try-error')) return(rep(NA, npars))
        structure <- mod2values(mod)
        longpars <- structure$value
        if(length(longpars) != npars) return(rep(NA, npars)) #in case intercepts dropped
        if(length(LR)){
            betas <- mod@pars[[length(mod@pars)]]@betas
            longpars <- c(longpars, as.numeric(betas))
        }
        return(longpars)
    }

    if(is(x, 'MixedClass'))
        stop('Bootstapped standard errors not supported for MixedClass objects')
    return.boot <- TRUE
    dat <- x@Data$data
    method <- x@method
    itemtype <- x@itemtype
    MG <- is(x, 'MultipleGroupClass')
    explor <- is(x, 'ExploratoryClass')
    discrete <- is(x, 'DiscreteClass')
    group <- if(MG) x@Data$group else NULL
    model <- x@model[[1L]]
    parprior <- x@parprior
    constrain <- x@constrain
    LR <- if(is(x, 'ConfirmatoryClass')) x@l.regress else list()
    if(length(parprior) == 0L) parprior <- NULL
    if(length(constrain) == 0L) constrain <- NULL
    prodlist <- x@prodlist
    ret <- x
    if(!require(boot)) require('boot')
    structure <- mod2values(x)
    longpars <- structure$value
    npars <- length(longpars)
    boots <- boot::boot(dat, boot.draws, R=R, npars=npars, constrain=constrain, discrete=discrete,
                  parprior=parprior, model=model, itemtype=itemtype, group=group, LR=LR, ...)
    if(explor) message('Note: bootstrapped standard errors for slope parameters for exploratory
                       models are not meaningful.')
    return(boots)
}





