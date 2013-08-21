#' Calculate bootstrapped standard errors for estimated models
#'
#' Given an internal mirt object estimate the bootstrapped standard errors. It may
#' be beneficial to run the computations using multicore architecture (e.g., the \code{parallel}
#' package).
#'
#' @aliases boot.mirt
#' @param x an estimated object from \code{mirt}, \code{bfactor}, or \code{multipleGroup}
#' @param R number of draws to use (passed to the \code{boot()} function)
#' @param return.boot logical; return the estimated object from the \code{boot} package? If \code{FALSE}
#' the estimated model is returned with the bootstrapped standard errors
#' @param ... additional arguments to be passed on to \code{boot(...)}
#' @keywords bootstrapped standard errors
#' @export boot.mirt
#' @seealso
#' \code{\link{PLCI.mirt}}
#' @examples
#'
#' \dontrun{
#' mod <- mirt(Science, 1)
#' booted <- boot.mirt(mod)
#' booted
#' 
#' #run in parallel using snow backend
#' modwithboot <- boot.mirt(mod, return.boot = FALSE, parallel = 'snow', ncpus = 4L)
#' coef(modwithboot)
#'
#' }
boot.mirt <- function(x, R = 100, return.boot = TRUE, ...){
    boot.draws <- function(orgdat, ind, npars, constrain, parprior, model, itemtype, group) {
        ngroup <- length(unique(group))
        dat <- orgdat[ind, ]
        g <- group[ind]
        if(length(unique(g)) != ngroup) return(rep(NA, npars))
        if(!is.null(group)){
            mod <- try(multipleGroup(data=dat, model=model, itemtype=itemtype, group=g,
                                 constrain=constrain, parprior=parprior, method='EM',
                                 calcNull=FALSE, verbose = FALSE))
        } else {
            mod <- try(mirt(data=dat, model=model, itemtype=itemtype, constrain=constrain,
                        parprior=parprior, calcNull=FALSE, verbose=FALSE))
        }
        if(is(mod, 'try-error')) return(rep(NA, npars))
        structure <- mod2values(mod)
        longpars <- structure$value
        if(length(longpars) != npars) return(rep(NA, npars)) #in case intercepts dropped
        return(longpars)
    }
    loadSE <- function(pars, SEs, nfact, MG, explor){
        ind1 <- 1L
        if(MG){
            for(g in 1L:length(pars)){
                for(i in 1L:length(pars[[g]]@pars)){
                    ind2 <- ind1 + length(pars[[g]]@pars[[i]]@par) - 1L
                    pars[[g]]@pars[[i]]@SEpar <- SEs[ind1:ind2]
                    ind1 <- ind2 + 1L
                }
            }
        } else {
            for(i in 1L:length(x@pars)){
                ind2 <- ind1 + length(pars[[i]]@par) - 1L
                pars[[i]]@SEpar <- SEs[ind1:ind2]
                if(explor) pars[[i]]@SEpar[1L:nfact] <- NA
                ind1 <- ind2 + 1L
            }
        }
        return(pars)
    }
    if(is(x, 'MixedClass'))
        stop('Bootstapped standard errors not supported for MixedClass objects')
    dat <- x@data
    method <- x@method
    itemtype <- x@itemtype
    MG <- is(x, 'MultipleGroupClass')
    explor <- is(x, 'ExploratoryClass')
    pars <- if(MG) x@cmods else x@pars
    group <- if(MG) x@group else NULL
    model <- x@model[[1L]]
    parprior <- x@parprior
    constrain <- x@constrain
    if(length(parprior) == 0L) parprior <- NULL
    if(length(constrain) == 0L) constrain <- NULL
    prodlist <- x@prodlist
    ret <- x
    if(!require(boot)) require(boot)
    structure <- mod2values(x)
    longpars <- structure$value
    npars <- length(longpars)
    boots <- boot(dat, boot.draws, R=R, npars=npars, constrain=constrain,
                  parprior=parprior, model=model, itemtype=itemtype, group=group, ...)
    if(return.boot){
        if(explor) message('Note: bootstrapped standard errors for slope parameters for exploratory
                           models are not meaningful.')
        return(boots)
    }
    ret@information <- matrix(0)
    SEs <- apply(boots$t,2, sd)
    SEs[SEs == 0] <- NA
    retpars <- loadSE(pars=pars, SEs=SEs, nfact=x@nfact, MG=MG, explor=explor)
    if(MG) ret@cmods <- retpars else ret@pars <- retpars
    return(ret)
}





