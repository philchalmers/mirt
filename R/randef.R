#' Compute random effects
#'
#' Stochastically compute random effects for \code{MixedClass} objects with Metropolis-Hastings
#' samplers and averaging over the draws. Returns a list of the estimated effects.
#'
#' @aliases randef
#' @param x an estimated model object from the \code{\link{mixedmirt}} function
#' @param ndraws total number of draws to perform. Default is 1000
#' @param thin amount of thinning to apply. Default is to use every 10th draw
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords random effects
#' @export randef
#' @examples
#' \dontrun{
#' effects <- randef(mod1, ndraws = 2000, thin = 20)
#'
#' }
randef <- function(x, ndraws = 1000, thin = 10){
    if(!is(x, 'MixedClass'))
        stop('Only applicable to MixedClass objects')
    div <- ndraws / thin
    if(!closeEnough(floor(ndraws/thin) == (ndraws/thin), -1e4, 1e4))
        stop('ndraws and thin are not the correct dimensions')
    random <- x@random
    if(length(random) > 0L){
        Random <- vector('list', length(random))
        for(i in 1L:length(Random))
            Random[[i]] <- matrix(0, nrow(x@random[[i]]@drawvals), ncol(x@random[[i]]@drawvals))
    }
    J <- ncol(x@data)
    N <- nrow(x@fulldata)
    Theta <- tmpTheta <- matrix(0, N, x@nfact)
    if(length(random) > 0L){
        OffTerm <- OffTerm(random, J=J, N=N)
    } else OffTerm <- matrix(0, 1, ncol(x@data))
    gstructgrouppars <- ExtractGroupPars(x@pars[[J+1L]])
    CUSTOM.IND <- x@CUSTOM.IND
    for(i in 1L:20L){
        tmpTheta <- draw.thetas(theta0=tmpTheta, pars=x@pars, fulldata=x@fulldata,
                                itemloc=x@itemloc, cand.t.var=x@cand.t.var,
                                prior.t.var=gstructgrouppars$gcov, OffTerm=OffTerm,
                                prior.mu=gstructgrouppars$gmeans, prodlist=list(),
                                CUSTOM.IND=CUSTOM.IND)
        if(length(random) > 0L){
            for(j in 1L:length(random))
                random[[j]]@drawvals <- DrawValues(random[[j]], Theta=tmpTheta, itemloc=x@itemloc,
                                                   pars=x@pars, fulldata=x@fulldata,
                                                   offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND)
            OffTerm <- OffTerm(random, J=J, N=N)
        }
    }
    for(i in 1L:ndraws){
        tmpTheta <- draw.thetas(theta0=tmpTheta, pars=x@pars, fulldata=x@fulldata,
                                itemloc=x@itemloc, cand.t.var=x@cand.t.var,
                                prior.t.var=gstructgrouppars$gcov, OffTerm=OffTerm,
                                prior.mu=gstructgrouppars$gmeans, prodlist=list(),
                                CUSTOM.IND=CUSTOM.IND)
        if(i %% thin == 0) Theta <- Theta + tmpTheta
        if(length(random) > 0L){
            for(j in 1L:length(random)){
                random[[j]]@drawvals <- DrawValues(random[[j]], Theta=tmpTheta, itemloc=x@itemloc,
                                                   pars=x@pars, fulldata=x@fulldata,
                                                   offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND)
                if(i %% thin == 0) Random[[j]] <- Random[[j]] + random[[j]]@drawvals
            }
            OffTerm <- OffTerm(random, J=J, N=N)
        }
    }
    Theta <- Theta / (ndraws/thin)
    attr(Theta, 'Proportion Accepted') <- attr(Theta, 'log.lik') <- NULL
    colnames(Theta) <- x@factorNames
    ret <- list(Theta)
    retnames <- 'Theta'
    if(length(random) > 0L){
        for(j in 1L:length(random)){
            Random[[j]] <- Random[[j]] / (ndraws/thin)
            attr(Random[[j]], 'Proportion Accepted') <- NULL
            colnames(Random[[j]]) <- colnames(x@random[[j]]@gdesign)
            ret[[length(ret) + 1L]] <- Random[[j]]
            retnames <- c(retnames, colnames(x@random[[j]]@gdesign)[1L])
        }
    }
    names(ret) <- retnames
    ret
}
