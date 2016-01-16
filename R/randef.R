#' Compute posterior estimates of random effect
#'
#' Stochastically compute random effects for \code{MixedClass} objects with Metropolis-Hastings
#' samplers and averaging over the draws. Returns a list of the estimated effects.
#'
#' @aliases randef
#' @param x an estimated model object from the \code{\link{mixedmirt}} function
#' @param ndraws total number of draws to perform. Default is 1000
#' @param thin amount of thinning to apply. Default is to use every 10th draw
#' @param return.draws logical; return a list containing the thinned draws of the posterior?
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords random effects
#' @export randef
#' @examples
#' \dontrun{
#' #make an arbitrary groups
#' covdat <- data.frame(group = rep(paste0('group', 1:49), each=nrow(Science)/49))
#'
#' #partial credit model
#' mod <- mixedmirt(Science, covdat, model=1, random = ~ 1|group)
#' summary(mod)
#'
#' effects <- randef(mod, ndraws = 2000, thin = 20)
#' head(effects$Theta)
#' head(effects$group)
#'
#' }
randef <- function(x, ndraws = 1000, thin = 10, return.draws=FALSE){
    if(missing(x)) missingMsg('x')
    if(!is(x, 'MixedClass'))
        stop('Only applicable to MixedClass objects', call.=FALSE)
    if(!closeEnough(floor(ndraws/thin) == (ndraws/thin), -1e4, 1e4))
        stop('ndraws and thin are not the correct dimensions', call.=FALSE)
    random <- x@ParObjects$random
    if(length(random) > 0L){
        Random <- vector('list', length(random))
        for(i in 1L:length(Random))
            Random[[i]] <- matrix(0, nrow(x@ParObjects$random[[i]]@drawvals), ncol(x@ParObjects$random[[i]]@drawvals))
    }
    J <- ncol(x@Data$data)
    N <- nrow(x@Data$fulldata[[1L]])
    Theta <- tmpTheta <- matrix(0, N, x@Model$nfact)
    if(length(random) > 0L){
        OffTerm <- OffTerm(random, J=J, N=N)
    } else OffTerm <- matrix(0, 1, ncol(x@Data$data))
    gstructgrouppars <- ExtractGroupPars(x@ParObjects$pars[[J+1L]])
    CUSTOM.IND <- x@Internals$CUSTOM.IND
    if(length(x@Model$lrPars))
        gstructgrouppars$gmeans <- fixef(x)
    prodlist <- attr(x@ParObjects$pars, 'prodlist')
    for(i in 1L:20L){
        tmpTheta <- draw.thetas(theta0=tmpTheta, pars=x@ParObjects$pars, fulldata=x@Data$fulldata[[1L]],
                                itemloc=x@Model$itemloc, cand.t.var=x@OptimInfo$cand.t.var,
                                prior.t.var=gstructgrouppars$gcov, OffTerm=OffTerm,
                                prior.mu=gstructgrouppars$gmeans, prodlist=prodlist,
                                CUSTOM.IND=CUSTOM.IND)
        if(length(random) > 0L){
            for(j in 1L:length(random))
                random[[j]]@drawvals <- DrawValues(random[[j]], Theta=tmpTheta, itemloc=x@Model$itemloc,
                                                   pars=x@ParObjects$pars, fulldata=x@Data$fulldata[[1L]],
                                                   offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND)
            OffTerm <- OffTerm(random, J=J, N=N)
        }
    }
    if(return.draws){
        DRAWS <- vector('list', 1 + length(random))
        retnames <- "Theta"
        if(length(random))
            retnames <- c('Theta', sapply(random, function(x) colnames(x@gdesign)[1L]))
        names(DRAWS) <- retnames
    }
    for(i in 1L:ndraws){
        tmpTheta <- draw.thetas(theta0=tmpTheta, pars=x@ParObjects$pars, fulldata=x@Data$fulldata[[1L]],
                                itemloc=x@Model$itemloc, cand.t.var=x@OptimInfo$cand.t.var,
                                prior.t.var=gstructgrouppars$gcov, OffTerm=OffTerm,
                                prior.mu=gstructgrouppars$gmeans, prodlist=prodlist,
                                CUSTOM.IND=CUSTOM.IND)
        if(i %% thin == 0){
            Theta <- Theta + tmpTheta
            if(return.draws) DRAWS[[1L]][[length(DRAWS[[1L]]) + 1L]] <-
                tmpTheta[,1L:ncol(Theta), drop=FALSE]
        }
        if(length(random) > 0L){
            for(j in 1L:length(random)){
                random[[j]]@drawvals <- DrawValues(random[[j]], Theta=tmpTheta, itemloc=x@Model$itemloc,
                                                   pars=x@ParObjects$pars, fulldata=x@Data$fulldata[[1L]],
                                                   offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND)
                if(i %% thin == 0){
                    Random[[j]] <- Random[[j]] + random[[j]]@drawvals
                    if(return.draws) DRAWS[[j+1L]][[length(DRAWS[[j+1L]]) + 1L]] <-
                        random[[j]]@drawvals[,1L:ncol(random[[j]]@drawvals), drop=FALSE]
                }
            }
            OffTerm <- OffTerm(random, J=J, N=N)
        }
    }
    if(return.draws) return(DRAWS)
    Theta <- Theta / (ndraws/thin)
    attr(Theta, 'Proportion Accepted') <- attr(Theta, 'log.lik') <- NULL
    nfact <- extract.mirt(x, 'nfact')
    colnames(Theta) <- x@Model$factorNames[1L:nfact]
    ret <- list(Theta)
    retnames <- 'Theta'
    if(length(random) > 0L){
        for(j in 1L:length(random)){
            Random[[j]] <- Random[[j]] / (ndraws/thin)
            attr(Random[[j]], 'Proportion Accepted') <- NULL
            colnames(Random[[j]]) <- colnames(x@ParObjects$random[[j]]@gdesign)
            ret[[length(ret) + 1L]] <- Random[[j]]
            retnames <- c(retnames, colnames(x@ParObjects$random[[j]]@gdesign)[1L])
        }
    }
    ret <- lapply(ret, function(x){attr(x, 'log.lik_full') <- NULL; x} )
    names(ret) <- retnames
    ret
}
