#' Calculate bootstrapped standard errors for estimated models
#'
#' Given an internal mirt object estimate the bootstrapped standard errors. If possible, it will 
#' be benifitial to run the computations using multicore architecher (e.g., using the \code{parallel} 
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
#' @examples 
#' 
#' \dontrun{
#' mod <- mirt(Science, 1)
#' booted <- boot.mirt(mod, R = 100)
#' booted
#' modwithboot <- boot.mirt(mod, R = 100, return.boot = FALSE)
#' coef(modwithboot)
#' 
#' } 
boot.mirt <- function(x, R = 1000, return.boot = TRUE, ...){    
    boot.draws <- function(orgdat, ind, npars, constrain, parprior, model, itemtype, group) { 
        ngroup <- length(unique(group))             
        while(TRUE){
            dat <- orgdat[ind, ]
            g <- group[ind]
            ind <- sample(1:nrow(orgdat), nrow(orgdat), TRUE)            
            if(length(unique(g)) != ngroup) next                
            if(!is.null(group)){
                mod <- try(multipleGroup(data=dat, model=model, itemtype=itemtype, group=g, 
                                     constrain=constrain, parprior=parprior, method='EM', 
                                     calcNull=FALSE, verbose = FALSE))
            } else {
                mod <- try(mirt(data=dat, model=model, itemtype=itemtype, constrain=constrain, 
                            parprior=parprior, calcNull=FALSE))
            }
            if(is(mod, 'try-error')) next
            if(MG){
                longpars <- c()
                tmp <- coef(mod)
                for(g in 1:length(tmp))
                    longpars <- c(longpars, do.call(c, tmp[[g]]))                
            } else longpars <- do.call(c, coef(mod))            
            if(length(longpars) != npars) next
            break
        }
        return(longpars)
    }    
    loadSE <- function(pars, SEs, nfact, MG, explor){        
        ind1 <- 1                
        if(MG){            
            for(g in 1:length(pars)){
                for(i in 1:length(pars[[g]]@pars)){
                    ind2 <- ind1 + length(pars[[g]]@pars[[i]]@par) - 1
                    pars[[g]]@pars[[i]]@SEpar <- SEs[ind1:ind2]
                    ind1 <- ind2 + 1
                }
            }            
        } else {
            for(i in 1:length(x@pars)){
                ind2 <- ind1 + length(pars[[i]]@par) - 1
                pars[[i]]@SEpar <- SEs[ind1:ind2]
                if(explor) pars[[i]]@SEpar[1:nfact] <- NA 
                ind1 <- ind2 + 1
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
    model <- x@model[[1]]
    parprior <- x@parprior    
    constrain <- x@constrain    
    if(length(parprior) == 0) parprior <- NULL
    if(length(constrain) == 0) constrain <- NULL
    prodlist <- x@prodlist    
    ret <- x
    if(!require(boot)) require(boot)        
    if(MG){
        longpars <- c()
        tmp <- coef(x)
        for(g in 1:length(tmp))
            longpars <- c(longpars, do.call(c, tmp[[g]]))                
    } else longpars <- do.call(c, coef(x))
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
    
    
    
    
    
