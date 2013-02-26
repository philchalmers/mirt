#' Calculate standard errors for estimated model
#'
#' Given an internal mirt object estimate the standard errors, and, if possible, 
#' the parameter information matrix.
#' 
#' @aliases calcSE
#' @param x an estimated object
#' @param SE.type a character indicating which class of standard errors to compute. Can be 
#' \code{'boot'} for boostrapped standard errors (does not compute the information matrix), 
#' \code{'BL'} for Bock and Leiberman computed standard errors, or \code{'MC'} for Monte Carlo
#' compution of the information matrix to obtain stanard errors. 
#' If \code{SE = TRUE} was used at estimation runtime (indicating that 
#' MH-RM standard errors were computed) then these will be overwritten in the returned object. 
#' @param R number of draws to use (also passed to the \code{boot()} function)
#' @param return.boot logical; return the estimated object from the \code{boot} package?
#' @param ... additional arguments to be passed
#' @keywords standard errors
#' @export calcSE
#' @examples 
#' 
#' \dontrun{
#' mod <- mirt(Science, 1)
#' modwithSE <- calcSE(mod)
#' coef(modwithSE)
#' } 
calcSE <- function(x, SE.type = 'boot', R = 1000, return.boot = FALSE, ...){
    BL.LL <- function(pars, constrain){}
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
    if(method == 'MHRM' || method == 'MIXED' && SE.type == 'MHRM')
        stop('MHRM standard errors already calculated during estimation.')
    ret <- x
    if(SE.type == 'boot'){
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
        if(return.boot) return(boots)                   
        ret@information <- matrix(0)
        SEs <- apply(boots$t,2, sd)
        SEs[SEs == 0] <- NA
        retpars <- loadSE(pars=pars, SEs=SEs, nfact=x@nfact, MG=MG, explor=explor)
        if(MG) ret@cmods <- retpars else ret@pars <- retpars    
    }
    if(SE.type == 'MHRM'){
        stop('MHRM standard errors must be calulated during the estimation runtime.')
    }
    if(SE.type == 'MC'){    
        
        
    }
    return(ret)
}
    
    
    
    
    
