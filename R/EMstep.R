# Estep for mirt
Estep.mirt <- function(pars, tabdata, Theta, prior, itemloc, debug) 
{   
    if(debug == 'Estep') browser()
    nfact <- ncol(Theta)
    nquad <- nrow(Theta)	
    r <- tabdata[ ,ncol(tabdata)]
    X <- tabdata[ ,1:(ncol(tabdata) - 1)]	
    itemtrace <- matrix(0, ncol=ncol(X), nrow=nrow(Theta))	
    for (i in 1:length(pars))
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)
    retlist <- .Call("Estep", itemtrace, prior, X, nfact, r)    
    return(retlist)
} 

# Estep for bfactor
Estep.bfactor <- function(pars, tabdata, Theta, prior, specific, sitems, itemloc, debug) 
{	    
    if(debug == 'Estep') browser()
    nfact <- pars[[1]]@nfact
    J <- length(pars)
    nquad <- nrow(Theta)		
    r <- tabdata[ ,ncol(tabdata)]
    X <- tabdata[ ,1:(ncol(tabdata) - 1)]	
    itemtrace <- matrix(0, ncol=ncol(X), nrow=nrow(Theta))	
    for (i in 1:J)
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)			
    retlist <- .Call("Estepbfactor", itemtrace, prior, X, r, sitems)	
    r1 <- matrix(0, nrow(Theta), ncol(X))	
    for (i in 1:J){
        if(is.na(specific[i])){
            for(j in 1:(nfact-1))
                r1[ ,itemloc[i]:(itemloc[i+1]-1)] <- r1[ ,itemloc[i]:(itemloc[i+1]-1)] + 	
                    retlist$r1[ ,itemloc[i]:(itemloc[i+1]-1) + (j - 1)*ncol(X) ]               
        } else {
            r1[ ,itemloc[i]:(itemloc[i+1]-1)] <- 		
                retlist$r1[ ,itemloc[i]:(itemloc[i+1]-1) + (specific[i] - 1)*ncol(X) ]		
        }
    }
    return(list(r1=r1, expected=retlist$expected))	
}      

Mstep.mirt <- function(par, obj, Theta, prior, constr = list(), debug){     
    if(debug == 'Mstep.mirt') browser()
    if(length(constr) < 1){
        obj@par[obj@est] <- par    
        ret <- LogLik(x=obj, Theta=Theta)                
    } else {        
        obj <- reloadConstr(par=par, constr=constr, obj=obj)        
        ret <- 0
        for(i in 1:length(obj))            
            ret <- ret + LogLik(x=obj[[i]], Theta=Theta)               
    }
    return(ret)
}

# Mstep for group pars
Mstep.group <- function(par, pars, gobj, Theta, tabdata, r, itemloc, constr = list(), debug)
{   
    if(debug == 'Mstep.group') browser()
    gobj@par[gobj@est] <- par    
    gpars <- ExtractGroupPars(gobj)
    mu <- gpars$gmeans
    sigma <- gpars$gcov    
    prior <- dmvnorm(Theta, mean=mu, sigma=sigma)
    prior <- prior/sum(prior)    
    rlist <- Estep.mirt(pars=pars, tabdata=tabdata, Theta=Theta, prior=prior, itemloc=itemloc, 
                        debug=debug)
    L <- sum(r*log(Pl))
    L   
}
