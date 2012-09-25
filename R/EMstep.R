EM <- function(pars, NCYCLES, MSTEPMAXIT, TOL, NULL.MODEL = FALSE, tabdata, tabdata2, Theta, 
               itemloc, debug, verbose, constrain, npars, data, sitems = NULL, specific = NULL, 
               theta = NULL, BFACTOR = FALSE, itemtype)
{
    if(debug == 'EM') browser()
    r <- tabdata[, ncol(tabdata)]    
    N <- sum(r)
    converge <- 1
    J <- length(itemloc) - 1
    listpars <- list()
    for(i in 1:length(pars))
        listpars[[i]] <- pars[[i]]@par
    lastpars2 <- lastpars1 <- listpars
    nfact <- pars[[length(pars)]]@nfact
    #initial priors
    if(BFACTOR){
        prior <- dnorm(theta)
        prior <- prior/sum(prior)
    } else {        
        gp <- ExtractGroupPars(pars[[length(pars)]])
        prior <- mvtnorm::dmvnorm(Theta,gp$gmeans,gp$gcov)
        prior <- prior/sum(prior)
    }
    prodlist <- attr(pars, 'prodlist')
    ThetaShort <- Theta
    if(length(prodlist) > 0)        
        Theta <- prodterms(Theta,prodlist)    
    
    #EM cycles
    for (cycles in 1:NCYCLES){        
        if(BFACTOR) 
            rlist <- Estep.bfactor(pars=pars, tabdata=tabdata, Theta=Theta, prior=prior,
                                   specific=specific, sitems=sitems, itemloc=itemloc, debug=debug)
        else
            rlist <- Estep.mirt(pars=pars, tabdata=tabdata, Theta=Theta, prior=prior, itemloc=itemloc, 
                                debug=debug)
        if(verbose){
            print(Pl <- sum(r*log(rlist$expected)))                            
            flush.console()
        }
        for(i in 1:J){
            tmp <- c(itemloc[i]:(itemloc[i+1] - 1))
            pars[[i]]@rs <- rlist$r1[, tmp]           
        }            
        lastpars2 <- lastpars1
        lastpars1 <- listpars
        #####
        #items without constraints
        for(i in 1:J){ 
            if(pars[[i]]@constr) next    	       
            estpar <- pars[[i]]@par[pars[[i]]@est]
            maxim <- try(optim(estpar, fn=Mstep.mirt, obj=pars[[i]], 
                               Theta=Theta, prior=prior, debug=debug,
                               method=pars[[i]]@method,
                               lower=pars[[i]]@lbound, 
                               upper=pars[[i]]@ubound,
                               control=list(maxit=MSTEPMAXIT)))
            if(class(maxim) == "try-error"){    			
                converge <- 0
                next
            }		  
            pars[[i]]@par[pars[[i]]@est] <- maxim$par            
        }
        #group      
        i = i + 1         
        estpar <- pars[[i]]@par[pars[[i]]@est]        
        if(length(estpar) > 0){
            maxim <- try(optim(estpar, fn=Mstep.group, 
                               pars=pars, gobj=pars[[i]], Theta=Theta, ThetaShort=ThetaShort,
                               tabdata=tabdata, r=r, itemloc=itemloc, constr=pars[[i]]@constr, 
                               debug=debug, 
                               method=pars[[i]]@method, 
                               lower=pars[[i]]@lbound, 
                               upper=pars[[i]]@ubound,
                               control=list(maxit=MSTEPMAXIT)))
            if(class(maxim) == "try-error"){
                converge <- 0
                next
            }		  
            pars[[i]]@par[pars[[i]]@est] <- maxim$par 
            gp <- ExtractGroupPars(pars[[i]])
            prior <- mvtnorm::dmvnorm(Theta,gp$gmeans,gp$gcov)
            prior <- prior/sum(prior)
        }
        ####    
        #items with constraints
        if(length(constrain) > 0){
            constrpars <- constrlist <- list()
            tmp <- 1
            for(i in 1:J){ 
                if(pars[[i]]@constr){
                    constrpars[[tmp]] <- pars[[i]] 
                    tmp <- tmp + 1
                }
            }
            tmp <- numpars <- c()
            for(i in 1:length(constrpars)){
                tmp <- c(tmp, constrpars[[i]]@par[constrpars[[i]]@est])
                numpars <- c(numpars, constrpars[[i]]@parnum[constrpars[[i]]@est])                
            }
            estpar <- c(rep(NA, length(constrain)), tmp[!(numpars %in% attr(pars, 'uniqueconstr'))])
            for(i in 1:length(constrain)){                
                constrlist[[i]] <- numpars %in% constrain[[i]]
                estpar[i] <- mean(tmp[constrlist[[i]]])
            }            
            maxim <- try(optim(estpar, fn=Mstep.mirt, obj=constrpars, debug=debug,
                               Theta=Theta, prior=prior, constr=constrlist,
                               method='Nelder-Mead',
                               lower=-Inf, 
                               upper=Inf,
                               control=list(maxit=MSTEPMAXIT)))            
            constrpars <- reloadConstr(maxim$par, constr=constrlist, obj=constrpars)
            tmp <- 1
            for(i in 1:J){ 
                if(pars[[i]]@constr){
                    pars[[i]] <- constrpars[[tmp]]
                    tmp <- tmp + 1
                }
            }
        }
        #apply sum(t) == 1 constraint for mcm
        if(is(pars[[i]], 'mcm')){
            tmp <- pars[[i]]@par
            tmp[length(tmp) - pars[[i]]@ncat + 1] <- 1 - sum(tmp[length(tmp):(length(tmp) - 
                pars[[i]]@ncat + 2)])
            pars[[i]]@par <- tmp
        }
        for(i in 1:J) listpars[[i]] <- pars[[i]]@par
        maxdif <- max(do.call(c,listpars) - do.call(c,lastpars1))
        if(maxdif < TOL && cycles > 10) break
        if(cycles %% 3 == 0 & cycles > 6)
            pars <- rateChange(pars=pars, listpars=listpars, lastpars1=lastpars1, 
                               lastpars2=lastpars2)
    }###END EM	    
    if(converge == 0) 
        warning("Parameter estimation reached unacceptable values. 
			Model probably did not converged.")  		
	lastchange <- do.call(c,listpars) - do.call(c,lastpars1)
    if(cycles == NCYCLES){
        converge <- 0  
        message("Estimation terminated after ", cycles, " EM loops and likely did not converge.")
    }      
    Pl <- rlist$expected  
    logLik <- sum(r*log(Pl))			
    logN <- 0
    npatmissing <- sum(is.na(rowSums(tabdata2)))
    logr <- rep(0,length(r))	
    for (i in 1:N) logN <- logN + log(i)
    for (i in 1:length(r)) 
        for (j in 1:r[i]) 
            logr[i] <- logr[i] + log(j)    	
    nconstr <- 0
    if(length(constrain) > 0)
        for(i in 1:length(constrain))
            nconstr <- nconstr + length(constrain[[i]]) - 1
    df <- (length(r) - 1) - npars + nfact*(nfact - 1)/2  - npatmissing + nconstr	
    if(NULL.MODEL) df <- (length(r) - 1) - npars - npatmissing
    G2 <- 2 * sum(r * log(r/(N*Pl)))	
    logLik <- logLik + logN/sum(logr)	
    p <- 1 - pchisq(G2,df)  
    AIC <- (-2) * logLik + 2 * npars
    BIC <- (-2) * logLik + npars*log(N)
    RMSEA <- ifelse((G2 - df) > 0, 
                    sqrt(G2 - df) / sqrt(df * (N-1)), 0)	    	
    null.mod <- unclass(new('ExploratoryClass'))
    if(!NULL.MODEL)        
        null.mod <- unclass(mirt(data, 1, itemtype=itemtype, technical = list(NULL.MODEL = TRUE), 
                                 SE = FALSE))    
    TLI <- NaN    
    if(!NULL.MODEL)
        TLI <- (null.mod@G2 / null.mod@df - G2/df) / (null.mod@G2 / null.mod@df - 1)
    if(npatmissing > 0) p <- RMSEA <- G2 <- TLI <- NaN
    ret <- list(pars=pars, converge=converge, cycles=cycles, maxdif=maxdif, G2=G2, df=df, p=p, TLI=TLI,
                AIC=AIC, BIC=BIC, RMSEA=RMSEA, null.mod=null.mod, logLik=logLik, Pl=Pl)
    ret
}

# Estep for mirt
Estep.mirt <- function(pars, tabdata, Theta, prior, itemloc, debug, deriv = FALSE) 
{   
    if(debug == 'Estep.mirt') browser()
    nfact <- ncol(Theta)
    nquad <- nrow(Theta)	
    J <- length(itemloc) - 1
    r <- tabdata[ ,ncol(tabdata)]
    X <- tabdata[ ,1:(ncol(tabdata) - 1)]	
    itemtrace <- matrix(0, ncol=ncol(X), nrow=nrow(Theta))	
    for (i in 1:J)
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)
    retlist <- .Call("Estep", itemtrace, prior, X, nfact, r)    
    if(deriv) retlist$itemtrace <- itemtrace        
    return(retlist)
} 

# Estep for bfactor
Estep.bfactor <- function(pars, tabdata, Theta, prior, specific, sitems, itemloc, debug) 
{	    
    if(debug == 'Estep.bfactor') browser()
    J <- length(itemloc) - 1
    r <- tabdata[ ,ncol(tabdata)]
    X <- tabdata[ ,1:(ncol(tabdata) - 1)]	
    itemtrace <- matrix(0, ncol=ncol(X), nrow=nrow(Theta))	
    for (i in 1:J)
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)			
    retlist <- .Call("Estepbfactor", itemtrace, prior, X, r, sitems)	
    r1 <- matrix(0, nrow(Theta), ncol(X))	
    for (i in 1:J){
        r1[ ,itemloc[i]:(itemloc[i+1]-1)] <- 		
            retlist$r1[ ,itemloc[i]:(itemloc[i+1]-1) + (specific[i] - 1)*ncol(X) ]		        
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
Mstep.group <- function(par, pars, gobj, Theta, ThetaShort, tabdata, r, itemloc, constr = list(), debug)
{   
    if(debug == 'Mstep.group') browser()
    gobj@par[gobj@est] <- par    
    gpars <- ExtractGroupPars(gobj)
    mu <- gpars$gmeans
    sigma <- gpars$gcov
    prior <- dmvnorm(ThetaShort, mean=mu, sigma=sigma)
    prior <- prior/sum(prior)
    rlist <- Estep.mirt(pars=pars, tabdata=tabdata, Theta=Theta, prior=prior, itemloc=itemloc, 
                        debug=debug)
    L <- (-1)*sum(r*log(rlist$expected))
    L   
}
