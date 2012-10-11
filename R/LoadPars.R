LoadPars <- function(itemtype, itemloc, lambdas, zetas, guess, upper, fulldata, J, K, nfact, 
                     constrain, startvalues, freepars, parprior, parnumber, 
                     estLambdas, BFACTOR = FALSE, debug)
    {       
    if(debug == 'LoadPars') browser() 
    if(any(itemtype[1] == c('Rasch', '1PL', 'rating') && nfact > 1)) 
        stop('Rasch, 1PL, and rating scale models can only be estimated for unidimensional models')
    pars <- list()           
    constr <- c()
    if(!is.null(constrain) && is.list(constrain)) 
        for(i in 1:length(constrain))
            constr <- c(constr, constrain[[i]])
    constr <- unique(constr)
    #startvalues
    if(is.null(startvalues) || startvalues =='index'){        
        startvalues <- list()
        for(i in 1:J){            
            if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] == 2){
                val <- c(1/1.702, zetas[[i]], guess[i], upper[i])
                names(val) <- c(paste('a', 1:nfact, sep=''), 'd', 'g','u')
            }
            if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] > 2){
                val <- c(1/1.702, zetas[[i]])
                names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 0:(K[i]-1), sep=''))
            }
            if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){
                val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i])
                names(val) <- c(paste('a', 1:nfact, sep=''), 'd', 'g','u')
            }
            if(itemtype[i] == 'graded'){
                val <- c(lambdas[i,], zetas[[i]])
                names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:(K[i]-1), sep=''))    
            }            
            if(itemtype[i] == 'grsm'){
                val <- c(lambdas[i,], zetas[[1]], 0) #first item intercepts
                names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:(K[i]-1), sep=''), 'c')
            }
            if(itemtype[i] == 'gpcm'){
                val <- c(lambdas[i,], 0, zetas[[i]])
                names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 0:(K[i]-1), sep=''))                
            }
            if(itemtype[i] == 'nominal'){
                val <- c(lambdas[i,], 0, rep(.5, K[i] - 2), K[i]-1, rep(0, K[i]))
                names(val) <- c(paste('a', 1:nfact, sep=''), paste('ak', 0:(K[i]-1), sep=''), 
                                paste('d', 0:(K[i]-1), sep=''))                
            }
            if(any(itemtype[i] == c('PC2PL','PC3PL'))){
                val <- c(lambdas[i,], rep(-1, nfact), 0, 1)
                names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:nfact, sep=''), 'g','u')
            }
            if(itemtype[i] == 'mcm'){
                val <- c(lambdas[i,], 0, rep(.5, K[i] - 2), K[i]-1, rep(0, K[i]), 
                         rep(1/K[i], K[i]))
                names(val) <- c(paste('a', 1:nfact, sep=''), paste('ak', 0:(K[i]-1), sep=''), 
                                paste('d', 0:(K[i]-1), sep=''), paste('t', 0:(K[i]-1), sep=''))                
            }            
            startvalues[[i]] <- val
        } 
    }  
    #freepars
    if(is.null(freepars) || freepars == 'index'){
        freepars <- list()
        for(i in 1:J){                        
            if(itemtype[i] == 'Rasch' && K[i] == 2)
                freepars[[i]] <- c(FALSE,TRUE,FALSE,FALSE)            
            if(any(itemtype[i] == c('1PL', '2PL', '3PL', '3PLu', '4PL'))){
                estpars <- c(estLambdas[i, ], TRUE, FALSE, FALSE) 
                if(any(itemtype[i] == c('3PL', '4PL'))) estpars[length(estpars)-1] <- TRUE
                if(any(itemtype[i] == c('3PLu', '4PL'))) estpars[length(estpars)] <- TRUE
                freepars[[i]] <- estpars
            }
            if(itemtype[i] == 'Rasch' && K[i] > 2)            
                freepars[[i]] <- c(FALSE, rep(TRUE, K[i]))
            if(itemtype[i] == '1PL' && K[i] > 2)            
                freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]))            
            if(itemtype[i] == 'grsm')
                freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]))
            if(itemtype[i] == 'graded')
                freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]-1))
            if(itemtype[i] == 'gpcm')            
                freepars[[i]] <- c(estLambdas[i, ], FALSE, rep(TRUE, K[i]-1))            
            if(itemtype[i] == 'nominal'){
                estpars <- c(estLambdas[i, ], rep(TRUE, K[i]*2))
                #identifiction constraints
                estpars[c(nfact+1, nfact + K[i], nfact + K[i] + 1)] <- FALSE
                freepars[[i]] <- estpars
            }
            if(any(itemtype[i] == c('PC2PL','PC3PL'))){
                estpars <- c(estLambdas[i, ], estLambdas[i, ], FALSE, FALSE)
                if(itemtype[i] == 'PC3PL') estpars[length(estpars) - 1] <- TRUE
                freepars[[i]] <- estpars
            }
            if(itemtype[i] == 'mcm'){
                estpars <- c(estLambdas[i, ], rep(TRUE, K[i]*3))
                #identifiction constraints
                estpars[c(nfact+1, nfact + K[i], nfact + K[i] + 1, length(estpars) - K[i] + 1)] <- FALSE
                freepars[[i]] <- estpars
            }
        }         
    }
    for(i in 1:J) names(freepars[[i]]) <- names(startvalues[[i]])    
    #load items
    for(i in 1:J){
        tmp <- c(itemloc[i]:(itemloc[i+1] - 1)) #item location         
        
        if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] == 2){ 
            pars[[i]] <- new('dich', par=startvalues[[i]], est=freepars[[i]],
                             nfact=nfact, 
                             dat=fulldata[ ,tmp], 
                             constr=TRUE,                              
                             ncat=2,
                             lbound=ifelse(itemtype[i] == 'Rasch', -25, -Inf),
                             ubound=ifelse(itemtype[i] == 'Rasch', 25, Inf),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))       
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])            
        }

        if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] > 2){ 
            pars[[i]] <- new('gpcm', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             constr=TRUE,                              
                             ncat=K[i],
                             lbound=-Inf,
                             ubound=Inf,                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            pars[[i]]@par[nfact+1] <- 0            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])            
        }
        
        if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){ 
            pars[[i]] <- new('dich', 
                             par=startvalues[[i]], 
                             est=freepars[[i]],
                             nfact=nfact, 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE,                              
                             ncat=2,
                             lbound=-Inf,                                           
                             ubound=Inf,                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))       
            if(itemtype[i] != '2PL'){       
                tmp2 <- c(rep(-Inf, length(startvalues[[i]]) - 2),0,0)
                tmp3 <- c(rep(Inf, length(startvalues[[i]]) - 2),1,1)
                pars[[i]]@lbound <- tmp2[freepars[[i]]]
                pars[[i]]@ubound <- tmp3[freepars[[i]]]
            }
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])            
        }
        
        if(any(itemtype[i] == 'grsm')){
            pars[[i]] <- new('rating', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             constr=TRUE,                              
                             lbound=-Inf,
                             ubound=Inf,                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
        }
        
        if(itemtype[i] == 'graded'){
            pars[[i]] <- new('graded', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE,                              
                             lbound=-Inf,
                             ubound=Inf,                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
        }
        
        if(itemtype[i] == 'gpcm'){            
            pars[[i]] <- new('gpcm', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE,                              
                             lbound=-Inf,
                             ubound=Inf,                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            pars[[i]]@par[nfact+1] <- 0            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
        }        
        
        if(itemtype[i] == 'nominal'){
            pars[[i]] <- new('nominal', 
                             par=startvalues[[i]], 
                             est=freepars[[i]], 
                             nfact=nfact, 
                             ncat=K[i], 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE,                              
                             lbound=-Inf,
                             ubound=Inf,                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))       
            pars[[i]]@par[c(nfact + 1, nfact + K[i] + 1)] <- 0
            pars[[i]]@par[nfact + K[i]] <- K[i] - 1            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
        } 
        
        if(any(itemtype[i] == c('PC2PL','PC3PL'))){
            pars[[i]] <- new('partcomp', 
                             par=startvalues[[i]], 
                             est=freepars[[i]],
                             nfact=nfact, 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE,                              
                             lbound=-Inf,
                             ubound=Inf,                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))
            if(itemtype[i] != 'PC2PL'){       
                tmp2 <- c(rep(-Inf, length(startvalues[[i]]) - 2),0,0)
                tmp3 <- c(rep(Inf, length(startvalues[[i]]) - 2),1,1)
                pars[[i]]@lbound <- tmp2[freepars[[i]]]
                pars[[i]]@ubound <- tmp3[freepars[[i]]]
            }
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
        }
        
        if(itemtype[i] == 'mcm'){
            pars[[i]] <- new('mcm', 
                             par=startvalues[[i]], 
                             est=freepars[[i]], 
                             nfact=nfact, 
                             ncat=K[i], 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE,                              
                             lbound=-Inf,
                             ubound=Inf,                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
        }
    }   
    #priors
    for(i in 1:J){
        names(pars[[i]]@parnum) <- names(startvalues[[i]])
        if(!is.null(parprior) && parprior != 'index'){
            for(j in 1:length(parprior)){
                tmp <- pars[[i]]@parnum %in% as.numeric(parprior[[j]][1])
                if(any(tmp)){
                    if(parprior[[j]][2] == 'norm'){
                        pars[[i]]@n.prior.mu[tmp] <- as.numeric(parprior[[j]][3])
                        pars[[i]]@n.prior.sd[tmp] <- as.numeric(parprior[[j]][4])
                    } else {
                        pars[[i]]@b.prior.alpha[tmp] <- as.numeric(parprior[[j]][3])
                        pars[[i]]@b.prior.beta[tmp] <- as.numeric(parprior[[j]][4])
                    }                
                }          
            }
        }
    } 
    attr(pars, 'uniqueconstr') <- constr     
    attr(pars, 'parnumber') <- attr(startvalues, 'parnumber') <- 
        attr(freepars, 'parnumber') <- parnumber     
    return(pars)
}

LoadGroupPars <- function(gmeans, gcov, estgmeans, estgcov, parnumber, constrain, parprior, startvalues,
                          freepars, debug){
    if(debug == 'LoadGroupPars') browser()
    nfact <- length(gmeans)
    fn <- paste('COV_', 1:nfact, sep='')
    FNCOV <- outer(fn, 1:nfact, FUN=paste, sep='')
    FNMEANS <- paste('MEAN_', 1:nfact, sep='')      
    tri <- lower.tri(gcov, diag=TRUE)
    par <- c(gmeans, gcov[tri])
    parnum <- parnumber:(parnumber + length(par) - 1)
    est <- c(estgmeans,estgcov[tri])
    names(parnum) <- names(par) <- names(est) <- c(FNMEANS,FNCOV[tri])     
    ret <- new('GroupPars', par=par, est=est, nfact=nfact, 
               parnum=parnum, lbound=-Inf, ubound=Inf)                      
    if(!is.null(startvalues)){
        if(is.list(startvalues)) ret@par <- startvalues[[length(startvalues)]]
        else return(ret@par)        
    }
    if(!is.null(freepars)){
        if(is.list(freepars)) ret@est <- freepars[[length(freepars)]]
        else return(ret@est)        
    }
    if(!is.null(parprior) && is.list(parprior)){
        for(j in 1:length(parprior)){
            tmp <- parnum %in% as.numeric(parprior[[j]][1])
            if(any(tmp)){
                if(parprior[[j]][2] == 'norm'){
                    ret@n.prior.mu[tmp] <- as.numeric(parprior[[j]][3])
                    ret@n.prior.sd[tmp] <- as.numeric(parprior[[j]][4])
                } else {
                    ret@b.prior.alpha[tmp] <- as.numeric(parprior[[j]][3])
                    ret@b.prior.beta[tmp] <- as.numeric(parprior[[j]][4])
                }                
            }          
        }    
    }
    return(ret)    
}
