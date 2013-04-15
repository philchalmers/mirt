LoadPars <- function(itemtype, itemloc, lambdas, zetas, guess, upper, fulldata, J, K, nfact, 
                     parprior, parnumber, D, estLambdas, BFACTOR = FALSE, mixedlist, customItems, 
                     key)
    {           
    valid.items <- c('Rasch', '1PL', '2PL', '3PL', '3PLu', '4PL', 'graded', 
                    'grsm', 'gpcm', 'rsm', 'nominal', 'mcm', 'PC2PL','PC3PL',
                    '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')
    pars <- vector('list', J)
    #startvalues
    startvalues <- vector('list', J)
    for(i in 1:J){            
        if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] == 2){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1/D                
            val <- c(tmpval, zetas[[i]], guess[i], upper[i])
            names(val) <- c(paste('a', 1:nfact, sep=''), 'd', 'g','u')
        }
        if(itemtype[i] == 'Rasch' && K[i] > 2){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1/D                
            val <- c(tmpval, 0, zetas[[i]])
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 0:(K[i]-1), sep=''))
        }
        if(itemtype[i] == '1PL' && K[i] > 2){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1/D                
            val <- c(tmpval, zetas[[i]])
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:(K[i]-1), sep=''))
        }
        if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){
            val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i])
            names(val) <- c(paste('a', 1:nfact, sep=''), 'd', 'g','u')
        }
        if(any(itemtype[i] == c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM'))){               
            val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i],
                     0, rep(.5, K[i] - 2), rep(0, K[i]-1))  
            names(val) <- c(paste('a', 1:nfact, sep=''), 'd', 'g','u', 
                            paste('ak', 0:(K[i]-2), sep=''), 
                            paste('d', 0:(K[i]-2), sep=''))
        }
        if(itemtype[i] == 'graded'){
            val <- c(lambdas[i,], zetas[[i]])
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:(K[i]-1), sep=''))    
        }            
        if(itemtype[i] == 'grsm'){
            val <- c(lambdas[i,], seq(2.5, -2.5, length.out = length(zetas[[i]])), 0) 
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:(K[i]-1), sep=''), 'c')
        }
        if(itemtype[i] == 'gpcm'){
            val <- c(lambdas[i,], 0, zetas[[i]])
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 0:(K[i]-1), sep=''))                
        }
        if(itemtype[i] == 'rsm'){                
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1/D                
            val <- c(tmpval, 0, seq(2.5, -2.5, length.out = length(zetas[[i]])), 0)
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 0:(K[i]-1), sep=''), 'c')                
        }
        if(itemtype[i] == 'nominal'){
            val <- c(lambdas[i,], 0, rep(.5, K[i] - 2), K[i]-1, rep(0, K[i]))
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('ak', 0:(K[i]-1), sep=''), 
                            paste('d', 0:(K[i]-1), sep=''))                
        }
        if(any(itemtype[i] == c('PC2PL','PC3PL'))){
            val <- c(lambdas[i,], rep(1, nfact), 0, 1)
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:nfact, sep=''), 'g','u')
        }
        if(itemtype[i] == 'mcm'){
            val <- c(lambdas[i,], 1, 0, rep(.5, K[i] - 2), K[i]-1, rep(0, K[i]+1), 
                     rep(1/K[i], K[i]))
            names(val) <- c(paste('a', 1:nfact, sep=''), paste('ak', 0:(K[i]), sep=''), 
                            paste('d', 0:(K[i]), sep=''), paste('t', 1:(K[i]), sep=''))                
        }            
        if(all(itemtype[i] != valid.items)) next            
        startvalues[[i]] <- val
    } 
    #freepars
    freepars <- vector('list', J)
    for(i in 1:J){                        
        if(itemtype[i] == 'Rasch' && K[i] == 2)
            freepars[[i]] <- c(rep(FALSE,nfact),TRUE,FALSE,FALSE)            
        if(any(itemtype[i] == c('1PL', '2PL', '3PL', '3PLu', '4PL'))){
            estpars <- c(estLambdas[i, ], TRUE, FALSE, FALSE) 
            if(any(itemtype[i] == c('3PL', '4PL'))) estpars[length(estpars)-1] <- TRUE
            if(any(itemtype[i] == c('3PLu', '4PL'))) estpars[length(estpars)] <- TRUE
            freepars[[i]] <- estpars
        }
        if(any(itemtype[i] == c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM'))){               
            estpars <- c(estLambdas[i, ], TRUE, FALSE, FALSE, rep(TRUE, (K[i]-1)*2))             
            estpars[c(nfact+4, length(estpars)-(K[i]-2) )] <- FALSE
            if(any(itemtype[i] == c('3PLNRM', '4PLNRM'))) estpars[nfact+2] <- TRUE
            if(any(itemtype[i] == c('3PLuNRM', '4PLNRM'))) estpars[nfact+3] <- TRUE
            freepars[[i]] <- estpars
        }
        if(itemtype[i] == 'Rasch' && K[i] > 2)            
            freepars[[i]] <- c(rep(FALSE,nfact), FALSE, rep(TRUE, K[i]-1))
        if(itemtype[i] == '1PL' && K[i] > 2)            
            freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]))            
        if(itemtype[i] == 'grsm')
            freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]))
        if(itemtype[i] == 'graded')
            freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]-1))
        if(itemtype[i] == 'gpcm')            
            freepars[[i]] <- c(estLambdas[i, ], FALSE, rep(TRUE, K[i]-1))  
        if(itemtype[i] == 'rsm')            
            freepars[[i]] <- c(rep(FALSE, nfact), FALSE, rep(TRUE, K[i]))
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
            estpars <- c(estLambdas[i, ], rep(TRUE, 2+K[i]*3))
            #identifiction constraints
            tmp <- names(startvalues[[i]])
            tmp2 <- 1:length(tmp)            
            estpars[tmp2[tmp %in% c('ak0', 'ak1', 
                                    paste('ak',K[i],sep=''), 'd0', 'd1', 't1')]] <- FALSE                
            freepars[[i]] <- estpars
        }
        if(all(itemtype[i] != valid.items)) next
        names(freepars[[i]]) <- names(startvalues[[i]])    
    }         

    #augment startvalues and fixedpars for mixed effects
    nfixedeffects <- 0    
    if(!is.null(mixedlist)){ 
        fixed.design.list <- designMats(covdata=mixedlist$covdata, fixed=mixedlist$fixed, 
                                        Thetas=mixedlist$Theta, nitems=J, 
                                        itemdesign=mixedlist$itemdesign, 
                                        fixed.identical=mixedlist$fixed.identical)
        betas <- rep(0, ncol(fixed.design.list[[1]]))
        estbetas <- rep(TRUE, length(betas))
        names(estbetas) <- names(betas) <- colnames(fixed.design.list[[1]])
        nfixedeffects <- length(betas)
        nfact <- nfact + nfixedeffects
        for(i in 1:J){
            freepars[[i]] <- c(estbetas, freepars[[i]])
            startvalues[[i]] <- c(betas, startvalues[[i]])            
        }        
    }    
    
    #load items
    for(i in 1:J){
        tmp <- c(itemloc[i]:(itemloc[i+1] - 1)) #item location         
        
        if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] == 2){ 
            pars[[i]] <- new('dich', par=startvalues[[i]], est=freepars[[i]],
                             nfact=nfact, 
                             dat=fulldata[ ,tmp], 
                             ncat=2,
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             lbound=c(rep(-Inf, length(startvalues[[i]]) - 2),0,.5),
                             ubound=c(rep(Inf, length(startvalues[[i]]) - 2),.5,1),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))       
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])  
            next
        }

        if(itemtype[i] == 'Rasch' && K[i] > 2){ 
            pars[[i]] <- new('gpcm', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            pars[[i]]@par[nfact+1] <- 0            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]]) 
            next
        }
        
        if(itemtype[i] == '1PL' && K[i] > 2){
            pars[[i]] <- new('graded', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }
        
        if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){ 
            pars[[i]] <- new('dich', 
                             par=startvalues[[i]], 
                             est=freepars[[i]],
                             nfact=nfact, 
                             nfixedeffects=nfixedeffects, 
                             dat=fulldata[ ,tmp], 
                             ncat=2,
                             D=D,
                             lbound=c(rep(-Inf, length(startvalues[[i]]) - 2),0,.5),                                           
                             ubound=c(rep(Inf, length(startvalues[[i]]) - 2),.5,1),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])  
            next
        }
        
        if(any(itemtype[i] == c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM'))){             
            pars[[i]] <- new('nestlogit', 
                             par=startvalues[[i]], 
                             est=freepars[[i]],
                             nfact=nfact, 
                             nfixedeffects=nfixedeffects, 
                             dat=fulldata[ ,tmp], 
                             ncat=K[i],
                             correctcat=key[i],
                             D=D,
                             lbound=c(rep(-Inf, nfact+1),0,.5, rep(-Inf, length(startvalues[[i]])-nfact-3)),
                             ubound=c(rep(Inf, nfact+1),.5,1, rep(Inf, length(startvalues[[i]])-nfact-3)),
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])  
            next
        }
        
        if(any(itemtype[i] == 'grsm')){
            pars[[i]] <- new('rating', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }
        
        if(itemtype[i] == 'graded'){
            pars[[i]] <- new('graded', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }
        
        if(itemtype[i] == 'gpcm'){            
            pars[[i]] <- new('gpcm', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            pars[[i]]@par[nfact+1] <- 0            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }    
        
        if(itemtype[i] == 'rsm'){            
            pars[[i]] <- new('rsm', 
                             par=startvalues[[i]], 
                             nfact=nfact, 
                             ncat=K[i],
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                   
            pars[[i]]@par[nfact+1] <- 0            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }   
        
        if(itemtype[i] == 'nominal'){
            pars[[i]] <- new('nominal', 
                             par=startvalues[[i]], 
                             est=freepars[[i]], 
                             nfact=nfact, 
                             ncat=K[i], 
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             dat=fulldata[ ,tmp], 
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))       
            pars[[i]]@par[c(nfact + 1, nfact + K[i] + 1)] <- 0
            pars[[i]]@par[nfact + K[i]] <- K[i] - 1            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        } 
        
        if(any(itemtype[i] == c('PC2PL','PC3PL'))){
            pars[[i]] <- new('partcomp', 
                             par=startvalues[[i]], 
                             est=freepars[[i]],
                             nfact=nfact, 
                             ncat=2, 
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             dat=fulldata[ ,tmp], 
                             lbound=c(rep(-Inf, length(startvalues[[i]]) - 2),0,.5),
                             ubound=c(rep(Inf, length(startvalues[[i]]) - 2),.5,1),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }
        
        if(itemtype[i] == 'mcm'){
            pars[[i]] <- new('mcm', 
                             par=startvalues[[i]], 
                             est=freepars[[i]], 
                             nfact=nfact, 
                             ncat=K[i], 
                             nfixedeffects=nfixedeffects, 
                             D=D,
                             dat=fulldata[ ,tmp], 
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),                             
                             n.prior.mu=rep(NaN,length(startvalues[[i]])),
                             n.prior.sd=rep(NaN,length(startvalues[[i]])),
                             b.prior.alpha=rep(NaN,length(startvalues[[i]])),
                             b.prior.beta=rep(NaN,length(startvalues[[i]])))                            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }
        
        if(all(itemtype[i] != valid.items)){                  
            pars[[i]] <- customItems[[itemtype[i] == names(customItems)]]            
            pars[[i]]@nfact <- nfact 
            pars[[i]]@ncat <- K[i] 
            pars[[i]]@nfixedeffects <- nfixedeffects 
            pars[[i]]@D <- D
            pars[[i]]@dat <- fulldata[ ,tmp]
            pars[[i]]@n.prior.mu <- rep(NaN,length(pars[[i]]@par))
            pars[[i]]@n.prior.sd <- rep(NaN,length(pars[[i]]@par))
            pars[[i]]@b.prior.alpha <- rep(NaN,length(pars[[i]]@par))
            pars[[i]]@b.prior.beta <- rep(NaN,length(pars[[i]]@par))
            tmp2 <- parnumber:(parnumber + length(pars[[i]]@est) - 1)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(pars[[i]]@est)
            next
        }
    }   
    #priors
    for(i in 1:J){
        names(pars[[i]]@parnum) <- names(startvalues[[i]])
        if(!is.null(parprior) && parprior != 'index'){
            for(j in 1:length(parprior)){                                
                tmp <- pars[[i]]@parnum %in% as.numeric(parprior[[j]][1:(length(parprior[[j]])-3)])
                if(any(tmp)){
                    if(parprior[[j]][length(parprior[[j]]) - 2] == 'norm'){
                        pars[[i]]@n.prior.mu[tmp] <- as.numeric(parprior[[j]][length(parprior[[j]])-1])
                        pars[[i]]@n.prior.sd[tmp] <- as.numeric(parprior[[j]][length(parprior[[j]])])
                    } else if(parprior[[j]][length(parprior[[j]]) - 2] == 'beta'){
                        pars[[i]]@b.prior.alpha[tmp] <- as.numeric(parprior[[j]][length(parprior[[j]])-1])
                        pars[[i]]@b.prior.beta[tmp] <- as.numeric(parprior[[j]][length(parprior[[j]])])
                    }                
                }          
            }
        }
    } 
    attr(pars, 'parnumber') <- parnumber     
    return(pars)
}

LoadGroupPars <- function(gmeans, gcov, estgmeans, estgcov, parnumber, parprior){    
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
               parnum=parnum, lbound=rep(-Inf, length(par)), ubound=rep(Inf, length(par)))                      
    if(!is.null(parprior) && is.list(parprior)){
        for(j in 1:length(parprior)){            
            tmp <- parnum %in% as.numeric(parprior[[j]][1:(length(parprior[[j]])-3)])
            if(any(tmp)){
                if(parprior[[j]][length(parprior[[j]]) - 2] == 'norm'){
                    ret@n.prior.mu[tmp] <- as.numeric(parprior[[j]][length(parprior[[j]])-1])
                    ret@n.prior.sd[tmp] <- as.numeric(parprior[[j]][length(parprior[[j]])])
                } else if(parprior[[j]][length(parprior[[j]]) - 2] == 'beta'){
                    ret@b.prior.alpha[tmp] <- as.numeric(parprior[[j]][length(parprior[[j]])-1])
                    ret@b.prior.beta[tmp] <- as.numeric(parprior[[j]][length(parprior[[j]])])
                }                
            }          
        }    
    }    
    return(ret)    
}
