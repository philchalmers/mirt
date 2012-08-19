LoadPars <- function(itemtype, itemloc, lambdas, zetas, guess, upper, fulldata, J, K, nfact, 
                     constrain, startvalues, freepars, parprior, parnumber, 
                     estLambdas, BFACTOR = FALSE, nfactNames = NULL, debug)
    {       
    if(debug == 'LoadPars') browser() 
    if(any(itemtype[1] == c('Rasch', '1PL') && nfact > 1)) 
        stop('Rasch and 1PL models can only be estimated for unidimensional models')
    pars <- list()       
    RETURNSTARTVALUES <- ifelse(!is.null(startvalues) && startvalues == 'index', TRUE, FALSE)
    RETURNFREEPARS <- ifelse(!is.null(freepars) && freepars == 'index', TRUE, FALSE)
    if(is.null(nfactNames)) nfactNames <- nfact        
    constr <- c()
    if(!is.null(constrain) && is.list(constrain)) 
        for(i in 1:length(constrain))
            constr <- c(constr, constrain[[i]])
    constr <- unique(constr)
    if(is.null(startvalues) || startvalues =='index'){        
        startvalues <- list()
        for(i in 1:J){
            if(itemtype[i] == 'NullModel' && K[i] == 2) val <- c(0,zetas[[i]],0,1)                                
            if(itemtype[i] == 'NullModel' && K[i] > 2) val <- c(0,zetas[[i]]) 
            if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] == 2){
                val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i])
                names(val) <- c(paste('a', 1:nfactNames, sep=''), 'd', 'g','u')
            }
            if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] > 2){
                val <- c(lambdas[i,], zetas[[i]])
                names(val) <- c(paste('a', 1:nfactNames, sep=''), paste('d', 0:(K[i]-1), sep=''))
            }
            if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){
                val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i])
                names(val) <- c(paste('a', 1:nfactNames, sep=''), 'd', 'g','u')
            }
            if(itemtype[i] == 'graded'){
                val <- c(lambdas[i,], zetas[[i]])
                names(val) <- c(paste('a', 1:nfactNames, sep=''), paste('d', 1:(K[i]-1), sep=''))    
            }
            if(itemtype[i] == 'gpcm'){
                val <- c(lambdas[i,], zetas[[i]])
                names(val) <- c(paste('a', 1:nfactNames, sep=''), paste('d', 0:(K[i]-1), sep=''))                
            }
            if(itemtype[i] == 'nominal'){
                val <- c(rep(.5, nfactNames), 0, rep(.5, K[i] - 2), K[i]-1, rep(0, K[i]))
                names(val) <- c(paste('a', 1:nfactNames, sep=''), paste('ak', 0:(K[i]-1), sep=''), 
                                paste('d', 0:(K[i]-1), sep=''))                
            }
            if(any(itemtype[i] == c('PC2PL','PC3PL'))){
                val <- c(rep(.5, nfact), rep(-1, nfact), 0, 1)
                names(val) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:nfact, sep=''), 'g','u')
            }
            startvalues[[i]] <- val
        } 
    }        
    if(is.null(freepars) || freepars == 'index'){
        freepars <- list()
        for(i in 1:J){
            if(itemtype[i] == 'NullModel' && K[i] == 2)
                freepars[[i]] <- c(FALSE,TRUE,FALSE,FALSE)
            if(itemtype[i] == 'NullModel' && K[i] > 2)    
                freepars[[i]] <- c(FALSE,rep(TRUE,K[i]-1))            
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
            if(itemtype[i] == 'graded')
                freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]-1))
            if(itemtype[i] == 'gpcm')            
                freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]))            
            if(itemtype[i] == 'nominal'){
                estpars <- c(estLambdas[i, ], rep(TRUE, length(pars[[i]]@par) - nfact))
                #identifiction constraints
                estpars[c(nfact+1, nfact+ K[i], nfact + K[i] + 1)] <- FALSE
                freepars[[i]] <- estpars
            }
            if(any(itemtype[i] == c('PC2PL','PC3PL'))){
                estpars <- c(estLambdas[i, ], estLambdas[i, ], FALSE, FALSE)
                if(itemtype[i] == 'PC3PL') estpars[length(estpars) - 1] <- TRUE
                freepars[[i]] <- estpars
            }
        }         
    }
    for(i in 1:J) names(freepars[[i]]) <- names(startvalues[[i]])
    if(itemtype[1] == 'Rasch') 
        for(i in 1:J)
            startvalues[[i]][1] <- 1/1.702            
    for(i in 1:J){
        tmp <- c(itemloc[i]:(itemloc[i+1] - 1)) #item location 
        if(itemtype[i] == 'NullModel' && K[i] == 2){ 
            pars[[i]] <- new('dich', 
                             par=startvalues[[i]], nfact=1, 
                             bfactor=BFACTOR,
                             dat=fulldata[ ,tmp], 
                             est=freepars[[i]], 
                             constr=FALSE,
                             lbound=-25,
                             ubound=25,
                             method='Brent')
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])            
        }
        
        if(itemtype[i] == 'NullModel' && K[i] > 2){ 
            pars[[i]] <- new('graded', 
                             par=startvalues[[i]], 
                             nfact=1, 
                             ncat=K[i], 
                             bfactor=BFACTOR,
                             dat=fulldata[ ,tmp], 
                             est=freepars[[i]], 
                             constr=FALSE,
                             lbound=-Inf,
                             ubound=Inf,
                             method='Nelder-Mead')
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
        }
        
        if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] == 2){ 
            pars[[i]] <- new('dich', par=startvalues[[i]], est=freepars[[i]],
                             nfact=nfactNames, 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE, 
                             bfactor=BFACTOR,
                             lbound=-Inf,
                             ubound=Inf,
                             method='Nelder-Mead')            
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])            
        }

        if(any(itemtype[i] == c('Rasch', '1PL')) && K[i] > 2){ 
            pars[[i]] <- new('gpcm', 
                             par=startvalues[[i]], 
                             nfact=nfactNames, 
                             ncat=K[i],
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE, 
                             bfactor=BFACTOR,
                             lbound=-Inf,
                             ubound=Inf,
                             method='Nelder-Mead')                        
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
                             nfact=nfactNames, 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE, 
                             bfactor=BFACTOR,
                             lbound=-Inf,                                           
                             ubound=Inf,
                             method=ifelse(itemtype[i] == '2PL', 'Nelder-Mead', 'L-BFGS-B'))            
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
        
        if(itemtype[i] == 'graded'){
            pars[[i]] <- new('graded', 
                             par=startvalues[[i]], 
                             nfact=nfactNames, 
                             ncat=K[i],
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE, 
                             bfactor=BFACTOR,
                             lbound=-Inf,
                             ubound=Inf,
                             method='Nelder-Mead')                        
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE            
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
        }
        
        if(itemtype[i] == 'gpcm'){            
            pars[[i]] <- new('gpcm', 
                             par=startvalues[[i]], 
                             nfact=nfactNames, 
                             ncat=K[i],
                             est=freepars[[i]], 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE, 
                             bfactor=BFACTOR,
                             lbound=-Inf,
                             ubound=Inf,
                             method='Nelder-Mead')                        
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
                             nfact=nfactNames, 
                             ncat=K[i], 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE, 
                             bfactor=BFACTOR,
                             lbound=-Inf,
                             ubound=Inf,
                             method='Nelder-Mead')            
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
                             nfact=nfactNames, 
                             dat=fulldata[ ,tmp], 
                             constr=FALSE, 
                             bfactor=BFACTOR,
                             lbound=-Inf,
                             ubound=Inf,
                             method=ifelse(itemtype[i] == 'PC2PL', 'Nelder-Mead', 'L-BFGS-B'))
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
    }
    for(i in 1:J) names(pars[[i]]@parnum) <- names(startvalues[[i]])
    attr(pars, 'uniqueconstr') <- constr     
    attr(pars, 'parnumber') <- attr(startvalues, 'parnumber') <- attr(freepars, 'parnumber') <- 
        parnumber - length(freepars[[length(pars)]])
    if(RETURNSTARTVALUES) return(startvalues)
    if(RETURNFREEPARS) return(freepars)
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
    names(parnum) <- names(par) <- c(FNMEANS,FNCOV[tri])
    ret <- new('GroupPars', par=par, est=c(estgmeans,estgcov[tri]), nfact=nfact, 
               parnum=parnum)    
    if(!is.null(startvalues)){
        if(startvalues == 'index')
            return(ret@par)
        else ret@par <- startvalues[[length(startvalues)]]
    }
    if(!is.null(freepars)){
        if(freepars == 'index')
            return(ret@est)
        else ret@est <- freepars[[length(freepars)]]
    }
    return(ret)    
}
