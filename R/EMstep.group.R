EM.group <- function(pars, constrain, PrepList, list, Theta, debug)
{    
    if(debug == 'EM') browser()
    verbose <- list$verbose        
    nfact <- list$nfact
    NCYCLES <- list$NCYCLES    
    MSTEPMAXIT <- list$MSTEPMAXIT
    TOL <- list$TOL    
    BFACTOR <- list$BFACTOR
    itemloc <- list$itemloc
    ngroups <- length(pars)
    specific <- list$specific
    sitems <- list$sitems
    theta <- list$theta
    J <- length(itemloc) - 1
    nfullpars <- 0
    estpars <- c()
    gfulldata <- gtheta0 <- gstructgrouppars <- vector('list', ngroups)
    for(g in 1:ngroups){
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1]])
        gfulldata[[g]] <- PrepList[[g]]$fulldata
        gtheta0[[g]] <- matrix(0, nrow(gfulldata[[g]]), nfact)
        for(i in 1:(J+1)){        
            nfullpars <- nfullpars + length(pars[[g]][[i]]@par)                
            estpars <- c(estpars, pars[[g]][[i]]@est)
        }
    }
    listpars <- vector('list', ngroups)
    for(g in 1:ngroups){
        listpars[[g]] <- list()
        for(i in 1:(J + 1)){            
            listpars[[g]][[i]] <- pars[[g]][[i]]@par
        }
    }
    lastpars2 <- lastpars1 <- listpars
    index <- 1:nfullpars    
    longpars <- rep(NA,nfullpars)
    ind1 <- 1
    for(g in 1:ngroups){
        for(i in 1:(J+1)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1
            longpars[ind1:ind2] <- pars[[g]][[i]]@par
            ind1 <- ind2 + 1
        }                  
    }
    stagecycle <- 1    
    converge <- 1
    LLwarn <- FALSE
    inverse_fail_count <- 1
    ##
    L <- c()    
    for(g in 1:ngroups)
        for(i in 1:(J+1))
            L <- c(L, pars[[g]][[i]]@est)    
    estindex <- index[estpars]
    L <- diag(as.numeric(L))
    redun_constr <- rep(FALSE, length(estpars)) 
    if(length(constrain) > 0){
        for(i in 1:length(constrain)){            
            L[constrain[[i]], constrain[[i]]] <- 1
            for(j in 2:length(constrain[[i]]))
                redun_constr[constrain[[i]][j]] <- TRUE
        }
    }
    estindex_unique <- index[estpars & !redun_constr]
    if(any(diag(L)[!estpars] > 0)){
        redindex <- index[!estpars]        
        stop('Constraint applied to fixed parameter(s) ', 
             paste(redindex[diag(L)[!estpars] > 0]), ' but should only be applied to 
                 estimated parameters. Please fix!')
    }   
    Prior <- prior <- gstructgrouppars <- rlist <- r <- list()        
    #make sure constrained pars are equal    
    tmp <- rowSums(L)
    tmp[tmp == 0] <- 1
    tmp <- matrix(1/tmp, length(longpars), length(longpars), byrow = TRUE)
    tmp2 <- abs(diag(L) - 1)
    longpars <- diag((tmp * L) * longpars) + tmp2 * longpars
    LL <- 0
    for(g in 1:ngroups)
        r[[g]] <- PrepList[[g]]$tabdata[, ncol(PrepList[[g]]$tabdata)]        
    LBOUND <- UBOUND <- c()
    for(g in 1:ngroups){
        for(i in 1:(J+1)){
            LBOUND <- c(LBOUND, pars[[g]][[i]]@lbound)    
            UBOUND <- c(UBOUND, pars[[g]][[i]]@ubound)    
        }
    }
    
    #EM     
    for (cycles in 1:NCYCLES){  
        #priors
        for(g in 1:ngroups){
            gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1]])        
            if(BFACTOR){
                prior[[g]] <- dnorm(theta, 0, 1)
                prior[[g]] <- prior[[g]]/sum(prior[[g]])
                Prior[[g]] <- mvtnorm::dmvnorm(Theta[,1:2])
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])                
                next  
            } 
            Prior[[g]] <- mvtnorm::dmvnorm(Theta,gstructgrouppars[[g]]$gmeans,
                                           gstructgrouppars[[g]]$gcov)
            Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
        }
        #Estep
        lastLL <- LL
        LL <- 0
        for(g in 1:ngroups){
            if(BFACTOR){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=PrepList[[g]]$tabdata, 
                                            Theta=Theta, prior=prior[[g]],
                                            specific=specific, sitems=sitems, 
                                            itemloc=itemloc, debug=debug)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=PrepList[[g]]$tabdata, 
                                         Theta=Theta, prior=Prior[[g]], itemloc=itemloc, 
                                         debug=debug)                      
            }
            LL <- LL + sum(r[[g]]*log(rlist[[g]]$expected))
        }
        if(LL < lastLL && cycles > 1) LLwarn <- TRUE
        for(g in 1:ngroups){
            for(i in 1:J){
                tmp <- c(itemloc[i]:(itemloc[i+1] - 1))
                pars[[g]][[i]]@rs <- rlist[[g]]$r1[, tmp]           
            }
        }
        if(verbose){
            print(LL)                            
            flush.console()
        }
        
        #Mstep        
        lastpars2 <- lastpars1
        lastpars1 <- listpars
        preMstep.longpars <- longpars
        lastgrad <- 0
        stepLimit <- .1        
        for(mstep in 1:MSTEPMAXIT){                        
            #Reload pars list
            ind1 <- 1
            for(g in 1:ngroups){                
                for(i in 1:(J+1)){
                    ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1
                    pars[[g]][[i]]@par <- longpars[ind1:ind2]
                    ind1 <- ind2 + 1
                }
                for(i in 1:(J+1)){                        
                    #apply sum(t) == 1 constraint for mcm
                    if(is(pars[[g]][[i]], 'mcm')){                
                        tmp <- pars[[g]][[i]]@par
                        cat <- pars[[g]][[i]]@ncat
                        tmp2 <- (length(tmp) - (cat-1)):length(tmp) 
                        Num <- exp(tmp[tmp2])
                        tmp <- Num/sum(Num)                                    
                        pars[[g]][[i]]@par[tmp2] <- tmp
                    }
                }
            }
            #reset longpars and gradient
            g <- rep(0, nfullpars)
            h <- matrix(0, nfullpars, nfullpars)
            ind1 <- 1                    
            for(group in 1:ngroups){
                for (i in 1:J){                    
                    deriv <- Deriv(x=pars[[group]][[i]], Theta=Theta, EM = TRUE, prior=Prior[[group]])
                    ind2 <- ind1 + length(deriv$grad) - 1
                    longpars[ind1:ind2] <- pars[[group]][[i]]@par
                    g[ind1:ind2] <- pars[[group]][[i]]@gradient <- deriv$grad
                    h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
                    ind1 <- ind2 + 1 
                }
                i <- i + 1                        
                deriv <- Deriv(x=pars[[group]][[i]], Theta=Theta, EM = TRUE, 
                               pars=pars[[group]], tabdata=PrepList[[group]]$tabdata,
                               itemloc=itemloc)
                ind2 <- ind1 + length(deriv$grad) - 1
                longpars[ind1:ind2] <- pars[[group]][[i]]@par
                g[ind1:ind2] <- pars[[group]][[i]]@gradient <- deriv$grad
                h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
                ind1 <- ind2 + 1
            }
            grad <- g %*% L 
            hess <- L %*% h %*% L 			                   
            grad <- grad[1, estpars & !redun_constr]            
            if(any(is.na(grad))) 
                stop('Model did not converge (unacceptable gradient caused by extreme parameter values)')            
            Hess <- Matrix(hess[estpars & !redun_constr, estpars & !redun_constr], sparse = TRUE)            
            inv.Hess <- try(solve(Hess), silent = TRUE)        	                        
            if(class(inv.Hess) == 'try-error'){             
                if(inverse_fail_count == 5) 
                    stop('Hessian is not invertable. Likelihood surface is likely too flat.') 
                inverse_fail_count <- inverse_fail_count + 1
                inv.Hess <- Hess
                tmp <- .1*diag(inv.Hess)
                tmp[tmp > -25] <- -.25
                diag(inv.Hess) <- diag(inv.Hess) + tmp
                inv.Hess <- try(solve(inv.Hess))                
            }
            correction <- as.vector(inv.Hess %*% grad)
            #keep steps smaller
            correction[correction > stepLimit] <- stepLimit
            correction[correction < -stepLimit] <- -stepLimit
            #prevent guessing/upper pars from moving more than .001 at all times
            names(correction) <- names(estpars[estpars & !redun_constr])
            if(stepLimit > .002){                
                tmp <- correction[names(correction) == 'g']
                tmp[abs(tmp) > .002] <- sign(tmp[abs(tmp) > .002]) * .002
                correction[names(correction) == 'g'] <- tmp
                tmp <- correction[names(correction) == 'u']
                tmp[abs(tmp) > .002] <- sign(tmp[abs(tmp) > .002]) * .002
                correction[names(correction) == 'u'] <- tmp
            }
            longpars[estindex_unique] <- longpars[estindex_unique] - correction             
            longpars[longpars < LBOUND] <- LBOUND[longpars < LBOUND]
            longpars[longpars > UBOUND] <- UBOUND[longpars > UBOUND]
            if(mstep > 1){
                if (any(grad*lastgrad < 0.0)){    				# any changed sign
                    newcorrection <- rep(0, length(correction))
                    newcorrection[grad*lastgrad < 0.0] <- .5*correction[grad*lastgrad < 0.0]
                    longpars[estindex_unique] <- longpars[estindex_unique] + newcorrection	# back up 1/2
                    stepLimit <- 0.5*stepLimit				# split the difference                    
                    lastgrad <- 0
                } else {
                    lastgrad <- grad
                }        
            }            
            if(length(constrain) > 0)
                for(i in 1:length(constrain))
                    longpars[index %in% constrain[[i]][-1]] <- longpars[constrain[[i]][1]]
            if(all(abs(correction) < .0001)) break            
            if(is.list(debug)) print(longpars[debug[[1]]])
        }#END MSTEP        
        if(all(abs(preMstep.longpars - longpars) < TOL) || abs(lastLL - LL) < .01 ) break 
        for(g in 1:ngroups)
            for(i in 1:J) 
                listpars[[g]][[i]] <- pars[[g]][[i]]@par         
    } #END EM          
    
    if(cycles == NCYCLES) converge <- 0
    ret <- list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                logLik=LL, rlist=rlist)
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
    X <- tabdata[ ,1:(ncol(tabdata) - 1), drop = FALSE]    
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
