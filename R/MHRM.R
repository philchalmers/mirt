#MHRM optimimization algorithm for confmirt

setGeneric('MHRM', function(pars, list) standardGeneric("MHRM"))

setMethod(
    f = "MHRM",
    signature = signature(pars = 'list', list = 'list'),
    definition = function(pars, list)
    {           
        #for(i in 1:length(pars)) class(pars[[i]]) <- classes[[i]]
        verbose <- list$verbose
        constrain <- list$constrain
        npars <- list$npars
        nfact <- list$nfact
        NCYCLES <- list$NCYCLES
        BURNIN <- list$BURNIN
        SEMCYCLES <- list$SEMCYCLES
        KDRAWS <- list$KDRAWS
        TOL <- list$TOL
        gain <- list$gain
        nfactNames <- list$nfactNames
        itemloc <- list$itemloc
        fulldata <- list$fulldata
        N <- nrow(fulldata)
        #Burn in          
        nfullpars <- 0
        estpars <- c()
        for(i in 1:length(pars)){        
            nfullpars <- nfullpars + length(pars[[i]]@par)                
            estpars <- c(estpars, pars[[i]]@est)
        }
        index <- 1:nfullpars
        structgrouppars <- ExtractGroupPars(pars[[length(pars)]])
        prodlist <- attr(pars, 'prodlist')    
        theta0 <- matrix(0, N, nfact)	    
        cand.t.var <- 1
        tmp <- .1            
        for(i in 1:30){			
            theta0 <- draw.thetas(theta0=theta0, pars=pars, fulldata=fulldata, itemloc=itemloc, 
                                  cand.t.var=cand.t.var, prior.t.var=structgrouppars$gcov, 
                                  prior.mu=structgrouppars$gmeans, prodlist=prodlist)
            if(i > 5){		
                if(attr(theta0,"Proportion Accepted") > .35) cand.t.var <- cand.t.var + 2*tmp 
                else if(attr(theta0,"Proportion Accepted") > .25 && nfact > 3) 
                    cand.t.var <- cand.t.var + tmp
                else if(attr(theta0,"Proportion Accepted") < .2 && nfact < 4) 
                    cand.t.var <- cand.t.var - tmp
                else if(attr(theta0,"Proportion Accepted") < .1) 
                    cand.t.var <- cand.t.var - 2*tmp
                if (cand.t.var < 0){
                    cand.t.var <- tmp		
                    tmp <- tmp / 2
                }		
            }
        }	
        m.thetas <- grouplist <- list()		
        SEM.stores <- SEM.stores2 <- list()     
        m.list <- list()	  
        conv <- 0
        k <- 1	
        gamma <- .25
        longpars <- rep(NA,nfullpars)
        ind1 <- 1
        for(i in 1:length(pars)){
            ind2 <- ind1 + length(pars[[i]]@par) - 1
            longpars[ind1:ind2] <- pars[[i]]@par
            ind1 <- ind2 + 1
        }                  
        stagecycle <- 1	
        converge <- 1    
        noninvcount <- 0    
        L <- c()    
        for(i in 1:(length(pars)))
            L <- c(L, pars[[i]]@est)    
        estindex <- index[estpars]
        L <- diag(as.numeric(L))
        redun_constr <- rep(FALSE, length(estpars)) 
        if(!is.null(constrain)){
            for(i in 1:length(constrain)){            
                L[constrain[[i]], constrain[[i]]] <- 1/length(constrain[[i]]) 
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
        ####Big MHRM loop 
        for(cycles in 1:(NCYCLES + BURNIN + SEMCYCLES))								
        { 
            if(cycles == BURNIN + 1) stagecycle <- 2			
            if(stagecycle == 3)
                gamma <- (gain[1] / (cycles - SEMCYCLES - BURNIN - 1))^(gain[2]) - gain[3]
            if(cycles == (BURNIN + SEMCYCLES + 1)){ 
                stagecycle <- 3		
                longpars <- SEM.stores[[1]]
                Tau <- SEM.stores2[[1]]
                for(i in 2:SEMCYCLES){
                    longpars <- longpars + SEM.stores[[i]]
                    Tau <- Tau + SEM.stores2[[i]]
                }	
                longpars <- longpars/SEMCYCLES	
                Tau <- Tau/SEMCYCLES	
                k <- KDRAWS	
                gamma <- .25
            }
            
            ind1 <- 1
            for(i in 1:length(pars)){
                ind2 <- ind1 + length(pars[[i]]@par) - 1
                pars[[i]]@par <- longpars[ind1:ind2]
                ind1 <- ind2 + 1            
            }        
            structgrouppars <- ExtractGroupPars(pars[[length(pars)]])
            
            #Step 1. Generate m_k datasets of theta 
            for(j in 1:4) 
                theta0 <- draw.thetas(theta0=theta0, pars=pars, fulldata=fulldata, itemloc=itemloc, 
                                      cand.t.var=cand.t.var, prior.t.var=structgrouppars$gcov, 
                                      prior.mu=structgrouppars$gmeans, prodlist=prodlist)
            for(i in 1:k) 
                m.thetas[[i]] <- draw.thetas(theta0=theta0, pars=pars, fulldata=fulldata, itemloc=itemloc, 
                                             cand.t.var=cand.t.var, prior.t.var=structgrouppars$gcov, 
                                             prior.mu=structgrouppars$gmeans, prodlist=prodlist)
            theta0 <- m.thetas[[1]]
            
            #Step 2. Find average of simulated data gradients and hessian 		
            g.m <- h.m <- group.m <- list()
            longpars <- g <- rep(0, nfullpars)
            h <- matrix(0, nfullpars, nfullpars)        
            for (j in 1:k) { 
                ind1 <- 1
                thetatemp <- m.thetas[[j]]
                if(!is.null(prodlist)) thetatemp <- prodterms(thetatemp,prodlist)	
                for (i in 1:length(pars)){	
                    deriv <- Deriv(x=pars[[i]], Theta=thetatemp)
                    ind2 <- ind1 + length(deriv$grad) - 1
                    longpars[ind1:ind2] <- pars[[i]]@par
                    g[ind1:ind2] <- pars[[i]]@gradient <- deriv$grad
                    h[ind1:ind2, ind1:ind2] <- pars[[i]]@hessian <- deriv$hess
                    ind1 <- ind2 + 1 
                }            
                g.m[[j]] <- g %*% L 
                h.m[[j]] <- L %*% h %*% L 			
            }		
            ave.g <- rep(0,length(g))
            ave.h <- matrix(0,length(g),length(g))		
            for(i in 1:k){
                ave.g <- ave.g + g.m[[i]]
                ave.h <- ave.h + h.m[[i]]
            } 		
            grad <- ave.g/k
            ave.h <- (-1)*ave.h/k
    #         if(length(IND$parpriors) > 0){
    #             for(i in 1:length(IND$parpriors)){
    #                 tmp <- IND$parpriors[[i]]
    #                 if(tmp[1] == 1){
    #                     grad[tmp[2]] <- grad[tmp[2]] - (pars[tmp[2]] - tmp[3])/ tmp[4]^2
    #                     ave.h[tmp[2],tmp[2]] <- ave.h[tmp[2],tmp[2]] +  1/tmp[4]^2
    #                 }				
    #                 else if(tmp[1] == 2){		
    #                     tmp2 <- betaprior(tmp[3],tmp[4],pars[tmp[2]])					
    #                     grad[tmp[2]] <- grad[tmp[2]] + tmp2$g
    #                     ave.h[tmp[2],tmp[2]] <- ave.h[tmp[2],tmp[2]] + tmp2$h
    #                 }				
    #             }
    #         }	                
            grad <- grad[1, estpars & !redun_constr]		
            ave.h <- ave.h[estpars & !redun_constr, estpars & !redun_constr] 
            if(is.na(attr(theta0,"log.lik"))) stop('Estimation halted. Model did not converge.')		
            if(verbose){
                if((cycles + 1) %% 10 == 0){
                    if(cycles < BURNIN)
                        cat("Stage 1: Cycle = ", cycles + 1, ", Log-Lik = ", 
                            sprintf("%.1f",attr(theta0,"log.lik")), sep="")
                    if(cycles > BURNIN && cycles < BURNIN + SEMCYCLES)
                        cat("Stage 2: Cycle = ", cycles-BURNIN+1, ", Log-Lik = ",
                            sprintf("%.1f",attr(theta0,"log.lik")), sep="")
                    if(cycles > BURNIN + SEMCYCLES)
                        cat("Stage 3: Cycle = ", cycles-BURNIN-SEMCYCLES+1, 
                            ", Log-Lik = ", sprintf("%.1f",attr(theta0,"log.lik")), sep="")					
                }
            }			
            if(stagecycle < 3){	
                ave.h <- as(ave.h,'sparseMatrix')
                inv.ave.h <- try(solve(ave.h))			
                if(class(inv.ave.h) == 'try-error'){
                    inv.ave.h <- try(qr.solve(ave.h + 2*diag(ncol(ave.h))))
                    noninvcount <- noninvcount + 1
                    if(noninvcount == 3) 
                        stop('\nEstimation halted during burn in stages, solution is unstable')
                }
                correction <-  as.numeric(inv.ave.h %*% grad)
                correction[correction > .5] <- 1
                correction[correction < -.5] <- -1
                longpars[estindex_unique] <- longpars[estindex_unique] + gamma*correction           
                if(!is.null(constrain))
                    for(i in 1:length(constrain))
                        longpars[index %in% constrain[[i]][-1]] <- longpars[constrain[[i]][1]]           
                if(verbose && (cycles + 1) %% 10 == 0){ 
                    cat(", Max Change =", sprintf("%.4f", max(abs(gamma*correction))), "\n")
                    flush.console()
                }			            
                if(stagecycle == 2){
                    SEM.stores[[cycles - BURNIN]] <- longpars
                    SEM.stores2[[cycles - BURNIN]] <- ave.h
                }	
                next
            }	 
            
            #Step 3. Update R-M step		
            Tau <- Tau + gamma*(ave.h - Tau)
            Tau <- as(Tau,'sparseMatrix')	
            inv.Tau <- solve(Tau)
            if(class(inv.Tau) == 'try-error'){
                inv.Tau <- try(qr.solve(Tau + 2 * diag(ncol(Tau))))
                noninvcount <- noninvcount + 1
                if(noninvcount == 3) 
                    stop('\nEstimation halted during stage 3, solution is unstable')
            }		
            correction <-  as.numeric(inv.Tau %*% grad)
            longpars[estindex_unique] <- longpars[estindex_unique] + gamma*correction           
            if(!is.null(constrain))
                for(i in 1:length(constrain))
                    longpars[index %in% constrain[[i]][-1]] <- longpars[constrain[[i]][1]]
            if(verbose && (cycles + 1) %% 10 == 0){ 
                cat(", gam = ",sprintf("%.3f",gamma),", Max Change = ", 
                    sprintf("%.4f",max(abs(gamma*correction))), "\n", sep = '')
                flush.console()		
            }	
            if(all(abs(gamma*correction) < TOL)) conv <- conv + 1
                else conv <- 0		
            if(conv == 3) break        
            
            #Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE
            if(gamma == .25){
                gamma <- 1	
                phi <- rep(0, length(grad))            
                info <- matrix(0, length(grad), length(grad))
            }
            phi <- phi + gamma*(grad - phi)
            info <- info + gamma*(Tau - phi %*% t(phi) - info)		
        } ###END BIG LOOP      
        
        SEtmp <- diag(solve(info))    	
        if(any(SEtmp < 0)){
            warning("Information matrix is not positive definite, negative SEs set to 'NA'.\n")
            SEtmp <- rep(NA, length(SEtmp))
        } else SEtmp <- sqrt(SEtmp)
        SE <- rep(NA, length(longpars))
        SE[estindex_unique] <- SEtmp
        if(!is.null(constrain))
            for(i in 1:length(constrain))
                SE[index %in% constrain[[i]][-1]] <- SE[constrain[[i]][1]]
        ind1 <- 1
        for(i in 1:length(pars)){
            ind2 <- ind1 + length(pars[[i]]@par) - 1
            pars[[i]]@SEpar <- SE[ind1:ind2]
            ind1 <- ind2 + 1            
        }    
        ret <- list(pars=pars, cycles = cycles - BURNIN - SEMCYCLES, info=info, converge=converge)
        ret        
    }
)
