#MHRM optimimization algorithm for confmirt

MHRM <- function(mod, NCYCLES, BURNIN, SEMCYCLES, KDRAWS, TOL, gain, nfactNames, itemloc, fulldata, 
                 nfact, N, K, J, verbose)
{       
    pars <- mod$val$pars    
    VAL <- mod$val
    EST <- mod$est
    IND <- mod$ind 
    sind <- IND$sind
    indlist <- IND$indlist
    prodlist <- IND$prodlist
    constvalues <- mod$val$constvalues
    npars <- mod$npars
    itemtype <- mod$itemtype    
    parind <- IND$parind    
    
    #Burn in 
    pars[constvalues[,1] == 1] <- constvalues[constvalues[,1] == 1,2]
    theta0 <- matrix(0, N, nfact)	    
    cand.t.var <- 1
    tmp <- .1
    for(i in 1:30){			
        theta0 <- draw.thetas(theta0=theta0, lambdas=VAL$lambdas, zetas=VAL$zetas, 
                              guess=VAL$guess, upper=VAL$upper, fulldata=fulldata, K=K, 
                              itemloc=itemloc, cand.t.var=cand.t.var, prior.t.var=VAL$gcov, 
                              prior.mu=VAL$gmeans, estComp=EST$estComp, 
                              prodlist=IND$prodlist, itemtype=itemtype)
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
    SEM.stores <- matrix(0, SEMCYCLES, npars)
    SEM.stores2 <- list()
    phi <- rep(0,sum(sind))	
    h <- matrix(0,npars,npars)		
    Tau <- info <- matrix(0,sum(sind),sum(sind))		
    m.list <- list()	  
    conv <- 0
    k <- 1	
    gamma <- .25
    startvalues <- pars	
    stagecycle <- 1	
    converge <- 1
    nconstvalues <- mod$nconstvalues
    noninvcount <- 0    
       
    ####Big MHRM loop 
    for(cycles in 1:(NCYCLES + BURNIN + SEMCYCLES))								
    { 
        if(cycles == BURNIN + 1) stagecycle <- 2			
        if(stagecycle == 3)
            gamma <- (gain[1] / (cycles - SEMCYCLES - BURNIN - 1))^(gain[2]) - gain[3]
        if(cycles == (BURNIN + SEMCYCLES + 1)){ 
            stagecycle <- 3		
            pars <- rep(0, npars)
            for(i in 1:SEMCYCLES){
                pars <- pars + SEM.stores[i,]
                Tau <- Tau + SEM.stores2[[i]]
            }	
            pars <- pars/SEMCYCLES	
            Tau <- Tau/SEMCYCLES	
            k <- KDRAWS	
            gamma <- .25
        }	
        
        normpars <- sortParsConfmirt(pars=pars, indlist=IND, nfact=nfact, nfactNames=nfactNames)		
        lambdas <- normpars$lambdas
        zetas <- normpars$zetas
        guess <- normpars$guess	
        upper <- normpars$upper
        grouplist$u <- mu <- normpars$mu					
        grouplist$sig <- sig <- normpars$sig			
        
        #Step 1. Generate m_k datasets of theta 
        for(j in 1:4) 
            theta0 <- draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, 
                                  upper=upper, fulldata=fulldata, K=K, itemloc=itemloc, 
                                  cand.t.var=cand.t.var, prior.t.var=sig, prior.mu=mu, 
                                  estComp=EST$estComp, prodlist=IND$prodlist, itemtype=itemtype)
        for(i in 1:k) 
            m.thetas[[i]] <- draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, 
                                         upper=upper, fulldata=fulldata, K=K, itemloc=itemloc, 
                                         cand.t.var=cand.t.var, prior.t.var=sig, prior.mu=mu, 
                                         estComp=EST$estComp, prodlist=IND$prodlist, 
                                         itemtype=itemtype)
        theta0 <- m.thetas[[1]]
        
        #Step 2. Find average of simulated data gradients and hessian 		
        g.m <- h.m <- group.m <- list()
        g <- rep(0, npars)
        h <- matrix(0, npars, npars)	
        for (j in 1:k) {
            g <- rep(NA, npars)            
            thetatemp <- m.thetas[[j]]
            if(!is.null(prodlist)) thetatemp <- prodterms(thetatemp,prodlist)	
            for (i in 1:J){			                
                if(itemtype[i] == 'N2PL' || itemtype[i] == 'N3PL') {
                    temp <- dpars.comp(lambdas[i,][EST$estlam[i,]], zetas[[i]], 
                                       guess[i], fulldata[, itemloc[i]], thetatemp, EST$estGuess[i])
                    ind <- parind[is.na(g)][1]
                    ind2 <- ind + length(temp$grad) - 1
                    g[ind:ind2] <- temp$grad
                    h[ind:ind2, ind:ind2] <- temp$hess						
                    if(!EST$estGuess[i]) g[ind2 + 1] <- 0 #zero for guess
                    if(!EST$estUpper[i]) g[ind2 + 2] <- 0 #zero for upper
                    next
                }                
                if(itemtype[i] == '2PL' || itemtype[i] == '3PL'){
                    temp <- dpars.dich(lambda=lambdas[i, ], zeta=zetas[[i]], g=guess[i], u=upper[i],
                                       dat=fulldata[ ,itemloc[i]], Thetas=thetatemp, 
                                       estGuess=EST$estGuess[i])
                    ind <- parind[is.na(g)][1]					
                    ind2 <- ind + length(temp$g) - 1		                    
                    g[ind:ind2] <- temp$grad
                    h[ind:ind2,ind:ind2] <- temp$hess		
                    if(!EST$estGuess[i]) g[ind2 + 1] <- 0
                    if(!EST$estUpper[i]) g[ind2 + 2] <- 0
                    next
                }
                if(itemtype[i] == 'ordinal'){
                    temp <- dpars.poly(lambdas[i, ],zetas[[i]],
                                       fulldata[ ,itemloc[i]:(itemloc[i+1]-1)],thetatemp)
                    ind <- parind[is.na(g)][1]	
                    ind2 <- ind + length(temp$g) - 1		
                    g[ind:ind2] <- temp$grad
                    h[ind:ind2,ind:ind2] <- temp$hess
                    g[ind2 + 1] <- g[ind2 + 2] <- 0	#zeros for guess + upper
                    next
                }
                if(itemtype[i] == '3PLu'){
                    
                    next
                }
                if(itemtype[i] == '4PL'){
                
                    next
                }            
            }            
            tmp <- d.group(grouplist,as.matrix(thetatemp[ ,1:nfact]))
            g[IND$groupind] <- tmp$g
            h[IND$groupind,IND$groupind] <- tmp$h
            g.m[[j]] <- g
            h.m[[j]] <- h			
        }		
        ave.g <- rep(0,length(g))
        ave.h <- matrix(0,length(g),length(g))		
        for(i in 1:k){
            ave.g <- ave.g + g.m[[i]]
            ave.h <- ave.h + h.m[[i]]
        } 		
        grad <- ave.g/k
        ave.h <- (-1)*ave.h/k
        if(length(IND$parpriors) > 0){
            for(i in 1:length(IND$parpriors)){
                tmp <- IND$parpriors[[i]]
                if(tmp[1] == 1){
                    grad[tmp[2]] <- grad[tmp[2]] - (pars[tmp[2]] - tmp[3])/ tmp[4]^2
                    ave.h[tmp[2],tmp[2]] <- ave.h[tmp[2],tmp[2]] +  1/tmp[4]^2
                }				
                else if(tmp[1] == 2){		
                    tmp2 <- betaprior(tmp[3],tmp[4],pars[tmp[2]])					
                    grad[tmp[2]] <- grad[tmp[2]] + tmp2$g
                    ave.h[tmp[2],tmp[2]] <- ave.h[tmp[2],tmp[2]] + tmp2$h
                }				
            }
        }		
        grad <- grad[parind[sind]]		
        ave.h <- ave.h[parind[sind],parind[sind]] 
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
            correction <-  inv.ave.h %*% grad	
            correction[correction > 1] <- 1
            correction[correction < -1] <- -1			
            parsold <- pars
            correct <- rep(0,npars)
            correct[sind] <- as.vector(correction)
            correct[constvalues[,1] == 1] <- 0
            if(length(IND$equalconstr) > 0)	
                for(i in 1:length(IND$equalconstr))
                    correct[IND$equalconstr[[i]]] <- mean(correct[IND$equalconstr[[i]]])			
            correct[correct[IND$guessind] > .05] <- .05		
            correct[correct[IND$guessind] < -.05] <- -.05
            correct[correct[IND$upperind] > .05] <- .05    	
            correct[correct[IND$upperind] < -.05] <- -.05
            pars <- pars + gamma*correct
            if(verbose && (cycles + 1) %% 10 == 0){ 
                cat(", Max Change =", sprintf("%.4f", max(abs(gamma*correction))), "\n")
                flush.console()
            }			
            pars[IND$covind][pars[IND$covind] > .95] <- parsold[IND$covind][pars[IND$covind] > .95]
            pars[IND$covind][pars[IND$covind] < -.95] <- parsold[IND$covind][pars[IND$covind] < -.95]
            pars[IND$guessind][pars[IND$guessind] < 0] <- parsold[IND$guessind][pars[IND$guessind] < 0]
            pars[IND$upperind][pars[IND$upperind] > 1] <- parsold[IND$upperind][pars[IND$upperind] > 1]
            if(stagecycle == 2){
                SEM.stores[cycles - BURNIN,] <- pars
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
        correction <-  inv.Tau %*% grad
        parsold <- pars
        correct <- rep(0,npars)
        correct[sind] <- as.vector(correction)
        correct[constvalues[,1] == 1] <- 0
        if(length(IND$equalconstr) > 0)		
            for(i in 1:length(IND$equalconstr))
                correct[IND$equalconstr[[i]]] <- mean(correct[IND$equalconstr[[i]]])	
        if(verbose && (cycles + 1) %% 10 == 0){ 
            cat(", gam = ",sprintf("%.3f",gamma),", Max Change = ", 
                sprintf("%.4f",max(abs(gamma*correction))), "\n", sep = '')
            flush.console()		
        }	
        if(all(gamma*correct < TOL)) conv <- conv + 1
        else conv <- 0		
        if(conv == 3) break
        correct[correct[IND$guessind] > .025] <- .025		
        correct[correct[IND$guessind] < -.025] <- -.025	
        correct[correct[IND$upperind] > .025] <- .025    	
        correct[correct[IND$upperind] < -.025] <- -.025
        pars <- pars + gamma*correct	
        pars[IND$covind][pars[IND$covind] > .95] <- parsold[IND$covind][pars[IND$covind] > .95]
        pars[IND$covind][pars[IND$covind] < -.95] <- parsold[IND$covind][pars[IND$covind] < -.95]
        pars[IND$guessind][pars[IND$guessind] < 0] <- parsold[IND$guessind][pars[IND$guessind] < 0]
        pars[IND$upperind][pars[IND$upperind] > 1] <- parsold[IND$upperind][pars[IND$upperind] > 1]
        
        #Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE
        if(gamma == .25) gamma <- 1	
        phi <- phi + gamma*(grad - phi)
        info <- info + gamma*(Tau - phi %*% t(phi) - info)		
    } ###END BIG LOOP
    
    normpars <- sortParsConfmirt(pars=pars, indlist=IND, nfact=nfact, nfactNames=nfactNames)
    ret <- list(pars=pars, info=info, normpars=normpars, theta0=theta0, 
                cycles=cycles - SEMCYCLES - BURNIN, converge=converge)
    ret    
}
