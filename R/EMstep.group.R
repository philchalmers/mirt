EM.group <- function(pars, constrain, Ls, Data, PrepList, list, Theta, DERIV, solnp_args, control)
{
    verbose <- list$verbose
    lrPars <- list$lrPars
    nfact <- list$nfact
    if(list$EH && nfact != 1L)
        stop('empirical histogram only available for unidimensional models', call.=FALSE)
    NCYCLES <- list$NCYCLES
    TOL <- list$TOL
    BFACTOR <- list$BFACTOR
    CUSTOM.IND <- list$CUSTOM.IND
    itemloc <- list$itemloc
    ngroups <- length(pars)
    specific <- list$specific
    sitems <- list$sitems
    theta <- list$theta
    full <- list$full
    J <- length(itemloc) - 1L
    nfullpars <- 0L
    estpars <- c()
    prodlist <- PrepList[[1L]]$prodlist
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            nfullpars <- nfullpars + length(pars[[g]][[i]]@par)
            estpars <- c(estpars, pars[[g]][[i]]@est)
        }
        if(length(lrPars)){
            nfullpars <- nfullpars + length(lrPars@par)
            estpars <- c(estpars, lrPars@est)
        }
    }
    listpars <- vector('list', ngroups)
    for(g in 1L:ngroups){
        listpars[[g]] <- list()
        for(i in 1L:(J + 1L)){
            listpars[[g]][[i]] <- pars[[g]][[i]]@par
        }
        if(length(lrPars))
            listpars[[g]][[i+1L]] <- lrPars@par
    }
    index <- 1L:nfullpars
    longpars <- rep(NA,nfullpars)
    latent_longpars <- logical(nfullpars)
    ind1 <- 1L
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1L
            longpars[ind1:ind2] <- pars[[g]][[i]]@par
            if(i == (J+1L)) latent_longpars[ind1:ind2] <- TRUE
            ind1 <- ind2 + 1L
        }
        if(length(lrPars)){
            ind2 <- ind1 + length(lrPars@par) - 1L
            longpars[ind1:ind2] <- lrPars@par
            ind1 <- ind2 + 1L
        }
    }
    converge <- TRUE
    L <- Ls$L
    redun_constr <- Ls$redun_constr
    estindex_unique <- index[estpars & !redun_constr]
    if(any(diag(L)[!estpars] > 0L)){
        redindex <- index[!estpars]
        stop('Constraint applied to fixed parameter(s) ',
             paste(paste0(redindex[diag(L)[!estpars] > 0L]), ''), ' but should only be applied to
                 estimated parameters. Please fix!', call.=FALSE)
    }
    prior <- rlist <- r <- vector('list', ngroups)
    #make sure constrained pars are equal
    tmp <- rowSums(L)
    tmp[tmp == 0L] <- 1L
    check <- as.numeric(L %*% longpars) / tmp
    longpars[estpars] <- check[estpars]
    LL <- 0
    LBOUND <- UBOUND <- c()
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            LBOUND <- c(LBOUND, pars[[g]][[i]]@lbound)
            UBOUND <- c(UBOUND, pars[[g]][[i]]@ubound)
        }
        if(length(lrPars)){
            LBOUND <- c(LBOUND, lrPars@lbound)
            UBOUND <- c(UBOUND, lrPars@ubound)
        }
    }
    est <- c()
    for(g in 1L:ngroups){
       for(j in 1L:(J+1L))
           est <- c(est, pars[[g]][[j]]@est)
       if(length(lrPars))
           est <- c(est, rep(FALSE, length(lrPars@est)))
    }
    if(length(constrain))
       for(i in 1L:length(constrain))
           est[constrain[[i]][-1L]] <- FALSE
    names(longpars) <- names(est)
    if(list$Moptim != 'BFGS') {
        Moptim <- list$Moptim
    } else {
        Moptim <- if(all(c(LBOUND[est], UBOUND[est]) %in% c(-Inf, Inf))) 'BFGS' else 'L-BFGS-B'
    }
    if(Moptim == 'NR' && sum(est) > 300L && list$message)
        message('NR optimizer should not be used for models with a large number of parameters.
                Use the optimizer = \'BFGS\' or \'nlminb\' instead.')
    EMhistory <- matrix(NA, NCYCLES+1L, length(longpars))
    EMhistory[1L,] <- longpars
    ANY.PRIOR <- rep(FALSE, ngroups)
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta, prodlist)
    if(BFACTOR){
        Thetabetween <- thetaComb(theta=theta, nfact=nfact-ncol(sitems))
        for(g in 1L:ngroups){
            prior[[g]] <- dnorm(theta, 0, 1)
            prior[[g]] <- prior[[g]]/sum(prior[[g]])
        }
    }
    gTheta <- vector('list', ngroups)
    for(g in 1L:ngroups){
        ANY.PRIOR[g] <- any(sapply(pars[[g]], function(x) x@any.prior))
        gTheta[[g]] <- Theta
    }
    preMstep.longpars2 <- preMstep.longpars <- longpars
    accel <- 0; Mrate <- ifelse(list$SEM, 1, .4)
    Estep.time <- Mstep.time <- 0
    collectLL <- rep(NA, NCYCLES)
    hess <- matrix(0)
    if(list$BL){
        start <- proc.time()[3L]
        lower <- LBOUND[est]; upper <- UBOUND[est]
        Moptim <- ifelse(any(is.finite(lower) | is.finite(upper)), 'L-BFGS-B', 'BFGS')
        if(Moptim == 'BFGS'){
            control <- list(fnscale=-1, reltol=TOL)
            lower <- -Inf; upper <- Inf
        } else {
            control <- list(fnscale=-1, pgtol=TOL)
        }
        opt <- try(optim(longpars[est], BL.LL, BL.grad, est=est, longpars=longpars,
                         pars=pars, ngroups=ngroups, J=J, itemloc=itemloc,
                         Theta=Theta, PrepList=PrepList, BFACTOR=BFACTOR,
                         specific=specific, sitems=sitems, CUSTOM.IND=CUSTOM.IND,
                         EH=list$EH, constrain=constrain, EHPrior=NULL, Data=Data, method=Moptim,
                         control=control, hessian=list$SE,
                         lower=lower, upper=upper), silent=TRUE)
        cycles <- as.integer(opt$counts[1L])
        longpars[est] <- opt$par
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
        converge <- opt$convergence == 0
        if(list$SE) hess <- opt$hessian
        tmp <- updatePrior(pars=pars, Theta=Theta, Thetabetween=Thetabetween,
                           list=list, ngroups=ngroups, nfact=nfact, prior=prior,
                           J=J, BFACTOR=BFACTOR, sitems=sitems, cycles=cycles, rlist=rlist)
        Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
        LL <- 0
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        for(g in 1L:ngroups){
            if(BFACTOR){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                            Theta=Theta, prior=prior[[g]], Prior=Prior[[g]],
                                            Priorbetween=Priorbetween[[g]], specific=specific,
                                            sitems=sitems, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                         CUSTOM.IND=CUSTOM.IND, Theta=Theta,
                                         prior=Prior[[g]], itemloc=itemloc)
            }
            LL <- LL + sum(Data$Freq[[g]]*log(rlist[[g]]$expected))
        }
        Estep.time <- Estep.time + proc.time()[3L] - start
    } else {
        #EM
        for (cycles in 1L:NCYCLES){
            #priors
            start <- proc.time()[3L]
            if(length(lrPars)) lrPars@mus <- lrPars@X %*% lrPars@beta
            tmp <- updatePrior(pars=pars, Theta=Theta, Thetabetween=Thetabetween,
                               list=list, ngroups=ngroups, nfact=nfact, prior=prior,
                               J=J, BFACTOR=BFACTOR, sitems=sitems, cycles=cycles,
                               rlist=rlist, full=full, lrPars=lrPars)
            Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
            Elist <- Estep(pars=pars, Data=Data, Theta=Theta, prior=prior, Prior=Prior,
                           Priorbetween=Priorbetween, specific=specific, sitems=sitems,
                           ngroups=ngroups, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                           BFACTOR=BFACTOR, rlist=rlist, full=full, Etable=list$Etable)
            rlist <- Elist$rlist; LL <- Elist$LL
            collectLL[cycles] <- LL
            if(is.nan(LL))
                stop('Optimization error: Could not compute observed log-likelihood. Try
                     estimating with different starting values by passing GenRandomPars = TRUE',
                     call.=FALSE)
            if(!list$SEM){
                if(cycles > 1L){
                    tmp <- collectLL[cycles-1L] - collectLL[cycles]
                    if(tmp < 0)
                        Mrate <- exp(tmp)
                    Mrate <- ifelse(is.finite(Mrate), Mrate, 1e-6)
                }
            }
            for(g in 1L:ngroups){
                for(i in 1L:J)
                    pars[[g]][[i]]@dat <- rlist[[g]]$r1[, c(itemloc[i]:(itemloc[i+1L] - 1L)),
                                                        drop=FALSE]
                pars[[g]][[J+1L]]@rr <- rowSums(rlist[[g]]$r1)
            }
            Estep.time <- Estep.time + proc.time()[3L] - start
            start <- proc.time()[3L]
            preMstep.longpars2 <- preMstep.longpars
            preMstep.longpars <- longpars
            if(all(!est)) break
            if(is.nan(TOL)) break
            longpars <- Mstep(pars=pars, est=est, longpars=longpars, ngroups=ngroups, J=J,
                              gTheta=gTheta, itemloc=itemloc, Prior=Prior, ANY.PRIOR=ANY.PRIOR,
                              CUSTOM.IND=CUSTOM.IND, SLOW.IND=list$SLOW.IND,
                              PrepList=PrepList, L=L, UBOUND=UBOUND, LBOUND=LBOUND, Moptim=Moptim,
                              BFACTOR=BFACTOR, nfact=nfact, Thetabetween=Thetabetween,
                              rlist=rlist, constrain=constrain, DERIV=DERIV, Mrate=Mrate,
                              TOL=list$MSTEPTOL, solnp_args=solnp_args, full=full, lrPars=lrPars,
                              control=control)
            EMhistory[cycles+1L,] <- longpars
            if(verbose)
                cat(sprintf('\rIteration: %d, Log-Lik: %.3f, Max-Change: %.5f',
                            cycles, LL, max(abs(preMstep.longpars - longpars))))
            if(all(abs(preMstep.longpars - longpars) < TOL)){
                pars <- reloadPars(longpars=longpars, pars=pars,
                                   ngroups=ngroups, J=J)
                if(length(lrPars)){
                    lrPars@par <- longpars[lrPars@parnum]
                    lrPars@beta[] <- matrix(lrPars@par, lrPars@nfixed, lrPars@nfact)
                }
                break
            }
            if(list$accelerate != 'none' && cycles %% 3 == 0L){
                if(list$accelerate == 'Ramsay'){
                    if(Mrate > .01){
                        dX2 <- preMstep.longpars - preMstep.longpars2
                        dX <- longpars - preMstep.longpars
                        d2X2 <- dX - dX2
                        ratio <- sqrt((dX %*% dX) / (d2X2 %*% d2X2))
                        accel <- 1 - ratio
                        if(accel < -5) accel <- -5
                        tmp <- (1 - accel) * longpars + accel * preMstep.longpars
                        longpars[!latent_longpars] <- tmp[!latent_longpars]
                    }
                } else if(list$accelerate == 'squarem'){
                    r <- preMstep.longpars - preMstep.longpars2
                    v <- (longpars - preMstep.longpars) - r
                    ratio <- sqrt((r %*% r) / (v %*% v))
                    accel <- -ratio
                    if(accel > -1){
                        accel <- -1
                    } else {
                        count <- 1L
                        while(TRUE){
                            tmp <- preMstep.longpars2 - 2 * accel * r  + accel^2 * v
                            longpars[!latent_longpars] <- tmp[!latent_longpars]
                            pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
                            Elist <- Estep(pars=pars, Data=Data, Theta=Theta, prior=prior, Prior=Prior,
                                           Priorbetween=Priorbetween, specific=specific, sitems=sitems,
                                           ngroups=ngroups, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                           BFACTOR=BFACTOR, rlist=rlist, full=full)
                            if(Elist$LL <= collectLL[cycles]){
                                accel <- (accel - 1) / 2
                                count <- count + 1L
                            } else break
                            if(count == 5L){
                                accel <- -1
                                break
                            }
                        }
                    }
                    tmp <- preMstep.longpars2 - 2 * accel * r  + accel^2 * v
                    longpars[!latent_longpars] <- tmp[!latent_longpars]
                } else stop('acceleration option not defined', call.=FALSE)
            }
            pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
            if(length(lrPars)){
                lrPars@par <- longpars[lrPars@parnum]
                lrPars@beta[] <- matrix(lrPars@par, lrPars@nfixed, lrPars@nfact)
            }
            Mstep.time <- Mstep.time + proc.time()[3L] - start
        } #END EM
        if(cycles == NCYCLES){
            if(list$message)
                message('EM cycles terminated after ', cycles, ' iterations.')
            converge <- FALSE
        } else if(cycles == 1L && !all(!est)){
            if(list$warn && !is.nan(TOL))
                warning('M-step optimimizer converged immediately. Solution is either at the ML or
                     starting values are causing issues and should be adjusted. ', call.=FALSE)
        }
        if(cycles > 1L && list$warn && !ANY.PRIOR){
            diff <- c(-Inf, na.omit(collectLL)) - c(na.omit(collectLL), Inf)
            if(any(diff[length(diff):ceiling(length(diff)*.9)] > 0))
                warning('Log-likelihood was decreasing near the ML solution. EM method may be unstable',
                        call.=FALSE)
        }
    }
    infological <- estpars & !redun_constr
    correction <- numeric(length(estpars[estpars & !redun_constr]))
    names(correction) <- names(estpars[estpars & !redun_constr])
    collectLL <- as.numeric(na.omit(collectLL))
    LP <- 0
    if(any(ANY.PRIOR)){
        if(length(lrPars)){
            if(lrPars@any.prior)
                LP <- LL.Priors(x=lrPars, LL=LP)
        }
        for(g in 1L:length(pars)){
            for(i in 1L:length(pars[[1L]]))
                if(pars[[g]][[i]]@any.prior)
                    LP <- LL.Priors(x=pars[[g]][[i]], LL=LP)
        }
    }
    LP <- unname(LP)
    if(list$SEM){
        h <- matrix(0, nfullpars, nfullpars)
        ind1 <- 1L
        for(group in 1L:ngroups){
            for (i in 1L:J){
                deriv <- Deriv(x=pars[[group]][[i]], Theta=Theta, estHess=TRUE)
                ind2 <- ind1 + length(deriv$grad) - 1L
                h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
                ind1 <- ind2 + 1L
            }
            i <- i + 1L
            deriv <- Deriv(x=pars[[group]][[i]], CUSTOM.IND=CUSTOM.IND,
                           Theta=Theta, EM = TRUE,
                           pars=pars[[group]], tabdata=Data$tabdatalong,
                           freq=Data$Freq[[group]], prior=Prior[[group]],
                           itemloc=itemloc, estHess=TRUE)
            ind2 <- ind1 + length(deriv$grad) - 1L
            h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
            ind1 <- ind2 + 1L
            if(length(lrPars)){
                gp <- ExtractGroupPars(pars[[group]][[J+1L]])
                tmp <- Mstep.LR(Theta=Theta, CUSTOM.IND=CUSTOM.IND, pars=pars[[group]],
                                itemloc=itemloc, fulldata=Data$fulldata[[1L]], prior=Prior[[group]],
                                lrPars=lrPars, retscores=TRUE)
                deriv <- Deriv(lrPars, cov=gp$gcov, theta=tmp)
                for(i in 0L:(ncol(deriv$grad)-1L))
                    h[lrPars@parnum[1L:nrow(deriv$grad) + nrow(deriv$grad)*i],
                      lrPars@parnum[1L:nrow(deriv$grad) + nrow(deriv$grad)*i]] <- deriv$hess
            }
        }
        hess <- updateHess(h=h, L=L)
        hess <- hess[estpars & !redun_constr, estpars & !redun_constr]
        ret <- list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                    logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological, Moptim=Moptim,
                    estindex_unique=estindex_unique, correction=correction, hess=hess, Prior=Prior,
                    estpars=estpars & !redun_constr, redun_constr=redun_constr, ngroups=ngroups,
                    LBOUND=LBOUND, UBOUND=UBOUND, EMhistory=na.omit(EMhistory), random=list(),
                    time=c(Estep=as.numeric(Estep.time), Mstep=as.numeric(Mstep.time)),
                    collectLL=collectLL, shortpars=longpars[estpars & !redun_constr],
                    lrPars=lrPars, logPrior=LP, fail_invert_info=FALSE)
    } else {
        ret <- list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                    logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological, Moptim=Moptim,
                    estindex_unique=estindex_unique, correction=correction, hess=hess, random=list(),
                    Prior=Prior, time=c(Estep=as.numeric(Estep.time), Mstep=as.numeric(Mstep.time)),
                    prior=prior, Priorbetween=Priorbetween, sitems=sitems, collectLL=collectLL,
                    shortpars=longpars[estpars & !redun_constr], lrPars=lrPars,
                    logPrior=LP, fail_invert_info=FALSE)
    }
    for(g in 1L:ngroups)
        for(i in 1L:J)
            ret$pars[[g]][[i]]@dat <- matrix(0)
    ret
}
