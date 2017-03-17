EM.group <- function(pars, constrain, Ls, Data, PrepList, list, Theta, DERIV, solnp_args, control)
{
    verbose <- list$verbose
    lrPars <- list$lrPars
    nfact <- list$nfact
    if(list$dentype == 'EH' && nfact != 1L)
        stop('empirical histogram only available for unidimensional models', call.=FALSE)
    NCYCLES <- list$NCYCLES
    TOL <- list$TOL
    CUSTOM.IND <- list$CUSTOM.IND
    dentype <- list$dentype
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
    if(any(diag(L)[!estpars] > 0L) && !list$PLCI){
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
    LL <- LP <- 0
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
    if(all(!est) && list$SE)
        stop('Computing ACOV matrix is meaningless when no parameters are estimated', call.=FALSE)
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
    if(dentype == 'bfactor'){
        Thetabetween <- thetaComb(theta=theta, nfact=nfact-ncol(sitems))
        for(g in 1L:ngroups){
            pars[[g]][[J+1L]]@BFACTOR <- TRUE
            pars[[g]][[J+1L]]@theta <- theta
            pars[[g]][[J+1L]]@Thetabetween <- Thetabetween
        }
        np <- ncol(Thetabetween)
        ns <- ncol(Theta) - np
        gp <- ExtractGroupPars(pars[[1L]][[J+1L]])
        gp$gmeans <- 1L:length(gp$gmeans) - 1L
        ind <- length(gp$gmeans)
        for(i in 1L:length(gp$gmeans)){
            for(j in i:length(gp$gmeans)){
                gp$gcov[j,i] <- ind
                ind <- ind + 1L
            }
        }
        tmp <- gp$gcov[1L:np, 1L:np]
        tmp <- tmp[lower.tri(tmp, TRUE)]
        tmpmat <- matrix(0, ns, 2L)
        for(i in 1L:ns)
            tmpmat[i, ] <- c(gp$gmean[np + i], gp$gcov[np+i, np+i])
        for(g in 1L:ngroups){
            pars[[g]][[J+1L]]@bindex <- as.integer(c(gp$gmeans[1L:np], tmp))
            pars[[g]][[J+1L]]@sindex = tmpmat
        }
    }
    gTheta <- vector('list', ngroups)
    for(g in 1L:ngroups){
        ANY.PRIOR[g] <- any(sapply(pars[[g]], function(x) x@any.prior))
        gTheta[[g]] <- Theta
    }
    preMstep.longpars2 <- preMstep.longpars <- longpars
    is_SEM <- list$SE.type == 'SEM'
    accel <- 0; Mrate <- ifelse(is_SEM, 1, .4)
    Estep.time <- Mstep.time <- 0
    collectLL <- rep(NA, NCYCLES)
    hess <- matrix(0)
    Elist <- list()
    startMrate <- ifelse(Moptim == 'L-BFGS-B', 5L, 1L)
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
                         Theta=Theta, PrepList=PrepList, dentype=dentype, lrPars=lrPars,
                         specific=specific, sitems=sitems, CUSTOM.IND=CUSTOM.IND,
                         constrain=constrain, EHPrior=NULL, Data=Data, method=Moptim,
                         control=control, hessian=list$SE,
                         lower=lower, upper=upper), silent=TRUE)
        cycles <- as.integer(opt$counts[1L])
        longpars[est] <- opt$par
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
        converge <- opt$convergence == 0
        if(list$SE) hess <- opt$hessian
        tmp <- updatePrior(pars=pars, Theta=Theta,
                           list=list, ngroups=ngroups, nfact=nfact,
                           J=J, dentype=dentype, sitems=sitems, cycles=cycles, rlist=rlist)
        prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
        LL <- LP <- 0
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        for(g in 1L:ngroups){
            if(dentype == 'bfactor'){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                            Theta=Theta, prior=prior[[g]],
                                            Priorbetween=Priorbetween[[g]], specific=specific,
                                            sitems=sitems, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                         CUSTOM.IND=CUSTOM.IND, Theta=Theta,
                                         prior=Prior[[g]], itemloc=itemloc)
            }
            LL <- LL + sum(Data$Freq[[g]]*log(rlist[[g]]$expected))
        }
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
        Estep.time <- Estep.time + proc.time()[3L] - start
    } else {
        #EM
        for (cycles in 1L:NCYCLES){
            #priors
            start <- proc.time()[3L]
            if(length(lrPars)) lrPars@mus <- lrPars@X %*% lrPars@beta
            tmp <- updatePrior(pars=pars, Theta=Theta,
                               list=list, ngroups=ngroups, nfact=nfact,
                               J=J, dentype=dentype, sitems=sitems, cycles=cycles,
                               rlist=rlist, full=full, lrPars=lrPars)
            prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
            if(is.na(TOL) && !is.nan(TOL)){
                for(g in 1L:ngroups) rlist[[g]]$expected <- 1
                break
            }
            Elist <- Estep(pars=pars, Data=Data, Theta=Theta, prior=prior, Prior=Prior,
                           Priorbetween=Priorbetween, specific=specific, sitems=sitems,
                           ngroups=ngroups, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                           dentype=dentype, rlist=rlist, full=full, Etable=list$Etable)
            rlist <- Elist$rlist; LL <- Elist$LL
            if(any(ANY.PRIOR)){
                LP <- 0
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
            collectLL[cycles] <- LL
            if(is.nan(LL))
                stop('Optimization error: Could not compute observed log-likelihood. Try
                     estimating with different starting values by passing GenRandomPars = TRUE',
                     call.=FALSE)
            if(!is_SEM){
                if(cycles > startMrate){
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
                if(dentype == 'bfactor'){
                    pars[[g]][[J+1L]]@rrb <- rlist[[g]]$r2
                    pars[[g]][[J+1L]]@rrs <- rlist[[g]]$r3
                } else pars[[g]][[J+1L]]@rr <- rowSums(rlist[[g]]$r1) / J
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
                              dentype=dentype, nfact=nfact, keep_vcov_PD=list$keep_vcov_PD,
                              rlist=rlist, constrain=constrain, DERIV=DERIV, Mrate=Mrate,
                              TOL=list$MSTEPTOL, solnp_args=solnp_args, full=full, lrPars=lrPars,
                              control=control)
            EMhistory[cycles+1L,] <- longpars
            if(verbose)
                cat(sprintf('\rIteration: %d, Log-Lik: %.3f, Max-Change: %.5f',
                            cycles, LL + LP, max(abs(preMstep.longpars - longpars))))
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
                                           dentype=dentype, rlist=rlist, full=full)
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
            for(g in 1L:ngroups){
                if(any(pars[[g]][[J+1L]]@est) && nfact > 1L){
                    gp <- ExtractGroupPars(pars[[g]][[J+1L]])
                    chl <- try(chol(gp$gcov), silent=TRUE)
                    if(is(chl, 'try-error')){
                        if(list$warn)
                            warning('Latent trait variance-covariance matrix became non-positive definite.',
                                    call.=FALSE)
                        converge <- FALSE
                    }
                }
            }
            if(!converge) break
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
            if(list$warn && !(is.nan(TOL) || is.na(TOL)) && !list$NULL.MODEL)
                warning('M-step optimizer converged immediately. Solution is either at the ML or
                     starting values are causing issues and should be adjusted. ', call.=FALSE)
        }
        if(Moptim == 'L-BFGS-B' && cycles <= 10L && !all(!est) && !list$NULL.MODEL){
            if(list$warn && !(is.nan(TOL) || is.na(TOL)) && all( abs(preMstep.longpars - longpars) < 1e-30 ))
                warning(paste0("L-BFGS-B optimizer did not change any values across successive EM cycles;",
                               " likely indicates a problem in the M-step. \nCheck with the more stable ",
                               "optimizer = \'nlminb\', or supply better starting values"), call.=FALSE)
        }
        if(cycles > 1L && list$warn && !ANY.PRIOR){
            diff <- c(-Inf, na.omit(collectLL)) - c(na.omit(collectLL), Inf)
            if(any(diff[length(diff):ceiling(length(diff)*.9)] > .001))
                warning('Log-likelihood was decreasing near the ML solution. EM method may be unstable',
                        call.=FALSE)
        }
    }
    if(dentype == 'custom'){
        if(pars[[1L]][[J + 1L]]@itemclass == -1L){
            for(g in 1L:length(pars)){
                gp <- pars[[g]][[J + 1L]]
                pars[[g]][[J + 1L]]@density <- gp@safe_den(gp, gTheta[[g]])
            }
        }
    }
    infological <- estpars & !redun_constr
    correction <- numeric(length(estpars[estpars & !redun_constr]))
    names(correction) <- names(estpars[estpars & !redun_constr])
    collectLL <- as.numeric(na.omit(collectLL))
    LP <- unname(LP)
    start.time.SE <- proc.time()[3L]
    if(list$SE.type %in% c('SEM', 'Oakes', 'complete') && list$SE){
        h <- matrix(0, nfullpars, nfullpars)
        ind1 <- 1L
        for(group in 1L:ngroups){
            for (i in 1L:J){
                deriv <- Deriv(x=pars[[group]][[i]], Theta=Theta, estHess=TRUE)
                ind2 <- ind1 + length(pars[[group]][[i]]@par) - 1L
                h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
                ind1 <- ind2 + 1L
            }
            deriv <- Deriv(x=pars[[group]][[i+1L]], CUSTOM.IND=CUSTOM.IND,
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
        if(list$SE.type == 'Oakes' && length(lrPars) && list$SE){
            warning('Oakes method not supported for models with latent regression effects', call.=FALSE)
        } else if(list$SE.type == 'Oakes' && list$SE){
            complete_info <- hess
            shortpars <- longpars[estpars & !redun_constr]
            tmp <- updatePrior(pars=pars, Theta=Theta,
                               list=list, ngroups=ngroups, nfact=nfact,
                               J=J, dentype=dentype, sitems=sitems, cycles=cycles,
                               rlist=rlist, full=full, lrPars=lrPars)
            prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
            if(list$Norder >= 2){
                missing_info <- mySapply(1L:length(shortpars), SE.Oakes,
                                       pars=pars, L=L, constrain=constrain, delta=list$delta,
                                       est=est, shortpars=shortpars, longpars=longpars,
                                       Theta=Theta, list=list, ngroups=ngroups, J=J,
                                       dentype=dentype, sitems=sitems, nfact=nfact,
                                       rlist=rlist, full=full, Data=Data,
                                       specific=specific, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                       prior=prior, Priorbetween=Priorbetween, Prior=Prior,
                                       PrepList=PrepList, ANY.PRIOR=ANY.PRIOR, DERIV=DERIV,
                                       SLOW.IND=list$SLOW.IND, Norder=list$Norder)
            } else {
                zero_g <- SE.Oakes(0L, pars=pars, L=L, constrain=constrain, delta=0,
                                   est=est, shortpars=shortpars, longpars=longpars,
                                   Theta=Theta, list=list, ngroups=ngroups, J=J,
                                   dentype=dentype, sitems=sitems, nfact=nfact,
                                   rlist=rlist, full=full, Data=Data,
                                   specific=specific, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                   prior=prior, Priorbetween=Priorbetween, Prior=Prior,
                                   PrepList=PrepList, ANY.PRIOR=ANY.PRIOR, DERIV=DERIV,
                                   SLOW.IND=list$SLOW.IND, Norder=1L)
                missing_info <- mySapply(1L:length(shortpars), SE.Oakes,
                                       pars=pars, L=L, constrain=constrain, delta=list$delta,
                                       est=est, shortpars=shortpars, longpars=longpars,
                                       Theta=Theta, list=list, ngroups=ngroups, J=J,
                                       dentype=dentype, sitems=sitems, nfact=nfact,
                                       rlist=rlist, full=full, Data=Data,
                                       specific=specific, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                       prior=prior, Priorbetween=Priorbetween, Prior=Prior,
                                       PrepList=PrepList, ANY.PRIOR=ANY.PRIOR, DERIV=DERIV,
                                       SLOW.IND=list$SLOW.IND, zero_g=zero_g, Norder=1L)
            }
            if(list$symmetric) missing_info <- (missing_info + t(missing_info))/2
            pars <- reloadPars(longpars=longpars, pars=pars,
                               ngroups=ngroups, J=J)
            is.latent <- grepl('MEAN_', names(shortpars)) | grepl('COV_', names(shortpars))
            missing_info[is.latent, is.latent] <- 0
            hess <- complete_info + missing_info
        }
        ret <- list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                    logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological, Moptim=Moptim,
                    estindex_unique=estindex_unique, correction=correction, hess=hess, Prior=Prior,
                    estpars=estpars & !redun_constr, redun_constr=redun_constr, ngroups=ngroups,
                    LBOUND=LBOUND, UBOUND=UBOUND, EMhistory=na.omit(EMhistory), random=list(),
                    time=c(Estep=as.numeric(Estep.time), Mstep=as.numeric(Mstep.time)),
                    collectLL=collectLL, shortpars=longpars[estpars & !redun_constr],
                    lrPars=lrPars, logPrior=LP, fail_invert_info=FALSE, Etable=Elist$rlist,
                    start.time.SE=start.time.SE)
    } else {
        ret <- list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                    logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological, Moptim=Moptim,
                    estindex_unique=estindex_unique, correction=correction, hess=hess, random=list(),
                    Prior=Prior, time=c(Estep=as.numeric(Estep.time), Mstep=as.numeric(Mstep.time)),
                    prior=prior, Priorbetween=Priorbetween, sitems=sitems, collectLL=collectLL,
                    shortpars=longpars[estpars & !redun_constr], lrPars=lrPars,
                    logPrior=LP, fail_invert_info=FALSE, Etable=Elist$rlist)
    }
    for(g in 1L:ngroups)
        for(i in 1L:J)
            ret$pars[[g]][[i]]@dat <- matrix(0)
    ret
}
