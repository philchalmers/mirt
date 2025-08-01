EM.group <- function(pars, constrain, Ls, Data, PrepList, list, Theta, DERIV, solnp_args, control,
                     nconstrain=NULL, fixedEtable=NULL)
{
    verbose <- list$verbose
    lrPars <- list$lrPars
    nfact <- list$nfact
    if(list$dentype %in% c('EH', 'EHW') && nfact != 1L)
        stop('empirical histogram only available for unidimensional models', call.=FALSE)
    if(list$dentype == 'Davidian'){
        if(nfact != 1L)
            stop('Davidian curves estimation only availablae for unidimensional models', call.=FALSE)
        if(list$SE)
            stop('No standard error method currently supported for Davidian curves', call.=FALSE)
        J <- length(pars[[1L]]) - 1L
        for(g in seq_len(length(pars))){ # throw error if latent mean/var estimated TODO
            if(any(pars[[g]][[J+1L]]@est[1L:2L]))
                stop('Estimating Davidian mean-variance hyper-parameters is not currently supported',
                     call.=FALSE)
        }
    }
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
    MC <- list$method %in% c('QMCEM', 'MCEM')
    QMC <- list$method == 'QMCEM'
    for(g in seq_len(ngroups)){
        for(i in seq_len(J+1L)){
            nfullpars <- nfullpars + length(pars[[g]][[i]]@par)
            estpars <- c(estpars, pars[[g]][[i]]@est)
        }
        if(length(lrPars)){
            nfullpars <- nfullpars + length(lrPars@par)
            estpars <- c(estpars, lrPars@est)
        }
    }
    listpars <- vector('list', ngroups)
    for(g in seq_len(ngroups)){
        listpars[[g]] <- list()
        for(i in seq_len(J + 1L)){
            listpars[[g]][[i]] <- pars[[g]][[i]]@par
        }
        if(length(lrPars))
            listpars[[g]][[i+1L]] <- lrPars@par
    }
    index <- seq_len(nfullpars)
    longpars <- rep(NA,nfullpars)
    latent_longpars <- logical(nfullpars)
    ind1 <- 1L
    for(g in seq_len(ngroups)){
        for(i in seq_len(J+1L)){
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
    if(any(attr(L, 'diag')[!estpars] > 0L) && !list$PLCI){
        redindex <- index[!estpars]
        stop('Constraint applied to fixed parameter(s) ',
             paste(paste0(redindex[attr(L, 'diag')[!estpars] > 0L]), ''), ' but should only be applied to
                 estimated parameters. Please fix!', call.=FALSE)
    }
    prior <- rlist <- r <- vector('list', ngroups)
    #make sure constrained pars are equal
    tmp <- Matrix::rowSums(L)
    tmp[tmp == 0L] <- 1L
    check <- as.numeric(L %*% longpars) / tmp
    longpars[estpars] <- check[estpars]
    LL <- LP <- 0
    LBOUND <- UBOUND <- c()
    for(g in seq_len(ngroups)){
        for(i in seq_len(J+1L)){
            LBOUND <- c(LBOUND, pars[[g]][[i]]@lbound)
            UBOUND <- c(UBOUND, pars[[g]][[i]]@ubound)
        }
        if(length(lrPars)){
            LBOUND <- c(LBOUND, lrPars@lbound)
            UBOUND <- c(UBOUND, lrPars@ubound)
        }
    }
    est <- c()
    for(g in seq_len(ngroups)){
       for(j in seq_len(J+1L))
           est <- c(est, pars[[g]][[j]]@est)
       if(length(lrPars)){
           tmp <- rep(FALSE, length(lrPars@est))
           names(tmp) <- names(lrPars@est)
           est <- c(est, tmp)
       }
    }
    for(i in seq_len(length(constrain)))
        est[constrain[[i]][-1L]] <- FALSE
    if(list$dentype %in% c('Davidian', 'EHW')) # dont estimate mean/var for extrapolation
        est[names(est) %in% c('MEAN_1', 'COV_11')] <- FALSE
    if(all(!est) && list$SE)
        stop('Computing ACOV matrix is meaningless when no parameters are estimated', call.=FALSE)
    names(longpars) <- names(est)
    if(list$Moptim != 'BFGS') {
        Moptim <- list$Moptim
    } else {
        Moptim <- if(all(c(LBOUND[est], UBOUND[est]) %in% c(-Inf, Inf))) 'BFGS' else 'nlminb'
    }
    if(Moptim == 'NR' && sum(est) > 300L && list$message)
        message('NR optimizer should not be used for models with a large number of parameters.
                Use the optimizer = \'BFGS\' or \'nlminb\' instead.')
    EMhistory <- matrix(NA, NCYCLES+1L, length(longpars))
    EMhistory[1L,] <- longpars
    ANY.PRIOR <- rep(FALSE, ngroups)
    if(length(prodlist))
        Theta <- prodterms(Theta, prodlist)
    if(dentype == 'bfactor'){
        Thetabetween <- thetaComb(theta=theta, nfact=nfact-ncol(sitems))
        for(g in seq_len(ngroups)){
            pars[[g]][[J+1L]]@BFACTOR <- TRUE
            pars[[g]][[J+1L]]@theta <- theta
            pars[[g]][[J+1L]]@Thetabetween <- Thetabetween
        }
        np <- ncol(Thetabetween)
        ns <- ncol(Theta) - np
        gp <- ExtractGroupPars(pars[[1L]][[J+1L]])
        gp$gmeans <- seq_len(length(gp$gmeans)) - 1L
        ind <- length(gp$gmeans)
        for(i in seq_len(length(gp$gmeans))){
            for(j in i:length(gp$gmeans)){
                gp$gcov[j,i] <- ind
                ind <- ind + 1L
            }
        }
        tmp <- gp$gcov[seq_len(np), seq_len(np)]
        tmp <- tmp[lower.tri(tmp, TRUE)]
        tmpmat <- matrix(0, ns, 2L)
        for(i in seq_len(ns))
            tmpmat[i, ] <- c(gp$gmeans[np + i], gp$gcov[np+i, np+i])
        for(g in seq_len(ngroups)){
            pars[[g]][[J+1L]]@bindex <- as.integer(c(gp$gmeans[seq_len(np)], tmp))
            pars[[g]][[J+1L]]@sindex = tmpmat
        }
    }
    gTheta <- vector('list', ngroups)
    for(g in seq_len(ngroups)){
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
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain,
                                   nconstrain=nconstrain)
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    if(list$method == 'BL'){
        start <- proc.time()[3L]
        lower <- LBOUND[est]; upper <- UBOUND[est]
        Moptim <- ifelse(any(is.finite(lower) | is.finite(upper)), 'nlminb', 'BFGS')
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
                         constrain=constrain, nconstrain=nconstrain,
                         EHPrior=NULL, Data=Data, omp_threads=list$omp_threads,
                         method=Moptim, control=control, hessian=list$SE,
                         lower=lower, upper=upper), silent=TRUE)
        cycles <- as.integer(opt$counts[1L])
        longpars[est] <- opt$par
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain,
                                       nconstrain=nconstrain)
        converge <- opt$convergence == 0
        if(list$SE) hess <- opt$hessian
        tmp <- updatePrior(pars=pars, gTheta=gTheta, MC=MC,
                           list=list, ngroups=ngroups, nfact=nfact,
                           J=J, dentype=dentype, sitems=sitems, cycles=cycles, rlist=rlist)
        prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
        LL <- LP <- 0
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        for(g in seq_len(ngroups)){
            if(dentype == 'bfactor'){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                            Theta=Theta, prior=prior[[g]], wmiss=Data$wmiss,
                                            Priorbetween=Priorbetween[[g]], specific=specific,
                                            sitems=sitems, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                            omp_threads=list$omp_threads)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                         CUSTOM.IND=CUSTOM.IND, Theta=Theta, wmiss=Data$wmiss,
                                         prior=Prior[[g]], itemloc=itemloc, omp_threads=list$omp_threads)
            }
            LL <- LL + sum(Data$Freq[[g]]*log(rlist[[g]]$expected))
        }
        if(any(ANY.PRIOR)){
            if(length(lrPars)){
                if(lrPars@any.prior)
                    LP <- LL.Priors(x=lrPars, LL=LP)
            }
            for(g in seq_len(length(pars))){
                for(i in seq_len(length(pars[[1L]])))
                    if(pars[[g]][[i]]@any.prior)
                        LP <- LL.Priors(x=pars[[g]][[i]], LL=LP)
            }
        }
        Estep.time <- Estep.time + proc.time()[3L] - start
    } else {
        #EM
        for (cycles in seq_len(NCYCLES)){
            #priors
            start <- proc.time()[3L]
            if(length(lrPars)) lrPars@mus <- lrPars@X %*% lrPars@beta
            if(MC)
                gTheta <- updateTheta(npts=if(QMC) list$quadpts else list$MCEM_draws(cycles),
                                      nfact=nfact, pars=pars, QMC=QMC)
            tmp <- updatePrior(pars=pars, gTheta=gTheta,
                               list=list, ngroups=ngroups, nfact=nfact,
                               J=J, dentype=dentype, sitems=sitems, cycles=cycles,
                               rlist=rlist, full=full, lrPars=lrPars, MC=MC)
            prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
            if(is.na(TOL) && !is.nan(TOL)){
                for(g in seq_len(ngroups)) rlist[[g]]$expected <- 1
                break
            }
            Elist <- Estep(pars=pars, Data=Data, gTheta=gTheta, prior=prior, Prior=Prior,
                           Priorbetween=Priorbetween, specific=specific, sitems=sitems,
                           ngroups=ngroups, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                           dentype=dentype, rlist=rlist, full=full, Etable=list$Etable,
                           omp_threads=list$omp_threads)
            if(!is.null(fixedEtable)){
                if(!(Moptim %in% c('BFGS', 'L-BFGS-B', 'nlminb')))
                    stop('Optimizer not supported when using fixedEtable input', call.=FALSE)
                Elist$rlist[[1]]$r1 <- fixedEtable
                if(Moptim %in% c('BFGS', 'L-BFGS-B'))
                    control$maxit <- 200
                else if(Moptim == 'nlminb'){
                    control$rel.tol <- 1e-10
                    control$iter.max <- 200
                }
            }
            rlist <- Elist$rlist; LL <- Elist$LL
            if(any(ANY.PRIOR)){
                LP <- 0
                if(length(lrPars)){
                    if(lrPars@any.prior)
                        LP <- LL.Priors(x=lrPars, LL=LP)
                }
                for(g in seq_len(length(pars))){
                    for(i in seq_len(length(pars[[1L]])))
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
            for(g in seq_len(ngroups)){
                for(i in seq_len(J))
                    pars[[g]][[i]]@dat <- rlist[[g]]$r1[, c(itemloc[i]:(itemloc[i+1L] - 1L)),
                                                        drop=FALSE]
                if(dentype == 'bfactor'){
                    pars[[g]][[J+1L]]@rrb <- rlist[[g]]$r2
                    pars[[g]][[J+1L]]@rrs <- rlist[[g]]$r3
                } else {
                    pars[[g]][[J+1L]]@rr <- rlist[[g]]$r1g
                    if(dentype == 'EHW' && g == 1L || dentype == "Davidian" ||
                       (dentype == 'custom' && pars[[g]][[J+1L]]@standardize))
                        pars[[g]][[J+1L]]@rr <- standardizeQuadrature(gTheta[[g]],
                                                        nq=pars[[g]][[J+1L]]@rr,
                                                        estmean=pars[[g]][[J+1L]]@est['MEAN_1'],
                                                        estsd=pars[[g]][[J+1L]]@est['COV_11'])
                }
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
                              rlist=rlist, constrain=constrain, nconstrain=nconstrain, DERIV=DERIV, Mrate=Mrate,
                              TOL=list$MSTEPTOL, solnp_args=solnp_args, full=full, lrPars=lrPars,
                              control=control)
            if(dentype == 'EHW' && cycles > 1L){
                for(g in seq_len(ngroups)){
                    tmpparnum <- pars[[g]][[J+1L]]@parnum
                    if(pars[[g]][[J+1L]]@est[1L])
                        longpars[tmpparnum[1L]] <- attr(Prior[[g]], 'mean_var')['mean']
                    if(pars[[g]][[J+1L]]@est[2L])
                        longpars[tmpparnum[2L]] <- attr(Prior[[g]], 'mean_var')['var']
                }
            }
            if(verbose && interactive())
                printf('\rIteration: %d, Log-Lik: %.3f, Max-Change: %.5f',
                            cycles, LL + LP, max(abs(preMstep.longpars - longpars)))

            if(hasConverged(preMstep.longpars, longpars, TOL) || !is.null(fixedEtable)){
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
                        ratio <- as.vector(sqrt((dX %*% dX) / (d2X2 %*% d2X2)))
                        accel <- 1 - ratio
                        if(accel < -5) accel <- -5
                        tmp <- (1 - accel) * longpars + accel * preMstep.longpars
                        longpars[!latent_longpars] <- tmp[!latent_longpars]
                    }
                } else if(list$accelerate == 'squarem'){
                    r <- preMstep.longpars - preMstep.longpars2
                    v <- (longpars - preMstep.longpars) - r
                    ratio <- as.vector(sqrt((r %*% r) / (v %*% v)))
                    accel <- -ratio
                    if(accel > -1){
                        accel <- -1
                    } else {
                        count <- 1L
                        while(TRUE){
                            tmp <- preMstep.longpars2 - 2 * accel * r  + accel^2 * v
                            longpars[!latent_longpars] <- tmp[!latent_longpars]
                            pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
                            Elist <- Estep(pars=pars, Data=Data, gTheta=gTheta, prior=prior, Prior=Prior,
                                           Priorbetween=Priorbetween, specific=specific, sitems=sitems,
                                           ngroups=ngroups, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                           dentype=dentype, rlist=rlist, full=full, omp_threads=list$omp_threads)
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
            EMhistory[cycles+1L,] <- longpars
            pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
            for(g in seq_len(ngroups)){
                if(any(pars[[g]][[J+1L]]@est) && nfact > 1L){
                    gp <- ExtractGroupPars(pars[[g]][[J+1L]])
                    ev <- eigen(gp$gcov)
                    if(any(ev$values <= sqrt(.Machine$double.eps))){
                        warning('Latent trait covariance matrix became non-positive definite. Model may be unstable', call.=FALSE)
                        sds <- diag(sqrt(diag(gp$gcov)))
                        smoothed <- cov2cor(smooth.cov(gp$gcov))
                        gcov <- sds %*% smoothed %*% sds
                        longpars[pars[[g]][[J+1]]@parnum] <- c(pars[[g]][[J+1L]]@par[1:nfact],
                                                               gcov[lower.tri(gcov, TRUE)])
                        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
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
        if(verbose && !list$SE) cat('\n')
        if(is.nan(TOL) || is.numeric(TOL)){
            if(length(lrPars)) lrPars@mus <- lrPars@X %*% lrPars@beta
            if(MC)
                gTheta <- updateTheta(npts=if(QMC) list$quadpts else list$MCEM_draws(cycles),
                                      nfact=nfact, pars=pars, QMC=QMC)
            tmp <- updatePrior(pars=pars, gTheta=gTheta,
                               list=list, ngroups=ngroups, nfact=nfact,
                               J=J, dentype=dentype, sitems=sitems, cycles=2L,
                               rlist=rlist, full=full, lrPars=lrPars, MC=MC)
            prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
            Elist <- Estep(pars=pars, Data=Data, gTheta=gTheta, prior=prior, Prior=Prior,
                           Priorbetween=Priorbetween, specific=specific, sitems=sitems,
                           ngroups=ngroups, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                           dentype=dentype, rlist=rlist, full=full, Etable=list$Etable,
                           omp_threads=list$omp_threads)
            rlist <- Elist$rlist; LL <- Elist$LL
        }
        if(cycles == NCYCLES){
            if(list$warn)
                warning('EM cycles terminated after ', cycles, ' iterations.', call.=FALSE)
            converge <- FALSE
        } else if(cycles == 1L && !all(!est)){
            if(list$warn && !(is.nan(TOL) || is.na(TOL)) && !list$NULL.MODEL && is.null(fixedEtable))
                warning('M-step optimizer converged immediately. Estimates are either at the ML or
                     starting values are causing issues and should be adjusted. ', call.=FALSE)
        }
        if(Moptim == 'L-BFGS-B' && cycles <= 10L && !all(!est) && !list$NULL.MODEL){
            if(list$warn && !(is.nan(TOL) || is.na(TOL)) && all( abs(preMstep.longpars - longpars) < 1e-30 ))
                warning(paste0("L-BFGS-B optimizer did not change any values across successive EM cycles;",
                               " likely indicates a problem in the M-step. \nCheck with the more stable ",
                               "optimizer = \'nlminb\', or supply better starting values"), call.=FALSE)
        }
        if(cycles > 1L && list$warn && !any(ANY.PRIOR) &&
           list$method != 'MCEM' && !dentype %in% c('EHW', 'Davidian')){
            diff <- c(-Inf, na.omit(collectLL)) - c(na.omit(collectLL), Inf)
            if(any(diff[length(diff):ceiling(length(diff)*.9)] > .001))
                warning('Log-likelihood was decreasing near the ML solution. EM method may be unstable',
                        call.=FALSE)
        }
    }
    if(dentype %in% c('custom', 'discrete')){
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
    if(list$SE.type %in% c('SEM', 'Oakes', 'complete', 'sandwich', 'Louis', 'sandwich.Louis') && list$SE){
        h <- matrix(0, nfullpars, nfullpars)
        ind1 <- 1L
        usefixed <- pars[[1L]][[1]]@nfixedeffects > 0
        for(group in seq_len(ngroups)){
            for (i in seq_len(J)){
                Thetas <- gTheta[[group]]
                if(usefixed)
                    Thetas <- cbind(pars[[group]][[i]]@fixed.design[rep(1, nrow(Thetas)),],
                                    Thetas)
                deriv <- Deriv(x=pars[[group]][[i]], Theta=Thetas, estHess=TRUE)
                ind2 <- ind1 + length(pars[[group]][[i]]@par) - 1L
                h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
                ind1 <- ind2 + 1L
            }
            deriv <- Deriv(x=pars[[group]][[i+1L]], CUSTOM.IND=CUSTOM.IND,
                           Theta=gTheta[[group]], EM = TRUE,
                           pars=pars[[group]], tabdata=Data$tabdatalong,
                           freq=Data$Freq[[group]], prior=Prior[[group]],
                           itemloc=itemloc,
                           bfactor_info=if(dentype == 'bfactor')
                               list(specific=specific, sitems=sitems, nfact=nfact) else NULL,
                           estHess=TRUE)
            ind2 <- ind1 + length(deriv$grad) - 1L
            h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
            ind1 <- ind2 + 1L
            if(length(lrPars)){
                gp <- ExtractGroupPars(pars[[group]][[J+1L]])
                tmp <- Mstep.LR(Theta=gTheta[[group]], CUSTOM.IND=CUSTOM.IND, pars=pars[[group]],
                                itemloc=itemloc, fulldata=Data$fulldata[[1L]], prior=Prior[[group]],
                                lrPars=lrPars, retscores=TRUE)
                deriv <- Deriv(lrPars, cov=gp$gcov, theta=tmp)
                for(i in 0L:(ncol(deriv$grad)-1L))
                    h[lrPars@parnum[1L:nrow(deriv$grad) + nrow(deriv$grad)*i],
                      lrPars@parnum[1L:nrow(deriv$grad) + nrow(deriv$grad)*i]] <- deriv$hess
                est[lrPars@parnum] <- lrPars@est
            }
        }
        if(dentype == 'mixture'){
            mixtype <- list(par=sapply(1L:length(pars),
                                       function(g) sum(pars[[g]][[J+1L]]@par[length(pars[[g]][[J+1L]]@par)])),
                            est=as.logical(sapply(1L:length(pars),
                                                  function(g) sum(pars[[g]][[J+1L]]@est[length(pars[[g]][[J+1L]]@est)]))),
                            parnum=sapply(1L:length(pars),
                                          function(g) sum(pars[[g]][[J+1L]]@parnum[length(pars[[g]][[J+1L]]@parnum)])),
                            any.prior=FALSE,
                            dat=matrix(sapply(1L:length(pars),
                                              function(g) sum(pars[[g]][[J+1L]]@rr)), 1L))
            deriv <- Deriv.mix(mixtype, estHess=TRUE)
            h[mixtype$parnum, mixtype$parnum] <- deriv$hess
        } else mixtype <- NULL
        hess <- updateHess(h=h, L=L)
        hess <- as.matrix(hess[estpars & !redun_constr, estpars & !redun_constr])
        # if(list$SE.type %in% c('Oakes', 'sandwich') && length(lrPars) && list$SE){
        #     warning('Oakes method not supported for models with latent regression effects', call.=FALSE)
        # } else
        if(list$SE.type %in% c('Oakes', 'sandwich') && list$SE){
            if(verbose)
                catf('\n\nCalculating information matrix...\n')
            complete_info <- hess
            shortpars <- longpars[estpars & !redun_constr]
            tmp <- updatePrior(pars=pars, gTheta=gTheta,
                               list=list, ngroups=ngroups, nfact=nfact,
                               J=J, dentype=dentype, sitems=sitems, cycles=cycles,
                               rlist=rlist, full=full, lrPars=lrPars, MC=MC)
            prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
            if(list$Norder >= 2){
                missing_info <- mySapply(seq_len(length(shortpars)), SE.Oakes,
                                       pars=pars, L=L, constrain=constrain, delta=list$delta,
                                       est=est, shortpars=shortpars, longpars=longpars,
                                       gTheta=gTheta, list=list, ngroups=ngroups, J=J,
                                       dentype=dentype, sitems=sitems, nfact=nfact, lrPars=lrPars,
                                       rlist=rlist, full=full, Data=Data, mixtype=mixtype,
                                       specific=specific, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                       prior=prior, Priorbetween=Priorbetween, Prior=Prior,
                                       PrepList=PrepList, ANY.PRIOR=ANY.PRIOR, DERIV=DERIV,
                                       SLOW.IND=list$SLOW.IND, Norder=list$Norder, omp_threads=list$omp_threads)
            } else {
                zero_g <- SE.Oakes(0L, pars=pars, L=L, constrain=constrain, delta=0,
                                   est=est, shortpars=shortpars, longpars=longpars,
                                   gTheta=gTheta, list=list, ngroups=ngroups, J=J,
                                   dentype=dentype, sitems=sitems, nfact=nfact, lrPars=lrPars,
                                   rlist=rlist, full=full, Data=Data, mixtype=mixtype,
                                   specific=specific, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                   prior=prior, Priorbetween=Priorbetween, Prior=Prior,
                                   PrepList=PrepList, ANY.PRIOR=ANY.PRIOR, DERIV=DERIV,
                                   SLOW.IND=list$SLOW.IND, Norder=1L, omp_threads=list$omp_threads)
                missing_info <- mySapply(seq_len(length(shortpars)), SE.Oakes,
                                       pars=pars, L=L, constrain=constrain, delta=list$delta,
                                       est=est, shortpars=shortpars, longpars=longpars,
                                       gTheta=gTheta, list=list, ngroups=ngroups, J=J,
                                       dentype=dentype, sitems=sitems, nfact=nfact,
                                       rlist=rlist, full=full, Data=Data, mixtype=mixtype,
                                       specific=specific, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                       prior=prior, Priorbetween=Priorbetween, Prior=Prior,
                                       PrepList=PrepList, ANY.PRIOR=ANY.PRIOR, DERIV=DERIV,
                                       SLOW.IND=list$SLOW.IND, zero_g=zero_g, Norder=1L, omp_threads=list$omp_threads)
            }
            if(list$symmetric) missing_info <- (missing_info + t(missing_info))/2
            pars <- reloadPars(longpars=longpars, pars=pars,
                               ngroups=ngroups, J=J)
            is.latent <- grepl('MEAN_', names(shortpars)) | grepl('COV_', names(shortpars))
            if(length(lrPars))
                is.latent[names(shortpars) %in% names(lrPars@est)] <- TRUE
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
                    start.time.SE=start.time.SE, Theta=gTheta[[1L]], sitems=sitems)
    } else {
        ret <- list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                    logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological, Moptim=Moptim,
                    estindex_unique=estindex_unique, correction=correction, hess=hess, random=list(),
                    Prior=Prior, time=c(Estep=as.numeric(Estep.time), Mstep=as.numeric(Mstep.time)),
                    prior=prior, Priorbetween=Priorbetween, sitems=sitems, collectLL=collectLL,
                    shortpars=longpars[estpars & !redun_constr], lrPars=lrPars, EMhistory=na.omit(EMhistory),
                    logPrior=LP, fail_invert_info=FALSE, Etable=Elist$rlist, Theta=gTheta[[1L]])
    }
    attr(ret$EMhistory, "na.action") <- NULL
    for(g in seq_len(ngroups))
        for(i in seq_len(J))
            ret$pars[[g]][[i]]@dat <- matrix(0)
    ret
}
