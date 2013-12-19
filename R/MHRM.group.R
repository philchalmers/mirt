MHRM.group <- function(pars, constrain, PrepList, list, random = list(), DERIV)
{
    if(is.null(random)) random <- list()
    RAND <- length(random) > 0L
    verbose <- list$verbose
    nfact <- list$nfact
    NCYCLES <- list$NCYCLES
    BURNIN <- list$BURNIN
    if(RAND) BURNIN <- BURNIN + 50L
    SEMCYCLES <- list$SEMCYCLES
    KDRAWS <- list$KDRAWS
    TOL <- list$TOL
    NO.CUSTOM <- !any(sapply(pars[[1L]], class) %in% 'custom')
    gain <- list$gain
    itemloc <- list$itemloc
    ngroups <- length(pars)
    J <- length(itemloc) - 1L
    prodlist <- PrepList[[1L]]$prodlist
    nfullpars <- 0L
    estpars <- c()
    gfulldata <- gtheta0 <- gstructgrouppars <- gitemtrace <- vector('list', ngroups)
    for(g in 1L:ngroups){
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
        gfulldata[[g]] <- PrepList[[g]]$fulldata
        gtheta0[[g]] <- matrix(0, nrow(gfulldata[[g]]), nfact)
        for(i in 1L:(J+1L)){
            nfullpars <- nfullpars + length(pars[[g]][[i]]@par)
            estpars <- c(estpars, pars[[g]][[i]]@est)
        }
    }
    if(RAND){
        for(i in 1L:length(random)){
            nfullpars <- nfullpars + length(random[[i]]@par)
            estpars <- c(estpars, random[[i]]@est)
        }
    }
    index <- 1L:nfullpars
    N <- nrow(gtheta0[[1L]])
    #Burn in
    cand.t.var <- 1
    tmp <- .1
    OffTerm <- matrix(0, 1, J)
    for(g in 1L:ngroups){
        for(i in 1L:30L){
            gtheta0[[g]] <- draw.thetas(theta0=gtheta0[[g]], pars=pars[[g]], fulldata=gfulldata[[g]],
                                        itemloc=itemloc, cand.t.var=cand.t.var, NO.CUSTOM=NO.CUSTOM,
                                        prior.t.var=gstructgrouppars[[g]]$gcov, OffTerm=OffTerm,
                                        prior.mu=gstructgrouppars[[g]]$gmeans, prodlist=prodlist)
            if(i > 5L){
                if(attr(gtheta0[[g]],"Proportion Accepted") > .35) cand.t.var <- cand.t.var + 2*tmp
                else if(attr(gtheta0[[g]],"Proportion Accepted") > .25 && nfact > 3L)
                    cand.t.var <- cand.t.var + tmp
                else if(attr(gtheta0[[g]],"Proportion Accepted") < .2 && nfact < 4L)
                    cand.t.var <- cand.t.var - tmp
                else if(attr(gtheta0[[g]],"Proportion Accepted") < .1)
                    cand.t.var <- cand.t.var - 2*tmp
                if (cand.t.var < 0){
                    cand.t.var <- tmp
                    tmp <- tmp / 2
                }
            }
        }
    }
    if(RAND) OffTerm <- OffTerm(random, J=J, N=N)
    m.thetas <- grouplist <- SEM.stores <- SEM.stores2 <- m.list <- list()
    conv <- 0L
    k <- 1L
    gamma <- .25
    longpars <- rep(NA,nfullpars)
    ind1 <- 1L
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1L
            longpars[ind1:ind2] <- pars[[g]][[i]]@par
            ind1 <- ind2 + 1L
        }
        if(RAND){
            for(i in 1L:length(random)){
                ind2 <- ind1 + length(random[[i]]@par) - 1L
                longpars[ind1:ind2] <- random[[i]]@par
                ind1 <- ind2 + 1L
            }
        }
    }
    names(longpars) <- names(estpars)
    stagecycle <- 1L
    converge <- 1L
    noninvcount <- 0L
    L <- c()
    for(g in 1L:ngroups)
        for(i in 1L:(J+1L))
            L <- c(L, pars[[g]][[i]]@est)
    if(RAND)
        for(i in 1L:length(random))
            L <- c(L, random[[i]]@est)
    estindex <- index[estpars]
    L <- diag(as.numeric(L))
    redun_constr <- rep(FALSE, length(estpars))
    if(length(constrain) > 0L){
        for(i in 1L:length(constrain)){
            L[constrain[[i]], constrain[[i]]] <- 1L/length(constrain[[i]])
            for(j in 2L:length(constrain[[i]]))
                redun_constr[constrain[[i]][j]] <- TRUE
        }
    }
    estindex_unique <- index[estpars & !redun_constr]
    if(any(diag(L)[!estpars] > 0L)){
        redindex <- index[!estpars]
        stop('Constraint applied to fixed parameter(s) ',
             paste(paste0(redindex[diag(L)[!estpars] > 0L]), ''), ' but should only be applied to
                 estimated parameters. Please fix!')
    }
    #make sure constrained pars are equal
    tmp <- L
    tmp2 <- diag(tmp)
    tmp2[tmp2 == 0L] <- 1L
    diag(tmp) <- tmp2
    longpars <- as.numeric(tmp %*% longpars)
    LBOUND <- UBOUND <- c()
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            LBOUND <- c(LBOUND, pars[[g]][[i]]@lbound)
            UBOUND <- c(UBOUND, pars[[g]][[i]]@ubound)
        }
    }
    if(RAND){
        for(i in 1L:length(random)){
            LBOUND <- c(LBOUND, random[[i]]@lbound)
            UBOUND <- c(UBOUND, random[[i]]@ubound)
        }
    }
    Draws.time <- Mstep.time <- 0

    ####Big MHRM loop
    for(cycles in 1L:(NCYCLES + BURNIN + SEMCYCLES))
    {
        if(cycles == BURNIN + 1L) stagecycle <- 2L
        if(stagecycle == 3L)
            gamma <- (gain[1L] / (cycles - SEMCYCLES - BURNIN - 1L))^(gain[2L])
        if(cycles == (BURNIN + SEMCYCLES + 1L)){
            stagecycle <- 3L
            longpars <- SEM.stores[[1L]]
            Tau <- SEM.stores2[[1L]]
            for(i in 2L:SEMCYCLES){
                longpars <- longpars + SEM.stores[[i]]
                Tau <- Tau + SEM.stores2[[i]]
            }
            longpars <- longpars/SEMCYCLES
            Tau <- Tau/SEMCYCLES
            k <- KDRAWS
            gamma <- .25
        }
        #Reload pars list
        if(list$USEEM) longpars <- list$startlongpars
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        for(g in 1L:ngroups)
            gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
        if(RAND && cycles > 100L) random <- reloadRandom(random=random, longpars=longpars,
                                        parstart=max(pars[[1L]][[J+1L]]@parnum) + 1L)

        start <- proc.time()[3L]
        if(RAND && cycles == 100L){
            for(g in 1L:ngroups) gtheta0[[g]] <- matrix(0, nrow(gfulldata[[g]]), nfact)
            OffTerm <- OffTerm(random, J=J, N=N)
            for(j in 1L:length(random)){
                tmp <- .1
                for(i in 1L:30L){
                    random[[j]]@drawvals <- DrawValues(random[[j]], Theta=gtheta0[[1L]], itemloc=itemloc,
                                                       pars=pars[[1L]], fulldata=gfulldata[[1L]],
                                                       offterm0=OffTerm)
                    OffTerm <- OffTerm(random, J=J, N=N)
                    if(i > 5L){
                        if(attr(random[[j]]@drawvals,"Proportion Accepted") > .4)
                            random[[j]]@cand.t.var <- random[[j]]@cand.t.var + 2*tmp
                        if(attr(random[[j]]@drawvals,"Proportion Accepted") < .2)
                            random[[j]]@cand.t.var <- random[[j]]@cand.t.var - 2*tmp
                        if(attr(random[[j]]@drawvals,"Proportion Accepted") < .05)
                            random[[j]]@cand.t.var <- random[[j]]@cand.t.var - 5*tmp
                        if (random[[j]]@cand.t.var < 0){
                            random[[j]]@cand.t.var <- tmp
                            tmp <- tmp / 10
                        }
                    }
                }
                #better start values
                tmp <- nrow(random[[j]]@drawvals)
                tmp <- cov(random[[j]]@drawvals) * (tmp / (tmp-1L))
                random[[j]]@par[random[[j]]@est] <- tmp[lower.tri(tmp, TRUE)][random[[j]]@est]
            }
            cand.t.var <- .5
            tmp <- .1
            for(i in 1L:30L){
                gtheta0[[1L]] <- draw.thetas(theta0=gtheta0[[1L]], pars=pars[[1L]], fulldata=gfulldata[[1L]],
                                             itemloc=itemloc, cand.t.var=cand.t.var, NO.CUSTOM=NO.CUSTOM,
                                             prior.t.var=gstructgrouppars[[1L]]$gcov, OffTerm=OffTerm,
                                             prior.mu=gstructgrouppars[[1L]]$gmeans, prodlist=prodlist)
                if(i > 5L){
                    if(attr(gtheta0[[g]],"Proportion Accepted") > .35) cand.t.var <- cand.t.var + 2*tmp
                    else if(attr(gtheta0[[g]],"Proportion Accepted") > .25 && nfact > 3L)
                        cand.t.var <- cand.t.var + tmp
                    else if(attr(gtheta0[[g]],"Proportion Accepted") < .2 && nfact < 4L)
                        cand.t.var <- cand.t.var - tmp
                    else if(attr(gtheta0[[g]],"Proportion Accepted") < .1)
                        cand.t.var <- cand.t.var - 2*tmp
                    if (cand.t.var < 0){
                        cand.t.var <- tmp
                        tmp <- tmp / 2
                    }
                }
            }
            tmp <- nrow(gtheta0[[1L]])
            tmp <- cov(gtheta0[[1L]]) * (tmp / (tmp-1L))
            tmp2 <- c(rep(0, ncol(tmp)), tmp[lower.tri(tmp, TRUE)])
            pars[[1L]][[length(pars[[1L]])]]@par[pars[[1L]][[length(pars[[1L]])]]@est] <-
                tmp2[pars[[1L]][[length(pars[[1L]])]]@est]
        }

        #Step 1. Generate m_k datasets of theta
        LL <- 0
        for(g in 1L:ngroups){
            for(i in 1L:5L)
                gtheta0[[g]] <- draw.thetas(theta0=gtheta0[[g]], pars=pars[[g]], fulldata=gfulldata[[g]],
                                      itemloc=itemloc, cand.t.var=cand.t.var, NO.CUSTOM=NO.CUSTOM,
                                      prior.t.var=gstructgrouppars[[g]]$gcov, OffTerm=OffTerm,
                                      prior.mu=gstructgrouppars[[g]]$gmeans, prodlist=prodlist)
            LL <- LL + attr(gtheta0[[g]], "log.lik")
        }
        if(RAND && cycles > 100L){
            for(j in 1:length(random)){
                for(i in 1L:5L){
                    random[[j]]@drawvals <- DrawValues(random[[j]], Theta=gtheta0[[1L]], itemloc=itemloc,
                                                       pars=pars[[1L]], fulldata=gfulldata[[1L]],
                                                       offterm0=OffTerm)
                    OffTerm <- OffTerm(random, J=J, N=N)
                }
            }
        }
        Draws.time <- Draws.time + proc.time()[3L] - start

        #Step 2. Find average of simulated data gradients and hessian
        start <- proc.time()[3L]
        g.m <- h.m <- group.m <- list()
        longpars <- g <- rep(0, nfullpars)
        h <- matrix(0, nfullpars, nfullpars)
        for(group in 1L:ngroups){
            thetatemp <- gtheta0[[group]]
            if(length(prodlist) > 0L) thetatemp <- prodterms(thetatemp,prodlist)
            gitemtrace[[group]] <- computeItemtrace(pars=pars[[group]], offterm=OffTerm,
                                                Theta=thetatemp, itemloc=itemloc, NO.CUSTOM=NO.CUSTOM)
            pars[[group]] <- assignItemtrace(pars=pars[[group]], itemtrace=gitemtrace[[group]],
                                         itemloc=itemloc)
            for (i in 1L:J){
                deriv <- DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=thetatemp,
                               estHess=TRUE, offterm=OffTerm[,i])
                g[pars[[group]][[i]]@parnum] <- deriv$grad
                h[pars[[group]][[i]]@parnum,pars[[group]][[i]]@parnum] <- deriv$hess
                longpars[pars[[group]][[i]]@parnum] <- pars[[group]][[i]]@par
            }
            i <- i + 1L
            deriv <- Deriv(x=pars[[group]][[i]], Theta=gtheta0[[group]])
            longpars[pars[[group]][[i]]@parnum] <- pars[[group]][[i]]@par
            g[pars[[group]][[i]]@parnum] <- deriv$grad
            h[pars[[group]][[i]]@parnum, pars[[group]][[i]]@parnum] <- deriv$hess
        }
        if(RAND){
            for(i in 1L:length(random)){
                deriv <- RandomDeriv(x=random[[i]])
                longpars[random[[i]]@parnum] <- random[[i]]@par
                g[random[[i]]@parnum] <- deriv$grad
                h[random[[i]]@parnum, random[[i]]@parnum] <- deriv$hess
            }
        }
        grad <- g %*% L
        ave.h <- (-1)*L %*% h %*% L
        grad <- grad[1L, estpars & !redun_constr]
        ave.h <- ave.h[estpars & !redun_constr, estpars & !redun_constr]
        if(any(is.na(grad)))
            stop('Model did not converge (unacceptable gradient caused by extreme parameter values)')
        if(is.na(attr(gtheta0[[1L]],"log.lik")))
            stop('Estimation halted. Model did not converge.')
        if(verbose){
            if(cycles <= BURNIN)
                printmsg <- sprintf("\rStage 1: Cycle = %i, Log-Lik = %.1f", cycles, LL)
            if(cycles > BURNIN && cycles <= BURNIN + SEMCYCLES)
                printmsg <- sprintf("\rStage 2: Cycle = %i, Log-Lik = %.1f", cycles-BURNIN, LL)
            if(cycles > BURNIN + SEMCYCLES)
                printmsg <- sprintf("\rStage 3: Cycle = %i, Log-Lik = %.1f", cycles-BURNIN-SEMCYCLES, LL)
        }
        if(stagecycle < 3L){
            if(qr(ave.h)$rank != ncol(ave.h)){
                ev <- eigen(ave.h)
                eval <- ev$values
                eval[eval < 0] <- 100*.Machine$double.eps
                eval <- eval / sum(eval) * sum(ev$values)
                ave.h <- ev$vectors %*% diag(eval) %*% t(ev$vectors)
                noninvcount <- noninvcount + 1L
                if(noninvcount == 3L)
                    stop('\nEstimation halted during burn in stages, solution is unstable')
            }
            correction <- solve(ave.h, grad)
            correction[correction > 1] <- 1
            correction[correction < -1] <- -1
            longpars[estindex_unique] <- longpars[estindex_unique] + gamma*correction
            longpars[longpars < LBOUND] <- LBOUND[longpars < LBOUND]
            longpars[longpars > UBOUND] <- UBOUND[longpars > UBOUND]
            if(length(constrain) > 0L)
                for(i in 1L:length(constrain))
                    longpars[index %in% constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
            if(verbose)
                cat(printmsg, sprintf(", Max Change = %.4f\r", max(abs(gamma*correction))), sep='')
            if(stagecycle == 2L){
                SEM.stores[[cycles - BURNIN]] <- longpars
                SEM.stores2[[cycles - BURNIN]] <- ave.h
            }
            Mstep.time <- Mstep.time + proc.time()[3L] - start
            next
        }

        #Step 3. Update R-M step
        Tau <- Tau + gamma*(ave.h - Tau)
        if(qr(Tau)$rank != ncol(Tau)){
            ev <- eigen(Tau)
            eval <- ev$values
            eval[eval < 0] <- 100*.Machine$double.eps
            eval <- eval / sum(eval) * sum(ev$values)
            Tau <- ev$vectors %*% diag(eval) %*% t(ev$vectors)
            noninvcount <- noninvcount + 1L
            if(noninvcount == 3L)
                stop('\nEstimation halted during burn in stages, solution is unstable')
        }
        correction <- solve(Tau, grad)
        correction[gamma*correction > .25] <- .25/gamma
        correction[gamma*correction < -.25] <- -.25/gamma
        longpars[estindex_unique] <- longpars[estindex_unique] + gamma*correction
        longpars[longpars < LBOUND] <- LBOUND[longpars < LBOUND]
        longpars[longpars > UBOUND] <- UBOUND[longpars > UBOUND]
        if(length(constrain) > 0L)
            for(i in 1L:length(constrain))
                longpars[index %in% constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
        if(verbose)
            cat(printmsg, sprintf(", gam = %.3f, Max Change = %.4f\r",
                                  gamma, max(abs(gamma*correction))), sep='')
        if(all(abs(gamma*correction) < TOL)) conv <- conv + 1L
        else conv <- 0L
        if(!list$SE && conv == 3L) break
        if(list$SE && cycles >= (400L + BURNIN + SEMCYCLES) && conv == 3L) break

        #Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE
        if(gamma == .25){
            gamma <- 0
            phi <- grad
            Phi <- Tau
        }
        phi <- phi + gamma*(grad - phi)
        Phi <- Phi + gamma*(ave.h - outer(grad,grad) - Phi)
        Mstep.time <- Mstep.time + proc.time()[3L] - start
    } ###END BIG LOOP
    if(verbose) cat('\r\n')
    info <- Phi - outer(phi,phi)
    diag(info) <- abs(diag(info)) #diag of latent variances neg sometimes, why?
    #Reload final pars list
    if(cycles == NCYCLES + BURNIN + SEMCYCLES && !list$USEEM){
        message('MHRM iterations terminated after ', NCYCLES, ' iterations.')
        converge <- 0L
    }
    if(list$USEEM) longpars <- list$startlongpars
    ind1 <- 1L
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1L
            pars[[g]][[i]]@par <- longpars[ind1:ind2]
            ind1 <- ind2 + 1L
        }
    }
    SEtmp <- abs(diag(qr.solve(info)))
    if(any(SEtmp < 0)){
        warning("Negative SEs set to NaN.\n")
        SEtmp[SEtmp < 0 ] <- NaN
    }
    SEtmp <- sqrt(SEtmp)
    SE <- rep(NA, length(longpars))
    SE[estindex_unique] <- SEtmp
    if(length(constrain) > 0L)
        for(i in 1L:length(constrain))
            SE[index %in% constrain[[i]][-1L]] <- SE[constrain[[i]][1L]]
    ind1 <- 1L
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1L
            pars[[g]][[i]]@SEpar <- SE[ind1:ind2]
            ind1 <- ind2 + 1L
        }
    }
    if(RAND){
        for(i in 1L:length(random)){
            ind2 <- ind1 + length(random[[i]]@par) - 1L
            random[[i]]@SEpar <- SE[ind1:ind2]
            ind1 <- ind2 + 1L
        }
    }
    names(correction) <- names(estpars)[estindex_unique]
    info <- nameInfoMatrix(info=info, correction=correction, L=L, npars=length(longpars))
    ret <- list(pars=pars, cycles = cycles - BURNIN - SEMCYCLES, info=as.matrix(info),
                longpars=longpars, converge=converge, SElogLik=0, cand.t.var=cand.t.var,
                random=random, time=c(MH_draws = as.numeric(Draws.time), Mstep=as.numeric(Mstep.time)))
    ret
}
