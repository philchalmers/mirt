MHRM.group <- function(pars, constrain, Ls, Data, PrepList, list, random = list(),
                       lrPars = list(), DERIV)
{
    if(is.null(random)) random <- list()
    itemtype <- sapply(pars[[1]], class)
    has_graded <- any(itemtype == 'graded')
    lr.random <- list() #FIXME overwrite later
    RAND <- length(random) > 0L; LR.RAND <- length(lr.random) > 0L
    RANDSTART <- list$RANDSTART
    LRPARS <- length(lrPars) > 0L
    verbose <- list$verbose
    nfact <- list$nfact
    NCYCLES <- list$NCYCLES
    BURNIN <- list$BURNIN
    MHDRAWS <- list$MHDRAWS
    stopifnot(BURNIN >= RANDSTART)
    SEMCYCLES <- list$SEMCYCLES
    KDRAWS <- list$KDRAWS
    TOL <- list$TOL
    CUSTOM.IND <- list$CUSTOM.IND
    USE.FIXED <- nrow(pars[[1L]][[1L]]@fixed.design) > 1L
    gain <- list$gain
    itemloc <- list$itemloc
    ngroups <- length(pars)
    J <- length(itemloc) - 1L
    prodlist <- PrepList[[1L]]$prodlist
    nfullpars <- 0L
    estpars <- c()
    gtheta0 <- gstructgrouppars <- vector('list', ngroups)
    for(g in 1L:ngroups){
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
        gtheta0[[g]] <- matrix(0, nrow(Data$fulldata[[g]]), nfact)
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
    if(LRPARS){
        nfullpars <- nfullpars + length(lrPars@par)
        estpars <- c(estpars, lrPars@est)
    }
    index <- 1L:nfullpars
    N <- nrow(gtheta0[[1L]])
    #Burn in
    if(LRPARS){
        lrPars@beta[] <- lrPars@par
        lrPars@mus <- lrPars@X %*% lrPars@beta
        gstructgrouppars[[1L]]$gmeans <- lrPars@mus
    }
    cand.t.var <- 1
    tmp <- .1
    OffTerm <- matrix(0, 1, J)
    for(g in 1L:ngroups){
        for(i in 1L:31L){
            gtheta0[[g]] <- draw.thetas(theta0=gtheta0[[g]], pars=pars[[g]], fulldata=Data$fulldata[[g]],
                                        itemloc=itemloc, cand.t.var=cand.t.var, CUSTOM.IND=CUSTOM.IND,
                                        prior.t.var=gstructgrouppars[[g]]$gcov, OffTerm=OffTerm,
                                        prior.mu=gstructgrouppars[[g]]$gmeans, prodlist=prodlist)
            if(is.null(list$cand.t.var)){
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
            } else {
                cand.t.var <- list$cand.t.var[1L]
            }
        }
    }
    if(RAND) OffTerm <- OffTerm(random, J=J, N=N)
    m.thetas <- grouplist <- SEM.stores <- SEM.stores2 <- m.list <- list()
    conv <- 0L
    k <- 1L
    gamma <- .25
    longpars <- rep(NA,nfullpars)
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L))
            longpars[pars[[g]][[i]]@parnum] <- pars[[g]][[i]]@par
        if(RAND){
            for(i in 1L:length(random))
                longpars[random[[i]]@parnum] <- random[[i]]@par
        }
        if(LRPARS)
            longpars[lrPars@parnum] <- lrPars@par
    }
    names(longpars) <- names(estpars)
    stagecycle <- 1L
    converge <- 1L
    noninvcount <- 0L
    estindex <- index[estpars]
    L <- Ls$L
    redun_constr <- Ls$redun_constr
    estindex_unique <- index[estpars & !redun_constr]
    if(any(diag(L)[!estpars] > 0L)){
        redindex <- index[!estpars]
        stop('Constraint applied to fixed parameter(s) ',
             paste(paste0(redindex[diag(L)[!estpars] > 0L]), ''), ' but should only be applied to
                 estimated parameters. Please fix!', call.=FALSE)
    }
    #make sure constrained pars are equal
    tmp <- rowSums(L)
    tmp[tmp == 0L] <- 1L
    check <- as.numeric(L %*% longpars) / tmp
    longpars[estpars] <- check[estpars]
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
    if(LRPARS){
        LBOUND <- c(LBOUND, lrPars@lbound)
        UBOUND <- c(UBOUND, lrPars@ubound)
    }
    for(g in 1L:ngroups)
        for(i in 1L:J)
            pars[[g]][[i]]@dat <- Data$fulldata[[g]][, c(itemloc[i]:(itemloc[i+1L] - 1L))]
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
        if(has_graded){
            for(g in 1L:length(pars)){
                pars[[g]][1L:J] <- lapply(pars[[g]][1L:J], function(x){
                    if(class(x) == 'graded'){
                        ds <- x@par[-(1L:x@nfact)]
                        x@par[-(1L:x@nfact)] <- sort(ds, decreasing = TRUE)
                        names(x@par) <- names(x@est)
                    }
                    return(x)
                })
            }
        }
        for(g in 1L:ngroups)
            gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
        if(LRPARS){
            lrPars@par <- lrPars@beta[] <- longpars[lrPars@parnum]
            lrPars@mus <- lrPars@X %*% lrPars@beta
            gstructgrouppars[[1L]]$gmeans <- lrPars@mus
        }
        if(RAND && cycles > RANDSTART) random <- reloadRandom(random=random, longpars=longpars,
                                        parstart=max(pars[[1L]][[J+1L]]@parnum) + 1L)

        start <- proc.time()[3L]
        if(RAND && cycles == RANDSTART){
            gtheta0[[1L]] <- matrix(0, nrow(gtheta0[[1L]]), ncol(gtheta0[[1L]]))
            OffTerm <- OffTerm(random, J=J, N=N)
            for(j in 1L:length(random)){
                tmp <- .1
                for(i in 1L:31L){
                    random[[j]]@drawvals <- DrawValues(random[[j]], itemloc=itemloc,
                                                       Theta=gtheta0[[1L]],
                                                       pars=pars[[1L]], fulldata=Data$fulldata[[1L]],
                                                       offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND)
                    OffTerm <- OffTerm(random, J=J, N=N)
                    if(is.null(list$cand.t.var)){
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
                    } else {
                        random[[j]]@cand.t.var <- list$cand.t.var[j + 1L]
                    }
                }
                #better start values
                tmp <- nrow(random[[j]]@drawvals)
                tmp <- cov(random[[j]]@drawvals) * (tmp / (tmp-1L))
                random[[j]]@par[random[[j]]@est] <- tmp[lower.tri(tmp, TRUE)][random[[j]]@est]
            }
            cand.t.var <- .5
            tmp <- .1
            for(i in 1L:31L){
                gtheta0[[1L]] <- draw.thetas(theta0=gtheta0[[1L]], pars=pars[[1L]], fulldata=Data$fulldata[[1L]],
                                             itemloc=itemloc, cand.t.var=cand.t.var, CUSTOM.IND=CUSTOM.IND,
                                             prior.t.var=gstructgrouppars[[1L]]$gcov, OffTerm=OffTerm,
                                             prior.mu=gstructgrouppars[[1L]]$gmeans, prodlist=prodlist)
                if(is.null(list$cand.t.var)){
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
                } else {
                    cand.t.var <- list$cand.t.var[1L]
                }
            }
            tmp <- nrow(gtheta0[[1L]])
            tmp <- cov(gtheta0[[1L]]) * (tmp / (tmp-1L))
            tmp2 <- c(rep(0, ncol(tmp)), tmp[lower.tri(tmp, TRUE)])
            pars[[1L]][[length(pars[[1L]])]]@par[pars[[1L]][[length(pars[[1L]])]]@est] <-
                tmp2[pars[[1L]][[length(pars[[1L]])]]@est]
        }

        # deal with latent regression random effects TODO
        if(LR.RAND){
            browser()
        }

        #Step 1. Generate m_k datasets of theta
        LL <- 0
        for(g in 1L:ngroups){
            for(i in 1L:MHDRAWS)
                gtheta0[[g]] <- draw.thetas(theta0=gtheta0[[g]], pars=pars[[g]], fulldata=Data$fulldata[[g]],
                                      itemloc=itemloc, cand.t.var=cand.t.var, CUSTOM.IND=CUSTOM.IND,
                                      prior.t.var=gstructgrouppars[[g]]$gcov, OffTerm=OffTerm,
                                      prior.mu=gstructgrouppars[[g]]$gmeans, prodlist=prodlist)
            LL <- LL + attr(gtheta0[[g]], "log.lik")
        }
        if(RAND && cycles > RANDSTART){
            for(j in 1L:length(random)){
                for(i in 1L:MHDRAWS){
                    random[[j]]@drawvals <- DrawValues(random[[j]], Theta=gtheta0[[1L]], itemloc=itemloc,
                                                       pars=pars[[1L]], fulldata=Data$fulldata[[1L]],
                                                       offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND)
                    OffTerm <- OffTerm(random, J=J, N=N)
                }
            }
        }
        #adjust cand.t.var
        PA <- sapply(gtheta0, function(x) attr(x, "Proportion Accepted"))
        NS <- sapply(gtheta0, function(x) nrow(x))
        cand.t.var <- controlCandVar(sum(PA * NS / sum(NS)), cand.t.var)
        if(RAND && cycles > (RANDSTART + 1L)){
            for(j in 1L:length(random)){
                random[[j]]@cand.t.var <- controlCandVar(
                               attr(random[[j]]@drawvals, "Proportion Accepted"),
                               random[[j]]@cand.t.var, min = .01, max = .5)

            }
        }
        Draws.time <- Draws.time + proc.time()[3L] - start

        #Step 2. Find average of simulated data gradients and hessian
        start <- proc.time()[3L]
        gthetatmp <- gtheta0
        if(length(prodlist))
            gthetatmp <- lapply(gtheta0, function(x, prodlist) prodterms(x, prodlist),
                              prodlist=prodlist)
        tmp <- .Call('computeDPars', pars, gthetatmp, OffTerm, length(longpars), TRUE,
                     USE.FIXED)
        g <- tmp$grad; h <- tmp$hess
        if(length(list$SLOW.IND)){
            for(group in 1L:ngroups){
                for (i in list$SLOW.IND){
                    deriv <- DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=gthetatmp[[group]],
                                                 estHess=TRUE)
                    g[pars[[group]][[i]]@parnum] <- deriv$grad
                    h[pars[[group]][[i]]@parnum, pars[[group]][[i]]@parnum] <- deriv$hess
                }
            }
        }
        for(group in 1L:ngroups){
            tmptheta <- gtheta0[[group]]
            if(is(gstructgrouppars[[1L]]$gmeans, 'matrix'))
                tmptheta <- tmptheta - gstructgrouppars[[1L]]$gmeans
            i <- J + 1L
            deriv <- Deriv(x=pars[[group]][[i]], Theta=tmptheta, CUSTOM.IND=CUSTOM.IND, estHess=TRUE)
            g[pars[[group]][[i]]@parnum] <- deriv$grad
            h[pars[[group]][[i]]@parnum, pars[[group]][[i]]@parnum] <- deriv$hess
        }
        if(RAND){
            if(cycles <= RANDSTART){
                for(i in 1L:length(random)){
                    g[random[[i]]@parnum] <- 0
                    h[random[[i]]@parnum, random[[i]]@parnum] <- diag(length(random[[i]]@parnum))
                }
            } else {
                for(i in 1L:length(random)){
                    deriv <- RandomDeriv(x=random[[i]])
                    g[random[[i]]@parnum] <- deriv$grad
                    h[random[[i]]@parnum, random[[i]]@parnum] <- deriv$hess
                }
            }
        }
        if(LRPARS){
            inv_sigma <- solve(gstructgrouppars[[1L]]$gcov)
            tmp <- t(inv_sigma %*% t(gtheta0[[1L]] - lrPars@mus) %*% lrPars@X)
            g[lrPars@parnum] <- as.numeric(tmp)
            tmp2 <- -det(inv_sigma) * lrPars@tXX
            for(i in 0L:(ncol(tmp)-1L))
                h[lrPars@parnum[1:nrow(tmp) + nrow(tmp)*i],
                  lrPars@parnum[1:nrow(tmp) + nrow(tmp)*i]] <- tmp2
        }
        if(length(constrain)){
            grad <- as.numeric(updateGrad(g, L))
            ave.h <- updateHess(-h, L)
        } else {
            grad <- g
            ave.h <- -h
        }
        grad <- grad[estpars & !redun_constr]
        ave.h <- ave.h[estpars & !redun_constr, estpars & !redun_constr]
        if(any(is.na(grad)))
            stop('Model did not converge (unacceptable gradient caused by extreme parameter values)',
                 call.=FALSE)
        if(is.na(attr(gtheta0[[1L]],"log.lik")))
            stop('Estimation halted. Model did not converge.', call.=FALSE)
        if(verbose){
            AR <- do.call(c, lapply(gtheta0, function(x) attr(x, "Proportion Accepted")))
            CTV <- cand.t.var
            if(RAND && cycles > RANDSTART){
                AR <- c(AR, do.call(c, lapply(random,
                                        function(x) attr(x@drawvals, "Proportion Accepted"))))
                CTV <- c(CTV, do.call(c, lapply(random,
                                                function(x) x@cand.t.var)))
            }
            AR <- paste0(sapply(AR, function(x) sprintf('%.2f', x)), collapse='; ')
            CTV <- paste0(sapply(CTV, function(x) sprintf('%.2f', x)), collapse='; ')
            if(cycles <= BURNIN)
                printmsg <- sprintf("\rStage 1 = %i, LL = %.1f, AR(%s) = [%s]",
                                    cycles, LL, CTV, AR)
            if(cycles > BURNIN && cycles <= BURNIN + SEMCYCLES)
                printmsg <- sprintf("\rStage 2 = %i, LL = %.1f, AR(%s) = [%s]",
                                    cycles-BURNIN, LL, CTV, AR)
            if(cycles > BURNIN + SEMCYCLES)
                printmsg <- sprintf("\rStage 3 = %i, LL = %.1f, AR(%s) = [%s]",
                                    cycles-BURNIN-SEMCYCLES, LL, CTV, AR)
        }
        if(stagecycle < 3L){
            correction <- try(solve(ave.h, grad), TRUE)
            if(is(correction, 'try-error')){
                ave.h.inv <- MPinv(ave.h)
                correction <- as.vector(grad %*% ave.h.inv)
            }
            correction[correction > 1] <- 1
            correction[correction < -1] <- -1
            longpars[estindex_unique] <- longpars[estindex_unique] + gamma*correction
            longpars[longpars < LBOUND] <- LBOUND[longpars < LBOUND]
            longpars[longpars > UBOUND] <- UBOUND[longpars > UBOUND]
            if(length(constrain))
                for(i in 1L:length(constrain))
                    longpars[index %in% constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
            if(verbose)
                cat(printmsg, sprintf(", Max-Change = %.4f", max(abs(gamma*correction))), sep='')
            if(stagecycle == 2L){
                SEM.stores[[cycles - BURNIN]] <- longpars
                SEM.stores2[[cycles - BURNIN]] <- ave.h
            }
            Mstep.time <- Mstep.time + proc.time()[3L] - start
            next
        }

        #Step 3. Update R-M step
        Tau <- Tau + gamma*(ave.h - Tau)
        correction <- try(solve(Tau, grad), TRUE)
        if(is(correction, 'try-error')){
            Tau.inv <- MPinv(Tau)
            correction <- as.vector(grad %*% Tau.inv)
        }
        correction[gamma*correction > .25] <- .25/gamma
        correction[gamma*correction < -.25] <- -.25/gamma
        longpars[estindex_unique] <- longpars[estindex_unique] + gamma*correction
        longpars[longpars < LBOUND] <- LBOUND[longpars < LBOUND]
        longpars[longpars > UBOUND] <- UBOUND[longpars > UBOUND]
        if(length(constrain))
            for(i in 1L:length(constrain))
                longpars[index %in% constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
        if(verbose)
            cat(printmsg, sprintf(", gam = %.4f, Max-Change = %.4f",
                                  gamma, max(abs(gamma*correction))), sep='')
        if(all(abs(gamma*correction) < TOL)) conv <- conv + 1L
        else conv <- 0L
        if(!list$SE && conv >= 3L) break
        if(list$SE && cycles >= (400L + BURNIN + SEMCYCLES) && conv >= 3L) break
        #Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE
        if(gamma == .25){
            gamma <- 0
            phi <- grad
            Phi <- ave.h
        }
        phi <- phi + gamma*(grad - phi)
        Phi <- Phi + gamma*(ave.h - outer(grad,grad) - Phi)
        Mstep.time <- Mstep.time + proc.time()[3L] - start
    } ###END BIG LOOP
    if(verbose) cat('\r\n')
    info <- Phi + outer(phi,phi)
    #Reload final pars list
    if(cycles == NCYCLES + BURNIN + SEMCYCLES && !list$USEEM){
        if(list$message)
            message('MHRM terminated after ', NCYCLES, ' iterations.')
        converge <- 0L
    }
    if(list$USEEM) longpars <- list$startlongpars
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L))
            pars[[g]][[i]]@par <- longpars[pars[[g]][[i]]@parnum]
    }
    fail_invert_info <- TRUE
    if(list$SE){
        acov <- try(solve(info), TRUE)
        if(is(acov, 'try-error')){
            if(list$warn)
                warning('Could not invert information matrix; model likely is not identified.',
                        call.=FALSE)
        } else {
            fail_invert_info <- FALSE
            SEtmp <- diag(acov)
            if(any(SEtmp < 0))
                SEtmp[SEtmp < 0 ] <- NaN
            SEtmp <- sqrt(SEtmp)
            SE <- rep(NA, length(longpars))
            SE[estindex_unique] <- SEtmp
            if(length(constrain) > 0L)
                for(i in 1L:length(constrain))
                    SE[index %in% constrain[[i]][-1L]] <- SE[constrain[[i]][1L]]
            for(g in 1L:ngroups){
                for(i in 1L:(J+1L))
                    pars[[g]][[i]]@SEpar <- SE[pars[[g]][[i]]@parnum]
            }
            if(RAND){
                for(i in 1L:length(random))
                    random[[i]]@SEpar <- SE[random[[i]]@parnum]
            }
            if(LRPARS){
                lrPars@SEpar <- SE[lrPars@parnum]
            }
        }
    }
    names(correction) <- names(estpars)[estindex_unique]
    info <- nameInfoMatrix(info=info, correction=correction, L=L, npars=ncol(L))
    if(!list$SE) info <- matrix(0, 1, 1)
    ret <- list(pars=pars, cycles = cycles - BURNIN - SEMCYCLES, info=if(list$expl) matrix(0) else info,
                correction=correction, longpars=longpars, converge=converge, SElogLik=0, cand.t.var=cand.t.var, L=L,
                random=random, lrPars=lrPars, time=c(MH_draws = as.numeric(Draws.time), Mstep=as.numeric(Mstep.time)),
                estindex_unique=estindex_unique, shortpars=longpars[estpars & !redun_constr],
                fail_invert_info=fail_invert_info, Prior=vector('list', ngroups), collectLL=NaN)
    ret
}
