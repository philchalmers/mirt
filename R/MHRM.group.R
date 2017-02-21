MHRM.group <- function(pars, constrain, Ls, Data, PrepList, list, random = list(),
                       lrPars = list(), lr.random = list(), DERIV, solnp_args, control)
{
    if(is.null(random)) random <- list()
    itemtype <- sapply(pars[[1]], class)
    has_graded <- any(itemtype == 'graded')
    RAND <- length(random) > 0L; LR.RAND <- length(lr.random) > 0L
    RANDSTART <- list$RANDSTART
    LRPARS <- length(lrPars) > 0L
    verbose <- list$verbose
    nfact <- list$nfact
    NCYCLES <- list$NCYCLES
    no_stage_3 <- FALSE
    if(is.na(NCYCLES)){
        NCYCLES <- 1L
        no_stage_3 <- TRUE
    }
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
    if(LR.RAND){
        for(i in 1L:length(lr.random)){
            nfullpars <- nfullpars + length(lr.random[[i]]@par)
            estpars <- c(estpars, lr.random[[i]]@est)
        }
    }
    if(all(!estpars) && list$SE)
        stop('Computing ACOV matrix is meaningless when no parameters are estimated', call.=FALSE)
    index <- 1L:nfullpars
    N <- nrow(gtheta0[[1L]])
    #Burn in
    if(LRPARS){
        lrPars@beta[] <- lrPars@par
        lrPars@mus <- lrPars@X %*% lrPars@beta
        gstructgrouppars[[1L]]$gmeans <- lrPars@mus
    }
    correction <- numeric(0L)
    cand.t.var <- if(is.null(list$cand.t.var)) 1 else list$cand.t.var[1L]
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
            }
        }
    }
    if(RAND) OffTerm <- OffTerm(random, J=J, N=N)
    if(list$plausible.draws > 0L){
        ret <- vector('list', list$plausible.draws)
        for(i in 1L:length(ret)){
            for(j in 1L:MHDRAWS)
                gtheta0[[1L]] <- draw.thetas(theta0=gtheta0[[1L]],
                        pars=pars[[1L]], fulldata=Data$fulldata[[1L]],
                        itemloc=itemloc, cand.t.var=cand.t.var, CUSTOM.IND=CUSTOM.IND,
                        prior.t.var=gstructgrouppars[[g]]$gcov, OffTerm=OffTerm,
                        prior.mu=gstructgrouppars[[g]]$gmeans, prodlist=prodlist)
            ret[[i]] <- gtheta0[[1L]]
        }
        ret <- lapply(ret, function(x){
            attr(x, "Proportion Accepted") <- attr(x, "log.lik") <- attr(x, "log.lik_full") <- NULL
            x
        })
        return(ret)
    }
    SEM.stores <- SEM.stores2 <- list()
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
        if(LR.RAND){
            for(i in 1L:length(lr.random))
                longpars[lr.random[[i]]@parnum] <- lr.random[[i]]@par
        }
    }
    names(longpars) <- names(estpars)
    stagecycle <- 1L
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
    if(LR.RAND){
        for(i in 1L:length(lr.random)){
            LBOUND <- c(LBOUND, lr.random[[i]]@lbound)
            UBOUND <- c(UBOUND, lr.random[[i]]@ubound)
        }
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
            gamma <- 0
        }
        #Reload pars list
        if(list$SE) longpars <- list$startlongpars
        tmp <- MHRM.reloadPars(longpars=longpars, pars=pars, gstructgrouppars=gstructgrouppars,
                               ngroups=ngroups, J=J, has_graded=has_graded, cycles=cycles,
                               LRPARS=LRPARS, LR.RAND=LR.RAND, RANDSTART=RANDSTART,
                               RAND=RAND, lrPars=lrPars, lr.random=lr.random, random=random)
        pars <- with(tmp, pars)
        gstructgrouppars <- with(tmp, gstructgrouppars)
        lr.random <- with(tmp, lr.random)
        random <- with(tmp, random)
        lrPars <- with(tmp, lrPars)
        if(cycles == (BURNIN + SEMCYCLES + 1L) && no_stage_3){
            cycles <- cycles + SEMCYCLES - 1L
            break
        }
        start <- proc.time()[3L]
        if((RAND || LR.RAND) && cycles == RANDSTART){
            gtheta0[[1L]] <- matrix(0, nrow(gtheta0[[1L]]), ncol(gtheta0[[1L]]))
            if(RAND){
                OffTerm <- OffTerm(random, J=J, N=N)
                if(!is.null(list$cand.t.var))
                    for(j in 1L:length(random)) random[[j]]@cand.t.var <- list$cand.t.var[j + 1L]
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
                        }
                    }
                    #better start values
                    tmp <- nrow(random[[j]]@drawvals)
                    tmp <- cov(random[[j]]@drawvals) * (tmp / (tmp-1L))
                    random[[j]]@par[random[[j]]@est] <- tmp[lower.tri(tmp, TRUE)][random[[j]]@est]
                }
            }
            if(LR.RAND){
                if(!is.null(list$cand.t.var)){
                    for(j in 1L:length(lr.random)) lr.random[[j]]@cand.t.var <-
                            list$cand.t.var[j + length(random) + 1L]
                }
                for(j in 1L:length(lr.random)){
                    tmp <- .1
                    for(i in 1L:31L){
                        lr.random[[j]]@drawvals <- DrawValues(lr.random[[j]], itemloc=itemloc,
                                                           Theta=gtheta0[[1L]], LR=TRUE,
                                                           pars=pars[[1L]], fulldata=Data$fulldata[[1L]],
                                                           offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND)
                        if(is.null(list$cand.t.var)){
                            if(i > 5L){
                                if(attr(lr.random[[j]]@drawvals,"Proportion Accepted") > .4)
                                    lr.random[[j]]@cand.t.var <- lr.random[[j]]@cand.t.var + 2*tmp
                                if(attr(lr.random[[j]]@drawvals,"Proportion Accepted") < .2)
                                    lr.random[[j]]@cand.t.var <- lr.random[[j]]@cand.t.var - 2*tmp
                                if(attr(lr.random[[j]]@drawvals,"Proportion Accepted") < .05)
                                    lr.random[[j]]@cand.t.var <- lr.random[[j]]@cand.t.var - 5*tmp
                                if (lr.random[[j]]@cand.t.var < 0){
                                    lr.random[[j]]@cand.t.var <- tmp
                                    tmp <- tmp / 10
                                }
                            }
                        }
                    }
                    #better start values
                    tmp <- nrow(lr.random[[j]]@drawvals)
                    tmp <- cov(lr.random[[j]]@drawvals) * (tmp / (tmp-1L))
                    lr.random[[j]]@par[lr.random[[j]]@est] <-
                        tmp[lower.tri(tmp, TRUE)][lr.random[[j]]@est]
                    for(j in 1L:length(lr.random))
                        gstructgrouppars[[1L]]$gmeans <- gstructgrouppars[[1L]]$gmeans +
                            lr.random[[j]]@drawvals[lr.random[[j]]@mtch]
                    tmp <- c(numeric(nfact),
                             as.vector(gstructgrouppars[[1L]]$gcov[lower.tri(gstructgrouppars[[1L]]$gcov, TRUE)]))
                    pars[[1L]][[J+1L]]@par[pars[[1L]][[J+1L]]@est] <- tmp[pars[[1L]][[J+1L]]@est]
                }
            }
            cand.t.var <- if(is.null(list$cand.t.var)) .5 else list$cand.t.var[1L]
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
                }
            }
            tmp <- nrow(gtheta0[[1L]])
            tmp <- cov(gtheta0[[1L]]) * (tmp / (tmp-1L))
            tmp2 <- c(rep(0, ncol(tmp)), tmp[lower.tri(tmp, TRUE)])
            pars[[1L]][[length(pars[[1L]])]]@par[pars[[1L]][[length(pars[[1L]])]]@est] <-
                tmp2[pars[[1L]][[length(pars[[1L]])]]@est]
            cand.t.var <- if(is.null(list$cand.t.var)) .5 else list$cand.t.var[1L]
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
        if(LR.RAND && cycles > RANDSTART){
            for(j in 1L:length(lr.random)){
                for(i in 1L:MHDRAWS){
                    lr.random[[j]]@drawvals <- DrawValues(lr.random[[j]], Theta=gtheta0[[1L]],
                                                          itemloc=itemloc, pars=pars[[1L]],
                                                          fulldata=Data$fulldata[[1L]],
                                                          offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND, LR=TRUE)
                }
            }
        }
        #adjust cand.t.var
        PA <- sapply(gtheta0, function(x) attr(x, "Proportion Accepted"))
        NS <- sapply(gtheta0, function(x) nrow(x))
        if(is.null(list$cand.t.var)){
            cand.t.var <- controlCandVar(sum(PA * NS / sum(NS)), cand.t.var)
            if(RAND && cycles > (RANDSTART + 1L)){
                for(j in 1L:length(random)){
                    random[[j]]@cand.t.var <- controlCandVar(
                                   attr(random[[j]]@drawvals, "Proportion Accepted"),
                                   random[[j]]@cand.t.var, min = .01, max = .5)
                }
            }
            if(LR.RAND && cycles > (RANDSTART + 1L)){
                for(j in 1L:length(lr.random)){
                    lr.random[[j]]@cand.t.var <- controlCandVar(
                        attr(lr.random[[j]]@drawvals, "Proportion Accepted"),
                        lr.random[[j]]@cand.t.var, min = .01, max = .5)
                }
            }
        }
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
            if(LR.RAND && cycles > RANDSTART){
                AR <- c(AR, do.call(c, lapply(lr.random,
                                              function(x) attr(x@drawvals, "Proportion Accepted"))))
                CTV <- c(CTV, do.call(c, lapply(lr.random,
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
        gthetatmp <- gtheta0
        if(length(prodlist))
            gthetatmp <- lapply(gtheta0, function(x, prodlist) prodterms(x, prodlist),
                                prodlist=prodlist)

        #Step 2. Find average of simulated data gradients and hessian
        Draws.time <- Draws.time + proc.time()[3L] - start
        start <- proc.time()[3L]
        tmp <- MHRM.deriv(pars=pars, gtheta=gthetatmp, lrPars=lrPars, OffTerm=OffTerm, longpars=longpars,
                          USE.FIXED=USE.FIXED, list=list, ngroups=ngroups, LR.RAND=LR.RAND,
                          DERIV=DERIV, gtheta0=gtheta0, gstructgrouppars=gstructgrouppars,
                          CUSTOM.IND=CUSTOM.IND, RAND=RAND, cycles=cycles, lr.random=lr.random,
                          RANDSTART=RANDSTART, random=random, J=J, LRPARS=LRPARS, L=L,
                          constrain=constrain, estpars=estpars, redun_constr=redun_constr)
        grad <- tmp$grad
        ave.h <- tmp$ave.h
        if(stagecycle < 3L){
            if(all(!estpars)) break
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
        if(list$SE.type == 'MHRM' && list$SE &&
           cycles >= (400L + BURNIN + SEMCYCLES) && conv >= 3L) break
        #Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE
        if(cycles == (BURNIN + SEMCYCLES + 1L)){
            if(list$SE.type == 'MHRM'){
                phi <- grad
                Phi <- ave.h
            } else Phi <- 0
        }
        if(list$SE.type != 'none' && list$SE){
            if(list$SE.type == 'MHRM'){
                phi <- phi + gamma*(grad - phi)
                Phi <- Phi + gamma*(ave.h - outer(grad,grad) - Phi)
            } else if(list$SE.type == 'FMHRM'){
                Phi <- Phi + 1/NCYCLES * (ave.h - outer(grad,grad))
            }
        }
        Mstep.time <- Mstep.time + proc.time()[3L] - start
    } ###END BIG LOOP
    if(verbose) cat('\r\n')
    #Reload final pars list
    if(cycles == NCYCLES + BURNIN + SEMCYCLES && !list$SE && !no_stage_3){
        if(list$message)
            message('MHRM terminated after ', NCYCLES, ' iterations.')
        converge <- FALSE
    }
    if(list$SE) longpars <- list$startlongpars
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L))
            pars[[g]][[i]]@par <- longpars[pars[[g]][[i]]@parnum]
    }
    fail_invert_info <- TRUE
    if(list$SE){
        if(list$SE.type == 'MHRM') info <- Phi + outer(phi,phi)
        else if(list$SE.type == 'FMHRM') info <- Phi
        acov <- try(solve(info), TRUE)
        if(is(acov, 'try-error')){
            if(list$warn)
                warning('Could not invert information matrix; may not be identified.',
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
            if(LR.RAND){
                for(i in 1L:length(lr.random))
                    lr.random[[i]]@SEpar <- SE[lr.random[[i]]@parnum]
            }
        }
    }
    if(any(estpars)){
        names(correction) <- names(estpars)[estindex_unique]
        if(list$SE) info <- nameInfoMatrix(info=info, correction=correction, L=L, npars=ncol(L))
        else info <- matrix(0, 1L, 1L)
    } else {
        info <- matrix(0, 1L, 1L)
        cycles <- BURNIN + SEMCYCLES
    }
    ret <- list(pars=pars, cycles = cycles - BURNIN - SEMCYCLES, info=if(list$expl) matrix(0) else info,
                correction=correction, longpars=longpars, converge=converge, SElogLik=0, cand.t.var=cand.t.var, L=L,
                random=random, lrPars=lrPars, lr.random=lr.random,
                time=c(MH_draws = as.numeric(Draws.time), Mstep=as.numeric(Mstep.time)),
                estindex_unique=estindex_unique, shortpars=longpars[estpars & !redun_constr],
                fail_invert_info=fail_invert_info, Prior=vector('list', ngroups), collectLL=NaN)
    ret
}

MHRM.deriv <- function(pars, gtheta, OffTerm, longpars, USE.FIXED, list, ngroups,
                      DERIV, gtheta0, gstructgrouppars, CUSTOM.IND, RAND,
                      cycles, RANDSTART, random, J, LRPARS, lrPars, LR.RAND, lr.random,
                      constrain, estpars, redun_constr, L, estHess = TRUE){
    tmp <- .Call('computeDPars', pars, gtheta, OffTerm, length(longpars), estHess,
                 USE.FIXED, 0L, FALSE)
    g <- tmp$grad
    h <- tmp$hess
    if(length(list$SLOW.IND)){
        for(group in 1L:ngroups){
            for (i in list$SLOW.IND){
                deriv <- DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=gtheta[[group]],
                                             estHess=estHess)
                g[pars[[group]][[i]]@parnum] <- deriv$grad
                if(estHess)
                    h[pars[[group]][[i]]@parnum, pars[[group]][[i]]@parnum] <- deriv$hess
            }
        }
    }
    for(group in 1L:ngroups){
        tmptheta <- gtheta0[[group]]
        if(is(gstructgrouppars[[1L]]$gmeans, 'matrix'))
            tmptheta <- tmptheta - gstructgrouppars[[1L]]$gmeans
        i <- J + 1L
        deriv <- Deriv(x=pars[[group]][[i]], Theta=tmptheta,
                       CUSTOM.IND=CUSTOM.IND, estHess=estHess)
        g[pars[[group]][[i]]@parnum] <- deriv$grad
        if(estHess)
            h[pars[[group]][[i]]@parnum, pars[[group]][[i]]@parnum] <- deriv$hess
    }
    if(RAND){
        if(cycles <= RANDSTART){
            for(i in 1L:length(random)){
                g[random[[i]]@parnum] <- 0
                h[random[[i]]@parnum, random[[i]]@parnum] <- -diag(length(random[[i]]@parnum))
            }
        } else {
            for(i in 1L:length(random)){
                deriv <- RandomDeriv(x=random[[i]], estHess=estHess)
                g[random[[i]]@parnum] <- deriv$grad
                if(estHess)
                    h[random[[i]]@parnum, random[[i]]@parnum] <- deriv$hess
            }
        }
    }
    if(LRPARS){
        deriv <- Deriv(lrPars, cov=gstructgrouppars[[1L]]$gcov, theta=gtheta0[[1L]],
                       estHess=estHess)
        g[lrPars@parnum] <- deriv$grad
        if(estHess)
            for(i in 0L:(ncol(deriv$grad)-1L))
                h[lrPars@parnum[1L:nrow(deriv$grad) + nrow(deriv$grad)*i],
                  lrPars@parnum[1L:nrow(deriv$grad) + nrow(deriv$grad)*i]] <- deriv$hess

    }
    if(LR.RAND){
        if(cycles <= RANDSTART){
            for(i in 1L:length(lr.random)){
                g[lr.random[[i]]@parnum] <- 0
                h[lr.random[[i]]@parnum, lr.random[[i]]@parnum] <- -diag(length(lr.random[[i]]@parnum))
            }
        } else {
            for(i in 1L:length(lr.random)){
                deriv <- RandomDeriv(x=lr.random[[i]], estHess=estHess)
                g[lr.random[[i]]@parnum] <- deriv$grad
                if(estHess)
                    h[lr.random[[i]]@parnum, lr.random[[i]]@parnum] <- deriv$hess
            }
        }
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
    if(any(is.na(grad))){
        stop('Model did not converge (unacceptable gradient caused by extreme parameter values)',
             call.=FALSE)
    }
    list(grad=grad, ave.h=ave.h)
}

MHRM.reloadPars <- function(longpars, pars, gstructgrouppars, ngroups, J, has_graded,
                            cycles, LRPARS, LR.RAND, RANDSTART, RAND, lrPars, lr.random, random){
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
    if(LR.RAND && cycles > RANDSTART){
        for(j in 1L:length(lr.random))
            gstructgrouppars[[1L]]$gmeans <- gstructgrouppars[[1L]]$gmeans +
                lr.random[[j]]@drawvals[lr.random[[j]]@mtch]
    }
    if(RAND && cycles > RANDSTART) random <- reloadRandom(random=random, longpars=longpars)
    if(LR.RAND && cycles > RANDSTART) lr.random <- reloadRandom(random=lr.random, longpars=longpars)
    list(pars=pars, gstructgrouppars=gstructgrouppars, lr.random=lr.random, random=random,
         lrPars=lrPars)
}
