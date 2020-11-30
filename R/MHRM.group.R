MHRM.group <- function(pars, constrain, Ls, Data, PrepList, list, random = list(),
                       lrPars = list(), lr.random = list(), DERIV, solnp_args, control)
{
    if(is.null(random)) random <- list()
    correction <- 0
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
    for(g in seq_len(ngroups)){
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
        gtheta0[[g]] <- matrix(0, nrow(Data$fulldata[[g]]), nfact)
        for(i in seq_len(J+1L)){
            nfullpars <- nfullpars + length(pars[[g]][[i]]@par)
            estpars <- c(estpars, pars[[g]][[i]]@est)
        }
    }
    if(RAND){
        for(i in seq_len(length(random))){
            nfullpars <- nfullpars + length(random[[i]]@par)
            estpars <- c(estpars, random[[i]]@est)
        }
    }
    if(LRPARS){
        nfullpars <- nfullpars + length(lrPars@par)
        estpars <- c(estpars, lrPars@est)
    }
    if(LR.RAND){
        for(i in seq_len(length(lr.random))){
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
    cand.t.var <- if(is.null(list$cand.t.var)) 1 else list$cand.t.var[1L]
    tmp <- .1
    OffTerm <- matrix(0, 1, J)
    for(g in seq_len(ngroups)){
        PAs <- CTVs <- rep(NA, 25L)
        for(i in seq_len(31L)){
            gtheta0[[g]] <- draw.thetas(theta0=gtheta0[[g]], pars=pars[[g]], fulldata=Data$fulldata[[g]],
                                        itemloc=itemloc, cand.t.var=cand.t.var, CUSTOM.IND=CUSTOM.IND,
                                        prior.t.var=gstructgrouppars[[g]]$gcov, OffTerm=OffTerm,
                                        prior.mu=gstructgrouppars[[g]]$gmeans, prodlist=prodlist)
            if(is.null(list$cand.t.var)){
                if(i > 5L){
                    pa <- attr(gtheta0[[g]],"Proportion Accepted")
                    PAs[i-5L] <- pa
                    CTVs[i-5L] <- cand.t.var
                    cand.t.var <- update_cand.var(PAs, CTVs)
                }
            }
        }
    }
    if(RAND) OffTerm <- OffTerm(random, J=J, N=N)
    if(list$plausible.draws > 0L){
        ret <- vector('list', list$plausible.draws)
        for(i in seq_len(length(ret))){
            for(j in seq_len(MHDRAWS))
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
    SEM.stores <- matrix(0, SEMCYCLES, nfullpars)
    SEM.stores2 <- vector('list', SEMCYCLES)
    conv <- 0L
    k <- 1L
    gamma <- .2
    longpars <- rep(NA,nfullpars)
    for(g in seq_len(ngroups)){
        for(i in seq_len(J+1L))
            longpars[pars[[g]][[i]]@parnum] <- pars[[g]][[i]]@par
        if(RAND){
            for(i in seq_len(length(random)))
                longpars[random[[i]]@parnum] <- random[[i]]@par
        }
        if(LRPARS)
            longpars[lrPars@parnum] <- lrPars@par
        if(LR.RAND){
            for(i in seq_len(length(lr.random)))
                longpars[lr.random[[i]]@parnum] <- lr.random[[i]]@par
        }
    }
    names(longpars) <- names(estpars)
    stagecycle <- 1L
    converge <- TRUE
    L <- Ls$L
    redun_constr <- Ls$redun_constr
    estindex_unique <- index[estpars & !redun_constr]
    if(any(attr(L, 'diag')[!estpars] > 0L)){
        redindex <- index[!estpars]
        stop('Constraint applied to fixed parameter(s) ',
             paste(paste0(redindex[attr(L, 'diag')[!estpars] > 0L]), ''), ' but should only be applied to
                 estimated parameters. Please fix!', call.=FALSE)
    }
    #make sure constrained pars are equal
    tmp <- Matrix::rowSums(L)
    tmp[tmp == 0L] <- 1L
    check <- as.numeric(L %*% longpars) / tmp
    longpars[estpars] <- check[estpars]
    LBOUND <- UBOUND <- c()
    for(g in seq_len(ngroups)){
        for(i in seq_len(J+1L)){
            LBOUND <- c(LBOUND, pars[[g]][[i]]@lbound)
            UBOUND <- c(UBOUND, pars[[g]][[i]]@ubound)
        }
    }
    if(RAND){
        for(i in seq_len(length(random))){
            LBOUND <- c(LBOUND, random[[i]]@lbound)
            UBOUND <- c(UBOUND, random[[i]]@ubound)
        }
    }
    if(LRPARS){
        LBOUND <- c(LBOUND, lrPars@lbound)
        UBOUND <- c(UBOUND, lrPars@ubound)
    }
    if(LR.RAND){
        for(i in seq_len(length(lr.random))){
            LBOUND <- c(LBOUND, lr.random[[i]]@lbound)
            UBOUND <- c(UBOUND, lr.random[[i]]@ubound)
        }
    }
    aveAR <- vector('list', NCYCLES)
    control$fnscale <- -1
    if(is.null(control$maxit)) control$maxit <- 15
    if(list$Moptim == 'BFGS')
        if(any(is.finite(LBOUND[estindex_unique]) | is.finite(UBOUND[estindex_unique])))
            list$Moptim <- 'L-BFGS-B'
    for(g in seq_len(ngroups))
        for(i in seq_len(J))
            pars[[g]][[i]]@dat <- Data$fulldata[[g]][, c(itemloc[i]:(itemloc[i+1L] - 1L))]
    Draws.time <- Mstep.time <- 0

    ####Big MHRM loop
    for(cycles in seq_len(NCYCLES + BURNIN + SEMCYCLES))
    {
        if(cycles == BURNIN + 1L) stagecycle <- 2L
        if(stagecycle == 3L)
            gamma <- gain_fun(gain, t = cycles - SEMCYCLES - BURNIN - 1L)
        if(cycles == (BURNIN + SEMCYCLES + 1L)){
            stagecycle <- 3L
            longpars <- colMeans(SEM.stores)
            if(!no_stage_3){
                Tau <- SEM.stores2[[1L]]/SEMCYCLES
                if(SEMCYCLES > 1L){
                    for(i in 2L:SEMCYCLES)
                        Tau <- Tau + SEM.stores2[[i]]/SEMCYCLES
                }
            }
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
        tmp <- MHRM.draws(pars=pars, lrPars=lrPars, lr.random=lr.random, random=random,
                          gstructgrouppars=gstructgrouppars, RAND=RAND, LR.RAND=LR.RAND,
                          RANDSTART=RANDSTART, gtheta0=gtheta0, OffTerm=OffTerm, J=J, N=N, cycles=cycles,
                          itemloc=itemloc, CUSTOM.IND=CUSTOM.IND, Data=Data, nfact=nfact,
                          prodlist=prodlist, ngroups=ngroups, MHDRAWS=MHDRAWS, BURNIN=BURNIN,
                          SEMCYCLES=SEMCYCLES, cand.t.var=cand.t.var, list=list, verbose=verbose)
        lr.random <- with(tmp, lr.random)
        random <- with(tmp, random)
        gtheta0 <- with(tmp, gtheta0)
        OffTerm <- with(tmp, OffTerm)
        printmsg <- with(tmp, printmsg)

        if(stagecycle == 3L)
            aveAR[[cycles - SEMCYCLES - BURNIN]] <- tmp$AR
        cand.t.var <- with(tmp, cand.t.var)
        gthetatmp <- gtheta0
        if(length(prodlist))
            gthetatmp <- lapply(gtheta0, function(x, prodlist) prodterms(x, prodlist),
                                prodlist=prodlist)
        if(all(!estpars)) break

        #Step 2. Find average of simulated data gradients and hessian
        Draws.time <- Draws.time + proc.time()[3L] - start
        start <- proc.time()[3L]
        longpars0 <- longpars
        tmp <- MHRM.Mstep(pars=pars, gtheta=gthetatmp, lrPars=lrPars, OffTerm=OffTerm, longpars=longpars,
                          USE.FIXED=USE.FIXED, list=list, ngroups=ngroups, LR.RAND=LR.RAND,
                          DERIV=DERIV, gstructgrouppars=gstructgrouppars, nfact=nfact,
                          CUSTOM.IND=CUSTOM.IND, RAND=RAND, cycles=cycles, lr.random=lr.random,
                          RANDSTART=RANDSTART, random=random, J=J, LRPARS=LRPARS, L=L, has_graded=has_graded,
                          constrain=constrain, estpars=estpars, redun_constr=redun_constr,
                          estindex_unique=estindex_unique,LBOUND=LBOUND, UBOUND=UBOUND,
                          gfulldata=Data$fulldata, itemloc=itemloc, control=control)
        grad <- tmp$grad
        ave.h <- tmp$hess
        correction <- tmp$correction
        if(stagecycle < 3L){
            correction[correction > 1] <- 1
            correction[correction < -1] <- -1
            longpars[estindex_unique] <- longpars[estindex_unique] + gamma*correction
            if(any(longpars < LBOUND))
                longpars[longpars < LBOUND] <- (longpars0[longpars < LBOUND] + LBOUND[longpars < LBOUND])/2
            if(any(longpars > UBOUND))
                longpars[longpars > UBOUND] <- (longpars0[longpars > UBOUND] + UBOUND[longpars > UBOUND])/2
            if(length(constrain))
                for(i in seq_len(length(constrain)))
                    longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
            if(verbose)
                cat(printmsg, sprintf(", Max-Change = %.4f", max(abs(gamma*correction))), sep='')
            if(stagecycle == 2L){
                SEM.stores[cycles - BURNIN, ] <- longpars
                if(!no_stage_3)
                    SEM.stores2[[cycles - BURNIN]] <- ave.h
            }
            Mstep.time <- Mstep.time + proc.time()[3L] - start
            next
        }

        #Step 3. Update R-M step
        Tau <- Tau + gamma*(ave.h - Tau)
        if(list$Moptim == 'NR1'){
            correction <- try(solve(Tau, grad), TRUE)
            if(is(correction, 'try-error')){
                Tau.inv <- MPinv(Tau)
                correction <- as.vector(grad %*% Tau.inv)
            }
            longpars[estindex_unique] <- longpars[estindex_unique] + gamma*correction
            if(any(longpars < LBOUND))
                longpars[longpars < LBOUND] <- (longpars0[longpars < LBOUND] + LBOUND[longpars < LBOUND])/2
            if(any(longpars > UBOUND))
                longpars[longpars > UBOUND] <- (longpars0[longpars > UBOUND] + UBOUND[longpars > UBOUND])/2
            if(length(constrain))
                for(i in seq_len(length(constrain)))
                    longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
        }
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
            } else phi <- Phi <- 0
        }
        if(list$SE.type != 'none' && list$SE){
            if(list$SE.type == 'MHRM'){
                phi <- phi + gamma*(grad - phi)
                Phi <- Phi + gamma*(ave.h - outer(grad,grad) - Phi)
            } else if(list$SE.type == 'FMHRM'){
                phi <- phi + 1/NCYCLES * grad
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
    tmp <- MHRM.reloadPars(longpars=longpars, pars=pars, gstructgrouppars=gstructgrouppars,
                           ngroups=ngroups, J=J, has_graded=has_graded, cycles=cycles,
                           LRPARS=LRPARS, LR.RAND=LR.RAND, RANDSTART=RANDSTART,
                           RAND=RAND, lrPars=lrPars, lr.random=lr.random, random=random)
    pars <- with(tmp, pars)
    gstructgrouppars <- with(tmp, gstructgrouppars)
    lr.random <- with(tmp, lr.random)
    random <- with(tmp, random)
    lrPars <- with(tmp, lrPars)
    fail_invert_info <- TRUE
    if(list$SE){
        info <- Phi + outer(phi,phi)
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
                for(i in seq_len(length(constrain)))
                    SE[constrain[[i]][-1L]] <- SE[constrain[[i]][1L]]
            for(g in seq_len(ngroups)){
                for(i in seq_len(J+1L))
                    pars[[g]][[i]]@SEpar <- SE[pars[[g]][[i]]@parnum]
            }
            if(RAND){
                for(i in seq_len(length(random)))
                    random[[i]]@SEpar <- SE[random[[i]]@parnum]
            }
            if(LRPARS){
                lrPars@SEpar <- SE[lrPars@parnum]
            }
            if(LR.RAND){
                for(i in seq_len(length(lr.random)))
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
    ret <- list(pars=pars, cycles = cycles-BURNIN-SEMCYCLES, info=if(list$expl) matrix(0) else as.matrix(info),
                correction=correction, longpars=longpars, converge=converge, SElogLik=0, cand.t.var=cand.t.var,
                L=L, random=random, lrPars=lrPars, lr.random=lr.random, aveAR=colMeans(do.call(rbind, aveAR)),
                time=c(MH_draws = as.numeric(Draws.time), Mstep=as.numeric(Mstep.time)),
                estindex_unique=estindex_unique, shortpars=longpars[estpars & !redun_constr],
                fail_invert_info=fail_invert_info, Prior=vector('list', ngroups), collectLL=NaN)
    ret
}
