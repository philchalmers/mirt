MHRM.deriv <- function(pars, gtheta, OffTerm, longpars, USE.FIXED, list, ngroups,
                      DERIV, gstructgrouppars, CUSTOM.IND, RAND, nfact,
                      cycles, RANDSTART, random, J, LRPARS, lrPars, LR.RAND, lr.random,
                      constrain, estpars, redun_constr, L, estHess = TRUE){
    tmp <- .Call('computeDPars', pars, gtheta, OffTerm, length(longpars), estHess,
                 USE.FIXED, 0L, FALSE)
    g <- tmp$grad
    h <- tmp$hess
    if(length(list$SLOW.IND)){
        for(group in seq_len(ngroups)){
            for (i in list$SLOW.IND){
                deriv <- DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=gtheta[[group]],
                                             estHess=estHess)
                g[pars[[group]][[i]]@parnum] <- deriv$grad
                if(estHess)
                    h[pars[[group]][[i]]@parnum, pars[[group]][[i]]@parnum] <- deriv$hess
            }
        }
    }
    for(group in seq_len(ngroups)){
        tmptheta <- gtheta[[group]][,1L:nfact, drop=FALSE]
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
            for(i in seq_len(length(random))){
                g[random[[i]]@parnum] <- 0
                h[random[[i]]@parnum, random[[i]]@parnum] <- -diag(length(random[[i]]@parnum))
            }
            if(cycles == (RANDSTART - 1L)) #hack for R 3.4.0
                deriv2 <- RandomDeriv(x=random[[i]], estHess=estHess)
        } else {
            for(i in seq_len(length(random))){
                deriv2 <- RandomDeriv(x=random[[i]], estHess=estHess)
                g[random[[i]]@parnum] <- deriv2$grad
                if(estHess)
                    h[random[[i]]@parnum, random[[i]]@parnum] <- deriv2$hess
            }
        }
    }
    if(LRPARS){
        deriv3 <- Deriv(lrPars, cov=gstructgrouppars[[1L]]$gcov,
                       theta=gtheta[[1L]][,1L:nfact, drop=FALSE], estHess=estHess)
        g[lrPars@parnum] <- deriv3$grad
        if(estHess)
            for(i in 0L:(ncol(deriv3$grad)-1L))
                h[lrPars@parnum[1L:nrow(deriv3$grad) + nrow(deriv3$grad)*i],
                  lrPars@parnum[1L:nrow(deriv3$grad) + nrow(deriv3$grad)*i]] <- deriv3$hess

    }
    if(LR.RAND){
        if(cycles <= RANDSTART){
            for(i in seq_len(length(lr.random))){
                g[lr.random[[i]]@parnum] <- 0
                h[lr.random[[i]]@parnum, lr.random[[i]]@parnum] <- -diag(length(lr.random[[i]]@parnum))
            }
            if(cycles == (RANDSTART - 1L)) #hack for R 3.4.0
                deriv4 <- RandomDeriv(x=lr.random[[i]], estHess=estHess)
        } else {
            for(i in seq_len(length(lr.random))){
                deriv4 <- RandomDeriv(x=lr.random[[i]], estHess=estHess)
                g[lr.random[[i]]@parnum] <- deriv4$grad
                if(estHess)
                    h[lr.random[[i]]@parnum, lr.random[[i]]@parnum] <- deriv4$hess
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
    if(any(is.na(grad) | is.infinite(grad))){
        pick <- which(is.na(grad) | is.infinite(grad))
        grad[pick] <- sample(seq(-.05, .05, length.out=100L), length(pick))
        ave.h[pick, pick] <- -diag(length(pick))
    }
    list(grad=grad, ave.h=ave.h)
}

MHRM.deriv_reload <- function(shortpars, longpars, pars, gtheta, lrPars, OffTerm, USE.FIXED,
                             list, ngroups, LR.RAND, DERIV, has_graded, nfact,
                             CUSTOM.IND, RAND, cycles, lr.random, RANDSTART, random, J, LRPARS, L,
                             constrain, estpars, redun_constr, estindex_unique, gfulldata, itemloc,
                             estHess = FALSE){
    longpars[estindex_unique] <- shortpars
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    tmp <- MHRM.reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J,
                           has_graded=has_graded, cycles=cycles, LRPARS=LRPARS,
                           LR.RAND=LR.RAND, RANDSTART=RANDSTART, gstructgrouppars=vector('list', ngroups),
                           RAND=RAND, lrPars=lrPars, lr.random=lr.random, random=random)
    pars <- with(tmp, pars)
    gstructgrouppars <- with(tmp, gstructgrouppars)
    lr.random <- with(tmp, lr.random)
    random <- with(tmp, random)
    lrPars <- with(tmp, lrPars)
    ret <- MHRM.deriv(pars=pars, gtheta=gtheta, lrPars=lrPars, OffTerm=OffTerm, longpars=longpars,
                      USE.FIXED=USE.FIXED, list=list, ngroups=ngroups, LR.RAND=LR.RAND,
                      DERIV=DERIV, gstructgrouppars=gstructgrouppars, nfact=nfact,
                      CUSTOM.IND=CUSTOM.IND, RAND=RAND, cycles=cycles, lr.random=lr.random,
                      RANDSTART=RANDSTART, random=random, J=J, LRPARS=LRPARS, L=L,
                      constrain=constrain, estpars=estpars, redun_constr=redun_constr,
                      estHess = FALSE)
    if(estHess) return(ret)
    else return(ret$grad)
}

MHRM.LL <- function(pars, gstructgrouppars, gtheta, lr.random, random, lrPars, OffTerm,
                    ngroups, nfact, J, USE.FIXED, LR.RAND, RAND, RANDSTART, LRPARS, CUSTOM.IND,
                    gfulldata, itemloc){
    LL <- 0
    for(g in seq_len(ngroups))
        LL <- LL + sum(complete.LL(theta=gtheta[[g]], pars=pars[[g]], nfact=nfact,
                               prior.mu=gstructgrouppars[[g]]$gmeans,
                               prior.t.var=gstructgrouppars[[g]]$gcov,
                               OffTerm=OffTerm, CUSTOM.IND=CUSTOM.IND,
                               itemloc=itemloc, fulldata=gfulldata[[g]]))
    LL
}

MHRM.LL_reload <- function(shortpars, longpars, pars, gtheta, lrPars, OffTerm, USE.FIXED,
                           list, ngroups, LR.RAND, DERIV, has_graded, nfact,
                           CUSTOM.IND, RAND, cycles, lr.random, RANDSTART, random, J, LRPARS, L,
                           constrain, estpars, redun_constr, estindex_unique, gfulldata, itemloc){
    longpars[estindex_unique] <- shortpars
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    tmp <- MHRM.reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J,
                           has_graded=has_graded, cycles=cycles, LRPARS=LRPARS,
                           LR.RAND=LR.RAND, RANDSTART=RANDSTART, gstructgrouppars=vector('list', ngroups),
                           RAND=RAND, lrPars=lrPars, lr.random=lr.random, random=random)
    pars <- with(tmp, pars)
    gstructgrouppars <- with(tmp, gstructgrouppars)
    lr.random <- with(tmp, lr.random)
    random <- with(tmp, random)
    lrPars <- with(tmp, lrPars)
    LL <- MHRM.LL(pars=pars, gstructgrouppars=gstructgrouppars, gtheta=gtheta, lr.random=lr.random,
                  random=random, lrPars=lrPars, ngroups=ngroups, nfact=nfact, J=J, OffTerm=OffTerm,
                  USE.FIXED=USE.FIXED, LR.RAND=LR.RAND, RAND=RAND, RANDSTART=RANDSTART, LRPARS=LRPARS,
                  CUSTOM.IND=CUSTOM.IND, gfulldata=gfulldata, itemloc=itemloc)
    LL
}

MHRM.NR <- function(shortpars, pars, gtheta, lrPars, OffTerm, longpars, USE.FIXED, list, ngroups,
                    LR.RAND, DERIV, has_graded, nfact, CUSTOM.IND, RAND, cycles, lr.random,
                    RANDSTART, random, J, LRPARS, L, constrain, estpars, redun_constr,
                    estindex_unique, gfulldata, itemloc, lbound, ubound, control){
    plast2 <- plast <- p <- shortpars
    lastchange <- 0
    gstructgrouppars <- vector('list', ngroups)
    if(is.null(control$maxit)) control$maxit <- 50L
    for(iter in seq_len(control$maxit)){
        longpars[estindex_unique] <- shortpars
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        for(g in seq_len(ngroups))
            gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
        tmp <- MHRM.deriv(pars=pars, gtheta=gtheta, lrPars=lrPars, OffTerm=OffTerm, longpars=longpars,
                         USE.FIXED=USE.FIXED, list=list, ngroups=ngroups, LR.RAND=LR.RAND,
                         DERIV=DERIV, nfact=nfact, CUSTOM.IND=CUSTOM.IND, RAND=RAND, cycles=cycles,
                         lr.random=lr.random, RANDSTART=RANDSTART, random=random, J=J, LRPARS=LRPARS, L=L,
                         constrain=constrain, estpars=estpars, redun_constr=redun_constr,
                         gstructgrouppars=gstructgrouppars, estHess=TRUE)
        h <- tmp$ave.h
        g <- tmp$grad
        inv <- MPinv(h) #TODO this could be avoided if no constrains present
        change <- as.vector(g %*% inv)
        change <- ifelse(change > .25, .25, change)
        change <- ifelse(change < -.25, -.25, change)
        plast2 <- plast
        plast <- p
        p <- p + change
        if(iter > 1L){
            flip <- (sign(lastchange) * sign(change)) == -1L
            p[flip] <- (plast[flip] + p[flip]) / 2
        }
        if(any(p < lbound))
            p[p < lbound] <- (plast[p < lbound] + lbound[p < lbound])/2
        if(any(p > ubound))
            p[p > ubound] <- (plast[p > ubound] + ubound[p > ubound])/2
        dif <- plast - p
        if(all(abs(dif) < control$tol)) break
        lastchange <- change
    }
    p - shortpars
}

MHRM.reloadPars <- function(longpars, pars, gstructgrouppars, ngroups, J, has_graded,
                            cycles, LRPARS, LR.RAND, RANDSTART, RAND, lrPars, lr.random, random){
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    if(has_graded){
        for(g in seq_len(length(pars))){
            pars[[g]][seq_len(J)] <- lapply(pars[[g]][seq_len(J)], function(x){
                if(class(x) == 'graded'){
                    ds <- x@par[-seq_len(x@nfact)]
                    x@par[-seq_len(x@nfact)] <- sort(ds, decreasing = TRUE)
                    names(x@par) <- x@parnames
                }
                return(x)
            })
        }
    }
    for(g in seq_len(ngroups))
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
    if(LRPARS){
        lrPars@par <- lrPars@beta[] <- longpars[lrPars@parnum]
        lrPars@mus <- lrPars@X %*% lrPars@beta
        gstructgrouppars[[1L]]$gmeans <- lrPars@mus
    }
    if(LR.RAND && cycles > RANDSTART){
        for(j in seq_len(length(lr.random)))
            gstructgrouppars[[1L]]$gmeans <- gstructgrouppars[[1L]]$gmeans +
                lr.random[[j]]@drawvals[lr.random[[j]]@mtch]
    }
    if(RAND && cycles > RANDSTART) random <- reloadRandom(random=random, longpars=longpars)
    if(LR.RAND && cycles > RANDSTART) lr.random <- reloadRandom(random=lr.random, longpars=longpars)
    list(pars=pars, gstructgrouppars=gstructgrouppars, lr.random=lr.random, random=random,
         lrPars=lrPars)
}

MHRM.draws <- function(pars, lrPars, lr.random, random, gstructgrouppars, OffTerm, RAND, LR.RAND, RANDSTART,
                       gtheta0, J, N, cycles, itemloc, CUSTOM.IND, Data, nfact, prodlist, ngroups,
                       MHDRAWS, BURNIN, SEMCYCLES, cand.t.var, list, verbose){
    if((RAND || LR.RAND) && cycles == RANDSTART){
        gtheta0[[1L]] <- matrix(0, nrow(gtheta0[[1L]]), ncol(gtheta0[[1L]]))
        if(RAND){
            OffTerm <- OffTerm(random, J=J, N=N)
            if(!is.null(list$cand.t.var))
                for(j in seq_len(length(random))) random[[j]]@cand.t.var <- list$cand.t.var[j + 1L]
                for(j in seq_len(length(random))){
                    tmp <- .1
                    for(i in seq_len(31L)){
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
                for(j in seq_len(length(lr.random))) lr.random[[j]]@cand.t.var <-
                        list$cand.t.var[j + length(random) + 1L]
            }
            for(j in seq_len(length(lr.random))){
                tmp <- .1
                for(i in seq_len(31L)){
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
                for(j in seq_len(length(lr.random)))
                    gstructgrouppars[[1L]]$gmeans <- gstructgrouppars[[1L]]$gmeans +
                    lr.random[[j]]@drawvals[lr.random[[j]]@mtch]
                tmp <- c(numeric(nfact),
                         as.vector(gstructgrouppars[[1L]]$gcov[lower.tri(gstructgrouppars[[1L]]$gcov, TRUE)]))
                pars[[1L]][[J+1L]]@par[pars[[1L]][[J+1L]]@est] <- tmp[pars[[1L]][[J+1L]]@est]
            }
        }
        cand.t.var <- if(is.null(list$cand.t.var)) .5 else list$cand.t.var[1L]
        tmp <- .1
        for(i in seq_len(31L)){
            gtheta0[[1L]] <- draw.thetas(theta0=gtheta0[[1L]], pars=pars[[1L]], fulldata=Data$fulldata[[1L]],
                                         itemloc=itemloc, cand.t.var=cand.t.var, CUSTOM.IND=CUSTOM.IND,
                                         prior.t.var=gstructgrouppars[[1L]]$gcov, OffTerm=OffTerm,
                                         prior.mu=gstructgrouppars[[1L]]$gmeans, prodlist=prodlist)
            if(is.null(list$cand.t.var)){
                if(i > 5L){
                    if(attr(gtheta0[[1L]],"Proportion Accepted") > .35) cand.t.var <- cand.t.var + 2*tmp
                    else if(attr(gtheta0[[1L]],"Proportion Accepted") > .25 && nfact > 3L)
                        cand.t.var <- cand.t.var + tmp
                    else if(attr(gtheta0[[1L]],"Proportion Accepted") < .2 && nfact < 4L)
                        cand.t.var <- cand.t.var - tmp
                    else if(attr(gtheta0[[1L]],"Proportion Accepted") < .1)
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
        for(i in seq_len(31L)){
            gtheta0[[1L]] <- draw.thetas(theta0=gtheta0[[1L]], pars=pars[[1L]], fulldata=Data$fulldata[[1L]],
                                         itemloc=itemloc, cand.t.var=cand.t.var, CUSTOM.IND=CUSTOM.IND,
                                         prior.t.var=gstructgrouppars[[1L]]$gcov, OffTerm=OffTerm,
                                         prior.mu=gstructgrouppars[[1L]]$gmeans, prodlist=prodlist)
            if(is.null(list$cand.t.var)){
                if(i > 5L){
                    if(attr(gtheta0[[1L]],"Proportion Accepted") > .35) cand.t.var <- cand.t.var + 2*tmp
                    else if(attr(gtheta0[[1L]],"Proportion Accepted") > .25 && nfact > 3L)
                        cand.t.var <- cand.t.var + tmp
                    else if(attr(gtheta0[[1L]],"Proportion Accepted") < .2 && nfact < 4L)
                        cand.t.var <- cand.t.var - tmp
                    else if(attr(gtheta0[[1L]],"Proportion Accepted") < .1)
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
    for(g in seq_len(ngroups)){
        for(i in seq_len(MHDRAWS))
            gtheta0[[g]] <- draw.thetas(theta0=gtheta0[[g]], pars=pars[[g]], fulldata=Data$fulldata[[g]],
                                        itemloc=itemloc, cand.t.var=cand.t.var, CUSTOM.IND=CUSTOM.IND,
                                        prior.t.var=gstructgrouppars[[g]]$gcov, OffTerm=OffTerm,
                                        prior.mu=gstructgrouppars[[g]]$gmeans, prodlist=prodlist)
        LL <- LL + attr(gtheta0[[g]], "log.lik")
    }
    if(RAND && cycles > RANDSTART){
        for(j in seq_len(length(random))){
            for(i in seq_len(MHDRAWS)){
                random[[j]]@drawvals <- DrawValues(random[[j]], Theta=gtheta0[[1L]], itemloc=itemloc,
                                                   pars=pars[[1L]], fulldata=Data$fulldata[[1L]],
                                                   offterm0=OffTerm, CUSTOM.IND=CUSTOM.IND)
                OffTerm <- OffTerm(random, J=J, N=N)
            }
        }
    }
    if(LR.RAND && cycles > RANDSTART){
        for(j in seq_len(length(lr.random))){
            for(i in seq_len(MHDRAWS)){
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
            for(j in seq_len(length(random))){
                random[[j]]@cand.t.var <- controlCandVar(
                    attr(random[[j]]@drawvals, "Proportion Accepted"),
                    random[[j]]@cand.t.var, min = .01, max = .5)
            }
        }
        if(LR.RAND && cycles > (RANDSTART + 1L)){
            for(j in seq_len(length(lr.random))){
                lr.random[[j]]@cand.t.var <- controlCandVar(
                    attr(lr.random[[j]]@drawvals, "Proportion Accepted"),
                    lr.random[[j]]@cand.t.var, min = .01, max = .5)
            }
        }
    }
    if(is.na(attr(gtheta0[[1L]],"log.lik")))
        stop('Estimation halted. Model did not converge.', call.=FALSE)
    printmsg <- ""
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
    list(pars=pars, gstructgrouppars=gstructgrouppars, lr.random=lr.random, random=random,
         lrPars=lrPars, gtheta0=gtheta0, OffTerm=OffTerm, printmsg=printmsg, cand.t.var=cand.t.var)
}

MHRM.Mstep <- function(pars, gtheta, OffTerm, longpars, USE.FIXED, list, ngroups, LBOUND, UBOUND,
                       DERIV, gstructgrouppars, CUSTOM.IND, RAND, has_graded, nfact,
                       cycles, RANDSTART, random, J, LRPARS, lrPars, LR.RAND, lr.random,
                       constrain, estpars, redun_constr, L, estindex_unique, gfulldata, itemloc, control){
    Moptim <- list$Moptim
    shortpars <- longpars[estindex_unique]
    lbound <- LBOUND[estindex_unique]
    ubound <- UBOUND[estindex_unique]
    grad <- numeric(1)
    hess <- matrix(0)
    if(Moptim == 'NR1'){
        tmp <- MHRM.deriv(pars=pars, gtheta=gtheta, lrPars=lrPars, OffTerm=OffTerm, longpars=longpars,
                          USE.FIXED=USE.FIXED, list=list, ngroups=ngroups, LR.RAND=LR.RAND,
                          DERIV=DERIV, gstructgrouppars=gstructgrouppars, nfact=nfact,
                          CUSTOM.IND=CUSTOM.IND, RAND=RAND, cycles=cycles, lr.random=lr.random,
                          RANDSTART=RANDSTART, random=random, J=J, LRPARS=LRPARS, L=L,
                          constrain=constrain, estpars=estpars, redun_constr=redun_constr)
        grad <- tmp$grad
        hess <- tmp$ave.h
        correction <- try(solve(hess, grad), TRUE)
        if(is(correction, 'try-error')){
            ave.h.inv <- MPinv(hess)
            correction <- as.vector(grad %*% ave.h.inv)
        }
    } else if(Moptim == 'BFGS'){
        opt <- optim(shortpars, MHRM.LL_reload, gr=MHRM.deriv_reload, pars=pars,
                     gtheta=gtheta, lrPars=lrPars, OffTerm=OffTerm, longpars=longpars,
                     USE.FIXED=USE.FIXED, list=list, ngroups=ngroups, LR.RAND=LR.RAND,
                     DERIV=DERIV, has_graded=has_graded, nfact=nfact,
                     CUSTOM.IND=CUSTOM.IND, RAND=RAND, cycles=cycles, lr.random=lr.random,
                     RANDSTART=RANDSTART, random=random, J=J, LRPARS=LRPARS, L=L,
                     constrain=constrain, estpars=estpars, redun_constr=redun_constr,
                     estindex_unique=estindex_unique, gfulldata=gfulldata, itemloc=itemloc,
                     method='BFGS', control=control)
        correction <- opt$par - shortpars
    } else if(Moptim == 'L-BFGS-B'){
        opt <- optim(shortpars, MHRM.LL_reload, gr=MHRM.deriv_reload, pars=pars,
                     gtheta=gtheta, lrPars=lrPars, OffTerm=OffTerm, longpars=longpars,
                     USE.FIXED=USE.FIXED, list=list, ngroups=ngroups, LR.RAND=LR.RAND,
                     DERIV=DERIV, has_graded=has_graded, nfact=nfact,
                     CUSTOM.IND=CUSTOM.IND, RAND=RAND, cycles=cycles, lr.random=lr.random,
                     RANDSTART=RANDSTART, random=random, J=J, LRPARS=LRPARS, L=L,
                     constrain=constrain, estpars=estpars, redun_constr=redun_constr,
                     estindex_unique=estindex_unique, gfulldata=gfulldata, itemloc=itemloc,
                     method='L-BFGS-B', lower=lbound, upper=ubound, control=control)
        correction <- opt$par - shortpars
    } else if(Moptim == 'NR'){
        correction <- MHRM.NR(shortpars=shortpars, pars=pars,
                              gtheta=gtheta, lrPars=lrPars, OffTerm=OffTerm, longpars=longpars,
                              USE.FIXED=USE.FIXED, list=list, ngroups=ngroups, LR.RAND=LR.RAND,
                              DERIV=DERIV, has_graded=has_graded, nfact=nfact,
                              CUSTOM.IND=CUSTOM.IND, RAND=RAND, cycles=cycles, lr.random=lr.random,
                              RANDSTART=RANDSTART, random=random, J=J, LRPARS=LRPARS, L=L,
                              constrain=constrain, estpars=estpars, redun_constr=redun_constr,
                              estindex_unique=estindex_unique, gfulldata=gfulldata, itemloc=itemloc,
                              lbound=lbound, ubound=ubound, control=control)
    } else {
        stop('Optimizer currently not supported for stochastic optimization method', call.=FALSE)
    }
    list(correction=correction, grad=grad, hess=hess)
}
