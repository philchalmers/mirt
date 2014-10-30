Estep <- function(pars, Data, Theta, prior, Prior, Priorbetween, specific, sitems,
                  itemloc, CUSTOM.IND, BFACTOR, ngroups, rlist, full){
    LL <- 0
    tabdata <- if(full) Data$fulldata[[1L]] else Data$tabdatalong
    for(g in 1L:ngroups){
        freq <- if(full) 1 else Data$Freq[[g]]
        if(BFACTOR){
            rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=tabdata, freq=Data$Freq[[g]],
                                        Theta=Theta, prior=prior[[g]], Prior=Prior[[g]],
                                        Priorbetween=Priorbetween[[g]], specific=specific,
                                        sitems=sitems, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
        } else {
            rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=tabdata, freq=Data$Freq[[g]],
                                     CUSTOM.IND=CUSTOM.IND, Theta=Theta,
                                     prior=Prior[[g]], itemloc=itemloc, full=full)
        }
        LL <- LL + sum(freq * log(rlist[[g]]$expected), na.rm = TRUE)
    }
    return(list(rlist=rlist, LL=LL))
}

# Estep for mirt
Estep.mirt <- function(pars, tabdata, freq, Theta, prior, itemloc, CUSTOM.IND, full,
                       itemtrace=NULL, deriv = FALSE)
{
    nquad <- nrow(Theta)
    J <- length(itemloc) - 1L
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
    retlist <- if(full) .Call("Estep2", itemtrace, prior, tabdata, mirtClusterEnv$ncores)
        else .Call("Estep", itemtrace, prior, tabdata, freq, mirtClusterEnv$ncores)
    if(deriv) retlist$itemtrace <- itemtrace
    return(retlist)
}

# Estep for bfactor
Estep.bfactor <- function(pars, tabdata, freq, Theta, prior, Prior, Priorbetween, specific,
                          CUSTOM.IND, sitems, itemloc, itemtrace=NULL)
{
    J <- length(itemloc) - 1L
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
    retlist <- .Call("Estepbfactor", itemtrace, prior, Priorbetween, tabdata, freq, sitems, Prior,
                     mirtClusterEnv$ncores)
    return(retlist)
}

Mstep <- function(pars, est, longpars, ngroups, J, gTheta, itemloc, PrepList, L, ANY.PRIOR,
                  UBOUND, LBOUND, constrain, DERIV, Prior, rlist, CUSTOM.IND, solnp_args,
                  SLOW.IND, groupest, BFACTOR, nfact, Thetabetween, Moptim, Mrate, TOL, full){
    p <- longpars[est]
    if(length(p)){
        if(Moptim == 'BFGS'){
            maxit <- max(ceiling(Mrate * 50), 15)
            opt <- try(optim(p, fn=Mstep.LL, gr=Mstep.grad, method='BFGS',
                             control=list(maxit=maxit),
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc),
                       silent=TRUE)
        } else if(Moptim == 'L-BFGS-B'){
            maxit <- max(ceiling(Mrate * 50), 15)
            opt <- try(optim(p, fn=Mstep.LL, gr=Mstep.grad, method='L-BFGS-B',
                             control=list(maxit=maxit),
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc, lower=LBOUND[est],
                             upper=UBOUND[est]),
                       silent=TRUE)
        } else if(Moptim == 'Nelder-Mead'){
            opt <- try(optim(p, fn=Mstep.LL, method='Nelder-Mead',
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc),
                       silent=TRUE)
        } else if(Moptim == 'SANN'){
            opt <- try(optim(p, fn=Mstep.LL, method='SANN',
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc),
                       silent=TRUE)
        } else if(Moptim == 'NR'){
            opt <- try(Mstep.NR(p=p, est=est, longpars=longpars, pars=pars, ngroups=ngroups,
                                J=J, gTheta=gTheta, PrepList=PrepList, L=L,  ANY.PRIOR=ANY.PRIOR,
                                constrain=constrain, LBOUND=LBOUND, UBOUND=UBOUND, SLOW.IND=SLOW.IND,
                                itemloc=itemloc, DERIV=DERIV, rlist=rlist, TOL=TOL))
        } else if(Moptim %in% c('solnp', 'alabama')){
            optim_args <- list(CUSTOM.IND=CUSTOM.IND, est=est, longpars=longpars, pars=pars,
                               ngroups=ngroups, J=J, gTheta=gTheta, PrepList=PrepList, L=L,
                               ANY.PRIOR=ANY.PRIOR, constrain=constrain, LBOUND=LBOUND,
                               UBOUND=UBOUND, SLOW.IND=SLOW.IND, itemloc=itemloc, DERIV=DERIV,
                               rlist=rlist)
            if(Moptim == 'solnp'){
                if(requireNamespace("Rsolnp", quietly = TRUE)){
                    opt <- try(Rsolnp::solnp(p, Mstep.LL_alt, eqfun = solnp_args$eqfun, eqB = solnp_args$eqB,
                                     ineqfun = solnp_args$ineqfun, ineqLB = solnp_args$ineqLB,
                                     ineqUB = solnp_args$ineqUB, LB = solnp_args$LB, UB = solnp_args$UB,
                                     control=solnp_args$control, optim_args=optim_args), silent=TRUE)
                    if(!is(opt, 'try-error')) opt$par <- opt$pars
                }
            } else {
                if(requireNamespace("alabama", quietly = TRUE)){
                    opt <- try(alabama::constrOptim.nl(p, Mstep.LL_alt, Mstep.grad_alt,
                                              hin = solnp_args$hin, hin.jac = solnp_args$hin.jac,
                                              heq = solnp_args$heq, heq.jac = solnp_args$heq.jac,
                                              control.outer = solnp_args$control.outer,
                                              control.optim = solnp_args$control.optim, optim_args=optim_args),
                               silent=TRUE)
                }
            }
        } else {
            stop('M-step optimizer not supported')
        }
        if(is(opt, 'try-error'))
            stop(opt)
        longpars[est] <- opt$par
    }
    if(!full){
        if(any(groupest)){
            p <- longpars[groupest]
            maxit <- max(ceiling(Mrate * 100), 35)
            res <- try(nlm(Mstep.LL2, p, pars=pars, Theta=gTheta[[1L]], nfact=nfact, BFACTOR=BFACTOR,
                           constrain=constrain, groupest=groupest, longpars=longpars, rlist=rlist,
                           Thetabetween=Thetabetween, iterlim=maxit),
                       silent=TRUE)
            if(is(res, 'try-error')) stop(res)
            longpars[groupest] <- res$estimate
        }
    } else {
        res <- Mstep.LR(Theta=gTheta[[1L]], CUSTOM.IND=CUSTOM.IND, pars=pars[[1L]],
                        itemloc=itemloc, fulldata=PrepList[[1L]]$fulldata, prior=Prior[[1L]])
        attr(longpars, 'beta') <- res$beta
        longpars[groupest] <- res$siglong
    }
    if(length(constrain))
        for(i in 1L:length(constrain))
            longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
    return(longpars)
}

Mstep.LL <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L, CUSTOM.IND,
                     SLOW.IND, constrain, LBOUND, UBOUND, itemloc, DERIV, rlist, ANY.PRIOR){
    longpars[est] <- p
    if(length(constrain))
        for(i in 1L:length(constrain))
            longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    LLs <- numeric(ngroups)
    for(g in 1L:ngroups)
        LLs[g] <- LogLikMstep(pars[[g]], Theta=gTheta[[g]], rs=rlist[[g]],
                              itemloc=itemloc, CUSTOM.IND=CUSTOM.IND, any.prior=ANY.PRIOR[g])
    return(-sum(LLs))
}

Mstep.LL_alt <- function(x0, optim_args){
    return(Mstep.LL(p=x0, est=optim_args$est, longpars=optim_args$longpars, pars=optim_args$pars,
                    ngroups=optim_args$ngroups, J=optim_args$J, gTheta=optim_args$gTheta,
                    PrepList=optim_args$PrepList, L=optim_args$L, CUSTOM.IND=optim_args$CUSTOM.IND,
                    SLOW.IND=optim_args$SLOW.IND,
                    constrain=optim_args$constrain, LBOUND=optim_args$LBOUND,
                    UBOUND=optim_args$UBOUND, itemloc=optim_args$itemloc,
                    DERIV=optim_args$DERIV, rlist=optim_args$rlist, ANY.PRIOR=optim_args$ANY.PRIOR))
}

Mstep.LL2 <- function(p, longpars, pars, Theta, BFACTOR, nfact, constrain, groupest, rlist,
                      Thetabetween){
    ngroups <- length(pars); J <- length(pars[[1L]]) - 1L
    longpars[groupest] <- p
    if(length(constrain))
        for(i in 1L:length(constrain))
            longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    LL <- 0
    ind <- 1L
    for(g in 1L:ngroups){
        if(BFACTOR){
            theta <- Thetabetween
            rr <- rlist[[g]]$r2
        } else {
            theta <- Theta[ ,1L:nfact,drop=FALSE]
            rr <- rlist[[g]]$r1
        }
        est <- pars[[g]][[J+1L]]@est
        nest <- sum(est)
        if(nest){
            pars[[g]][[J+1L]]@par[est] <- p[ind:(nest + ind - 1L)]
            ind <- ind + nest
        } else next
        gp <- ExtractGroupPars(pars[[g]][[J+1L]])
        chl <- try(chol(gp$gcov), silent=TRUE)
        if(is(chl, 'try-error')) return(1e100)
        tmp <- outer(diag(gp$gcov), diag(gp$gcov))
        if(any(gp$gcov[lower.tri(tmp)] >= tmp[lower.tri(tmp)])) return(1e100)
        tmp <- rr * mirt_dmvnorm(theta, gp$gmeans[1L:ncol(theta)],
                                 gp$gcov[1L:ncol(theta),1L:ncol(theta), drop=FALSE], log=TRUE)
        LL <- LL + sum(tmp)
        if(pars[[g]][[J+1L]]@any.prior)
            LL <- LL.Priors(x=pars[[g]][[J+1L]], LL=LL)
    }
    return(ifelse(is.nan(LL), 1e100, -LL))
}

LogLikMstep <- function(x, Theta, itemloc, rs, any.prior, CUSTOM.IND){
    log_itemtrace <- log(computeItemtrace(pars=x, Theta=Theta,
                                          itemloc=itemloc, CUSTOM.IND=CUSTOM.IND))
    LL <- sum(rs$r1 * log_itemtrace)
    if(any.prior){
        for(i in 1L:(length(x)-1L))
            if(x[[i]]@any.prior)
                LL <- LL.Priors(x=x[[i]], LL=LL)
    }
    return(ifelse(is.finite(LL), LL, -1e100))
}

Mstep.grad <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L, ANY.PRIOR,
                       constrain, LBOUND, UBOUND, itemloc, DERIV, rlist, CUSTOM.IND,
                       SLOW.IND){
    longpars[est] <- p
    if(length(constrain) > 0L)
        for(i in 1L:length(constrain))
            longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    g <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, J), length(est), 0L, 0L)$grad
    if(length(SLOW.IND)){
        for(group in 1L:ngroups){
            for (i in SLOW.IND){
                deriv <- DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=gTheta[[group]])
                g[pars[[group]][[i]]@parnum] <- deriv$grad
            }
        }
    }
    if(length(constrain)){
        grad <- g %*% L
    } else {
        grad <- g
    }
    return(-grad[est])
}

Mstep.grad_alt <- function(x0, optim_args){
    return(Mstep.grad(p=x0, est=optim_args$est, longpars=optim_args$longpars, pars=optim_args$pars,
                      ngroups=optim_args$ngroups, J=optim_args$J, gTheta=optim_args$gTheta,
                      PrepList=optim_args$PrepList, L=optim_args$L, CUSTOM.IND=optim_args$CUSTOM.IND,
                      SLOW.IND=optim_args$SLOW.IND,
                      constrain=optim_args$constrain, LBOUND=optim_args$LBOUND,
                      UBOUND=optim_args$UBOUND, itemloc=optim_args$itemloc,
                      DERIV=optim_args$DERIV, rlist=optim_args$rlist, ANY.PRIOR=optim_args$ANY.PRIOR))
}

Mstep.NR <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L,  ANY.PRIOR,
                     constrain, LBOUND, UBOUND, itemloc, DERIV, rlist, NO.CUSTOM, SLOW.IND,
                     TOL)
{
    plast2 <- plast <- p
    for(iter in 1L:50L){
        longpars[est] <- p
        if(length(constrain) > 0L)
            for(i in 1L:length(constrain))
                longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        dd <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, J), length(est), 1L, 0L)
        if(length(SLOW.IND)){
            for(group in 1L:ngroups){
                for (i in SLOW.IND){
                    deriv <- DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=gTheta[[group]],
                                                 estHess=TRUE)
                    tmp <- pars[[group]][[i]]@parnum
                    dd$grad[tmp] <- deriv$grad
                    dd$hess[tmp, tmp] <- deriv$hess
                }
            }
        }
        if(length(constrain)){
            grad <- updateGrad(dd$grad, L)
            hess <- updateHess(-dd$hess, L)
        } else {
            grad <- dd$grad
            hess <- -dd$hess
        }
        g <- grad[est]
        h <- hess[est, est]
        trychol <- try(chol(h), silent=TRUE)
        if(is(trychol, 'try-error')){
            ev <- eigen(h)
            ev$values[ev$values <= 0] <- .Machine$double.eps*100
            h <- t(ev$vector) %*% diag(ev$values) %*% ev$vector
            trychol <- chol(h)
        }
        change <- as.numeric(g %*% chol2inv(trychol))
        change <- ifelse(change > .25, .25, change)
        change <- ifelse(change < -.25, -.25, change)
        plast2 <- plast
        plast <- p
        p <- p + change
        if(iter > 1L){
            flip <- (sign(lastchange) * sign(change)) == -1L
            p[flip] <- (plast[flip] + p[flip]) / 2
        }
        dif <- plast - p
        if(all(abs(dif) < TOL)) break
        lastchange <- change
    }
    return(list(par=p))
}

BL.grad <- function(x, ...) numDeriv::grad(BL.LL, x=x, ...)

Mstep.LR <- function(Theta, CUSTOM.IND, pars, itemloc, fulldata, prior){
    x <- pars[[length(pars)]]
    X <- x@X
    J <- length(pars) - 1L
    N <- nrow(fulldata)
    nfact <- ncol(Theta)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
    mu <- X %*% x@betas
    ret <- .Call('EAPgroup', itemtrace, fulldata, Theta, prior, mu, mirtClusterEnv$ncores)
    scores <- ret[[1L]]; vars <- ret[[2L]]
    beta <- rep(0, length(x@betas))
    for(i in 1:length(beta))
        beta[i] <- solve(X[,i] %*% X[,i]) %*% X[,i] %*% scores
    siglong <- colMeans(vars)
    siglong <- siglong[x@est[-(1L:nfact)]]
    return(list(beta=beta, siglong=siglong))
}
