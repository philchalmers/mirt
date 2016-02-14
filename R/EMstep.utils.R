Estep <- function(pars, Data, Theta, prior, Prior, Priorbetween, specific, sitems,
                  itemloc, CUSTOM.IND, BFACTOR, ngroups, rlist, full, Etable = TRUE){
    LL <- 0
    tabdata <- if(full) Data$fulldata[[1L]] else Data$tabdatalong
    for(g in 1L:ngroups){
        freq <- if(full) 1 else Data$Freq[[g]]
        if(BFACTOR){
            rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=tabdata, freq=Data$Freq[[g]],
                                        Theta=Theta, prior=prior[[g]], Prior=Prior[[g]],
                                        Priorbetween=Priorbetween[[g]], specific=specific,
                                        sitems=sitems, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                        Etable=Etable)
        } else {
            rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=tabdata, freq=Data$Freq[[g]],
                                     CUSTOM.IND=CUSTOM.IND, Theta=Theta,
                                     prior=Prior[[g]], itemloc=itemloc, full=full, Etable=Etable)
        }
        LL <- LL + sum(freq * log(rlist[[g]]$expected), na.rm = TRUE)
        rlist[[g]]$r1[is.nan(rlist[[g]]$r1)] <- 0
    }
    return(list(rlist=rlist, LL=LL))
}

# Estep for mirt
Estep.mirt <- function(pars, tabdata, freq, Theta, prior, itemloc, CUSTOM.IND, full = FALSE,
                       itemtrace=NULL, deriv = FALSE, Etable = TRUE)
{
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
    retlist <- if(full) .Call("Estep2", itemtrace, prior, tabdata, Etable)
        else .Call("Estep", itemtrace, prior, tabdata, freq, Etable)
    if(deriv) retlist$itemtrace <- itemtrace
    return(retlist)
}

# Estep for bfactor
Estep.bfactor <- function(pars, tabdata, freq, Theta, prior, Prior, Priorbetween, specific,
                          CUSTOM.IND, sitems, itemloc, itemtrace=NULL, Etable = TRUE)
{
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
    retlist <- .Call("Estepbfactor", itemtrace, prior, Priorbetween, tabdata, freq, sitems, Prior, Etable)
    return(retlist)
}

Mstep <- function(pars, est, longpars, ngroups, J, gTheta, itemloc, PrepList, L, ANY.PRIOR,
                  UBOUND, LBOUND, constrain, DERIV, Prior, rlist, CUSTOM.IND, solnp_args,
                  SLOW.IND, BFACTOR, nfact, Thetabetween, Moptim, Mrate, TOL, full,
                  lrPars, control){
    p <- longpars[est]
    if(length(p)){
        if(Moptim == 'BFGS'){
            if(is.null(control$maxit))
                control$maxit <- max(ceiling(Mrate * 50), 15)
            opt <- try(optim(p, fn=Mstep.LL, gr=Mstep.grad, method='BFGS', control=control,
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc),
                       silent=TRUE)
        } else if(Moptim == 'L-BFGS-B'){
            if(is.null(control$maxit))
                control$maxit <- max(ceiling(Mrate * 50), 15)
            opt <- try(optim(p, fn=Mstep.LL, gr=Mstep.grad, method='L-BFGS-B', control=control,
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc, lower=LBOUND[est],
                             upper=UBOUND[est]),
                       silent=TRUE)
        } else if(Moptim == 'Nelder-Mead'){
            opt <- try(optim(p, fn=Mstep.LL, method='Nelder-Mead', control=control,
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc),
                       silent=TRUE)
        } else if(Moptim == 'SANN'){
            opt <- try(optim(p, fn=Mstep.LL, method='SANN', control=control,
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc),
                       silent=TRUE)
        } else if(Moptim == 'NR'){
            opt <- try(Mstep.NR(p=p, est=est, longpars=longpars, pars=pars, ngroups=ngroups,
                                J=J, gTheta=gTheta, PrepList=PrepList, L=L,  ANY.PRIOR=ANY.PRIOR,
                                constrain=constrain, LBOUND=LBOUND, UBOUND=UBOUND, SLOW.IND=SLOW.IND,
                                itemloc=itemloc, DERIV=DERIV, rlist=rlist, TOL=TOL, control=control))
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
                                     control = control, optim_args=optim_args), silent=TRUE)
                    if(!is(opt, 'try-error')) opt$par <- opt$pars
                } else {
                    stop('Rsolnp package is not available. Please install.', call.=FALSE)
                }
            } else {
                if(requireNamespace("alabama", quietly = TRUE)){
                    opt <- try(alabama::constrOptim.nl(p, Mstep.LL_alt, Mstep.grad_alt,
                                              hin = solnp_args$hin, hin.jac = solnp_args$hin.jac,
                                              heq = solnp_args$heq, heq.jac = solnp_args$heq.jac,
                                              control.outer = solnp_args$control.outer,
                                              control.optim = solnp_args$control.optim, optim_args=optim_args),
                               silent=TRUE)
                } else {
                    stop('alabama package is not available. Please install.', call.=FALSE)
                }
            }
        } else if(Moptim == 'nlminb'){
            opt <- try(nlminb(p, Mstep.LL, Mstep.grad,
                              DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                              est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                              PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                              UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc, lower=LBOUND[est],
                              upper=UBOUND[est], control=control),
                       silent=TRUE)
        } else {
            stop('M-step optimizer not supported', call.=FALSE)
        }
        if(is(opt, 'try-error'))
            stop(opt, call.=FALSE)
        longpars[est] <- opt$par
    }
    if(full){
        res <- Mstep.LR(Theta=gTheta[[1L]], CUSTOM.IND=CUSTOM.IND, pars=pars[[1L]], lrPars=lrPars,
                        itemloc=itemloc, fulldata=PrepList[[1L]]$fulldata, prior=Prior[[1L]])
        longpars[lrPars@parnum] <- res$beta
        longpars[pars[[1L]][[J+1L]]@parnum[pars[[1L]][[J+1L]]@est]] <-
            res$siglong[pars[[1L]][[J+1L]]@est]
    }
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    return(longpars)
}

Mstep.LL <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L, CUSTOM.IND,
                     SLOW.IND, constrain, LBOUND, UBOUND, itemloc, DERIV, rlist, ANY.PRIOR){
    longpars[est] <- p
    if(any(longpars > UBOUND | longpars < LBOUND)) return(-1e100)
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    LLs <- numeric(ngroups)
    for(g in 1L:ngroups){
        LLs[g] <- LogLikMstep(pars[[g]], Theta=gTheta[[g]], rs=rlist[[g]],
                              itemloc=itemloc, CUSTOM.IND=CUSTOM.IND, any.prior=ANY.PRIOR[g])
        if(any(pars[[g]][[J+1L]]@est))
            LLs[g] <- LLs[g] + Mstep.LL.group(pars=pars[[g]], Theta=gTheta[[g]],
                                              rr=rowSums(rlist[[g]]$r1))
    }
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

Mstep.LL.group <- function(pars, Theta, rr){
    theta <- Theta
    pick <- length(pars)
    gp <- ExtractGroupPars(pars[[pick]])
    chl <- try(chol(gp$gcov), silent=TRUE)
    if(is(chl, 'try-error')) return(-1e100)
    tmp <- outer(diag(gp$gcov), diag(gp$gcov))
    if(any(gp$gcov[lower.tri(tmp)] >= tmp[lower.tri(tmp)])) return(-1e100)
    tmp <- rr * mirt_dmvnorm(theta, gp$gmeans, gp$gcov, log=TRUE)
    LL <- sum(tmp)
    if(pars[[pick]]@any.prior)
        LL <- LL.Priors(x=pars[[pick]], LL=LL)
    return(ifelse(is.nan(LL), -1e100, LL))
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
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    g <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, J), length(est), 0L, 0L, 1L, TRUE)$grad
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
                     constrain, LBOUND, UBOUND, itemloc, DERIV, rlist, SLOW.IND, TOL, control)
{
    plast2 <- plast <- p
    ubound <- UBOUND[est]
    lbound <- LBOUND[est]
    lastchange <- 0
    if(is.null(control$maxit)) control$maxit <- 50L
    for(iter in 1L:control$maxit){
        longpars[est] <- p
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        dd <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, J), length(est), 1L, 0L, 1L, TRUE)
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
        p[p > ubound] <- ubound[p > ubound]
        p[p < lbound] <- lbound[p < lbound]
        dif <- plast - p
        if(all(abs(dif) < TOL)) break
        lastchange <- change
    }
    return(list(par=p))
}

BL.grad <- function(x, ...){
    numerical_deriv(x, BL.LL, ...)
}

Mstep.LR <- function(Theta, CUSTOM.IND, pars, itemloc, fulldata, prior, lrPars, retscores=FALSE){
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
    mu <- lrPars@mus
    X <- lrPars@X
    ret <- .Call('EAPgroup', itemtrace, fulldata, Theta, prior, mu)
    scores <- ret[[1L]]; vars <- ret[[2L]]
    if(retscores) return(scores)
    beta <- lrPars@inv_tXX %*% t(X) %*% scores
    siglong <- colMeans(vars)
    beta[!lrPars@est] <- lrPars@par[!lrPars@est]
    return(list(beta=beta, siglong=c(rep(0, ncol(Theta)), siglong)))
}
