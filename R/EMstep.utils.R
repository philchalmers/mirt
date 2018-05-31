Estep <- function(pars, Data, gTheta, prior, Prior, Priorbetween, specific, sitems,
                  itemloc, CUSTOM.IND, dentype, ngroups, rlist, full, Etable = TRUE){
    LL <- 0
    tabdata <- if(full) Data$fulldata[[1L]] else Data$tabdatalong
    if(dentype == 'mixture'){
        rlist <- Estep.mixture(pars=pars, tabdata=tabdata, freq=Data$Freq[[1L]],
                               CUSTOM.IND=CUSTOM.IND, Theta=gTheta[[1L]],
                               prior=Prior, itemloc=itemloc, full=full, Etable=Etable)
        LL <- sum(Data$Freq[[1L]] * log(rlist[[1L]]$expected))
    } else {
        for(g in seq_len(ngroups)){
            freq <- if(full) 1 else Data$Freq[[g]]
            if(dentype == 'bfactor'){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=tabdata, freq=Data$Freq[[g]],
                                            Theta=gTheta[[g]], prior=prior[[g]],
                                            Priorbetween=Priorbetween[[g]], specific=specific,
                                            sitems=sitems, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                            Etable=Etable)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=tabdata, freq=Data$Freq[[g]],
                                         CUSTOM.IND=CUSTOM.IND, Theta=gTheta[[g]],
                                         prior=Prior[[g]], itemloc=itemloc, full=full, Etable=Etable)
            }
            LL <- LL + sum(freq * log(rlist[[g]]$expected), na.rm = TRUE)
            rlist[[g]]$r1[is.nan(rlist[[g]]$r1)] <- 0
        }
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
Estep.bfactor <- function(pars, tabdata, freq, Theta, prior, Priorbetween, specific,
                          CUSTOM.IND, sitems, itemloc, itemtrace=NULL, Etable = TRUE)
{
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
    retlist <- .Call("Estepbfactor", itemtrace, prior, Priorbetween, tabdata, freq, sitems, Etable)
    return(retlist)
}

# Estep for mixture Gaussian
Estep.mixture <- function(pars, tabdata, freq, Theta, prior, itemloc, CUSTOM.IND, full = FALSE,
                       itemtrace=NULL, deriv = FALSE, Etable = TRUE)
{
    ngroups <- length(pars)
    if(is.null(itemtrace)){
        itemtrace <- vector('list', ngroups)
        for(g in seq_len(ngroups))
            itemtrace[[g]] <- computeItemtrace(pars=pars[[g]], Theta=Theta,
                                               itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
    }
    tmp <- .Call("Estep", do.call(rbind, itemtrace),
                 do.call(c, prior), tabdata, freq, Etable)
    retlist <- vector('list', ngroups)
    nrows <- nrow(itemtrace[[1L]])
    for(g in seq_len(ngroups))
        retlist[[g]] <- list(r1=tmp$r1[1L:nrows + (g-1)*nrows, ],
                             expected= if(g == 1) tmp$expected else NA)
    return(retlist)
}

Mstep <- function(pars, est, longpars, ngroups, J, gTheta, itemloc, PrepList, L, ANY.PRIOR,
                  UBOUND, LBOUND, constrain, DERIV, Prior, rlist, CUSTOM.IND, solnp_args,
                  SLOW.IND, dentype, nfact, Moptim, Mrate, TOL, full,
                  lrPars, keep_vcov_PD, control){
    p <- longpars[est]
    if(length(p)){
        if(Moptim == 'BFGS'){
            if(is.null(control$maxit))
                control$maxit <- max(ceiling(Mrate * 50), 15)
            opt <- try(optim(p, fn=Mstep.LL, gr=Mstep.grad, method='BFGS', control=control,
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             itemloc=itemloc, keep_vcov_PD=keep_vcov_PD),
                       silent=TRUE)
        } else if(Moptim == 'L-BFGS-B'){
            if(is.null(control$maxit))
                control$maxit <- max(ceiling(Mrate * 50), 15)
            opt <- try(optim(p, fn=Mstep.LL, gr=Mstep.grad, method='L-BFGS-B', control=control,
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             itemloc=itemloc, keep_vcov_PD=keep_vcov_PD, lower=LBOUND[est], upper=UBOUND[est]),
                       silent=TRUE)
        } else if(Moptim == 'Nelder-Mead'){
            opt <- try(optim(p, fn=Mstep.LL, method='Nelder-Mead', control=control,
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             itemloc=itemloc, keep_vcov_PD=keep_vcov_PD),
                       silent=TRUE)
        } else if(Moptim == 'SANN'){
            opt <- try(optim(p, fn=Mstep.LL, method='SANN', control=control,
                             DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                             est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                             PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                             itemloc=itemloc, keep_vcov_PD=keep_vcov_PD),
                       silent=TRUE)
        } else if(Moptim == 'NR'){
            opt <- try(Mstep.NR(p=p, est=est, longpars=longpars, pars=pars, ngroups=ngroups,
                                J=J, gTheta=gTheta, PrepList=PrepList, L=L,  ANY.PRIOR=ANY.PRIOR,
                                constrain=constrain, LBOUND=LBOUND, UBOUND=UBOUND, SLOW.IND=SLOW.IND,
                                itemloc=itemloc, DERIV=DERIV, rlist=rlist, control=control), TRUE)
        } else if(Moptim %in% c('solnp', 'nloptr')){
            optim_args <- list(CUSTOM.IND=CUSTOM.IND, est=est, longpars=longpars, pars=pars,
                               ngroups=ngroups, J=J, gTheta=gTheta, PrepList=PrepList, L=L,
                               ANY.PRIOR=ANY.PRIOR, constrain=constrain, LBOUND=LBOUND,
                               UBOUND=UBOUND, SLOW.IND=SLOW.IND, itemloc=itemloc, DERIV=DERIV,
                               rlist=rlist, keep_vcov_PD=keep_vcov_PD)
            if(Moptim == 'solnp'){
                if(requireNamespace("Rsolnp", quietly = TRUE)){
                    opt <- try(Rsolnp::solnp(p, Mstep.LL_alt, eqfun = solnp_args$eqfun, eqB = solnp_args$eqB,
                                     ineqfun = solnp_args$ineqfun, ineqLB = solnp_args$ineqLB,
                                     ineqUB = solnp_args$ineqUB, LB = solnp_args$LB, UB = solnp_args$UB,
                                     control = solnp_args$control, optim_args=optim_args), silent=TRUE)
                    if(!is(opt, 'try-error')) opt$par <- opt$pars
                } else {
                    stop('Rsolnp package is not available. Please install.', call.=FALSE)
                }
            } else {
                if(requireNamespace("nloptr", quietly = TRUE)){
                    opt <- try(nloptr::nloptr(p, Mstep.LL_alt, Mstep.grad_alt,
                                              lb=solnp_args$lb,
                                              ub=solnp_args$ub,
                                              eval_g_ineq=solnp_args$eval_g_ineq,
                                              eval_jac_g_ineq=solnp_args$eval_jac_g_ineq,
                                              eval_g_eq=solnp_args$eval_g_eq,
                                              eval_jac_g_eq=solnp_args$eval_jac_g_eq,
                                              opts=solnp_args$opts,
                                              optim_args=optim_args),
                               silent=TRUE)
                    if(!is(opt, 'try-error')) opt$par <- opt$solution
                } else {
                    stop('nloptr package is not available. Please install.', call.=FALSE)
                }
            }
        } else if(Moptim == 'nlminb'){
            opt <- try(nlminb(p, Mstep.LL, Mstep.grad,
                              DERIV=DERIV, rlist=rlist, CUSTOM.IND=CUSTOM.IND, SLOW.IND=SLOW.IND,
                              est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                              PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                              itemloc=itemloc, keep_vcov_PD=keep_vcov_PD,
                              lower=LBOUND[est], upper=UBOUND[est], control=control),
                       silent=TRUE)
        } else {
            stop('M-step optimizer not supported', call.=FALSE)
        }
        if(is(opt, 'try-error'))
            stop(opt, call.=FALSE)

        longpars[est] <- opt$par
    }
    if (dentype == "Davidian") {
        optDC <- Mstep.DC.optim(pars=pars, J=J, gTheta=gTheta, constrain=constrain)
        for(g in seq_len(length(pars))){
            tmp <- pars[[g]][[J+1L]]@est
            pick <- pars[[g]][[J+1L]]@parnum[-c(1L:2L)]
            pick2 <- 1L:length(pick) + length(pick) * (g-1L)
            pick3 <- pars[[g]][[J+1L]]@parnum
            longpars[pick3[tmp]] <- c(optDC$mean_var[g, ], optDC$par[pick2])[tmp]
        }
    }
    if (dentype == "mixture") {
        optMix <- Mstep.mixture(pars=pars, Prior=Prior,
                                J=J, gTheta=gTheta, constrain=constrain)
        for(g in seq_len(length(pars))){
            tmp <- pars[[g]][[J+1L]]@est
            pick <- pars[[g]][[J+1L]]@parnum[length(tmp)]
            if(tmp[length(tmp)])
                longpars[pick] <- optMix[g]
        }
    }
    if(full){
        res <- Mstep.LR(Theta=gTheta[[1L]], CUSTOM.IND=CUSTOM.IND, pars=pars[[1L]], lrPars=lrPars,
                        itemloc=itemloc, fulldata=PrepList[[1L]]$fulldata, prior=Prior[[1L]])
        longpars[lrPars@parnum] <- res$beta
        if(dentype != 'discrete')
            longpars[pars[[1L]][[J+1L]]@parnum[pars[[1L]][[J+1L]]@est]] <-
                res$siglong[pars[[1L]][[J+1L]]@est]
    }
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    longpars
}

Mstep.LL <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L, CUSTOM.IND,
                     SLOW.IND, constrain, itemloc, DERIV, rlist, ANY.PRIOR, keep_vcov_PD){
    longpars[est] <- p
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    LLs <- numeric(ngroups)
    for(g in seq_len(ngroups)){
        LLs[g] <- LogLikMstep(pars[[g]], Theta=gTheta[[g]], rs=rlist[[g]],
                              itemloc=itemloc, CUSTOM.IND=CUSTOM.IND, any.prior=ANY.PRIOR[g])
        if(any(pars[[g]][[J+1L]]@est))
            LLs[g] <- LLs[g] + Mstep.LL.group(pars=pars[[g]], Theta=gTheta[[g]],
                                              keep_vcov_PD=keep_vcov_PD)
    }
    return(-sum(LLs))
}

Mstep.LL_alt <- function(x0, optim_args){
    return(Mstep.LL(p=x0, est=optim_args$est, longpars=optim_args$longpars, pars=optim_args$pars,
                    ngroups=optim_args$ngroups, J=optim_args$J, gTheta=optim_args$gTheta,
                    PrepList=optim_args$PrepList, L=optim_args$L, CUSTOM.IND=optim_args$CUSTOM.IND,
                    SLOW.IND=optim_args$SLOW.IND, keep_vcov_PD=optim_args$keep_vcov_PD,
                    constrain=optim_args$constrain, itemloc=optim_args$itemloc,
                    DERIV=optim_args$DERIV, rlist=optim_args$rlist, ANY.PRIOR=optim_args$ANY.PRIOR))
}

Mstep.LL.group <- function(pars, Theta, keep_vcov_PD){
    pick <- length(pars)
    rr <- pars[[pick]]@rr
    if(pars[[pick]]@itemclass < 0L){
        gp <- pars[[pick]]
        d <- gp@safe_den(gp, Theta)
        LL <- sum(rr * log(d))
    } else {
        gp <- ExtractGroupPars(pars[[pick]])
        chl <- try(chol(gp$gcov), silent=TRUE)
        if(is(chl, 'try-error')){
            if(keep_vcov_PD){
                sds <- diag(sqrt(diag(gp$gcov)))
                smoothed <- cov2cor(smooth.cov(gp$gcov))
                gp$gcov <- sds %*% smoothed %*% sds
            } else return(-1e100)
        }
        if (pars[[pick]]@BFACTOR) {
            theta <- pars[[pick]]@theta
            Thetabetween <- pars[[pick]]@Thetabetween
            nsfact <- ncol(pars[[pick]]@rrs)
            npfact <- ncol(Theta) - nsfact
            pick2 <- 1L:npfact
            LL <- sum(pars[[pick]]@rrb * mirt_dmvnorm(Thetabetween, gp$gmeans[pick2],
                                                      gp$gcov[pick2,pick2,drop=FALSE], log=TRUE))
            for(i in 1L:nsfact){
                pick2 <- npfact + i
                LL <- LL + sum(pars[[pick]]@rrs[,i] * dnorm(theta, gp$gmeans[pick2],
                                                            sqrt(gp$gcov[pick2,pick2]), log=TRUE))
            }
        } else {
            LL <- sum(rr * mirt_dmvnorm(Theta, gp$gmeans, gp$gcov, log=TRUE))
        }
    }
    if(pars[[pick]]@any.prior)
        LL <- LL.Priors(x=pars[[pick]], LL=LL)
    LL
}

Mstep.DC.optim <- function(pars, J, gTheta, constrain) {
    ngroup <- length(pars)
    orgphi <- std_rr <- vector("list", ngroup)
    for(g in seq_len(ngroup)){
        orgphi[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])$phi
        std_rr[[g]] <- standardizeQuadrature(gTheta[[g]], nq=pars[[g]][[J+1L]]@rr,
                                             estmean=pars[[g]][[J+1L]]@est[1L],
                                             estsd=pars[[g]][[J+1L]]@est[2L])
    }
    phi <- do.call(c, orgphi) # adjust based on constrain TODO
    optDC <- nlminb(phi, function (...) -Mstep.DC.LL_group(...),
                    gradient = function(...) -Mstep.DC.grad_group(...),
                    gTheta=gTheta, rr = std_rr, constrain = constrain,
                    orgphi=orgphi)
    # find equivalent set equal within [-pi/2, pi/2]
    # there's likely an analytical way to do this, but it's not important
    if(!all(optDC$par >= -pi/2 & optDC$par <= pi/2)){
        phi_s <- phi_bound(optDC$par)
        tmp_LL <- 0
        all_signs <- lapply(1L:length(phi_s), function(x) c(-1,1))
        all_signs <- as.matrix(expand.grid(all_signs))
        for(i in seq_len(nrow(all_signs))){
            signs <- all_signs[i, ]
            tmp_LL <- -Mstep.DC.LL_group(phi_s * signs, gTheta=gTheta,
                                         rr = std_rr, constrain = constrain,
                                         orgphi=orgphi)
            if(isTRUE(all.equal(tmp_LL, optDC$objective))) break
        }
        optDC$par <- phi_s * signs
    }
    optDC$value <- optDC$objective
    optDC$mean_var <- do.call(rbind, lapply(std_rr, function(x) attr(x, 'mean_var')))
    optDC # return as long parameter vector for easier transition TODO
}

Mstep.DC.LL <- function (phi, Theta, rr) {
    LL <- sum(rr * log(dcurver::ddc(Theta, phi)))
    LL
}

Mstep.DC.LL_group <- function (phi, gTheta, rr, constrain, orgphi) {
    ind1 <- 1L
    for(g in seq_len(length(gTheta))){
        ind2 <- length(orgphi[[g]]) + ind1 - 1L
        orgphi[[g]] <- phi[ind1:ind2]
        ind1 <- ind2 + 1L
    }
    # constrain.... TODO
    LL <- 0
    for(g in seq_len(length(gTheta)))
        LL <- LL + Mstep.DC.LL(orgphi[[g]], gTheta[[g]], rr[[g]])
    LL
}

Mstep.DC.grad <- function (phi, Theta, rr) {
    colSums(dcurver::dc_grad(Theta, phi) * rr)
}

Mstep.DC.grad_group <- function (phi, gTheta, rr, constrain, orgphi) {
    ind1 <- 1L
    for(g in seq_len(length(gTheta))){
        ind2 <- length(orgphi[[g]]) + ind1 - 1L
        orgphi[[g]] <- phi[ind1:ind2]
        ind1 <- ind2 + 1L
    }
    # constrain.... TODO
    grad <- c()
    for(g in seq_len(length(gTheta)))
        grad <- c(grad, Mstep.DC.grad(orgphi[[g]], gTheta[[g]], rr[[g]]))
    grad
}

Mstep.mixture <- function(pars, Prior, gTheta, J, constrain){
    rrs <- sapply(1L:length(pars),
                  function(g) sum(pars[[g]][[J+1L]]@rr))
    total <- sum(rrs)
    lps <- log(rrs/total)
    lps <- lps - lps[1L]
    lps
}

LogLikMstep <- function(x, Theta, itemloc, rs, any.prior, CUSTOM.IND){
    log_itemtrace <- log(computeItemtrace(pars=x, Theta=Theta,
                                          itemloc=itemloc, CUSTOM.IND=CUSTOM.IND))
    LL <- sum(rs$r1 * log_itemtrace)
    if(any.prior){
        for(i in seq_len(length(x)-1L))
            if(x[[i]]@any.prior)
                LL <- LL.Priors(x=x[[i]], LL=LL)
    }
    return(ifelse(is.finite(LL), LL, -1e100))
}

Mstep.grad <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L, ANY.PRIOR,
                       constrain, itemloc, DERIV, rlist, CUSTOM.IND, SLOW.IND, keep_vcov_PD){
    longpars[est] <- p
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    if(pars[[1L]][[J + 1L]]@itemclass == -1L){
        for(g in seq_len(length(pars))){
            gp <- pars[[g]][[J + 1L]]
            pars[[g]][[J + 1L]]@density <- gp@safe_den(gp, gTheta[[g]])
        }
    }
    g <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, J), length(est), 0L, 0L, 1L, TRUE)$grad
    if(length(SLOW.IND)){
        for(group in seq_len(ngroups)){
            for (i in SLOW.IND){
                deriv <- if(i == (J + 1L)){
                    Deriv(pars[[group]][[i]], Theta=gTheta[[group]])
                } else {
                    DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=gTheta[[group]])
                }
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
                      SLOW.IND=optim_args$SLOW.IND, keep_vcov_PD=optim_args$keep_vcov_PD,
                      constrain=optim_args$constrain, itemloc=optim_args$itemloc,
                      DERIV=optim_args$DERIV, rlist=optim_args$rlist, ANY.PRIOR=optim_args$ANY.PRIOR))
}

Mstep.NR <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L, ANY.PRIOR,
                     constrain, LBOUND, UBOUND, itemloc, DERIV, rlist, SLOW.IND, control)
{
    plast2 <- plast <- p
    ubound <- UBOUND[est]
    lbound <- LBOUND[est]
    lastchange <- 0
    for(iter in seq_len(control$maxit)){
        longpars[est] <- p
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        dd <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, J), length(est), 1L, 0L, 1L, TRUE)
        if(length(SLOW.IND)){
            for(group in seq_len(ngroups)){
                for (i in SLOW.IND){
                    deriv <- if(i == (J + 1L)){
                        Deriv(pars[[group]][[i]], Theta=gTheta[[group]], estHess=TRUE)
                    } else {
                        DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=gTheta[[group]], estHess=TRUE)
                    }
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
        inv <- try(MPinv(h), TRUE) #TODO this could be avoided if no constrains present
        if(is(inv, 'try-error'))
            stop('Newton-Raphson Hessian become non-positive definite in the M-step', call. = FALSE)
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
    return(list(par=p))
}

BL.grad <- function(x, ...){
    numerical_deriv(BL.LL, x, ...)
}

Mstep.LR <- function(Theta, CUSTOM.IND, pars, itemloc, fulldata, prior,
                     lrPars, retscores=FALSE){
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
