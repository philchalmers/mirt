SE.Numerical <- function(pars, Theta, theta, dentype, itemloc, PrepList, ESTIMATE, constrain, Ls,
                  CUSTOM.IND, specific=NULL, sitems=NULL, EHPrior = NULL, warn, Data, type,
                  delta, lrPars){
    longpars <- ESTIMATE$longpars
    rlist <- ESTIMATE$rlist
    infological=ESTIMATE$infological
    ngroups <- length(pars)
    J <- length(pars[[1L]]) - 1L
    est <- c()
    for(g in 1L:ngroups)
        for(j in 1L:(J+1L))
            est <- c(est, pars[[g]][[j]]@est)
    shortpars <- longpars[est]
    gstructgrouppars <- vector('list', ngroups)
    for(g in 1L:ngroups)
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
    for(g in 1L:ngroups){
        for(i in 1L:J){
            tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
            pars[[g]][[i]]@dat <- rlist[[g]]$r1[, tmp]
        }
    }
    if(length(lrPars)){
        est <- c(est, lrPars@est)
        shortpars <- longpars[est]
    }
    hess <- numerical_deriv(shortpars, BL.LL, est=est, longpars=longpars, lrPars=lrPars,
                            pars=pars, ngroups=ngroups, J=J, itemloc=itemloc,
                            Theta=Theta, PrepList=PrepList, dentype=dentype,
                            constrain=constrain,
                            specific=specific, sitems=sitems, CUSTOM.IND=CUSTOM.IND,
                            EHPrior=EHPrior, Data=Data, theta=theta, type=type,
                            delta = delta, gradient = FALSE)
    Hess <- matrix(0, length(longpars), length(longpars))
    Hess[est, est] <- -hess
    Hess <- updateHess(h=Hess, L=Ls$L)
    info <- Hess[infological, infological]
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain, warn=warn)
    return(ESTIMATE)
}

SE.SEM <- function(index, estmat, pars, constrain, Ls, PrepList, list, Theta, theta, dentype, ESTIMATE, DERIV,
                   collectLL, from, to, is.latent, Data, solnp_args, control){
    est <- estmat[,index]
    lrPars <- list$lrPars
    full <- list$full
    TOL <- list$TOL
    dentype <- list$dentype
    itemloc <- list$itemloc
    J <- length(itemloc) - 1L
    L <- Ls$L
    Moptim <- list$Moptim
    sitems <- list$sitems
    specific <- list$specific
    ngroups <- ESTIMATE$ngroups
    NCYCLES <- ESTIMATE$cycles
    if(NCYCLES <= 5L) stop('SEM can not be computed due to short EM history', call.=FALSE)
    prior <- rlist <- vector('list', ngroups)
    estpars <- ESTIMATE$estpars
    redun_constr <- ESTIMATE$redun_constr
    EMhistory <- ESTIMATE$EMhistory
    MLestimates <- EMhistory[nrow(EMhistory), ]
    UBOUND <- ESTIMATE$UBOUND
    LBOUND <- ESTIMATE$LBOUND
    estindex <- estpars
    estindex[estpars] <- est
    rij <- 1
    prodlist <- PrepList[[1L]]$prodlist
    nfact <- ncol(Theta)
    ANY.PRIOR <- rep(FALSE, ngroups)
    converged <- logical(sum(estpars & !redun_constr))
    rijfull <- rep(NA, length(converged))
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta, prodlist)
    gTheta <- vector('list', ngroups)
    for(g in 1L:ngroups){
        ANY.PRIOR[g] <- any(sapply(pars[[g]], function(x) x@any.prior))
        gTheta[[g]] <- Theta
    }
    converged <- FALSE

    for (cycles in from:to){

        longpars <- MLestimates
        longpars[estindex] <- EMhistory[cycles, estindex]
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        tmp <- updatePrior(pars=pars, Theta=Theta, list=list, ngroups=ngroups,
                           nfact=nfact, lrPars=lrPars, J=J, dentype=dentype,
                           sitems=sitems, cycles=cycles, rlist=rlist, full=full)
        prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
        #Estep
        for(g in 1L:ngroups){
            if(dentype == 'bfactor'){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                            Theta=Theta, prior=prior[[g]],
                                            Priorbetween=Priorbetween[[g]], specific=specific, sitems=sitems,
                                            itemloc=itemloc, CUSTOM.IND=list$CUSTOM.IND)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                         CUSTOM.IND=list$CUSTOM.IND, Theta=Theta,
                                         prior=Prior[[g]], itemloc=itemloc, full=full)
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
        longpars <- Mstep(pars=pars, est=estpars, longpars=longpars, ngroups=ngroups, J=J,
                          gTheta=gTheta, itemloc=itemloc, Prior=Prior, ANY.PRIOR=ANY.PRIOR,
                          PrepList=PrepList, L=L, UBOUND=UBOUND, LBOUND=LBOUND, nfact=nfact,
                          rlist=rlist, constrain=constrain, DERIV=DERIV, keep_vcov_PD=list$keep_vcov_PD,
                          CUSTOM.IND=list$CUSTOM.IND, SLOW.IND=list$SLOW.IND, dentype=dentype,
                          Moptim=Moptim, Mrate=1, TOL=list$MSTEPTOL, solnp_args=solnp_args, full=full,
                          lrPars=lrPars, control=control)
        rijlast <- rij
        denom <- (EMhistory[cycles, estindex] - MLestimates[estindex])
        if(denom == 0){
            converged <- FALSE
            break
        }
        rij <- (longpars[estpars & !redun_constr] - MLestimates[estpars & !redun_constr]) / denom
        diff <- abs(rij - rijlast) < TOL
        converged <- diff | converged
        which <- is.na(rijfull) & converged
        rijfull[which] <- rij[which]
        if(all(!is.na(rijfull))){
            converged <- TRUE
            break
        }
    } #END EM

    rijfull[is.na(rijfull)] <- rij[is.na(rijfull)]
    attr(rijfull, 'converged') <- converged
    return(rijfull)
}

SE.simple <- function(PrepList, ESTIMATE, Theta, constrain, Ls, N, type,
                      CUSTOM.IND, SLOW.IND, warn, message, Data, complete){
    pars <- ESTIMATE$pars
    itemloc <- PrepList[[1L]]$itemloc
    ngroups <- length(pars)
    nitems <- length(pars[[1L]]) - 1L
    L <- ESTIMATE$L
    Prior <- ESTIMATE$Prior
    prior <- ESTIMATE$prior
    Priorbetween <- ESTIMATE$Priorbetween
    isbifactor <- length(Priorbetween[[1L]]) > 0L
    if(!isbifactor)
        prior <- Priorbetween <- list(matrix(0))
    sitems <- ESTIMATE$sitems
    iscross <- ifelse(type == 'crossprod', TRUE, FALSE)
    gitemtrace <- vector('list', ngroups)
    for(g in 1L:ngroups){
        gitemtrace[[g]] <- computeItemtrace(pars=pars[[g]], Theta=Theta,
                                            itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
        gp <- ExtractGroupPars(pars[[g]][[nitems+1L]])
        pars[[g]][[nitems+1L]]@mu <- gp$gmeans
        pars[[g]][[nitems+1L]]@sig <- gp$gcov
        pars[[g]][[nitems+1L]]@invsig <- solve(gp$gcov)
        pars[[g]][[nitems+1L]]@meanTheta <- colMeans(Theta)
    }
    npars <- ncol(L)
    gPrior <- t(do.call(rbind, Prior))
    rs <- do.call(rbind, Data$Freq)
    whichitems <- unique(c(CUSTOM.IND, SLOW.IND))
    infolist <- .Call("computeInfo", pars, Theta, gPrior, prior, do.call(rbind, Priorbetween),
                      Data$tabdatalong, rs, sitems, itemloc, gitemtrace, npars, isbifactor, iscross)
    Igrad <- infolist[["Igrad"]]; IgradP <- infolist[["IgradP"]]; Ihess <- infolist[["Ihess"]]
    if(length(whichitems)){
        warning('Internal information matrix computations currently not supported for at
        least one of the supplied items. Information matrix/standard errors not computed', call.=FALSE)
        return(ESTIMATE)
    }
    Igrad <- updateHess(Igrad, L=Ls$L)
    IgradP <- updateHess(IgradP, L=Ls$L)
    Ihess <- updateHess(Ihess, L=Ls$L)
    Igrad <- Igrad[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    IgradP <- IgradP[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    Ihess <- Ihess[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    if(type == 'Louis'){
        info <- -Ihess - IgradP + Igrad
    } else if(type == 'crossprod'){
        info <- Igrad
    } else if(type == 'sandwich'){
        tmp <- solve(-Ihess - IgradP + Igrad)
        info <- solve(tmp %*% Igrad %*% tmp)
    }
    colnames(info) <- rownames(info) <- names(ESTIMATE$correction)
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain, warn=warn)
    return(ESTIMATE)
}

SE.Oakes <- function(pick, pars, L, constrain, est, shortpars, longpars,
                     Theta, list, ngroups, J, dentype, sitems,
                     rlist, full, Data, specific, itemloc, CUSTOM.IND,
                     delta, prior, Prior, Priorbetween, nfact,
                     PrepList, ANY.PRIOR, DERIV, SLOW.IND, Norder, zero_g = NULL){
    row <- 1L
    if(is.null(zero_g)){
        grad <- matrix(0, Norder, length(shortpars))
        signs <- seq(-Norder/2, Norder/2, by = 1L)
        signs <- signs[signs != 0L]
    } else {
        grad <- matrix(0, 2L, length(shortpars))
        grad[1L, ] <- zero_g
        signs <- 1
        row <- 2L
    }
    for(sign in signs){
        if(pick != 0){
            longpars_old <- longpars
            d <- sign * delta
            longpars[which(est)[pick]] <- shortpars[pick] + d
            longpars <- longpars_constrain(longpars, constrain)
            pars <- reloadPars(longpars=longpars, pars=pars,
                               ngroups=ngroups, J=J)
            nms <- names(shortpars)[pick]
            if(grepl('MEAN_', nms) || grepl('COV_', nms)){
                tmp <- updatePrior(pars=pars, Theta=Theta,
                                   list=list, ngroups=ngroups, nfact=nfact,
                                   J=J, dentype=dentype, sitems=sitems, cycles=100L,
                                   rlist=rlist, full=full)
                prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
            }
            Elist <- Estep(pars=pars, Data=Data, Theta=Theta, prior=prior, Prior=Prior,
                           Priorbetween=Priorbetween, specific=specific, sitems=sitems,
                           ngroups=ngroups, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                           dentype=dentype, rlist=rlist, full=full, Etable=list$Etable)
            rlist <- Elist$rlist
            longpars <- longpars_old
            pars <- reloadPars(longpars=longpars, pars=pars,
                               ngroups=ngroups, J=J)
        }

        gTheta <- vector('list', ngroups)
        for(g in 1L:ngroups) gTheta[[g]] <- Theta
        if(pars[[1L]][[J + 1L]]@itemclass == -1L){
            for(g in 1L:length(pars)){
                gp <- pars[[g]][[J + 1L]]
                pars[[g]][[J + 1L]]@density <- gp@safe_den(gp, gTheta[[g]])
            }
        }
        for(g in 1L:ngroups){
            for(i in 1L:J)
                pars[[g]][[i]]@dat <- rlist[[g]]$r1[ , c(itemloc[i]:(itemloc[i+1L] - 1L)),
                                                     drop=FALSE]
            if(dentype == 'bfactor'){
                pars[[g]][[J+1L]]@rrb <- rlist[[g]]$r2
                pars[[g]][[J+1L]]@rrs <- rlist[[g]]$r3
            } else pars[[g]][[J+1L]]@rr <- rowSums(rlist[[g]]$r1) / J
        }
        g <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, J), length(est), 0L, 0L, 1L, TRUE)$grad
        if(length(SLOW.IND)){
            for(group in 1L:ngroups){
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
        tmp <- g %*% L
        if(pick == 0L) return(tmp[est])
        grad[row, ] <- tmp[est]
        row <- row + 1L
    }
    cfs <- get_deriv_coefs(order=Norder)
    ret <- colSums(cfs * grad) / delta
    ret
}

SE.Fisher <- function(PrepList, ESTIMATE, Theta, constrain, Ls, N, CUSTOM.IND, SLOW.IND,
                      warn, Data, full){
    pars <- ESTIMATE$pars
    itemloc <- PrepList[[1L]]$itemloc
    ngroups <- length(pars)
    nitems <- length(pars[[1L]]) - 1L
    L <- Ls$L
    DX <- numeric(ncol(L))
    Prior <- ESTIMATE$Prior
    tabdata <- Data$tabdatalong
    K <- Data$K
    resp <- vector('list', nitems)
    for(i in 1L:nitems)
        resp[[i]] <- 0L:(K[i]-1L)
    resp <- expand.grid(resp)
    stringfulldata <- apply(resp, 1L, paste, sep='', collapse = '/')
    stringtabdata <- unique(stringfulldata)
    tabdata2 <- lapply(strsplit(stringtabdata, split='/'), as.integer)
    tabdata2 <- do.call(rbind, tabdata2)
    tabdata2[tabdata2 == 99999L] <- NA
    tabdata <- matrix(0L, nrow(tabdata2), sum(K))
    for(i in 1L:nitems){
        uniq <- sort(na.omit(unique(tabdata2[,i])))
        if(length(uniq) < K[i]) uniq <- 0L:(K[i]-1L)
        for(j in 1L:length(uniq))
            tabdata[,itemloc[i] + j - 1L] <- as.integer(tabdata2[,i] == uniq[j])
    }
    collectL <- numeric(nrow(tabdata))
    collectgrad <- matrix(0, nrow(tabdata), length(DX))
    gTheta <- vector('list', ngroups)
    for(g in 1L:ngroups){
        PrepList[[g]]$tabdata <- tabdata
        gTheta[[g]] <- Theta
    }
    for(pat in 1L:nrow(tabdata)){
        for(g in 1L:ngroups){
            gtabdata <- PrepList[[g]]$tabdata[pat, , drop=FALSE]
            rlist <- Estep.mirt(pars=pars[[g]], tabdata=gtabdata, freq=1L, CUSTOM.IND=CUSTOM.IND, full=full,
                                Theta=Theta, prior=Prior[[g]], itemloc=itemloc, deriv=TRUE)
            for(i in 1L:nitems){
                tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                pars[[g]][[i]]@dat <- rlist$r1[, tmp]
                pars[[g]][[i]]@itemtrace <- rlist$itemtrace[, tmp]
            }
            pars[[g]][[nitems + 1L]]@rr <- rowSums(rlist$r1)
        }
        DX <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, nitems), ncol(L), 0L, 0L, 1L, TRUE)$grad
        collectL[pat] <- rlist$expected
        DX[DX != 0] <-  rlist$expected - exp(log(rlist$expected) - DX[DX != 0])
        collectgrad[pat, ] <- DX
    }
    info <- N * t(collectgrad) %*% diag(1/collectL) %*% collectgrad
    info <- updateHess(info, L=Ls$L)[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    colnames(info) <- rownames(info) <- names(ESTIMATE$correction)
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain, warn=warn)
    return(ESTIMATE)
}
