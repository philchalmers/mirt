SE.Numerical <- function(pars, Theta, theta, dentype, itemloc, PrepList, ESTIMATE, constrain, Ls,
                  CUSTOM.IND, specific=NULL, sitems=NULL, EHPrior = NULL, warn, Data, type,
                  delta, lrPars, omp_threads){
    longpars <- ESTIMATE$longpars
    rlist <- ESTIMATE$rlist
    infological=ESTIMATE$infological
    ngroups <- length(pars)
    J <- length(pars[[1L]]) - 1L
    est <- c()
    for(g in seq_len(ngroups))
        for(j in seq_len(J+1L))
            est <- c(est, pars[[g]][[j]]@est)
    shortpars <- longpars[est]
    gstructgrouppars <- vector('list', ngroups)
    for(g in seq_len(ngroups))
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
    for(g in seq_len(ngroups)){
        for(i in seq_len(J)){
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
                            delta = delta, omp_threads=omp_threads, gradient = FALSE)
    Hess <- matrix(0, length(longpars), length(longpars))
    Hess[est, est] <- -hess
    Hess <- updateHess(h=Hess, L=Ls$L)
    info <- as.matrix(Hess[infological, infological])
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
    if(length(prodlist))
        Theta <- prodterms(Theta, prodlist)
    gTheta <- vector('list', ngroups)
    for(g in seq_len(ngroups)){
        ANY.PRIOR[g] <- any(sapply(pars[[g]], function(x) x@any.prior))
        gTheta[[g]] <- Theta
    }
    converged <- FALSE

    for (cycles in from:to){

        longpars <- MLestimates
        longpars[estindex] <- EMhistory[cycles, estindex]
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        tmp <- updatePrior(pars=pars, gTheta=gTheta, list=list, ngroups=ngroups,
                           nfact=nfact, lrPars=lrPars, J=J, dentype=dentype,
                           sitems=sitems, cycles=cycles, rlist=rlist, full=full, MC=list$method == 'QMC')
        prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
        #Estep
        for(g in seq_len(ngroups)){
            if(dentype == 'bfactor'){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                            Theta=Theta, prior=prior[[g]],
                                            Priorbetween=Priorbetween[[g]], specific=specific, sitems=sitems,
                                            itemloc=itemloc, CUSTOM.IND=list$CUSTOM.IND, omp_threads=1L)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                         CUSTOM.IND=list$CUSTOM.IND, Theta=Theta,
                                         prior=Prior[[g]], itemloc=itemloc, full=full, omp_threads=1L)
            }
        }
        for(g in seq_len(ngroups)){
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
    iscross <- ifelse(type %in% c('crossprod', 'sandwich'), TRUE, FALSE)
    gitemtrace <- vector('list', ngroups)
    for(g in seq_len(ngroups)){
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
    Igrad <- infolist[["Igrad"]]; IgradP <- infolist[["IgradP"]]
    if(length(whichitems)){
        warning('Internal information matrix computations currently not supported for at
        least one of the supplied items. Information matrix/standard errors not computed', call.=FALSE)
        return(ESTIMATE)
    }
    if(type != 'Louis'){
        Iprior <- matrix(0, nrow(Igrad), ncol(Igrad))
        ind1 <- 1L
        for(group in seq_len(length(pars))){
            for (i in seq_len(length(pars[[group]]))){
                np <- length(pars[[group]][[i]]@par)
                deriv <- DerivativePriors(x=pars[[group]][[i]], grad=numeric(np),
                                          hess=matrix(0, np, np))
                ind2 <- ind1 + np - 1L
                Iprior[ind1:ind2, ind1:ind2] <- deriv$hess
                ind1 <- ind2 + 1L
            }
        }
        Iprior <- -updateHess(Iprior, L=Ls$L)
        Iprior <- Iprior[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    }
    Igrad <- updateHess(Igrad, L=Ls$L)
    IgradP <- updateHess(IgradP, L=Ls$L)
    Igrad <- Igrad[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    IgradP <- IgradP[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    Ihess <- complete
    if(type == 'Louis'){
        info <- -Ihess - IgradP + Igrad
    } else if(type == 'crossprod'){
        info <- Igrad + Iprior
    } else if(type == 'sandwich.Louis'){
        tmp <- solve(-Ihess - IgradP + Igrad)
        info <- solve(tmp %*% (Igrad + Iprior) %*% tmp)
    } else if(type == 'sandwich'){
        tmp <- -solve(ESTIMATE$hess)
        info <- solve(tmp %*% (Igrad + Iprior) %*% tmp)
    }
    info <- as.matrix(info)
    colnames(info) <- rownames(info) <- names(ESTIMATE$correction)
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain, warn=warn)
    return(ESTIMATE)
}

SE.Oakes <- function(pick, pars, L, constrain, est, shortpars, longpars,
                     gTheta, list, ngroups, J, dentype, sitems,
                     rlist, full, Data, specific, itemloc, CUSTOM.IND,
                     delta, prior, Prior, Priorbetween, nfact, mixtype,
                     PrepList, ANY.PRIOR, DERIV, SLOW.IND, Norder, omp_threads,
                     zero_g = NULL){
    r <- 1L
    Richardson <- if(Norder > 2L) TRUE else FALSE
    if(Richardson){
        delta <- delta * 10
        r <- Norder
        Norder <- 2L
        R <- array(0, dim = c(length(shortpars), r, 2L))
    }
    if(is.null(zero_g)){
        grad <- matrix(0, 2L, length(shortpars))
        signs <- c(-1, 1)
    } else {
        grad <- matrix(0, 2L, length(shortpars))
        grad[1L, ] <- zero_g
        signs <- 1
    }
    for(rr in 0L:(r-1L)){
        row <- ifelse(is.null(zero_g), 1L, 2L)
        for(sign in signs){
            longpars_old <- longpars
            d <- sign * delta
            longpars[which(est)[pick]] <- shortpars[pick] + d
            longpars <- longpars_constrain(longpars, constrain)
            pars <- reloadPars(longpars=longpars, pars=pars,
                               ngroups=ngroups, J=J)
            tmp <- updatePrior(pars=pars, gTheta=gTheta,
                               list=list, ngroups=ngroups, nfact=nfact,
                               J=J, dentype=dentype, sitems=sitems, cycles=100L,
                               rlist=rlist, full=full, MC=list$method == 'QMCEM')
            prior <- tmp$prior; Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
            Elist <- Estep(pars=pars, Data=Data, gTheta=gTheta, prior=prior, Prior=Prior,
                           Priorbetween=Priorbetween, specific=specific, sitems=sitems,
                           ngroups=ngroups, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                           dentype=dentype, rlist=rlist, full=full, Etable=list$Etable,
                           omp_threads=omp_threads)
            rlist <- Elist$rlist
            longpars <- longpars_old
            pars <- reloadPars(longpars=longpars, pars=pars,
                               ngroups=ngroups, J=J)
            if(pars[[1L]][[J + 1L]]@itemclass == -1L){
                for(g in seq_len(length(pars))){
                    gp <- pars[[g]][[J + 1L]]
                    pars[[g]][[J + 1L]]@density <- gp@safe_den(gp, gTheta[[g]])
                }
            }
            for(g in seq_len(ngroups)){
                for(i in 1L:J)
                    pars[[g]][[i]]@dat <- rlist[[g]]$r1[ , c(itemloc[i]:(itemloc[i+1L] - 1L)),
                                                         drop=FALSE]
                if(dentype == 'bfactor'){
                    pars[[g]][[J+1L]]@rrb <- rlist[[g]]$r2
                    pars[[g]][[J+1L]]@rrs <- rlist[[g]]$r3
                } else pars[[g]][[J+1L]]@rr <- rowSums(rlist[[g]]$r1) / J
            }
            browser()
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
            if(dentype == 'mixture'){
                mixtype$par <- sapply(1L:length(pars),
                                      function(g) sum(pars[[g]][[J+1L]]@par[length(pars[[g]][[J+1L]]@par)]))
                mixtype$dat <- matrix(sapply(1L:length(pars),
                                             function(g) sum(pars[[g]][[J+1L]]@rr)), 1L)
                deriv <- Deriv.mix(mixtype)
                g[mixtype$parnum] <- deriv$grad
            }
            tmp <- g %*% L
            if(pick == 0L) return(tmp[est])
            grad[row, ] <- tmp[est]
            row <- row + 1L
        }
        cfs <- get_deriv_coefs(order=Norder)
        ret <- colSums(cfs * grad) / delta
        if(!Richardson) return(ret)
        delta <- delta / 2
        if(rr == 0L){
            R[, 1L, 1L] <- ret
        } else {
            R[ ,1L, 2L] <- ret
            for (j in 1L:rr)
                R[ ,j + 1L, 2L] <- (4^j*R[ , j, 2L] - R[, j, 1L])/(4^j - 1)
            R[ , , 1L] <- R[ , , 2L]
        }
    }
    R[ , r, 2L]
}

SE.Fisher <- function(PrepList, ESTIMATE, Theta, constrain, Ls, CUSTOM.IND,
                      SLOW.IND, warn, Data, full, omp_threads){
    pars <- ESTIMATE$pars
    Ns <- table(Data$group)
    itemloc <- PrepList[[1L]]$itemloc
    ngroups <- length(pars)
    nitems <- length(pars[[1L]]) - 1L
    L <- Ls$L
    DX <- numeric(ncol(L))
    Prior <- ESTIMATE$Prior
    tabdata <- Data$tabdatalong
    K <- Data$K
    resp <- vector('list', nitems)
    for(i in seq_len(nitems))
        resp[[i]] <- 0L:(K[i]-1L)
    resp <- expand.grid(resp)
    stringfulldata <- apply(resp, 1L, paste, sep='', collapse = '/')
    stringtabdata <- unique(stringfulldata)
    tabdata2 <- lapply(strsplit(stringtabdata, split='/'), as.integer)
    tabdata2 <- do.call(rbind, tabdata2)
    tabdata2[tabdata2 == 99999L] <- NA
    tabdata <- matrix(0L, nrow(tabdata2), sum(K))
    for(i in seq_len(nitems)){
        uniq <- sort(na.omit(unique(tabdata2[,i])))
        if(length(uniq) < K[i]) uniq <- 0L:(K[i]-1L)
        for(j in seq_len(length(uniq)))
            tabdata[,itemloc[i] + j - 1L] <- as.integer(tabdata2[,i] == uniq[j])
    }
    gTheta <- vector('list', ngroups)
    for(g in seq_len(ngroups)){
        PrepList[[g]]$tabdata <- tabdata
        gTheta[[g]] <- Theta
    }
    info_list <- vector('list', ngroups)
    for(g in seq_len(ngroups)){
        info <- 0
        # reset
        for(gg in seq_len(ngroups)){
            for(i in seq_len(nitems)){
                tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                pars[[gg]][[i]]@dat <- pars[[gg]][[i]]@itemtrace <-
                    matrix(0, nrow(Theta), length(tmp))
            }
            pars[[gg]][[nitems + 1L]]@rr <-
                numeric(length(pars[[gg]][[nitems + 1L]]@rr))
        }
        for(pat in seq_len(nrow(tabdata))){
            gtabdata <- PrepList[[g]]$tabdata[pat, , drop=FALSE]
            browser()
            rlist <- Estep.mirt(pars=pars[[g]], tabdata=gtabdata, freq=1L,
                                CUSTOM.IND=CUSTOM.IND, full=full,
                                Theta=Theta, prior=Prior[[g]], itemloc=itemloc,
                                deriv=TRUE, omp_threads=omp_threads)
            for(i in seq_len(nitems)){
                tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                pars[[g]][[i]]@dat <- rlist$r1[, tmp]
                pars[[g]][[i]]@itemtrace <- rlist$itemtrace[, tmp]
            }
            pars[[g]][[nitems + 1L]]@rr <- rowSums(rlist$r1)
            DX <- .Call('computeDPars', pars, gTheta, matrix(0L, 1L, nitems),
                        ncol(L), 0L, 0L, 1L, TRUE)$grad
            info <- info + DX %*% t(DX) * rlist$expected
        }
        info_list[[g]] <- Ns[g] * info
    }
    info <- Reduce("+", info_list)
    info <- as.matrix(updateHess(info, L=Ls$L)[ESTIMATE$estindex_unique,
                                               ESTIMATE$estindex_unique])
    colnames(info) <- rownames(info) <- names(ESTIMATE$correction)
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE,
                                 constrain=constrain, warn=warn)
    return(ESTIMATE)
}
