SE.BL <- function(pars, Theta, theta, prior, BFACTOR, itemloc, PrepList, ESTIMATE, constrain, Ls,
                  specific=NULL, sitems=NULL, EH = FALSE, EHPrior = NULL){
    LL <- function(p, est, longpars, pars, ngroups, J, Theta, PrepList, specific, sitems,
                   NO.CUSTOM, EH, EHPrior){
        longpars[est] <- p
        pars2 <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        gstructgrouppars <- prior <- Prior <- vector('list', ngroups)
        if(EH){
            Prior[[1L]] <- EHPrior[[1L]]
        } else {
            for(g in 1L:ngroups){
                gstructgrouppars[[g]] <- ExtractGroupPars(pars2[[g]][[J+1L]])
                if(BFACTOR){
                    prior[[g]] <- dnorm(theta, 0, 1)
                    prior[[g]] <- prior[[g]]/sum(prior[[g]])
                    Prior[[g]] <- apply(expand.grid(prior[[g]], prior[[g]]), 1L, prod)
                    next
                }
                Prior[[g]] <- mvtnorm::dmvnorm(Theta,gstructgrouppars[[g]]$gmeans,
                                               gstructgrouppars[[g]]$gcov)
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
            }
        }
        LL <- 0
        for(g in 1L:ngroups){
            if(BFACTOR){
                expected <- Estep.bfactor(pars=pars2[[g]], tabdata=PrepList[[g]]$tabdata,
                                          Theta=Theta, prior=prior[[g]],
                                          specific=specific, sitems=sitems,
                                          itemloc=itemloc, NO.CUSTOM=NO.CUSTOM)$expected
            } else {
                expected <- Estep.mirt(pars=pars2[[g]], tabdata=PrepList[[g]]$tabdata,
                                       Theta=Theta, prior=Prior[[g]], itemloc=itemloc,
                                       NO.CUSTOM=NO.CUSTOM)$expected
            }
            LL <- LL + sum(PrepList[[g]]$tabdata[,ncol(PrepList[[g]]$tabdata)] * log(expected))
        }
        LL
    }
    
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
    NO.CUSTOM <- !any(sapply(pars, class) %in% 'custom')
    hess <- numDeriv::hessian(LL, x=shortpars, est=est, longpars=longpars,
                              pars=pars, ngroups=ngroups, J=J,
                              Theta=Theta, PrepList=PrepList,
                              specific=specific, sitems=sitems, NO.CUSTOM=NO.CUSTOM,
                              EH=EH, EHPrior=EHPrior)
    Hess <- matrix(0, length(longpars), length(longpars))
    Hess[est, est] <- -hess
    Hess <- updateHess(h=Hess, L2=Ls$L2, L3=Ls$L3)
    info <- Hess[infological, infological]
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain)
    return(ESTIMATE)
}

SE.SEM <- function(est, pars, constrain, Ls, PrepList, list, Theta, theta, BFACTOR, ESTIMATE, DERIV,
                   collectLL, from, to, is.latent){
    TOL <- list$TOL
    itemloc <- list$itemloc
    J <- length(itemloc) - 1L
    L <- Ls$L
    MSTEPTOL <- list$MSTEPTOL
    sitems <- list$sitems
    specific <- list$specific
    ngroups <- ESTIMATE$ngroups
    NCYCLES <- ESTIMATE$cycles
    if(NCYCLES <= 5L) stop('SEM can not be computed due to short EM history')
    BFACTOR <- list$BFACTOR
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
    gTheta <- vector('list', ngroups)
    ANY.PRIOR <- rep(FALSE, ngroups)
    NO.CUSTOM <- !any(sapply(pars, class) %in% 'custom')
    converged <- logical(sum(estpars & !redun_constr))
    if(!is.latent[est]) converged[is.latent] <- TRUE
    rijfull <- rep(NA, length(converged))
    for(g in 1L:ngroups){
        gTheta[[g]] <- Theta
        if(length(prodlist) > 0L)
            gTheta[[g]] <- prodterms(gTheta[[g]],prodlist)
        if(BFACTOR){
            Thetabetween <- thetaComb(theta=theta, nfact=nfact-ncol(sitems))
            prior[[g]] <- dnorm(theta, 0, 1)
            prior[[g]] <- prior[[g]]/sum(prior[[g]])
        }
        ANY.PRIOR[g] <- any(sapply(pars[[g]], function(x) x@any.prior))
    }
    groupest <- FALSE
    for(g in 1L:ngroups)
        groupest <- any(groupest, pars[[g]][[J+1]]@est)
    
    for (cycles in from:to){
        
        longpars <- MLestimates
        longpars[estindex] <- EMhistory[cycles, estindex]
        if(length(constrain) > 0L)
            for(i in 1L:length(constrain))
                longpars[constrain[[i]][-1L]] <- longpars[[constrain[[i]][1L]]]
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        tmp <- updatePrior(pars=pars, gTheta=gTheta, Thetabetween=Thetabetween,
                           list=list, ngroups=ngroups, nfact=nfact, prior=prior,
                           J=J, BFACTOR=BFACTOR, sitems=sitems, cycles=cycles, rlist=rlist)
        Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
        #Estep
        for(g in 1L:ngroups){
            if(BFACTOR){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=PrepList[[g]]$tabdata,
                                            Theta=gTheta[[g]], prior=prior[[g]], Prior=Prior[[g]],
                                            Priorbetween=Priorbetween[[g]], specific=specific, sitems=sitems,
                                            itemloc=itemloc, NO.CUSTOM=NO.CUSTOM)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=PrepList[[g]]$tabdata, NO.CUSTOM=NO.CUSTOM,
                                         Theta=gTheta[[g]], prior=Prior[[g]], itemloc=itemloc)
            }
        }
        for(g in 1L:ngroups){
            for(i in 1L:J){
                tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                pars[[g]][[i]]@dat <- rlist[[g]]$r1[, tmp]
            }
        }
        longpars <- Mstep(pars=pars, est=estpars, longpars=longpars, ngroups=ngroups, J=J,
                          gTheta=gTheta, itemloc=itemloc, Prior=Prior, ANY.PRIOR=ANY.PRIOR,
                          NO.CUSTOM=NO.CUSTOM, PrepList=PrepList, L=L, UBOUND=UBOUND, LBOUND=LBOUND,
                          rlist=rlist, constrain=constrain, cycle=cycles, DERIV=DERIV, groupest=groupest)
        rijlast <- rij
        denom <- (EMhistory[cycles, estindex] - MLestimates[estindex])
        rij <- (longpars[estpars & !redun_constr] - MLestimates[estpars & !redun_constr]) / denom
        diff <- abs(rij - rijlast) < TOL
        converged <- diff | converged
        which <- is.na(rijfull) & converged
        rijfull[which] <- rij[which]
        if(all(!is.na(rijfull))) break
    } #END EM
    rijfull[is.na(rijfull)] <- rij[is.na(rijfull)]
    return(rijfull)
}

SE.simple <- function(PrepList, ESTIMATE, Theta, constrain, Ls, N, type){
    pars <- ESTIMATE$pars
    itemloc <- PrepList[[1L]]$itemloc
    ngroups <- length(pars)
    nitems <- length(pars[[1L]]) - 1L
    L <- ESTIMATE$L
    DX <- numeric(ncol(L))
    Prior <- ESTIMATE$Prior
    Igrad <- IgradP <- Ihess <- matrix(0, length(DX), length(DX))
    fast <- ifelse(type == 'crossprod', TRUE, FALSE)
    if(fast){
        for(pat in 1L:nrow(PrepList[[1L]]$tabdata)){
            for(g in 1L:ngroups){
                gtabdata <- PrepList[[g]]$tabdata[pat, , drop=FALSE]
                r <- gtabdata[,ncol(gtabdata)]
                gtabdata[,ncol(gtabdata)] <- 1L
                pick <- min(which(gtabdata == 1L))
                rlist <- Estep.mirt(pars=pars[[g]], tabdata=gtabdata,
                                    Theta=Theta, prior=Prior[[g]], itemloc=itemloc, deriv=TRUE)
                w <- rlist$r1[,pick]
                tmpderiv <- matrix(0, nrow(Theta), length(DX))
                for(i in 1L:nitems){
                    tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                    pars[[g]][[i]]@itemtrace <- rlist$itemtrace[, tmp]
                    pars[[g]][[i]]@dat <- matrix(gtabdata[, tmp, drop=FALSE], nrow(Theta), 
                                                 length(tmp), byrow = TRUE) * w
                    tmp <- Deriv(pars[[g]][[i]], Theta=Theta, EM=TRUE, estHess=FALSE)
                    DX[pars[[g]][[i]]@parnum] <- tmp$grad
                }
            }
            Igrad <- Igrad + updateHess(outer(DX, DX), L2=Ls$L2, L3=Ls$L3) * r
        }
    } else {
        for(pat in 1L:nrow(PrepList[[1L]]$tabdata)){
            for(g in 1L:ngroups){
                gtabdata <- PrepList[[g]]$tabdata[pat, , drop=FALSE]
                r <- gtabdata[,ncol(gtabdata)]
                gtabdata[,ncol(gtabdata)] <- 1L
                pick <- min(which(gtabdata == 1L))
                rlist <- Estep.mirt(pars=pars[[g]], tabdata=gtabdata,
                                    Theta=Theta, prior=Prior[[g]], itemloc=itemloc, deriv=TRUE)
                w <- rlist$r1[,pick]
                tmpderiv <- matrix(0, nrow(Theta), length(DX))
                tmphess <- matrix(0, nrow(Ihess), ncol(Ihess))
                for(i in 1L:nitems){
                    tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                    pars[[g]][[i]]@itemtrace <- rlist$itemtrace[, tmp]
                    pars[[g]][[i]]@dat <- gtabdata[, tmp, drop=FALSE]
                    for(j in 1L:nrow(Theta)){
                        tmp <- Deriv(pars[[g]][[i]], Theta=Theta[j, , drop=FALSE], EM=TRUE, estHess=TRUE)
                        tmpderiv[j,pars[[g]][[i]]@parnum] <- tmp$grad
                        tmphess[pars[[g]][[i]]@parnum, pars[[g]][[i]]@parnum] <- 
                            tmphess[pars[[g]][[i]]@parnum, pars[[g]][[i]]@parnum] + tmp$hess * w[j] * r
                    }
                }
            }
            Ihess <- Ihess + updateHess(tmphess, L2=Ls$L2, L3=Ls$L3)
            tmpderiv <- tmpderiv * w 
            for(j in 1L:nrow(Theta))
                IgradP <- IgradP + updateHess(outer(tmpderiv[j,], tmpderiv[j,]), L2=Ls$L2, L3=Ls$L3) * r
            DX <- colSums(tmpderiv)
            Igrad <- Igrad + updateHess(outer(DX, DX), L2=Ls$L2, L3=Ls$L3) * r
        }
    }
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
    lengthsplit <- do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'COV_'), length))
    lengthsplit <- lengthsplit + do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'MEAN_'), length))
    info[lengthsplit > 2L, lengthsplit > 2L] <- 1
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain)
    if(any(lengthsplit > 2L)){
        for(g in 1L:ngroups){
            tmp <- ESTIMATE$pars[[g]][[nitems+1L]]@SEpar
            tmp[!is.na(tmp)] <- NaN
            ESTIMATE$pars[[g]][[nitems+1L]]@SEpar <- tmp
        }
    }
    return(ESTIMATE)
}

SE.Fisher <- function(PrepList, ESTIMATE, Theta, constrain, Ls, N){
    pars <- ESTIMATE$pars
    itemloc <- PrepList[[1L]]$itemloc
    ngroups <- length(pars)
    nitems <- length(pars[[1L]]) - 1L
    L <- Ls$L
    DX <- numeric(ncol(L))
    Prior <- ESTIMATE$Prior
    Igrad <- Ihess <- matrix(0, length(DX), length(DX))
    tabdata <- PrepList[[1L]]$tabdata
    K <- PrepList[[1L]]$K
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
    tabdata <- cbind(tabdata, 1L)
    collectL <- numeric(nrow(tabdata))
    collectgrad <- matrix(0, nrow(tabdata), length(DX))
    for(g in 1L:ngroups)
        PrepList[[g]]$tabdata <- tabdata
    for(pat in 1L:nrow(tabdata)){
        for(g in 1L:ngroups){
            gtabdata <- PrepList[[g]]$tabdata[pat, , drop=FALSE]
            rlist <- Estep.mirt(pars=pars[[g]], tabdata=gtabdata,
                                Theta=Theta, prior=Prior[[g]], itemloc=itemloc, deriv=TRUE)
            for(i in 1L:nitems){
                tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                pars[[g]][[i]]@dat <- rlist$r1[, tmp]
                pars[[g]][[i]]@itemtrace <- rlist$itemtrace[, tmp]
                tmp <- Deriv(pars[[g]][[i]], Theta=Theta, EM = TRUE, estHess=FALSE)
                dx <- tmp$grad
                DX[pars[[g]][[i]]@parnum] <- dx
            }
        }        
        collectL[pat] <- rlist$expected
        DX <- as.numeric(L %*% DX)
        DX[DX != 0] <-  rlist$expected - exp(log(rlist$expected) - DX[DX != 0])
        collectgrad[pat, ] <- DX
    }
    Igrad <- Igrad[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    colnames(Igrad) <- rownames(Igrad) <- names(ESTIMATE$correction)
    collectgrad <- collectgrad[, ESTIMATE$estindex_unique]
    info <- N * t(collectgrad) %*% diag(1/collectL) %*% collectgrad
    colnames(info) <- rownames(info) <- names(ESTIMATE$correction)
    lengthsplit <- do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'COV_'), length))
    lengthsplit <- lengthsplit + do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'MEAN_'), length))
    info[lengthsplit > 2L, lengthsplit > 2L] <- 1
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain)
    if(any(lengthsplit > 2L)){
        for(g in 1L:ngroups){
            tmp <- ESTIMATE$pars[[g]][[nitems+1L]]@SEpar
            tmp[!is.na(tmp)] <- NaN
            ESTIMATE$pars[[g]][[nitems+1L]]@SEpar <- tmp
        }
    }
    return(ESTIMATE)
}
