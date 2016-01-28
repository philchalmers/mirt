SE.Numerical <- function(pars, Theta, theta, BFACTOR, itemloc, PrepList, ESTIMATE, constrain, Ls,
                  CUSTOM.IND, specific=NULL, sitems=NULL, EH = FALSE, EHPrior = NULL, warn, Data, type,
                  delta){
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
    if(type == 'Richardson'){
        hess <- numDeriv::hessian(BL.LL, x=shortpars, est=est, longpars=longpars,
                                  pars=pars, ngroups=ngroups, J=J, itemloc=itemloc,
                                  Theta=Theta, PrepList=PrepList, BFACTOR=BFACTOR,
                                  specific=specific, sitems=sitems, CUSTOM.IND=CUSTOM.IND,
                                  EH=EH, EHPrior=EHPrior, Data=Data, theta=theta)
    } else {
        hess <- numerical_deriv(shortpars, BL.LL, est=est, longpars=longpars,
                                pars=pars, ngroups=ngroups, J=J, itemloc=itemloc,
                                Theta=Theta, PrepList=PrepList, BFACTOR=BFACTOR,
                                specific=specific, sitems=sitems, CUSTOM.IND=CUSTOM.IND,
                                EH=EH, EHPrior=EHPrior, Data=Data, theta=theta, type=type,
                                h = delta, gradient = FALSE)
    }
    Hess <- matrix(0, length(longpars), length(longpars))
    Hess[est, est] <- -hess
    Hess <- updateHess(h=Hess, L=Ls$L)
    info <- Hess[infological, infological]
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain, warn=warn)
    return(ESTIMATE)
}

SE.SEM <- function(index, estmat, pars, constrain, Ls, PrepList, list, Theta, theta, BFACTOR, ESTIMATE, DERIV,
                   collectLL, from, to, is.latent, Data, solnp_args, control){
    est <- estmat[,index]
    lrPars <- list$lrPars
    full <- list$full
    TOL <- list$TOL
    itemloc <- list$itemloc
    J <- length(itemloc) - 1L
    L <- Ls$L
    Moptim <- list$Moptim
    sitems <- list$sitems
    specific <- list$specific
    ngroups <- ESTIMATE$ngroups
    NCYCLES <- ESTIMATE$cycles
    if(NCYCLES <= 5L) stop('SEM can not be computed due to short EM history', call.=FALSE)
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
    ANY.PRIOR <- rep(FALSE, ngroups)
    converged <- logical(sum(estpars & !redun_constr))
    rijfull <- rep(NA, length(converged))
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta, prodlist)
    if(BFACTOR){
        Thetabetween <- thetaComb(theta=theta, nfact=nfact-ncol(sitems))
        for(g in 1L:ngroups){
            prior[[g]] <- dnorm(theta, 0, 1)
            prior[[g]] <- prior[[g]]/sum(prior[[g]])
        }
    }
    gTheta <- vector('list', ngroups)
    for(g in 1L:ngroups){
        ANY.PRIOR[g] <- any(sapply(pars[[g]], function(x) x@any.prior))
        gTheta[[g]] <- Theta
    }
    converged <- FALSE

    for (cycles in from:to){

        longpars <- MLestimates
        longpars[estindex] <- EMhistory[cycles, estindex]
        if(length(constrain) > 0L)
            for(i in 1L:length(constrain))
                longpars[constrain[[i]][-1L]] <- longpars[[constrain[[i]][1L]]]
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        tmp <- updatePrior(pars=pars, Theta=Theta, Thetabetween=Thetabetween,
                           list=list, ngroups=ngroups, nfact=nfact, prior=prior, lrPars=lrPars,
                           J=J, BFACTOR=BFACTOR, sitems=sitems, cycles=cycles, rlist=rlist,
                           full=full)
        Prior <- tmp$Prior; Priorbetween <- tmp$Priorbetween
        #Estep
        for(g in 1L:ngroups){
            if(BFACTOR){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                            Theta=Theta, prior=prior[[g]], Prior=Prior[[g]],
                                            Priorbetween=Priorbetween[[g]], specific=specific, sitems=sitems,
                                            itemloc=itemloc, CUSTOM.IND=list$CUSTOM.IND)
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                         CUSTOM.IND=list$CUSTOM.IND, Theta=Theta,
                                         prior=Prior[[g]], itemloc=itemloc, full=full)
            }
        }
        for(g in 1L:ngroups){
            for(i in 1L:J){
                tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                pars[[g]][[i]]@dat <- rlist[[g]]$r1[, tmp]
            }
        }
        longpars <- Mstep(pars=pars, est=estpars & !ESTIMATE$groupest, longpars=longpars, ngroups=ngroups, J=J,
                          gTheta=gTheta, itemloc=itemloc, Prior=Prior, ANY.PRIOR=ANY.PRIOR,
                          PrepList=PrepList, L=L, UBOUND=UBOUND, LBOUND=LBOUND, nfact=nfact,
                          rlist=rlist, constrain=constrain, DERIV=DERIV, groupest=ESTIMATE$groupest,
                          CUSTOM.IND=list$CUSTOM.IND, SLOW.IND=list$SLOW.IND, BFACTOR=list$BFACTOR,
                          Moptim=Moptim, Mrate=1, TOL=list$MSTEPTOL, solnp_args=solnp_args, full=full,
                          Thetabetween=Thetabetween, lrPars=lrPars, control=control)
        rijlast <- rij
        denom <- (EMhistory[cycles, estindex] - MLestimates[estindex])
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
                      CUSTOM.IND, SLOW.IND, warn, message, Data){

    fn <- function(which, PrepList, ngroups, pars, Theta, Prior, itemloc, Igrad, Igrad2, Ihess,
                   CUSTOM.IND, SLOW.IND, whichitems, iscross, npars, Data){
        Igrad <- matrix(0, npars, npars)
        gtabdatafull <- Data$tabdatalong
        J <- Data$nitems
        if(iscross){
            for(pat in which){
                for(g in 1L:ngroups){
                    gtabdata <- gtabdatafull[pat, , drop=FALSE]
                    r <- Data$Freq[[g]][pat]
                    rlist <- Estep.mirt(pars=pars[[g]], tabdata=gtabdata, CUSTOM.IND=CUSTOM.IND, full=FALSE,
                                        freq=1L, Theta=Theta, prior=Prior[[g]], itemloc=itemloc, deriv=TRUE)
                    for(i in whichitems){
                        tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                        pars[[g]][[i]]@itemtrace <- rlist$itemtrace[, tmp]
                        pars[[g]][[i]]@dat <- rlist[[1L]][,tmp]
                        tmp <- Deriv(pars[[g]][[i]], Theta=Theta, estHess=FALSE)
                        DX[pars[[g]][[i]]@parnum] <- tmp$grad
                    }
                    if(any(pars[[g]][[J+1L]]@est)){
                       out <- Deriv(pars[[g]][[J+1L]], Theta=Theta, CUSTOM.IND=list(), EM = TRUE,
                                    pars = pars[[g]], itemloc = itemloc, tabdata = gtabdata, freq = 1L,
                                    estHess=FALSE, prior = Prior[[g]])
                       DX[pars[[g]][[J+1L]]@parnum] <- out$grad
                   }
                }
                Igrad <- Igrad + outer(DX, DX) * r
            }
            return(list(Igrad=Igrad, Ihess=matrix(0, ncol(Igrad), ncol(Igrad)),
                        IgradP=matrix(0, ncol(Igrad), ncol(Igrad))))
        } else {
            IgradP <- Ihess <- matrix(0, npars, npars)
            for(pat in 1L:nrow(gtabdatafull)){
                for(g in 1L:ngroups){
                    gtabdata <- gtabdatafull[pat, , drop=FALSE]
                    r <- Data$Freq[[g]][pat]
                    pick <- min(which(gtabdata == 1L))
                    rlist <- Estep.mirt(pars=pars[[g]], tabdata=gtabdata, CUSTOM.IND=CUSTOM.IND, full=FALSE,
                                        freq=1L, Theta=Theta, prior=Prior[[g]], itemloc=itemloc, deriv=TRUE)
                    w <- rlist$r1[,pick]
                    tmpderiv <- matrix(0, nrow(Theta), length(DX))
                    tmphess <- matrix(0, nrow(Ihess), ncol(Ihess))
                    for(i in whichitems){
                        tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                        pars[[g]][[i]]@itemtrace <- rlist$itemtrace[, tmp]
                        pars[[g]][[i]]@dat <- gtabdata[, tmp, drop=FALSE]
                        for(j in 1L:nrow(Theta)){
                            tmp <- Deriv(pars[[g]][[i]], Theta=Theta[j, , drop=FALSE], estHess=TRUE)
                            tmpderiv[j,pars[[g]][[i]]@parnum] <- tmp$grad
                            tmphess[pars[[g]][[i]]@parnum, pars[[g]][[i]]@parnum] <-
                                tmphess[pars[[g]][[i]]@parnum, pars[[g]][[i]]@parnum] + tmp$hess * w[j] * r
                        }
                    }
                }
                Ihess <- Ihess + tmphess
                tmpderiv2 <- tmpderiv * sqrt(w)
                for(j in 1L:nrow(Theta))
                    IgradP <- IgradP + outer(tmpderiv2[j,], tmpderiv2[j,]) * r
                DX <- colSums(tmpderiv * w)
                Igrad <- Igrad + outer(DX, DX) * r
            }
        }
        return(list(Igrad=Igrad, Ihess=Ihess, IgradP=IgradP))
    }

    pars <- ESTIMATE$pars
    itemloc <- PrepList[[1L]]$itemloc
    ngroups <- length(pars)
    nitems <- length(pars[[1L]]) - 1L
    L <- ESTIMATE$L
    DX <- numeric(ncol(L))
    Prior <- ESTIMATE$Prior
    prior <- ESTIMATE$prior
    Priorbetween <- ESTIMATE$Priorbetween
    isbifactor <- length(Priorbetween[[1L]]) > 0L
    if(!isbifactor)
        prior <- Priorbetween <- list(matrix(0))
    sitems <- ESTIMATE$sitems
    iscross <- ifelse(type == 'crossprod', TRUE, FALSE)
    gitemtrace <- rs <- vector('list', ngroups)
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
#     infolist <- fn(1L:ncol(rs), PrepList, ngroups, pars, Theta, Prior, itemloc, Igrad, Igrad2, Ihess,
#                     CUSTOM.IND, SLOW.IND, whichitems=1:length(Data$K), iscross, npars, Data)
    whichitems <- unique(c(CUSTOM.IND, SLOW.IND))
    infolist <- .Call("computeInfo", pars, Theta, gPrior, prior[[1L]], Priorbetween[[1L]],
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
    lengthsplit <- do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'COV_'), length))
    lengthsplit <- lengthsplit + do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'MEAN_'), length))
    if(!iscross){
        is.latent <- lengthsplit > 2L
        Ihess <- Ihess[!is.latent, !is.latent]; Igrad <- Igrad[!is.latent, !is.latent]
        IgradP <- IgradP[!is.latent, !is.latent]
    } else is.latent <- logical(length(lengthsplit))
    if(type == 'Louis'){
        info <- -Ihess - IgradP + Igrad
    } else if(type == 'crossprod'){
        info <- Igrad
    } else if(type == 'sandwich'){
        tmp <- solve(-Ihess - IgradP + Igrad)
        info <- solve(tmp %*% Igrad %*% tmp)
    }
    colnames(info) <- rownames(info) <- names(ESTIMATE$correction)[!is.latent]
    tmp <- matrix(NA, length(is.latent), length(is.latent))
    tmp[!is.latent, !is.latent] <- info
    ESTIMATE <- loadESTIMATEinfo(info=tmp, ESTIMATE=ESTIMATE, constrain=constrain, warn=warn)
    if(!iscross && any(lengthsplit > 2L)){
        for(g in 1L:ngroups){
            tmp <- ESTIMATE$pars[[g]][[nitems+1L]]@SEpar
            tmp[!is.na(tmp)] <- NaN
            ESTIMATE$pars[[g]][[nitems+1L]]@SEpar <- tmp
        }
    }
    return(ESTIMATE)
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
    Igrad <- Ihess <- matrix(0, length(DX), length(DX))
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
    for(g in 1L:ngroups)
        PrepList[[g]]$tabdata <- tabdata
    for(pat in 1L:nrow(tabdata)){
        for(g in 1L:ngroups){
            gtabdata <- PrepList[[g]]$tabdata[pat, , drop=FALSE]
            rlist <- Estep.mirt(pars=pars[[g]], tabdata=gtabdata, freq=1L, CUSTOM.IND=CUSTOM.IND, full=full,
                                Theta=Theta, prior=Prior[[g]], itemloc=itemloc, deriv=TRUE)
            for(i in 1L:nitems){
                tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                pars[[g]][[i]]@dat <- rlist$r1[, tmp]
                pars[[g]][[i]]@itemtrace <- rlist$itemtrace[, tmp]
                tmp <- Deriv(pars[[g]][[i]], Theta=Theta, estHess=FALSE)
                dx <- tmp$grad
                DX[pars[[g]][[i]]@parnum] <- dx
            }
        }
        collectL[pat] <- rlist$expected
        DX[DX != 0] <-  rlist$expected - exp(log(rlist$expected) - DX[DX != 0])
        collectgrad[pat, ] <- DX
    }
    info <- N * t(collectgrad) %*% diag(1/collectL) %*% collectgrad
    info <- updateHess(info, L=Ls$L)[ESTIMATE$estindex_unique, ESTIMATE$estindex_unique]
    colnames(info) <- rownames(info) <- names(ESTIMATE$correction)
    lengthsplit <- do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'COV_'), length))
    lengthsplit <- lengthsplit + do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'MEAN_'), length))
    info[lengthsplit > 2L, lengthsplit > 2L] <- NA
    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain, warn=warn)
    if(any(lengthsplit > 2L)){
        for(g in 1L:ngroups){
            tmp <- ESTIMATE$pars[[g]][[nitems+1L]]@SEpar
            tmp[!is.na(tmp)] <- NaN
            ESTIMATE$pars[[g]][[nitems+1L]]@SEpar <- tmp
        }
    }
    return(ESTIMATE)
}
