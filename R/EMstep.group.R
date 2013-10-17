EM.group <- function(pars, constrain, PrepList, list, Theta, DERIV)
{
    verbose <- list$verbose
    nfact <- list$nfact
    NCYCLES <- list$NCYCLES
    TOL <- list$TOL
    MSTEPTOL <- list$MSTEPTOL
    BFACTOR <- list$BFACTOR
    itemloc <- list$itemloc
    ngroups <- length(pars)
    specific <- list$specific
    sitems <- list$sitems
    theta <- list$theta
    J <- length(itemloc) - 1L
    nfullpars <- 0L
    estpars <- c()
    prodlist <- PrepList[[1L]]$prodlist
    gfulldata <- gtheta0 <- gstructgrouppars <- vector('list', ngroups)
    for(g in 1L:ngroups){
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])
        gfulldata[[g]] <- PrepList[[g]]$fulldata
        gtheta0[[g]] <- matrix(0, nrow(gfulldata[[g]]), nfact)
        for(i in 1L:(J+1L)){
            nfullpars <- nfullpars + length(pars[[g]][[i]]@par)
            estpars <- c(estpars, pars[[g]][[i]]@est)
        }
    }
    listpars <- vector('list', ngroups)
    for(g in 1L:ngroups){
        listpars[[g]] <- list()
        for(i in 1L:(J + 1L)){
            listpars[[g]][[i]] <- pars[[g]][[i]]@par
        }
    }
    lastpars2 <- lastpars1 <- listpars
    index <- 1L:nfullpars
    longpars <- rep(NA,nfullpars)
    ind1 <- 1L
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1L
            longpars[ind1:ind2] <- pars[[g]][[i]]@par
            ind1 <- ind2 + 1L
        }
    }
    stagecycle <- 1L
    converge <- 1L    
    inverse_fail_count <- 1L
    ##
    L <- c()
    for(g in 1L:ngroups)
        for(i in 1L:(J+1L))
            L <- c(L, pars[[g]][[i]]@est)
    estindex <- index[estpars]
    L <- diag(as.numeric(L))
    redun_constr <- rep(FALSE, length(estpars))
    if(length(constrain) > 0L){
        for(i in 1L:length(constrain)){
            L[constrain[[i]], constrain[[i]]] <- 1L
            for(j in 2L:length(constrain[[i]]))
                redun_constr[constrain[[i]][j]] <- TRUE
        }
    }
    estindex_unique <- index[estpars & !redun_constr]
    if(any(diag(L)[!estpars] > 0L)){
        redindex <- index[!estpars]
        stop('Constraint applied to fixed parameter(s) ',
             paste(paste0(redindex[diag(L)[!estpars] > 0L]), ''), ' but should only be applied to
                 estimated parameters. Please fix!')
    }
    Prior <- prior <- Priorbetween <- gstructgrouppars <- rlist <- r <- list()
    #make sure constrained pars are equal
    tmp <- rowSums(L)
    tmp[tmp == 0] <- 1L
    tmp <- matrix(1/tmp, length(longpars), length(longpars), byrow = TRUE)
    tmp2 <- abs(diag(L) - 1L)
    longpars <- diag((tmp * L) * longpars) + tmp2 * longpars
    LL <- 0
    for(g in 1L:ngroups)
        r[[g]] <- PrepList[[g]]$tabdata[, ncol(PrepList[[g]]$tabdata)]
    LBOUND <- UBOUND <- c()
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            LBOUND <- c(LBOUND, pars[[g]][[i]]@lbound)
            UBOUND <- c(UBOUND, pars[[g]][[i]]@ubound)
        }
    }
    est <- c()
    for(g in 1L:ngroups){
       for(j in 1L:(J+1L)){
           if(!is(pars[[g]][[j]], 'GroupPars')) est <- c(est, pars[[g]][[j]]@est)
           else est <- c(est, rep(FALSE, length(pars[[g]][[j]]@est)))
       }
    }
    if(length(constrain) > 0L)
       for(i in 1L:length(constrain))
           est[constrain[[i]][-1L]] <- FALSE
    names(longpars) <- names(est)
    EMhistory <- matrix(NA, NCYCLES+1L, length(longpars))
    EMhistory[1L,] <- longpars    
    gTheta <- vector('list', ngroups)
    ANY.PRIOR <- rep(FALSE, ngroups)
    NO.CUSTOM <- !any(sapply(pars[[1L]], class) %in% 'custom')
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
    preMstep.longpars2 <- preMstep.longpars <- longpars
    accel <- 0
    
    #EM
    for (cycles in 1L:NCYCLES){
        #priors
        for(g in 1L:ngroups){            
            gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1L]])            
            if(BFACTOR){
                sel <- 1L:(nfact-ncol(sitems))
                Priorbetween[[g]] <- mvtnorm::dmvnorm(Thetabetween,
                                        gstructgrouppars[[g]]$gmeans[sel],
                                        gstructgrouppars[[g]]$gcov[sel,sel,drop=FALSE])                
                Priorbetween[[g]] <- Priorbetween[[g]]/sum(Priorbetween[[g]]) 
                Prior[[g]] <- apply(expand.grid(Priorbetween[[g]], prior[[g]]), 1, prod)                
                next
            }
            Prior[[g]] <- mvtnorm::dmvnorm(gTheta[[g]][ ,1L:nfact,drop=FALSE],gstructgrouppars[[g]]$gmeans,
                                           gstructgrouppars[[g]]$gcov)
            Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])                                              
        }
        #Estep
        lastLL <- LL
        LL <- 0
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
            LL <- LL + sum(r[[g]]*log(rlist[[g]]$expected))
        }
        for(g in 1L:ngroups){
            for(i in 1L:J){
                tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
                pars[[g]][[i]]@rs <- rlist[[g]]$r1[, tmp]
            }
        }
        preMstep.longpars2 <- preMstep.longpars
        preMstep.longpars <- longpars
        if(all(!est)) break
        longpars <- Mstep(pars=pars, est=est, longpars=longpars, ngroups=ngroups, J=J,
                          gTheta=gTheta, itemloc=itemloc, Prior=Prior, ANY.PRIOR=ANY.PRIOR, 
                          NO.CUSTOM=NO.CUSTOM, PrepList=PrepList, L=L, UBOUND=UBOUND, LBOUND=LBOUND,
                          rlist=rlist, constrain=constrain, cycle=cycles, DERIV=DERIV)                   
        if(list$accelerate && cycles > 10L && cycles %% 3 == 0L){            
            dX2 <- preMstep.longpars - preMstep.longpars2
            dX <- longpars - preMstep.longpars
            d2X2 <- dX - dX2
            accel <- 1 - sqrt((dX2 %*% dX2) / (d2X2 %*% d2X2))
            if(accel < -4) accel <- -4
            longpars <- (1 - accel) * longpars + accel * preMstep.longpars
        }        
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        EMhistory[cycles+1L,] <- longpars
        if(verbose)
            if(cycles > 1L)
                cat(sprintf('\rIteration: %d, Log-Lik: %.3f, Max-Change: %.5f',
                            cycles, LL, max(abs(preMstep.longpars - longpars))))
        if(cycles > 3L && all(abs(preMstep.longpars - longpars) < TOL))  break
    } #END EM
    if(cycles == NCYCLES){
        message('EM iterations terminated after ', cycles, ' iterations.')        
        converge <- 0L
    }
    infological <- estpars & !redun_constr
    correction <- numeric(length(estpars[estpars & !redun_constr]))
    names(correction) <- names(estpars[estpars & !redun_constr])
    hess <- matrix(0)
    if(list$SEM){
        h <- matrix(0, nfullpars, nfullpars)
        ind1 <- 1L
        for(group in 1L:ngroups){
            for (i in 1L:J){
                deriv <- Deriv(x=pars[[group]][[i]], Theta=gTheta[[g]], EM = TRUE, 
                               estHess=TRUE)
                ind2 <- ind1 + length(deriv$grad) - 1L
                h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
                ind1 <- ind2 + 1L
            }
            i <- i + 1L
            deriv <- Deriv(x=pars[[group]][[i]], Theta=gTheta[[g]], EM = TRUE,
                           pars=pars[[group]], tabdata=PrepList[[group]]$tabdata,
                           itemloc=itemloc, estHess=TRUE)
            ind2 <- ind1 + length(deriv$grad) - 1L
            h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
            ind1 <- ind2 + 1L
        }
        hess <- L %*% h %*% L
        hess <- hess[estpars & !redun_constr, estpars & !redun_constr]
        #return internal functions for SEM.SE
        return(list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                    logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological,
                    estindex_unique=estindex_unique, correction=correction, hess=hess,
                    estpars=estpars & !redun_constr, redun_constr=redun_constr, ngroups=ngroups,
                    LBOUND=LBOUND, UBOUND=UBOUND, EMhistory=na.omit(EMhistory), random=list()))
    }
    ret <- list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological,
                estindex_unique=estindex_unique, correction=correction, hess=hess, random=list())
    ret
}

# Estep for mirt
Estep.mirt <- function(pars, tabdata, Theta, prior, itemloc, NO.CUSTOM=FALSE, 
                       itemtrace=NULL, deriv = FALSE)
{
    nquad <- nrow(Theta)
    J <- length(itemloc) - 1L
    r <- tabdata[ ,ncol(tabdata)]
    X <- tabdata[ ,1L:(ncol(tabdata) - 1L), drop = FALSE]
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, NO.CUSTOM=NO.CUSTOM)
    retlist <- .Call("Estep", itemtrace, prior, X, r)
    if(deriv) retlist$itemtrace <- itemtrace
    return(retlist)
}

# Estep for bfactor
Estep.bfactor <- function(pars, tabdata, Theta, prior, Prior, Priorbetween, specific, 
                          NO.CUSTOM=FALSE, sitems, itemloc, itemtrace=NULL)
{
    J <- length(itemloc) - 1L
    r <- tabdata[ ,ncol(tabdata)]
    X <- tabdata[ ,1L:(ncol(tabdata) - 1L)]
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, NO.CUSTOM=NO.CUSTOM)
    retlist <- .Call("Estepbfactor", itemtrace, prior, Priorbetween, X, r, sitems)
    r1 <- matrix(0, nrow(Theta), ncol(X))
    for (i in 1L:J){
        r1[ ,itemloc[i]:(itemloc[i+1L]-1L)] <-
            retlist$r1[ ,itemloc[i]:(itemloc[i+1L]-1L) + (specific[i] - 1L)*ncol(X) ]
    }
    r1 <- r1 * Prior 
    return(list(r1=r1, expected=retlist$expected))
}

Mstep <- function(pars, est, longpars, ngroups, J, gTheta, itemloc, PrepList, L, ANY.PRIOR,
                  UBOUND, LBOUND, constrain, cycle, DERIV, Prior, rlist, NO.CUSTOM){
    p <- longpars[est]
    opt <- try(optim(p, fn=Mstep.LL, gr=Mstep.grad, method='L-BFGS-B', 
                     control=list(maxit=ifelse(cycle > 10L, 10L, 5L), fnscale = -1L),  
                     DERIV=DERIV, rlist=rlist, NO.CUSTOM=NO.CUSTOM,
                     est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, gTheta=gTheta,
                     PrepList=PrepList, L=L, constrain=constrain, ANY.PRIOR=ANY.PRIOR,
                     UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc, lower=LBOUND[est], upper=UBOUND[est]),
            silent=TRUE)
    if(is(opt, 'try-error'))
        stop(opt)
    longpars[est] <- opt$par    
    if(length(constrain) > 0L)
        for(i in 1L:length(constrain))
            longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
    i = J + 1L
    for(group in 1L:ngroups){
        if(any(pars[[group]][[i]]@est)){
            newpars <- Deriv(x=pars[[group]][[i]], Theta=gTheta[[group]], EM = TRUE,
                           pars=pars[[group]], tabdata=PrepList[[group]]$tabdata,
                           itemloc=itemloc, prior=Prior[[group]])
            longpars[pars[[group]][[i]]@parnum[pars[[group]][[i]]@est]] <- newpars
        }
    }
    return(longpars)
}

Mstep.LL <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L, NO.CUSTOM,
                     constrain, LBOUND, UBOUND, itemloc, DERIV, rlist, ANY.PRIOR){
    longpars[est] <- p
    if(length(constrain) > 0L)
       for(i in 1L:length(constrain))
           longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    LLs <- numeric(ngroups)
    for(g in 1L:ngroups)
        LLs[g] <- LogLikMstep(pars[[g]], Theta=gTheta[[g]], rs=rlist[[g]], 
                              itemloc=itemloc, NO.CUSTOM=NO.CUSTOM, any.prior=ANY.PRIOR[g])        
    return(sum(LLs))
}

LogLikMstep <- function(x, Theta, itemloc, rs, any.prior, NO.CUSTOM){
    log_itemtrace <- log(computeItemtrace(pars=x, Theta=Theta, itemloc=itemloc, NO.CUSTOM=NO.CUSTOM))
    LL <- sum(rs$r1 * log_itemtrace)
    if(any.prior){
        for(i in 1L:(length(x)-1L))
            if(x[[i]]@any.prior)
                LL <- LL.Priors(x=x[[i]], LL=LL)
    }
    return(LL)
}

Mstep.grad <- function(p, est, longpars, pars, ngroups, J, gTheta, PrepList, L,  ANY.PRIOR,
                       constrain, LBOUND, UBOUND, itemloc, DERIV, rlist, NO.CUSTOM){
    longpars[est] <- p
    if(length(constrain) > 0L)
        for(i in 1L:length(constrain))
            longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    g <- rep(0, ncol(L))
    ind1 <- 1L
    for(group in 1L:ngroups){
        for (i in 1L:J){
            deriv <- DERIV[[group]][[i]](x=pars[[group]][[i]], Theta=gTheta[[group]], EM=TRUE)
            ind2 <- ind1 + length(deriv$grad) - 1L
            g[ind1:ind2] <- deriv$grad
            ind1 <- ind2 + 1L
        }
        i <- i + 1L
        ind2 <- ind1 + length(pars[[group]][[i]]@par) - 1L
        g[ind1:ind2] <- 0
        ind1 <- ind2 + 1L
    }
    grad <- g %*% L
    return(grad[est])
}
