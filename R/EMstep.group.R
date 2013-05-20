EM.group <- function(pars, constrain, PrepList, list, Theta)
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
    J <- length(itemloc) - 1
    nfullpars <- 0
    estpars <- c()
    gfulldata <- gtheta0 <- gstructgrouppars <- gitemtrace <- vector('list', ngroups)
    for(g in 1:ngroups){
        gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1]])
        gfulldata[[g]] <- PrepList[[g]]$fulldata
        gtheta0[[g]] <- matrix(0, nrow(gfulldata[[g]]), nfact)
        for(i in 1:(J+1)){
            nfullpars <- nfullpars + length(pars[[g]][[i]]@par)
            estpars <- c(estpars, pars[[g]][[i]]@est)
        }
    }
    listpars <- vector('list', ngroups)
    for(g in 1:ngroups){
        listpars[[g]] <- list()
        for(i in 1:(J + 1)){
            listpars[[g]][[i]] <- pars[[g]][[i]]@par
        }
    }
    lastpars2 <- lastpars1 <- listpars
    index <- 1:nfullpars
    longpars <- rep(NA,nfullpars)
    ind1 <- 1
    for(g in 1:ngroups){
        for(i in 1:(J+1)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1
            longpars[ind1:ind2] <- pars[[g]][[i]]@par
            ind1 <- ind2 + 1
        }
    }
    stagecycle <- 1
    converge <- 1
    LLwarn <- FALSE
    inverse_fail_count <- 1
    ##
    L <- c()
    for(g in 1:ngroups)
        for(i in 1:(J+1))
            L <- c(L, pars[[g]][[i]]@est)
    estindex <- index[estpars]
    L <- diag(as.numeric(L))
    redun_constr <- rep(FALSE, length(estpars))
    if(length(constrain) > 0){
        for(i in 1:length(constrain)){
            L[constrain[[i]], constrain[[i]]] <- 1
            for(j in 2:length(constrain[[i]]))
                redun_constr[constrain[[i]][j]] <- TRUE
        }
    }
    estindex_unique <- index[estpars & !redun_constr]
    if(any(diag(L)[!estpars] > 0)){
        redindex <- index[!estpars]
        stop('Constraint applied to fixed parameter(s) ',
             paste(redindex[diag(L)[!estpars] > 0]), ' but should only be applied to
                 estimated parameters. Please fix!')
    }
    Prior <- prior <- gstructgrouppars <- rlist <- r <- list()
    #make sure constrained pars are equal
    tmp <- rowSums(L)
    tmp[tmp == 0] <- 1
    tmp <- matrix(1/tmp, length(longpars), length(longpars), byrow = TRUE)
    tmp2 <- abs(diag(L) - 1)
    longpars <- diag((tmp * L) * longpars) + tmp2 * longpars
    LL <- 0
    for(g in 1:ngroups)
        r[[g]] <- PrepList[[g]]$tabdata[, ncol(PrepList[[g]]$tabdata)]
    LBOUND <- UBOUND <- c()
    for(g in 1:ngroups){
        for(i in 1:(J+1)){
            LBOUND <- c(LBOUND, pars[[g]][[i]]@lbound)
            UBOUND <- c(UBOUND, pars[[g]][[i]]@ubound)
        }
    }
    est <- c()
    for(g in 1:ngroups){
       for(j in 1:(J+1)){
           if(!is(pars[[g]][[j]], 'GroupPars')) est <- c(est, pars[[g]][[j]]@est)
           else est <- c(est, rep(FALSE, length(pars[[g]][[j]]@est)))
       }
    }
    if(length(constrain) > 0)
       for(i in 1:length(constrain))
           est[constrain[[i]][-1]] <- FALSE
    EMhistory <- matrix(NA, NCYCLES+1, length(longpars))
    EMhistory[1,] <- longpars
    bump <- 1000
    #EM
    for (cycles in 1:NCYCLES){
        #priors
        for(g in 1:ngroups){
            gitemtrace[[g]] <- computeItemtrace(pars=pars[[g]], Theta=Theta, itemloc=itemloc)
            pars[[g]] <- assignItemtrace(pars=pars[[g]], itemtrace=gitemtrace[[g]], itemloc=itemloc)
            gstructgrouppars[[g]] <- ExtractGroupPars(pars[[g]][[J+1]])
            if(BFACTOR){
                prior[[g]] <- dnorm(theta, 0, 1)
                prior[[g]] <- prior[[g]]/sum(prior[[g]])                
                Prior[[g]] <- apply(expand.grid(prior[[g]], prior[[g]]), 1, prod)                
                next
            }
            Prior[[g]] <- mvtnorm::dmvnorm(Theta,gstructgrouppars[[g]]$gmeans,
                                           gstructgrouppars[[g]]$gcov)
            Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
        }
        #Estep
        lastLL <- LL
        LL <- 0
        for(g in 1:ngroups){
            if(BFACTOR){
                rlist[[g]] <- Estep.bfactor(pars=pars[[g]], tabdata=PrepList[[g]]$tabdata,
                                            Theta=Theta, prior=prior[[g]],
                                            specific=specific, sitems=sitems,
                                            itemloc=itemloc, itemtrace=gitemtrace[[g]])
            } else {
                rlist[[g]] <- Estep.mirt(pars=pars[[g]], tabdata=PrepList[[g]]$tabdata,
                                         Theta=Theta, prior=Prior[[g]], itemloc=itemloc,
                                         itemtrace=gitemtrace[[g]])
            }
            LL <- LL + sum(r[[g]]*log(rlist[[g]]$expected))
        }
        if(LL < lastLL && cycles > 1) LLwarn <- TRUE
        for(g in 1:ngroups){
            for(i in 1:J){
                tmp <- c(itemloc[i]:(itemloc[i+1] - 1))
                pars[[g]][[i]]@rs <- rlist[[g]]$r1[, tmp]
            }
        }
        preMstep.longpars <- longpars
        if(all(!est)) break
        if(cycles == 2) bump <- 1
        longpars <- Mstep(pars=pars, est=est, longpars=longpars, ngroups=ngroups, J=J,
                      Theta=Theta, Prior=Prior, BFACTOR=BFACTOR, itemloc=itemloc,
                      PrepList=PrepList, L=L, UBOUND=UBOUND, LBOUND=LBOUND,
                      constrain=constrain, TOL=MSTEPTOL, bump=bump)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
        EMhistory[cycles+1,] <- longpars
        if(verbose)
            if(cycles > 1)
                cat(sprintf('\rIteration: %d, Log-Lik: %.3f, Max-Change: %.5f',
                            cycles, LL, max(abs(preMstep.longpars - longpars))))
        if(cycles > 3 && all(abs(preMstep.longpars - longpars) < TOL))  break
    } #END EM
    if(cycles == NCYCLES)
        message('EM iterations terminated after ', cycles, ' iterations.')
    infological <- estpars & !redun_constr
    correction <- numeric(length(estpars[estpars & !redun_constr]))
    names(correction) <- names(estpars[estpars & !redun_constr])
    hess <- matrix(0)
    if(list$SEM){
        h <- matrix(0, nfullpars, nfullpars)
        ind1 <- 1
        for(group in 1:ngroups){
            for (i in 1:J){
                deriv <- Deriv(x=pars[[group]][[i]], Theta=Theta, EM = TRUE, prior=Prior[[group]],
                               estHess=TRUE, BFACTOR=BFACTOR)
                ind2 <- ind1 + length(deriv$grad) - 1
                h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
                ind1 <- ind2 + 1
            }
            i <- i + 1
            deriv <- Deriv(x=pars[[group]][[i]], Theta=Theta, EM = TRUE,
                           pars=pars[[group]], tabdata=PrepList[[group]]$tabdata,
                           itemloc=itemloc, estHess=TRUE)
            ind2 <- ind1 + length(deriv$grad) - 1
            h[ind1:ind2, ind1:ind2] <- pars[[group]][[i]]@hessian <- deriv$hess
            ind1 <- ind2 + 1
        }
        hess <- L %*% h %*% L
        hess <- hess[estpars & !redun_constr, estpars & !redun_constr]
        #return internal functions for SEM.SE
        return(list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                    logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological,
                    estindex_unique=estindex_unique, correction=correction, hess=hess,
                    estpars=estpars & !redun_constr, redun_constr=redun_constr, ngroups=ngroups,
                    LBOUND=LBOUND, UBOUND=UBOUND, EMhistory=na.omit(EMhistory)))
    }
    ret <- list(pars=pars, cycles = cycles, info=matrix(0), longpars=longpars, converge=converge,
                logLik=LL, rlist=rlist, SElogLik=0, L=L, infological=infological,
                estindex_unique=estindex_unique, correction=correction, hess=hess)
    ret
}

# Estep for mirt
Estep.mirt <- function(pars, tabdata, Theta, prior, itemloc, itemtrace=NULL, deriv = FALSE)
{
    nfact <- ncol(Theta)
    nquad <- nrow(Theta)
    J <- length(itemloc) - 1
    r <- tabdata[ ,ncol(tabdata)]
    X <- tabdata[ ,1:(ncol(tabdata) - 1), drop = FALSE]
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
    retlist <- .Call("Estep", itemtrace, prior, X, nfact, r)
    if(deriv) retlist$itemtrace <- itemtrace
    return(retlist)
}

# Estep for bfactor
Estep.bfactor <- function(pars, tabdata, Theta, prior, specific, sitems, itemloc, itemtrace=NULL)
{
    J <- length(itemloc) - 1
    r <- tabdata[ ,ncol(tabdata)]
    X <- tabdata[ ,1:(ncol(tabdata) - 1)]
    if(is.null(itemtrace))
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
    retlist <- .Call("Estepbfactor", itemtrace, prior, X, r, sitems)
    r1 <- matrix(0, nrow(Theta), ncol(X))
    for (i in 1:J){
        r1[ ,itemloc[i]:(itemloc[i+1]-1)] <-
            retlist$r1[ ,itemloc[i]:(itemloc[i+1]-1) + (specific[i] - 1)*ncol(X) ]
    }
    return(list(r1=r1, expected=retlist$expected))
}

Mstep <- function(pars, est, longpars, ngroups, J, Theta, Prior, BFACTOR, itemloc, PrepList, L,
                  UBOUND, LBOUND, constrain, TOL, bump=1){
    p <- longpars[est]
    opt <- optim(p, fn=Mstep.LL, gr=Mstep.grad, method='BFGS', control=list(reltol = TOL*bump, abstol=TOL*bump),
                 est=est, longpars=longpars, pars=pars, ngroups=ngroups, J=J, Theta=Theta,
                 PrepList=PrepList, Prior=Prior, L=L, BFACTOR=BFACTOR, constrain=constrain,
                 UBOUND=UBOUND, LBOUND=LBOUND, itemloc=itemloc)
    longpars[est] <- opt$par
    if(length(constrain) > 0)
        for(i in 1:length(constrain))
            longpars[constrain[[i]][-1]] <- longpars[constrain[[i]][1]]
    i = J + 1
    for(group in 1:ngroups){
        if(any(pars[[group]][[i]]@est)){
            newpars <- Deriv(x=pars[[group]][[i]], Theta=Theta, EM = TRUE,
                           pars=pars[[group]], tabdata=PrepList[[group]]$tabdata,
                           itemloc=itemloc, prior=Prior[[group]])
            longpars[pars[[group]][[i]]@parnum[pars[[group]][[i]]@est]] <- newpars
        }
    }
    return(longpars)
}

Mstep.LL <- function(p, est, longpars, pars, ngroups, J, Theta, PrepList, Prior, L, BFACTOR,
                     constrain, LBOUND, UBOUND, itemloc){
    longpars[est] <- p
    if(length(constrain) > 0)
       for(i in 1:length(constrain))
           longpars[constrain[[i]][-1]] <- longpars[constrain[[i]][1]]
    if(any(longpars < LBOUND | longpars > UBOUND )) return(1e10)
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    LL <- 0
    for(g in 1:ngroups)
        for(i in 1:J)
            LL <- LL + LogLik(pars[[g]][[i]], Theta=Theta, EM=BFACTOR, prior=Prior[[g]])
    if(is.nan(LL)) return(1e10)
    LL
}

Mstep.grad <- function(p, est, longpars, pars, ngroups, J, Theta, PrepList, Prior, L, BFACTOR,
                       constrain, LBOUND, UBOUND, itemloc){
    longpars[est] <- p
    if(length(constrain) > 0)
        for(i in 1:length(constrain))
            longpars[constrain[[i]][-1]] <- longpars[constrain[[i]][1]]
    pars <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    g <- rep(0, ncol(L))
    ind1 <- 1
    for(group in 1:ngroups){
        for (i in 1:J){
            deriv <- Deriv(x=pars[[group]][[i]], Theta=Theta, EM=TRUE, BFACTOR=BFACTOR,
                           prior=Prior[[group]])
            ind2 <- ind1 + length(deriv$grad) - 1
            g[ind1:ind2] <- deriv$grad
            ind1 <- ind2 + 1
        }
        i <- i + 1
        ind2 <- ind1 + length(pars[[group]][[i]]@par) - 1
        g[ind1:ind2] <- 0
        ind1 <- ind2 + 1
    }
    grad <- g %*% L
    grad <- grad[est]
    -grad
}
