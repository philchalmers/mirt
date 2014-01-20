LoadPars <- function(itemtype, itemloc, lambdas, zetas, guess, upper, fulldata, J, K, nfact,
                     parprior, parnumber, D, estLambdas, BFACTOR = FALSE, mixed.design, customItems,
                     key, nominal.highlow)
{
    customItemNames <- unique(names(customItems))
    if(is.null(customItemNames)) customItemNames <- 'UsElEsSiNtErNaLNaMe'
    valid.items <- c('Rasch', '2PL', '3PL', '3PLu', '4PL', 'graded',
                    'grsm', 'gpcm', 'rsm', 'nominal', 'PC2PL','PC3PL',
                    '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')
    invalid.items <- is.na(match(itemtype, valid.items))
    if (any(invalid.items & !(itemtype %in% customItemNames))) {
      stop(paste("Unknown itemtype", paste(itemtype[invalid.items], collapse=" ")))
    }
    pars <- vector('list', J)
    guess <- logit(guess)
    upper <- logit(upper)
    #startvalues
    startvalues <- vector('list', J)
    for(i in 1L:J){
        if(any(itemtype[i] == c('Rasch')) && K[i] == 2L){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1
            val <- c(tmpval, zetas[[i]], guess[i], upper[i])
            names(val) <- c(paste('a', 1L:nfact, sep=''), 'd', 'g','u')
        }
        if(itemtype[i] == 'Rasch' && K[i] > 2L){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1
            val <- c(tmpval, 0, zetas[[i]])
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 0L:(K[i]-1L), sep=''))
        }
        if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){
            val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i])
            names(val) <- c(paste('a', 1L:nfact, sep=''), 'd', 'g','u')
        }
        if(any(itemtype[i] == c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM'))){
            val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i],
                     0, rep(.5, K[i] - 2L), rep(0, K[i]-1L))
            names(val) <- c(paste('a', 1L:nfact, sep=''), 'd', 'g','u',
                            paste('ak', 0L:(K[i]-2L), sep=''),
                            paste('d', 0L:(K[i]-2L), sep=''))
        }
        if(itemtype[i] == 'graded'){
            val <- c(lambdas[i,], zetas[[i]])
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 1L:(K[i]-1L), sep=''))
        }
        if(itemtype[i] == 'grsm'){
            val <- c(lambdas[i,], seq(2.5, -2.5, length.out = length(zetas[[i]])), 0)
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 1L:(K[i]-1L), sep=''), 'c')
        }
        if(itemtype[i] == 'gpcm'){
            val <- c(lambdas[i,], 0, zetas[[i]])
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 0L:(K[i]-1L), sep=''))
        }
        if(itemtype[i] == 'rsm'){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1
            val <- c(tmpval, 0, seq(2.5, -2.5, length.out = length(zetas[[i]])), 0)
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 0L:(K[i]-1L), sep=''), 'c')
        }
        if(itemtype[i] == 'nominal'){
            val <- c(lambdas[i,], rep(.5, K[i]), rep(0, K[i]))
            if(is.null(nominal.highlow)){
                val[nfact + 1L] <- 0
                val[nfact + K[i]] <- K[[i]] - 1
            } else {
                val[nfact + nominal.highlow[2L, i]] <- 0
                val[nfact + nominal.highlow[1L, i]] <- K[i] - 1
            }
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('ak', 0L:(K[i]-1L), sep=''),
                            paste('d', 0L:(K[i]-1L), sep=''))
        }
        if(any(itemtype[i] == c('PC2PL','PC3PL'))){
            val <- c(lambdas[i,], rep(1, nfact), guess[i], 999)
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 1L:nfact, sep=''), 'g','u')
        }
        if(all(itemtype[i] != valid.items)) next
        startvalues[[i]] <- val
    }
    #freepars
    freepars <- vector('list', J)
    for(i in 1L:J){
        if(itemtype[i] == 'Rasch' && K[i] == 2L)
            freepars[[i]] <- c(rep(FALSE,nfact),TRUE,FALSE,FALSE)
        if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){
            if(K[i] != 2L) stop(paste0('Item ', i, ' requires exactly 2 unique categories'))
            estpars <- c(estLambdas[i, ], TRUE, FALSE, FALSE)
            if(any(itemtype[i] == c('3PL', '4PL'))) estpars[length(estpars)-1L] <- TRUE
            if(any(itemtype[i] == c('3PLu', '4PL'))) estpars[length(estpars)] <- TRUE
            freepars[[i]] <- estpars
        }
        if(any(itemtype[i] == c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM'))){
            estpars <- c(estLambdas[i, ], TRUE, FALSE, FALSE, rep(TRUE, (K[i]-1L)*2))
            estpars[c(nfact+4L, length(estpars)-(K[i]-2L) )] <- FALSE
            if(any(itemtype[i] == c('3PLNRM', '4PLNRM'))) estpars[nfact+2] <- TRUE
            if(any(itemtype[i] == c('3PLuNRM', '4PLNRM'))) estpars[nfact+3] <- TRUE
            freepars[[i]] <- estpars
        }
        if(itemtype[i] == 'Rasch' && K[i] > 2L)
            freepars[[i]] <- c(rep(FALSE,nfact), FALSE, rep(TRUE, K[i]-1L))
        if(itemtype[i] == 'grsm')
            freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]))
        if(itemtype[i] == 'graded')
            freepars[[i]] <- c(estLambdas[i, ], rep(TRUE, K[i]-1L))
        if(itemtype[i] == 'gpcm')
            freepars[[i]] <- c(estLambdas[i, ], FALSE, rep(TRUE, K[i]-1L))
        if(itemtype[i] == 'rsm')
            freepars[[i]] <- c(rep(FALSE, nfact), FALSE, rep(TRUE, K[i]))
        if(itemtype[i] == 'nominal'){
            estpars <- c(estLambdas[i, ], rep(TRUE, K[i]*2))
            if(is.null(nominal.highlow)){
                estpars[nfact + 1L] <- FALSE
                estpars[nfact + K[i]] <- FALSE
            } else {
                estpars[nfact + nominal.highlow[2L, i]] <- FALSE
                estpars[nfact + nominal.highlow[1L, i]] <- FALSE
            }
            estpars[c(nfact + K[i] + 1L)] <- FALSE
            freepars[[i]] <- estpars
        }
        if(any(itemtype[i] == c('PC2PL','PC3PL'))){
            if(K[i] != 2L) stop(paste0('Item ', i, ' requires exactly 2 unique categories'))
            estpars <- c(estLambdas[i, ], estLambdas[i, ], FALSE, FALSE)
            if(itemtype[i] == 'PC3PL') estpars[length(estpars) - 1L] <- TRUE
            freepars[[i]] <- estpars
        }
        if(all(itemtype[i] != valid.items)) next
        names(freepars[[i]]) <- names(startvalues[[i]])
    }

    #augment startvalues and fixedpars for mixed effects
    nfixedeffects <- 0
    fixed.design.list <- vector('list', J)
    for(i in 1L:J) fixed.design.list[[i]] <- matrix(0)
    if(!is.null(mixed.design)){
        fixed.design <- mixed.design$fixed
        betas <- rep(0, ncol(fixed.design))
        estbetas <- rep(TRUE, length(betas))
        names(estbetas) <- names(betas) <- colnames(fixed.design)
        nfixedeffects <- length(betas)
        nfact <- nfact + nfixedeffects
        for(i in 1L:J){
            freepars[[i]] <- c(estbetas, freepars[[i]])
            startvalues[[i]] <- c(betas, startvalues[[i]])
        }
        valid.ints <- ifelse(any(K > 2), '', 'd')
        freepars <- lapply(freepars, function(x, valid){
            x[names(x) %in% valid] <- FALSE
            return(x)}, valid=valid.ints)
        startvalues <- lapply(startvalues, function(x, valid){
            x[names(x) %in% valid] <- 0
            return(x)}, valid=valid.ints)
        N <- nrow(mixed.design$fixed) / J
        for(i in 1L:J)
            fixed.design.list[[i]] <- mixed.design$fixed[1L:N + N*(i-1L), , drop = FALSE]
    }

    #load items
    nfact <- as.integer(nfact)
    K <- as.integer(K)
    for(i in 1L:J){
        tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L)) #item location

        if(any(itemtype[i] == c('Rasch')) && K[i] == 2L){
            pars[[i]] <- new('dich', par=startvalues[[i]], est=freepars[[i]],
                             nfact=nfact,
                             ncat=2L,
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=1L,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(itemtype[i] == 'Rasch' && K[i] > 2L){
            pars[[i]] <- new('gpcm',
                             par=startvalues[[i]],
                             nfact=nfact,
                             ncat=K[i],
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=3L,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             est=freepars[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            pars[[i]]@par[nfact+1L] <- 0
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){
            pars[[i]] <- new('dich',
                             par=startvalues[[i]],
                             est=freepars[[i]],
                             nfact=nfact,
                             itemclass=1L,
                             nfixedeffects=nfixedeffects,
                             ncat=2L,
                             any.prior=FALSE,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(any(itemtype[i] == c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM'))){
            pars[[i]] <- new('nestlogit',
                             par=startvalues[[i]],
                             est=freepars[[i]],
                             nfact=nfact,
                             itemclass=8L,
                             nfixedeffects=nfixedeffects,
                             ncat=K[i],
                             correctcat=key[i],
                             any.prior=FALSE,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(any(itemtype[i] == 'grsm')){
            pars[[i]] <- new('rating',
                             par=startvalues[[i]],
                             nfact=nfact,
                             ncat=K[i],
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=5L,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             est=freepars[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(itemtype[i] == 'graded'){
            pars[[i]] <- new('graded',
                             par=startvalues[[i]],
                             nfact=nfact,
                             ncat=K[i],
                             itemclass=2L,
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             est=freepars[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(itemtype[i] == 'gpcm'){
            pars[[i]] <- new('gpcm',
                             par=startvalues[[i]],
                             nfact=nfact,
                             ncat=K[i],
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=3L,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             est=freepars[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            pars[[i]]@par[nfact+1L] <- 0
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(itemtype[i] == 'rsm'){
            pars[[i]] <- new('rsm',
                             par=startvalues[[i]],
                             nfact=nfact,
                             ncat=K[i],
                             itemclass=6L,
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             est=freepars[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            pars[[i]]@par[nfact+1L] <- 0
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(itemtype[i] == 'nominal'){
            pars[[i]] <- new('nominal',
                             par=startvalues[[i]],
                             est=freepars[[i]],
                             nfact=nfact,
                             ncat=K[i],
                             itemclass=4L,
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(any(itemtype[i] == c('PC2PL','PC3PL'))){
            pars[[i]] <- new('partcomp',
                             par=startvalues[[i]],
                             est=freepars[[i]],
                             nfact=nfact,
                             ncat=2L,
                             itemclass=7L,
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(all(itemtype[i] != valid.items)){
            pars[[i]] <- customItems[[itemtype[i] == names(customItems)]]
            pars[[i]]@nfact <- nfact
            pars[[i]]@ncat <- K[i]
            pars[[i]]@nfixedeffects <- nfixedeffects
            pars[[i]]@any.prior <- FALSE
            pars[[i]]@itemclass <- 9L
            pars[[i]]@prior.type <- rep(0L, length(pars[[i]]@par))
            pars[[i]]@prior_1 <- rep(NaN,length(pars[[i]]@par))
            pars[[i]]@prior_2 <- rep(NaN,length(pars[[i]]@par))
            tmp2 <- parnumber:(parnumber + length(pars[[i]]@est) - 1L)
            pars[[i]]@parnum <- tmp2
            pars[[i]]@fixed.design <- fixed.design.list[[i]]
            parnumber <- parnumber + length(pars[[i]]@est)
            next
        }
    }
    #priors
    for(i in 1L:J){
        names(pars[[i]]@parnum) <- names(startvalues[[i]])
        if(!is.null(parprior) && parprior != 'index'){
            for(j in 1L:length(parprior)){
                tmp <- pars[[i]]@parnum %in% as.numeric(parprior[[j]][1L:(length(parprior[[j]])-3L)])
                if(any(tmp)){
                    pars[[i]]@any.prior <- TRUE
                    p2 <- as.numeric(parprior[[j]][length(parprior[[j]])])
                    p1 <- as.numeric(parprior[[j]][length(parprior[[j]])-1L])
                    type <- parprior[[j]][length(parprior[[j]])-2L]
                    pars[[i]]@prior.type[tmp] <- switch(type, norm=1L, lnorm=2L, beta=3L, 0L)
                    pars[[i]]@prior_1[tmp] <- p1
                    pars[[i]]@prior_2[tmp] <- p2
                }
            }
        }
    }
    attr(pars, 'parnumber') <- parnumber
    return(pars)
}

LoadGroupPars <- function(gmeans, gcov, estgmeans, estgcov, parnumber, parprior, Rasch = FALSE){
    nfact <- length(gmeans)
    fn <- paste('COV_', 1L:nfact, sep='')
    FNCOV <- outer(fn, 1L:nfact, FUN=paste, sep='')
    FNMEANS <- paste('MEAN_', 1L:nfact, sep='')
    tri <- lower.tri(gcov, diag=TRUE)
    par <- c(gmeans, gcov[tri])
    parnum <- parnumber:(parnumber + length(par) - 1L)
    if(Rasch) diag(estgcov) <- TRUE
    est <- c(estgmeans,estgcov[tri])
    names(parnum) <- names(par) <- names(est) <- c(FNMEANS,FNCOV[tri])
    tmp <- matrix(-Inf, nfact, nfact)
    diag(tmp) <- 1e-4
    lbound <- c(rep(-Inf, nfact), tmp[tri])
    Nans <- rep(NaN,length(par))
    ret <- new('GroupPars', par=par, est=est, nfact=nfact, any.prior=FALSE,
               parnum=parnum, lbound=lbound, ubound=rep(Inf, length(par)),
               prior.type=rep(0L, length(par)), prior_1=Nans, prior_2=Nans)
    return(ret)
}
