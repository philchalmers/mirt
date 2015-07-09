LoadPars <- function(itemtype, itemloc, lambdas, zetas, guess, upper, fulldata, J, K, nfact,
                     parprior, parnumber, estLambdas, BFACTOR = FALSE, mixed.design, customItems,
                     key, nominal.highlow, gpcm_mats)
{
    customItemNames <- unique(names(customItems))
    if(is.null(customItemNames)) customItemNames <- 'UsElEsSiNtErNaLNaMe'
    valid.items <- Valid_iteminputs()
    invalid.items <- is.na(match(itemtype, valid.items))
    if (any(invalid.items & !(itemtype %in% customItemNames)))
        stop(paste("Unknown itemtype", paste(itemtype[invalid.items], collapse=" ")), call.=FALSE)
    if(length(gpcm_mats)){
        tmp <- sapply(gpcm_mats, ncol)
        if(!all(tmp == nfact))
            stop(paste0('Matricies in gpcm_mats should only have ', nfact, ' column(s)'), call.=FALSE)
        use_gpcm_mats <- as.logical(sapply(gpcm_mats, length))
    } else use_gpcm_mats <- rep(FALSE, length(itemtype))
    pars <- vector('list', J)
    guess <- logit(guess)
    upper <- logit(upper)

    #start values and free parameters
    startvalues <- freepars <- vector('list', J)
    for(i in 1L:J){
        if(any(itemtype[i] == c('Rasch')) && K[i] == 2L){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1
            val <- c(tmpval, zetas[[i]], guess[i], upper[i])
            fp <- c(rep(FALSE,nfact),TRUE,FALSE,FALSE)
            names(val) <- c(paste('a', 1L:nfact, sep=''), 'd', 'g','u')
        } else if(itemtype[i] == 'Rasch' && K[i] > 2L){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1
            if(use_gpcm_mats[i]){
                val <- c(tmpval, as.vector(gpcm_mats[[i]]), 0, zetas[[i]])
                fp <- c(rep(FALSE,nfact), rep(FALSE, K[i]*nfact),
                        FALSE, rep(TRUE, K[i]-1L))
                names(val) <- c(paste('a', 1L:nfact, sep=''),
                                outer(paste0('ak', 0L:(K[i]-1L)), paste0('_', 1L:nfact), FUN=paste0),
                                paste('d', 0L:(K[i]-1L), sep=''))
            } else {
                val <- c(tmpval, 0:(K[i]-1L), 0, zetas[[i]])
                fp <- c(rep(FALSE,nfact), rep(FALSE, K[i]),
                        FALSE, rep(TRUE, K[i]-1L))
                names(val) <- c(paste('a', 1L:nfact, sep=''),
                                paste('ak', 0L:(K[i]-1L), sep=''),
                                paste('d', 0L:(K[i]-1L), sep=''))
            }
        } else if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){
            if(K[i] != 2L)
                stop(paste0('Item ', i, ' requires exactly 2 unique categories'), call.=FALSE)
            val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i])
            fp <- c(estLambdas[i, ], TRUE, FALSE, FALSE)
            if(any(itemtype[i] == c('3PL', '4PL'))) fp[length(fp)-1L] <- TRUE
            if(any(itemtype[i] == c('3PLu', '4PL'))) fp[length(fp)] <- TRUE
            names(val) <- c(paste('a', 1L:nfact, sep=''), 'd', 'g','u')
        } else if(any(itemtype[i] == c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM'))){
            val <- c(lambdas[i,], zetas[[i]], guess[i], upper[i],
                     0, rep(.5, K[i] - 2L), rep(0, K[i]-1L))
            fp <- c(estLambdas[i, ], TRUE, FALSE, FALSE, rep(TRUE, (K[i]-1L)*2))
            fp[c(nfact+4L, length(fp)-(K[i]-2L) )] <- FALSE
            if(any(itemtype[i] == c('3PLNRM', '4PLNRM'))) fp[nfact+2] <- TRUE
            if(any(itemtype[i] == c('3PLuNRM', '4PLNRM'))) fp[nfact+3] <- TRUE
            names(val) <- c(paste('a', 1L:nfact, sep=''), 'd', 'g','u',
                            paste('ak', 0L:(K[i]-2L), sep=''),
                            paste('d', 0L:(K[i]-2L), sep=''))
        } else if(itemtype[i] == 'graded'){
            val <- c(lambdas[i,], zetas[[i]])
            fp <- c(estLambdas[i, ], rep(TRUE, K[i]-1L))
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 1L:(K[i]-1L), sep=''))
        } else if(itemtype[i] == 'grsm'){
            tmp <- zetas[[min(which(K[i] == K))]]
            val <- c(lambdas[i,], tmp, ifelse(min(which(K[i] == K)) == i, 0, tmp[1L] + zetas[[i]][1L]))
            fp <- c(estLambdas[i, ], rep(TRUE, K[i]))
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 1L:(K[i]-1L), sep=''), 'c')
        } else if(itemtype[i] == 'gpcm'){
            if(use_gpcm_mats[i]){
                val <- c(lambdas[i,], as.vector(gpcm_mats[[i]]), 0, zetas[[i]])
                fp <- c(estLambdas[i, ], rep(FALSE, K[i]*nfact),
                        FALSE, rep(TRUE, K[i]-1L))
                names(val) <- c(paste('a', 1L:nfact, sep=''),
                                outer(paste0('ak', 0L:(K[i]-1L)), paste0('_', 1:nfact), FUN=paste0),
                                paste('d', 0L:(K[i]-1L), sep=''))
            } else {
                val <- c(lambdas[i,], 0:(K[i]-1), 0, zetas[[i]])
                fp <- c(estLambdas[i, ], rep(FALSE, K[i]), FALSE, rep(TRUE, K[i]-1L))
                names(val) <- c(paste('a', 1L:nfact, sep=''), paste('ak', 0L:(K[i]-1L), sep=''),
                                paste('d', 0L:(K[i]-1L), sep=''))
            }
#         } else if(itemtype[i] == 'rsm'){
#             tmpval <- rep(0, nfact)
#             tmpval[lambdas[i,] != 0] <- 1
#             val <- c(tmpval, 0:(K[i]-1), 0, seq(2.5, -2.5, length.out = length(zetas[[i]])), 0)
#             fp <- c(rep(FALSE, nfact), rep(FALSE, K[i]), FALSE, rep(TRUE, K[i]))
#             names(val) <- c(paste('a', 1L:nfact, sep=''), paste('ak', 0L:(K[i]-1L), sep=''),
#                             paste('d', 0L:(K[i]-1L), sep=''), 'c')
        } else if(itemtype[i] == 'nominal'){
            val <- c(lambdas[i,], rep(.5, K[i]), rep(0, K[i]))
            fp <- c(estLambdas[i, ], rep(TRUE, K[i]*2))
            if(is.null(nominal.highlow)){
                val[nfact + 1L] <- 0
                val[nfact + K[i]] <- K[[i]] - 1
                fp[nfact + 1L] <- FALSE
                fp[nfact + K[i]] <- FALSE
            } else {
                val[nfact + nominal.highlow[2L, i]] <- 0
                val[nfact + nominal.highlow[1L, i]] <- K[i] - 1
                fp[nfact + nominal.highlow[2L, i]] <- FALSE
                fp[nfact + nominal.highlow[1L, i]] <- FALSE
            }
            fp[c(nfact + K[i] + 1L)] <- FALSE
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('ak', 0L:(K[i]-1L), sep=''),
                            paste('d', 0L:(K[i]-1L), sep=''))
        } else if(any(itemtype[i] == c('PC2PL','PC3PL'))){
            if(K[i] != 2L)
                stop(paste0('Item ', i, ' requires exactly 2 unique categories'), call.=FALSE)
            val <- c(lambdas[i,], rep(1, nfact), guess[i], 999)
            fp <- c(estLambdas[i, ], estLambdas[i, ], FALSE, FALSE)
            if(itemtype[i] == 'PC3PL') fp[length(fp) - 1L] <- TRUE
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 1L:nfact, sep=''), 'g','u')
        } else if(itemtype[i] == 'ideal'){
            if(K[i] != 2L)
                stop(paste0('Item ', i, ' requires exactly 2 unique categories'), call.=FALSE)
            val <- c(lambdas[i,]/2, -0.5)
            fp <- c(estLambdas[i, ], TRUE)
            names(val) <- c(paste('a', 1L:nfact, sep=''), 'd')
        } else if (itemtype[i] %in% c('lca', 'nlca')){
            val <- rep(lambdas[i,], K[i]-1L)
            fp <- rep(TRUE, length(val))
            names(val) <- paste('a', 1L:length(val), sep='')
        }
        if(all(itemtype[i] != valid.items) || itemtype[i] %in% Experimental_itemtypes()) next
        names(fp) <- names(val)
        startvalues[[i]] <- val
        freepars[[i]] <- fp
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
                             mat=use_gpcm_mats[i],
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
                             mat=use_gpcm_mats[i],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=rep(Inf, length(startvalues[[i]])),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

#         if(itemtype[i] == 'rsm'){
#             pars[[i]] <- new('rsm',
#                              par=startvalues[[i]],
#                              nfact=nfact,
#                              ncat=K[i],
#                              itemclass=6L,
#                              nfixedeffects=nfixedeffects,
#                              any.prior=FALSE,
#                              prior.type=rep(0L, length(startvalues[[i]])),
#                              fixed.design=fixed.design.list[[i]],
#                              est=freepars[[i]],
#                              lbound=rep(-Inf, length(startvalues[[i]])),
#                              ubound=rep(Inf, length(startvalues[[i]])),
#                              prior_1=rep(NaN,length(startvalues[[i]])),
#                              prior_2=rep(NaN,length(startvalues[[i]])))
#             pars[[i]]@par[nfact+1L] <- 0
#             tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
#             pars[[i]]@parnum <- tmp2
#             parnumber <- parnumber + length(freepars[[i]])
#             next
#         }

        if(itemtype[i] == 'nominal'){
            pars[[i]] <- new('nominal',
                             par=startvalues[[i]],
                             est=freepars[[i]],
                             mat=FALSE,
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

        if(any(itemtype[i] == c('ideal'))){
            pars[[i]] <- new('ideal', par=startvalues[[i]], est=freepars[[i]],
                             nfact=nfact,
                             ncat=2L,
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=9L,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             lbound=rep(-Inf, length(startvalues[[i]])),
                             ubound=c(rep(Inf, length(startvalues[[i]])-1L), 0),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(any(itemtype[i] %in% c('lca', 'nlca'))){
            pars[[i]] <- new('lca', par=startvalues[[i]], est=freepars[[i]],
                             nfact=nfact,
                             ncat=K[i],
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             score=if(itemtype[i] == 'lca') 0:(K[i]-1) else rep(1, K[i]),
                             itemclass=10L,
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

        if(all(itemtype[i] %in% Experimental_itemtypes())){
            pars[[i]] <- new(itemtype[i], nfact=nfact, ncat=K[i])
            names(pars[[i]]@est) <- names(pars[[i]]@par)
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

        if(all(itemtype[i] != valid.items)){
            #Allowing multiple customItems
            #pars[[i]] <- customItems[[itemtype[i] == names(customItems)]]
            pars[[i]] <- customItems[[which(itemtype[i] == names(customItems))]]

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
               prior.type=rep(0L, length(par)), prior_1=Nans, prior_2=Nans, itemclass=0L)
    return(ret)
}
