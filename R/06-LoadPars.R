LoadPars <- function(itemtype, itemloc, lambdas, zetas, guess, upper, fulldata, J, K, nfact, data,
                     parprior, parnumber, estLambdas, BFACTOR = FALSE, mixed.design, customItems,
                     key, gpcm_mats, spline_args, itemnames, monopoly.k, customItemsData, item.Q)
{
    customItemNames <- unique(names(customItems))
    if(is.null(customItemNames)) customItemNames <- 'UsElEsSiNtErNaLNaMe'
    valid.items <- Valid_iteminputs()
    invalid.items <- is.na(match(itemtype, valid.items))
    if(any(invalid.items & !(itemtype %in% customItemNames)))
        stop(paste("Unknown itemtype:", paste(itemtype[invalid.items], collapse=" ")), call.=FALSE)
    if(any(itemtype %in% c('gpcmIRT', 'monopoly', 'grsmIRT', 'crm')) && nfact > 1L)
        stop('Multidimensional model not supported for select itemtype(s)', call.=FALSE)
    if(any(itemtype == 'spline' & K > 2L))
        stop('spline itemtype only supported for dichotomous items', call.=FALSE)
    if(length(gpcm_mats)){
        tmp <- sapply(gpcm_mats, ncol)
        if(!all(tmp == nfact))
            stop(paste0('Matricies in gpcm_mats should only have ', nfact, ' column(s)'), call.=FALSE)
        use_gpcm_mats <- as.logical(sapply(gpcm_mats, length))
    } else use_gpcm_mats <- rep(FALSE, length(itemtype))
    pars <- vector('list', J)
    guess <- logit(guess)
    upper <- logit(upper)
    ggum.start.values <- vector('list', length(K))
    if(any(itemtype == 'ggum')){
        dca <- try(vegan::decorana(fulldata), TRUE)
        if(is(dca, 'try-error') || nfact > 4){
            for (i in 1L:length(K))
                ggum.start.values[[i]] <- c(rep(1, nfact), numeric(nfact),
                                            seq(3, -3, length.out = K[i]-1))
        } else {
            tmp <- utils::capture.output(a <- as.list(summary(dca, digits=5, origin=TRUE,
                                                       display="species")))
            data.dca <- as.data.frame(a$spec.scores)
            dca.mat <- as.matrix(data.dca[,1:nfact])

            for (i in 1L:length(K)) {
                if(itemtype[i] != 'ggum') next

                dca.dist <- 0
                tmppar <- numeric(2*nfact + K[i]-1)
                for (d in 1:nfact) {
                    tmppar[d] <- 1   #alphas
                    tmppar[nfact+d] <- dca.mat[i,d]   #deltas
                    dca.dist <- dca.mat[i,d]^2 + dca.dist
                }

                for (k in 1:(K[i]-1)) {
                    origin <- 1.002+.449*sqrt(dca.dist) - .093*K[i]
                    delta <- .921+.058*sqrt(dca.dist) - .129*K[i]
                    tmppar[2*nfact+k] <- origin + delta*(K[i]-k) #taus
                }
                ggum.start.values[[i]] <- tmppar
            }
        }
    }


    #start values and free parameters
    startvalues <- freepars <- vector('list', J)
    for(i in seq_len(J)){
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
        } else if(itemtype[i] %in% c('sequential', 'Tutz')){
            val <- c(lambdas[i,], zetas[[i]])
            fp <- c(estLambdas[i, ], rep(TRUE, K[i]-1L))
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('d', 1L:(K[i]-1L), sep=''))
            if(itemtype[i] == 'Tutz'){
                val[1L] <- 1
                fp[1L] <- FALSE
            }
        } else if(itemtype[i] == 'monopoly'){
            val <- c(suppressWarnings(log(lambdas[i,])),
                     zetas[[i]], rep(0, monopoly.k[i]*2))
            if(!is.finite(val[1L])) val[1L] <- -4
            fp <- c(estLambdas[i, ], rep(TRUE, K[i]-1L + monopoly.k[i]*2))
            names(val) <- c('omega', paste('xi', 1L:(K[i]-1L), sep=''),
                            paste(c('alpha', 'tau'), 1L:(monopoly.k[i]*2), sep=''))
        } else if(itemtype[i] %in% c('gpcmIRT', 'rsm')){
            if(itemtype[i] == 'rsm'){
                val <- c(1, seq(-2, 2, length.out=K[i]-1), 0)
                fp <- c(FALSE, rep(TRUE, K[i]-1L), TRUE)
            } else {
                val <- c(lambdas[i,], -zetas[[i]], 0)
                fp <- c(estLambdas[i, ], rep(TRUE, K[i]-1L), FALSE)
            }
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('b', 1L:(K[i]-1L), sep=''), 'c')
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
        } else if(itemtype[i] == 'rsm'){
            tmpval <- rep(0, nfact)
            tmpval[lambdas[i,] != 0] <- 1
            val <- c(tmpval, 0:(K[i]-1), 0, seq(2.5, -2.5, length.out = length(zetas[[i]])), 0)
            fp <- c(rep(FALSE, nfact), rep(FALSE, K[i]), FALSE, rep(TRUE, K[i]))
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('ak', 0L:(K[i]-1L), sep=''),
                            paste('d', 0L:(K[i]-1L), sep=''), 'c')
        } else if(itemtype[i] == 'nominal'){
            val <- c(lambdas[i,], rep(.5, K[i]), rep(0, K[i]))
            fp <- c(estLambdas[i, ], rep(TRUE, K[i]*2))
            val[nfact + 1L] <- 0
            val[nfact + K[i]] <- K[[i]] - 1
            fp[nfact + 1L] <- FALSE
            fp[nfact + K[i]] <- FALSE
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
        } else if (itemtype[i] == 'lca'){
            if(K[i] == 2L){
                tmp <- length(estLambdas[i, ])
                val <- seq(-1 - log(tmp), 1 + log(tmp), length.out = tmp)
                fp <- estLambdas[i, ]
                names(val) <- paste('a', 1L:length(val), sep='')
            } else {
                tmp <- length(estLambdas[i, ])
                val <- rep(seq(-1 - log(tmp), 1 + log(tmp), length.out = tmp), K[i]-1L)
                fp <- rep(TRUE, length(val))
                names(val) <- paste('a', 1L:length(val), sep='')
            }
            val[!fp] <- 0
        } else if (itemtype[i] == 'ggum'){
            val <- ggum.start.values[[i]]
            val[!estLambdas[i,]] <- 0
            fp <- c(estLambdas[i,], estLambdas[i,], rep(TRUE,length(val)-(ncol(estLambdas)*2)))
            names(val) <- c(paste('a', 1L:nfact, sep=''), paste('b', 1L:nfact, sep=''),
                paste('t', 1L:(max(K[i]) - 1), sep=''))
            names(fp) <- names(val)
            startvalues[[i]] <- val
            freepars[[i]] <- fp
        } else if (itemtype[i] == 'crm'){
            val <- c(lambdas[i,], -mean(qlogis(data[,i]), na.rm=TRUE), 1)
            names(val) <- c(paste0('a', 1L:nfact), 'b', 'alpha')
            fp <- c(estLambdas[i,], rep(TRUE, 2))
            names(fp) <- names(val)
            startvalues[[i]] <- val
            freepars[[i]] <- fp
        } else if (itemtype[i] == 'spline') next
        if(all(itemtype[i] != valid.items) || itemtype[i] %in% Experimental_itemtypes()) next
        names(fp) <- names(val)
        startvalues[[i]] <- val
        freepars[[i]] <- fp
    }

    #augment startvalues and fixedpars for mixed effects
    nfixedeffects <- 0
    fixed.design.list <- vector('list', J)
    for(i in seq_len(J)) fixed.design.list[[i]] <- matrix(0)
    if(!is.null(mixed.design)){
        fixed.design <- mixed.design$fixed
        betas <- rep(0, ncol(fixed.design))
        estbetas <- rep(TRUE, length(betas))
        names(estbetas) <- names(betas) <- colnames(fixed.design)
        nfixedeffects <- length(betas)
        nfact <- nfact + nfixedeffects
        for(i in seq_len(J)){
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
        for(i in seq_len(J))
            fixed.design.list[[i]] <- mixed.design$fixed[1L:N + N*(i-1L), , drop = FALSE]
    }

    #load items
    nfact <- as.integer(nfact)
    K <- as.integer(K)
    for(i in seq_len(J)){
        tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L)) #item location

        if(any(itemtype[i] == c('Rasch')) && K[i] == 2L){
            pars[[i]] <- new('dich', par=startvalues[[i]], est=freepars[[i]],
                             parnames=names(freepars[[i]]),
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
                             parnames=names(freepars[[i]]),
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
                             parnames=names(freepars[[i]]),
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
                             parnames=names(freepars[[i]]),
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
                             parnames=names(freepars[[i]]),
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
                             parnames=names(freepars[[i]]),
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

        if(itemtype[i] == 'crm'){
            pars[[i]] <- new('crm',
                             par=startvalues[[i]],
                             parnames=names(freepars[[i]]),
                             nfact=nfact,
                             ncat=1L,
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=13L,
                             orgdat=data[,i,drop=FALSE],
                             transdat=qlogis(data[,i,drop=FALSE]),
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             est=freepars[[i]],
                             lbound=c(0,-Inf, 0),
                             ubound=c(Inf, Inf, Inf),
                             prior_1=rep(NaN,length(startvalues[[i]])),
                             prior_2=rep(NaN,length(startvalues[[i]])))
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(itemtype[i] %in% c('sequential', 'Tutz')){
            pars[[i]] <- new('sequential',
                             par=startvalues[[i]],
                             parnames=names(freepars[[i]]),
                             nfact=nfact,
                             ncat=K[i],
                             itemclass=9L,
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
                             parnames=names(freepars[[i]]),
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

        if(itemtype[i] %in% c('gpcmIRT', 'rsm')){
            pars[[i]] <- new('gpcmIRT',
                             par=startvalues[[i]],
                             parnames=names(freepars[[i]]),
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
            tmp2 <- parnumber:(parnumber + length(freepars[[i]]) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(freepars[[i]])
            next
        }

        if(itemtype[i] == 'monopoly'){
            pars[[i]] <- new('monopoly',
                             par=startvalues[[i]],
                             parnames=names(freepars[[i]]),
                             nfact=nfact,
                             ncat=K[i],
                             k=as.integer(monopoly.k[i]),
                             itemclass=12L,
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

        if(itemtype[i] == 'nominal'){
            pars[[i]] <- new('nominal',
                             par=startvalues[[i]],
                             parnames=names(freepars[[i]]),
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
                             parnames=names(freepars[[i]]),
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
                             parnames=names(freepars[[i]]),
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

        if(any(itemtype[i] == 'lca')){
            pars[[i]] <- new('lca', par=startvalues[[i]], est=freepars[[i]],
                             parnames=names(freepars[[i]]),
                             nfact=nfact,
                             ncat=K[i],
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=10L,
                             item.Q=item.Q[[i]],
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

        if(itemtype[i] == 'spline'){
            stype <- 'bs'
            intercept <- TRUE
            df <- knots <- NULL
            degree <- 3
            if(any(names(spline_args) == itemnames[i])){
                sargs <- spline_args[[which(names(spline_args) == itemnames[i])]]
                if(!is.null(sargs$fun)) stype <- sargs$fun
                if(!is.null(sargs$intercept)) intercept <- sargs$intercept
                if(!is.null(sargs$df)) df <- sargs$df
                if(!is.null(sargs$knots)) knots <- sargs$knots
                if(!is.null(sargs$degree)) degree <- sargs$degree
            }
            sargs <- list(stype=stype, intercept=intercept, df=df, knots=knots, degree=degree)
            Theta_prime <- if(stype == 'bs'){
                splines::bs(c(-2,2), df=df, knots=knots, degree=degree, intercept=intercept)
            } else if(stype == 'ns'){
                splines::ns(c(-2,2), df=df, knots=knots, intercept=intercept)
            } else stop('splines function not supported', call.=FALSE)
            p <- seq(-10, 10, length.out=ncol(Theta_prime))
            est <- rep(TRUE, ncol(Theta_prime))
            names(est) <- paste0('s', 1L:length(p))
            pars[[i]] <- new('spline', par=p,
                             parnames=names(est),
                             est=est,
                             nfact=nfact,
                             ncat=K[i],
                             stype=stype,
                             item.Q=matrix(1, K[i], length(p)),
                             Theta_prime=matrix(0),
                             sargs=sargs,
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=11L,
                             prior.type=rep(0L, length(p)),
                             fixed.design=fixed.design.list[[i]],
                             lbound=rep(-Inf, length(p)),
                             ubound=rep(Inf, length(p)),
                             prior_1=rep(NaN,length(p)),
                             prior_2=rep(NaN,length(p)))
            pars[[i]]@item.Q[1L, ] <- 0
            tmp2 <- parnumber:(parnumber + length(p) - 1L)
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(p)
            next
        }

        if(itemtype[i] == 'ggum'){
            pars[[i]] <- new('ggum',
                             par=startvalues[[i]],
                             est=freepars[[i]],
                             parnames=names(freepars[[i]]),
                             nfact=nfact,
                             ncat=K[i],
                             nfixedeffects=nfixedeffects,
                             any.prior=FALSE,
                             itemclass=11L,
                             prior.type=rep(0L, length(startvalues[[i]])),
                             fixed.design=fixed.design.list[[i]],
                             lbound=c(rep(1e-4, nfact),
                                      rep(-Inf, length(startvalues[[i]])-nfact)),
                             ubound=c(rep(Inf, length(startvalues[[i]]))),
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
            pars[[i]]@parnames <- names(pars[[i]]@est)
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
            if(pars[[i]]@useuserdata) pars[[i]]@userdata <- list(customItemsData[[i]])
            parnumber <- parnumber + length(pars[[i]]@est)
            next
        }

        # continuous or discrete
        .object@discrete <- if(Continuous_itemtypes() %in% itemtype[i]) TRUE else FALSE

    }
    #priors
    for(i in seq_len(J)){
        names(pars[[i]]@parnum) <- names(pars[[i]]@par)
        if(!is.null(parprior)){
            for(j in seq_len(length(parprior))){
                tmp <- pars[[i]]@parnum %in% as.numeric(parprior[[j]][1L:(length(parprior[[j]])-3L)])
                if(any(tmp)){
                    pars[[i]]@any.prior <- TRUE
                    p2 <- as.numeric(parprior[[j]][length(parprior[[j]])])
                    p1 <- as.numeric(parprior[[j]][length(parprior[[j]])-1L])
                    type <- parprior[[j]][length(parprior[[j]])-2L]
                    pars[[i]]@prior.type[tmp] <- switch(type, norm=1L, lnorm=2L, beta=3L, expbeta=4L, 0L)
                    whc <- which(tmp)
                    for(w in 1:length(whc)){
                        pars[[i]]@par[whc[w]] <- switch(type[w],
                                                     'norm'=pars[[i]]@par[whc[w]],
                                                     'lnorm'=exp(p1[w]),
                                                     'beta'=(p1[w]-1)/(p1[w] + p2[w] - 2),
                                                     'expbeta'=expbeta_sv(p1[w], p2[w]))
                        if(type[w] == 'lnorm')
                            pars[[i]]@lbound[whc[w]] <- 0
                        if(type[w] == 'beta'){
                            pars[[i]]@lbound[whc[w]] <- 0
                            pars[[i]]@ubound[whc[w]] <- 1
                        }
                    }
                    pars[[i]]@prior_1[tmp] <- p1
                    pars[[i]]@prior_2[tmp] <- p2
                }
            }
        }
    }
    attr(pars, 'parnumber') <- parnumber
    return(pars)
}

LoadGroupPars <- function(gmeans, gcov, estgmeans, estgcov, parnumber, parprior, dentype, Rasch = FALSE,
                          customGroup = NULL, dcIRT_nphi = NULL, ...){
    if (dentype != 'Davidian') {
        if(!is.null(customGroup)){
            par <- customGroup@par
            parnum <- parnumber:(parnumber + length(par) - 1L)
            customGroup@parnum <- parnum
            return(customGroup)
        } else {
            nfact <- length(gmeans)
            den <- Theta_mvtnorm_den
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
            if(dentype == 'mixture'){
                par <- c(par, "PI" = 0)
                est <- c(est, "PI" = TRUE)
                parnum <- c(parnum, max(parnum)+1L)
                Nans <- c(Nans, NaN)
                lbound <- c(lbound, -Inf)
            }
            ubound <- rep(Inf, length(par))
            ret <- new('GroupPars', par=par, est=est, parnames=names(est),
                       nfact=nfact, any.prior=FALSE, den=den,
                       safe_den=den, parnum=parnum, lbound=lbound, ubound=ubound,
                       prior.type=rep(0L, length(par)), prior_1=Nans, prior_2=Nans, rrb=0, rrs=matrix(0),
                       BFACTOR=FALSE, itemclass=0L, dentype=dentype)
            return(ret)
        }
    } else {
        nfact <- length(gmeans)
        if (nfact > 1) stop("Multidimensional DC-IRT models are not supported.")
        # DC Density:
        den <- Theta_DC_den
        fn <- paste('COV_', 1L:nfact, sep='')
        FNCOV <- outer(fn, 1L:nfact, FUN=paste, sep='')
        FNMEANS <- paste('MEAN_', 1L:nfact, sep='')
        FNPHI <- paste('PHI_', 1L:dcIRT_nphi, sep='')
        tri <- lower.tri(gcov, diag=TRUE)
        phi <- rep(1.570789, dcIRT_nphi)
        estphi <- rep(T, dcIRT_nphi)
        par <- c(gmeans, gcov[tri], phi)
        parnum <- parnumber:(parnumber + length(par) - 1L)
        if(Rasch) diag(estgcov) <- TRUE
        est <- c(estgmeans,estgcov[tri], estphi)
        names(parnum) <- names(par) <- names(est) <- c(FNMEANS,FNCOV[tri],FNPHI)
        tmp <- matrix(-Inf, nfact, nfact)
        diag(tmp) <- 1e-4
        lbound <- c(rep(-Inf, nfact), tmp[tri], rep(-pi/2, dcIRT_nphi))
        ubound <- c(rep(Inf, length(par)-dcIRT_nphi), rep(pi/2, dcIRT_nphi))
        Nans <- rep(NaN,length(par))
        ret <- new('GroupPars', par=par, est=est, parnames=names(est),
               nfact=nfact, any.prior=FALSE, den=den,
               safe_den=den, parnum=parnum, lbound=lbound, ubound=ubound,
               prior.type=rep(0L, length(par)), prior_1=Nans, prior_2=Nans, rrb=0, rrs=matrix(0),
               BFACTOR=FALSE, itemclass=0L, dentype = 'Davidian')
        return(ret)
      }
}

Theta_mvtnorm_den <- function(obj, Theta){
    gpars <- ExtractGroupPars(obj)
    mu <- gpars$gmeans
    sigma <- gpars$gcov
    d <- mirt_dmvnorm(Theta, mean=mu, sigma=sigma)
    d <- ifelse(d < 1e-300, 1e-300, d)
    d
}

Theta_DC_den <- function(obj, Theta) {
    gpars <- ExtractGroupPars(obj)
    phi <- gpars$phi
    d <- dcurver::ddc(Theta, phi = gpars$phi)
    d
}

Theta_discrete_den <- function(obj, Theta, mus = 0){
    if(length(Theta) == 1) return(1)
    par <- obj@par
    if(length(mus) > 1L){
        ret <- t(apply(mus, 1L, function(x)
            c(exp(par + x[-ncol(mus)]), 1)))
        ret <- ret / rowSums(ret)
    } else if(length(obj@structure)){
        d <- exp(obj@structure %*% par)
        d[length(d)] <- 1
        ret <- as.vector(d / sum(d))
    } else {
        d <- c(exp(par), 1)
        ret <- d / sum(d)
    }
    ret
}