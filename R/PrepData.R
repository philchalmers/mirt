PrepData <- function(data, model, itemtype, guess, upper, gpcm_mats, opts,
                     parprior, verbose, technical, parnumber = 1, BFACTOR = FALSE,
                     grsm.block = NULL, rsm.block = NULL, mixed.design, customItems,
                     customGroup, fulldata = NULL, key, spline_args, internal_constraints)
{
    if(is.null(grsm.block)) grsm.block <- rep(1, ncol(data))
    if(is.null(rsm.block)) rsm.block <- rep(1, ncol(data))
    itemnames <- colnames(data)
    keywords <- c('COV', 'CONSTRAIN', 'CONSTRAINB', 'PRIOR', 'MEAN', 'START', 'LBOUND', 'UBOUND',
                  'FIXED', 'FREE', 'NEXPLORE')
    data <- as.matrix(data)
    colnames(data) <- itemnames
    J <- ncol(data)
    N <- nrow(data)
    exploratory <- attr(model, 'exploratory')
    if(length(guess) == 1L) guess <- rep(guess,J)
    if(length(guess) > J || length(guess) < J)
        stop("The number of guessing parameters is incorrect.", call.=FALSE)
    if(length(upper) == 1L) upper <- rep(upper,J)
    if(length(upper) > J || length(upper) < J)
        stop("The number of upper bound parameters is incorrect.", call.=FALSE)
    if(is.null(key) && any(itemtype %in% c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')))
        stop('When using nested logit items a scoring key must be provided with key = c(...)',
             call.=FALSE)
    if(is.null(key))  key <- rep(1L, J)
    if(length(key) != J)
        stop("The number of elements in the key input is incorrect.", call.=FALSE)
    key[is.na(key) | is.nan(key)] <- 1
    key <- as.integer(key)
    uniques <- list()
    for(i in 1L:J){
        uniques[[i]] <- sort(unique(data[,i]))
        if(any(key[i] == uniques[[i]]))
            key[i] <- which(key[i] == uniques[[i]])
    }
    K <- rep(0L,J)
    for(i in 1L:J) K[i] <- length(uniques[[i]])
    if(any(K > 30L) && opts$warn)
        warning(paste0('The following items have a large number of categories which may cause estimation issues: ',
                       paste0(as.character(which(K > 30L)), collapse = " ")), call. = FALSE)
    if(any(itemtype %in% c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM') & K < 3))
        stop('Nested-logit models must have 3 or more categories', call.=FALSE)
    if(!is.null(technical$customK)){
        K <- technical$customK
        for(i in 1L:J)
            uniques[[i]] <- 0L:(K[i]-1L)
    }
    if(any(K < 2L))
        stop('The following items have only one response category and cannot be estimated: ',
             paste(itemnames[K < 2L], ''), call.=FALSE)
    if(is.null(itemtype)) {
        itemtype <- rep('', J)
        for(i in 1L:J){
            if(K[i] > 2L) itemtype[i] <- 'graded'
            if(K[i] == 2L) itemtype[i] <- '2PL'
        }
    }
    if(length(itemtype) == 1L) itemtype <- rep(itemtype, J)
    if(length(itemtype) != J) stop('itemtype specification is not the correct length', call.=FALSE)
    guess[guess == 0 & itemtype %in% c('3PL', '4PL', 'PC3PL', '3PLNRM', '4PLNRM')] <- .15
    upper[upper == 1 & itemtype %in% c('4PL', '3PLu', '3PLuNRM', '4PLNRM')] <- .85
    if(length(gpcm_mats)){
        if(length(gpcm_mats) != ncol(data))
            stop('gpcm_mats list does not correspond to columns in data', call.=FALSE)
        pick <- !sapply(gpcm_mats, is.null) & itemtype %in% c('gpcm', 'Rasch')
        tmp <- gpcm_mats[pick]
        if(!all(sapply(tmp, is.matrix)))
            stop('Matricies must be used in gpcm_mats', call.=FALSE)
        if(!all(sapply(tmp, nrow) == K[pick])){
            nrows <- sapply(tmp, nrow)
            out <- !sapply(tmp, nrow) == K[pick]
            stop(sprintf("Item %i should have %i rows in gpcm_mats, but instead has %i \n  ",
                         seq(1L, J)[out], K[out], nrows[out]), call.=FALSE)
        }
    }
    itemloc <- cumsum(c(1L,K))
    factorNames <- setdiff(model$x[,1L],keywords)
    nfactNames <- length(factorNames)
    nfact <- sum(!grepl('\\(',factorNames))
    index <- 1L:J
    Names <- NULL
    for(i in 1L:J)
        Names <- c(Names, paste("Item.",i,"_",1L:K[i],sep=""))
    if(is.null(fulldata)){
        fulldata <- matrix(0L,N,sum(K))
        colnames(fulldata) <- Names
        for(i in 1L:J){
            ind <- index[i]
            dummy <- matrix(0L,N,K[ind])
            for (j in 0L:(K[ind]-1L))
                dummy[,j+1L] <- as.integer(data[,ind] == uniques[[ind]][j+1L])
            fulldata[ ,itemloc[ind]:(itemloc[ind+1L]-1L)] <- dummy
        }
        fulldata[is.na(fulldata)] <- 0L
    }
    pars <- model.elements(model=model$x, itemtype=itemtype, factorNames=factorNames,
                           nfactNames=nfactNames, nfact=nfact, J=J, K=K, fulldata=fulldata,
                           itemloc=itemloc, data=data, N=N, guess=guess, upper=upper,
                           itemnames=itemnames, exploratory=exploratory, parprior=parprior,
                           parnumber=parnumber, BFACTOR=BFACTOR, mixed.design=mixed.design,
                           customItems=customItems, customGroup=customGroup, key=key,
                           gpcm_mats=gpcm_mats, spline_args=spline_args)
    prodlist <- attr(pars, 'prodlist')
    exploratory <- attr(pars, 'exploratory')
    if(is(pars[[1L]], 'numeric') || is(pars[[1L]], 'logical')){
        names(pars) <- c(itemnames, 'Group_Parameters')
        attr(pars, 'parnumber') <- NULL
        return(pars)
    }
    #within group constraints
    constrain <- list()
    if(internal_constraints){
        if(any(itemtype == 'grsm')){
            unique.grsmgroups <- unique(na.omit(grsm.block))
            for(group in unique.grsmgroups){
                Kk <- unique(K[grsm.block == unique.grsmgroups[group]])
                if(length(Kk) > 1L) stop(c('Rating scale models require items to have the ',
                                           'same number of categories'), call.=FALSE)
                for(k in 1L:(Kk-1L)){
                    grsmConstraint <- c()
                    for(i in 1L:J){
                        if(grsm.block[i] == unique.grsmgroups[group] && itemtype[i] == 'grsm'){
                            if(length(grsmConstraint) == 0L){
                                pars[[i]]@est[length(pars[[i]]@est)] <- FALSE
                                grsmConstraint <- c(grsmConstraint, pars[[i]]@parnum[length(pars[[i]]@parnum)-k])
                            } else grsmConstraint <- c(grsmConstraint, pars[[i]]@parnum[length(pars[[i]]@parnum)-k])
                        }
                    }
                    constrain[[length(constrain) + 1L]] <- grsmConstraint
                }
            }
        }
        if(any(itemtype == 'grsmIRT')){
            unique.grsmgroups <- unique(na.omit(grsm.block))
            for(group in unique.grsmgroups){
                Kk <- unique(K[grsm.block == unique.grsmgroups[group]])
                if(length(Kk) > 1L) stop(c('Rating scale models require items to have the ',
                                           'same number of categories'), call.=FALSE)
                for(k in 1L:(Kk-1L)){
                    grsmConstraint <- c()
                    for(i in 1L:J){
                        if(grsm.block[i] == unique.grsmgroups[group] && itemtype[i] == 'grsmIRT'){
                            if(length(grsmConstraint) == 0L){
                                pars[[i]]@est[length(pars[[i]]@est)] <- FALSE
                                grsmConstraint <- c(grsmConstraint, pars[[i]]@parnum[length(pars[[i]]@parnum)-k])
                            } else grsmConstraint <- c(grsmConstraint, pars[[i]]@parnum[length(pars[[i]]@parnum)-k])
                        }
                    }
                    constrain[[length(constrain) + 1L]] <- grsmConstraint
                }
            }
        }
        if(any(itemtype == 'rsm')){
            unique.rsmgroups <- unique(na.omit(rsm.block))
            for(group in unique.rsmgroups){
                Kk <- unique(K[rsm.block[rsm.block == unique.rsmgroups[group]]])
                if(length(Kk) > 1L) stop('Rating scale models require that items to have the
                                        same number of categories', call.=FALSE)
                for(k in 1L:(Kk-1L)){
                    rsmConstraint <- c()
                    for(i in 1L:J){
                        if(rsm.block[i] == unique.rsmgroups[group] && itemtype[i] == 'rsm'){
                            if(length(rsmConstraint) == 0L){
                                pars[[i]]@est[length(pars[[i]]@est)] <- FALSE
                                rsmConstraint <- c(rsmConstraint, pars[[i]]@parnum[length(pars[[i]]@parnum)-k])
                            } else rsmConstraint <- c(rsmConstraint, pars[[i]]@parnum[length(pars[[i]]@parnum)-k])
                        }
                    }
                    constrain[[length(constrain) + 1L]] <- rsmConstraint
                }
            }
        }
    }
    npars <- sum(sapply(pars, function(x) sum(x@est)))
    if(is.null(prodlist)) prodlist <- list()
    ret <- list(pars=pars, npars=npars, constrain=constrain, prodlist=prodlist, itemnames=itemnames,
                K=K, fulldata=fulldata, nfactNames=nfactNames, nfact=nfact, npars=npars,
                exploratory=exploratory, J=J, itemloc=itemloc, factorNames=factorNames, Names=Names,
                itemtype=itemtype, nLambdas=nfact+length(prodlist))
    ret
}
