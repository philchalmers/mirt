model.elements <- function(model, factorNames, itemtype, nfactNames, nfact, J, K, fulldata,
                           itemloc, data, N, guess, upper, itemnames, exploratory, parprior,
                           parnumber, BFACTOR = FALSE, mixed.design, customItems, customItemsData,
                           dentype, item.Q, customGroup, key, gpcm_mats, spline_args,
                           monopoly.k, dcIRT_nphi = NULL)
{
    hasProdTerms <- ifelse(nfact == nfactNames, FALSE, TRUE)
    prodlist <- NULL
    if(hasProdTerms){
        tmp <- factorNames[grepl('\\(',factorNames)]
        tmp2 <- factorNames[!grepl('\\(',factorNames)]
        tmp <- gsub("\\(","",tmp)
        tmp <- gsub("\\)","",tmp)
        tmp <- gsub(" ","",tmp)
        prodlist <- strsplit(tmp,"\\*")
        for(j in seq_len(length(prodlist))){
            for(i in seq_len(nfact))
                prodlist[[j]][prodlist[[j]] == tmp2[[i]]] <- i
            prodlist[[j]] <- as.numeric(prodlist[[j]])
        }
    }
    #slopes specification
    estlam <- matrix(FALSE, ncol = nfactNames, nrow = J)
    for(i in seq_len(nfactNames)){
        tmp <- model[model[ ,1L] == factorNames[i],2L]
        if(any(regexpr(",",tmp)))
            tmp <- strsplit(tmp,",")[[1L]]
        popout <- c()
        for(j in seq_len(length(tmp))){
            if(regexpr("-",tmp[j]) > 1L){
                popout <- c(popout,j)
                tmp2 <- as.numeric(strsplit(tmp[j],"-")[[1L]])
                tmp2 <- as.character(tmp2[1L]:tmp2[2L])
                tmp <- c(tmp,tmp2)
            }
        }
        if(length(popout != 0L))
            estlam[as.numeric(tmp[-popout]),i] <- TRUE
        else
            estlam[as.numeric(tmp),i] <- TRUE
    }
    lambdas <- ifelse(estlam, .5, 0)
    #INT
    cs <- sqrt(abs(1-rowSums(lambdas^2)))
    lambdas <- lambdas * 1.702
    zetas <- list()
    for(i in seq_len(J)){
        div <- ifelse(cs[i] > .25, cs[i], .25) / 1.702
        if(K[i] == 2L){
            zetas[[i]] <- (-1)*qnorm(mean(fulldata[,itemloc[i]]))/div
        } else if(itemtype[i] %in% c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')){
            zetas[[i]] <- qnorm(mean(fulldata[,itemloc[i] + key[i]-1L]))/div
        } else if(itemtype[i] %in% c('gpcm', 'nominal', 'Rasch')){
            temp <- qnorm(table(data[,i])[1L:(K[i])]/N)
            temp <- (temp - temp[1L])/div
            zetas[[i]] <- temp[-1L]
        } else {
            temp <- table(data[,i])[1L:(K[i]-1L)]/N
            temp <- cumsum(temp)
            zetas[[i]] <- qnorm(1 - temp)/div
        }
    }
    estzetas <- list()
    for(i in seq_len(J))
        estzetas[[i]] <- length(zetas[[i]])
    #COV
    find <- seq_len(nfact)
    estgcov <- matrix(FALSE,nfact,nfact)
    if(any(model[,1L] == 'COV')){
        tmp <- model[model[,1L] == 'COV',2L]
        tmp <- strsplit(tmp,",")[[1L]]
        tmp <- gsub(" ","",tmp)
        for(i in seq_len(length(tmp))){
            tmp2 <- strsplit(tmp[i],"*",fixed=TRUE)[[1L]]
            for(j in seq_len(length(tmp2))){
                for(k in seq_len(length(tmp2))){
                    if(j > k){
                        ind1 <- find[tmp2[k] == factorNames]
                        ind2 <- find[tmp2[j] == factorNames]
                        estgcov[ind2,ind1] <- TRUE
                    }
                }
            }
        }
    }
    gcov <- ifelse(estgcov,.25,0)
    diag(gcov) <- 1
    #MEAN
    gmeans <- rep(0, nfact)
    estgmeans <- rep(FALSE, nfact)
    if(any(model[,1L] == 'MEAN')){
        tmp <- model[model[,1L] == 'MEAN',2L]
        tmp <- strsplit(tmp,",")[[1L]]
        tmp <- gsub(" ","",tmp)
        for(i in seq_len(length(tmp)))
            estgmeans[find[tmp[i] == factorNames]] <- TRUE
    }

    if(exploratory){
        Rpoly <- cormod(data, K, guess)
        loads <- eigen(Rpoly)$vector[,seq_len(nfact), drop = FALSE]
        u <- 1 - rowSums(loads^2)
        u[u < .001 ] <- .2
        cs <- sqrt(u)
        lambdas <- loads/cs * 1.702
        if(!all(itemtype == 'lca'))
            lambdas[!estlam] <- 0
    }
    if(exploratory && any(itemtype %in% c('PC2PL', 'PC3PL')))
        stop('Partially compensatory models can only be estimated within a confirmatory model',
             call.=FALSE)
    if(is.null(item.Q) && any(itemtype == 'lca')){
        item.Q <- vector('list', J)
        for(i in seq_len(J)){
            item.Q[[i]] <- matrix(1, K[i], nfact)
            item.Q[[i]][1L, ] <- 0
        }
    }
    if(!is.null(item.Q)){
        if(length(item.Q) != J)
            stop('item.Q list does not have the correct length', call.=FALSE)
        for(i in seq_len(J)){
            if(is.null(item.Q[[i]])){
                item.Q[[i]] <- matrix(1, K[i], nfact)
                item.Q[[i]][1L, ] <- 0
            }
        }
        tmptest <- sapply(item.Q, nrow) == K & sapply(item.Q, ncol) == nfact & sapply(item.Q, is.matrix)
        if(!all(tmptest))
            stop(sprintf("item.Q list has incorrect matrix structures for the following item(s): %s",
                         paste0(which(!tmptest), collapse=', ')), call.=FALSE)
        if(any(sapply(item.Q, function(x) all(x[1L,] != 0))))
            stop('The first row of ever item.Q matrix must consist only of 0\'s for proper identification',
                 call.=FALSE)
    }
    ret <- LoadPars(itemtype=itemtype, itemloc=itemloc, lambdas=lambdas, zetas=zetas,
                    guess=guess, upper=upper, fulldata=fulldata, J=J, K=K, customItemsData=customItemsData,
                    nfact=nfact+length(prodlist), parprior=parprior, monopoly.k=monopoly.k,
                    parnumber=parnumber, estLambdas=estlam, BFACTOR=BFACTOR, item.Q=item.Q,
                    mixed.design=mixed.design, customItems=customItems, key=key,
                    gpcm_mats=gpcm_mats, spline_args=spline_args, itemnames=itemnames)
    ret[[length(ret) + 1L]] <- LoadGroupPars(gmeans=gmeans, gcov=gcov, estgmeans=estgmeans,
                                             estgcov=estgcov, parnumber=attr(ret, 'parnumber'),
                                             parprior=parprior, Rasch=all(itemtype %in% c('Rasch', 'rsm', 'Tutz')),
                                             customGroup=customGroup, dcIRT_nphi=dcIRT_nphi, dentype=dentype)
    attr(ret, 'prodlist') <- prodlist
    attr(ret, 'exploratory') <- exploratory
    return(ret)
}
