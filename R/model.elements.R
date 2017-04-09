model.elements <- function(model, factorNames, itemtype, nfactNames, nfact, J, K, fulldata,
                           itemloc, data, N, guess, upper, itemnames, exploratory, parprior,
                           parnumber, BFACTOR = FALSE, mixed.design, customItems,
                           customGroup, key, gpcm_mats, spline_args)
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
    ret <- LoadPars(itemtype=itemtype, itemloc=itemloc, lambdas=lambdas, zetas=zetas,
                    guess=guess, upper=upper, fulldata=fulldata, J=J, K=K,
                    nfact=nfact+length(prodlist), parprior=parprior,
                    parnumber=parnumber, estLambdas=estlam, BFACTOR=BFACTOR,
                    mixed.design=mixed.design, customItems=customItems, key=key,
                    gpcm_mats=gpcm_mats, spline_args=spline_args, itemnames=itemnames)
    if(any(model[,1L] == 'START')){
        input <- gsub(" ","", model[model[,1L] == 'START', 2L])
        elements <- strsplit(input, '\\),\\(')[[1L]]
        elements <- gsub('\\(', replacement='', x=elements)
        elements <- gsub('\\)', replacement='', x=elements)
        esplit <- strsplit(elements, ',')
        esplit <- lapply(esplit, function(x){
            newx <- c()
            if(length(x) < 3L)
                stop('START = ... has not been supplied enough arguments', call.=FALSE)
            for(i in seq_len(length(x)-2L)){
                if(grepl('-', x[i])){
                    tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                    newx <- c(newx, tmp[1L]:tmp[2L])
                } else newx <- c(newx, x[i])
            }
            x <- c(newx, x[length(x)-1L], x[length(x)])
            x
        })
        picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-2)]))
        for(i in seq_len(length(picks))){
            tmp <- ret[picks[[i]]]
            len <- length(esplit[[i]])
            tmp <- lapply(tmp, function(x, which, val){
                if(which %in% c('g', 'u')) val <- qlogis(val)
                x@par[names(x@parnum) == which] <- val
                x
            }, which=esplit[[i]][len-1L], val = as.numeric(esplit[[i]][len]))
            ret[picks[[i]]] <- tmp
        }
    }
    if(any(model[,1L] == 'FIXED')){
        input <- gsub(" ","", model[model[,1L] == 'FIXED', 2L])
        elements <- strsplit(input, '\\),\\(')[[1L]]
        elements <- gsub('\\(', replacement='', x=elements)
        elements <- gsub('\\)', replacement='', x=elements)
        esplit <- strsplit(elements, ',')
        esplit <- lapply(esplit, function(x){
            newx <- c()
            if(length(x) < 2L)
                stop('FIXED = ... has not been supplied enough arguments', call.=FALSE)
            for(i in seq_len(length(x)-1L)){
                if(grepl('-', x[i])){
                    tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                    newx <- c(newx, tmp[1L]:tmp[2L])
                } else newx <- c(newx, x[i])
            }
            x <- c(newx, x[length(x)])
            x
        })
        picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-1L)]))
        for(i in seq_len(length(picks))){
            tmp <- ret[picks[[i]]]
            len <- length(esplit[[i]])
            tmp <- lapply(tmp, function(x, which){
                x@est[names(x@parnum) == which] <- FALSE
                x
            }, which=esplit[[i]][len])
            ret[picks[[i]]] <- tmp
        }
    }
    if(any(model[,1L] == 'FREE')){
        input <- gsub(" ","", model[model[,1L] == 'FREE', 2L])
        elements <- strsplit(input, '\\),\\(')[[1L]]
        elements <- gsub('\\(', replacement='', x=elements)
        elements <- gsub('\\)', replacement='', x=elements)
        esplit <- strsplit(elements, ',')
        esplit <- lapply(esplit, function(x){
            newx <- c()
            if(length(x) < 2L)
                stop('FREE = ... has not been supplied enough arguments', call.=FALSE)
            for(i in seq_len(length(x)-1L)){
                if(grepl('-', x[i])){
                    tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                    newx <- c(newx, tmp[1L]:tmp[2L])
                } else newx <- c(newx, x[i])
            }
            x <- c(newx, x[length(x)])
            x
        })
        picks <- lapply(esplit, function(x) as.integer(x[seq_len(length(x)-1L)]))
        for(i in seq_len(length(picks))){
            tmp <- ret[picks[[i]]]
            len <- length(esplit[[i]])
            tmp <- lapply(tmp, function(x, which){
                x@est[names(x@parnum) == which] <- TRUE
                x
            }, which=esplit[[i]][len])
            ret[picks[[i]]] <- tmp
        }
    }
    ret[[length(ret) + 1L]] <- LoadGroupPars(gmeans=gmeans, gcov=gcov, estgmeans=estgmeans,
                                            estgcov=estgcov, parnumber=attr(ret, 'parnumber'),
                                            parprior=parprior, Rasch=all(itemtype %in% c('Rasch', 'rsm')),
                                            customGroup=customGroup)
    if(any(model[,1L] == 'LBOUND')){
        input <- gsub(" ","", model[model[,1L] == 'LBOUND', 2L])
        elements <- strsplit(input, '\\),\\(')[[1L]]
        elements <- gsub('\\(', replacement='', x=elements)
        elements <- gsub('\\)', replacement='', x=elements)
        esplit <- strsplit(elements, ',')
        esplit <- lapply(esplit, function(x){
            newx <- c()
            if(length(x) < 3L)
                stop('LBOUND = ... has not been supplied enough arguments', call.=FALSE)
            for(i in seq_len(length(x)-2L)){
                if(grepl('-', x[i])){
                    tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                    newx <- c(newx, tmp[1L]:tmp[2L])
                } else newx <- c(newx, x[i])
            }
            x <- c(newx, x[length(x)-1L], x[length(x)])
            x
        })
        picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-2)]))
        for(i in seq_len(length(picks))){
            tmp <- ret[picks[[i]]]
            len <- length(esplit[[i]])
            tmp <- lapply(tmp, function(x, which, val){
                if(which %in% c('g', 'u')) val <- qlogis(val)
                x@lbound[names(x@parnum) == which] <- val
                x
            }, which=esplit[[i]][len-1L], val = as.numeric(esplit[[i]][len]))
            ret[picks[[i]]] <- tmp
        }
    }
    if(any(model[,1L] == 'UBOUND')){
        input <- gsub(" ","", model[model[,1L] == 'UBOUND', 2L])
        elements <- strsplit(input, '\\),\\(')[[1L]]
        elements <- gsub('\\(', replacement='', x=elements)
        elements <- gsub('\\)', replacement='', x=elements)
        esplit <- strsplit(elements, ',')
        esplit <- lapply(esplit, function(x){
            newx <- c()
            if(length(x) < 3L)
                stop('UBOUND = ... has not been supplied enough arguments', call.=FALSE)
            for(i in seq_len(length(x)-2L)){
                if(grepl('-', x[i])){
                    tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                    newx <- c(newx, tmp[1L]:tmp[2L])
                } else newx <- c(newx, x[i])
            }
            x <- c(newx, x[length(x)-1L], x[length(x)])
            x
        })
        picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-2)]))
        for(i in seq_len(length(picks))){
            tmp <- ret[picks[[i]]]
            len <- length(esplit[[i]])
            tmp <- lapply(tmp, function(x, which, val){
                if(which %in% c('g', 'u')) val <- qlogis(val)
                x@ubound[names(x@parnum) == which] <- val
                x
            }, which=esplit[[i]][len-1L], val = as.numeric(esplit[[i]][len]))
            ret[picks[[i]]] <- tmp
        }
    }
    attr(ret, 'prodlist') <- prodlist
    attr(ret, 'exploratory') <- exploratory
    return(ret)
}
