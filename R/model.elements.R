model.elements <- function(model, factorNames, itemtype, nfactNames, nfact, J, K, fulldata,
                           itemloc, data, N, guess, upper, itemnames, exploratory, parprior,
                           parnumber, BFACTOR = FALSE, D, mixed.design, customItems, key,
                           nominal.highlow)
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
        for(j in 1L:length(prodlist)){
            for(i in 1L:nfact)
                prodlist[[j]][prodlist[[j]] == tmp2[[i]]] <- i
            prodlist[[j]] <- as.numeric(prodlist[[j]])
        }
    }
    #slopes specification
    estlam <- matrix(FALSE, ncol = nfactNames, nrow = J)
    for(i in 1L:nfactNames){
        tmp <- model[model[ ,1L] == factorNames[i],2L]
        if(any(regexpr(",",tmp)))
            tmp <- strsplit(tmp,",")[[1L]]
        popout <- c()
        for(j in 1L:length(tmp)){
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
    lambdas <- lambdas * 1.702/D
    zetas <- list()
    loc <- 1L
    for(i in 1L:J){
        if(K[i] == 2L){
            div <- ifelse(cs[i] > .25, cs[i], .25)
            zetas[[i]] <- (-1)*qnorm(mean(fulldata[,itemloc[i]]))/div * 1.702/D
        } else if(itemtype[i] %in% c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')){
            div <- ifelse(cs[i] > .25, cs[i], .25)
            zetas[[i]] <- qnorm(mean(fulldata[,itemloc[i] + key[i]-1L]))/div * 1.702/D
        } else {
            temp <- table(data[,i])[1L:(K[i]-1L)]/N
            temp <- cumsum(temp)
            div <- ifelse(cs[i] > .25, cs[i], .25)
            zetas[[i]] <- qnorm(1 - temp)/div * 1.702/D
        }
    }
    estzetas <- list()
    for(i in 1L:J)
        estzetas[[i]] <- length(zetas[[i]])
    #COV
    find <- 1L:nfact
    estgcov <- matrix(FALSE,nfact,nfact)
    if(any(model[,1L] == 'COV')){
        tmp <- model[model[,1L] == 'COV',2L]
        tmp <- strsplit(tmp,",")[[1L]]
        tmp <- gsub(" ","",tmp)
        for(i in 1L:length(tmp)){
            tmp2 <- strsplit(tmp[i],"*",fixed=TRUE)[[1L]]
            for(j in 1L:length(tmp2)){
                for(k in 1L:length(tmp2)){
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
        for(i in 1L:length(tmp))
            estgmeans[find[tmp[i] == factorNames]] <- TRUE
    }
    if(exploratory){
        Rpoly <- cormod(data, K, guess)
        loads <- eigen(Rpoly)$vector[,1L:nfact, drop = FALSE]
        u <- 1 - rowSums(loads^2)
        u[u < .001 ] <- .2
        cs <- sqrt(u)
        lambdas <- loads/cs * (1.702/D)
        if(!all(itemtype %in% c('lca', 'nlca')))
            lambdas[!estlam] <- 0
    }
    if(exploratory && any(itemtype %in% c('PC2PL', 'PC3PL')))
        stop('Partially compensatory models can only be estimated within a confirmatory model')
    ret <- LoadPars(itemtype=itemtype, itemloc=itemloc, lambdas=lambdas, zetas=zetas,
                    guess=guess, upper=upper, fulldata=fulldata, J=J, K=K,
                    nfact=nfact+length(prodlist), parprior=parprior,
                    parnumber=parnumber, estLambdas=estlam, BFACTOR=BFACTOR, D=D,
                    mixed.design=mixed.design, customItems=customItems, key=key,
                    nominal.highlow=nominal.highlow)
    if(any(model[,1L] == 'START')){
        start <- gsub(" ","", model[model[,1] == 'START', ])
        start <- strsplit(start[2], '),')[[1]]
        start <- sapply(start, function(x) gsub("\\(","", x))
        start <- sapply(start, function(x) gsub("\\)","", x))
        for(i in 1L:length(start)){
            tmp <- strsplit(start[i], ',')[[1L]]
            ret[[as.integer(tmp[1L])]]@par[tmp[2L]] <- as.numeric(tmp[3L])
        }
    }
    ret[[length(ret) + 1L]] <- LoadGroupPars(gmeans=gmeans, gcov=gcov, estgmeans=estgmeans,
                                            estgcov=estgcov, parnumber=attr(ret, 'parnumber'),
                                            parprior=parprior, Rasch=all(itemtype %in% c('Rasch', 'rsm')))
    attr(ret, 'prodlist') <- prodlist
    attr(ret, 'exploratory') <- exploratory
    return(ret)
}
