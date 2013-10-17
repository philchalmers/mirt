ESTIMATION <- function(data, model, group, itemtype = NULL, guess = 0, upper = 1,
                       invariance = '', pars = NULL, constrain = NULL, key = NULL,
                       parprior = NULL, mixed.design = NULL, customItems = NULL, 
                       nominal.highlow = NULL, ...)
{
    opts <- makeopts(...)
    if(!is.null(customItems)) opts$calcNull <- FALSE
    opts$start.time <- proc.time()[3L]
    # on exit, reset the seed to override internal
    if(opts$method == 'MHRM' || opts$method == 'MIXED')
        on.exit(set.seed((as.numeric(Sys.time()) - floor(as.numeric(Sys.time()))) * 1e8))
    #change itemtypes if NULL.MODEL
    if(opts$NULL.MODEL){
        constrain <- NULL
        if(!is.null(itemtype)){
            itemtype[itemtype == 'grsm'] <- 'graded'
            itemtype[itemtype == 'rsm'] <- 'gpcm'
            itemtype[itemtype == '3PL' | itemtype == '3PLu' | itemtype == '4PL'] <- '2PL'
            itemtype[itemtype == '3PLNRM' | itemtype == '3PLuNRM' | itemtype == '4PLNRM'] <- '2PLNRM'
        }
    }
    if(length(group) != nrow(data))
        stop('length of group not equal to number of rows in data.')
    if(any(is.na(group))){
        stop('Unknown group memberships cannot be estimated. Please remove the NA values in group
                and the associated rows in the data input.')
    }
    
    ##
    Data <- list()
    data <- as.matrix(data)
    if(!all(apply(data, 2L, class) %in% c('integer', 'numeric')))
        stop('Input data must be integer or numeric values only')
    if(is.numeric(data))
        data <- matrix(as.integer(data), nrow(data), dimnames=list(rownames(data), colnames(data)))
    rownames(data) <- 1L:nrow(data)
    if(is.null(colnames(data)))
        colnames(data) <- paste0('Item.', 1L:ncol(data))
    Data$data <- data
    if(is.null(opts$grsm.block)) Data$grsm.block <- rep(1L, ncol(data))
    if(is.null(opts$rsm.block)) Data$rsm.block <- rep(1L, ncol(data))
    Data$group <- factor(group)
    Data$groupNames <- unique(Data$group)
    if(any(grepl('-', Data$groupNames)))
        stop('Group names cannot contain a dash (-) character')
    Data$ngroups <- length(Data$groupNames)
    Data$nitems <- ncol(Data$data)
    Data$N <- nrow(Data$data)
    oldmodel <- model    
    if(length(model) == 1L){
        newmodel <- list()
        for(g in 1L:Data$ngroups)
            newmodel[[g]] <- model
        names(newmodel) <- Data$groupNames
        model <- newmodel
    }
    Data$model <- model
    PrepList <- vector('list', Data$ngroups)
    names(PrepList) <- Data$groupNames
    tmp <- 1L:Data$ngroups
    selectmod <- Data$model[[tmp[names(Data$model) == Data$groupNames[1L]]]]    
    PrepListFull <- PrepList[[1L]] <-
        PrepData(data=Data$data, model=selectmod, itemtype=itemtype, guess=guess,
                 upper=upper, parprior=parprior, verbose=opts$verbose,
                 technical=opts$technical, parnumber=1L, BFACTOR=opts$BFACTOR,
                 grsm.block=Data$grsm.block, rsm.block=Data$rsm.block,
                 D=opts$D, mixed.design=mixed.design, customItems=customItems,
                 fulldata=opts$PrepList[[1L]]$fulldata, key=key, nominal.highlow=nominal.highlow)
    if(any(PrepListFull$itemtype == 'nominal') && is.null(nominal.highlow) && !opts$NULL.MODEL
       && opts$verbose)
        message('\'nominal.highlow\' matrix not specified, highest and lowest categories are used by default')    
    parnumber <- max(PrepList[[1L]]$pars[[Data$nitems+1L]]@parnum) + 1L
    for(g in 1L:Data$ngroups){
        if(g != 1L){
            PrepList[[g]] <- PrepList[[1L]]
            for(i in 1L:length(PrepList[[g]]$pars)){
                PrepList[[g]]$pars[[i]]@parnum <- parnumber:(parnumber + length(PrepList[[g]]$pars[[i]]@parnum) - 1L)
                parnumber <- max(PrepList[[g]]$pars[[i]]@parnum) + 1L
            }
        }
    }
    if(length(mixed.design$random) > 0L){
        for(i in 1L:length(mixed.design$random)){
            mixed.design$random[[i]]@parnum <- parnumber:(parnumber - 1L +  
                                                              length(mixed.design$random[[i]]@par))
            parnumber <- max(mixed.design$random[[i]]@parnum) + 1L            
        }
    }
    if(!is.null(opts$PrepList)){
        PrepList <- opts$PrepList
    } else {
        tmpdata <- Data$data
        tmpdata[is.na(tmpdata)] <- 99999L
        stringfulldata <- apply(tmpdata, 1L, paste, sep='', collapse = '/')
        stringtabdata <- unique(stringfulldata)
        tmptabdata <- maketabData(stringfulldata=stringfulldata, stringtabdata=stringtabdata,
                                  group=Data$group, groupNames=Data$groupNames, nitem=Data$nitems,
                                  K=PrepListFull$K, itemloc=PrepListFull$itemloc,
                                  Names=PrepListFull$Names, itemnames=PrepListFull$itemnames)
        for(g in 1L:Data$ngroups){
            select <- Data$group == Data$groupNames[g]
            for(i in 1L:Data$nitems)
                PrepList[[g]]$pars[[i]]@dat <- 
                    PrepListFull$fulldata[select, PrepListFull$itemloc[i]:(PrepListFull$itemloc[i+1L]-1L), drop = FALSE]
            PrepList[[g]]$fulldata <- PrepListFull$fulldata[select, ]
            PrepList[[g]]$tabdata <- tmptabdata$tabdata[[g]]
            PrepList[[g]]$tabdata2 <- tmptabdata$tabdata2[[g]]
        }
        rm(tmpdata, tmptabdata, stringfulldata, stringtabdata, select, parnumber,
           PrepListFull, selectmod)
    }
    if(opts$returnPrepList) return(PrepList)
    if(opts$BFACTOR){
        #better start values
        if((PrepList[[1L]]$nfact - attr(model[[1L]], 'nspec')) == 1L){
            Rpoly <- cormod(Data$data, PrepList[[1L]]$K, guess, use=opts$use)
            loads <- abs(eigen(Rpoly)$vector[,1L, drop = FALSE])
            u <- 1 - rowSums(loads^2)
            u[u < .001] <- .2
            cs <- sqrt(u)
            astart <- loads/cs
            astart <- cbind(astart,astart/2L) * 1.702
            nfact <- PrepList[[1L]]$pars[[1L]]@nfact
            for(g in 1L:Data$ngroups)
                for(i in 1L:Data$nitems)
                    PrepList[[g]]$pars[[i]]@par[PrepList[[g]]$pars[[i]]@est][1L:2L] <- astart[i, ]
            rm(Rpoly, loads, u, cs, astart)
        }
    }
    if(!is.null(pars)){        
        if(is(pars, 'data.frame')){
            PrepList <- UpdatePrepList(PrepList, pars, random=mixed.design$random, MG = TRUE)            
            mixed.design$random <- attr(PrepList, 'random')
            attr(PrepList, 'random') <- NULL
        }          
        if(!is.null(attr(pars, 'values'))){
            return(ReturnPars(PrepList, PrepList[[1L]]$itemnames, 
                              random=mixed.design$random, MG = TRUE))
        } else if(is.character(pars)){
            if(pars == 'values') 
                return(ReturnPars(PrepList, PrepList[[1L]]$itemnames, 
                                  random=mixed.design$random, MG = TRUE))
        } 
    }
    pars <- vector('list', Data$ngroups)
    for(g in 1L:Data$ngroups)
        pars[[g]] <- PrepList[[g]]$pars
    nitems <- Data$nitems
    K <- PrepList[[1L]]$K
    Data$nfact <- nfact <- PrepList[[1L]]$pars[[nitems+1L]]@nfact
    if(nfact != 1L && any(c('Rasch') %in% itemtype ) && PrepList[[1L]]$exploratory)
       stop('Rasch itemtypes are for confimatory models only.')
    nLambdas <- PrepList[[1L]]$pars[[1L]]@nfact
    if(is.null(constrain)) constrain <- list()
    #default MG uses configural model (independent groups but each identified)
    if('free_means' %in% invariance ){ #Free factor means (means 0 for ref)
        for(g in 2L:Data$ngroups)
            pars[[g]][[nitems + 1L]]@est[1L:nfact] <- TRUE
    }
    dummymat <- matrix(FALSE, pars[[1L]][[nitems + 1L]]@nfact, pars[[1L]][[nitems + 1L]]@nfact)
    if(any(c('free_cov', 'free_varcov') %in% invariance)){ #Free factor covs (vars 1 for ref)
        dummymat <- matrix(TRUE, pars[[1L]][[nitems + 1L]]@nfact, pars[[1L]][[nitems + 1L]]@nfact)
        diag(dummymat) <- FALSE        
    }
    if(any(c('free_var', 'free_varcov') %in% invariance)) #Free factor vars (vars 1 for ref)
        diag(dummymat) <- TRUE
    if(any(c('free_var', 'free_varcov', 'free_cov') %in% invariance)){
        tmp <- dummymat[lower.tri(dummymat, TRUE)]
        for(g in 2L:Data$ngroups){
            pars[[g]][[nitems + 1L]]@est <- c(pars[[g]][[nitems + 1L]]@est[1L:pars[[g]][[nitems + 1L]]@nfact], tmp)
            names(pars[[g]][[nitems + 1L]]@est) <- names(pars[[g]][[nitems + 1L]]@par)
        }
    }
    constrain <- UpdateConstrain(pars=pars, constrain=constrain, invariance=invariance, nfact=Data$nfact,
                                 nLambdas=nLambdas, J=nitems, ngroups=Data$ngroups, PrepList=PrepList,
                                 method=opts$method, itemnames=PrepList[[1L]]$itemnames, model=model,
                                 groupNames=Data$groupNames)
    startlongpars <- c()
    if(opts$NULL.MODEL){
        constrain <- list()
        for(i in 1L:nitems){
            pars[[1L]][[i]]@par[1L] <- 0
            pars[[1L]][[i]]@est[1L] <- FALSE
            if(is(pars[[1L]][[i]], 'nominal'))
                pars[[1L]][[i]]@est[(nfact+1L):(nfact + K[i])] <- FALSE            
            if(is(pars[[1L]][[i]], 'nestlogit'))
                pars[[1L]][[i]]@est[(nfact+5L):(nfact + K[i] + 1L)] <- FALSE
        }
    }
    df <- rr <- 0L
    for(g in 1L:Data$ngroups){
        r <- PrepList[[g]]$tabdata
        r <- r[, ncol(r)]
        rr <- rr + r
        df <- df + sum(r != 0) - 1L
    }
    if(!is.null(opts$technical$override.df)){
        df <- 1e10
        opts$calcNull <- FALSE        
    }
    nestpars <- nconstr <- 0L
    for(g in 1L:Data$ngroups)
        for(i in 1L:(nitems+1L))
            nestpars <- nestpars + sum(pars[[g]][[i]]@est)
    if(!is.null(mixed.design$random)){
        for(i in 1L:length(mixed.design$random))
            nestpars <- nestpars + sum(mixed.design$random[[i]]@est)        
    }
    if(length(constrain) > 0L)
        for(i in 1L:length(constrain))
            nconstr <- nconstr + length(constrain[[i]]) - 1L
    nmissingtabdata <- sum(is.na(rowSums(PrepList[[1L]]$tabdata2)))
    dfsubtr <- nestpars - nconstr
    if(PrepList[[1L]]$exploratory) dfsubtr <- dfsubtr - nfact*(nfact - 1L)/2L
    if(df <= dfsubtr) 
        stop('Too few degrees of freedom. There are only ', df, ' degrees of freedom but ', 
             dfsubtr, ' parameters were freely estimated.')
    df <- df - dfsubtr
    if(!is.null(customItems)){
        for(g in 1L:Data$ngroups)
            PrepList[[g]]$exploratory <- FALSE
    }
    G2group <- X2group <- numeric(Data$ngroups)
    #avoid some S4 overhead by caching functions....as if this works
    DERIV <- vector('list', Data$ngroups)
    for(g in 1L:Data$ngroups){
        DERIV[[g]] <- vector('list', Data$nitems)
        for(i in 1L:Data$nitems)            
            DERIV[[g]][[i]] <- selectMethod(Deriv, c(class(pars[[g]][[i]]), 'matrix'))
    }    
    
    #EM estimation    
    if(opts$method == 'EM'){
        if(is.null(opts$quadpts))
            opts$quadpts <- switch(as.character(nfact), '1'=41, '2'=21, '3'=11, '4'=7, '5'=5, 3)
        if(opts$quadpts %% 2 == 0) stop('Must use an odd number for quadpts')
        if(opts$quadpts < 3) stop('Must use more than 2 quadpts')
        Theta <- theta <- as.matrix(seq(-(.8 * sqrt(opts$quadpts)), .8 * sqrt(opts$quadpts),
                                        length.out = opts$quadpts))
        nspec <- ifelse(!is.null(attr(model[[1L]], 'nspec')), attr(model[[1L]], 'nspec'), 1L)
        temp <- matrix(0L,nrow=nitems,ncol=nspec)
        sitems <- matrix(0L, nrow=sum(PrepList[[1L]]$K), ncol=nspec)
        specific <- NULL
        if(opts$BFACTOR){
            specific <- attr(oldmodel, 'specific')
            specific[is.na(specific)] <- 1L
            for(i in 1L:nitems) temp[i, specific[i]] <- 1L
            ind <- 1L
            for(i in 1L:nitems){
                for(j in 1L:PrepList[[1L]]$K[i]){
                    sitems[ind, ] <- temp[i, ]
                    ind <- ind + 1L
                }
            }
            theta <- seq(-4, 4, length.out = opts$quadpts)
            tmp <- PrepList[[1L]]$nfact - attr(model[[1L]], 'nspec') + 1L
            Theta <- thetaComb(theta, tmp)
            Theta <- cbind(Theta[,1L:(tmp-1L),drop=FALSE], 
                           matrix(Theta[,tmp], nrow=nrow(Theta), ncol=ncol(sitems)))
        } else {
            if(opts$quadpts^nfact <= opts$MAXQUAD){
                Theta <- thetaComb(theta, nfact)
            } else stop('Greater than ', opts$MAXQUAD, ' quadrature points.')
        }
        ESTIMATE <- EM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                             list = list(NCYCLES=opts$NCYCLES, TOL=opts$TOL, MSTEPTOL=opts$MSTEPTOL,
                                         nfactNames=PrepList[[1L]]$nfactNames, theta=theta,
                                         itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                         sitems=sitems, specific=specific, NULL.MODEL=opts$NULL.MODEL,
                                         nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                         SEM=opts$SE.type == 'SEM' && opts$SE, accelerate=opts$accelerate),
                             Theta=Theta, DERIV=DERIV)
        startlongpars <- ESTIMATE$longpars
        rlist <- ESTIMATE$rlist
        logLik <- G2 <- X2 <- SElogLik <- 0
        for(g in 1L:Data$ngroups){
            Pl <- rlist[[g]]$expected
            rg <- PrepList[[g]]$tabdata[,ncol(PrepList[[g]]$tabdata)]
            Pl <- Pl[rg != 0]
            rg <- rg[rg != 0]
            Ng <- sum(rg)
            G2group[g] <- 2 * sum(rg * log(rg/(Ng*Pl)))
            G2 <- G2 + G2group[g]
            X2group[g] <- sum((rg - Ng*Pl)^2 / (Ng*Pl))
            X2 <- X2 + X2group[g]
            logLik <- logLik + sum(rg*log(Pl))
        }
        Pl <- list(Pl)
        if(!opts$NULL.MODEL && opts$SE){
            tmp <- ESTIMATE
            if(opts$verbose) cat('\n\nCalculating information matrix...\n')
            if(opts$SE.type == 'SEM'){
                dontrun <- FALSE
                if(ESTIMATE$cycles <= 5L){
                    dontrun <- TRUE
                    warning('Too few EM interations to compute SEM information matrix')
                }
                if(ESTIMATE$cycles == opts$NCYCLES) dontrun <- TRUE
                lengthsplit <- do.call(c, lapply(strsplit(names(ESTIMATE$correct), 'COV_'), length))
                if(any(lengthsplit > 1L)){
                    warning('Group covariance parameters are not well estimated by S-EM algorithm. 
                            Better to use SE.type = \'MHRM\', SE.type = \'BL\', or resort to the 
                            boot.mirt() and PLCI.mirt() functions')
                }
                if(!dontrun){
                    if(ESTIMATE$cycles <= 10L)
                        message('Very few EM cycles performed. Consider decreasing TOL further to
                                increase EM iteration count')
                    estmat <- matrix(FALSE, length(ESTIMATE$correction), length(ESTIMATE$correction))
                    DM <- estmat + 0
                    diag(estmat) <- TRUE
                    if(!is.null(globalenv()$MIRTCLUSTER)){
                        DM <- t(parallel::parApply(cl=globalenv()$MIRTCLUSTER, estmat, MARGIN=1L, FUN=SEM.SE,
                                pars=ESTIMATE$pars, constrain=constrain, PrepList=PrepList,
                                list = list(NCYCLES=opts$NCYCLES, TOL=opts$TOL, MSTEPTOL=opts$MSTEPTOL,
                                           nfactNames=PrepList[[1L]]$nfactNames, theta=theta,
                                           itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                           sitems=sitems, specific=oldmodel, NULL.MODEL=opts$NULL.MODEL,
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose),
                                                   Theta=Theta, theta=theta, ESTIMATE=ESTIMATE, 
                                                   DERIV=DERIV))
                    } else {
                        for(i in 1L:ncol(DM))
                            DM[i, ] <- SEM.SE(est=estmat[i,], pars=ESTIMATE$pars, constrain=constrain,
                                              PrepList=PrepList,
                                          list = list(NCYCLES=opts$NCYCLES, TOL=opts$TOL, MSTEPTOL=opts$MSTEPTOL,
                                                      nfactNames=PrepList[[1L]]$nfactNames, theta=theta,
                                                      itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                                      sitems=sitems, specific=oldmodel, NULL.MODEL=opts$NULL.MODEL,
                                                      nfact=nfact, constrain=constrain, verbose=opts$verbose),
                                              Theta=Theta, theta=theta, ESTIMATE=ESTIMATE, 
                                              DERIV=DERIV)
                    }
                    ESTIMATE$pars <- reloadPars(longpars=ESTIMATE$longpars, pars=ESTIMATE$pars, 
                                                ngroups=Data$ngroups, J=Data$nitems)
                    info <- solve(-solve(ESTIMATE$hess) %*% solve(diag(ncol(DM)) - DM))
                    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain)
                }
            }
            if(opts$SE.type == 'BL')
                ESTIMATE <- BL.SE(pars=ESTIMATE$pars, Theta=Theta, theta=theta, PrepList=PrepList,
                                  BFACTOR=opts$BFACTOR, itemloc=PrepList[[1L]]$itemloc, ESTIMATE=ESTIMATE,
                                  constrain=constrain, specific=oldmodel, sitems=sitems)
            if(opts$SE.type == 'MHRM')
                ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=1L, SEMCYCLES=5L,
                                           KDRAWS=opts$KDRAWS, TOL=opts$SEtol, USEEM=opts$USEEM,
                                           gain=opts$gain,
                                           nfactNames=PrepList[[1L]]$nfactNames,
                                           itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                           nfact=nfact, constrain=constrain, verbose=FALSE,
                                           startlongpars=startlongpars),
                                       DERIV=DERIV)
            ESTIMATE$cycles <- tmp$cycles
        }
    } else if(opts$method == 'MHRM'){ #MHRM estimation
        Theta <- matrix(0, Data$N, nitems)
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=opts$BURNIN,
                                           SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                           KDRAWS=opts$KDRAWS, TOL=opts$TOL, USEEM=FALSE,
                                           nfactNames=PrepList[[1L]]$nfactNames,
                                           itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                           startlongpars=startlongpars),
                               DERIV=DERIV)
        rlist <- vector('list', Data$ngroups)
        for(g in 1L:Data$ngroups)
            rlist[[g]]$expected = numeric(1L)
    } else if(opts$method == 'MIXED'){
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain,
                                    PrepList=PrepList, random=mixed.design$random,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=opts$BURNIN,
                                           SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                           KDRAWS=opts$KDRAWS, TOL=opts$TOL, USEEM=FALSE,
                                           nfactNames=PrepList[[1L]]$nfactNames,
                                           itemloc=PrepList[[1L]]$itemloc, BFACTOR=opts$BFACTOR,
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                           startlongpars=startlongpars),
                               DERIV=DERIV)
        rlist <- vector('list', Data$ngroups)
        for(g in 1L:Data$ngroups)
            rlist[[g]]$expected = numeric(1L)
    }
    cmods <- list()
    for(g in 1L:Data$ngroups){
        lambdas <- Lambdas(ESTIMATE$pars[[g]]) * opts$D/1.702
        if (ncol(lambdas) > 1L) norm <- sqrt(1 + rowSums(lambdas^2))
        else norm <- as.matrix(sqrt(1 + lambdas[ ,1L]^2))
        alp <- as.matrix(lambdas/norm)
        F <- alp
        if(opts$method != 'MIXED')
            colnames(F) <- PrepList[[g]]$factorNames
        h2 <- rowSums(F^2)
        cmods[[g]] <- new('ConfirmatoryClass', pars=ESTIMATE$pars[[g]], itemloc=PrepList[[g]]$itemloc,
                          tabdata=PrepList[[g]]$tabdata2, data=Data$data[group == Data$groupNames[[g]], ],
                          converge=ESTIMATE$converge, esttype='MHRM', F=F, h2=h2, prodlist=PrepList[[g]]$prodlist,
                          K=PrepList[[g]]$K, tabdatalong=PrepList[[g]]$tabdata, nfact=nfact,
                          constrain=constrain, G2=G2group[g], X2=X2group[g], Pl = rlist[[g]]$expected,                          
                          fulldata=PrepList[[g]]$fulldata, factorNames=PrepList[[g]]$factorNames, 
                          random=ESTIMATE$random)
    }
    rm(lambdas, norm)
    #missing stats for MHRM
    if(opts$method =='MHRM' || opts$method == 'MIXED'){
        if(opts$verbose) cat("\nCalculating log-likelihood...\n")
        flush.console()
        logLik <- G2 <- X2 <- SElogLik <- 0
        Pl <- list()
        for(g in 1L:Data$ngroups){
            cmods[[g]] <- calcLogLik(cmods[[g]], opts$draws, G2 = 'return')
            logLik <- logLik + cmods[[g]]@logLik
            SElogLik <- SElogLik + cmods[[g]]@SElogLik
            G2 <- G2 + cmods[[g]]@G2
            Pl[[g]] <- cmods[[g]]@Pl
            X2 <- X2 + cmods[[g]]@X2
        }
    }

    ####post estimation stats    
    r <- rr
    N <- sum(r)
    logN <- 0
    logr <- rep(0,length(r))
    for (i in 1L:N) logN <- logN + log(i)
    for (i in 1L:length(r))
        for (j in 1L:r[i])
            logr[i] <- logr[i] + log(j)
    if(sum(logr) != 0)
        logLik <- logLik + logN/sum(logr)
    tmp <- (length(r) - df - 1L)
    AIC <- (-2) * logLik + 2 * tmp
    AICc <- AIC + 2 * tmp * (tmp + 1) / (length(r) - tmp - 1L)
    BIC <- (-2) * logLik + tmp*log(N)
    SABIC <- (-2) * logLik + tmp*log((N+2)/24)
    p.G2 <- 1 - pchisq(G2,df)
    p.X2 <- 1 - pchisq(X2,df)
    RMSEA.G2 <- ifelse((G2 - df) > 0,
                    sqrt(G2 - df) / sqrt(df * (N-1)), 0)
    RMSEA.X2 <- ifelse((X2 - df) > 0,
                       sqrt(X2 - df) / sqrt(df * (N-1)), 0)
    null.mod <- unclass(new('ConfirmatoryClass'))
    TLI.G2 <- TLI.X2 <- CFI.G2 <- CFI.X2 <- NaN
    if(!opts$NULL.MODEL && opts$method != 'MIXED' && opts$calcNull && nmissingtabdata == 0L){
        null.mod <- try(unclass(mirt(data, 1, itemtype=itemtype, technical=list(NULL.MODEL=TRUE, TOL=1e-3),
                                     large=opts$PrepList, key=key, verbose=FALSE)))
        if(is(null.mod, 'try-error')){
            message('Null model calculation did not converge.')
            null.mod <- unclass(new('ConfirmatoryClass'))
        } else {
            TLI.G2 <- (null.mod@G2 / null.mod@df - G2/df) / (null.mod@G2 / null.mod@df - 1)
            TLI.X2 <- (null.mod@X2 / null.mod@df - X2/df) / (null.mod@X2 / null.mod@df - 1)
            CFI.G2 <- 1 - (G2 - df) / (null.mod@G2 - null.mod@df)
            CFI.X2 <- 1 - (X2 - df) / (null.mod@X2 - null.mod@df)
            CFI.G2 <- ifelse(CFI.G2 > 1, 1, CFI.G2)
            CFI.G2 <- ifelse(CFI.G2 < 0, 0, CFI.G2)
            CFI.X2 <- ifelse(CFI.X2 > 1, 1, CFI.X2)
            CFI.X2 <- ifelse(CFI.X2 < 0, 0, CFI.X2)
        }
    }
    if(nmissingtabdata > 0L)
        p.G2 <- p.X2 <- RMSEA.G2 <- RMSEA.X2 <- G2 <- X2 <- TLI.G2 <-
            TLI.X2 <- CFI.G2 <- CFI.X2 <- NaN
    if(!is.nan(G2) && !opts$NULL.MODEL)
        if(X2/G2 > 10) TLI.X2 <- CFI.X2 <- X2 <- p.X2 <- RMSEA.X2 <- NaN
    if(is.null(parprior)) parprior <- list()
    if(is.null(opts$quadpts)) opts$quadpts <- NaN
    if(Data$ngroups == 1L){
        if(opts$method == 'MIXED'){            
            mod <- new('MixedClass',
                       iter=ESTIMATE$cycles,
                       pars=cmods[[1L]]@pars,
                       model=list(oldmodel),
                       df=df,
                       itemloc=PrepList[[1L]]$itemloc,
                       method=opts$method,
                       AIC=AIC,
                       AICc=AICc,
                       BIC=BIC,
                       SABIC=SABIC,
                       logLik=logLik,
                       F=F,
                       h2=h2,
                       tabdata=PrepList[[1L]]$tabdata2,
                       Pl=Pl[[1L]],
                       data=Data$data,
                       converge=ESTIMATE$converge,
                       nfact=nfact,
                       K=PrepList[[1L]]$K,                       
                       tabdatalong=PrepList[[1L]]$tabdata,                       
                       factorNames=PrepList[[1L]]$factorNames,
                       constrain=constrain,
                       parprior=parprior,
                       fulldata=PrepList[[1L]]$fulldata,
                       itemtype=PrepList[[1L]]$itemtype,
                       random=ESTIMATE$random,
                       cand.t.var=ESTIMATE$cand.t.var,
                       information=ESTIMATE$info)
        } else if(PrepList[[1L]]$exploratory){
            FF <- alp %*% t(alp)
            V <- eigen(FF)$vector[ ,1L:nfact]
            L <- eigen(FF)$values[1L:nfact]
            if (nfact == 1L) F <- as.matrix(V * sqrt(L))
            else F <- V %*% sqrt(diag(L))
            if (sum(F[ ,1L] < 0)) F <- (-1) * F
            colnames(F) <- paste("F_", 1:ncol(F),sep="")
            h2 <- rowSums(F^2)
            mod <- new('ExploratoryClass', iter=ESTIMATE$cycles,
                       pars=cmods[[1L]]@pars,
                       model=list(oldmodel),
                       G2=G2,
                       X2=X2,
                       p=p.G2,
                       p.X2=p.X2,
                       TLI=TLI.G2,
                       TLI.X2=TLI.X2,
                       CFI=CFI.G2,
                       CFI.X2=CFI.X2,
                       RMSEA=RMSEA.G2,
                       RMSEA.X2=RMSEA.X2,
                       df=df,
                       itemloc=PrepList[[1L]]$itemloc,
                       method=opts$method,
                       AIC=AIC,
                       AICc=AICc,
                       BIC=BIC,
                       SABIC=SABIC,
                       logLik=logLik,
                       F=F,
                       h2=h2,
                       tabdata=PrepList[[1L]]$tabdata2,
                       Theta=Theta,
                       Pl=Pl[[1L]],
                       data=Data$data,
                       converge=ESTIMATE$converge,
                       nfact=nfact,
                       quadpts=opts$quadpts,
                       K=PrepList[[1L]]$K,
                       tabdatalong=PrepList[[1L]]$tabdata,
                       rotate=opts$rotate,
                       null.mod=null.mod,
                       Target=opts$Target,
                       factorNames=PrepList[[1L]]$factorNames,
                       constrain=constrain,
                       parprior=parprior,
                       fulldata=PrepList[[1L]]$fulldata,
                       itemtype=PrepList[[1L]]$itemtype,
                       information=ESTIMATE$info)
        } else {
            mod <- new('ConfirmatoryClass', iter=ESTIMATE$cycles,
                       pars=cmods[[1L]]@pars,
                       model=list(oldmodel),
                       G2=G2,
                       X2=X2,
                       p=p.G2,
                       p.X2=p.X2,
                       TLI=TLI.G2,
                       TLI.X2=TLI.X2,
                       CFI=CFI.G2,
                       CFI.X2=CFI.X2,
                       RMSEA=RMSEA.G2,
                       RMSEA.X2=RMSEA.X2,
                       df=df,
                       itemloc=PrepList[[1L]]$itemloc,
                       AIC=AIC,
                       AICc=AICc,
                       BIC=BIC,
                       SABIC=SABIC,
                       logLik=logLik,
                       F=F,
                       h2=h2,
                       tabdata=PrepList[[1L]]$tabdata2,
                       Theta=Theta,
                       method=opts$method,
                       Pl=Pl[[1L]],
                       data=Data$data,
                       converge=ESTIMATE$converge,
                       nfact=nfact,
                       quadpts=opts$quadpts,
                       K=PrepList[[1L]]$K,
                       tabdatalong=PrepList[[1L]]$tabdata,
                       null.mod=null.mod,
                       factorNames=PrepList[[1L]]$factorNames,
                       constrain=constrain,
                       parprior=parprior,
                       fulldata=PrepList[[1L]]$fulldata,
                       itemtype=PrepList[[1L]]$itemtype,
                       information=ESTIMATE$info)
        }
    } else {
        tabdatalong <- PrepList[[1L]]$tabdata
        tabdata <- PrepList[[1L]]$tabdata2
        tabdata[,ncol(tabdata)] <- tabdatalong[,ncol(tabdatalong)] <- r
        if(is.nan(X2))
            for(g in 1L:length(cmods))
                cmods[[g]]@X2 <- NaN
        mod <- new('MultipleGroupClass', iter=ESTIMATE$cycles,
                   cmods=cmods,
                   model=list(oldmodel),
                   itemloc=PrepList[[1L]]$itemloc,
                   tabdata=tabdata,
                   data=Data$data,
                   converge=ESTIMATE$converge,
                   esttype=opts$method,
                   K=PrepList[[1L]]$K,
                   tabdatalong=tabdatalong,
                   constrain=constrain,
                   parprior=parprior,
                   group=Data$group,
                   groupNames=Data$groupNames,
                   invariance=invariance,
                   df=as.integer(df),
                   logLik=logLik,
                   method=opts$method,
                   SElogLik=SElogLik,
                   AIC=AIC,
                   AICc=AICc,
                   BIC=BIC,
                   SABIC=SABIC,
                   nfact=nfact,
                   G2=G2,
                   X2=X2,
                   p=p.G2,
                   p.X2=p.X2,
                   TLI=TLI.G2,
                   TLI.X2=TLI.X2,
                   CFI=CFI.G2,
                   CFI.X2=CFI.X2,
                   RMSEA=RMSEA.G2,
                   RMSEA.X2=RMSEA.X2,
                   Theta=Theta,
                   Pl=Pl,
                   itemtype=PrepList[[1L]]$itemtype,
                   information=ESTIMATE$info)
    }
    mod@time <- proc.time()[3L] - opts$start.time
    return(mod)
}
