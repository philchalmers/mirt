ESTIMATION <- function(data, model, group, itemtype = NULL, guess = 0, upper = 1, 
                       invariance = '', pars = NULL, constrain = NULL, key = NULL,
                       parprior = NULL, mixedlist = NULL, customItems = NULL, ...)
{       
    opts <- makeopts(...)
    if(!is.null(customItems)) opts$calcNull <- FALSE
    opts$start.time <- proc.time()[3]                 
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
    ##	    
    Data <- list()
    data <- as.matrix(data)
    rownames(data) <- 1:nrow(data)
    Data$data <- data
    if(is.null(opts$grsm.block)) Data$grsm.block <- rep(1, ncol(data))    
    if(is.null(opts$rsm.block)) Data$rsm.block <- rep(1, ncol(data))
    Data$group <- factor(group)
    Data$groupNames <- unique(Data$group)
    Data$ngroups <- length(Data$groupNames)
    Data$nitems <- ncol(Data$data)
    Data$N <- nrow(Data$data)
    oldmodel <- model
    if(is(model, 'numeric') && length(model) > 1)
        model <- bfactor2mod(model, ncol(data))
    if(length(model) == 1){
        newmodel <- list()
        for(g in 1:Data$ngroups)
            newmodel[[g]] <- model
        names(newmodel) <- Data$groupNames
        model <- newmodel
    }
    Data$model <- model    
    PrepList <- vector('list', Data$ngroups)        
    names(PrepList) <- Data$groupNames    
    tmp <- 1:Data$ngroups
    selectmod <- Data$model[[tmp[names(Data$model) == Data$groupNames[1]]]]         
    PrepListFull <- PrepList[[1]] <- 
        PrepData(data=Data$data, model=selectmod, itemtype=itemtype, guess=guess, 
                 upper=upper, parprior=parprior, verbose=opts$verbose,  
                 technical=opts$technical, parnumber=1, BFACTOR=opts$BFACTOR,
                 grsm.block=Data$grsm.block, rsm.block=Data$rsm.block, 
                 D=opts$D, mixedlist=mixedlist, customItems=customItems,
                 fulldata=opts$PrepList[[1]]$fulldata, key=key)                    
    parnumber <- 1
    for(g in 1:Data$ngroups){                    
        tmp <- 1:Data$ngroups
        selectmod <- Data$model[[tmp[names(Data$model) == Data$groupNames[g]]]]
        if(g != 1)
            PrepList[[g]] <- PrepData(data=Data$data, model=selectmod, itemtype=itemtype, guess=guess, 
                                      upper=upper, parprior=parprior, verbose=opts$verbose, 
                                      technical=opts$technical, parnumber=parnumber, BFACTOR=opts$BFACTOR,
                                      grsm.block=opts$grsm.block, D=opts$D, mixedlist=mixedlist, 
                                      customItems=customItems, fulldata=PrepList[[1]]$fulldata, key=key)        
        tmp <- PrepList[[g]]$pars[[length(PrepList[[g]]$pars)]]
        parnumber <- tmp@parnum[length(tmp@parnum)] + 1
    }
    if(!is.null(opts$PrepList)){
        PrepList <- opts$PrepList
    } else {                
        tmpdata <- Data$data
        tmpdata[is.na(tmpdata)] <- 99999
        stringfulldata <- apply(tmpdata, 1, paste, sep='', collapse = '/')
        stringtabdata <- unique(stringfulldata)          
        tmptabdata <- maketabData(stringfulldata=stringfulldata, stringtabdata=stringtabdata,                               
                                  group=Data$group, groupNames=Data$groupNames, nitem=Data$nitems, 
                                  K=PrepListFull$K, itemloc=PrepListFull$itemloc, 
                                  Names=PrepListFull$Names, itemnames=PrepListFull$itemnames)
        for(g in 1:Data$ngroups){              
            select <- Data$group == Data$groupNames[g]
            for(i in 1:Data$nitems)
                PrepList[[g]]$pars[[i]]@dat <- PrepList[[g]]$pars[[i]]@dat[select, , drop = FALSE]
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
        Rpoly <- cormod(Data$data, PrepList[[1]]$K, guess, use=opts$use)        
        loads <- abs(eigen(Rpoly)$vector[,1, drop = FALSE])
        u <- 1 - rowSums(loads^2)       
        u[u < .001 ] <- .2
        cs <- sqrt(u)              
        astart <- loads/cs
        astart <- cbind(astart,astart/2)* (1.702 / opts$D) #reweight due to D
        nfact <- PrepList[[1]]$pars[[1]]@nfact
        for(g in 1:Data$ngroups)
            for(i in 1:Data$nitems)
                PrepList[[g]]$pars[[i]]@par[PrepList[[g]]$pars[[i]]@est][1:2] <- astart[i, ] 
        rm(Rpoly, loads, u, cs, astart) 
    }    
    if(!is.null(pars)){
        if(is(pars, 'data.frame')){
            PrepList <- UpdatePrepList(PrepList, pars, MG = TRUE)
        } else if(pars == 'values'){            
            return(ReturnPars(PrepList, PrepList[[1]]$itemnames, MG = TRUE))            
        }                
    }    
    pars <- vector('list', Data$ngroups)
    for(g in 1:Data$ngroups)
        pars[[g]] <- PrepList[[g]]$pars
    nitems <- Data$nitems
    K <- PrepList[[1]]$K
    Data$nfact <- nfact <- PrepList[[1]]$pars[[nitems+1]]@nfact      
    if(nfact != 1 && any(c('1PL') %in% itemtype )) 
        stop('1PL itemtype for multidimenional models is ambiguous. Please specify the 
             appropriate constraints manually using the 2PL model and the constrain argument.')
    if(nfact != 1 && any(c('Rasch') %in% itemtype ) && PrepList[[1]]$exploratory) 
       stop('Rasch itemtypes are for confimatory models only.')
    nLambdas <- PrepList[[1]]$pars[[1]]@nfact
    if(is.null(constrain)) constrain <- list()         
    #default MG uses configural model (independent groups but each identified)        
    if('free_means' %in% invariance ){ #Free factor means (means 0 for ref)
        for(g in 2:Data$ngroups)
            pars[[g]][[nitems + 1]]@est[1:nfact] <- TRUE             
    }    
    if('free_varcov' %in% invariance){ #Free factor vars and covs (vars 1 for ref)        
        for(g in 2:Data$ngroups)
            pars[[g]][[nitems + 1]]@est[(nfact+1):length(pars[[g]][[nitems + 1]]@est)] <- TRUE            
    }   
    constrain <- UpdateConstrain(pars=pars, constrain=constrain, invariance=invariance, nfact=Data$nfact, 
                                 nLambdas=nLambdas, J=nitems, ngroups=Data$ngroups, PrepList=PrepList, 
                                 mixedlist=mixedlist, method=opts$method, itemnames=PrepList[[1]]$itemnames)     
    startlongpars <- c()    
    if(opts$NULL.MODEL){        
        constrain <- list()
        for(i in 1:nitems){
            pars[[1]][[i]]@par[1] <- 0
            pars[[1]][[i]]@est[1] <- FALSE            
            if(is(pars[[1]][[i]], 'nominal'))
                pars[[1]][[i]]@est[(nfact+1):(nfact + K[i])] <- FALSE             
            if(is(pars[[1]][[i]], 'mcm')) 
                pars[[1]][[i]]@est[c((nfact+1):(nfact + K[i]+1), 
                    length(pars[[1]][[i]]@est):(length(pars[[1]][[i]]@est)-K[i]+1))] <- FALSE                           
            if(is(pars[[1]][[i]], 'nestlogit'))
                pars[[1]][[i]]@est[(nfact+5):(nfact + K[i] + 1)] <- FALSE             
        }       
    }
    #EM estimation
    G2group <- X2group <- numeric(Data$ngroups)        
    if(opts$method == 'EM'){        
        if(opts$method == 'EM' && nLambdas > nfact) 
            stop('Polynominals and product terms not supported for EM method')        
        if (is.null(opts$quadpts)) opts$quadpts <- ceiling(40/(nfact^1.5))
        Theta <- theta <- as.matrix(seq(-(.8 * sqrt(opts$quadpts)), .8 * sqrt(opts$quadpts),
                                        length.out = opts$quadpts))        
        temp <- matrix(0,nrow=nitems,ncol=(nfact-1))
        sitems <- matrix(0, nrow=sum(PrepList[[1]]$K), ncol=(nfact-1))
        if(opts$BFACTOR){                
            for(i in 1:nitems) temp[i, oldmodel[i]] <- 1
            ind <- 1
            for(i in 1:nitems){
                for(j in 1:PrepList[[1]]$K[i]){
                    sitems[ind, ] <- temp[i, ]
                    ind <- ind + 1
                }		
            }    
            theta <- seq(-4, 4, length.out = opts$quadpts)
            Theta <- thetaComb(theta, 2)
            Theta <- cbind(Theta[,1], matrix(Theta[,2], nrow=nrow(Theta), ncol=ncol(sitems)))            
        } else {
            if(opts$quadpts^nfact <= opts$MAXQUAD){
                Theta <- thetaComb(theta, nfact)    	
            } else stop('Greater than ', opts$MAXQUAD, ' quadrature points.')
        }
        ESTIMATE <- EM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                             list = list(NCYCLES=opts$NCYCLES, TOL=opts$TOL, MSTEPMAXIT=opts$MSTEPMAXIT,
                                         nfactNames=PrepList[[1]]$nfactNames, theta=theta,
                                         itemloc=PrepList[[1]]$itemloc, BFACTOR=opts$BFACTOR,
                                         sitems=sitems, specific=oldmodel, NULL.MODEL=opts$NULL.MODEL,
                                         nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                         SEM=opts$SE.type == 'SEM' && opts$SE), 
                             Theta=Theta)         
        startlongpars <- ESTIMATE$longpars                     
        rlist <- ESTIMATE$rlist
        logLik <- G2 <- X2 <- SElogLik <- 0          
        for(g in 1:Data$ngroups){
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
            if(opts$verbose) cat('\nCalculating information matrix...\n')
            if(opts$SE.type == 'SEM'){                                
                dontrun <- FALSE
                if(ESTIMATE$cycles <= 5 ){
                    dontrun <- TRUE
                    warning('Too few EM interations to compute SEM information matrix')               
                }
                if(!dontrun){
                    if(ESTIMATE$cycles <= 10) 
                        message('Very few EM cycles performed. Consider decreasing TOL further to 
                                increase EM iteration count')
                    estmat <- matrix(FALSE, length(ESTIMATE$correction), length(ESTIMATE$correction))
                    DM <- estmat + 0
                    diag(estmat) <- TRUE                    
                    if(!is.null(list(...)$cl)){                        
                        DM <- t(parallel::parApply(cl=list(...)$cl, estmat, MARGIN=1, FUN=SEM.SE, 
                                pars=ESTIMATE$pars, constrain=constrain, PrepList=PrepList,
                                list = list(NCYCLES=opts$NCYCLES, TOL=opts$TOL, MSTEPMAXIT=opts$MSTEPMAXIT,
                                           nfactNames=PrepList[[1]]$nfactNames, theta=theta,
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=opts$BFACTOR,
                                           sitems=sitems, specific=oldmodel, NULL.MODEL=opts$NULL.MODEL,
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose), 
                                Theta=Theta, theta=theta, ESTIMATE=ESTIMATE))
                    } else {
                        for(i in 1:ncol(DM))
                            DM[i, ] <- SEM.SE(est=estmat[i,], pars=ESTIMATE$pars, constrain=constrain,
                                              PrepList=PrepList,
                                          list = list(NCYCLES=opts$NCYCLES, TOL=opts$TOL, MSTEPMAXIT=opts$MSTEPMAXIT,
                                                      nfactNames=PrepList[[1]]$nfactNames, theta=theta,
                                                      itemloc=PrepList[[1]]$itemloc, BFACTOR=opts$BFACTOR,
                                                      sitems=sitems, specific=oldmodel, NULL.MODEL=opts$NULL.MODEL,
                                                      nfact=nfact, constrain=constrain, verbose=opts$verbose), 
                                        Theta=Theta, theta=theta, ESTIMATE=ESTIMATE)
                    }
                    info <- solve(-solve(ESTIMATE$hess) %*% solve(diag(ncol(DM)) - DM))
                    ESTIMATE <- loadESTIMATEinfo(info=info, ESTIMATE=ESTIMATE, constrain=constrain)
                }
            }
            if(opts$SE.type == 'BL')
                ESTIMATE <- BL.SE(pars=ESTIMATE$pars, Theta=Theta, theta=theta, PrepList=PrepList, 
                                  BFACTOR=opts$BFACTOR, itemloc=PrepList[[1]]$itemloc, ESTIMATE=ESTIMATE, 
                                  constrain=constrain, specific=oldmodel, sitems=sitems)
            if(opts$SE.type == 'MHRM')            
                ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=1, SEMCYCLES=5,
                                           KDRAWS=opts$KDRAWS, TOL=opts$SEtol, USEEM=opts$USEEM, 
                                           gain=opts$gain, 
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=opts$BFACTOR,  
                                           nfact=nfact, constrain=constrain, verbose=FALSE,
                                           startlongpars=startlongpars))                 
            ESTIMATE$cycles <- tmp$cycles
        }
    } else if(opts$method == 'MHRM'){ #MHRM estimation
        Theta <- matrix(0, Data$N, nitems)             
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=opts$BURNIN, 
                                           SEMCYCLES=opts$SEMCYCLES, gain=opts$gain, 
                                           KDRAWS=opts$KDRAWS, TOL=opts$TOL, USEEM=FALSE, 
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=opts$BFACTOR, 
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                           startlongpars=startlongpars))                
        rlist <- vector('list', Data$ngroups)        
        for(g in 1:Data$ngroups)
            rlist[[g]]$expected = numeric(1)
    } else if(opts$method == 'MIXED'){  
        ESTIMATE <- MHRM.mixed(pars=pars, constrain=constrain, 
                                    PrepList=PrepList, mixedlist=mixedlist,                                    
                               list = list(NCYCLES=opts$NCYCLES, BURNIN=opts$BURNIN, 
                                           SEMCYCLES=opts$SEMCYCLES, gain=opts$gain,
                                           KDRAWS=opts$KDRAWS, TOL=opts$TOL, USEEM=FALSE,  
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=opts$BFACTOR, 
                                           nfact=nfact, constrain=constrain, verbose=opts$verbose,
                                           startlongpars=startlongpars))                
        rlist <- vector('list', Data$ngroups)        
        for(g in 1:Data$ngroups)
            rlist[[g]]$expected = numeric(1)
    }   
    cmods <- list()
    for(g in 1:Data$ngroups){        
        lambdas <- Lambdas(ESTIMATE$pars[[g]]) * opts$D/1.702
        if (ncol(lambdas) > 1) norm <- sqrt(1 + rowSums(lambdas^2))
        else norm <- as.matrix(sqrt(1 + lambdas[ ,1]^2))  
        alp <- as.matrix(lambdas/norm)
        F <- alp
        if(opts$method != 'MIXED')
            colnames(F) <- PrepList[[g]]$factorNames
        h2 <- rowSums(F^2)        
        cmods[[g]] <- new('ConfirmatoryClass', pars=ESTIMATE$pars[[g]], itemloc=PrepList[[g]]$itemloc, 
                          tabdata=PrepList[[g]]$tabdata2, data=Data$data[group == Data$groupNames[[g]], ], 
                          converge=ESTIMATE$converge, esttype='MHRM', F=F, h2=h2,                
                          K=PrepList[[g]]$K, tabdatalong=PrepList[[g]]$tabdata, nfact=nfact, 
                          constrain=constrain, G2=G2group[g], X2=X2group[g], Pl = rlist[[g]]$expected,
                          mixedlist=if(opts$method == 'MIXED') mixedlist else list(),
                          fulldata=PrepList[[g]]$fulldata, factorNames=PrepList[[g]]$factorNames)        
    }
    rm(lambdas, norm)
    #missing stats for MHRM
    if(opts$method =='MHRM' || opts$method == 'MIXED'){        
        if(opts$verbose) cat("\nCalculating log-likelihood...\n")
        flush.console()      
        logLik <- G2 <- X2 <- SElogLik <- 0        
        Pl <- list()        
        for(g in 1:Data$ngroups){
            cmods[[g]] <- calcLogLik(cmods[[g]], opts$draws, G2 = 'return', cl=opts$cl)                
            logLik <- logLik + cmods[[g]]@logLik
            SElogLik <- SElogLik + cmods[[g]]@SElogLik
            G2 <- G2 + cmods[[g]]@G2            
            Pl[[g]] <- cmods[[g]]@Pl            
            X2 <- X2 + cmods[[g]]@X2
        }            
    } 
    #constraint for grsm blocks
    if(any(itemtype == 'grsm')){
        for(g in 1:Data$ngroups){
            for(u in unique(na.omit(Data$grsm.block))){
                tmp <- Data$grsm.block == u
                tmp2 <- rep(0, nitems)                
                for(i in 1:nitems){
                    if(tmp[i]){
                        ind <- length(ESTIMATE$pars[[g]][[i]]@par)
                        tmp2[i] <- ESTIMATE$pars[[g]][[i]]@par[ind]
                    }
                }
                tmp2 <- tmp2 - mean(tmp2)
                for(i in 1:nitems){
                    if(tmp[i]){
                        ind <- length(ESTIMATE$pars[[g]][[i]]@par)
                        ESTIMATE$pars[[g]][[i]]@par[ind] <- tmp2[i]
                    }               
                }
            }
        }       
    }
    
    ####post estimation stats    
    df <- rr <- 0    
    for(g in 1:Data$ngroups){
        r <- PrepList[[g]]$tabdata
        r <- r[, ncol(r)]
        rr <- rr + r
        df <- df + sum(r != 0) - 1 
    }
    r <- rr    
    N <- sum(r)
    logN <- 0
    logr <- rep(0,length(r))
    for (i in 1:N) logN <- logN + log(i)
    for (i in 1:length(r)) 
        for (j in 1:r[i]) 
            logr[i] <- logr[i] + log(j)    		
    if(sum(logr) != 0)		
        logLik <- logLik + logN/sum(logr)							
    nestpars <- nconstr <- 0
    for(g in 1:Data$ngroups)
        for(i in 1:(nitems+1))
            nestpars <- nestpars + sum(pars[[g]][[i]]@est)
    if(length(constrain) > 0)
        for(i in 1:length(constrain))
            nconstr <- nconstr + length(constrain[[i]]) - 1     
    nmissingtabdata <- sum(is.na(rowSums(PrepList[[1]]$tabdata2)))
    df <- df - nestpars + nconstr + nfact*(nfact - 1)/2 - nmissingtabdata	
    tmp <- (length(r) - df - 1)
    AIC <- (-2) * logLik + 2 * tmp
    AICc <- AIC + 2 * tmp * (tmp + 1) / (length(r) - tmp - 1)
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
    if(!opts$NULL.MODEL && opts$method != 'MIXED' && opts$calcNull){        
        null.mod <- try(unclass(mirt(data, 1, itemtype=itemtype, technical=list(NULL.MODEL=TRUE),
                                     large=opts$PrepList, key=key)))
        if(is(null.mod, 'try-error')){
            message('Null model calculation did not converge.')
            null.mod <- unclass(new('ConfirmatoryClass'))            
        } else {
            TLI.G2 <- (null.mod@G2 / null.mod@df - G2/df) / (null.mod@G2 / null.mod@df - 1)
            TLI.X2 <- (null.mod@X2 / null.mod@df - X2/df) / (null.mod@X2 / null.mod@df - 1)
            CFI.G2 <- 1 - (G2 - df) / (null.mod@G2 - null.mod@df)
            CFI.X2 <- 1 - (X2 - df) / (null.mod@X2 - null.mod@df)
        }
    }
    if(X2/G2 > 4) TLI.X2 <- CFI.X2 <- X2 <- p.X2 <- RMSEA.X2 <- NaN
    if(nmissingtabdata > 0) 
        p.G2 <- p.X2 <- RMSEA.G2 <- RMSEA.X2 <- G2 <- X2 <- TLI.G2 <- 
            TLI.X2 <- CFI.G2 <- CFI.X2 <- NaN
    if(is.null(parprior)) parprior <- list()
    if(Data$ngroups == 1){        
        if(opts$method == 'MIXED'){                        
            mixedlist$betas <- cmods[[1]]@pars[[1]]@par[1:ncol(mixedlist$FDL[[1]])]
            mixedlist$SEbetas <- cmods[[1]]@pars[[1]]@SEpar[1:ncol(mixedlist$FDL[[1]])]
            names(mixedlist$SEbetas) <- names(mixedlist$betas) <- 
                colnames(mixedlist$FDL[[1]])
            mod <- new('MixedClass', 
                       iter=ESTIMATE$cycles, 
                       pars=cmods[[1]]@pars,                         
                       model=list(oldmodel),                       
                       df=df,                        
                       itemloc=PrepList[[1]]$itemloc, 
                       method=opts$method,
                       AIC=AIC, 
                       AICc=AICc,
                       BIC=BIC, 
                       SABIC=SABIC,
                       logLik=logLik, 
                       F=F, 
                       h2=h2, 
                       tabdata=PrepList[[1]]$tabdata2, 
                       Pl=Pl[[1]], 
                       data=Data$data, 
                       converge=ESTIMATE$converge, 
                       nfact=nfact,                                       
                       K=PrepList[[1]]$K, 
                       tabdatalong=PrepList[[1]]$tabdata,                        
                       mixedlist=mixedlist,                       
                       factorNames=PrepList[[1]]$factorNames, 
                       constrain=PrepList[[1]]$constrain, 
                       parprior=parprior,
                       fulldata=PrepList[[1]]$fulldata,
                       itemtype=PrepList[[1]]$itemtype,
                       information=ESTIMATE$info)
        } else if(PrepList[[1]]$exploratory){
            FF <- alp %*% t(alp)
            V <- eigen(FF)$vector[ ,1:nfact]
            L <- eigen(FF)$values[1:nfact]
            if (nfact == 1) F <- as.matrix(V * sqrt(L))
            else F <- V %*% sqrt(diag(L))  
            if (sum(F[ ,1] < 0)) F <- (-1) * F 
            colnames(F) <- paste("F_", 1:ncol(F),sep="")
            h2 <- rowSums(F^2)
            mod <- new('ExploratoryClass', iter=ESTIMATE$cycles, 
                       pars=cmods[[1]]@pars, 
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
                       itemloc=PrepList[[1]]$itemloc, 
                       method=opts$method,
                       AIC=AIC, 
                       AICc=AICc,
                       BIC=BIC,
                       SABIC=SABIC,
                       logLik=logLik, 
                       F=F, 
                       h2=h2, 
                       tabdata=PrepList[[1]]$tabdata2, 
                       Theta=Theta, 
                       Pl=Pl[[1]], 
                       data=Data$data, 
                       converge=ESTIMATE$converge, 
                       nfact=nfact,               
                       quadpts=opts$quadpts,                        
                       K=PrepList[[1]]$K, 
                       tabdatalong=PrepList[[1]]$tabdata, 
                       rotate=opts$rotate,
                       null.mod=null.mod,                        
                       Target=opts$Target,
                       factorNames=PrepList[[1]]$factorNames, 
                       constrain=PrepList[[1]]$constrain, 
                       parprior=parprior,
                       fulldata=PrepList[[1]]$fulldata,
                       itemtype=PrepList[[1]]$itemtype,
                       information=ESTIMATE$info)
        } else {                        
            mod <- new('ConfirmatoryClass', iter=ESTIMATE$cycles, 
                       pars=cmods[[1]]@pars, 
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
                       itemloc=PrepList[[1]]$itemloc, 
                       AIC=AIC, 
                       AICc=AICc,
                       BIC=BIC, 
                       SABIC=SABIC,
                       logLik=logLik, 
                       F=F, 
                       h2=h2, 
                       tabdata=PrepList[[1]]$tabdata2, 
                       Theta=Theta, 
                       method=opts$method,
                       Pl=Pl[[1]], 
                       data=Data$data, 
                       converge=ESTIMATE$converge, 
                       nfact=nfact,               
                       quadpts=opts$quadpts,                        
                       K=PrepList[[1]]$K, 
                       tabdatalong=PrepList[[1]]$tabdata, 
                       null.mod=null.mod,                       
                       factorNames=PrepList[[1]]$factorNames, 
                       constrain=PrepList[[1]]$constrain, 
                       parprior=parprior,
                       fulldata=PrepList[[1]]$fulldata,
                       itemtype=PrepList[[1]]$itemtype,
                       information=ESTIMATE$info)
        }        
    } else {
        tabdatalong <- PrepList[[1]]$tabdata
        tabdata <- PrepList[[1]]$tabdata2
        tabdata[,ncol(tabdata)] <- tabdatalong[,ncol(tabdatalong)] <- r
        mod <- new('MultipleGroupClass', iter=ESTIMATE$cycles, 
                   cmods=cmods, 
                   model=list(oldmodel),
                   itemloc=PrepList[[1]]$itemloc, 
                   tabdata=tabdata, 
                   data=Data$data, 
                   converge=ESTIMATE$converge, 
                   esttype=opts$method,                
                   K=PrepList[[1]]$K, 
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
                   itemtype=PrepList[[1]]$itemtype,
                   information=ESTIMATE$info)  
    }
    mod@time <- proc.time()[3] - opts$start.time  
    return(mod)
}
