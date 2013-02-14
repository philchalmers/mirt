ESTIMATION <- function(data, model, group, itemtype = NULL, guess = 0, upper = 1, 
                       invariance = '', pars = NULL, method = 'MHRM', constrain = NULL, 
                       parprior = NULL, draws = 2000, calcLL = TRUE,
                       quadpts = NaN, rotate = 'varimax', Target = NaN, SE = TRUE,
                       technical = list(), debug = FALSE, verbose = TRUE, BFACTOR = FALSE,
                       SEtol = .01, nested.mod = NULL, grsm.block = NULL, D = 1.702, 
                       rsm.block = NULL, mixedlist=NULL, calcNull=TRUE, ...)
{    
    start.time <- proc.time()[3]
    if(debug == 'ESTIMATION') browser()    
    set.seed(12345)       
    MAXQUAD <- ifelse(is.null(technical$MAXQUAD), 10000, technical$MAXQUAD)
    MSTEPMAXIT <- ifelse(is.null(technical$MSTEPMAXIT), 15, technical$MSTEPMAXIT)        
    NCYCLES <- ifelse(is.null(technical$NCYCLES), 2000, technical$NCYCLES)
    if(method == 'EM')
        NCYCLES <- ifelse(is.null(technical$NCYCLES), 300, technical$NCYCLES)
    BURNIN <- ifelse(is.null(technical$BURNIN), 150, technical$BURNIN)
    SEMCYCLES <- ifelse(is.null(technical$SEMCYCLES), 50, technical$SEMCYCLES)
    KDRAWS  <- ifelse(is.null(technical$KDRAWS), 1, technical$KDRAWS)
    TOL <- ifelse(is.null(technical$TOL), .001, technical$TOL)      
    if(!is.null(technical$set.seed)) set.seed(technical$set.seed)	
    gain <- c(0.05,0.5,0.004)
    if(!is.null(technical$gain)){
        if(length(technical$gain) == 3 && is.numeric(technical$gain))
            gain <- technical$gain
    }	 
    NULL.MODEL <- ifelse(is.null(technical$NULL.MODEL), FALSE, TRUE)
    USEEM <- ifelse(method == 'EM', TRUE, FALSE)
    #change itemtypes if NULL.MODEL 
    if(NULL.MODEL){
        constrain <- NULL        
        if(!is.null(itemtype)){
            itemtype[itemtype == 'grsm'] <- 'graded'
            itemtype[itemtype == 'rsm'] <- 'gpcm'        
        }
    }
    ##	            
    data <- as.matrix(data)
    if(is.null(grsm.block)) grsm.block <- rep(1, ncol(data))
    rownames(data) <- 1:nrow(data)
    group <- factor(group)
    groupNames <- unique(group)
    ngroups <- length(groupNames)
    oldmodel <- model
    if(is(model, 'numeric') && length(model) > 1)
        model <- bfactor2mod(model, ncol(data))
    if(length(model) == 1){
        newmodel <- list()
        for(g in 1:ngroups)
            newmodel[[g]] <- model
        names(newmodel) <- groupNames
        model <- newmodel
    }     
    parnumber <- 1
    PrepList <- vector('list', ngroups)        
    names(PrepList) <- groupNames    
    tmp <- 1:ngroups
    selectmod <- model[[tmp[names(model) == groupNames[1]]]]
    PrepListFull <- PrepList[[1]] <- 
        PrepData(data=data, model=selectmod, itemtype=itemtype, guess=guess, 
                             upper=upper, startvalues=NULL, constrain=NULL, freepars=NULL, 
                             parprior=parprior, verbose=verbose, debug=debug, free.start=NULL,
                             technical=technical, parnumber=parnumber, BFACTOR=BFACTOR,
                             grsm.block=grsm.block, D=D, mixedlist=mixedlist, ...)            
    stringtabdata <- apply(PrepListFull$tabdata2[, -ncol(PrepListFull$tabdata2)], 
                           1, paste, sep='', collapse = '/')
    stringfulldata <- apply(data, 1, paste, sep='', collapse = '/')     
    for(g in 1:ngroups){                    
        tmp <- 1:ngroups
        selectmod <- model[[tmp[names(model) == groupNames[g]]]]
        if(g != 1)
            PrepList[[g]] <- PrepData(data=data, model=selectmod, itemtype=itemtype, guess=guess, 
                                  upper=upper, startvalues=NULL, constrain=NULL, freepars=NULL, 
                                  parprior=parprior, verbose=verbose, debug=debug, free.start=NULL,
                                  technical=technical, parnumber=parnumber, BFACTOR=BFACTOR,
                                  grsm.block=grsm.block, D=D, mixedlist=mixedlist, ...)        
        tmp <- PrepList[[g]]$pars[[length(PrepList[[g]]$pars)]]
        parnumber <- tmp@parnum[length(tmp@parnum)] + 1
    }    
    if(ngroups > 1) {
        tmprs <- makerData(stringfulldata=stringfulldata, stringtabdata=stringtabdata,                                
                               r=PrepListFull$tabdata[,ncol(PrepListFull$tabdata)], 
                               group=group, groupNames=groupNames)
        for(g in 1:ngroups){              
            select <- group == groupNames[g]
            for(i in 1:ncol(data))
                PrepList[[g]]$pars[[i]]@dat <- PrepList[[g]]$pars[[i]]@dat[select, , drop = FALSE]
            PrepList[[g]]$fulldata <- PrepListFull$fulldata[select, ]
            PrepList[[g]]$tabdata[,ncol(PrepListFull$tabdata)] <- tmprs[,g]
            PrepList[[g]]$tabdata2[,ncol(PrepListFull$tabdata2)] <- tmprs[,g]
        }        
    }    
    if(BFACTOR){
        #better start values        
        J <- length(PrepList[[1]]$pars) - 1
        Rpoly <- cormod(na.omit(data), PrepListFull$K, guess = rep(0,J))
        suppressWarnings(FA <- psych::fa(Rpoly, 1, warnings = FALSE))
        loads <- unclass(FA$load)
        cs <- sqrt(abs(FA$u))              
        astart <- loads/cs
        astart <- cbind(astart,astart/2)* (1.702 / D) #reweight due to D
        nfact <- PrepList[[1]]$pars[[1]]@nfact
        for(g in 1:ngroups)
            for(i in 1:J)
                PrepList[[g]]$pars[[i]]@par[PrepList[[g]]$pars[[i]]@est][1:2] <- astart[i, ]                
    }
    if(!is.null(nested.mod) && is(nested.mod, 'MultipleGroupClass')){          
        for(g in 1:ngroups){
            for(i in 1:length(PrepList[[g]]$pars)){
                tmp <- nested.mod@cmods[[g]]@pars
                tmp2 <- PrepList[[g]]$pars[[i]]@est
                PrepList[[g]]$pars[[i]]@par[tmp2] <- tmp[[i]]@par[tmp2]                
            }            
        }       
    }
    if(!is.null(pars)){
        if(is(pars, 'matrix') || is(pars, 'data.frame')){
            PrepList <- UpdatePrepList(PrepList, pars, MG = TRUE)
        } else if(pars == 'values'){            
            return(ReturnPars(PrepList, PrepList[[1]]$itemnames, MG = TRUE))            
        }                
    }    
    pars <- vector('list', ngroups)
    for(g in 1:ngroups)
        pars[[g]] <- PrepList[[g]]$pars
    J <- length(PrepList[[1]]$itemtype)
    K <- PrepListFull$K
    nfact <- PrepList[[1]]$pars[[J+1]]@nfact    
    if(nfact != 1 && any(c('1PL') %in% itemtype )) 
        stop('1PL itemtype for multidimenional models is ambiguous. Please specify the 
             appropriate constraints manually using the 2PL model and the constrain argument.')
    if(nfact != 1 && any(c('Rasch') %in% itemtype ) && PrepList[[1]]$exploratory) 
       stop('Rasch itemtypes are for confimatory models only.')
    nLambdas <- PrepList[[1]]$pars[[1]]@nfact
    if(is.null(constrain)) constrain <- list()         
    #default MG uses configural model (independent groups but each identified)        
    if('free_means' %in% invariance ){ #Free factor means (means 0 for ref)
        for(g in 2:ngroups)
            pars[[g]][[J + 1]]@est[1:nfact] <- TRUE             
    }    
    if('free_varcov' %in% invariance){ #Free factor vars and covs (vars 1 for ref)        
        for(g in 2:ngroups)
            pars[[g]][[J + 1]]@est[(nfact+1):length(pars[[g]][[J + 1]]@est)] <- TRUE            
    }   
    constrain <- UpdateConstrain(pars=pars, constrain=constrain, invariance=invariance, nfact=nfact, 
                                 nLambdas=nLambdas, J=J, ngroups=ngroups, PrepList=PrepList, 
                                 mixedlist=mixedlist, method=method)        
    if(!is.null(technical$return_newconstrain)) return(constrain)    
    startlongpars <- c()
    if(NULL.MODEL){
        constrain <- list()
        for(i in 1:J){
            pars[[1]][[i]]@par[1] <- 0
            pars[[1]][[i]]@est[1] <- FALSE            
            if(is(pars[[1]][[i]], 'nominal'))
                pars[[1]][[i]]@est[(nfact+1):(nfact + K[i])] <- FALSE             
            if(is(pars[[1]][[i]], 'mcm')) 
                pars[[1]][[i]]@est[c((nfact+1):(nfact + K[i]+1), 
                    length(pars[[1]][[i]]@est):(length(pars[[1]][[i]]@est)-K[i]+1))] <- FALSE                           
        }       
    }
    #EM estimation
    G2group <- numeric(ngroups)
    if(method == 'EM'){
        esttype <- 'EM'
        if(method == 'EM' && nLambdas > nfact) 
            stop('Polynominals and product terms not supported for EM method')        
        if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))
        Theta <- theta <- as.matrix(seq(-(.8 * sqrt(quadpts)), .8 * sqrt(quadpts),
                                        length.out = quadpts))        
        temp <- matrix(0,nrow=J,ncol=(nfact-1))
        sitems <- matrix(0, nrow=sum(PrepListFull$K), ncol=(nfact-1))
        if(BFACTOR){                
            for(i in 1:J) temp[i,oldmodel[i]] <- 1
            ind <- 1
            for(i in 1:J){
                for(j in 1:PrepListFull$K[i]){
                    sitems[ind, ] <- temp[i, ]
                    ind <- ind + 1
                }		
            }    
            theta <- seq(-4,4,length.out = quadpts)
            Theta <- thetaComb(theta, 2)
            Theta <- cbind(Theta[,1], matrix(Theta[,2], nrow=nrow(Theta), ncol=ncol(sitems)))            
        } else {
            if(quadpts^nfact <= MAXQUAD){
                Theta <- thetaComb(theta,nfact)    	
            } else stop('Greater than ', MAXQUAD, ' quadrature points.')
        }
        ESTIMATE <- EM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                             list = list(NCYCLES=NCYCLES, TOL=TOL, MSTEPMAXIT=MSTEPMAXIT,
                                         nfactNames=PrepList[[1]]$nfactNames, theta=theta,
                                         itemloc=PrepList[[1]]$itemloc, BFACTOR=BFACTOR,
                                         sitems=sitems, specific=oldmodel, NULL.MODEL=NULL.MODEL,
                                         nfact=nfact, constrain=constrain, verbose=verbose), 
                             Theta=Theta, debug=debug)                
        startlongpars <- ESTIMATE$longpars
        logLik <- ESTIMATE$logLik
        SElogLik <- 0     
        rlist <- ESTIMATE$rlist
        logLik <- G2 <- 0
        iter <- ESTIMATE$cycles
        for(g in 1:ngroups){
            Pl <- rlist[[g]]$expected
            rg <- PrepList[[g]]$tabdata[,ncol(PrepList[[g]]$tabdata)]
            Pl <- Pl[rg != 0]
            rg <- rg[rg != 0]            
            Ng <- sum(rg) 
            G2group[g] <- 2 * sum(rg * log(rg/(Ng*Pl)))
            G2 <- G2 + G2group[g]
            logLik <- logLik + sum(rg*log(Pl))
        }
        Pl <- list(Pl)
        if(!NULL.MODEL && SE){
            tmp <- ESTIMATE
            ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=NCYCLES, BURNIN=1, SEMCYCLES=5,
                                           KDRAWS=KDRAWS, TOL=SEtol, USEEM=USEEM, gain=gain, 
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=BFACTOR,  
                                           nfact=nfact, constrain=constrain, verbose=FALSE,
                                           startlongpars=startlongpars), 
                               debug=debug)                 
            ESTIMATE$cycles <- tmp$cycles
        }
    } else if(method == 'MHRM'){ #MHRM estimation
        Theta <- matrix(0, nrow(data), J)
        esttype <- 'MHRM'
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=NCYCLES, BURNIN=BURNIN, SEMCYCLES=SEMCYCLES,
                                           KDRAWS=KDRAWS, TOL=TOL, USEEM=FALSE, gain=gain, 
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=BFACTOR, 
                                           nfact=nfact, constrain=constrain, verbose=verbose,
                                           startlongpars=startlongpars), 
                               debug=debug)        
        iter <- ESTIMATE$cycles
        rlist <- vector('list', ngroups)        
        for(g in 1:ngroups)
            rlist[[g]]$expected = numeric(1)
    } else if(method == 'MIXED'){  
        ESTIMATE <- MHRM.mixed(pars=pars, constrain=constrain, 
                                    PrepList=PrepList, mixedlist=mixedlist,                                    
                               list = list(NCYCLES=NCYCLES, BURNIN=BURNIN, SEMCYCLES=SEMCYCLES,
                                           KDRAWS=KDRAWS, TOL=TOL, USEEM=FALSE, gain=gain, 
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=BFACTOR, 
                                           nfact=nfact, constrain=constrain, verbose=verbose,
                                           startlongpars=startlongpars), 
                               debug=debug, ...)        
        iter <- ESTIMATE$cycles
        rlist <- vector('list', ngroups)        
        for(g in 1:ngroups)
            rlist[[g]]$expected = numeric(1)
    }   
    cmods <- list()
    for(g in 1:ngroups){
        lambdas <- Lambdas(ESTIMATE$pars[[g]]) * D/1.702
        if (ncol(lambdas) > 1) norm <- sqrt(1 + rowSums(lambdas^2))
        else norm <- as.matrix(sqrt(1 + lambdas[ ,1]^2))  
        alp <- as.matrix(lambdas/norm)
        F <- alp
        if(method != 'MIXED')
            colnames(F) <- PrepList[[g]]$factorNames
        h2 <- rowSums(F^2)        
        cmods[[g]] <- new('ConfirmatoryClass', pars=ESTIMATE$pars[[g]], itemloc=PrepList[[g]]$itemloc, 
                          tabdata=PrepList[[g]]$tabdata2, data=data[group == groupNames[[g]], ], 
                          converge=ESTIMATE$converge, esttype='MHRM', F=F, h2=h2,                
                          K=PrepList[[g]]$K, tabdatalong=PrepList[[g]]$tabdata, nfact=nfact, 
                          constrain=constrain, G2=G2group[g], Pl = rlist[[g]]$expected,
                          mixedlist=if(method == 'MIXED') mixedlist else list(),
                          fulldata=PrepList[[g]]$fulldata, factorNames=PrepList[[g]]$factorNames)        
    }
    #missing stats for MHRM
    if(method =='MHRM' || method == 'MIXED'){
        if(verbose) cat("\nCalculating log-likelihood...\n")
        flush.console()      
        logLik <- G2 <- SElogLik <- 0        
        Pl <- list()
        for(g in 1:ngroups){
            cmods[[g]] <- calcLogLik(cmods[[g]], draws, G2 = 'return')                
            logLik <- logLik + cmods[[g]]@logLik
            SElogLik <- SElogLik + cmods[[g]]@SElogLik
            G2 <- G2 + cmods[[g]]@G2
            Pl[[g]] <- cmods[[g]]@Pl
        }            
    } 
    #constraint for grsm blocks
    if(any(itemtype == 'grsm')){
        for(g in 1:ngroups){
            for(u in unique(na.omit(grsm.block))){
                tmp <- grsm.block == u
                tmp2 <- rep(0, J)                
                for(i in 1:J){
                    if(tmp[i]){
                        ind <- length(ESTIMATE$pars[[g]][[i]]@par)
                        tmp2[i] <- ESTIMATE$pars[[g]][[i]]@par[ind]
                    }
                }
                tmp2 <- tmp2 - mean(tmp2)
                for(i in 1:J){
                    if(tmp[i]){
                        ind <- length(ESTIMATE$pars[[g]][[i]]@par)
                        ESTIMATE$pars[[g]][[i]]@par[ind] <- tmp2[i]
                    }               
                }
            }
        }       
    }
    
    ####post estimation stats
    df <- 0    
    for(g in 1:ngroups){
        r <- PrepList[[g]]$tabdata
        r <- r[, ncol(r)]
        df <- df + sum(r != 0) - 1 
    }
    r <- PrepListFull$tabdata
    r <- r[, ncol(r)]
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
    for(g in 1:ngroups)
        for(i in 1:(J+1))
            nestpars <- nestpars + sum(pars[[g]][[i]]@est)
    if(length(constrain) > 0)
        for(i in 1:length(constrain))
            nconstr <- nconstr + length(constrain[[i]]) - 1     
    nmissingtabdata <- sum(is.na(rowSums(PrepListFull$tabdata2)))
    df <- df - nestpars + nconstr + nfact*(nfact - 1)/2 - nmissingtabdata	
    tmp <- (length(r) - df - 1)
    AIC <- (-2) * logLik + 2 * tmp
    AICc <- AIC + 2 * tmp * (tmp + 1) / (length(r) - tmp - 1)
    BIC <- (-2) * logLik + tmp*log(N) 
    SABIC <- (-2) * logLik + tmp*log((N+2)/24)
    p <- 1 - pchisq(G2,df)    
    RMSEA <- ifelse((G2 - df) > 0, 
                    sqrt(G2 - df) / sqrt(df * (N-1)), 0)
    null.mod <- unclass(new('ConfirmatoryClass'))
    TLI <- NaN    
    if(!NULL.MODEL && method != 'MIXED' && calcNull){
        null.mod <- try(unclass(mirt(data, 1, itemtype=itemtype, technical=list(NULL.MODEL=TRUE))))
        if(is(null.mod, 'try-error')){
            message('Null model calculation did not converge.')
            null.mod <- unclass(new('ConfirmatoryClass'))            
        } else {
            TLI <- (null.mod@G2 / null.mod@df - G2/df) / (null.mod@G2 / null.mod@df - 1)
        }
    }
    if(nmissingtabdata > 0) p <- RMSEA <- G2 <- TLI <- NaN
    if(ngroups == 1){        
        if(method == 'MIXED'){                        
            mixedlist$betas <- cmods[[1]]@pars[[1]]@par[1:ncol(mixedlist$FDL[[1]])]
            mixedlist$SEbetas <- cmods[[1]]@pars[[1]]@SEpar[1:ncol(mixedlist$FDL[[1]])]
            names(mixedlist$SEbetas) <- names(mixedlist$betas) <- 
                colnames(mixedlist$FDL[[1]])
            mod <- new('MixedClass', 
                       iter=ESTIMATE$cycles, 
                       pars=cmods[[1]]@pars,                         
                       df=df,                        
                       itemloc=PrepListFull$itemloc, 
                       method=method,
                       AIC=AIC, 
                       AICc=AICc,
                       BIC=BIC, 
                       SABIC=SABIC,
                       logLik=logLik, 
                       F=F, 
                       h2=h2, 
                       tabdata=PrepListFull$tabdata2, 
                       Pl=Pl[[1]], 
                       data=data, 
                       converge=ESTIMATE$converge, 
                       nfact=nfact,                                       
                       K=PrepListFull$K, 
                       tabdatalong=PrepListFull$tabdata,                        
                       mixedlist=mixedlist,                       
                       factorNames=PrepListFull$factorNames, 
                       constrain=PrepList[[1]]$constrain, 
                       fulldata=PrepListFull$fulldata,
                       itemtype=PrepListFull$itemtype,
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
                       G2=G2, 
                       df=df, 
                       p=p, 
                       itemloc=PrepListFull$itemloc, 
                       method=method,
                       AIC=AIC, 
                       AICc=AICc,
                       BIC=BIC,
                       SABIC=SABIC,
                       logLik=logLik, 
                       F=F, 
                       h2=h2, 
                       tabdata=PrepListFull$tabdata2, 
                       Theta=Theta, 
                       Pl=Pl[[1]], 
                       data=data, 
                       converge=ESTIMATE$converge, 
                       nfact=nfact,               
                       quadpts=quadpts, 
                       RMSEA=RMSEA, 
                       K=PrepListFull$K, 
                       tabdatalong=PrepListFull$tabdata, 
                       rotate=rotate,  #missing
                       null.mod=null.mod, 
                       TLI=TLI, 
                       Target=Target, #missing
                       factorNames=PrepListFull$factorNames, 
                       constrain=PrepList[[1]]$constrain, 
                       fulldata=PrepListFull$fulldata,
                       itemtype=PrepListFull$itemtype,
                       information=ESTIMATE$info)
        } else {                        
            mod <- new('ConfirmatoryClass', iter=ESTIMATE$cycles, 
                       pars=cmods[[1]]@pars, 
                       G2=G2, 
                       df=df, 
                       p=p, 
                       itemloc=PrepListFull$itemloc, 
                       AIC=AIC, 
                       AICc=AICc,
                       BIC=BIC, 
                       SABIC=SABIC,
                       logLik=logLik, 
                       F=F, 
                       h2=h2, 
                       tabdata=PrepListFull$tabdata2, 
                       Theta=Theta, 
                       method=method,
                       Pl=Pl[[1]], 
                       data=data, 
                       converge=ESTIMATE$converge, 
                       nfact=nfact,               
                       quadpts=quadpts, 
                       RMSEA=RMSEA, 
                       K=PrepListFull$K, 
                       tabdatalong=PrepListFull$tabdata, 
                       null.mod=null.mod, 
                       TLI=TLI, 
                       factorNames=PrepListFull$factorNames, 
                       constrain=PrepList[[1]]$constrain, 
                       fulldata=PrepListFull$fulldata,
                       itemtype=PrepListFull$itemtype,
                       information=ESTIMATE$info)
        }        
    } else {
        mod <- new('MultipleGroupClass', iter=ESTIMATE$cycles, 
                   cmods=cmods, 
                   itemloc=PrepListFull$itemloc, 
                   tabdata=PrepListFull$tabdata2, 
                   data=data, 
                   converge=ESTIMATE$converge, 
                   esttype=esttype,                
                   K=PrepListFull$K, 
                   tabdatalong=PrepListFull$tabdata, 
                   constrain=constrain,               
                   group=group, 
                   groupNames=groupNames, 
                   invariance=invariance, 
                   df=as.integer(df),
                   logLik=logLik, 
                   method=method,
                   SElogLik=SElogLik, 
                   AIC=AIC, 
                   AICc=AICc,
                   BIC=BIC,
                   SABIC=SABIC,
                   nfact=nfact,
                   G2=G2,
                   RMSEA=RMSEA,
                   TLI=TLI,
                   p=p,
                   Theta=Theta,
                   Pl=Pl,
                   itemtype=PrepListFull$itemtype,
                   information=ESTIMATE$info)  
    }
    mod@time <- proc.time()[3] - start.time  
    return(mod)
}
