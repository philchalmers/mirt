ESTIMATION <- function(data, model, group, itemtype = NULL, guess = 0, upper = 1, 
                       invariance = '', pars = NULL, method = 'MHRM', constrain = NULL, 
                       parprior = NULL, draws = 2000, calcLL = TRUE,
                       quadpts = NaN, rotate = 'varimax', Target = NaN, SE = TRUE,
                       technical = list(), debug = FALSE, verbose = TRUE, BFACTOR = FALSE,
                       SEtol = .01, nested.mod = NULL, grsm.block = NULL, D = 1.702, 
                       multilevel=NULL, ...)
{    
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
    PrepListFull <- PrepData(data=data, model=model[[1]], itemtype=itemtype, guess=guess, upper=upper, 
                             startvalues=NULL, constrain=NULL, freepars=NULL, 
                             parprior=NULL, verbose=verbose, debug=debug, free.start=NULL, D=D,
                             technical=technical) #just a dummy model to collect fulldata stuff
    for(g in 1:ngroups){    
        select <- group == groupNames[g]        
        tmp <- 1:ngroups
        selectmod <- model[[tmp[names(model) == groupNames[g]]]]
        PrepList[[g]] <- PrepData(data=data[select,], model=selectmod, itemtype=itemtype, guess=guess, 
                                  upper=upper, startvalues=NULL, constrain=constrain, freepars=NULL, 
                                  parprior=parprior, verbose=verbose, debug=debug, free.start=NULL,
                                  technical=technical, parnumber=parnumber, BFACTOR=BFACTOR,
                                  grsm.block=grsm.block, D=D)        
        tmp <- PrepList[[g]]$pars[[length(PrepList[[g]]$pars)]]
        parnumber <- tmp@parnum[length(tmp@parnum)] + 1
    }    
    names(PrepList) <- groupNames    
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
    if(nfact != 1 && any(c('1PL', 'Rasch') %in% itemtype )) 
       stop('1PL and Rasch itemtypes are for unidimensional models only.')
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
                                 nLambdas=nLambdas, J=J, ngroups=ngroups, PrepList=PrepList)    
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
        Theta <- theta <- as.matrix(seq(-4,4,length.out = quadpts))
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
                                         sitems=sitems, specific=oldmodel,
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
            Ng <- sum(rg) 
            G2group[g] <- 2 * sum(rg * log(rg/(Ng*Pl)))
            G2 <- G2 + G2group[g]
            logLik <- logLik + sum(rg*log(Pl))
        }
        Pl <- list(Pl)
        if(!NULL.MODEL && SE){
            ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=NCYCLES, BURNIN=1, SEMCYCLES=5,
                                           KDRAWS=KDRAWS, TOL=SEtol, USEEM=USEEM, gain=gain, 
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=BFACTOR,  
                                           nfact=nfact, constrain=constrain, verbose=FALSE,
                                           startlongpars=startlongpars), 
                               debug=debug)                 
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
    } else if(method == 'MULTILEVEL'){        
        if(is.null(pars)) {
            burninmod <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                                   list = list(NCYCLES=2, BURNIN=BURNIN, SEMCYCLES=2,
                                               KDRAWS=KDRAWS, TOL=TOL, USEEM=FALSE, gain=gain, 
                                               nfactNames=PrepList[[1]]$nfactNames, 
                                               itemloc=PrepList[[1]]$itemloc, BFACTOR=BFACTOR, 
                                               nfact=nfact, constrain=constrain, verbose=verbose,
                                               startlongpars=startlongpars), 
                                   debug=debug)        
            pars <- burninmod@pars        
        }
        ESTIMATE <- MHRM.multilevel(pars=pars, constrain=constrain, 
                                    PrepList=PrepList, multilevel=multilevel,                                    
                               list = list(NCYCLES=NCYCLES, BURNIN=BURNIN, SEMCYCLES=SEMCYCLES,
                                           KDRAWS=KDRAWS, TOL=TOL, USEEM=FALSE, gain=gain, 
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc, BFACTOR=BFACTOR, 
                                           nfact=nfact, constrain=constrain, verbose=verbose,
                                           startlongpars=startlongpars), 
                               debug=debug, ...)        
        iter <- ESTIMATE$cycles
    }   
    cmods <- list()
    for(g in 1:ngroups){
        lambdas <- Lambdas(ESTIMATE$pars[[g]]) * 1.702/D
        if (ncol(lambdas) > 1) norm <- sqrt(1 + rowSums(lambdas^2))
        else norm <- as.matrix(sqrt(1 + lambdas[ ,1]^2))  
        alp <- as.matrix(lambdas/norm)
        F <- alp
        colnames(F) <- PrepList[[g]]$factorNames
        h2 <- rowSums(F^2)
        cmods[[g]] <- new('ConfirmatoryClass', pars=ESTIMATE$pars[[g]], itemloc=PrepList[[g]]$itemloc, 
                          tabdata=PrepList[[g]]$tabdata2, data=data[group == groupNames[[g]], ], 
                          converge=ESTIMATE$converge, esttype='MHRM', F=F, h2=h2,                
                          K=PrepList[[g]]$K, tabdatalong=PrepList[[g]]$tabdata, nfact=nfact, 
                          constrain=constrain, G2=G2group[g], 
                          fulldata=PrepList[[g]]$fulldata, factorNames=PrepList[[g]]$factorNames)        
    }
    #missing stats for MHRM
    if(method =='MHRM'){
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
        df <- df + length(r) - 1 
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
    AIC <- (-2) * logLik + 2 * (length(r) - df - 1)
    BIC <- (-2) * logLik + (length(r) - df - 1)*log(N) 
    p <- 1 - pchisq(G2,df)
    RMSEA <- ifelse((G2 - df) > 0, 
                    sqrt(G2 - df) / sqrt(df * (N-1)), 0)
    null.mod <- unclass(new('ConfirmatoryClass'))
    TLI <- NaN
    if(!NULL.MODEL){
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
        if(method == 'MULTILEVEL'){
            mod <- new('MultilevelClass', 
                       iter=ESTIMATE$cycles, 
                       pars=cmods[[1]]@pars, 
                       G2=G2, 
                       df=df, 
                       p=p, 
                       itemloc=PrepListFull$itemloc, 
                       method=method,
                       AIC=AIC, 
                       BIC=BIC, 
                       logLik=logLik, 
                       F=F, 
                       h2=h2, 
                       tabdata=PrepListFull$tabdata2, 
                       Pl=Pl[[1]], 
                       data=data, 
                       converge=ESTIMATE$converge, 
                       nfact=nfact,               
                       RMSEA=RMSEA, 
                       K=PrepListFull$K, 
                       tabdatalong=PrepListFull$tabdata, 
                       null.mod=null.mod, #? what to do about this....
                       TLI=TLI, 
                       factorNames=PrepListFull$factorNames, 
                       constrain=PrepList[[1]]$constrain, 
                       fulldata=PrepListFull$fulldata,
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
                       BIC=BIC, 
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
                       information=ESTIMATE$info)
        } else {                        
            mod <- new('ConfirmatoryClass', iter=ESTIMATE$cycles, 
                       pars=cmods[[1]]@pars, 
                       G2=G2, 
                       df=df, 
                       p=p, 
                       itemloc=PrepListFull$itemloc, 
                       AIC=AIC, 
                       BIC=BIC, 
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
                   BIC=BIC,
                   nfact=nfact,
                   G2=G2,
                   RMSEA=RMSEA,
                   TLI=TLI,
                   p=p,
                   Theta=Theta,
                   Pl=Pl,
                   information=ESTIMATE$info)  
    }
    return(mod)
}
