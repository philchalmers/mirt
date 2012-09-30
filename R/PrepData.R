PrepData <- function(data, model, itemtype, guess, upper, startvalues, constrain, freepars, 
                     free.start, parprior, verbose, debug, technical, parnumber = 1, BFACTOR = FALSE,
                     rsm.group = NULL)
{
    if(debug == 'PrepData') browser()
    if(is.null(rsm.group)) rsm.group <- rep(1, ncol(data))
    itemnames <- colnames(data)
    keywords <- c('COV')
    data <- as.matrix(data)    	
    colnames(data) <- itemnames	
    J <- ncol(data)
    N <- nrow(data)
    exploratory <- FALSE
    if(is(model, 'numeric') && length(model) == 1){
        if(model != 1) exploratory <- TRUE        
        tmp <- tempfile('tempfile')
        cat(paste('F',1:model,' = 1-', J, "\n", sep=''), file=tmp)
        model <- confmirt.model(tmp, quiet = TRUE)        
        unlink(tmp)
    }
    if(exploratory && any(itemtype == c('PC2PL', 'PC3PL'))) 
        stop('Partially compensatory models can only be estimated within a confirmatory model')
    if(is(model, 'numeric') && length(model) > 1)
        model <- bfactor2mod(model, J)
    if(length(guess) == 1) guess <- rep(guess,J)
    if(length(guess) > J || length(guess) < J) 
        stop("The number of guessing parameters is incorrect.")					
    if(length(upper) == 1) upper <- rep(upper,J)
    if(length(upper) > J || length(upper) < J) 
        stop("The number of upper bound parameters is incorrect.")
    uniques <- list()
    for(i in 1:J)
        uniques[[i]] <- sort(unique(data[,i]))
    K <- rep(0,J)
    for(i in 1:J) K[i] <- length(uniques[[i]])	
    guess[K > 2] <- 0
    upper[K > 2] <- 1		
    if(is.null(itemtype)) {
        itemtype <- rep('', J)
        for(i in 1:J){
            if(K[i] > 2) itemtype[i] <- 'graded'
            if(K[i] == 2) itemtype[i] <- '2PL'                            
        }        
    } 
    if(length(itemtype) == 1) itemtype <- rep(itemtype, J)
    if(length(itemtype) != J) stop('itemtype specification is not the correct length')    
    itemloc <- cumsum(c(1,K))	
    model <- matrix(model$x,ncol=2)
    factorNames <- setdiff(model[,1],keywords)
    nfactNames <- length(factorNames)
    nfact <- sum(!grepl('\\(',factorNames))
    index <- 1:J	
    fulldata <- matrix(0,N,sum(K))
    Names <- NULL
    for(i in 1:J)
        Names <- c(Names, paste("Item.",i,"_",1:K[i],sep=""))				
    colnames(fulldata) <- Names			
    for(i in 1:J){
        ind <- index[i]		
        dummy <- matrix(0,N,K[ind])
        for (j in 0:(K[ind]-1))  
            dummy[,j+1] <- as.integer(data[,ind] == uniques[[ind]][j+1])  		
        fulldata[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- dummy		
    }	
    fulldata[is.na(fulldata)] <- 0    
    pats <- apply(fulldata, 1, paste, collapse = "/") 
    freqs <- rev(table(pats))
    nfreqs <- length(freqs)
    r <- as.vector(freqs)	
    tabdata <- unlist(strsplit(cbind(names(freqs)), "/"))
    tabdata <- matrix(as.numeric(tabdata), nfreqs, sum(K), TRUE)	
    tabdata2 <- matrix(NA, nfreqs, J)
    tmp <- c()
    for(i in 1:J){ 
        if(K[i] == 2) tmp <- c(tmp,0,1)
        else tmp <- c(tmp, 1:K[i])
    }
    for(i in 1:nfreqs){
        if(sum(tabdata[i, ]) < J){
            tmp2 <- rep(NA,J)
            ind <- tmp[as.logical(tabdata[i, ])]
            logicalind <- as.logical(tabdata[i, ])
            k <- 1
            for(j in 1:J){
                if(sum(logicalind[itemloc[j]:(itemloc[j+1]-1)]) != 0){
                    tmp2[j] <- ind[k]
                    k <- k + 1
                }
            }
            tabdata2[i, ] <- tmp2
        } else tabdata2[i, ] <- tmp[as.logical(tabdata[i, ])]
    }
    tabdata <- cbind(tabdata,r) 
    tabdata2 <- cbind(tabdata2,r)
    colnames(tabdata) <- c(Names,'Freq')	
    colnames(tabdata2) <- c(itemnames, 'Freq')             
    pars <- model.elements(model=model, itemtype=itemtype, factorNames=factorNames, 
                           nfactNames=nfactNames, nfact=nfact, J=J, K=K, fulldata=fulldata, 
                           itemloc=itemloc, data=data, N=N, guess=guess, upper=upper,  
                           itemnames=itemnames, exploratory=exploratory, constrain=constrain,
                           startvalues=startvalues, freepars=freepars, parprior=parprior, 
                           parnumber=parnumber, BFACTOR=BFACTOR, debug=debug)   
    prodlist <- attr(pars, 'prodlist')
    if(is(pars[[1]], 'numeric') || is(pars[[1]], 'logical')){
        names(pars) <- c(itemnames, 'Group_Parameters')
        attr(pars, 'parnumber') <- NULL
        return(pars)  
    }     
    onePLconstraint <- c()
    if(itemtype[1] == '1PL'){
        constrain <- list()
        for(i in 1:J)
            onePLconstraint <- c(onePLconstraint, pars[[i]]@parnum[1])    
        constrain[[length(constrain) + 1]] <- onePLconstraint
        pars <- model.elements(model=model, itemtype=itemtype, factorNames=factorNames, 
                               nfactNames=nfactNames, nfact=nfact, J=J, K=K, fulldata=fulldata, 
                               itemloc=itemloc, data=data, N=N, guess=guess, upper=upper,  
                               itemnames=itemnames, exploratory=exploratory, constrain=constrain,
                               startvalues=startvalues, freepars=freepars, parprior=parprior, 
                               parnumber=parnumber, BFACTOR=BFACTOR, debug=debug)
    }        
    if(any(itemtype == c('rsm','grsm'))){          
        unique.rsmgroups <- unique(na.omit(rsm.group))        
        for(group in unique.rsmgroups){                
            Kk <- unique(K[rsm.group[rsm.group == unique.rsmgroups[group]]])
            if(length(Kk) > 1) stop('Rating scale models require that items to have the 
                                       same number of categories')
            for(k in 2:(Kk-1)){
                rsmConstraint <- c()    
                for(i in 1:J){
                    if(rsm.group[i] == unique.rsmgroups[group]){
                        if(length(rsmConstraint) == 0){ 
                            pars[[i]]@est[length(pars[[i]]@est)] <- FALSE
                            rsmConstraint <- c(rsmConstraint, pars[[i]]@parnum[nfact+k])
                        } else rsmConstraint <- c(rsmConstraint, pars[[i]]@parnum[nfact+k])    
                    }
                }
                constrain[[length(constrain) + 1]] <- rsmConstraint
            }
        }            
    }
    npars <- 0
    for(i in 1:length(pars))
        npars <- npars + sum(pars[[i]]@est)    
    if(is.null(constrain)) constrain <- list()
    if(is.null(prodlist)) prodlist <- list()
    ret <- list(pars=pars, npars=npars, constrain=constrain, prodlist=prodlist, itemnames=itemnames,
                K=K, fulldata=fulldata, nfactNames=nfactNames, nfact=nfact, npars=npars, 
                exploratory=exploratory, J=J, itemloc=itemloc, factorNames=factorNames, 
                itemtype=itemtype, tabdata=tabdata, tabdata2=tabdata2, nLambdas=nfact+length(prodlist))
    ret
}
