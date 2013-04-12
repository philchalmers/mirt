PrepData <- function(data, model, itemtype, guess, upper, 
                     parprior, verbose, technical, parnumber = 1, BFACTOR = FALSE,
                     grsm.block = NULL, rsm.block = NULL, D, mixedlist, customItems, 
                     fulldata = NULL, key)
{    
    if(is.null(grsm.block)) grsm.block <- rep(1, ncol(data))
    if(is.null(rsm.block)) rsm.block <- rep(1, ncol(data))
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
        model <- confmirt.model(file=tmp, quiet = TRUE)        
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
    if(is.null(key) && any(itemtype %in% c('2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')))
        stop('When using nested logit items a scoring key must be provided with key = c(...)')
    if(is.null(key))  key <- rep(1L, J)
    if(length(key) != J)
        stop("The number of elements in the key input is incorrect.")
    key <- as.integer(key)
    uniques <- list()
    for(i in 1:J)
        uniques[[i]] <- sort(unique(data[,i]))
    K <- rep(0,J)
    for(i in 1:J) K[i] <- length(uniques[[i]])	
    guess[K > 2] <- 0
    upper[K > 2] <- 1
    if(any(K < 2))
        stop('The following items have only one response category and cannot be estimated: ', 
             paste(itemnames[K < 2], ''))    
    if(is.null(itemtype)) {
        itemtype <- rep('', J)
        for(i in 1:J){
            if(K[i] > 2) itemtype[i] <- 'graded'
            if(K[i] == 2) itemtype[i] <- '2PL'                            
        }        
    } 
    if(length(itemtype) == 1) itemtype <- rep(itemtype, J)
    if(length(itemtype) != J) stop('itemtype specification is not the correct length')    
    guess[guess == 0 & itemtype %in% c('3PL', '4PL', 'PC3PL')] <- .15
    upper[upper == 1 & itemtype %in% c('4PL', '3PLu')] <- .85
    itemloc <- cumsum(c(1,K))	
    model <- matrix(model$x,ncol=2)
    factorNames <- setdiff(model[,1],keywords)
    nfactNames <- length(factorNames)
    nfact <- sum(!grepl('\\(',factorNames))
    index <- 1:J
    Names <- NULL
    for(i in 1:J)
        Names <- c(Names, paste("Item.",i,"_",1:K[i],sep=""))    			
    if(is.null(fulldata)){
        fulldata <- matrix(0,N,sum(K))    
        colnames(fulldata) <- Names			
        for(i in 1:J){
            ind <- index[i]		
            dummy <- matrix(0,N,K[ind])
            for (j in 0:(K[ind]-1))  
                dummy[,j+1] <- as.integer(data[,ind] == uniques[[ind]][j+1])  		
            fulldata[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- dummy		
        }	
        fulldata[is.na(fulldata)] <- 0        
    }
    pars <- model.elements(model=model, itemtype=itemtype, factorNames=factorNames, 
                           nfactNames=nfactNames, nfact=nfact, J=J, K=K, fulldata=fulldata, 
                           itemloc=itemloc, data=data, N=N, guess=guess, upper=upper,  
                           itemnames=itemnames, exploratory=exploratory, parprior=parprior, 
                           parnumber=parnumber, BFACTOR=BFACTOR, D=D, mixedlist=mixedlist, 
                           customItems=customItems, key=key)   
    prodlist <- attr(pars, 'prodlist')
    if(is(pars[[1]], 'numeric') || is(pars[[1]], 'logical')){
        names(pars) <- c(itemnames, 'Group_Parameters')
        attr(pars, 'parnumber') <- NULL
        return(pars)  
    }   
    #within group constraints
    constrain <- list()
    onePLconstraint <- c()    
    if(itemtype[1] == '1PL'){        
        for(i in 1:J)
            onePLconstraint <- c(onePLconstraint, pars[[i]]@parnum[1 + pars[[1]]@nfixedeffects])    
        constrain[[length(constrain) + 1]] <- onePLconstraint
    }        
    if(any(itemtype == 'grsm')){           
        unique.grsmgroups <- unique(na.omit(grsm.block))        
        for(group in unique.grsmgroups){                
            Kk <- unique(K[grsm.block[grsm.block == unique.grsmgroups[group]]])
            if(length(Kk) > 1) stop('Rating scale models require that items to have the 
                                       same number of categories')
            for(k in 1:(Kk-1)){
                grsmConstraint <- c()    
                for(i in 1:J){
                    if(grsm.block[i] == unique.grsmgroups[group]){
                        if(length(grsmConstraint) == 0){ 
                            pars[[i]]@est[length(pars[[i]]@est)] <- FALSE
                            grsmConstraint <- c(grsmConstraint, pars[[i]]@parnum[nfact+k])
                        } else grsmConstraint <- c(grsmConstraint, pars[[i]]@parnum[nfact+k])    
                    }
                }
                constrain[[length(constrain) + 1]] <- grsmConstraint
            }
        }            
    }  
    if(any(itemtype == 'rsm')){         
        unique.rsmgroups <- unique(na.omit(rsm.block))        
        for(group in unique.rsmgroups){                
            Kk <- unique(K[rsm.block[rsm.block == unique.rsmgroups[group]]])
            if(length(Kk) > 1) stop('Rating scale models require that items to have the 
                                    same number of categories')
            for(k in 1:(Kk-1)){
                rsmConstraint <- c()    
                for(i in 1:J){
                    if(rsm.block[i] == unique.rsmgroups[group]){
                        if(length(rsmConstraint) == 0){ 
                            pars[[i]]@est[length(pars[[i]]@est)] <- FALSE
                            rsmConstraint <- c(rsmConstraint, pars[[i]]@parnum[nfact+k+1])
                        } else rsmConstraint <- c(rsmConstraint, pars[[i]]@parnum[nfact+k+1])    
                    }
                }
                constrain[[length(constrain) + 1]] <- rsmConstraint
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
