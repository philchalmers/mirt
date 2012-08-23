PrepData <- function(data, model, itemtype, guess, upper, startvalues, constrain, freepars, 
                     parprior, verbose, calcLL, debug, technical)
{
    if(debug == 'PrepData') browser()
    itemnames <- colnames(data)
    keywords <- c('COV')
    data <- as.matrix(data)    	
    colnames(data) <- itemnames	
    J <- ncol(data)
    N <- nrow(data)
    exploratory <- FALSE
    if(is(model, 'numeric')){
        tmp <- tempfile('tempfile')
        cat(paste('F',1:model,' = 1-', J, "\n", sep=''), file=tmp)
        model <- confmirt.model(tmp, quiet = TRUE)
        exploratory <- TRUE
        unlink(tmp)
    }    
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
    if(length(itemtype) != J) stop('itemtype specification is not the correct length')
    if(length(itemtype) == 1) itemtype <- rep(itemtype, J)
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
    parnumber <- 1 #to be used later when looping over more than 1 group       
    pars <- model.elements(model=model, itemtype=itemtype, factorNames=factorNames, 
                           nfactNames=nfactNames, nfact=nfact, J=J, K=K, fulldata=fulldata, 
                           itemloc=itemloc, data=data, N=N, guess=guess, upper=upper,  
                           itemnames=itemnames, exploratory=exploratory, constrain=constrain,
                           startvalues=startvalues, freepars=freepars, parprior=parprior, 
                           parnumber=parnumber, debug=debug)   
    prodlist <- attr(pars, 'prodlist')
    if(is(pars[[1]], 'numeric') || is(pars[[1]], 'logical')){
        names(pars) <- c(itemnames, 'Group_Parameters')
        attr(pars, 'parnumber') <- NULL
        return(pars)  
    }
    if(!is.null(constrain) || !is.null(parprior)){
        if(any(constrain == 'index', parprior == 'index')){
            returnedlist <- list()                        
            for(i in 1:length(pars))
                returnedlist[[i]] <- pars[[i]]@parnum 
            names(returnedlist) <- c(itemnames, 'Group_Parameters')            
            return(returnedlist)
        }
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
                               parnumber=parnumber, debug=debug)
    }
    npars <- 0
    for(i in 1:length(pars))
        npars <- npars + sum(pars[[i]]@est)    
    if(is.null(constrain)) constrain <- list()
    if(is.null(prodlist)) prodlist <- list()
    ret <- list(pars=pars, npars=npars, constrain=constrain, prodlist=prodlist, itemnames=itemnames,
                K=K, fulldata=fulldata, nfactNames=nfactNames, nfact=nfact, npars=npars, 
                exploratory=exploratory, J=J, itemloc=itemloc, factorNames=factorNames)
    ret
}
