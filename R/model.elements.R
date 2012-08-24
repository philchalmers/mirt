model.elements <- function(model, factorNames, itemtype, nfactNames, nfact, J, K, fulldata, 
                           itemloc, data, N, guess, upper, itemnames, exploratory, constrain, 
                           startvalues, freepars, parprior, parnumber, BFACTOR = FALSE, debug)
{       
    if(debug == 'model.elements') browser()
    hasProdTerms <- ifelse(nfact == nfactNames, FALSE, TRUE)
    prodlist <- NULL
    if(hasProdTerms){
        tmp <- factorNames[grepl('\\(',factorNames)]
        tmp2 <- factorNames[!grepl('\\(',factorNames)] 
        tmp <- gsub("\\(","",tmp)    
        tmp <- gsub("\\)","",tmp)
        tmp <- gsub(" ","",tmp)
        prodlist <- strsplit(tmp,"\\*")
        for(j in 1:length(prodlist)){
            for(i in 1:nfact)
                prodlist[[j]][prodlist[[j]] == tmp2[[i]]] <- i    	
            prodlist[[j]] <- as.numeric(prodlist[[j]])	
        }		
    }  
    #slopes specification
    estlam <- matrix(FALSE, ncol = nfactNames, nrow = J)	
    for(i in 1:nfactNames){
        tmp <- model[model[ ,1] == factorNames[i],2]
        if(any(regexpr(",",tmp)))
            tmp <- strsplit(tmp,",")[[1]]
        popout <- c()	
        for(j in 1:length(tmp)){
            if(regexpr("-",tmp[j]) > 1){
                popout <- c(popout,j)
                tmp2 <- as.numeric(strsplit(tmp[j],"-")[[1]])
                tmp2 <- as.character(tmp2[1]:tmp2[2])
                tmp <- c(tmp,tmp2)
            }
        }
        if(length(popout != 0))	
            estlam[as.numeric(tmp[-popout]),i] <- TRUE
        else 
            estlam[as.numeric(tmp),i] <- TRUE
    }
    lambdas <- ifelse(estlam, .5, 0)	  
    #INT
    cs <- sqrt(abs(1-rowSums(lambdas^2)))	
    zetas <- list()
    loc <- 1	
    for(i in 1:J){        
        if(K[i] == 2){
            div <- ifelse(cs[i] > .25, cs[i], .25)		
            zetas[[i]] <- qnorm(mean(fulldata[,itemloc[i]]))/div            
        } else {			
            temp <- table(data[,i])[1:(K[i]-1)]/N
            temp <- cumsum(temp)
            div <- ifelse(cs[i] > .25, cs[i], .25)		
            zetas[[i]] <- qnorm(1 - temp)/div	            
        }		
    }    
    estzetas <- list()        
    for(i in 1:J)
        estzetas[[i]] <- length(zetas[[i]])                
    #COV
    find <- 1:nfact
    estgcov <- matrix(FALSE,nfact,nfact)    
    if(any(model[,1] == 'COV')){
        tmp <- model[model[,1] == 'COV',2]		
        tmp <- strsplit(tmp,",")[[1]]
        tmp <- gsub(" ","",tmp)
        for(i in 1:length(tmp)){            
            tmp2 <- strsplit(tmp[i],"*",fixed=TRUE)[[1]]				
            ind1 <- find[tmp2[1] == factorNames]
            ind2 <- find[tmp2[2] == factorNames]
            estgcov[ind2,ind1] <- TRUE            	
        }
    }
    gcov <- ifelse(estgcov,.25,0) 
    diag(gcov) <- 1	  
    #MEAN
    gmeans <- rep(0, nfact)
    estgmeans <- rep(FALSE, nfact)    
    if(exploratory){        
        Rpoly <- cormod(na.omit(data),K,guess)
        FA <- psych::fa(Rpoly, nfact, rotate = 'none', warnings= FALSE, fm="minres")    
        loads <- unclass(loadings(FA))
        u <- FA$unique
        u[u < .001 ] <- .2
        cs <- sqrt(u)
        lambdas <- loads/cs                
    }
    ret <- LoadPars(itemtype=itemtype, itemloc=itemloc, lambdas=lambdas, zetas=zetas, guess=guess, upper=upper,
                    fulldata=fulldata, J=J, K=K, nfact=nfact, constrain=constrain, nfactNames=nfactNames,
                    startvalues=startvalues, freepars=freepars, parprior=parprior, parnumber=parnumber,
                    estLambdas=estlam, BFACTOR=BFACTOR, debug=debug)      
    ret[[length(ret) + 1]] <- LoadGroupPars(gmeans=gmeans, gcov=gcov, estgmeans=estgmeans, 
                                            estgcov=estgcov, parnumber=attr(ret, 'parnumber')+1,
                                            startvalues=startvalues, freepars=freepars, parprior=parprior,
                                            constrain=constrain, debug=debug)
    attr(ret, 'prodlist') <- prodlist     
    return(ret)    
}