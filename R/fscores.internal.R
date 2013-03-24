setMethod(
	f = "fscores.internal",
	signature = 'ExploratoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", 
                          quadpts = NULL, response.vector = NULL, degrees = NULL, verbose = TRUE)
	{          
        if(!is.null(response.vector)){            
            if(!is.matrix(response.vector)) response.vector <- matrix(response.vector, nrow = 1)
            v <- response.vector
            newdata <- rbind(object@data, v)
            nfact <- object@nfact 
            newmod <- mirt(newdata, nfact, technical = list(TOL = 10))
            newmod@pars <- object@pars
            tabdata <- newmod@tabdata
            index <- rep(FALSE, nrow(tabdata))
            for(i in 1:nrow(v)){
                vfull <-  matrix(v[i, ], nrow(tabdata), ncol(v), byrow = TRUE)
                index[rowSums(tabdata[,1:ncol(v)] == vfull) == ncol(v)] <- TRUE                
            }            
            newmod@tabdata <- newmod@tabdata[index, , drop = FALSE]
            newmod@tabdatalong <- newmod@tabdatalong[index, , drop = FALSE]
            ret <- fscores(newmod, rotate=rotate, full.scores=FALSE, method=method, 
                           quadpts=quadpts, verbose=FALSE)
            ret <- ret[, -(ncol(ret) - nfact*2)]            
            return(ret)
        }
        if(method == 'EAPsum') return(EAPsum(object, full.scores=full.scores, 
                                             quadpts=quadpts))
        pars <- object@pars        
		K <- object@K        
        J <- length(K)        
        prodlist <- attr(pars, 'prodlist')
        nfact <- object@nfact
        nLambdas <- object@nfact
        itemloc <- object@itemloc
        gp <- ExtractGroupPars(object@pars[[length(itemloc)]])
        if(rotate != 'CONFIRMATORY'){
            so <- summary(object, rotate = rotate, verbose = FALSE)            
            a <- rotateLambdas(so)
            for(i in 1:J)
                object@pars[[i]]@par[1:nfact] <- a[i, ]            
            gp$gmeans <- rep(0, nfact)
            gp$gcov <- so$fcor
        }               			        
        if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))
		theta <- as.matrix(seq(-4,4,length.out = quadpts))
		ThetaShort <- Theta <- thetaComb(theta,nfact)         
        if(length(prodlist) > 0)        
            Theta <- prodterms(Theta,prodlist)
		fulldata <- object@data 
		tabdata <- object@tabdatalong
		tabdata <- tabdata[ ,-ncol(tabdata), drop = FALSE]
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)			                
		W <- mvtnorm::dmvnorm(ThetaShort,gp$gmeans,gp$gcov) 
		W <- W/sum(W)                		
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
		for (i in 1:nrow(tabdata)){				
			L <- rowSums(log(itemtrace)[ ,as.logical(tabdata[i,]), drop = FALSE])			
			thetas <- colSums(ThetaShort * exp(L) * W / sum(exp(L) * W))
			SE <- sqrt(colSums(t((t(ThetaShort) - thetas))^2 * exp(L) * W / sum(exp(L) * W)))	
			scores[i, ] <- thetas
			SEscores[i, ] <- SE
		}		        
		if(method == "MAP"){ 
			for (i in 1:nrow(scores)){       
				tmp <- scores[i, ]	  
                estimate <- try(nlm(MAP.mirt,tmp,pars=pars, patdata=tabdata[i, ],
                                itemloc=itemloc, gp=gp, prodlist=prodlist, hessian=TRUE))
                if(is(estimate, 'try-error')) {
                    scores[i, ] <- SEscores[i, ] <- NA
                    next
                }
				scores[i, ] <- estimate$estimate
				SEest <- try(sqrt(diag(solve(estimate$hessian))))
				if(is(SEest, 'try-error')) SEest <- rep(NA, nfact)
				SEscores[i, ] <- SEest
			}  
		}
		if(method == "ML"){            		            
			tmp2 <- tabdata[,itemloc[-1] - 1, drop = FALSE]	
            tmp2[is.na(rowSums(object@tabdata))] <- 0
			scores[rowSums(tmp2) == J,] <- Inf
            tmp2 <- tabdata[,itemloc[-length(itemloc)], drop = FALSE]
            tmp2[is.na(rowSums(object@tabdata))] <- 1
            scores[rowSums(tmp2) == J,] <- -Inf
			SEscores[is.na(scores[,1]), ] <- rep(NA, nfact)
			for (i in 1:nrow(scores)){
				if(any(scores[i, ] %in% c(-Inf, Inf))) next 
				Theta <- scores[i, ]	  
				estimate <- try(nlm(MAP.mirt,Theta,pars=pars,patdata=tabdata[i, ],
				    itemloc=itemloc, gp=gp, prodlist=prodlist, ML=TRUE, hessian = TRUE))
				if(is(estimate, 'try-error')) {
				    scores[i, ] <- SEscores[i, ] <- NA
				    next
				}
				scores[i, ] <- estimate$estimate                
                SEest <- try(sqrt(diag(solve(estimate$hessian))))
                if(is(SEest, 'try-error')) SEest <- rep(NA, nfact)
				SEscores[i, ] <- SEest
			}  			
		}
        if(method == 'WLE'){                            
            tmp2 <- tabdata[,itemloc[-1] - 1, drop = FALSE] 
            tmp2[is.na(rowSums(object@tabdata))] <- 0
            scores[rowSums(tmp2) == J,] <- Inf
            tmp2 <- tabdata[,itemloc[-length(itemloc)], drop = FALSE]
            tmp2[is.na(rowSums(object@tabdata))] <- 1
            scores[rowSums(tmp2) == J,] <- -Inf            
            SEscores <- matrix(NA, nrow(SEscores), ncol(SEscores))
            for (i in 1:nrow(scores)){
                if(any(scores[i, ] %in% c(-Inf, Inf))) next
                Theta <- scores[i, ]	  
                estimate <- try(nlm(gradnorm.WLE,Theta,pars=pars,patdata=tabdata[i, ],
                                itemloc=itemloc, gp=gp, prodlist=prodlist, degrees=degrees,
                                hessian = TRUE))
                if(is(estimate, 'try-error')) {
                    scores[i, ] <- NA
                    next
                }
                scores[i, ] <- estimate$estimate                                
            }  
        }
		colnames(scores) <- paste('F', 1:ncol(scores), sep='')          
		if (full.scores){                                       
            tabdata2 <- object@tabdatalong
            tabdata2 <- tabdata2[,-ncol(tabdata2)]
            stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
            sfulldata <- apply(object@fulldata, 1, paste, sep='', collapse = '/')
            scoremat <- scores[match(sfulldata, stabdata2), , drop = FALSE]                                    
			colnames(scoremat) <- colnames(scores)	
			return(cbind(fulldata,scoremat))
		} else {            
            r <- object@tabdata[,ncol(object@tabdata)]            
            T <- E <- matrix(NA, 1, ncol(scores))
            for(i in 1:nrow(scores)){
                if(any(scores[i, ] %in% c(Inf, -Inf))) next
                T <- rbind(T, matrix(rep(scores[i, ], r[i]), ncol=ncol(scores), byrow = TRUE))
                E <- rbind(E, matrix(rep(SEscores[i, ], r[i]), ncol=ncol(scores), byrow = TRUE))
            }            
            T <- na.omit(T)
            E <- na.omit(E)
            reliability <- diag(var(T)) / (diag(var(T)) + colMeans(E^2))
            names(reliability) <- colnames(scores)
			if(verbose){
                cat("\nMethod: ", method)                
                cat("\n\nEmpirical Reliability:\n")
                print(round(reliability, 4))                
			}
			colnames(SEscores) <- paste('SE_', colnames(scores), sep='')
			return(cbind(object@tabdata,scores,SEscores))
		}   
	}  
)

#------------------------------------------------------------------------------
setMethod(
	f = "fscores.internal",
	signature = 'ConfirmatoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", 
	                      quadpts = NULL, response.vector = NULL, degrees = NULL, verbose = TRUE)
	{ 	        
        class(object) <- 'ExploratoryClass'
        ret <- fscores(object, rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, quadpts=quadpts, 
                       response.vector=response.vector, degrees=degrees, verbose=verbose)
        return(ret)
	}	
)

#------------------------------------------------------------------------------
setMethod(
    f = "fscores.internal",
    signature = 'MultipleGroupClass',
    definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", 
                          quadpts = NULL, response.vector = NULL, degrees = NULL, verbose = TRUE)
    { 	        
        cmods <- object@cmods
        ngroups <- length(cmods)
        for(g in 1:ngroups)
            class(cmods[[g]]) <- 'ConfirmatoryClass'
        ret <- vector('list', length(cmods))
        for(g in 1:ngroups)
            ret[[g]] <- fscores(cmods[[g]], rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, 
                           quadpts=quadpts, degrees=degrees, verbose=verbose)
        names(ret) <- object@groupNames
        if(full.scores){
            id <- c()
            fulldata <- matrix(NA, 1, ncol(ret[[1]]))
            for(g in 1:ngroups){
                id <- c(id, rownames(ret[[g]]))
                fulldata <- rbind(fulldata, ret[[g]])
            }
            fulldata <- fulldata[-1, ]
            fulldata <- data.frame(id=as.numeric(id), fulldata)
            ret <- fulldata[order(fulldata$id), ]
            ret <- ret[ ,-1]        
        }         
        return(ret)
    }	
)

# MAP scoring for mirt
MAP.mirt <- function(Theta, pars, patdata, itemloc, gp, prodlist, ML=FALSE)
{       
    ThetaShort <- Theta
    Theta <- matrix(Theta, nrow=1)
    if(length(prodlist) > 0)
        Theta <- prodterms(Theta,prodlist)
    itemtrace <- matrix(0, ncol=length(patdata), nrow=nrow(Theta))        
    for (i in 1:(length(itemloc)-1))
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)		
    itemtrace[itemtrace < 1e-8] <- 1e-8
    L <- sum(log(itemtrace)[as.logical(patdata)])    
    prior <- mvtnorm::dmvnorm(ThetaShort, gp$gmeans, gp$gcov)
    L <- ifelse(ML, -L, (-1)*(L + log(prior)))
    L  
}

gradnorm.WLE <- function(Theta, pars, patdata, itemloc, gp, prodlist, degrees){    
    ThetaShort <- Theta
    Theta <- matrix(Theta, nrow=1)
    if(length(prodlist) > 0)
        Theta <- prodterms(Theta,prodlist)
    nfact <- ncol(Theta)
    itemtrace <- matrix(0, ncol=length(patdata), nrow=nrow(Theta))  
    dP <- d2P <- vector('list', nfact)
    for(i in 1:nfact)
        dP[[i]] <- d2P[[i]] <- itemtrace 
    I <- numeric(1)
    dW <- dL <- numeric(nfact)
    for (i in 1:(length(itemloc)-1)){
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)		
        deriv <- DerivTheta(x=pars[[i]], Theta=Theta)
        for(k in 1:nfact){
            dPitem <- d2Pitem <- matrix(0, 1, length(deriv[[1]]))
            for(j in 1:length(deriv[[1]])){
                dPitem[1, j] <- deriv$grad[[j]][ ,k]
                d2Pitem[1, j] <- deriv$hess[[j]][ ,k]           
            }
            dP[[k]][ ,itemloc[i]:(itemloc[i+1] - 1)] <- dPitem
            d2P[[k]][ ,itemloc[i]:(itemloc[i+1] - 1)] <- d2Pitem            
            dW[k] <- dW[k] + sum(dPitem * d2Pitem / itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)])
        }        
        I <- I + iteminfo(x=pars[[i]], Theta=Theta, degrees=degrees)       
    }
    dW <- 1/(2*I^2) * dW
    for(i in 1:nfact)
        dL[i] <- sum(patdata * dP[[i]] / itemtrace)
    grad <- dL - dW*I    
    MIN <- sum(grad^2)
    MIN
}

EAPsum <- function(x, full.scores = FALSE, quadpts = NULL, S_X2 = FALSE){    
    calcL1 <- function(itemtrace, K){        
        J <- length(K)
        L0 <- L1 <- matrix(1, sum(K-1) + 1, ncol(itemtrace))                
        L0[1:K[1], ] <- itemtrace[1:K[1], ]
        nstar <- K[1] + K[2] - 3    
        Sum.Scores <- 1:nrow(L0)-1    
        MAX.Scores <- cumsum(K-1)            
        for(i in 1:(J-1)){
            T <- itemtrace[itemloc[i+1]:(itemloc[i+2] - 1), ]
            L1[1, ] <- L0[1, ] * T[1, ]        
            #recursive rule for internal values (gets a little ugly at polytomous data edges)
            for(j in 1:nstar+1){        
                sums <- 0                            
                for(k in 1:K[i+1]-1)
                    if(Sum.Scores[j] >= k && (MAX.Scores[i] + k) >= Sum.Scores[j])
                        sums <- sums + L0[j - k, ] * T[1 + k, ]
                L1[j, ] <- sums                   
            }        
            L1[j+1, ] <- L0[j - k + 1, ] * T[nrow(T), ]
            L0 <- L1        
            nstar <- nstar + K[i+1] - 1
        }        
        list(L1=L1, Sum.Scores=Sum.Scores) 
    }
    if(x@nfact > 1) stop('EAP sum score method only is applicable to unidimensional models')             
    if(is.null(quadpts)) quadpts <- 40
    Theta <- as.matrix(seq(-4,4,length.out = quadpts))    
    prior <- dnorm(Theta)
    prior <- prior/sum(prior)
    pars <- x@pars
    K <- x@K
    J <- length(K)
    itemloc <- x@itemloc    
    itemtrace <- matrix(0, ncol=ncol(x@tabdatalong)-1, nrow=nrow(Theta))        
    for (i in 1:J)
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)    
    itemtrace <- t(itemtrace)    
    tmp <- calcL1(itemtrace=itemtrace, K=K)    
    L1 <- tmp$L1
    Sum.Scores <- tmp$Sum.Scores            
    if(S_X2){        
        L1total <- L1 %*% prior         
        Elist <- vector('list', J)        
        for(i in 1:J){
            KK <- K[-i]
            T <- itemtrace[c(itemloc[i]:(itemloc[i+1]-1)), ]
            itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i+1]-1)), ]
            tmp <- calcL1(itemtrace=itemtrace2, K=KK)    
            E <- matrix(NA, nrow(L1total), nrow(T))            
            for(j in 1:(nrow(T)))                
                E[1:nrow(tmp$L1)+j-1,j] <- tmp$L1 %*% (T[j,] * prior) / 
                    L1total[1:nrow(tmp$L1)+j-1, ]                                        
            Elist[[i]] <- E[-c(1, nrow(E)), ]
        }
        return(Elist)
    }
    thetas <- SEthetas <- numeric(nrow(L1))    
    for(i in 1:length(thetas)){        
        thetas[i] <- sum(Theta * L1[i, ] * prior / sum(L1[i,] * prior))
        SEthetas[i] <- sqrt(sum((Theta - thetas[i])^2 * L1[i, ] * prior / sum(L1[i,] * prior)))
    }
    ret <- data.frame(Sum.Scores=Sum.Scores, Theta=thetas, SE.Theta=SEthetas)        
    rownames(ret) <- ret$Sum.Scores
    if(full.scores){               
        if(any(is.na(x@data))) stop('Full scores requires a complete dataset (no N\'s)')
        dat <- x@data
        adj <- apply(dat, 2, min)
        if(any(adj > 0)) message('Data adjusted so that every item has a lowest score of 0')
        dat <- t(t(dat) - adj)
        scores <- rowSums(dat)
        EAPscores <- ret$Theta[match(scores, ret$Sum.Scores)]
        ret <- data.frame(Sum.Scores=scores, Theta=EAPscores)
    }
    ret   
}
