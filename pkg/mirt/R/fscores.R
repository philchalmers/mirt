setMethod(
	f = "fscores",
	signature = 'mirtClass',
	definition = function(object, full.scores = FALSE, method = "EAP", ...)
	{    
		rotF <- object@F
		cs <- sqrt(1 - object@h2)
		a <- as.matrix(rotF / cs)
		d <- qnorm(object@facility) / cs  
		g <- object@guess  
		nfact <- ncol(a)
		theta <- as.matrix(seq(-4,4,length.out = 15))
		Theta <- thetaComb(theta,nfact)
		fulldata <- object@fulldata  
		tabdata <- object@tabdata[ ,1:ncol(fulldata)]
		colnames(tabdata) <- colnames(fulldata) 
		scores <- matrix(0,ncol=ncol(Theta),nrow=nrow(tabdata))
		thetas <- rep(0,nfact)
		W <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
		for (i in 1:nrow(scores)) {
			L <- 0  
			for (j in 1:nrow(a)){
				if(tabdata[i,j] == 1) L <- log(P.mirt(a[j, ],d[j],Theta,g[j])) + L
				else L <- log(1 - P.mirt(a[j, ],d[j],Theta,g[j])) + L	
			}	
			for (k in 1:ncol(Theta))
				thetas[k] <- sum(Theta[ ,k] * exp(L) * W / sum(exp(L) * W))
			scores[i, ] <- thetas
		}
		if(method == "MAP"){
			for (i in 1:nrow(scores)) {       
				Theta <- scores[i, ]	  
				thetas <- nlm(MAP.mirt,Theta,a=a,d=d,guess=g,patdata=tabdata[i, ])@estimate 
				scores[i, ] <- thetas
			}  
		}
		colnames(scores) <- paste("F",1:ncol(scores),sep="")  
		if (full.scores){      
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=ncol(Theta)) 
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tabdata[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, ] <- scores[j, ]
			}              
			return(cbind(fulldata,scoremat))
		} else {
			r <- as.matrix(object@tabdata[ ,ncol(tabdata)+1])
			colnames(r) <- 'Freq'			
			cat("\nMethod: ", method,"\n\n")	
			return(cbind(tabdata,r,scores))
		}   
	}  
)


setMethod(
	f = "fscores",
	signature = 'bfactorClass',
	definition = function(object, full.scores = FALSE, method = "EAP", ...){  
		g <- object@guess
		d <- object@pars[ ,ncol(object@pars)]
		a <- as.matrix(object@pars[ ,1:(ncol(object@pars) - 1)])
		nfact <- 2
		theta <- as.matrix(seq(-4,4,length.out = 15))
		Theta <- thetaComb(theta,nfact)
		fulldata <- object@fulldata  
		tabdata <- object@tabdata[ ,1:ncol(fulldata)]
		colnames(tabdata) <- colnames(fulldata) <- object@itemnames
		scores <- matrix(0,ncol=ncol(Theta),nrow=nrow(tabdata))
		thetas <- rep(0,nfact)
		logicalfact <- object@logicalfact
		W <- dmvnorm(Theta,rep(0,nfact),diag(nfact))  
		for (i in 1:nrow(scores)){
			L <- 0  
			for (j in 1:nrow(a)){
				if(tabdata[i,j] == 1) 
					L <- log(P.bfactor(a[j, ],d[j],Theta,g[j],logicalfact[j, ])) + L
				else 
					L <- log(1 - P.bfactor(a[j, ],d[j],Theta,g[j],logicalfact[j, ])) + L	
			}	
			for (k in 1:ncol(Theta))
				thetas[k] <- sum(Theta[ ,k] * exp(L) * W / sum(exp(L) * W))
			scores[i, ] <- thetas
		}
		if(method == "MAP"){
			for (i in 1:nrow(scores)) {       
				Theta <- scores[i, ]	  
				thetas <- nlm(MAP.bfactor,Theta,a=a,d=d,guess=g,
				patdata=tabdata[i, ],logicalfact=logicalfact)@estimate 
				scores[i, ] <- thetas
			}  
		}    
		colnames(scores) <- c("g","s")
		if (full.scores){
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=2) 
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tabdata[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, ] <- scores[j, ]
			} 
			return(cbind(fulldata,scoremat))
		} else {  
			cat("Method: ", method,"\n")
			return(cbind(tabdata,scores))
		}   
	}  
)