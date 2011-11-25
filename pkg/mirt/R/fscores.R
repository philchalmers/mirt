#' Methods for Function fscores
#' 
#' Save tabulated or full data factor scores for \code{mirt} or \code{bfactor}
#' using EAP or MAP scoring.
#' 
#' 
#' @name fscores-methods
#' @aliases fscores-methods fscores,bfactorClass-method
#' fscores,mirtClass-method fscores,polymirtClass-method
#' fscores,confmirtClass-method
#' @docType methods
#' @section Methods: \describe{ \item{fscores}{\code{signature(object =
#' "bfactorClass")}} \item{fscores}{\code{signature(object = "mirtClass")}}
#' \item{fscores}{\code{signature(object = "polymirtClass")}}
#' \item{fscores}{\code{signature(object = "confmirtClass")}} }
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportMethod fscores
#' @keywords methods
setGeneric("fscores", 
	def = function(object, ...) standardGeneric("fscores")
)

#' Compute factor scores
#' 
#' Computes MAP or EAP factor scores for \code{mirt} and \code{bfactor} models,
#' or stocastic approximations for \code{polymirt} and \code{confmirt}. Note
#' that only the general factor scores are computed for bifactor models.
#' 
#' 
#' @aliases fscores fscores,mirt-method fscores,bfactor-method
#' fscores,polymirt-method fscores,confmirt-method fscores
#' @param object a model of class \code{mirtClass} or \code{bfactorClass}
#' @param full.scores if \code{FALSE} (default) then a summary table with
#' factor scores for each unique pattern is displayed. Otherwise the original
#' data matrix is returned with the computed factor scores appended to the
#' rightmost column
#' @param method type of factor score estimation method. Can be expected
#' a-posteriori (\code{"EAP"}) or Bayes modal (\code{"MAP"})
#' @param ndraws number of MH samplers to draw for each response pattern
#' @param thin controls how much the chain should be thinned by, default
#' collects every 5th draw. Note that \code{ndraws/thin} must be a whole number
#' @param ... additional arguments to be passed
#' @return Returns either a summary table with the response patterns and
#' expected factor scores, or a complete data matrix with factor scores
#' appended to the last column.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords factor.scores
#' @exportMethod fscores
#' @export fscores
#' @examples
#' 
#' \dontrun{
#' 
#' tabscores <- fscores(mod)
#' fullscores <- fscores(mod, full.scores = TRUE)
#' 
#' 
#'   }
#' 
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
		SEscores <- scores <- matrix(0,ncol=ncol(Theta),nrow=nrow(tabdata))
		SE <- thetas <- rep(0,nfact)
		W <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
		W <- W/sum(W)
		for (i in 1:nrow(scores)){
			L <- 0  
			for (j in 1:nrow(a)){
				if(tabdata[i,j] == 1) L <- log(P.mirt(a[j, ],d[j],Theta,g[j])) + L
				else L <- log(1 - P.mirt(a[j, ],d[j],Theta,g[j])) + L	
			}	
			for (k in 1:ncol(Theta))
				thetas[k] <- sum(Theta[ ,k] * exp(L) * W / sum(exp(L) * W))
			for (k in 1:ncol(Theta))	
				SE[k] <- sqrt(sum((Theta[,k] - thetas[k])^2 * exp(L) * W / sum(exp(L) * W)))
			scores[i, ] <- thetas
			SEscores[i, ] <- SE
		}
		if(method == "MAP"){
			for (i in 1:nrow(scores)){       
				Theta <- scores[i, ]	  
				thetas <- nlm(MAP.mirt,Theta,a=a,d=d,guess=g,patdata=tabdata[i, ])$estimate 
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
			r <- matrix(object@tabdata[ ,ncol(tabdata)+1])
			colnames(r) <- 'Freq'			
			colnames(SEscores) <- paste("SE_",colnames(scores),sep='')
			cat("\nMethod: ", method,"\n\n")
			if(method == "EAP")	
				return(cbind(tabdata,r,scores,SEscores))
			else
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
		SEscores <- scores <- rep(0,nrow(tabdata))
		SE <- thetas <- 0
		logicalfact <- object@logicalfact
		W <- dmvnorm(Theta,rep(0,nfact),diag(nfact)) 
		W <- W/sum(W)	
		for (i in 1:nrow(tabdata)){
			L <- 0  
			for (j in 1:nrow(a)){
				if(tabdata[i,j] == 1) 
					L <- log(P.bfactor(a[j, ],d[j],Theta,g[j],logicalfact[j, ])) + L
				else 
					L <- log(1 - P.bfactor(a[j, ],d[j],Theta,g[j],logicalfact[j, ])) + L	
			}	
			thetas <- sum(Theta[ ,1] * exp(L) * W / sum(exp(L) * W))
			SE <- sqrt(sum((Theta[ ,1] - thetas)^2 * exp(L) * W / sum(exp(L) * W)))	
			scores[i] <- thetas
			SEscores[i] <- SE
		}
		if(method == "MAP"){
			for (i in 1:length(scores)) {       
				Theta <- scores[i]	  
				thetas <- nlm(MAP.bfactor,Theta,a=a,d=d,guess=g,
					patdata=tabdata[i, ],logicalfact=logicalfact)$estimate 
				scores[i] <- thetas
			}  
		}
		if(method == 'EAP'){	
			scores <- cbind(scores,SEscores)			
			colnames(scores) <- c("g","SE_g")
		} else {
			scores <- matrix(scores)
			colnames(scores) <- 'g'
		}	
		if(full.scores){
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=1) 
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tabdata[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, ] <- scores[j,1]
			} 
			return(cbind(fulldata,scoremat))
		} else {  	
			r <- matrix(object@tabdata[,ncol(tabdata)+1])
			colnames(r) <- 'Freq'
			cat("\nMethod: ", method,"\n\n")
			if(method == "EAP")	
				return(cbind(tabdata,r,scores))
			else
				return(cbind(tabdata,r,scores))	
		}   
	}  
)

setMethod(
	f = "fscores",
	signature = 'polymirtClass',
	definition = function(object, full.scores = FALSE, ndraws = 3000, thin = 5, ...)
	{ 	
		cand.t.var <- 1
		theta0 <- object@Theta
		K <- object@K
		nfact <- ncol(theta0)
		lambdas <- matrix(object@pars[,1:nfact],ncol=nfact)
		zetas <- na.omit(as.numeric(t(object@pars[,(nfact+1):ncol(object@pars)])))
		guess <- object@guess
		guess[is.na(guess)] <- 0
		data <- cbind(object@data,object@fulldata)
		Names <- c(colnames(object@data[,1:length(K)]),paste("F",1:nfact,sep=''),paste("SE_F",1:nfact,sep=''))
		tabdata <- unique(data)[,-c(1:length(K))]			
		itemloc <- object@itemloc
		Theta <- list()
		for(i in 1:nfact)
			Theta[[i]] <- matrix(0,ncol=ndraws/thin,nrow=nrow(tabdata))		
		theta0 <- matrix(0,nrow(tabdata),nfact)
		for(i in 1:30){			
			theta0 <- draw.thetas(theta0,lambdas,zetas,guess,tabdata,K,itemloc,cand.t.var)
			if(attr(theta0,'Proportion Accepted') > .4) cand.t.var <- cand.t.var + .2
			if(attr(theta0,'Proportion Accepted') < .3) cand.t.var <- cand.t.var - .2
		}
		ind <- 1
		for(i in 1:ndraws){			
			theta0 <- draw.thetas(theta0,lambdas,zetas,guess,tabdata,K,itemloc,cand.t.var)
			if(i %% thin == 0){
				for(j in 1:nfact)
					Theta[[j]][,ind] <- theta0[,j]									
				ind <- ind + 1
			}			
		}

		expscores <- matrix(0,ncol=nfact,nrow=nrow(tabdata))
		sdscores <- matrix(0,ncol=nfact,nrow=nrow(tabdata))
		for(i in 1:nfact){
			expscores[,i] <- rowMeans(Theta[[i]])
			sdscores[,i] <- apply(Theta[[i]],1,sd)
		}
				
		ret <- cbind(unique(data)[,1:length(K)],expscores,sdscores)
		colnames(ret) <- Names
		
		if(!full.scores){ 
			ret <- ret[order(expscores[,1]),]
			rownames(ret) <- NULL
			return(ret)
		} else {
			fulldata <- object@data
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=nfact)
			colnames(scoremat) <- paste("F",1:nfact,sep='')
			tmp <- unique(data)[,1:length(K)]
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tmp[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, ] <- expscores[j, ]
			}              
			return(cbind(object@data,scoremat))
		}	
	}	
)

setMethod(
	f = "fscores",
	signature = 'confmirtClass',
	definition = function(object, full.scores = FALSE, ndraws = 3000, thin = 5, ...)
	{ 	
		cand.t.var <- 1
		estComp <- object@estComp
		sig <- object@gpars$sig
		mu <- object@gpars$u
		theta0 <- object@Theta
		K <- object@K
		nfact <- ncol(theta0)
		lambdas <- matrix(object@pars[,1:nfact],ncol=nfact)
		lambdas[is.na(lambdas)] <- 0
		zetas <- na.omit(as.numeric(t(object@pars[,(nfact+1):ncol(object@pars)])))
		guess <- object@guess
		guess[is.na(guess)] <- 0
		data <- cbind(object@data,object@fulldata)		
		Names <- c(colnames(object@data[,1:length(K)]),paste("F",1:nfact,sep=''),paste("SE_F",1:nfact,sep=''))
		tabdata <- unique(data)[,-c(1:length(K))]			
		itemloc <- object@itemloc
		Theta <- list()
		for(i in 1:nfact)
			Theta[[i]] <- matrix(0,ncol=ndraws/thin,nrow=nrow(tabdata))		
		theta0 <- matrix(0,nrow(tabdata),nfact)
		for(i in 1:30){			
			theta0 <- draw.thetas(theta0,lambdas,zetas,guess,tabdata,K,itemloc,cand.t.var,sig,mu,estComp)
			if(attr(theta0,'Proportion Accepted') > .4) cand.t.var <- cand.t.var + .2
			if(attr(theta0,'Proportion Accepted') < .3) cand.t.var <- cand.t.var - .2
		}
		ind <- 1
		for(i in 1:ndraws){			
			theta0 <- draw.thetas(theta0,lambdas,zetas,guess,tabdata,K,itemloc,cand.t.var,sig,mu,estComp)
			if(i %% thin == 0){
				for(j in 1:nfact)
					Theta[[j]][,ind] <- theta0[,j]									
				ind <- ind + 1
			}			
		}

		expscores <- matrix(0,ncol=nfact,nrow=nrow(tabdata))
		sdscores <- matrix(0,ncol=nfact,nrow=nrow(tabdata))
		for(i in 1:nfact){
			expscores[,i] <- rowMeans(Theta[[i]])
			sdscores[,i] <- apply(Theta[[i]],1,sd)
		}
				
		ret <- cbind(unique(data)[,1:length(K)],expscores,sdscores)
		colnames(ret) <- Names
		
		if(!full.scores){ 
			ret <- ret[order(expscores[,1]),]
			rownames(ret) <- NULL
			return(ret)
		} else {
			fulldata <- object@data
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=nfact)
			colnames(scoremat) <- paste("F",1:nfact,sep='')
			tmp <- unique(data)[,1:length(K)]
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tmp[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, ] <- expscores[j, ]
			}              
			return(cbind(object@data,scoremat))
		}	
	}	
)


