#' Methods for Function fscores
#' 
#' Computes MAP, EAP, or ML factor scores for \code{mirt} and \code{bfactor} models,
#' or a stochastic approximation with a multivariate normal prior for \code{polymirt} and 
#' \code{confmirt}. Note that only the general factor scores are computed for bifactor 
#' models.
#'
#' 
#' @usage 
#' fscores(object, ...)
#' 
#' @aliases fscores-method fscores,bfactorClass-method
#' fscores,mirtClass-method fscores,polymirtClass-method
#' fscores,confmirtClass-method
#' @docType methods
#' @section Methods: \describe{ \item{fscores}{\code{signature(object =
#' "bfactorClass")}} \item{fscores}{\code{signature(object = "mirtClass")}}
#' \item{fscores}{\code{signature(object = "polymirtClass")}}
#' \item{fscores}{\code{signature(object = "confmirtClass")}} }
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @rdname fscores-methods   
#' @exportMethod fscores
#' @keywords methods
setGeneric("fscores", 
           def = function(object, ...) standardGeneric("fscores")
)

#------------------------------------------------------------------------------
# Methods for Function fscores
#
#' @name fscores
#' @param object a model of class \code{mirtClass}, \code{bfactorClass}, \code{polymirtClass},
#' or \code{confmirtClass}
#' @param full.scores if \code{FALSE} (default) then a summary table with
#' factor scores for each unique pattern is displayed. Otherwise the original
#' data matrix is returned with the computed factor scores appended to the
#' rightmost column
#' @param rotate rotation declaration to be used when estimating the factor scores. If \code{""} then the 
#' \code{object@@rotate} default value is used
#' @param method type of factor score estimation method. Can be expected
#' a-posteriori (\code{"EAP"}), Bayes modal (\code{"MAP"}), or maximum likelihood 
#' (\code{"ML"}). Only applicable to \code{mirtClass} and \code{bfactorClass} objects 
#' @param ndraws number of MH samples to draw for each response pattern for \code{polymirtClass}
#' or \code{confmirtClass} objects
#' @param thin controls how much the chain should be thinned by, default
#' collects every 5th draw (\code{thin = 5}). Note that \code{ndraws/thin} must be a whole number.
#' for \code{polymirtClass} or \code{confmirtClass} objects only
#' @param verbose logical; print verbose output messages?
#' @return Returns either a summary table with the response patterns and
#' expected factor scores, or a complete data matrix with factor scores
#' appended to the last column.
#' @keywords factor.scores
#' @rdname fscores-methods   
#' @export fscores
#' @examples
#' 
#' \dontrun{
#'
#' mod <- mirt(Science, 1)
#' tabscores <- fscores(mod)
#' fullscores <- fscores(mod, full.scores = TRUE)
#' fullscores <- fscores(mod, full.scores = TRUE, method='MAP')
#' 
#' mod2 <- polymirt(Science, 1)
#' tabscores2 <- fscores(mod2, ndraws = 5000)
#' 
#'   }
#' 
setMethod(
	f = "fscores",
	signature = 'mirtClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", verbose = TRUE)
	{        
		K <- object@K				
        so <- summary(object, rotate = rotate, print = FALSE)
        a <- rotateLambdas(so)
		d <- object@pars$zetas		
		g <- object@guess				
		u <- object@upper
		itemloc <- object@itemloc
		J <- nrow(a)
		nfact <- ncol(a)
		theta <- as.matrix(seq(-4,4,length.out = 15))
		Theta <- thetaComb(theta,nfact)
		fulldata <- object@data 
		tabdata <- object@tabdatalong
		tabdata <- tabdata[,-ncol(tabdata)]
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)			
		W <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact)) 
		W <- W/sum(W)
		itemtrace <- matrix(0, ncol=ncol(tabdata), nrow=nrow(Theta))
		for (i in 1:J){
			if(length(d[[i]]) == 1){
				itemtrace[ ,itemloc[i] + 1] <- P.mirt(a[i, ], d[[i]], Theta, g[i], u[i]) 
				itemtrace[ ,itemloc[i]] <- 1.0 - itemtrace[ ,itemloc[i] + 1]
			} else {
				itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- 
					P.poly(a[i, ], d[[i]], Theta, TRUE)	
			}
		}	        
		for (i in 1:nrow(tabdata)){				
			L <- rowSums(log(itemtrace)[ ,as.logical(tabdata[i,])])			
			thetas <- colSums(Theta * exp(L) * W / sum(exp(L) * W))
			SE <- sqrt(colSums(t((t(Theta) - thetas))^2 * exp(L) * W / sum(exp(L) * W)))	
			scores[i, ] <- thetas
			SEscores[i, ] <- SE
		}		
		if(method == "MAP"){ 
			for (i in 1:nrow(scores)){       
				tmp <- scores[i, ]	  
				thetas <- nlm(MAP.mirt,tmp,a=a,d=d,guess=g,upper=u,patdata=tabdata[i, ],
                              itemloc=itemloc)$estimate 
				scores[i, ] <- thetas
			}  
		}
		if(method == "ML"){
			tmp <- tabdata[,itemloc[-length(itemloc)]]			 
			tmp2 <- tabdata[,itemloc[-1] - 1]			 
            scores[rowSums(tmp) == J, ] <- Inf
			scores[rowSums(tmp2) == J,] <- -Inf
			for (i in 1:nrow(scores)){
				if(any((scores[i, ]) == -Inf | scores[i, ] == Inf)) next 
				Theta <- scores[i, ]	  
				thetas <- nlm(MAP.mirt,Theta,a=a,d=d,guess=g,upper=u,patdata=tabdata[i, ],
                              itemloc=itemloc, ML=TRUE)$estimate 
				scores[i, ] <- thetas
			}  
		}
		colnames(scores) <- colnames(object@F)
		if (full.scores){      
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=ncol(Theta))
			tabdata2 <- object@tabdata[,-(ncol(fulldata)+1)]	
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tabdata2[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, ] <- scores[j, ]
			} 
			colnames(scoremat) <- colnames(object@F)	
			return(cbind(fulldata,scoremat))
		} else {						
			if(verbose) cat("\nMethod: ", method,"\n\n")
			colnames(SEscores) <- paste('SE_', colnames(scores), sep='')
			return(cbind(object@tabdata,scores,SEscores))
				
		}   
	}  
)

#------------------------------------------------------------------------------
#' @rdname fscores-methods
setMethod(
	f = "fscores",
	signature = 'bfactorClass',
	definition = function(object, full.scores = FALSE, method = "EAP", verbose = TRUE)
	{  
		K <- object@K		
		a <- object@pars$lambdas		
		d <- object@pars$zetas		
		g <- object@guess				
		u <- object@upper
		itemloc <- object@itemloc
		J <- nrow(a)
		nfact <- 2
		theta <- as.matrix(seq(-4,4,length.out = 15))
		Theta <- thetaComb(theta,nfact)
		fulldata <- object@data 
		tabdata <- object@tabdatalong
		tabdata <- tabdata[,-ncol(tabdata)]
		SEscores <- scores <- rep(0,nrow(tabdata))
		SE <- thetas <- 0
		logicalfact <- object@logicalfact
		W <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact)) 
		W <- W/sum(W)
		itemtrace <- matrix(0, ncol=ncol(tabdata), nrow=nrow(Theta))
		for (i in 1:J){
			if(length(d[[i]]) == 1){
				itemtrace[ ,itemloc[i] + 1] <- P.bfactor(a[i, ], d[[i]], Theta, g[i], u[i], 
                                                         logicalfact[i, ]) 
				itemtrace[ ,itemloc[i]] <- 1.0 - itemtrace[ ,itemloc[i] + 1]
			} else {
				itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- 
					P.bfactor(a[i, ], d[[i]], Theta, 0, 1, logicalfact[i, ])	
			}
		}		        
		for (i in 1:nrow(tabdata)){							
			L <- rowSums(log(itemtrace)[ ,as.logical(tabdata[i,])])				
			thetas <- sum(Theta[ ,1] * exp(L) * W / sum(exp(L) * W))
			SE <- sqrt(sum((Theta[ ,1] - thetas)^2 * exp(L) * W / sum(exp(L) * W)))	
			scores[i] <- thetas
			SEscores[i] <- SE
		}
		if(method == "MAP"){ 
			for (i in 1:length(scores)) {       
				Theta <- scores[i]	  
				thetas <- nlm(MAP.bfactor,Theta,a=a,d=d,guess=g,upper=u,
					patdata=tabdata[i, ],logicalfact=logicalfact,itemloc=itemloc)$estimate 
				scores[i] <- thetas
			}  
		}
		if(method == "ML"){			
		    tmp <- tabdata[,itemloc[-length(itemloc)]]			 
			tmp2 <- tabdata[,itemloc[-1] - 1]
			scores[rowSums(tmp) == J, ] <- Inf
			scores[rowSums(tmp2) == J,] <- -Inf
		    for (i in 1:length(scores)) { 
		        if(any(scores[i] == -Inf | scores[i] == Inf)) next
		        Theta <- scores[i]	  
		        thetas <- nlm(MAP.bfactor,Theta,a=a,d=d,guess=g,upper=u,
		                      patdata=tabdata[i, ],logicalfact=logicalfact,itemloc=itemloc,
                              ML=TRUE)$estimate 
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
			tabdata2 <- object@tabdata[,-(ncol(fulldata)+1)]	
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tabdata2[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, 1] <- scores[j, 1]
			} 
			colnames(scoremat) <- 'G'	
			return(cbind(fulldata,scoremat))
		} else {  				
		    if(verbose) cat("\nMethod: ", method,"\n\n")			
			return(cbind(object@tabdata,scores))			
		}   
	}  
)

#------------------------------------------------------------------------------
#' @rdname fscores-methods  
setMethod(
	f = "fscores",
	signature = 'polymirtClass',
	definition = function(object, rotate = '', full.scores = FALSE, ndraws = 3000, thin = 5)
	{ 	        
		cand.t.var <- 1
		theta0 <- object@Theta
		K <- object@K
		nfact <- ncol(theta0)		
		so <- summary(object, rotate = rotate, print = FALSE)
		lambdas <- rotateLambdas(so)
		zetas <- object@pars$zetas
		guess <- object@guess	   
		guess[is.na(guess)] <- 0
	    upper <- object@upper       
	    upper[is.na(upper)] <- 1
		data <- cbind(object@data,object@fulldata)
		Names <- c(colnames(object@data[,1:length(K)]),colnames(object@F),
			paste("SE_F",1:nfact,sep=''))
		tabdata <- unique(data)[ ,-c(1:length(K))]			
		itemloc <- object@itemloc            		
		Theta <- list()
		for(i in 1:nfact)
			Theta[[i]] <- matrix(0,ncol=ndraws/thin,nrow=nrow(tabdata))		
		theta0 <- matrix(0,nrow(tabdata),nfact)        
		for(i in 1:30){			
			theta0 <- draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, 
			                      upper=upper, fulldata=tabdata, K=K, itemloc=itemloc, 
			                      cand.t.var=cand.t.var)
			if(attr(theta0,'Proportion Accepted') > .4) cand.t.var <- cand.t.var + .2
			if(attr(theta0,'Proportion Accepted') < .3) cand.t.var <- cand.t.var - .2
		}
		ind <- 1
		for(i in 1:ndraws){			
		    theta0 <- draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, 
		                upper=upper, fulldata=tabdata, K=K, itemloc=itemloc, 
		                cand.t.var=cand.t.var)			
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

#------------------------------------------------------------------------------
#' @rdname fscores-methods  	
setMethod(
	f = "fscores",
	signature = 'confmirtClass',
	definition = function(object, full.scores = FALSE, ndraws = 3000, thin = 5)
	{ 	
		cand.t.var <- 1
		estComp <- object@estComp
		sig <- object@pars$sig
		mu <- object@pars$mu
		theta0 <- object@Theta
		K <- object@K
		nfact <- ncol(theta0)
		nfactNames <- ncol(object@F)
		factorNames <- colnames(object@F)
		lambdas <- object@pars$lambdas		
		zetas <- object@pars$zetas
		guess <- object@guess
		guess[is.na(guess)] <- 0
		upper <- object@upper       
		upper[is.na(upper)] <- 1
		data <- cbind(object@data,object@fulldata)		
		Names <- c(colnames(object@data[,1:length(K)]),paste('F_',1:nfact,sep=''),
			paste("SE_",1:nfact,sep=''))
		tabdata <- unique(data)[,-c(1:length(K))]			
		fulldata <- object@fulldata
		itemloc <- object@itemloc
		Theta <- list()
		prodlist <- object@prodlist
		if(length(prodlist) == 0) prodlist <- NULL	
		for(i in 1:nfact)
			Theta[[i]] <- matrix(0,ncol=ndraws/thin,nrow=nrow(tabdata))		
		theta0 <- matrix(0,nrow(tabdata),nfact)		
		for(i in 1:30){
		    theta0 <- draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, 
		                upper=upper, fulldata=tabdata, K=K, itemloc=itemloc, 
		                cand.t.var=cand.t.var, prior.t.var=sig, prior.mu=mu, estComp=estComp, 
		                prodlist=prodlist)			
			if(attr(theta0,'Proportion Accepted') > .4) cand.t.var <- cand.t.var + .2
			if(attr(theta0,'Proportion Accepted') < .3) cand.t.var <- cand.t.var - .2
		}
		ind <- 1
		for(i in 1:ndraws){			
		    theta0 <- draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, 
		                upper=upper, fulldata=tabdata, K=K, itemloc=itemloc, 
		                cand.t.var=cand.t.var, prior.t.var=sig, prior.mu=mu, estComp=estComp, 
		                prodlist=prodlist)
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
			colnames(scoremat) <- factorNames
			tmp <- unique(data)[,1:length(K)]
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tmp[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, ] <- expscores[j, ]
			}              
			return(cbind(object@data,scoremat))
		}	
	}	
)

