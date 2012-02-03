#' Methods for Function fscores
#' 
#' Computes MAP or EAP factor scores for \code{mirt} and \code{bfactor} models,
#' or stochastic approximations for \code{polymirt} and \code{confmirt}. Note
#' that only the general factor scores are computed for bifactor models.
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

#' @name fscores
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
#' @keywords factor.scores
#' @rdname fscores-methods   
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

#' @rdname fscores-methods
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


