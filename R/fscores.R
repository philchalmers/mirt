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
#' @aliases fscores-method fscores,ExploratoryClass-method
#' fscores,ConfirmatoryClass-method
#' @docType methods
#' @section Methods: 
#' \describe{ \item{fscores}{\code{signature(object =
#' "ExploratoryClass")}} \item{fscores}{\code{signature(object = "ConfirmatoryClass")}}}
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
#' @param object a model of class \code{ExploratoryClass} or \code{ConfirmatoryClass}
#' @param full.scores if \code{FALSE} (default) then a summary table with
#' factor scores for each unique pattern is displayed. Otherwise the original
#' data matrix is returned with the computed factor scores appended to the
#' rightmost column
#' @param rotate rotation declaration to be used when estimating the factor scores. If \code{""} then the 
#' \code{object@@rotate} default value is used
#' @param method type of factor score estimation method. Can be expected
#' a-posteriori (\code{"EAP"}), Bayes modal (\code{"MAP"}), or maximum likelihood 
#' (\code{"ML"})
#' @param quadpts number of quadratures to use per dimension
#' @param ndraws number of MH samples to draw for each response pattern 
#' @param thin controls how much the chain should be thinned by, default
#' collects every 5th draw (\code{thin = 5}). Note that \code{ndraws/thin} must be a whole number.
#' for \code{confmirtClass} objects only
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
#' mod2 <- confmirt(Science, 1)
#' tabscores2 <- fscores(mod2, ndraws = 5000)
#' 
#'   }
#' 
setMethod(
	f = "fscores",
	signature = 'ExploratoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", 
                          quadpts = NULL, verbose = TRUE)
	{           
        pars <- object@pars        
		K <- object@K        
        J <- length(K)
        nfact <- pars[[1]]@nfact
        if(!pars[[1]]@bfactor){
            so <- summary(object, rotate = rotate, print = FALSE)
            a <- rotateLambdas(so)		
        } else {
            a <- matrix(0, J, nfact)
            for(i in 1:J){
                a[i, ] <- ExtractLambdas(pars[[i]])            
                pars[[i]]@bfactor <- FALSE
            }
        }
		itemloc <- object@itemloc	
        if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))
		theta <- as.matrix(seq(-4,4,length.out = quadpts))
		Theta <- thetaComb(theta,nfact)
		fulldata <- object@data 
		tabdata <- object@tabdatalong
		tabdata <- tabdata[ ,-ncol(tabdata)]
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)			
		W <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact)) 
		W <- W/sum(W)
		itemtrace <- matrix(0, ncol=ncol(tabdata), nrow=nrow(Theta))        
        for (i in 1:J)
            itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)                    
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
                estimate <- nlm(MAP.mirt,tmp,pars=pars,patdata=tabdata[i, ],
                                itemloc=itemloc, hessian=TRUE)				
				scores[i, ] <- estimate$estimate
				SEest <- try(sqrt(diag(solve(estimate$hessian))))
				if(is(SEest, 'try-error')) SEest <- rep(NA, nfact)
				SEscores[i, ] <- SEest
			}  
		}
		if(method == "ML"){
			tmp <- tabdata[,itemloc[-length(itemloc)]]			 
			tmp2 <- tabdata[,itemloc[-1] - 1]			 
            scores[rowSums(tmp) == J, ] <- Inf
			scores[rowSums(tmp2) == J,] <- -Inf
			SEscores[rowSums(tmp) == J, ] <- SEscores[rowSums(tmp2) == J, ] <- rep(NA, nfact)
			for (i in 1:nrow(scores)){
				if(any((scores[i, ]) == -Inf | scores[i, ] == Inf)) next 
				Theta <- scores[i, ]	  
				estimate <- nlm(MAP.mirt,Theta,pars=pars,patdata=tabdata[i, ],
				    itemloc=itemloc, ML=TRUE, hessian = TRUE)
				scores[i, ] <- estimate$estimate                
                SEest <- try(sqrt(diag(solve(estimate$hessian))))
                if(is(SEest, 'try-error')) SEest <- rep(NA, nfact)
				SEscores[i, ] <- SEest
			}  			
		}
		colnames(scores) <- colnames(object@F)
		if (full.scores){      
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=ncol(scores))
			tabdata2 <- object@tabdata[,-(ncol(fulldata)+1)]	
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tabdata2[j, ],1,0)) == ncol(fulldata)
                tmp <- matrix(rep(scores[j, ], sum(TFvec)), nrow=sum(TFvec), byrow=TRUE)
                scoremat[TFvec, ] <- tmp
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
	signature = 'ConfirmatoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, ndraws = 3000, thin = 5)
	{ 	        
        cand.t.var <- 1    	
        pars <- object@pars
        K <- object@K
        nfact <- pars[[length(pars)]]@nfact
        nfactNames <- ncol(object@F)
        factorNames <- colnames(object@F)		
        data <- cbind(object@data,object@fulldata)		
        Names <- c(colnames(object@data[,1:length(K)]),paste('F_',1:nfact,sep=''),
                   paste("SE_",1:nfact,sep=''))
        tabdata <- unique(data)[,-c(1:length(K))]			
        fulldata <- object@fulldata
        itemloc <- object@itemloc
        Theta <- list()
        for(i in 1:nfact)
            Theta[[i]] <- matrix(0,ncol=ndraws/thin,nrow=nrow(tabdata))
        structgrouppars <- ExtractGroupPars(pars[[length(pars)]])
        prodlist <- attr(pars, 'prodlist')
        theta0 <- matrix(0, nrow(tabdata), nfact)        
        for(i in 1:30){
            theta0 <- draw.thetas(theta0=theta0, pars=pars, fulldata=tabdata, itemloc=itemloc, 
                                  cand.t.var=cand.t.var, prior.t.var=structgrouppars$gcov, 
                                  prior.mu=structgrouppars$gmeans, prodlist=prodlist, debug=FALSE)			
            if(attr(theta0,'Proportion Accepted') > .4) cand.t.var <- cand.t.var + .2
            if(attr(theta0,'Proportion Accepted') < .3) cand.t.var <- cand.t.var - .2
        }
        ind <- 1
        for(i in 1:ndraws){			
            theta0 <- draw.thetas(theta0=theta0, pars=pars, fulldata=tabdata, itemloc=itemloc, 
                                  cand.t.var=cand.t.var, prior.t.var=structgrouppars$gcov, 
                                  prior.mu=structgrouppars$gmeans, prodlist=prodlist, debug=FALSE)
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
			    tmp2 <- matrix(rep(expscores[j, ], sum(TFvec)), nrow=sum(TFvec), byrow=TRUE)
			    scoremat[TFvec, ] <- tmp2
			}              
			return(cbind(object@data,scoremat))
		}	
	}	
)

# MAP scoring for mirt
MAP.mirt <- function(Theta, pars, patdata, itemloc, ML=FALSE)
{    
    itemtrace <- rep(0, ncol=length(patdata))
    Theta <- matrix(Theta, nrow=1)    
    itemtrace <- matrix(0, ncol=length(patdata), nrow=nrow(Theta))        
    for (i in 1:(length(pars)))
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)		
    L <- sum(log(itemtrace)[as.logical(patdata)])
    mu <- 0
    sigma <- 1
    L <- ifelse(ML, -L, (-1)*(L + sum(log(exp(-0.5*((Theta - mu)/sigma)^2)))))
    L  
}