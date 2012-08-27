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
        nfact <- object@nfact
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
        gp <- ExtractGroupPars(object@pars[[length(itemloc)]])        
		W <- mvtnorm::dmvnorm(Theta,gp$gmeans,gp$gcov) 
		W <- W/sum(W)
        prodlist <- attr(pars, 'prodlist')
        ThetaShort <- Theta
        if(length(prodlist) > 0)        
            Theta <- prodterms(Theta,prodlist)
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
                estimate <- nlm(MAP.mirt,tmp,pars=pars, patdata=tabdata[i, ],
                                itemloc=itemloc, gp=gp, prodlist=prodlist, hessian=TRUE)				
				scores[i, ] <- estimate$estimate
				SEest <- try(sqrt(diag(solve(estimate$hessian))))
				if(is(SEest, 'try-error')) SEest <- rep(NA, nfact)
				SEscores[i, ] <- SEest
			}  
		}
		if(method == "ML"){
			tmp <- tabdata[,itemloc[-length(itemloc)]]			 
			tmp2 <- tabdata[,itemloc[-1] - 1]			 
            scores[rowSums(tmp) == J, ] <- NA
			scores[rowSums(tmp2) == J,] <- NA
			SEscores[rowSums(tmp) == J, ] <- SEscores[rowSums(tmp2) == J, ] <- rep(NA, nfact)
			for (i in 1:nrow(scores)){
				if(any(is.na(scores[i, ]))) next 
				Theta <- scores[i, ]	  
				estimate <- nlm(MAP.mirt,Theta,pars=pars,patdata=tabdata[i, ],
				    itemloc=itemloc, gp=gp, prodlist=prodlist, ML=TRUE, hessian = TRUE)
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
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", 
	                      quadpts = NULL, verbose = TRUE)
	{ 	        
        class(object) <- 'ExploratoryClass'
        ret <- fscores(object, rotate = 'none', full.scores=full.scores, method=method, quadpts=quadpts, 
                       verbose=verbose)
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
