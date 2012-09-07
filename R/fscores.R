#' Methods for Function fscores
#' 
#' Computes MAP, EAP, or ML factor scores with a multivariate normal prior distribution. 
#'
#' 
#' @usage 
#' fscores(object, ...)
#' 
#' @aliases fscores-method fscores,ExploratoryClass-method
#' fscores,ConfirmatoryClass-method fscores,MultipleGroupClass-method
#' @docType methods
#' @section Methods: 
#' \describe{ \item{fscores}{\code{signature(object =
#' "ExploratoryClass")}} \item{fscores}{\code{signature(object = "ConfirmatoryClass")}}
#' \item{fscores}{\code{signature(object = "MultipleGroupClass")}}}
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
#' @param verbose logical; print verbose output messages?
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
        prodlist <- attr(pars, 'prodlist')
        nfact <- object@nfact - length(prodlist)
        nLambdas <- object@nfact
        itemloc <- object@itemloc
        gp <- ExtractGroupPars(object@pars[[length(itemloc)]])
        if(!pars[[1]]@bfactor && rotate != 'CONFIRMATORY'){
            so <- summary(object, rotate = rotate, verbose = FALSE)            
            a <- rotateLambdas(so)
            gp$gmeans <- rep(0, nfact)
            gp$gcov <- so$fcor
        } else {
            a <- matrix(0, J, nLambdas)
            for(i in 1:J){
                a[i, ] <- ExtractLambdas(pars[[i]])            
                pars[[i]]@bfactor <- FALSE
            }
        }        			
        if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))
		theta <- as.matrix(seq(-4,4,length.out = quadpts))
		ThetaShort <- Theta <- thetaComb(theta,nfact)         
        if(length(prodlist) > 0)        
            Theta <- prodterms(Theta,prodlist)
		fulldata <- object@data 
		tabdata <- object@tabdatalong
		tabdata <- tabdata[ ,-ncol(tabdata)]
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)			                
		W <- mvtnorm::dmvnorm(ThetaShort,gp$gmeans,gp$gcov) 
		W <- W/sum(W)                
		itemtrace <- matrix(0, ncol=ncol(tabdata), nrow=nrow(Theta))        
        for (i in 1:J)
            itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)                    
		for (i in 1:nrow(tabdata)){				
			L <- rowSums(log(itemtrace)[ ,as.logical(tabdata[i,])])			
			thetas <- colSums(ThetaShort * exp(L) * W / sum(exp(L) * W))
			SE <- sqrt(colSums(t((t(ThetaShort) - thetas))^2 * exp(L) * W / sum(exp(L) * W)))	
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
			tmp2 <- tabdata[,itemloc[-1] - 1]			             
			scores[rowSums(tmp2) == J,] <- scores[rowSums(tmp2) == 0,] <- NA
			SEscores[is.na(scores[,1]), ] <- rep(NA, nfact)
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
		colnames(scores) <- paste('F', 1:ncol(scores), sep='')          
		if (full.scores){               
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=ncol(scores))
			tabdata2 <- object@tabdata
            tabdata2 <- tabdata2[,-ncol(tabdata2)]
            fulldata2 <- fulldata
            fulldata2[is.na(fulldata2)] <- tabdata2[is.na(tabdata2)] <- 999
			for (j in 1:nrow(tabdata2)){          
				TFvec <- colSums(ifelse(t(fulldata2) == tabdata2[j, ],1,0)) == ncol(fulldata2)
                tmp <- matrix(rep(scores[j, ], sum(TFvec)), nrow=sum(TFvec), byrow=TRUE)
                scoremat[TFvec, ] <- tmp
			} 
			colnames(scoremat) <- colnames(scores)	
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
        ret <- fscores(object, rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, quadpts=quadpts, 
                       verbose=verbose)
        return(ret)
	}	
)

#------------------------------------------------------------------------------
#' @rdname fscores-methods      
setMethod(
    f = "fscores",
    signature = 'MultipleGroupClass',
    definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", 
                          quadpts = NULL, verbose = TRUE)
    { 	        
        cmods <- object@cmods
        ngroups <- length(cmods)
        for(g in 1:ngroups)
            class(cmods[[g]]) <- 'ConfirmatoryClass'
        ret <- vector('list', length(cmods))
        for(g in 1:ngroups)
            ret[[g]] <- fscores(cmods[[g]], rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, 
                           quadpts=quadpts, verbose=verbose)
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
