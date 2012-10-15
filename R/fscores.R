#' Methods for Function fscores
#' 
#' Computes MAP, EAP, or ML factor scores with a multivariate normal prior distribution.
#' 
#'
#' @aliases fscores
#' @param object a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or 
#' \code{MultipleGroupClass}
#' @param full.scores if \code{FALSE} (default) then a summary table with
#' factor scores for each unique pattern is displayed. Otherwise the original
#' data matrix is returned with the computed factor scores appended to the
#' rightmost column
#' @param rotate rotation declaration to be used when estimating the factor scores. If \code{""} then the 
#' \code{object@@rotate} default value is used (only applicable to \code{ExploratoryClass} objects)
#' @param method type of factor score estimation method. Can be expected
#' a-posteriori (\code{"EAP"}), Bayes modal (\code{"MAP"}), or maximum likelihood 
#' (\code{"ML"})
#' @param quadpts number of quadratures to use per dimension
#' @param response.vector an optional argument used to calculate the factor scores and standard errors
#' for a given response vector that may or may not have been in the original dataset 
#' @param verbose logical; print verbose output messages?
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords factor.scores
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
#' #calculate MAP for a given response vector
#' fscores(mod, method='MAP', response.vector = c(1,2,3,4))
#'   }
#'
fscores <- function(object, rotate = '', full.scores = FALSE, method = "EAP", 
                    quadpts = NULL, response.vector = NULL, verbose = TRUE)
{
    ret <- fscores.internal(object=object, rotate=rotate, full.scores=full.scores, method=method, 
                            quadpts=quadpts, response.vector=response.vector, verbose=verbose)
    ret    
}

 
#------------------------------------------------------------------------------
setGeneric("fscores.internal", 
           def = function(object, ...) standardGeneric("fscores.internal")
)

setMethod(
	f = "fscores.internal",
	signature = 'ExploratoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", 
                          quadpts = NULL, response.vector = NULL, verbose = TRUE)
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
			tmp2 <- tabdata[,itemloc[-1] - 1, drop = FALSE]			             
			scores[rowSums(tmp2) == J,] <- NA
            tmp2 <- tabdata[,itemloc[-length(itemloc)], drop = FALSE]
            scores[rowSums(tmp2) == J,] <- NA
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
            tabdata2 <- object@tabdatalong
            tabdata2 <- tabdata2[,-ncol(tabdata2)]
            scoremat <- .Call("fullScores", object@fulldata, tabdata2, scores)			 
			colnames(scoremat) <- colnames(scores)	
			return(cbind(fulldata,scoremat))
		} else {						            
            r <- object@tabdata[,ncol(object@tabdata)]            
            T <- E <- matrix(NA, 1, ncol(scores))
            for(i in 1:nrow(scores)){
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
	                      quadpts = NULL, response.vector = NULL, verbose = TRUE)
	{ 	        
        class(object) <- 'ExploratoryClass'
        ret <- fscores(object, rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, quadpts=quadpts, 
                       response.vector=response.vector, verbose=verbose)
        return(ret)
	}	
)

#------------------------------------------------------------------------------
setMethod(
    f = "fscores.internal",
    signature = 'MultipleGroupClass',
    definition = function(object, rotate = '', full.scores = FALSE, method = "EAP", 
                          quadpts = NULL, response.vector = NULL, verbose = TRUE)
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
