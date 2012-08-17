#' Monte Carlo Log-Likelihood Calculation
#' 
#' Calculates a new object that contain the Monte Carlo estimated observed
#' log-likelihood values for \code{polymirt} and \code{confmirt} objects
#' 
#' @name calcLogLik
#' @usage 
#' calcLogLik(object, ...)
#'
#' \S4method{calcLogLik}{confmirtClass}(object,
#'    draws = 2000, G2 = TRUE)
#' @aliases calcLogLik-method calcLogLik,confmirtClass-method
#' @param object a model of class \code{confmirtClass}
#' @param draws the number of Monte Carlo draws
#' @param G2 logical; estimate the G2 model fit statistic?
#' @param ... parameters that are passed
#' @section Methods: \describe{ \item{calcLogLik}{\code{signature(object = "confmirtClass")}} }
#' @return Returns an object of class \code{confmirtClass} with the log-likelihood, standard error,
#' and (possibly) the G^2 model fit statistic if there is no missing data.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @docType methods
#' @rdname calcLogLik-methods  
#' @seealso
#' \code{\link{confmirt}}
#' @keywords calcLogLik
#' @examples
#' 
#' \dontrun{
#' 
#' mod1withLogLik <- calcLogLik(mod1, draws = 5000)
#' 
#'   }
#'
setMethod(
	f = "calcLogLik",
	signature = signature(object = 'confmirtClass'),
	definition = function(object, draws = 2000, G2 = TRUE)
	{	        
        pars <- object@pars
	    tol <- 1e-8	    
        fulldata <- object@fulldata
        prodlist <- object@prodlist
        itemloc <- object@itemloc
        N <- nrow(fulldata)
	    J <- length(pars)-1
	    nfact <- length(ExtractLambdas(pars[[1]])) - length(object@prodlist)	    
        LL <- matrix(0, N, draws)
        grp <- ExtractGroupPars(pars[[length(pars)]])
        for(draw in 1:draws){
	        if(nfact > 1)		
	            theta <-  mvtnorm::rmvnorm(N,grp$gmeans, grp$gcov) 
	        else
	            theta <- as.matrix(rnorm(N,grp$gmeans, grp$gcov))
	        if(length(prodlist) > 0)
	            theta <- prodterms(theta,prodlist)	        	    	
	        itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)    
	        for (i in 1:J)
	            itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=theta)	            	        
	        tmp <- itemtrace*fulldata	        
	        tmp[tmp < tol] <- 1    
	        LL[ ,draw] <- exp(rowSums(log(tmp)))
        }
        LL[is.nan(LL)] <- 0 
        rwmeans <- rowMeans(LL) 
        logLik <- sum(log(rwmeans))
		data <- object@data
		pats <- apply(data,1,paste,collapse = "/")			
		freqs <- table(pats)
		nfreqs <- length(freqs)		
		r <- as.vector(freqs)
		ncolfull <- ncol(data)
		tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
		tabdata <- suppressWarnings(matrix(as.numeric(tabdata),nfreqs,ncolfull,TRUE))
		tabdata <- cbind(tabdata,r)	
		expected <- rep(0,nrow(tabdata))
		for (j in 1:nrow(tabdata)){          
			TFvec <- colSums(ifelse(t(data) == tabdata[j,1:ncolfull],1,0)) == ncolfull 
			TFvec[is.na(TFvec)] <- FALSE	
			expected[j] <- mean(rwmeans[TFvec])			
			rwmeans[TFvec] <- rwmeans[TFvec]/r[j]
		}
		expected[is.nan(expected)] <- NA
		tabdata <- cbind(tabdata,expected*N)
		object@tabdata <- tabdata		
		logN <- 0
		logr <- rep(0,length(r))
		for (i in 1:N) logN <- logN + log(i)
		for (i in 1:length(r)) 
			for (j in 1:r[i]) 
				logr[i] <- logr[i] + log(j)    		
		if(sum(logr) != 0)		
			logLik <- logLik + logN/sum(logr)			
		SElogLik <- sqrt(var(log(rowMeans(LL))) / draws)		
        nestpars <- nconstr <- 0
        for(i in 1:length(pars))
            nestpars <- nestpars + sum(pars[[i]]@est)
        if(length(object@constrain) > 0)
            for(i in 1:length(object@constrain))
                nconstr <- nconstr + length(object@constrain[[i]]) - 1 
        nfact <- object@nfact
        df <- length(r) - nestpars + nconstr + nfact*(nfact - 1)/2 - 1 #-nmissingtabdata					
		AIC <- (-2) * logLik + 2 * (length(r) - df - 1)
		BIC <- (-2) * logLik + (length(r) - df - 1)*log(N)				
		if(G2){						
			if(any(is.na(data))){
			    object@G2 <- object@p <- object@RMSEA <- object@TLI <- NaN
			} else {				
				G2 <- 2 * sum(log(1/(N*rwmeans)))
				p <- 1 - pchisq(G2,df) 
				object@G2 <- G2	
				object@p <- p				
				object@RMSEA <- ifelse((G2 - df) > 0, 
				    sqrt(G2 - df) / sqrt(df * (N-1)), 0)
				null.mod <- object@null.mod
				object@TLI <- (null.mod@X2 / null.mod@df - G2/df) / (null.mod@X2 / null.mod@df - 1)
			}	            
		}        
		object@tabdata <- tabdata
		object@logLik <- logLik
		object@SElogLik <- SElogLik		
		object@AIC <- AIC
		object@BIC <- BIC
		object@df <- as.integer(df)
		return(object)
	} 	
)
