#' Monte Carlo Log-Likelihood Calculation
#' 
#' Calculates a new object that contain the Monte Carlo estimated observed
#' log-likelihood values for mirt objects estimated with the MH-RM algorithm
#' 
#' @name calcLogLik
#' @usage 
#' calcLogLik(object, ...)
#'
#' \S4method{calcLogLik}{ExploratoryClass}(object,
#'    draws = 3000, G2 = TRUE, cl = NULL)
#' \S4method{calcLogLik}{ConfirmatoryClass}(object,
#'    draws = 3000, G2 = TRUE, cl = NULL)
#' @aliases calcLogLik-method calcLogLik,ExploratoryClass-method
#' calcLogLik,ConfirmatoryClass-method
#' @param object a model of class \code{ConfirmatoryClass} or \code{ExploratoryClass}
#' @param draws the number of Monte Carlo draws
#' @param G2 logical; estimate the G2 model fit statistic?
#' @param cl a cluster object from the \code{parallel} package (set from using \code{makeCluster(ncores)})
#' @param ... parameters that are passed
#' @section Methods: 
#' \describe{ \item{calcLogLik}{\code{signature(object = "ConfirmatoryClass")}, 
#' \code{signature(object = "ExploratoryClass")} }
#' } 
#' @return Returns an object with the log-likelihood, standard errors, information matrix, 
#' and (possibly) the G^2 and other model fit statistic if there is no missing data.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @docType methods
#' @rdname calcLogLik-methods  
#' @export calcLogLik
#' @seealso
#' \code{\link{confmirt}}, \code{\link{multipleGroup}}
#' @keywords calcLogLik
#' @examples
#' 
#' \dontrun{
#' 
#' # no parallel
#' mod1withLogLik <- calcLogLik(mod1, draws=5000)
#' 
#' #with parallel using detected number of cores
#' library(parallel)
#' cl <- makeCluster(detectCores())
#' mod1withLogLik <- calcLogLik(mod1, draws=5000, cl=cl)
#' 
#'   }
#'
setMethod(
	f = "calcLogLik",
	signature = signature(object = 'ExploratoryClass'),
	definition = function(object, draws = 3000, G2 = TRUE, cl = NULL)
	{
        LLdraws <- function(LLDUMMY=NULL, nfact, N, grp, prodlist, fulldata, object, J){
            if(nfact > 1) theta <-  mvtnorm::rmvnorm(N,grp$gmeans, grp$gcov)
            else theta <- as.matrix(rnorm(N,grp$gmeans, grp$gcov))
            if(length(prodlist) > 0)
                theta <- prodterms(theta,prodlist)	        	    	
            itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)    
            if(length(object@mixedlist) > 1){ 
                colnames(theta) <- object@mixedlist$factorNames
                fixed.design.list <- designMats(covdata=object@mixedlist$covdata, 
                                                fixed=object@mixedlist$fixed, 
                                                Thetas=theta, nitems=J, 
                                                itemdesign=object@mixedlist$itemdesign, 
                                                fixed.identical=object@mixedlist$fixed.identical)	            
            }
            for (i in 1:J) itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- 
                ProbTrace(x=pars[[i]], Theta=theta, fixed.design=fixed.design.list[[i]])	            	        
            tmp <- itemtrace*fulldata	        
            tmp[tmp < tol] <- 1    
            return(exp(rowSums(log(tmp))))           
        }
        pars <- object@pars
	    tol <- .Machine$double.eps	    
        fulldata <- object@fulldata
        prodlist <- object@prodlist
        itemloc <- object@itemloc
        N <- nrow(fulldata)
	    J <- length(pars)-1
	    nfact <- length(ExtractLambdas(pars[[1]])) - length(object@prodlist) - pars[[1]]@nfixedeffects	    
        LL <- matrix(0, N, draws)
        grp <- ExtractGroupPars(pars[[length(pars)]]) 
        fixed.design.list <- vector('list', J)         
        if(!is.null(cl))            
            LL <- parallel::parApply(cl=cl, LL, MARGIN=1, FUN=LLdraws, nfact=nfact, N=N, grp=grp, prodlist=prodlist, 
                           fulldata=fulldata, object=object, J=J)
        else
            for(draw in 1:draws)
                LL[ ,draw] <- LLdraws(nfact=nfact, N=N, grp=grp, prodlist=prodlist, 
                                      fulldata=fulldata, object=object, J=J)    
        LL[is.nan(LL)] <- 0 
        rwmeans <- rowMeans(LL) 
        logLik <- sum(log(rwmeans))
        SElogLik <- sqrt(var(log(rowMeans(LL))) / draws) 
        logLikpre <- 0
        if(G2 == 'return'){
            logLikpre <- logLik                       
            G2 <- TRUE
        }
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
		}
		expected[is.nan(expected)] <- NA
		tabdata <- cbind(tabdata,expected*N)        
        object@Pl <- expected
		logN <- 0
		logr <- rep(0,length(r))
		for (i in 1:N) logN <- logN + log(i)
		for (i in 1:length(r)) 
			for (j in 1:r[i]) 
				logr[i] <- logr[i] + log(j)    		
		if(sum(logr) != 0)		
			logLik <- logLik + logN/sum(logr)							
        nestpars <- nconstr <- 0
        for(i in 1:length(pars))
            nestpars <- nestpars + sum(pars[[i]]@est)
        if(length(object@constrain) > 0)
            for(i in 1:length(object@constrain))
                nconstr <- nconstr + length(object@constrain[[i]]) - 1 
        nfact <- object@nfact - length(prodlist)        
        nmissingtabdata <- sum(is.na(rowSums(object@tabdata)))
        df <- length(r) - nestpars + nconstr + nfact*(nfact - 1)/2 - 1 - nmissingtabdata	
		AIC <- (-2) * logLik + 2 * (length(r) - df - 1)
		BIC <- (-2) * logLik + (length(r) - df - 1)*log(N)        
		if(G2){						
			if(any(is.na(data))){
			    object@G2 <- object@p <- object@RMSEA <- object@TLI <- NaN
			} else {
				G2 <- 2 * sum(r*log(r/(N*expected)))
				p <- 1 - pchisq(G2,df) 
				object@G2 <- G2	
				object@p <- p				
				object@RMSEA <- ifelse((G2 - df) > 0, 
				    sqrt(G2 - df) / sqrt(df * (N-1)), 0)
                if(logLikpre == 0){
    				null.mod <- object@null.mod
    				object@TLI <- (null.mod@G2 / null.mod@df - G2/df) / (null.mod@G2 / null.mod@df - 1)
                }
			}	            
		}                
		object@logLik <- logLik
        if(logLikpre < 0)
            object@logLik <- logLikpre
		object@SElogLik <- SElogLik		
		object@AIC <- AIC
		object@BIC <- BIC
		object@df <- as.integer(df)
		return(object)
	} 	
)

setMethod(
    f = "calcLogLik",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, draws = 2000, G2 = TRUE, cl = NULL)
    {	        
        class(object) <- 'ExploratoryClass'
        ret <- calcLogLik(object, draws=draws, G2=G2, cl=cl)
        class(ret) <- 'ConfirmatoryClass'
        return(ret)
    } 	
)

