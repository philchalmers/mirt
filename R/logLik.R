#' Monte Carlo Log-Likelihood Calculation
#' 
#' Calculates a new object that contain the Monte Carlo estimated observed
#' log-likelihood values for \code{polymirt} and \code{confmirt} objects
#' 
#' @name logLik
#' @usage 
#' logLik(object, ...)
#'
#' \S4method{logLik}{polymirtClass}(object,
#'    draws = 2000, G2 = TRUE)
#'
#' \S4method{logLik}{confmirtClass}(object,
#'    draws = 2000, G2 = TRUE)
#' @aliases logLik-method logLik,polymirtClass-method 
#' logLik,confmirtClass-method
#' @param object a model of class \code{mirtClass} or \code{bfactorClass}
#' @param draws the number of Monte Carlo draws
#' @param G2 logical; estimate the G2 model fit statistic?
#' @param ... parameters that are passed
#' @section Methods: \describe{ \item{logLik}{\code{signature(object =
#' "polymirtClass")}} \item{logLik}{\code{signature(object = "confmirtClass")}} }
#' @return Returns an object of class \code{polymirtClass} or
#' \code{confmirtClass} with the log-likelihood, standard error, and (possibly)
#' the G^2 model fit statistic.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @docType methods
#' @rdname logLik-methods  
#' @seealso
#' \code{\link{polymirt}}, \code{\link{confmirt}}
#' @keywords logLik
#' @examples
#' 
#' \dontrun{
#' 
#' mod1withLogLik <- logLik(mod1, draws = 5000)
#' 
#'   }
#'
setMethod(
	f = "logLik",
	signature = signature(object = 'confmirtClass'),
	definition = function(object, draws = 2000, G2 = TRUE)
	{	
		nfact <- ncol(object@Theta)
		nfactNames <- ifelse(length(object@prodlist) > 0, 
			length(object@prodlist) + nfact, nfact)		
		N <- nrow(object@Theta)
		J <- length(object@K)
		pars <- object@pars
		lambdas <- pars[ ,1:nfactNames]
		lambdas[is.na(lambdas)] <- 0
		zetas <- pars[,(nfactNames+1):ncol(pars)]
		zetas <- t(zetas)[!is.na(t(zetas))]		
		mu <- object@gpars$u
		sigma <- object@gpars$sig		
		LL <- matrix(0,N,draws)		
		guess <- object@guess
		guess[is.na(guess)] <- 0
		K <- object@K	
		fulldata <- object@fulldata	
		for(i in 1:draws){
			theta <- rmvnorm(N,mu,sigma)	
			if(nfact < nfactNames) 
				theta <- prodterms(theta, object@prodlist)	
			LL[,i] <- .Call('logLik', 					
						as.numeric(lambdas),
						as.numeric(zetas),
						as.numeric(guess),
						as.numeric(theta),
						as.integer(fulldata),
						as.integer(object@itemloc-1),
						as.integer(object@K),
						as.integer(J),
						as.integer(N),
						as.integer(nfactNames),
						as.integer(object@estComp))		
		}		
		rwmeans <- rowMeans(LL) 
		logLik <- sum(log(rwmeans))				
		pats <- apply(fulldata,1,paste,collapse = "/")
		freqs <- table(pats)
		nfreqs <- length(freqs)		
		r <- as.vector(freqs)
		ncolfull <- ncol(fulldata)
		tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
		tabdata <- matrix(as.numeric(tabdata),nfreqs,ncolfull,TRUE)
		tabdata <- cbind(tabdata,r)		 		 				
		pats <- apply(fulldata,1,paste,collapse = "/")
		freqs <- table(pats)			
		r <- as.vector(freqs)
		logN <- 0
		logr <- rep(0,length(r))
		for (i in 1:N) logN <- logN + log(i)
		for (i in 1:length(r)) 
			for (j in 1:r[i]) 
				logr[i] <- logr[i] + log(j)    		
		if(sum(logr) != 0)		
			logLik <- logLik + logN/sum(logr)			
		SElogLik <- sqrt(var(log(rowMeans(LL))) / draws)
		x <- object@estpars	
		df <- as.integer(length(r) - sum(x$estlam) - sum(x$estgcov) - 
			sum(x$estgmeans) - length(zetas) + object@nconstvalues + 
			nfact*(nfact - 1)/2 - sum(x$estGuess) - 1)			
		AIC <- (-2) * logLik + 2 * (length(r) - df - 1)
		BIC <- (-2) * logLik + (length(r) - df - 1)*log(N)
		if(G2){			
			data <- object@data
			if(any(is.na(data))){
				object@G2 <- 0	
				object@p <- 1					
			} else {			
				pats <- apply(data,1,paste,collapse = "/")			
				freqs <- table(pats)
				nfreqs <- length(freqs)		
				r <- as.vector(freqs)
				ncolfull <- ncol(data)
				tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
				tabdata <- matrix(as.numeric(tabdata),nfreqs,ncolfull,TRUE)
				tabdata <- cbind(tabdata,r)	
				expected <- rep(0,nrow(tabdata))	
				for (j in 1:nrow(tabdata)){          
					TFvec <- colSums(ifelse(t(data) == tabdata[j,1:ncolfull],1,0)) == ncolfull        
					expected[j] <- mean(rwmeans[TFvec])
					rwmeans[TFvec] <- rwmeans[TFvec]/r[j]
				}
				tabdata <- cbind(tabdata,expected*N)
				G2 <- 2 * sum(log(1/(N*rwmeans)))
				p <- 1 - pchisq(G2,df) 
				object@G2 <- G2	
				object@p <- p
				object@tabdata <- tabdata
			}	
		}	
		object@logLik <- logLik
		object@SElogLik <- SElogLik		
		object@AIC <- AIC
		object@BIC <- BIC
		object@df <- df		
		return(object)
	} 	
)

# @rdname logLik-methods  
setMethod(
	f = "logLik",
	signature = signature(object = 'polymirtClass'),
	definition = function(object, draws = 2000, G2 = TRUE)
	{	
		nfact <- ncol(object@Theta)
		N <- nrow(object@Theta)
		J <- length(object@K)
		pars <- object@pars
		lambdas <- pars[,1:nfact]
		zetas <- pars[,(nfact+1):ncol(pars)]
		zetas <- t(zetas)[!is.na(t(zetas))]		
		mu <- rep(0,nfact)
		sigma <- diag(nfact)		
		LL <- matrix(0,N,draws)		
		guess <- object@guess
		guess[is.na(guess)] <- 0
		K <- object@K		
		fulldata <- object@fulldata
		estComp <- rep(FALSE,J)
		for(i in 1:draws){
			theta <- rmvnorm(N,mu,sigma)				
			LL[,i] <- .Call('logLik', 					
						as.numeric(lambdas),
						as.numeric(zetas),
						as.numeric(guess),
						as.numeric(theta),
						as.integer(fulldata),
						as.integer(object@itemloc-1),
						as.integer(object@K),
						as.integer(J),
						as.integer(N),
						as.integer(nfact),
						as.integer(estComp))		
		}		
		rwmeans <- rowMeans(LL)
		logLik <- sum(log(rwmeans))		
		pats <- apply(fulldata,1,paste,collapse = "/")
		freqs <- table(pats)
		nfreqs <- length(freqs)		
		r <- as.vector(freqs)
		ncolfull <- ncol(fulldata)
		tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
		tabdata <- matrix(as.numeric(tabdata),nfreqs,ncolfull,TRUE)
		tabdata <- cbind(tabdata,r)		
		logN <- 0
		logr <- rep(0,length(r))
		for (i in 1:N) logN <- logN + log(i)
		for (i in 1:length(r)) 
			for (j in 1:r[i]) 
				logr[i] <- logr[i] + log(j) 
		if(sum(logr) != 0)								
			logLik <- logLik + logN/sum(logr)		
		SElogLik <- sqrt(var(log(rwmeans)) / draws)
		df <- (length(r) - 1) - nfact*J - sum(K - 1) + nfact*(nfact - 1)/2 - sum(object@estGuess)
		AIC <- (-2) * logLik + 2 * (length(r) - df - 1)
		BIC <- (-2) * logLik + (length(r) - df - 1)*log(N)
		if(G2){				
			data <- object@data
			if(any(is.na(data))){
				object@G2 <- 0	
				object@p <- 1					
			} else {
				pats <- apply(data,1,paste,collapse = "/")			
				freqs <- table(pats)
				nfreqs <- length(freqs)		
				r <- as.vector(freqs)
				ncolfull <- ncol(data)
				tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
				tabdata <- matrix(as.numeric(tabdata),nfreqs,ncolfull,TRUE)
				tabdata <- cbind(tabdata,r)	
				expected <- rep(0,nrow(tabdata))	
				for (j in 1:nrow(tabdata)){          
					TFvec <- colSums(ifelse(t(data) == tabdata[j,1:ncolfull],1,0)) == ncolfull        
					expected[j] <- mean(rwmeans[TFvec])
					rwmeans[TFvec] <- rwmeans[TFvec]/r[j]					
				}
				tabdata <- cbind(tabdata,expected*N)
				G2 <- 2 * sum(log(1/(N*rwmeans)))
				p <- 1 - pchisq(G2,df) 
				object@G2 <- G2	
				object@p <- p
				object@tabdata <- tabdata
			}	
		}		
		object@logLik <- logLik
		object@SElogLik <- SElogLik		
		object@AIC <- AIC
		object@BIC <- BIC
		object@df <- as.integer(df)
		return(object)
	} 	
)

