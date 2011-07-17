########################################## 

setMethod(
	f = "print",
	signature = signature(x = 'bfactorClass'),
	definition = function(x, ...){  
		cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")		
		cat("Full-information bifactor analysis with ", 
		length(unique(x@specific)), " specific factors \n", sep='')
		if(x@converge == 1)	
			cat("Converged in ", x@EMiter, " iterations using ",x@quadpts,
			" quadrature points. \n", sep="")
		else 	
			cat("Estimation stopped after ", x@EMiter, " iterations using ",x@quadpts,
			" quadrature points. \n", sep="")
		cat("Log-likelihood = ", x@log.lik, "\n")
		cat("AIC = ", x@AIC, "\n")
		cat("G^2 = ", round(x@X2,2), ", df = ", 
		x@df, ", p = ", round(x@p,4), "\n")
	}
)

setMethod(
	f = "show",
	signature = signature(object = 'bfactorClass'),
	definition = function(object){  
		cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")		
		cat("Full-information bifactor analysis with ", 
		length(unique(object@specific)), " specific factors \n", sep='')
		if(object@converge == 1)	
			cat("Converged in ", object@EMiter, " iterations using ", object@quadpts,
				" quadrature points.\n", sep="")
		else 	
			cat("Estimation stopped after ", object@EMiter, " iterations using ", 
				object@quadpts,	" quadrature points.\n", sep="")
		cat("Log-likelihood = ", object@log.lik, "\n")
		cat("AIC = ", object@AIC, "\n")
		cat("G^2 = ", round(object@X2,2), ", df = ", 
		object@df, ", p = ", round(object@p,4), "\n")
	}
)

setMethod(
	f = "summary",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, digits = 3, ...){
		F <- round(object@F,digits)
		h2 <- round(object@h2,digits)
		SS <- colSums(F^2)		
		colnames(F) <- c('G',paste("F_", 1:(ncol(F)-1),sep=""))
		names(h2) <- "h2"		
		loads <- round(cbind(F,h2),digits)
		rownames(loads) <- object@itemnames  
		cat("\nFactor loadings: \n\n")
		print(loads)
		cat("\nSS loadings: ",round(SS,digits), "\n")
		cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
		if(any(h2 > 1)) 
			warning("Solution has heywood cases. Interpret with caution.")
	}
)

setMethod(
	f = "coef",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, digits = 3, ...){
		a <- as.matrix(object@pars[ ,1:(ncol(object@pars)-1)])
		d <- object@pars[ ,ncol(object@pars)]
		A <- sqrt(apply(a^2,1,sum))
		B <- -d/A 
		fac <- object@facility  
		parameters <- round(cbind(object@pars,object@guess,fac,A,B),digits)
		colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),"d", "guess", 
			"facility","mvdisc", "mvint")  
		cat("\nParameters with multivariate discrimination and intercept: \n\n")
		print(parameters)	    
		invisible(parameters)
	}
)

setMethod(
	f = "residuals",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, restype = 'LD', digits = 3, p = FALSE, ...)
	{       
		Theta <- object@Theta
		fulldata <- object@fulldata	
		N <- nrow(fulldata)	
		J <- ncol(fulldata)
		nfact <- ncol(object@F)
		lambdas <- matrix(object@pars[,1:nfact], J)
		zetas <- object@pars[,(nfact+1)]
		guess <- object@guess
		guess[is.na(guess)] <- 0
		logicalfact <- object@logicalfact
		if(restype == 'LD'){
			res <- matrix(0,J,J)
			diag(res) <- NA
			colnames(res) <- rownames(res) <- colnames(fulldata)
			prior <- dmvnorm(Theta,rep(0,2),diag(2))
			prior <- prior/sum(prior)
			for(i in 1:J){			
				for(j in 1:J){
					if(i < j){
						P1 <- P.bfactor(lambdas[i,],zetas[i], Theta, guess[i],logicalfact[i,])
						P2 <- P.bfactor(lambdas[j,],zetas[j], Theta, guess[j],logicalfact[j,])
						E22 <- N * sum(P1 * P2 * prior)
						E12 <- N * sum(P1 * (1-P2) * prior)
						E21 <- N * sum((1-P1) * P2 * prior)
						E11 <- N * sum((1-P1) * (1-P2) * prior)
						tab <- table(fulldata[,i],fulldata[,j])
						Etab <- matrix(c(E11,E12,E21,E22),2)
						s <- phi(tab) - phi(Etab)
						if(s == 0) s <- 1
						res[j,i] <- sum(((tab - Etab)^2)/Etab) * sign(s)
						res[i,j] <- sqrt( abs(res[j,i]) / N ) 
					}
				}
			}
			cat("\nLD matrix:\n\n")			
			res <- round(res,digits)	
			return(res)	
		}
		if(restype == 'exp'){
			r <- object@tabdata[ ,ncol(object@tabdata)]
			res <- round((r - object@Pl * nrow(object@fulldata)) / 
				sqrt(object@Pl * nrow(object@fulldata)),digits)
			expected <- round(object@sampsize * object@Pl/sum(object@Pl),digits)  
			tabdata <- cbind(object@tabdata,expected,res)
			colnames(tabdata) <- c(object@itemnames, "freq", "exp", "std_res")								
			return(tabdata)
		}				
	}
)

setMethod(
	f = "fitted",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, digits = 3, ...){  
		expected <- round(object@sampsize * object@Pl/sum(object@Pl),digits)  
		tabdata <- cbind(object@tabdata,expected)
		colnames(tabdata) <- c(object@itemnames, "freq", "exp")	
		print(tabdata)
		invisible(tabdata)
	}
)


########################################## 
#Main function

bfactor <- function(fulldata, specific, guess = 0, prev.cor = NULL, 
  par.prior = FALSE, startvalues = NULL, quadpts = NULL, ncycles = 300, 
  EMtol = .001, nowarn = TRUE, debug = FALSE, ...)
{ 
	#local functions
	fn <- function(pars, r1, N, guess, Theta, prior, parprior){
		a <- pars[1:(length(pars)-1)]
		d <- pars[length(pars)]
		r1 <- r1 * prior	
		N <- N * prior
		result <- .Call("loglik", 	                
						as.double(a),				
						as.double(d),
						as.double(r1),
						as.double(N),
						as.double(guess),
						as.double(as.matrix(Theta)),
						as.integer(parprior))					
	}  
	gr <- function(pars, r1, N, guess, Theta, prior, parprior){
		a <- pars[1:(length(pars)-1)]
		d <- pars[length(pars)]			
		result <- .Call("grad", 	                
						as.double(a),				
						as.double(d),
						as.double(r1),
						as.double(N),
						as.double(guess),
						as.double(as.matrix(Theta)),
						as.double(prior),
						as.integer(parprior))	    				
	} 
	
	#Main
	Call <- match.call() 
	rotate <- 'oblimin'
	itemnames <- colnames(fulldata) 
	fulldata <- as.matrix(fulldata)
	fulldata[is.na(fulldata)] <- 0
	if (!any(fulldata %in% c(0,1,NA))) stop("Data must contain only 0, 1, or NA.")
	if (length(specific) != ncol(fulldata)) 
	stop("Specific factor loadings have been declared incorrectly")  
	nfact <- length(unique(specific)) + 1  
	nitems <- ncol(fulldata)
	if(2*nfact >= nitems) stop('Model is not identified.')
	if (length(guess) == 1) guess <- rep(guess,nitems)
		else if (length(guess) > nitems || length(guess) < nitems) 
			stop("The number of guessing parameters is incorrect.")
	pats <- apply(fulldata, 1, paste, collapse = "/") 
	freqs <- table(pats)
	nfreqs <- length(freqs)
	r <- as.vector(freqs)
	sampsize <- nrow(fulldata)
	K <- rep(2,nitems)  
	tabdata <- unlist(strsplit(cbind(names(freqs)), "/"))
	tabdata <- matrix(as.numeric(tabdata), nfreqs, nitems, TRUE)
	tabdata <- cbind(tabdata,r)  
	logicalfact <- matrix(FALSE,nitems,nfact - 1)
	is.na(specific) <- FALSE
	for (i in 1:nitems) logicalfact[i,specific[i]] <- TRUE
	logicalfact <- cbind(rep(TRUE,nitems),logicalfact)     
	if (is.null(quadpts)) quadpts <- 15
	theta <- as.matrix(seq(-4,4,length.out = quadpts))
	Theta <- as.matrix(expand.grid(theta,theta))
	facility <- colMeans(fulldata)
	selvec <- 2:(nfact)    
	suppressAutoPrior <- TRUE
	if(is.logical(par.prior)) 
		if(par.prior) suppressAutoPrior <- FALSE  
			temp <- matrix(c(1,0,0),ncol = 3, nrow=nitems, byrow=TRUE)
	if(!is.logical(par.prior)){
		if(!is.null(par.prior$slope.items))
		for(i in 1:length(par.prior$slope.items))
			temp[par.prior$slope.items[i],1] <- par.prior$slope		
		if(!is.null(par.prior$int.items))
			for(i in 1:length(par.prior$int.items))
				temp[par.prior$int.items[i],2:3] <- par.prior$int		 
	}   
	par.prior <- temp  
	if (any(class(prev.cor) == c('mirt','bfactor'))) Rpoly <- prev.cor$cormat
		else if(!is.null(prev.cor)) {
			if (ncol(prev.cor) == nrow(prev.cor)) Rpoly <- prev.cor
				else stop("Correlation matrix is not square.\n")
		} else Rpoly <- cormod(fulldata,K,guess)       
	pars <- matrix(0,nrow=nitems, ncol=nfact + 1)
	if (is.null(startvalues)){  
		suppressMessages(startvalues <- start.values(fulldata,guess,Rpoly,bfactor=TRUE))
		pars[logicalfact] <- startvalues    
		sload <- startvalues[(nitems+1):(2*nitems)]
		for(i in 1:nitems){
			temp <- selvec[pars[i,2:nfact] > 0]
			pars[i,temp] <- sload[i]
		}  
	} else {
		if (ncol(startvalues) == 3){
			pars[logicalfact] <- startvalues	  	
			sload <- startvalues[(nitems+1):(2*nitems)]
			for(i in 1:nitems){
				temp <- selvec[pars[i,2:nfact] > 0]
				pars[i,temp] <- sload[i]
			} 
		} else pars <- startvalues
	}
	diag(Rpoly) <- 1
	item <- 1
	lastpars2 <- lastpars1 <- rate <- matrix(0,nrow=nitems,ncol=ncol(pars))  
	prior <- dnorm(theta) 
	Prior <- dmvnorm(Theta,rep(0,2),diag(2))  
	startvalues <- pars  
	converge <- 1
	problemitems <- c()
	index <- 1:nitems
	sitems <- matrix(0,ncol=nitems,nrow=(nfact-1))
	for(i in 1:nitems) sitems[specific[i],i] <- 1      
	if(debug){
		print(startvalues)
		print(sitems)	 
	} 

	#EM  loop  
	for (cycles in 1:ncycles) 
	{    
		rlist <- Estep.bfactor(pars, tabdata, Theta, prior, guess, logicalfact, specific, sitems)
		if(debug) print(sum(r*log(rlist[[3]])))
		lastpars2 <- lastpars1
		lastpars1 <- pars	
		mpars <- matrix(pars[logicalfact], ncol=3)    
		temp <- rowSums(pars[,2:nfact])
		mpars[ ,2] <- temp	
		for(i in 1:nitems){ 
			if(guess[i] == 0)	
				maxim <- try(optim(mpars[i, ],fn=fn,gr=gr,r1=rlist[[1]][i, ],N=rlist[[2]][i, ],
					guess=guess[i],Theta=Theta,prior=Prior,parprior=par.prior[i, ],method="BFGS"))
			else	  
				maxim <- try(optim(mpars[i, ],fn=fn,r1=rlist[[1]][i, ],N=rlist[[2]][i, ],
					guess=guess[i],Theta=Theta,prior=Prior,parprior=par.prior[i, ],method="BFGS"))
			if(class(maxim) == "try-error") {
				problemitems <- c(problemitems, i)	  
				converge <- 0
				next
			}	
			mpars[i, ] <- maxim$par	  
		}			
		sload <- mpars[(nitems+1):(2*nitems)]		
		for(i in 1:nitems){
			temp <- selvec[pars[i,2:nfact] > 0]
			pars[i,temp] <- sload[i]
		}  
		pars[ ,1] <- mpars[,1]
		pars[ ,nfact+1] <- mpars[ ,3]		
		pars[is.na(pars)] <- lastpars1[is.na(pars)]
		if (max(abs(lastpars1 - pars)) < EMtol) break
		if(!suppressAutoPrior){
			if(any(abs(pars[ ,nfact+1]) > 4)){
				ints <- index[abs(pars[ ,nfact+1]) > 4] 	
				par.prior[ints,3] <- 2
				if(any(abs(pars[ ,nfact+1]) > 5.5)){
					ints <- index[abs(pars[ ,nfact+1]) > 5.5] 	
					par.prior[ints,3] <- 1
				}
			}
			norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2))
			alp <- as.matrix(pars[ ,1:nfact]/norm)
			FF <- alp %*% t(alp)
			V <- eigen(FF)$vector[ ,1:nfact]
			L <- eigen(FF)$values[1:nfact]
			F <- as.matrix(V * sqrt(L))
			F <- V %*% sqrt(diag(L))
			h2 <- rowSums(F^2)
			if(any(h2 > .95)){
				if(any(h2 > .95)){
					ind <- index[h2 > .95]
					par.prior[ind,1] <- 1.2
				} 
				if(any(h2 > .98)){  
					ind <- index[h2 > .98]
					par.prior[ind,1] <- 1.5		
				}
			}
		}	
	# apply rate acceleration every third cycle    
		if (cycles %% 3 == 0 & cycles > 6){
			d1 <- lastpars1 - pars
			d2 <- lastpars2 - pars      
			for (i in 1:nitems) {
				for(j in 1:ncol(pars)){      
					if((abs(d1[i,j]) > 0.001) & (d1[i,j]*d2[i,j] > 0.0) & (d1[i,j]/d2[i,j] < 1.0))
						rate[i,j] <- (1 - (1 - rate[i,j]) * (d1[i,j]/d2[i,j]))
						else rate[i,j] <- 0
				}        
			}      
		}  
		rate[pars > 4] <- 0
		rate[pars < -4] <- 0
		pars <- lastpars1*rate*(-2) + (1 - rate*(-2))*pars          
	}

	if(any(par.prior[,1] != 1)) cat("Slope prior for item(s):",
		as.character(index[par.prior[,1] > 1]), "\n")
	if(any(par.prior[,3] != 0)) cat("Intercept prior for item(s):",
		as.character(index[par.prior[,3] > 0]), "\n")  
	if(converge == 0) 
		warning("Parameter estimation reached unacceptable values. 
		Model probably did not converge.")
	if(length(problemitems) > 0) warning("Problem with the M-step for item(s): ", 
		paste(unique(problemitems), " "))	
	lastchange <- abs(lastpars1 - pars)
	if (cycles == ncycles){ 
		converge <- 0
		message("Estimation terminated after ", cycles, " EM loops. Maximum changes: 
			\n slopes = ", round(max(abs(lastchange[,1:nfact])),4), ", intercepts = ", 
			round(max(abs(lastchange[,ncol(pars)])),4) ,"\n")
	}	
	rlist <- Estep.bfactor(pars, tabdata, Theta, prior, guess, logicalfact, specific, sitems)
	Pl <- rlist[[3]]
	log.lik <- sum(r * log(Pl))
	logN <- 0
	logr <- rep(0,length(r))
	for (i in 1:sampsize) logN <- logN + log(i)
	for (i in 1:length(r)) 
	for (j in 1:r[i]) 
	logr[i] <- logr[i] + log(j)	
	log.lik <- log.lik + logN/sum(logr)
	AIC <- (-2) * log.lik + 3 * length(specific)
	X2 <- 2 * sum(r * log(r / (sampsize*Pl)))  
	df <- length(r) - 1 - 2*nitems - length(specific)
	p <- 1 - pchisq(X2,df)

	#from last EM cycle pars to FA
	norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2))
	gam <- (-1)*pars[ ,nfact + 1]/norm  
	F <- matrix(0,ncol = nfact, nrow = nitems)
	for (i in 1:nitems) F[i,1:nfact] <- pars[i,1:nfact]/norm[i]  
	h2 <- rowSums(F^2)  

	mod <- new('bfactorClass',EMiter=cycles, pars=pars, guess=guess, AIC=AIC, X2=X2, 
		df=df, log.lik=log.lik, p=p, F=F, h2=h2, itemnames=itemnames, 
		tabdata=tabdata, sampsize=sampsize, Pl=Pl, Theta=Theta, fulldata=fulldata, 
		logicalfact=logicalfact, facility=facility, specific=specific,
		cormat=Rpoly, converge=converge, par.prior=par.prior, quadpts=quadpts,Call=Call)  
	return(mod)  
}
