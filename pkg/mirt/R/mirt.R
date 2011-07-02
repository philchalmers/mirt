########################################## 
#Methods 

setMethod(
	f = "print",
	signature = signature(x = 'mirtClass'),
	definition = function(x, ...){  
		cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information factor analysis with ", ncol(x@F), " factor",
			if(ncol(x@F)>1) "s", "\n", sep="")
			if(x@converge == 1)	
				cat("Converged in ", x@EMiter, " iterations using ", x@quadpts,
				" quadrature points.\n", sep="")
		else 	
			cat("Estimation stopped after ", x@EMiter, " iterations using ", 
				x@quadpts, " quadrature points.\n", sep="")
		cat("Log-likelihood =", x@log.lik, "\n")
		cat("AIC =", x@AIC, "\n")
		cat("G^2 = ", round(x@X2,2), ", df = ", 
			x@df, ", p = ", round(x@p,4), "\n", sep="")
	}
)

setMethod(
	f = "show",
	signature = signature(object = 'mirtClass'),
	definition = function(object){  
		cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information factor analysis with ", ncol(object@F), " factor",
			if(ncol(object@F)>1) "s", "\n", sep="")
			if(object@converge == 1)	
				cat("Converged in ", object@EMiter, " iterations using ", object@quadpts,
				" quadrature points.\n", sep="")
		else 	
			cat("Estimation stopped after ", object@EMiter, " iterations using ", 
				object@quadpts,	" quadrature points.\n", sep="")
		cat("Log-likelihood =", object@log.lik, "\n")
		cat("AIC =", object@AIC, "\n")
		cat("G^2 = ", round(object@X2,2), ", df = ", 
			object@df, ", p = ", round(object@p,4), "\n", sep="")
	}
)

setMethod(
	f = "summary",
	signature = 'mirtClass',
	definition = function(object, rotate = 'varimax', suppress = 0, digits = 3, ...){
		nfact <- ncol(object@F)
		if (rotate == 'none' || nfact == 1) {
			F <- object@F
			F[abs(F) < suppress] <- NA
			h2 <- as.matrix(object@h2)
			fac <- as.matrix(object@facility)	
			SS <- apply(F^2,2,sum)
			colnames(h2) <- "h2"
			colnames(fac) <- "facility"
			colnames(F) <- names(SS) <- paste("F_", 1:ncol(F),sep="")
			cat("\nUnrotated factor loadings: \n\n")
			loads <- round(cbind(F,h2,fac),digits)
			print(loads)	    	 
			cat("\nSS loadings: ",round(SS,digits), "\n")
			cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
			invisible(list(F,h2))
		} else {	
			F <- object@F
			h2 <- as.matrix(object@h2)
			fac <- as.matrix(object@facility)
			colnames(F) <- paste("F_", 1:ncol(F),sep="")
			colnames(h2) <- "h2"				
			cat("\nRotation: ", rotate, "\n")
			rotF <- Rotate(F,rotate)
			SS <- apply(rotF$loadings^2,2,sum)
			L <- rotF$loadings
			L[abs(L) < suppress] <- NA	
			loads <- round(cbind(L,h2),digits)			
			cat("\nRotated factor loadings: \n\n")
			print(loads,digits)		
			if(attr(rotF, "oblique")){
				cat("\nFactor correlations: \n")
				Phi <- rotF$Phi	  
				Phi <- round(Phi, digits)
				colnames(Phi) <- rownames(Phi) <- colnames(F)
				print(Phi)            
			}	
			cat("\nRoateted SS loadings: ",round(SS,digits), "\n")		
			if(any(h2 > 1)) 
				warning("Solution has heywood cases. Interpret with caution.") 
			invisible(list(rotF$loadings,h2))  
		}  
	}
)

setMethod(
	f = "coef",
	signature = 'mirtClass',
	definition = function(object, digits = 3, ...){  
		a <- as.matrix(object@pars[ ,1:(ncol(object@pars)-1)])
		d <- object@pars[ ,ncol(object@pars)]
		A <- sqrt(apply(a^2,1,sum))
		B <- -d/A  
		if (ncol(a) > 1){  
			parameters <- cbind(object@pars,object@guess,object@facility,A,B)    
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),"d","guess", 
			"facility","mvdisc","mvint")	  
			cat("\nUnrotated parameters, multivariate discrimination and intercept: \n\n")
			print(round(parameters, digits))  	
		} else {
			parameters <- cbind(object@pars,object@guess,object@facility) 
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),"d","guess","facility")    
			cat("\nParameter slopes and intercepts: \n\n")	
			print(round(parameters, digits))	  
		}
		invisible(parameters)
	}
)

setMethod(
	f = "anova",
	signature = signature(object = 'mirtClass'),
	definition = function(object, object2, ...){
		dots <- list(...)		
		df <- object@df - object2@df  
		X2 <- 2*object2@log.lik - 2*object@log.lik 
		AICdiff <- object@AIC - object2@AIC    
		cat("\nChi-squared difference: \n\nX2 = ", round(X2,3), ", df = ",
		df, ", p = ", round(1 - pchisq(X2,df),4), "\n", sep="")
		cat("AIC difference = ", round(AICdiff,3), "\n")  
	}
)

setMethod(
	f = "residuals",
	signature = signature(object = 'mirtClass'),
	definition = function(object, type = 'LD', digits = 3, ...){   	
		Theta <- object@Theta
		fulldata <- object@fulldata	
		N <- nrow(fulldata)	
		J <- ncol(fulldata)
		nfact <- ncol(object@F)
		lambdas <- matrix(object@pars[,1:nfact], J)
		zetas <- object@pars[,(nfact+1)]
		guess <- object@guess
		guess[is.na(guess)] <- 0		
		if(type == 'LD'){
			res <- matrix(0,J,J)
			diag(res) <- NA
			colnames(res) <- rownames(res) <- colnames(fulldata)
			prior <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
			prior <- prior/sum(prior)
			for(i in 1:J){			
				for(j in 1:J){
					if(i < j){
						P1 <- P.mirt(lambdas[i,],zetas[i], Theta, guess[i])
						P2 <- P.mirt(lambdas[j,],zetas[j], Theta, guess[j])
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
			print(res,digits)			
		}		
		if(type == 'exp'){	
			r <- object@tabdata[ ,ncol(object@tabdata)]
			res <- (r - object@Pl * nrow(object@fulldata)) / 
				sqrt(object@Pl * nrow(object@fulldata))	
			cat("\nExpected values:\n\n")
			print(res,digits)		
		}			
		return(invisible(res))	
	}
)

setMethod(
	f = "plot",
	signature = signature(x = 'mirtClass', y = "missing"),
	definition = function(x, y, type = 'info', npts = 50, 
		rot = list(xaxis = -70, yaxis = 30, zaxis = 10))
	{  
		if (!type %in% c('curve','info')) stop(type, " is not a valid plot type.")
		rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
		a <- as.matrix(x@pars[ ,1:(ncol(x@pars) - 1)])
		d <- x@pars[ ,ncol(x@pars)]
		g <- x@guess
		A <- as.matrix(sqrt(apply(a^2,1,sum)))
		B <- -d/A
		if(ncol(a) > 2 ) stop("Can't plot high dimentional solutions.\n")
		theta <- seq(-4,4,length.out=npts)
		Theta <- thetaComb(theta, ncol(a))
		P <- matrix(0, ncol=length(g), nrow = nrow(as.matrix(Theta)))
		for(i in 1:nrow(a)) P[ ,i] <- P.mirt(a[i, ],d[i],as.matrix(Theta),g[i])  
		Ptot <- rowSums(P) 
		if(ncol(a) == 2){
			require(lattice)
			if(type == 'info'){
				I <- (P * (1 - P)) %*% A^2 
				plt <- data.frame(cbind(I,Theta))
				colnames(plt) <- c("I", "Theta1", "Theta2")
				wireframe(I ~ Theta1 + Theta2, data = plt, main = "Test Information", 
					zlab = "I", xlab = "Theta 1", ylab = "Theta 2", 
					scales = list(arrows = FALSE), screen = rot)
			} else {  
				plt <- data.frame(cbind(Ptot,Theta))
				colnames(plt) <- c("Ptot", "Theta1", "Theta2")
				wireframe(Ptot ~ Theta1 + Theta2, data = plt, main = "Test score surface", 
					zlab = "Test \nScore", xlab = "Theta 1", ylab = "Theta 2", 
					scales = list(arrows = FALSE), screen = rot)
			}	
		} else {
			if(type == 'curve')  
				plot(Theta, Ptot, type='l', main = 'Test score plot', xlab = 'Theta', ylab='Test Score')
			 else {
				I <- (P * (1 - P)) %*% a^2 
				plot(Theta, I, type='l', main = 'Test Information', xlab = 'Theta', ylab='Information')
			} 	
		}  
	}		
)	

setMethod(
	f = "fitted",
	signature = signature(object = 'mirtClass'),
	definition = function(object, digits = 3, ...){  
		expected <- round(nrow(object@fulldata) * object@Pl,digits)  
		tabdata <- cbind(object@tabdata,expected)
		colnames(tabdata) <- c(colnames(object@fulldata),"freq","exp")	
		print(tabdata)
		invisible(expected)
	}
)
    
############################################
#Main function

mirt <- function(fulldata, nfact, guess = 0, prev.cor = NULL, par.prior = FALSE, 
	startvalues = NULL, quadpts = NULL, ncycles = 300, tol = .001, nowarn = TRUE, 
	debug = FALSE, ...)
{ 
	fn <- function(pars, r1, N, guess, Theta, prior, parprior){
		a <- pars[1:(length(pars)-1)]
		d <- pars[length(pars)]		
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
  
	Call <- match.call()    
	itemnames <- colnames(fulldata)
	fulldata <- as.matrix(fulldata)
	if (!any(fulldata %in% c(0,1,NA))) stop("Data must contain only 0, 1, or NA.")
	fulldata[is.na(fulldata)] <- 0  
	nitems <- ncol(fulldata)  
	colnames(fulldata) <- itemnames
	if (length(guess) == 1) guess <- rep(guess,nitems)
		else if (length(guess) > nitems || length(guess) < nitems) 
	stop("The number of guessing parameters is incorrect.")	
	pats <- apply(fulldata,1,paste,collapse = "/")
	freqs <- table(pats)
	nfreqs <- length(freqs)
	K <- rep(2,nitems)
	r <- as.vector(freqs)
	sampsize <- nrow(fulldata) 
	tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
	tabdata <- matrix(as.numeric(tabdata),nfreqs,nitems,TRUE)
	tabdata <- cbind(tabdata,r)    
	if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))  
	theta <- as.matrix(seq(-4,4,length.out = quadpts))
	if(quadpts^nfact < 10000){
		Theta <- thetaComb(theta,nfact)
		prior <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
		prior <- prior/sum(prior)
	} else stop('Greater than 10000 quadrature points, reduce number.')
	facility <- colMeans(fulldata)
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
	if (any(class(prev.cor) == c('mirt','bmirt'))) Rpoly <- prev.cor$cormat
		else if(!is.null(prev.cor)) {
		if (ncol(prev.cor) == nrow(prev.cor)) Rpoly <- prev.cor
			else stop("Correlation matrix is not square.\n")
	} else 
		Rpoly <- cormod(fulldata,K,guess)   
	if (is.null(startvalues)){ 
		suppressMessages(pars <- start.values(fulldata,guess,Rpoly,
			nfact=nfact,nowarn=nowarn))
		pars[pars > 3] <- 3
		pars[pars < -3] <- -3	
	} else {
		if ((ncol(startvalues) != (nfact + 1)) || (nrow(startvalues) != nitems))
			stop("Startvalues are declared incorrectly.")  
		pars <- startvalues  
	} 
	diag(Rpoly) <- 1
	item <- 1
	lastpars2 <- lastpars1 <- rate <- matrix(0,nrow=nitems,ncol=ncol(pars))    
	startvalues <- pars
	converge <- 1  
	problemitems <- c()
	index <- 1:nitems  
	if(debug) print(startvalues)      

	# EM loop
	for (cycles in 1:ncycles)
	{       
		rlist <- Estep.mirt(pars,tabdata,Theta,prior,guess)
		prior <- rlist[[4]]
		if (debug) print(sum(r*log(rlist[[3]])))
		lastpars2 <- lastpars1
		lastpars1 <- pars	
		for(i in 1:nitems){
			if(guess[i] == 0)	
				maxim <- try(optim(pars[i, ],fn=fn,gr=gr,r1=rlist[[1]][i, ],N=rlist[[2]][i, ],
				guess=guess[i],Theta=Theta,prior=prior,parprior=par.prior[i, ],method="BFGS"))
			else 
				maxim <- try(optim(pars[i, ],fn=fn,r1=rlist[[1]][i, ],N=rlist[[2]][i, ],
				guess=guess[i],Theta=Theta,prior=prior,parprior=par.prior[i, ],method="BFGS"))
			if(class(maxim) == "try-error"){
				problemitems <- c(problemitems, i)
				converge <- 0
				next
			}		  
			pars[i, ] <- maxim$par	  
		}	
		if(!suppressAutoPrior){
			if(any(abs(pars[ ,nfact+1]) > 4)){
				ints <- index[abs(pars[ ,nfact+1]) > 4] 	
				par.prior[ints,3] <- 2
				if(any(abs(pars[ ,nfact+1]) > 5.5)){
					ints <- index[abs(pars[ ,nfact+1]) > 5.5] 	
					par.prior[ints,3] <- 1
				} 
			}
			if(nfact > 1){ 
				norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2))
				alp <- as.matrix(pars[ ,1:nfact]/norm)
				FF <- alp %*% t(alp)
				V <- eigen(FF)$vector[ ,1:nfact]
				L <- eigen(FF)$values[1:nfact]      
				F <- V %*% sqrt(diag(L))
				h2 <- rowSums(F^2)
				if(any(h2 > .95)){
					ind <- index[h2 > .95]
					par.prior[ind,1] <- 1.2
					if(any(h2 > .98)){
						ind <- index[h2 > .98]
						par.prior[ind,1] <- 1.5
					} 
				}
			}
		}  
		maxdif <- max(abs(lastpars1 - pars))	
		if (maxdif < tol) break    
		# rate acceleration adjusted every third cycle
		if (cycles %% 3 == 0 & cycles > 6) 
		{
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
		Model probably did not converged.")  
	if(length(problemitems) > 0) warning("Problem with the M-step for item(s): ", 
		paste(unique(problemitems), " "))	
	lastchange <- lastpars1 - pars
	if (cycles == ncycles){
		converge <- 0  
		message("Estimation terminated after ", cycles, " EM loops. Maximum changes:") 
		message("\n slopes = ", round(max(abs(lastchange[ ,1:nfact])),4), ", intercepts = ", 
			round(max(abs(lastchange[ ,ncol(pars)])),4) ,"\n", sep="")
	}	    
	prior <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
	prior <- prior/sum(prior)  
	rlist <- Estep.mirt(pars,tabdata,Theta,prior,guess)      	  
	Pl <- rlist[[3]]  
	log.lik <- sum(r*log(Pl))
	logN <- 0
	logr <- rep(0,length(r))
	for (i in 1:sampsize) logN <- logN + log(i)
	for (i in 1:length(r)) 
		for (j in 1:r[i]) 
			logr[i] <- logr[i] + log(j)    
	df <- (length(r) - 1) - nitems*(nfact + 1) + nfact*(nfact - 1)/2 
	X2 <- 2 * sum(r * log(r/(sampsize*Pl)))
	log.lik <- log.lik + logN/sum(logr)	
	p <- 1 - pchisq(X2,df)  
	AIC <- (-2) * log.lik + 2 * length(pars)

	# pars to FA loadings
	if (nfact > 1) norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2))
		else norm <- as.matrix(sqrt(1 + pars[ ,1]^2))  
	alp <- as.matrix(pars[ ,1:nfact]/norm)
	FF <- alp %*% t(alp)
	V <- eigen(FF)$vector[ ,1:nfact]
	L <- eigen(FF)$values[1:nfact]
	if (nfact == 1) F <- as.matrix(V * sqrt(L))
		else F <- V %*% sqrt(diag(L))  
	if (sum(F[ ,1] < 0)) F <- (-1) * F  
	h2 <- rowSums(F^2) 

	mod <- new('mirtClass', EMiter=cycles, pars=pars, guess=guess, X2=X2, df=df, p=p,
		AIC=AIC, log.lik=log.lik, F=F, h2=h2, tabdata=tabdata, Theta=Theta, Pl=Pl, 
		fulldata=fulldata, cormat=Rpoly, facility=facility, converge=converge, quadpts=quadpts,
		Call=Call)	  
	return(mod)    
}


