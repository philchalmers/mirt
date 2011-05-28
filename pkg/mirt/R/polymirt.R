plot.polymirt <- function(x, npts = 50,
  rot = list(x = -70, y = 30, z = 10), ...)
{  
	type = 'curve'
	K <- x$K		
	nfact <- ncol(x$Theta)
	a <- as.matrix(x$pars[ ,1:nfact])
    d <- as.matrix(x$pars[ ,(nfact+1):ncol(x$pars)])	
	guess <- x$guess
	guess[is.na(guess)] <- 0
	A <- as.matrix(sqrt(apply(a^2,1,sum)))	
	theta <- seq(-4,4,length.out=npts)
	Theta <- thetaComb(theta, nfact)
	info <- rep(0,nrow(Theta))
	for(j in 1:length(K)){
		if(K[j] > 2){
			P <- P.poly(a[j,], d[j,],Theta, itemexp = FALSE)		
			for(i in 1:K[j]){
				w1 <- P[,i]*(1-P[,i])*A[j]
				w2 <- P[,i+1]*(1-P[,i+1])*A[j]
				I <- ((w1 - w2)^2) / (P[,i] - P[,i+1]) * P[,i]
				info <- info + I
			}
		} else {
			P <- P.mirt(a[j,], d[j,],Theta, guess[j])
			Pstar <- P.mirt(a[j,], d[j,],Theta, 0)
			info <- info + A[j]^2 * P * (1-P) * Pstar/P
		}			
	}		
	plt <- cbind(info,Theta)
	if(nfact > 1){
		require(lattice)
		wireframe(info ~ Theta[ ,1] + Theta[ ,2], data = plt, main = "Item Information", 
			zlab = "I", xlab = "Theta 1", ylab = "Theta 2", scales = list(arrows = FALSE),
			screen = rot)	
	} else 
		plot(Theta, info, type='l',main = 'Item Information', xlab = 'Theta', ylab='Information')
}	  

coef.polymirt <- function(object, SE = TRUE, digits = 3, ...)
{  
	nfact <- ncol(object$Theta)	
	a <- matrix(object$pars[ ,1:nfact],ncol=nfact)
	d <- matrix(object$pars[,(nfact+1):ncol(object$pars)],
		ncol = ncol(object$pars)-nfact)    
	A <- sqrt(apply(a^2,1,sum))
	B <- -d/A  
	if (nfact > 1){  
		parameters <- cbind(object$pars,object$guess,A,B)
		SEs <- object$SEpars	
		colnames(parameters) <- c(paste("a_",1:nfact,sep=""),paste("d_",1:(ncol(object$pars)-nfact),sep=""),"guess","mvdisc",paste("mvint_",1:(ncol(object$pars)-nfact),sep=""))	
		colnames(SEs) <- c(paste("a_",1:nfact,sep=""),paste("d_",1:(ncol(object$pars)-nfact),sep=""),"guess")		
		cat("Unrotated parameters, multivariate discrimination and intercept: \n")
		print(round(parameters, digits))
		if(SE){
			cat("\nStd. Errors: \n")	
			print(round(SEs, digits))
		}				
	} else {
		parameters <- cbind(object$pars,object$guess)
		SEs <- object$SEpars	
		colnames(parameters) <- colnames(SEs) <- c(paste("a_",1:nfact,sep=""),paste("d_",1:(ncol(object$pars)-nfact),sep=""),"guess")			
		cat("Parameter slopes and intercepts: \n")	
		print(round(parameters, digits))
		if(SE){
			cat("\nStd. Errors: \n")	
			print(round(SEs, digits))
		}
	}
	invisible(parameters)
}

summary.polymirt <- function(object, rotate = 'varimax', digits = 3, ...)
{
	nfact <- ncol(object$F)
	if (rotate == 'none' || nfact == 1) {
		F <- object$F
		h2 <- as.matrix(object$h2)    	
		SS <- apply(F^2,2,sum)
		colnames(h2) <- "h2"	
		colnames(F) <- names(SS) <- paste("F_", 1:ncol(F),sep="")
		cat("\nUnrotated factor loadings: \n")
		loads <- round(cbind(F,h2),digits)
		print(loads)	    	 
		cat("\nSS loadings: ",round(SS,digits), "\n")
		cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
		invisible(list(F,h2))
	} else {	
		F <- object$F
		h2 <- as.matrix(object$h2)		
		colnames(F) <- paste("F_", 1:ncol(F),sep="")
		colnames(h2) <- "h2"		
		cat("Rotation: ", rotate, "\n")
		rotF <- Rotate(F,rotate)
		SS <- apply(rotF$loadings^2,2,sum)
		loads <- round(cbind(rotF$loadings,h2),digits)		
		cat("\nRotated factor loadings: \n")
		print(loads)		
		if(attr(rotF, "oblique")){
			cat("\nFactor correlations: \n")
			Phi <- rotF$Phi	  
			Phi <- round(Phi, digits)
			colnames(Phi) <- rownames(Phi) <- colnames(F)
			print(Phi)
			cat("\nRotated Sums of Squares: \n")
			round(colSums(rotF$loadings %*% Phi), digits)      
		} else {	
			cat("\nSS loadings: ",round(SS,digits), "\n")
			cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")	
		}
		if(any(h2 > 1)) 
			warning("Solution has heywood cases. Interpret with caution.") 
		invisible(list(loadings,h2))  
	}  
}

print.polymirt <- function(x, ...){
	cat("Call: ")
	print(x$Call)
	cat("\nFull-information factor analysis with ", ncol(x$F), " factor",
		if(ncol(x$F)>1) "s", "\n", sep="")
	if(x$converge == 1)	
		cat("Converged in ", x$cycles, " iterations.\n", sep="")
	else 	
		cat("Estimation stopped after ", x$cycles, " iterations.\n", sep="")	
} 


polymirt <- function(data, nfact, guess = 0, prev.cor = NULL, 
	ncycles = 2000, burnin = 200, SEM.cycles = 100, kdraws = 1, 
	tol = .001, debug = FALSE, ...){
		
	Call <- match.call()   
	itemnames <- colnames(data)
	data <- as.matrix(data)	
	if(any(is.na(data))) stop("polymirt function can't handle missing data.\n")  
	J <- ncol(data)
	N <- nrow(data)	
	colnames(data) <- itemnames
	if(length(guess) == 1) guess <- rep(guess,J)
	if(length(guess) > J || length(guess) < J) 
		stop("The number of guessing parameters is incorrect.")
	estGuess <- guess > 0					
	uniques <- list()
	for(i in 1:J)
		uniques[[i]] <- sort(unique(data[,i]))
	K <- rep(0,J)
	for(i in 1:J) K[i] <- length(uniques[[i]])
	guess[K > 2] <- 0
	estGuess[K > 2] <- FALSE	
	itemloc <- cumsum(c(1,K))
	index <- 1:J	
	fulldata <- matrix(0,N,sum(K))
	Names <- NULL
	for(i in 1:J)
        Names <- c(Names, paste("Item.",i,"_",1:K[i],sep=""))				
	colnames(fulldata) <- Names			
	for(i in 1:J){
		ind <- index[i]
		if(setequal(uniques[[i]], c(0,1))){
			fulldata[ ,itemloc[ind]:(itemloc[ind]+1)] <- cbind(data[,ind],abs(1-data[,ind]))
			next
		}
		dummy <- matrix(0,N,K[ind])
		for (j in 0:(K[ind]-1))  
			dummy[,j+1] <- as.integer(data[,ind] == uniques[[ind]][j+1])  		
		fulldata[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- dummy		
	}			
	if(!is.null(prev.cor)){
		if (ncol(prev.cor) == nrow(prev.cor)) Rpoly <- prev.cor
			else stop("Correlation matrix is not square.\n")
	} 	else Rpoly <- cormod(data,K,guess)
	FA <- factor.minres(Rpoly,nfact,rotate = 'none', warnings= FALSE)	
	loads <- unclass(loadings(FA))
	u <- FA$unique
	u[u < .001 ] <- .2
	cs <- sqrt(u)
	lambdas <- loads/cs
	zetas <- rep(0,ncol(fulldata) - J)
	loc <- 1	
	for(i in 1:J){
		if(K[i] == 2){
			zetas[loc] <- qnorm(mean(fulldata[,itemloc[i]]))/cs[i]
			loc <- loc + 1
		} else {			
			temp <- table(data[,i])[1:(K[i]-1)]/N
			temp <- cumsum(temp)			
			zetas[loc:(loc+K[i]-2)] <- qnorm(1 - temp)/cs[i]	
			loc <- loc + K[i] - 1	
		}		
	}	
	npars <- length(c(lambdas,zetas)) + sum(estGuess) #start here
	parind <- 1:npars
	pars <- rep(NA,npars)
	Ksum <- cumsum(K + nfact - 1 + estGuess) - (nfact-1)
	lamind	<- gind <- c()	 
	for(i in 1:J){
		pars[Ksum[i]:(Ksum[i] + nfact - 1)] <- lambdas[i,]
		lamind <- c(lamind,Ksum[i]:(Ksum[i] + nfact - 1))
		if(estGuess[i]){
			pars[Ksum[i] - 2] <- guess[i]
			gind <- c(gind,Ksum[i] - 2)
		}	
	}	
	zetaind <- parind[is.na(pars)]			
	pars[is.na(pars)] <- zetas
	diag(Rpoly) <- 1	
	converge <- 1    	
	if(debug){
		print(lambdas)
		print(zetas)
	}	
	
    #preamble for MRHM algorithm			
	theta0 <- matrix(0,N,nfact)	    
	cand.t.var <- 1	
	for(i in 1:20){
		theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)
		if(attr(theta0,"Proportion Accepted") > .5) cand.t.var <- cand.t.var + .1 
		else if(attr(theta0,"Proportion Accepted") < .35) cand.t.var <- cand.t.var - .1     	
	}	
	m.thetas <- list()		
	SEM.stores <- matrix(0,SEM.cycles,npars)
	phi <- rep(0,npars)
	Tau <- info <- h <- matrix(0,npars,npars)
	m.list <- list()	  
	conv <- 0
	k <- 1	
	gamma <- 0.1
	startvalues <- pars
	stagecycle <- 1	
	
	for(cycles in 1:(ncycles + burnin + SEM.cycles))
	{ 
		if(cycles == burnin + 1) stagecycle <- 2
		if(stagecycle == 3)
			gamma <- sqrt(0.1/(2.5*(cycles - SEM.cycles - burnin - 1)))
		if(cycles == (burnin + SEM.cycles + 1)){ 
			stagecycle <- 3		
		    pars <- rep(0,npars)
			for(i in 1:SEM.cycles) pars <- pars + SEM.stores[i,]
			pars <- pars/SEM.cycles	
			k <- kdraws	
			gamma <- 1
		}		
		lambdas <- matrix(pars[lamind],ncol=nfact,byrow=TRUE)
		zetas <- pars[zetaind]
		guess <- rep(0,J)
		guess[estGuess] <- pars[gind]
		for(j in 1:5) theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)
		
		#Step 1. Generate m_k datasets of theta 
		for(i in 1:k)			
			m.thetas[[i]] <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)		
		theta0 <- m.thetas[[1]]
		
		#Step 2. Find average of simulated data gradients and hessian 
		g.m <- h.m <- list()					
		for(j in 1:k){
			g <- rep(NA,npars)
			loc <- 1
			for(i in 0:(J - 1)){
				if(K[i+1]==2){
					temp <- dpars.dich(lambdas[i+1,],zetas[loc],guess[i+1],
						fulldata[,itemloc[i+1]],m.thetas[[j]], estGuess[i+1])
					ind <- parind[is.na(g)][1]
					ind2 <- ind+nfact+ estGuess[i+1]		
					g[ind:ind2] <- temp$grad
					h[ind:ind2,ind:ind2] <- temp$hess
					loc <- loc + 1
				} else {
					loc2 <- loc + K[i+1] - 2
					temp <- dpars.poly(lambdas[i+1,],zetas[loc:loc2],
						fulldata[,itemloc[i+1]:(itemloc[i+2]-1)],m.thetas[[j]])
					ind <- parind[is.na(g)][1]	
					ind2 <- ind+nfact+K[i+1]-2
					g[ind:ind2] <- temp$grad
					h[ind:ind2,ind:ind2] <- temp$hess
					loc <- loc + K[i+1] - 1				
				}
			} 
			g.m[[j]] <- g
			h.m[[j]] <- h
		}
		ave.g <- rep(0,length(g))
		ave.h <- matrix(0,length(g),length(g))		
		for(i in 1:k){
		  ave.g <- ave.g + g.m[[i]]
		  ave.h <- ave.h + h.m[[i]]
		}
		grad <- ave.g/k
		ave.h <- (-1)*ave.h/k 
		if(stagecycle < 3){
		    correction <- solve(ave.h) %*% grad
			correction[correction > .5] <- .5
			correction[correction < -.5] <- -.5			
			parsold <- pars
			pars <- pars + gamma*correction
			pars[pars[gind] < 0] <- parsold[pars[gind] < 0]
			pars[pars[gind] > .4] <- parsold[pars[gind] > .4]	
			if(stagecycle == 2) SEM.stores[cycles - burnin,] <- pars
			next
		}	
		
		#Step 3. Update R-M step		
		Tau <- Tau + gamma*(ave.h - Tau)		
		correction <- (solve(Tau) %*% grad)	
		correction[correction > .5] <- .5
		correction[correction < -.5] <- -.5	
		if(any(estGuess)){
			correction[correction[gind] > .05] <- .05
			correction[correction[gind] < -.05] <- -.05
		}				
		if(conv == 3) break				
		pars <- pars + gamma*correction
		if(all(abs(parsold - pars) < tol)) conv <- conv + 1
			else conv <- 0	
		parsold <- pars	
		pars[pars[gind] < 0] <- parsold[pars[gind] < 0]
		pars[pars[gind] > .4] <- parsold[pars[gind] > .4]
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE		
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	}
	
	SE <- diag(solve(info))
	if(any(SE < 0)){
		warning("Information matrix is not positive definite.\n")
		SE <- rep(0,npars)
	}
	if(any(guess < 0)) warning("Negative lower asymptote parameter(s). \n")		
	SE <- sqrt(SE)	
	lambdas <- matrix(pars[lamind],ncol=nfact,byrow=TRUE)
	SElam <- matrix(SE[lamind],ncol=nfact,byrow=TRUE)
	SEg <- guess <- rep(NA,J)
	guess[estGuess] <- pars[gind]
	SEg[estGuess] <- SE[gind]
	zetas <- SEzeta <- matrix(NA,J,(max(K)-1))
	temp <- pars[zetaind]
	temp1 <- SE[zetaind]	
	k <- 1
	for(i in 1:J){
		for(j in 1:(K[i]-1)){
			zetas[i,j] <- temp[k] 
			SEzeta[i,j] <- temp1[k]
			k <- k + 1
		}
	}	 
	guess[K > 2] <- NA
	pars <- cbind(lambdas,zetas)
	SEpars <- cbind(SElam,SEzeta,SEg)
	
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
	
	
	mod <- list(pars=pars, guess=guess, SEpars=SEpars, cycles=cycles - SEM.cycles,
		Theta=theta0, fulldata=fulldata,K=K, F=F, h2=h2, fulldata=fulldata, 
		itemloc=itemloc, converge = converge, Call=Call)	 
	class(mod) <- 'polymirt'
	mod	
}
