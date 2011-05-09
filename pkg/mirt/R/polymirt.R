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
			for(i in 1:K[i]){
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

coef.polymirt <- function(object, digits = 3, SE = TRUE, ...)
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
		cat("Parameters with multivariate discrimination and intercept: \n")	
		print(round(parameters, digits))
		if(SE){
			cat("\nStd. Errors: \n")	
			print(round(SEs, digits))
		}
	}
	invisible(parameters)
}

summary.polymirt <- function(object, digits = 3, rotate = 'varimax', ...)
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
	ncycles = 2000, SEM.cycles = 100, kdraws = 1, tol = .0005, debug = FALSE, ...){
	
	draw.thetas <- function(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var) { 		
		N <- nrow(fulldata)
		J <- length(K)
		nfact <- ncol(theta0)
		prior.t.var <- diag(nfact)
		locz <- 1		
		P0 <- P1 <- matrix(0,N,J)		
        if(nfact > 1)		
		  theta1 <- theta0 + rmvnorm(N,rep(0,nfact), diag(rep(sqrt(cand.t.var),nfact))) 
        else
          theta1 <- theta0 + rnorm(N,0,sqrt(cand.t.var))		
		for(i in 1:J){
			if(K[i]==2){
				tmp <- P.mirt(lambdas[i,],zetas[locz],theta0,guess[i])
				tmp[tmp < 1e-7] <- 1e-7	
				P0[,i] <- ifelse(fulldata[,itemloc[i]],tmp,1-tmp)
				tmp <- P.mirt(lambdas[i,],zetas[locz],theta1,guess[i])
				tmp[tmp < 1e-7] <- 1e-7	
				P1[,i] <- ifelse(fulldata[,itemloc[i]],tmp,1-tmp)						
				locz <- locz + 1				
			} else {
				upz <- locz + (K[i]-2)
				tmp <- P.poly(lambdas[i,],zetas[locz:upz],theta0,itemexp=TRUE)
				tmp[tmp < 1e-7] <- 1e-7	
				P0[,i] <- rowSums(tmp * fulldata[,itemloc[i]:(itemloc[i+1]-1)])
				tmp <- P.poly(lambdas[i,],zetas[locz:upz],theta1,itemexp=TRUE)
				tmp[tmp < 1e-7] <- 1e-7	
				P1[,i] <- rowSums(tmp * fulldata[,itemloc[i]:(itemloc[i+1]-1)])							
				locz <- locz + (K[i]-1)
			}
		}			
		irt0 <- rowSums(log(P0)) + dmvnorm(theta0,rep(0,nfact),prior.t.var,log=TRUE)		
		irt1 <- rowSums(log(P1)) + dmvnorm(theta1,rep(0,nfact),prior.t.var,log=TRUE)				
		accept <- irt1 - irt0
		accept <- ifelse(accept>0,0,accept)
		accept <- ifelse(runif(N) < exp(accept),TRUE,FALSE) 		
		theta1[!accept,] <- theta0[!accept,]	
		attr(theta1, "Proportion Accepted") <- sum(accept)/N 		
		return(theta1) 
	}
	dpars.dich <- function(lambda,zeta,g,dat,Thetas,estGuess) {
		nfact <- length(lambda)
		P <- P.mirt(lambda, zeta, Thetas, g)						
		if(estGuess){
			Pstar <- P.mirt(lambda, zeta, Thetas, 0)
			Q <- 1 - P		
			PQ <- P*Q	
			L1 <- sum((dat-P)/(P-g) * (Pstar/P))
			L2 <- sum((dat-P) * (Pstar/P))
			L3 <- colSums((dat-P) * Thetas * (Pstar/P))			
			dL <- c(L1,L2,L3)
			d2L <- matrix(0,nfact+2,nfact+2)			
			d2L[1,1] <- -sum(Q/(1-g) * 1/(P-g) *(Pstar/P))			
			d2L[2,2] <- -sum(PQ*(Pstar/P))			
			d2L[1,2] <- d2L[2,1] <- -sum(Q/(1-g)*(Pstar/P))	
			const <- PQ * (Pstar/P)^2
			d2L[3:(3+nfact-1),3:(3+nfact-1)] <- (-1)* .Call("dichOuter",Thetas,
				const,nfact,nrow(Thetas))
			d2L[1,3:(3+nfact-1)] <- d2L[3:(3+nfact-1),1] <- -colSums(Thetas*(Q/(1-g)*(Pstar/P)))	
			d2L[2,3:(3+nfact-1)] <- d2L[3:(3+nfact-1),2] <-	-colSums(Thetas*PQ*(Pstar/P))			
		} else {
			PQ <- P*(1-P)
			L1 <- sum(dat-P)
			L2 <- colSums((dat-P)*Thetas)
			dL <- c(L1,L2)		
			d2L <- matrix(0,nfact+1,nfact+1)						
			L11 <- .Call("dichOuter",Thetas,PQ,nfact,nrow(Thetas))
			if(nfact > 1) d2L[1:nfact+1,1:nfact+1] <- -L11
				 else d2L[nfact+1,nfact+1] <- -L11 				
			d2L[1,1]<- (-1)*sum(PQ)		
			d2L[1,1:nfact+1] <- d2L[1:nfact+1,1] <- (-1)*colSums(PQ * Thetas)
		}	
		list(grad = dL, hess = d2L)
	} 
	dpars.poly <- function(lambda,zeta,dat,Thetas){  
		nzeta <- length(zeta)
		ncat <- nzeta + 1		
		nfact <- length(lambda)
		factind <- ncat:(ncat+nfact-1)
		N <- nrow(Thetas)		
		P <- P.poly(lambda,zeta,Thetas)	
		dL <- rep(0,nzeta + nfact)
		d2L <- matrix(0,nzeta + nfact, nzeta + nfact)	
		PQfull <- P * (1 - P)			 
		for(i in 1:ncat){
			if(i < ncat){
				Pk_1 <- P[,i]	
				Pk <- P[,i+1]
				Pk.1 <- P[,i+2]
				PQ_1 <- PQfull[,i]			
				PQ <- PQfull[,i+1]
				PQ.1 <- PQfull[,i+2]		
				dif1 <- dat[,i]/(Pk_1 - Pk)
				dif1sq <- dat[,i]/(Pk_1 - Pk)^2		
				dif2 <- dat[,i+1]/(Pk - Pk.1)	  
				dif2sq <- dat[,i+1]/(Pk - Pk.1)^2		
				
				dL[i] <- sum(-1 * PQ * (dif1 - dif2)) 
				d2L[i,i] <- sum(-1 * (PQ^2) * (dif1sq + dif2sq) -
					(dif1 - dif2) * (Pk * (1 - Pk) * (1 - 2*Pk)))
				if(i < nzeta) d2L[i,i+1] <- d2L[i+1,i] <- sum(dif2sq * PQ.1 * PQ)        		
				d2L[factind,i] <- d2L[i,factind] <- 
					colSums(-(dif2sq * PQ * (PQ - PQ.1) * Thetas) +
					(dif1sq * PQ * (PQ_1 - PQ) * Thetas) - 
					((dif1 - dif2) * (Pk * (1 - Pk) * (1 - 2*Pk) * Thetas)))
			}		
			Pk_1 <- P[,i]	
			Pk <- P[,i+1]		
			PQ_1 <- PQfull[,i]			
			PQ <- PQfull[,i+1]		
			dif1 <- dat[,i]/(Pk_1 - Pk)
			dif1sq <- dat[,i]/(Pk_1 - Pk)^2
			dL[factind] <- dL[factind] + colSums(dif1 * (PQ_1 - PQ) * Thetas)			
			d2L[factind,factind] <- d2L[factind,factind] + .Call("polyOuter",Thetas,Pk,Pk_1,PQ_1,PQ,
				dif1sq,dif1,nfact,nrow(Thetas))					
		} 
		return(list(grad=dL, hess=d2L))	
	}
	
	############################
	
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
			zetas[i] <- qnorm(mean(fulldata[,itemloc[i]]))/cs[i]
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
	for(i in 1:20) theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)
	for(i in 1:20){
		theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)
		if(attr(theta0,"Proportion Accepted") > .5 && nfact < 5) cand.t.var <- cand.t.var + .1 
		else if(attr(theta0,"Proportion Accepted") > .3) cand.t.var <- cand.t.var + .1
     	if(attr(theta0,"Proportion Accepted") < .2)	 cand.t.var <- cand.t.var - .1
	}	
	m.thetas <- list()		
	SEM.stores <- matrix(0,SEM.cycles,npars)
	phi <- rep(0,npars)
	Tau <- info <- h <- matrix(0,npars,npars)
	m.list <- list()	  
	conv <- 0
	gamma <- k <- 1	
	startvalues <- pars	
	
	for(cycles in 1:ncycles)
	{ 
		if(cycles == (SEM.cycles + 1)){
		    pars <- rep(0,npars)
			for(i in 1:SEM.cycles) pars <- pars + SEM.stores[i,]
			pars <- pars/SEM.cycles	
			k <- kdraws	
		}
		if(cycles > SEM.cycles)
			gamma <- 1/(cycles - SEM.cycles)        					
		
		lambdas <- matrix(pars[lamind],ncol=nfact,byrow=TRUE)
		zetas <- pars[zetaind]
		guess <- rep(0,J)
		guess[estGuess] <- pars[gind]
		
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
		if(cycles <= SEM.cycles){
		    correction <- solve(ave.h) %*% grad
			correction[correction > .5] <- .5
			correction[correction < -.5] <- -.5
			if(any(estGuess)){
				correction[correction[gind] > .05] <- .05
				correction[correction[gind] < -.05] <- -.05
			}	
			SEM.stores[cycles,] <- pars <- pars + correction				
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
		if(all(gamma*correction < tol)) conv <- conv + 1
			else conv <- 0		
		if(conv == 3) break		
		pars <- pars + gamma*correction
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE 	
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	}
	SE <- diag(solve(info))
	if(any(SE < 0)){
		warning("Solution is not proper, information matrix is not positive definite.\n")
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
	if (sum(F[ ,1] < 0)) F[ ,1] <- (-1)*F[ ,1]  
	h2 <- rowSums(F^2) 
	
	mod <- list(pars=pars, guess=guess, SEpars=SEpars, cycles=cycles - SEM.cycles,
		Theta=theta0, fulldata=fulldata,K=K, F=F, h2=h2, fulldata=fulldata, 
		itemloc=itemloc, converge = converge,Call=Call)	 
	class(mod) <- 'polymirt'
	mod	
}
