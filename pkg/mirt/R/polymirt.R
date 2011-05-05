polymirt <- function(data, nfact, guess = 0, prev.cor = NULL, 
	ncycles = 2000, SEM.cycles = 30, kdraws = 5, tol = .001, debug = FALSE){
	
	draw.thetas <- function(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var) {    
		N <- nrow(fulldata)
		J <- length(K)
		nfact <- ncol(theta0)
		prior.t.var <- diag(nfact)
		locz <- locl <- 1
		if(nfact > 1) locl <- 1:(nfact - 1)		
		P0 <- P1 <- matrix(0,N,J)		
        if(nfact > 1)		
		  theta1 <- theta0 + rmvnorm(N,rep(0,nfact), diag(rep(sqrt(cand.t.var),nfact))) 
        else
          theta1 <- theta0 + rnorm(N,0,sqrt(cand.t.var))
		for(i in 1:J){
			if(K[i]==2){
				tmp <- P.mirt(lambdas[locl],zetas[locz],theta0,guess[i])
				tmp[tmp < 1e-7] <- 1e-7	
				P0[,i] <- ifelse(fulldata[,itemloc[i]],tmp,1-tmp)
				tmp <- P.mirt(lambdas[locl],zetas[locz],theta1,guess[i])
				tmp[tmp < 1e-7] <- 1e-7	
				P1[,i] <- ifelse(fulldata[,itemloc[i]],tmp,1-tmp)		
				locl <- locl + nfact
				locz <- locz + 1				
			} else {
				upz <- locz + (K[i]-2)
				tmp <- P.poly(lambdas[locl],zetas[locz:upz],theta0,itemexp=TRUE)
				tmp[tmp < 1e-7] <- 1e-7	
				P0[,i] <- rowSums(tmp * fulldata[,itemloc[i]:(itemloc[i+1]-1)])
				tmp <- P.poly(lambdas[locl],zetas[locz:upz],theta1,itemexp=TRUE)
				tmp[tmp < 1e-7] <- 1e-7	
				P1[,i] <- rowSums(tmp * fulldata[,itemloc[i]:(itemloc[i+1]-1)])			
				locl <- locl + nfact
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
	dpars.dich <- function(lambda,zeta,dat,Thetas) {
		nfact <- length(lambda)
		P <- P.mirt(lambda, zeta, Thetas, guess)								
		PQ <- P*(1-P)		
		L2 <- colSums((dat-P)*Thetas)						
		L1 <- sum(dat-P)
        dL <- c(L1,L2) 		
		d2L <- matrix(0,nfact+1,nfact+1)		
		L11 <- matrix(0,nfact,nfact)		
        if(nfact > 1){		
			for(i in 1:N){
			  temp <- outer(Thetas[i,], Thetas[i,])
			  L11 <- L11 + temp * PQ[i] 
			}  
			d2L[1:nfact+1,1:nfact+1] <- -L11
		} else d2L[nfact+1,nfact+1] <- (-1)*colSums(PQ * Thetas^2)
		d2L[1,1]<- (-1)*sum(PQ)		
		d2L[1,1:nfact+1] <- d2L[1:nfact+1,1] <- (-1)*colSums(PQ * Thetas)	
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
			if(nfact == 1)
				d2L[ncat,ncat] <- d2L[ncat,ncat] + sum(-1*(dif1sq * ((PQ_1 - PQ)^2) * Thetas^2) +  
					(dif1 * (Pk_1*(1-Pk_1)*(1-2*Pk_1) - Pk*(1-Pk)*(1-2*Pk)) * Thetas^2 ))		
			else {
				L11 <- matrix(0,nfact,nfact)
				for(j in 1:N){
					temp <- outer(Thetas[j,], Thetas[j,])
					L11 <- L11 + (-1*(dif1sq[j] * (PQ_1[j] * Thetas[j,] - PQ[j] * Thetas[j,]) %*% 
					t(PQ_1[j] * Thetas[j,] - PQ[j] * Thetas[j,])) + (dif1[j] * (Pk_1[j]*(1-Pk_1[j])*(1-2*Pk_1[j]) 
					* temp - Pk[j]*(1-Pk[j])*(1-2*Pk[j])* temp))) 
				}  
				d2L[factind,factind] <- L11 + d2L[factind,factind]
			}	
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
	uniques <- list()
	for(i in 1:J)
		uniques[[i]] <- sort(unique(data[,i]))
	K <- rep(0,J)
	for(i in 1:J) K[i] <- length(uniques[[i]])    
	itemloc <- cumsum(c(1,K))
	index <- 1:J	
	fulldata <- matrix(0,N,sum(K))
	Names <- NULL
	for(i in 1:J)
        Names <- c(Names, paste("Item.",i,"_",1:K[i],sep=""))				
	colnames(fulldata) <- Names			
	for(i in 1:J){
		if(setequal(uniques[[i]], c(0,1)))
			fulldata[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- cbind(data[,ind],abs(1-data[,ind]))
		ind <- index[i]
		dummy <- matrix(0,N,K[ind])
		for (j in 0:(K[ind]-1))  
			dummy[,j+1] <- as.integer(data[,ind] == uniques[[ind]][j+1])  		
		fulldata[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- dummy		
	}			
	if(!is.null(prev.cor)) {
		if (ncol(prev.cor) == nrow(prev.cor)) Rpoly <- prev.cor
			else stop("Correlation matrix is not square.\n")
	} 	else Rpoly <- abs(cor(data))^(1/1.15) * sign(cor(data))	  			
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
			temp <- colMeans(fulldata[,itemloc[i]:(itemloc[i+1]-2)])
			temp <- cumsum(temp)			
			zetas[loc:(loc+K[i]-2)] <- qnorm(1 - temp)/cs[i]	
			loc <- loc + K[i] - 1	
		}		
	}	
	npars <- length(c(lambdas,zetas))
	parind <- 1:npars
	pars <- rep(NA,npars)
	Ksum <- cumsum(K + nfact - 1) - (nfact-1)
	for(i in 1:J)
		pars[Ksum[i]:(Ksum[i] + nfact - 1)] <- lambdas[i,]
	lamind <- parind[!is.na(pars)]
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
		}
		if(cycles > SEM.cycles){
			gamma <- 1/(cycles - SEM.cycles)        		
			k <- kdraws
		}		
		lambdas <- matrix(pars[lamind],ncol=nfact)
		zetas <- pars[zetaind]
		
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
					temp <- dpars.dich(lambdas[i+1,],zetas[loc],fulldata[,itemloc[i+1]],
						m.thetas[[j]])
					ind <- parind[is.na(g)][1]
					ind2 <- ind+nfact		
					g[ind:ind2] <- temp$grad
					h[ind:ind2,ind:ind2] <- temp$hess
					loc <- loc + 1
				} else {
					loc2 <- loc + K[i+1] - 2
					temp <- dpars.poly(lambdas[i+1,],zetas[loc:loc2],fulldata[,itemloc[i+1]:(itemloc[i+2]-1)],
						m.thetas[[j]])
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
			SEM.stores[cycles,] <- pars <- pars + correction				
			next
		}	
		
		#Step 3. Update R-M step
		Tau <- Tau + gamma*(ave.h - Tau)		
		correction <- (solve(Tau) %*% grad)	
		correction[correction > .5] <- .5
		correction[correction < -.5] <- -.5	
		if(all(gamma*correction < tol)) conv <- conv + 1
			else conv <- 0		
		if(conv == 3) break		
		pars <- pars + gamma*correction
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE 	
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	}
	SE <- sqrt(diag(solve(info)))	
	lambdas <- matrix(pars[lamind],ncol=nfact)
	SElam <- SE[lamind]
	zetas <- pars[zetaind]	
	SEzeta <- SE[zetaind]
	
	list(lambdas = lambdas, zetas = zetas, SElam = SElam, SEzeta = SEzeta, 
		cycles = cycles - SEM.cycles, Theta = theta0, K=K, data = data, fulldata=fulldata)    
}

