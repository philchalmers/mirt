# theta combinations
thetaComb <- function(theta, nfact)
{
	if (nfact == 1) Theta <- matrix(theta)
		else if (nfact == 2) Theta <- expand.grid(theta,theta)   
		else if (nfact == 3) Theta <- expand.grid(theta,theta,theta)  
		else if (nfact == 4) Theta <- expand.grid(theta,theta,theta,theta)
		else if (nfact == 5) Theta <- expand.grid(theta,theta,theta,theta,theta)        
		else if (nfact == 6) Theta <- expand.grid(theta,theta,theta,theta,theta,theta)        
		else if (nfact == 7) Theta <- expand.grid(theta,theta,theta,theta,theta,theta,theta)        
		else if (nfact == 8) Theta <- expand.grid(theta,theta,theta,theta,theta,theta,theta,theta) 
	Theta <- as.matrix(Theta)	
	return(Theta)     
}

# start values
start.values <- function(fulldata, guess, Rpoly, nfact=2, bfactor = FALSE, nowarn = TRUE)
{
	if (length(guess) == 1) guess <- rep(guess,ncol(fulldata))
		else if (length(guess) > ncol(fulldata) || length(guess) < ncol(fulldata)) 
			stop("The number of guessing parameters is incorrect.")  	
	if (bfactor){ 
		FA <- fa(Rpoly,1, warnings= !nowarn)
		loads <- unclass(FA$load)
		cs <- sqrt(abs(FA$u))      
		dstart <- qnorm(colMeans(fulldata))/cs
		astart <- loads/cs
		startvalues <- cbind(astart,astart/2,dstart)
	} else {    
		FA <- fa(Rpoly,nfact,rotate = 'none', warnings= !nowarn)	
		loads <- unclass(loadings(FA))
		u <- FA$unique
		u[u < .001 ] <- .2
		cs <- sqrt(u)
		dstart <- qnorm(colMeans(fulldata))/cs
		astart <- loads/cs
		startvalues <- cbind(astart,dstart)
	}  	
	startvalues
}

# Rotation function
Rotate <- function(F, rotate)
{
	orthogonal <- c("varimax", "quartimax", "tandemI", "tandemII", "entropy", "mccammon")
	oblique <- c("promax", "oblimin", "quartimin", "oblimax", "simplimax")
	if (!any(rotate %in% c(orthogonal,oblique))) stop("Unknown rotation specified.")
	if(any(rotate %in% orthogonal)){
		oblique <- FALSE
		rotF <- GPForth(F, method = rotate)
	}
	if(any(rotate %in% oblique)){
		oblique <- TRUE
		if(rotate == 'promax') rotF <- Promax(F) 
			else rotF <- GPFoblq(F, method = rotate)
	}
	attr(rotF,"oblique") <- oblique 
	return(rotF)
}  

# MAP scoring for mirt
MAP.mirt <- function(Theta,a,d,guess,patdata)
{
	Theta <- t(as.matrix(Theta))
	L <- 0
	for (j in 1:length(patdata)){
		if(patdata[j] == 1) L <- log(P.mirt(a[j, ],d[j],Theta,guess[j])) + L
			else L <- log(1 - P.mirt(a[j, ],d[j],Theta,guess[j])) + L	
	}
	mu <- 0
	sigma <- 1
	L <- (-1)*(L + sum(log(exp(-0.5*((Theta - mu)/sigma)^2))))
	L  
}  

# MAP scoring for bfactor
MAP.bfactor <- function(Theta,a,d,guess,patdata,logicalfact)
{
	Theta <- t(as.matrix(Theta))
	L <- 0
	for (j in 1:length(patdata)){
		if(patdata[j] == 1) L <- log(P.bfactor(a[j, ],d[j],Theta,guess[j],logicalfact[j, ])) + L
			else L <- log(1 - P.bfactor(a[j, ],d[j],Theta,guess[j],logicalfact[j, ])) + L	
	}
	mu <- 0
	sigma <- 1
	L <- (-1)*(L + sum(log(exp(-0.5*((Theta - mu)/sigma)^2))))
	L  
}  

#trace lines for polymirt
P.poly <- function(lambda, zetas, Thetas, itemexp = FALSE)
{	
	ncat <- length(zetas) + 1
	nfact <- length(lambda)
	Pk <- matrix(0,nrow(Thetas),ncat+1)
	Pk[,1] <- 1	
	for(i in 1:(ncat-1))			
		Pk[ ,i+1] <- P.mirt(lambda,zetas[i],Thetas,0)		
	if(itemexp){
		P <- matrix(0,nrow(Thetas),ncat)		
		for(i in ncat:1)
			P[,i] <- Pk[,i] - Pk[,i+1]						
		Pk <- P
	}	
	return(Pk)
}

# Trace lines for mirt models
P.mirt <- function(a, d, Theta, g)
{ 
	nfact <- length(a)
	nquad <- nrow(Theta)
	traces <- .Call("traceLinePts",
					as.double(a), 
					as.double(d),
					as.double(g),  
					as.double(as.matrix(Theta)), 
					as.integer(nquad), 
					as.integer(nfact))
	return(traces)
}

P.comp <- function(a,d,thetas,c = 0){
	nfact <- length(a)
	P <- rep(1,nrow(thetas))
	for(i in 1:nfact)
		P <- P * P.mirt(a[i],d[i],matrix(thetas[,i]),0)
	P <- c + (1-c) * P
	P	
} 

# Estep
Estep.mirt <- function(pars, tabdata, Theta, prior, guess) 
{
	a <- as.matrix(pars[ ,1:(ncol(pars) - 1)])
	nfact <- ncol(a)
	nitems <- nrow(a)
	nquad <- nrow(Theta)
	d <- pars[ ,ncol(pars)]    
	r <- tabdata[ ,ncol(tabdata)]
	X <- tabdata[ ,1:(ncol(tabdata) - 1)]     

	itemtrace <- r1 <- r0 <- matrix(0,nrow=nitems,ncol=nrow(Theta))
	for (i in 1:nitems) itemtrace[i, ] <- 
	  P.mirt(a[i, ],d[i],Theta,guess[i])    
	  
	retlist <- .Call("Estep",                     	
					 as.double(itemtrace),
					 as.double(prior),
					 as.integer(X), 
					 as.integer(nfact),      
					 as.integer(r))   

	N <- retlist$r1 + retlist$r0    
	empprior <- colSums(N)/sum(N)    
	rlist <- list(retlist$r1, N, retlist$expected, empprior)
	return(rlist)
} 

P.bfactor <- function(a, d, Theta, g, patload)
{ 
	a <- a[patload]
	nfact <- length(a)
	nquad <- nrow(Theta)
	traces <- .Call("traceLinePts",                    
					as.double(a), 
					as.double(d),
					as.double(g),  
					as.double(as.matrix(Theta)), 
					as.integer(nquad), 
					as.integer(nfact))
	return(traces)
}

# Estep
Estep.bfactor <- function(pars, tabdata, Theta, prior, guess, logicalfact, specific, sitems) 
{
	a <- as.matrix(pars[ ,1:(ncol(pars) - 1)])
	nfact <- ncol(a)
	nitems <- nrow(a)
	nquad <- nrow(Theta)	
	d <- pars[ ,ncol(pars)]    
	r <- tabdata[ ,ncol(tabdata)]
	X <- tabdata[ ,1:(ncol(tabdata) - 1)]
	as <- rowSums(pars[,2:nfact])

	itemtrace <- r1 <- r0 <- matrix(0,nrow=nitems,ncol=nrow(Theta))
	for (i in 1:nitems) itemtrace[i, ] <- 
	  P.bfactor(a[i, ],d[i],Theta,guess[i],logicalfact[i,])

	retlist <- .Call("Estepbfactor",
					as.double(itemtrace),
					as.double(prior), 					
					as.integer(X), 
					as.integer(nfact),
					as.integer(r),
					as.double(sitems))      

	r1 <- N <- matrix(0, nitems, nrow(Theta))
	for (i in 1:nitems){
	  r1[i, ] <- retlist$r1[(specific[i] - 1)*nitems + i, ]
	  N[i, ] <- retlist$r0[(specific[i] - 1)*nitems + i, ] + r1[i, ]
	}
		
	rlist <- list(r1, N, retlist$expected)
	return(rlist)
}      

draw.thetas <- function(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var,
	prior.t.var = diag(ncol(theta0)), prior.mu = rep(0,ncol(theta0)), estComp = rep(FALSE,length(K))) 
{ 			
	N <- nrow(fulldata)
	J <- length(K)
	nfact <- ncol(theta0)				
	P0 <- P1 <- matrix(0,N,J)		
	unif <- runif(N)
	if(nfact > 1)		
		theta1 <- theta0 + rmvnorm(N,prior.mu, 
			diag(rep(cand.t.var,ncol(theta0)))) 
	else
		theta1 <- theta0 + rnorm(N,prior.mu,sqrt(cand.t.var))							
	den0 <- dmvnorm(theta0,prior.mu,prior.t.var)
	den1 <- dmvnorm(theta1,prior.mu,prior.t.var)		
	accept <- .Call("drawThetas",
					as.numeric(unif),
					as.numeric(den0),
					as.numeric(den1),
					as.numeric(lambdas),
					as.numeric(zetas),
					as.numeric(guess),
					as.numeric(theta0),
					as.numeric(theta1),
					as.integer(fulldata),
					as.integer(itemloc - 1),
					as.integer(K),
					as.integer(J),
					as.integer(N),
					as.integer(ncol(lambdas)),
					as.integer(estComp))
	log.lik <- accept[N+1]			
	accept <- as.logical(accept[-(N+1)])				
	theta1[!accept,] <- theta0[!accept,]	
	attr(theta1, "Proportion Accepted") <- sum(accept)/N 				
	attr(theta1, "log.lik") <- log.lik		
	return(theta1) 
}	

d.group <- function(grouplist,theta)
{		
	tr <- function(x) sum(diag(x))
	x <- theta
	u <- grouplist$u	
	sig <- grouplist$sig
	N <- nrow(x)
	nfact <- length(u)
	selcov <- matrix(FALSE,nfact,nfact)
	selcov <- lower.tri(selcov) 
	diag(selcov) <- TRUE
	npars <- length(sig) + nfact	
	g <- rep(0,nfact + nfact*(nfact+1)/2)	
	invSig <- solve(sig)	
	Z <- t(x-u) %*% (x-u)
	g[1:nfact] <- N * invSig %*% (colMeans(x) - u) 		
	tmp <- .5 * invSig %*% (Z - N * sig) %*% invSig  
	g[(nfact+1):length(g)] <- tmp[selcov]
	h <- matrix(0,npars,npars)
	sel <- 1:npars		
	cMeans <- N*(colMeans(x) - u)
	Zdif <- (Z - N * sig)		
	h <- .Call("dgroup",
				as.numeric(sig),
				as.numeric(invSig),
				as.numeric(cMeans),					
				as.numeric(Z),
				as.numeric(Zdif),
				as.integer(N),
				as.integer(nfact),
				as.integer(npars))				
	sel <- sel[c(rep(TRUE,nfact),as.logical(selcov))]	
	h <- h[sel,sel] 
	list(h=h,g=g) 
}

dpars.dich <- function(lambda,zeta,g,dat,Thetas,estGuess)
{
	nfact <- length(lambda)
	P <- P.mirt(lambda, zeta, Thetas, g)						
	if(estGuess){ 
		r <- dat
		f <- 1		
		c <- g
		a <- lambda
		d <- zeta			
		thetas <- Thetas
		Pstar <- P.mirt(lambda,zeta,Thetas,0)		
		Qstar <- 1 - Pstar
		Q <- 1 - P
		da <- rep(0,nfact)	
		dd <- sum((1-g)*Pstar*Qstar*(r/P - (f-r)/Q))
		dc <- sum(Qstar*(r/P - (f-r)/Q))
		for(i in 1:nfact){
			da[i] <- sum(Thetas[,i]*Pstar*Qstar*(1-g)*(r/P - (f-r)/Q))
		}
		dL <- c(dd,da,dc)
				
		hess <- matrix(0,nfact + 2,nfact + 2)	
		aNames <- paste("a",1:nfact,sep='_')
		Names <- c('d',paste("a",1:nfact,sep='_'),'c')
		colnames(hess) <- rownames(hess) <- Names
		hsize <- nfact+2
		const1 <- (r/P - (f-r)/Q)*(Qstar-Pstar)
		const2 <- (r/P^2 + (f-r)/Q^2)	
		hess[1,1] <- sum((1-c)*Pstar*Qstar*(const1 - 
			Pstar*Qstar*(1-c)*const2))		
		hess[hsize,hsize] <- -sum(Qstar^2 *(r/P^2 + (f-r)/Q^2))
		hess[hsize,1] <- hess[1,hsize] <- sum(-Pstar*Qstar*((r/P - (f-r)/Q) + Qstar*(1-c)*const2))  				
		for(i in 1:nfact){
			hess[1,1+i] <- hess[1+i,1] <- sum((1-c)*thetas[,i]*Pstar*Qstar*(const1 - 
				Pstar*Qstar*(1-c)*const2))			
			hess[hsize,1+i] <- hess[1+i,hsize] <- sum(-thetas[,i]*Pstar*Qstar*((r/P - (f-r)/Q) + Qstar*(1-c)*const2))		
			for(j in 1:nfact){
				if(i == j)
					hess[1+i,1+i] <- sum(thetas[,i]^2 *Pstar*Qstar*(1-c)*(const1 - 
						(1-c)*Pstar*Qstar*const2))
				if(i < j)
					hess[1+i,1+j] <- hess[1+j,1+i] <- sum(thetas[,i]*thetas[,j] *Pstar*Qstar*(1-c)*
						(const1 - (1-c)*Pstar*Qstar*const2))					
			}
		}	
		d2L <- hess			
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

dpars.comp <- function(lambda,zeta,g,dat,Thetas,estg = FALSE)
{	
	nfact <- length(lambda)	
	pars <- c(zeta,lambda,g)
	if(estg){
		grad <- function(pars, r, thetas){
			f <- 1
			nfact <- ncol(thetas)
			d <- pars[1:nfact]	
			a <- pars[(nfact+1):(length(pars)-1)]
			c <- pars[length(pars)]
			P <- P.comp(a,d,thetas,c)		
			Pstar <- P.comp(a,d,thetas,0)		
			Qstar <- 1 - Pstar
			Q <- 1 - P
			const1 <- (r/P - (f-r)/Q)
			dd <- da <- rep(0,nfact)		
			dc <- sum(Qstar*const1)
			for(i in 1:nfact){
				Pk <- P.mirt(a[i],d[i],matrix(thetas[,i]),0)
				Qk <- 1 - Pk
				dd[i] <- sum((1-c)*Pstar*Qk*const1)
				da[i] <- sum((1-c)*Pstar*Qk*thetas[,i]*const1)
			}
			return(c(dd,da,dc))
		}		
		hess <- function(pars, r, thetas){
			f <- 1
			nfact <- ncol(thetas)
			d <- pars[1:nfact]	
			a <- pars[(nfact+1):(length(pars)-1)]
			c <- pars[length(pars)]
			P <- P.comp(a,d,thetas,c)		
			Pstar <- P.comp(a,d,thetas,0)		
			Qstar <- 1 - Pstar
			Q <- 1 - P	
			const1 <- (r/P - (f-r)/Q)
			const2 <- (r/P^2 + (f-r)/Q^2)	
			hess <- matrix(0,nfact*2+1,nfact*2+1)
			dNames <- paste("d",1:nfact,sep='_')
			aNames <- paste("a",1:nfact,sep='_')
			Names <- c(paste("d",1:nfact,sep='_'),paste("a",1:nfact,sep='_'),'c_0')
			for(i in 1:(nfact*2+1)){		
				for(j in 1:(nfact*2+1)){
					if(i <= j){
						d1 <- strsplit(Names[c(i,j)],"_")[[1]]
						d2 <- strsplit(Names[c(i,j)],"_")[[2]]
						k <- as.numeric(d1[2])
						m <- as.numeric(d2[2])
						Pk <- P.mirt(a[k],d[k],matrix(thetas[,k]),0)
						Qk <- 1 - Pk	
						Pm <- P.mirt(a[m],d[m],matrix(thetas[,m]),0)
						Qm <- 1 - Pm									
						if(i == j && d1[1] == 'd'){
							hess[i,i] <- sum((1-c)*Pstar*Qk*(const1*((1-c)*Qk - Pk) - Pstar*Qk*(1-c)*const2))
							next
						}
						if(i == j && d1[1] == 'a'){
							hess[i,i] <- sum((1-c)*thetas[,k]^2*Pstar*Qk*(const1*((1-c)*Qk - Pk) - Pstar*Qk*(1-c)*const2))
							next		
						}
						if(i == j && d1[1] == 'c'){
							hess[i,i] <- -sum(Qstar^2 * const2)
							next		
						}	
						if(d1[1] == 'a' && d2[1] == 'a'){
							hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*thetas[,m]*Qk*Pstar*Qm*(const1 - Pstar*(1-c)*const2))
							next
						}
						if(d1[1] == 'd' && d2[1] == 'd'){
							hess[i,j] <- hess[j,i] <- sum((1-c)*Qk*Pstar*Qm*(const1 - Pstar*(1-c)*const2))
							next
						}
						if(d1[1] == 'a' && d2[1] == 'c'){
							hess[i,j] <- hess[j,i] <- -sum(thetas[,k]*Pstar*Qk*(const1 + Qstar*(1-c)*const2))
							next
						}
						if(d1[1] == 'd' && d2[1] == 'c'){
							hess[i,j] <- hess[j,i] <- -sum(Pstar*Qk*(const1 + Qstar*(1-c)*const2))
							next
						}
						if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
							hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*Pstar*Qk*(const1*((1-c)*Qk - Pk) - Pstar*Qk*(1-c)*const2))
							next	
						}
						if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
							hess[i,j] <- hess[j,i] <- sum((1-c)*Qk*thetas[,m]*Pstar*Qm*(const1 - Pstar*(1-c)*const2))
							next
						}						
					}
				}
			}	
			return(hess)
		}		
		return(list(grad = grad(pars, dat, Thetas), hess = hess(pars, dat, Thetas)))
	} else {			
		P <- P.comp(lambda,zeta,Thetas)	
		Q <- 1 - P	
		da <- dd <- rep(0,nfact)	
		for(i in 1:nfact){
			Pk <- P.mirt(lambda[i],zeta[i],matrix(Thetas[,i]),0)
			Qk <- 1 - Pk
			const <- (1 - dat)*P/Q
			dd[i] <- sum(Qk*(dat - const))
			da[i] <- sum(Thetas[,i]*Qk*(dat - const))
		}
		hess <- matrix(0,nfact*2,nfact*2)
		dNames <- paste("d",1:nfact,sep='_')
		aNames <- paste("a",1:nfact,sep='_')
		Names <- c(paste("d",1:nfact,sep='_'),paste("a",1:nfact,sep='_'))
		f <- 1
		r <- dat
		for(i in 1:(nfact*2)){		
			for(j in 1:(nfact*2)){
				if(i <= j){
					d1 <- strsplit(Names[c(i,j)],"_")[[1]]
					d2 <- strsplit(Names[c(i,j)],"_")[[2]]
					k <- as.numeric(d1[2])
					m <- as.numeric(d2[2])
					Pk <- P.mirt(lambda[k],zeta[k],matrix(Thetas[,k]),0)
					Qk <- 1 - Pk	
					Pm <- P.mirt(lambda[m],zeta[m],matrix(Thetas[,m]),0)
					Qm <- 1 - Pm									
					if(i == j && d1[1] == 'd'){
						hess[k,k] <- sum(-Pk*Qk*(r - (f-r)*P/Q) - Qk^2 * (f-r)*P/Q^2)
						next
					}
					if(i == j && d1[1] == 'a'){
						hess[k+nfact,k+nfact] <- sum(Thetas[,k]^2 *
							(-Pk*Qk*(r - (f-r)*P/Q) - Qk^2 * (f-r)*P/Q^2))
						next		
					}				
					if(d1[1] == 'a' && d2[1] == 'a'){
						hess[i,j] <- hess[j,i] <- -sum(Thetas[,k]*Thetas[,m]*Qk*Qm*(f-r)*P/Q^2) 
						next
					}
					if(d1[1] == 'd' && d2[1] == 'd'){
						hess[i,j] <- hess[j,i] <- -sum(Qk*Qm*(f-r)*P/Q^2)
						next
					}	
					if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
						hess[i,j] <- hess[j,i] <- sum(Thetas[,k]*Qk*(-Pk*(r - (f-r)*P/Q) - 
							Qk*(f-r)*P/Q^2))
						next	
					}
					if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
						hess[i,j] <- hess[j,i] <- -sum(Qk*Qm*Thetas[,m]*(f-r)*P/Q^2)
						next
					}						
				}
			}
		}		
	}	
	return(list(grad = c(dd,da), hess = hess))
}

dpars.poly <- function(lambda,zeta,dat,Thetas)
{  
	nzeta <- length(zeta)			
	nfact <- length(lambda)				
	N <- nrow(Thetas)		
	P <- P.poly(lambda,zeta,Thetas)			
	ret <- .Call("dparsPoly",
				as.numeric(P), 
				as.numeric(Thetas), 
				as.integer(dat),
				as.integer(nzeta),
				as.integer(nfact),
				as.integer(N)) 				 
	return(ret)	
}



gamma.cor <- function(x)
{ 
	concordant <- function(x){ 	  
			mat.lr <- function(r, c){ 
				lr <- x[(r.x > r) & (c.x > c)] 
				sum(lr) 
			} 	  
		r.x <- row(x) 
		c.x <- col(x) 	  
		sum(x * mapply(mat.lr, r = r.x, c = c.x)) 
	} 	
	discordant <- function(x){ 	  
		mat.ll <- function(r, c){ 
			ll <- x[(r.x > r) & (c.x < c)] 
			sum(ll) 
		} 	  
		r.x <- row(x) 
		c.x <- col(x) 	  
		sum(x * mapply(mat.ll, r = r.x, c = c.x)) 
	} 
	c <- concordant(x) 
	d <- discordant(x) 
	gamma <- (c - d) / (c + d) 
	gamma 
} 

betaprior <- function(a,b,g,W=20){
	a <- a + (1-g)*W
	b <- b + g*W
	grad <- ((a-1) * g^(a-1) * (1-g)^(b-1) - (b-1)*g^(a-1)*(1-g)^(b-1))/ 
		(g^(a-1) * (1-g)^(b-1))
	hess <- -((g^(a-1)*(a-1)^2*(1-g)^(b-1)/g^2 - g^(a-1)*(a-1)*(1-g)^(b-1)/g^2 
		- 2*g^(a-1)*(a-1)*(1-g)^(b-1)*(b-1)/(g*(1-g)) + g^(a-1)*(1-g)^(b-1)*(b-1)^2/(1-g)^2 
		- g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g)^2)/(g^(a-1)*(1-g)^(b-1))	
		- ((g^(a-1)*(a-1)*(1-g)^(b-1)/g-g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g))*(a-1)/(g^(a-1)*(1-g)^(b-1)*g))
		+ ((g^(a-1)*(a-1)*(1-g)^(b-1)/g-g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g))*(b-1)/(g^(a-1)*(1-g)^(b-1)*(1-g))))
	return(list(g=grad, h=hess))	
}

cormod <- function(fulldata, K, guess, smooth = TRUE) 
{  
	fulldata <- as.matrix(fulldata) 
	nitems <- ncol(fulldata)         
	cormat <- cor(fulldata)      
	if (any(guess > 0)){
		for (i in 1:nitems){
			for (j in 1:nitems){
				if (i < j & K[i] == 2 & K[j] == 2 & guess[i]!= 0 ){         
					g1 <- guess[i]
					g2 <- guess[j]
					tabp <- tab <- table(fulldata[ ,i],fulldata[ ,j])/length(fulldata[ ,i])
					w1 <- (1 - g1)
					w2 <- (1 - g2)
					tabp[1,1] <- tab[1,1]/(w1*w2)
					tabp[1,2] <- (w2*tab[1,2] - g2*tab[1,1])/(w1*w2)
					tabp[2,1] <- (w1*tab[2,1] - g1*tab[1,1])/(w1*w2)
					tabp[2,2] <- 1 - tabp[1,1] - tabp[1,2] - tabp[2,1]
					tabp <- round(tabp*length(fulldata[ ,i]))
					if(any(tabp < 0)) next	
					cormat[i,j] <- cormat[j,i] <- 
						abs(phi(tabp,6))^(1/1.15)*sign(phi(tabp,6))          		  
				} 
				if(i < j & K[i] == 2 & K[j] > 2 & guess[i]!= 0) 
					cormat[i,j] <- cormat[j,i] <- abs(cormat[i,j])^(1/1.15) * sign(cormat[i,j])
			}	
		}      
	} 
	cormat <- abs(cormat)^(1/1.15) * sign(cormat)  
	if(smooth){  
		eig <- eigen(cormat)
		negvalues <- eig$values < 0 
		if (any(negvalues)) {
			negeig <- sum(abs(eig$value[eig$value < 0])) 
			eig$value[eig$value < 0] <- 0
			L <- nitems/sum(eig$value)*eig$value[!negvalues]
			V <- eig$vector[ ,!negvalues] 
			cormat <- V %*% diag(L) %*% t(V)    
		}      
	}	
	cormat
}  


