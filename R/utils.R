# theta combinations
thetaComb <- function(theta, nfact)
{
	if (nfact == 1) Theta <- matrix(theta)
	else if (nfact == 2) Theta <- expand.grid(theta,theta)   
	else if (nfact == 3) Theta <- expand.grid(theta,theta,theta)  
	else if (nfact == 4) Theta <- expand.grid(theta,theta,theta,theta)
	else if (nfact == 5) Theta <- expand.grid(theta,theta,theta,theta,theta)        
	else if (nfact == 6) Theta <- expand.grid(theta,theta,theta,theta,theta,theta)        	
	if(nfact > 6) stop('Are you crazy?!?!? That\'s way too many factors for this quandrature method.')
	Theta <- as.matrix(Theta)	
	return(Theta)     
}

# start values
start.values <- function(fulldata, guess, Rpoly, nfact=2, bfactor = FALSE, nowarn = TRUE)
{	  	
	if (bfactor){ 
		suppressWarnings(FA <- psych::fa(Rpoly, 1, warnings = !nowarn))
		loads <- unclass(FA$load)
		cs <- sqrt(abs(FA$u))      
		dstart <- qnorm(colMeans(fulldata))/cs
		astart <- loads/cs
		startvalues <- cbind(astart,astart/2,dstart)
	} else {    
		suppressWarnings(FA <- psych::fa(Rpoly,nfact,rotate = 'none', warnings= !nowarn))	
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
Rotate <- function(F, rotate, Target = NULL, ...)
{	
    if(ncol(F) == 1) rotF <- list()
	if(rotate == 'promax'){
        rotF <- psych::Promax(F)
        rotF$orthogonal <- FALSE
	}    
    if(rotate == 'oblimin') rotF <- GPArotation::oblimin(F, ...)     
	if(rotate == 'quartimin') rotF <- GPArotation::quartimin(F, ...)
	if(rotate == 'targetT') rotF <- GPArotation::targetT(F, Target = Target, ...)
	if(rotate == 'targetQ') rotF <- GPArotation::targetQ(F, Target = Target, ...)
	if(rotate == 'pstT') rotF <- GPArotation::pstT(F, Target = Target, ...)
	if(rotate == 'pstQ') rotF <- GPArotation::pstQ(F, Target = Target, ...)
	if(rotate == 'oblimax') rotF <- GPArotation::oblimax(F, ...)
	if(rotate == 'entropy') rotF <- GPArotation::entropy(F, ...)
	if(rotate == 'quartimax') rotF <- GPArotation::quartimax(F, ...)
	if(rotate == 'varimax') rotF <- GPArotation::Varimax(F, ...)
	if(rotate == 'simplimax') rotF <- GPArotation::simplimax(F, ...)
	if(rotate == 'bentlerT') rotF <- GPArotation::bentlerT(F, ...)
	if(rotate == 'bentlerQ') rotF <- GPArotation::bentlerQ(F, ...)
	if(rotate == 'tandemI') rotF <- GPArotation::tandemI(F, ...)
	if(rotate == 'tandemII') rotF <- GPArotation::tandemII(F, ...)
	if(rotate == 'geominT') rotF <- GPArotation::geominT(F, ...)
	if(rotate == 'geominQ') rotF <- GPArotation::geominQ(F, ...)
	if(rotate == 'cfT') rotF <- GPArotation::cfT(F, ...)
	if(rotate == 'cfQ') rotF <- GPArotation::cfQ(F, ...)
	if(rotate == 'infomaxT') rotF <- GPArotation::infomaxT(F, ...)
	if(rotate == 'infomaxQ') rotF <- GPArotation::infomaxQ(F, ...)
	if(rotate == 'mccammon') rotF <- GPArotation::mccammon(F, ...)
	if(rotate == 'bifactorT') rotF <- GPArotation::bifactorT(F, ...)
	if(rotate == 'bifactorQ') rotF <- GPArotation::bifactorQ(F, ...)    		
	
	return(unclass(rotF))
}  

# MAP scoring for mirt
MAP.mirt <- function(Theta, a, d, guess, upper, patdata, itemloc, ML=FALSE)
{	
	itemtrace <- rep(0, ncol=length(patdata))
	Theta <- matrix(Theta, 1)
	for (i in 1:length(guess)){
		if(length(d[[i]]) == 1){
			itemtrace[itemloc[i]] <- P.mirt(a[i, ], d[[i]], Theta, guess[i], upper[i]) 
			itemtrace[itemloc[i] + 1] <- 1.0 - itemtrace[itemloc[i]]
		} else {
			itemtrace[itemloc[i]:(itemloc[i+1] - 1)] <- 
				P.poly(a[i, ], d[[i]], Theta, TRUE)	
		}
	}		
	L <- sum(log(itemtrace)[as.logical(patdata)])
	mu <- 0
	sigma <- 1
    L <- ifelse(ML, -L, (-1)*(L + sum(log(exp(-0.5*((Theta - mu)/sigma)^2)))))
	L  
}  

# MAP scoring for bfactor
MAP.bfactor <- function(Theta, a, d, guess, upper, patdata, logicalfact, itemloc, ML=FALSE)
{	
	itemtrace <- rep(0, ncol=length(patdata))
	Theta <- matrix(Theta, 1)
	for (i in 1:length(guess)){
		if(length(d[[i]]) == 1){
			itemtrace[itemloc[i]] <- P.mirt(a[i, logicalfact[i, ]], d[[i]], Theta, guess[i], upper[i]) 
			itemtrace[itemloc[i] + 1] <- 1.0 - itemtrace[itemloc[i]]
		} else {
			itemtrace[itemloc[i]:(itemloc[i+1] - 1)] <- 
				P.poly(a[i, logicalfact[i, ]], d[[i]], Theta, TRUE)	
		}
	}		
	L <- sum(log(itemtrace)[as.logical(patdata)])
	mu <- 0
	sigma <- 1
    L <- ifelse(ML, -L, (-1)*(L + sum(log(exp(-0.5*((Theta - mu)/sigma)^2)))))
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
		Pk[ ,i+1] <- P.mirt(lambda, zetas[i], Thetas, 0)		
	if(itemexp){
		P <- matrix(0,nrow(Thetas),ncat)		
		for(i in ncat:1)
			P[ ,i] <- Pk[ ,i] - Pk[ ,i+1]						
		Pk <- P
	}	
	return(Pk)
}

# Trace lines for mirt models
P.mirt <- function(a, d, Theta, g, u = 1)
{ 		
	traces <- .Call("traceLinePts", a, d, g, u, Theta)
	return(traces)
}

# Trace lines for partially compensetory models
P.comp <- function(a, d, thetas, c = 0, u = 1)
{
	nfact <- length(a)
	P <- rep(1,nrow(thetas))
	for(i in 1:nfact)
		P <- P * P.mirt(a[i], d[i], thetas[ ,i, drop=FALSE],0)
	P <- c + (u - c) * P
	P	
} 

# Estep for mirt
Estep.mirt <- function(pars, tabdata, Theta, prior, guess, upper, itemloc) 
{
	a <- pars$lambdas
	J <- nrow(a)
	nfact <- ncol(a)	
	nquad <- nrow(Theta)
	d <- pars$zetas    
	r <- tabdata[ ,ncol(tabdata)]
	X <- tabdata[ ,1:(ncol(tabdata) - 1)]	
	itemtrace <- matrix(0, ncol=ncol(X), nrow=nrow(Theta))
	r1 <- r0 <- matrix(0, ncol=length(guess), nrow=nrow(Theta))
	for (i in 1:J){
		if(length(d[[i]]) == 1){
			itemtrace[ ,itemloc[i] + 1] <- P.mirt(a[i, ], d[[i]], Theta, guess[i], upper[i]) 
			itemtrace[ ,itemloc[i]] <- 1.0 - itemtrace[ ,itemloc[i] + 1]
		} else {
			itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- 
				P.poly(a[i, ], d[[i]], Theta, TRUE)	
		}
	}		
    retlist <- .Call("Estep", itemtrace, prior, X, nfact, r)	    		
	return(retlist)
} 

# Trace lines for bfactor
P.bfactor <- function(a, d, Theta, g, u, patload)
{ 
	a <- a[patload]	
	if(length(d) > 1){
		ncat <- length(d) + 1
		nfact <- length(a)
		Pk <- matrix(0,nrow(Theta),ncat+1)
		Pk[,1] <- 1	
		for(i in 1:(ncat-1))			
			Pk[ ,i+1] <- P.mirt(a, d[i], Theta, 0)				
		P <- matrix(0,nrow(Theta),ncat)		
		for(i in ncat:1)
			P[ ,i] <- Pk[ ,i] - Pk[ ,i+1]		
	} else P <- .Call("traceLinePts", a, d, g, u, Theta)		
	return(P)
}

# Estep for bfactor
Estep.bfactor <- function(pars, tabdata, Theta, prior, guess, upper, specific, sitems, itemloc) 
{	
	a <- pars$lambdas
	logicalfact <- attr(pars, 'lamsel')
	nfact <- ncol(a)
	J <- nrow(a)
	nquad <- nrow(Theta)	
	d <- pars$zetas    
	r <- tabdata[ ,ncol(tabdata)]
	X <- tabdata[ ,1:(ncol(tabdata) - 1)]	
	itemtrace <- matrix(0, ncol=ncol(X), nrow=nrow(Theta))
	r1 <- r0 <- matrix(0, ncol=length(guess), nrow=nrow(Theta))
	for (i in 1:J){
		atmp <- a[i, logicalfact[i, ]]
		if(length(d[[i]]) == 1){
			itemtrace[ ,itemloc[i] + 1] <- P.mirt(atmp, d[[i]], Theta, guess[i], upper[i]) 
			itemtrace[ ,itemloc[i]] <- 1.0 - itemtrace[ ,itemloc[i] + 1]
		} else {
			itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- 
				P.poly(atmp, d[[i]], Theta, TRUE) 
		}
	}			
	retlist <- .Call("Estepbfactor", itemtrace, prior, X, r, sitems)	
	r1 <- matrix(0, nrow(Theta), ncol(X))	
	for (i in 1:J)
	    r1[ ,itemloc[i]:(itemloc[i+1]-1)] <- 		
			retlist$r1[ ,itemloc[i]:(itemloc[i+1]-1) + (specific[i] - 1)*ncol(X) ]	    
		
	return(list(r1=r1, expected=retlist$expected))	
}      

# MH sampler for theta values
draw.thetas <- function(theta0,lambdas,zetas,guess,upper = rep(1,length(K)),fulldata,K,itemloc,cand.t.var,
	prior.t.var = diag(ncol(theta0)), prior.mu = rep(0,ncol(theta0)), estComp = rep(FALSE,length(K)),
	prodlist = NULL, itemtype = NULL) 
{ 			
	N <- nrow(fulldata)
	J <- length(K)
	nfact <- ncol(theta0)					
	unif <- runif(N)
	if(nfact > 1)		
		theta1 <- theta0 + mvtnorm::rmvnorm(N,prior.mu, 
			diag(rep(cand.t.var,ncol(theta0)))) 
	else
		theta1 <- theta0 + rnorm(N,prior.mu,sqrt(cand.t.var))							
	den0 <- mvtnorm::dmvnorm(theta0,prior.mu,prior.t.var)
	den1 <- mvtnorm::dmvnorm(theta1,prior.mu,prior.t.var)		
	if(!is.null(prodlist)){
		theta0 <- prodterms(theta0,prodlist)
		theta1 <- prodterms(theta1,prodlist)	
	}	
	ThetaDraws <- .Call("drawThetas", unif, den0, den1, lambdas, zetas, guess, upper,
					theta0, theta1,	fulldata, (itemloc-1), as.numeric(estComp))
	log.lik <- ThetaDraws$cdloglik
	accept <- as.logical(ThetaDraws$accept)				
	theta1[!accept,] <- theta0[!accept,]	
	if(!is.null(prodlist)) 
		theta1 <- theta1[ ,1:(ncol(lambdas) - length(prodlist)), drop=FALSE]
	attr(theta1, "Proportion Accepted") <- sum(accept)/N 				
	attr(theta1, "log.lik") <- log.lik	
	return(theta1) 
}	

# Analytical derivatives for covariances parameters
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
				as.numeric(invSig),
				as.numeric(cMeans),				
				as.numeric(Zdif),
				as.integer(N),
				as.integer(nfact),
				as.integer(npars))				
	sel <- sel[c(rep(TRUE,nfact),as.logical(selcov))]	
	h <- h[sel,sel] 
	list(h=h,g=g) 
}

# Analytical derivatives for dichotomous items
dpars.dich <- function(lambda, zeta, g, u, dat, Thetas, estGuess)
{
	nfact <- length(lambda)
	P <- P.mirt(lambda, zeta, Thetas, g, u)						
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
			hess[hsize,1+i] <- hess[1+i,hsize] <- sum(-thetas[,i]*Pstar*Qstar*((r/P - (f-r)/Q) + 
                Qstar*(1-c)*const2))		
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
		L2 <- colSums((dat-P) * Thetas)
		dL <- c(L1,L2)		
		d2L <- matrix(0,nfact+1, nfact+1)						
		L11 <- .Call("dichOuter", Thetas, PQ, nrow(Thetas))
		if(nfact > 1) d2L[1:nfact+1, 1:nfact+1] <- -L11
			 else d2L[nfact+1, nfact+1] <- -L11 				
		d2L[1, 1] <- (-1)*sum(PQ)		
		d2L[1, 1:nfact+1] <- d2L[1:nfact+1, 1] <- (-1)*colSums(PQ * Thetas)
	}	
	list(grad = dL, hess = d2L)
}

# Analytical derivatives for partially compensatory items
dpars.comp <- function(lambda,zeta,g,dat,Thetas,estg = FALSE)
{	    
	nfact <- length(lambda)	
	pars <- c(zeta,lambda,g)
	if(estg){
		grad <- function(pars, r, thetas){
			f <- 1			
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
				Pk <- P.mirt(a[i],d[i],thetas[ , i, drop=FALSE],0)
				Qk <- 1 - Pk
				dd[i] <- sum((1-c)*Pstar*Qk*const1)
				da[i] <- sum((1-c)*Pstar*Qk*thetas[,i]*const1)
			}
			return(c(dd,da,dc))
		}		
		hess <- function(pars, r, thetas){
			f <- 1			
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
						Pk <- P.mirt(a[k],d[k],thetas[ , k, drop=FALSE],0)
						Qk <- 1 - Pk	
						Pm <- P.mirt(a[m],d[m],thetas[ , m, drop=FALSE],0)
						Qm <- 1 - Pm									
						if(i == j && d1[1] == 'd'){
							hess[i,i] <- sum((1-c)*Pstar*Qk*(const1*((1-c)*Qk - Pk) - Pstar*Qk*(1-c)*const2))
							next
						}
						if(i == j && d1[1] == 'a'){
							hess[i,i] <- sum((1-c)*thetas[,k]^2*Pstar*Qk*(const1*((1-c)*Qk - Pk) - Pstar*Qk*
                                (1-c)*const2))
							next		
						}
						if(i == j && d1[1] == 'c'){
							hess[i,i] <- -sum(Qstar^2 * const2)
							next		
						}	
						if(d1[1] == 'a' && d2[1] == 'a'){
							hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*thetas[,m]*Qk*Pstar*Qm*(const1 - 
                                Pstar*(1-c)*const2))
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
							hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*Pstar*Qk*(const1*((1-c)*Qk - Pk) - 
                                Pstar*Qk*(1-c)*const2))
							next	
						}
						if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
							hess[i,j] <- hess[j,i] <- sum((1-c)*Qk*thetas[,m]*Pstar*Qm*(const1 - 
                                Pstar*(1-c)*const2))
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
			Pk <- P.mirt(lambda[i],zeta[i],Thetas[ , i, drop=FALSE],0)
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
					Pk <- P.mirt(lambda[k],zeta[k],Thetas[ , k, drop=FALSE],0)
					Qk <- 1 - Pk	
					Pm <- P.mirt(lambda[m],zeta[m],Thetas[ , m, drop=FALSE],0)
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

# Analytical derivatives for polychotomous items (ordinal)
dpars.poly <- function(lambda,zeta,dat,Thetas)
{  
	nzeta <- length(zeta)				
	P <- P.poly(lambda,zeta,Thetas)			    	
	ret <- .Call("dparsPoly", P, Thetas, dat, nzeta)	
	return(ret)	
}

# Gamma correlation, mainly for obtaining a sign
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

# Beta prior for grad and hess
betaprior <- function(a,b,g,W=20)
{
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

# Approximation to polychoric matrix for initial values
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
						abs(psych::phi(tabp,6))^(1/1.15)*sign(psych::phi(tabp,6))          		  
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

# Product terms in confmirt
prodterms <- function(theta0, prodlist)
{
	products <- matrix(1, ncol = length(prodlist), nrow = nrow(theta0))
	for(i in 1:length(prodlist)){
		tmp <- prodlist[[i]]
		for(j in 1:length(tmp)) 
			products[ ,i] <- products[ ,i] * theta0[ ,tmp[j]]	
	}	
	ret <- cbind(theta0,products)
	ret
}

# Extract model matricies and values for user specified confmirt.model()
model.elements <- function(model, factorNames, nfactNames, nfact, J, K, fulldata, itemloc, data, N,
  estGuess, guess, upper, estUpper, guess.prior.n, itemnames, exploratory)
{    
  hasProdTerms <- ifelse(nfact == nfactNames, FALSE, TRUE)
  prodlist <- NULL
  if(hasProdTerms){
    tmp <- factorNames[grepl('\\(',factorNames)]
    tmp2 <- factorNames[!grepl('\\(',factorNames)] 
    tmp <- gsub("\\(","",tmp)	
    tmp <- gsub("\\)","",tmp)
    tmp <- gsub(" ","",tmp)
    prodlist <- strsplit(tmp,"\\*")
    for(j in 1:length(prodlist)){
      for(i in 1:nfact)
        prodlist[[j]][prodlist[[j]] == tmp2[[i]]] <- i		
      prodlist[[j]] <- as.numeric(prodlist[[j]])	
    }		
  } 
  
  #slopes specification
  estlam <- matrix(FALSE, ncol = nfactNames, nrow = J)	
  for(i in 1:nfactNames){
    tmp <- model[model[ ,1] == factorNames[i],2]
    if(any(regexpr(",",tmp)))
      tmp <- strsplit(tmp,",")[[1]]
    popout <- c()	
    for(j in 1:length(tmp)){
      if(regexpr("-",tmp[j]) > 1){
        popout <- c(popout,j)
        tmp2 <- as.numeric(strsplit(tmp[j],"-")[[1]])
        tmp2 <- as.character(tmp2[1]:tmp2[2])
        tmp <- c(tmp,tmp2)
      }
    }
    if(length(popout != 0))	
      estlam[as.numeric(tmp[-popout]),i] <- TRUE
    else 
      estlam[as.numeric(tmp),i] <- TRUE
  }
  lambdas <- ifelse(estlam, .5, 0)	
  
  #PARTCOMP
  estComp <- rep(FALSE,J)
  if(any(model[,1] == 'PARTCOMP')){
    tmp <- model[model[,1] == 'PARTCOMP',2]		
    tmp <- strsplit(tmp,",")[[1]]
    tmp <- gsub(" ","",tmp)		
    for(j in 1:length(tmp)){
      if(regexpr("-",tmp[j]) > 1){				
        tmp2 <- as.numeric(strsplit(tmp[j],"-")[[1]])				
        estComp[tmp2[1]:tmp2[2]] <- TRUE
      }
    }
    if(any(is.numeric(suppressWarnings(as.numeric(tmp)))))
      for(i in 1:length(tmp))
        estComp[suppressWarnings(as.numeric(tmp))] <- TRUE				
  }
  if(nfact == 1) estComp <- rep(FALSE,J)	
  
  #INT
  cs <- sqrt(abs(1-rowSums(lambdas^2)))	
  zetas <- rep(NA,200)	
  loc <- 1	
  for(i in 1:J){
    if(estComp[i]){ 
      div <- ifelse(cs[i] > .25, cs[i], .25)
      tmp <- rep(qnorm(mean(fulldata[,itemloc[i]]))/div, sum(estlam[i,]))
      zetas[loc:(loc+length(tmp)-1)] <- tmp
      loc <- loc + length(tmp)
      next
    }
    if(K[i] == 2){
      div <- ifelse(cs[i] > .25, cs[i], .25)		
      zetas[loc] <- qnorm(mean(fulldata[,itemloc[i]]))/div
      loc <- loc + 1
    } else {			
      temp <- table(data[,i])[1:(K[i]-1)]/N
      temp <- cumsum(temp)
      div <- ifelse(cs[i] > .25, cs[i], .25)		
      zetas[loc:(loc+K[i]-2)] <- qnorm(1 - temp)/div	
      loc <- loc + K[i] - 1	
    }		
  }
  zetas <- zetas[!is.na(zetas)]
  estzetas <- list()
  estzetas2 <- c()
  ind1 <- 1
  for(i in 1:J){
    if(estComp[i]){
      estzetas[[i]] <- rep(TRUE,sum(estlam[i,]))		
      estzetas2 <- c(estzetas2,estzetas[[i]])
      ind1 <- ind1 + sum(estlam[i,]) - 1
    } else {
      estzetas[[i]] <- rep(TRUE,length((ind1):(K[i] + ind1 - 2)))		
      estzetas2 <- c(estzetas2,estzetas[[i]])
      ind1 <- ind1 + K[i] - 1
    }	
  }
    
  #MEANS
  find <- 1:nfact
  gmeans <- rep(0,nfact)
  estgmeans <- rep(FALSE,nfact)	
  if(any(model[,1] == 'MEAN')){
    tmp <- model[model[,1] == 'MEAN',2]		
    tmp <- strsplit(tmp,",")[[1]]
    tmp <- gsub(" ","",tmp)
    for(i in 1:length(tmp)){
      tmp2 <- strsplit(tmp[i],"eq",fixed=TRUE)[[1]]
      ind1 <- find[tmp2[1] == factorNames]			
      gmeans[ind1] <- as.numeric(tmp2[2])
    }
  }
  
  #COV
  estgcov <- constgcov <- matrix(FALSE,nfact,nfact)
  equalcov <- list()
  equalcovind <- 1
  if(any(model[,1] == 'COV')){
    tmp <- model[model[,1] == 'COV',2]		
    tmp <- strsplit(tmp,",")[[1]]
    tmp <- gsub(" ","",tmp)
    for(i in 1:length(tmp)){
      if(regexpr("eq",tmp[i]) > 1){
        tmp2 <- strsplit(tmp[i],"eq",fixed=TRUE)[[1]]
        suppressWarnings(value <- as.numeric(tmp2[length(tmp2)]))
        if(!is.na(value)){
          tmp2 <- strsplit(tmp2[1],"*",fixed=TRUE)[[1]]
          ind1 <- find[tmp2[1] == factorNames]
          ind2 <- find[tmp2[2] == factorNames]
          constgcov[ind1,ind2] <- constgcov[ind2,ind1] <- value
        } else {
          tmp2 <- strsplit(tmp2,"*",fixed=TRUE)
          equalcov[[equalcovind]] <- matrix(FALSE,nfact,nfact)
          for(j in 1:length(tmp2)){
            ind1 <- find[tmp2[[j]][1] == factorNames]
            ind2 <- find[tmp2[[j]][2] == factorNames]
            estgcov[ind1,ind2] <- estgcov[ind2,ind1] <- TRUE						
            equalcov[[equalcovind]][ind1,ind2] <- equalcov[[equalcovind]][ind2,ind1] <- TRUE
          }
          equalcovind <- equalcovind + 1
        }
      } else {
        tmp2 <- strsplit(tmp[i],"*",fixed=TRUE)[[1]]				
        ind1 <- find[tmp2[1] == factorNames]
        ind2 <- find[tmp2[2] == factorNames]
        estgcov[ind1,ind2] <- estgcov[ind2,ind1] <- TRUE
      }	
    }
  }
  gcov <- ifelse(estgcov,.1,0) + constgcov
  diag(gcov) <- 1	
  tmp <- matrix(FALSE,nfact,nfact)
  tmp[lower.tri(tmp,diag=TRUE)] <- estgcov[lower.tri(tmp,diag=TRUE)]
  selgcov <- lower.tri(tmp,diag = TRUE)
  estgcov <- tmp
  
  #Housework
  loc1 <- 1
  lamind <- zetaind <- guessind <- upperind <- sind <- c()	
  for(i in 1:J){
    if(estComp[i])
      zetaind <- c(zetaind, loc1:(loc1+(length(estzetas[[i]])-1)))		
    else zetaind <- c(zetaind,loc1:(loc1 + K[i] - 2))
    lamind <- c(lamind,max(zetaind + 1):(max(zetaind) + nfactNames))		
    guessind <- c(guessind,max(lamind + 1):max(lamind + 1 ))
    upperind <- c(upperind,max(lamind + 2):max(lamind + 2 ))
    sind <- c(sind, estzetas[[i]], estlam[i,], estGuess[i], estUpper[i])
    loc1 <- loc1 + nfactNames + sum(estzetas[[i]]) + 2	
  }	
  sind <- c(sind, estgmeans, estgcov[lower.tri(estgcov,diag=TRUE)])
  npars <- length(sind)
  pars <- rep(0,npars)
  groupind <- (npars - length(c(gmeans,gcov[lower.tri(gcov,diag=TRUE)]))+1):npars
  meanind <- groupind[1:nfact]
  covind <- groupind[-(1:nfact)]	
  pars[lamind] <- t(lambdas)
  pars[zetaind] <- zetas
  pars[guessind] <- guess
  pars[upperind] <- upper
  pars[groupind] <- c(gmeans,gcov[lower.tri(gcov,diag=TRUE)])
  parnames <- rep('',npars)
  parnames[lamind] <- paste('lam',1:length(lamind),sep='')  
  parnames[zetaind] <- paste('zeta',1:length(zetaind),sep='')
  parnames[guessind] <- paste('guess',1:length(guessind),sep='')  
  parnames[upperind] <- paste('upper',1:length(upperind),sep='')
  parnames[groupind] <- c(paste('gmeans',1:length(gmeans),sep=''), 
	paste('gcov',1:length(gcov[lower.tri(gcov,diag=TRUE)]),sep=''))
  names(pars) <- names(sind) <- parnames
  
  parcount <- list(lam=estlam, zeta=estzetas2, guess=estGuess, upper=estUpper, cov=estgcov, mean=estgmeans)
  parind <- 1:npars
  loc1 <- 1
  parcount$lam <- matrix(lamind,J,byrow=TRUE)	
  zetaind2 <- estzetas
  k <- 1
  for(i in 1:J){
    for(j in 1:length(zetaind2[[i]])){
      zetaind2[[i]][j] <- zetaind[k]			
      k <- k + 1
    }
  }
  
  zetas <- list()  
  ind1 <- 1
  for(i in 1:J){ 
    zetas[[i]] <- pars[zetaind][ind1:(ind1+sum(estzetas[[i]])-1)]	
    ind1 <- ind1 + sum(estzetas[[i]])
  }
  names(zetas) <- names(zetaind2) <- itemnames   
  parcount$zeta <- zetaind2
  parcount$guess <- guessind
  parcount$upper <- upperind
  parcount$mean <- meanind
  parcount$cov <- matrix(0,nfact,nfact)
  parcount$cov[selgcov] <- covind		
  constvalues <- matrix(0,ncol = 2, npars)	
  
  #ADDITIONAL SPECS
  constvalues[parcount$cov[constgcov[selgcov] != 0], ] <- c(1,constgcov[constgcov[selgcov] != 0])
  equalconstr <- list()
  equalind <- 1
  if(length(equalcov) > 0){
    for(i in 1:length(equalcov)){
      equalconstr[[equalind]] <- parcount$cov[equalcov[[i]][lower.tri(estgcov,diag=TRUE)]]
      equalind <- equalind + 1
    }	
  }
  if(any(model[ ,1] == 'SLOPE')){
    tmp <- model[model[ ,1] == "SLOPE",2]
    if(any(regexpr(",",tmp)))
      tmp <- strsplit(tmp,",")[[1]]
    tmp <- gsub('\\s+','', tmp, perl = TRUE)	
    for(i in 1:length(tmp)){
      tmp2 <- strsplit(tmp[i],'eq')[[1]]
      suppressWarnings(attempt <- as.numeric(tmp2))
      if(any(!is.na(attempt))){
        value <- attempt[!is.na(attempt)]
        tmp3 <- tmp2[is.na(attempt)]
        tmp3 <- strsplit(tmp3,"@")								
        for(j in 1:length(tmp3)){					
          loc1 <- tmp3[[j]][1] == factorNames
          loc2 <- as.numeric(tmp3[[j]][2])
          constvalues[parcount$lam[loc2,loc1], ] <- c(1,value)
        }					
      } else {
        tmp3 <- strsplit(tmp2,"@")				
        equalconstr[[equalind]] <- rep(0,length(tmp3)) 
        for(j in 1:length(tmp3)){					
          loc1 <- tmp3[[j]][1] == factorNames
          loc2 <- as.numeric(tmp3[[j]][2])
          equalconstr[[equalind]][j] <- parcount$lam[loc2,loc1] 					
        }
        if(any(equalconstr[[equalind]] == 0)) stop("Improper constrainst specification.")
        equalind <- equalind + 1					
      }
    }
  }	
  zetaind2 <- estzetas
  k <- 1
  for(i in 1:J){
    for(j in 1:length(zetaind2[[i]])){
      zetaind2[[i]][j] <- zetaind[k]
      k <- k + 1
    }
  }
  if(any(model[,1] == 'INT')){
    tmp <- model[model[,1] == 'INT',2]
    if(any(regexpr(",",tmp)))
      tmp <- strsplit(tmp,",")[[1]]
    tmp <- gsub('\\s+','', tmp, perl = TRUE)	
    for(i in 1:length(tmp)){
      tmp2 <- strsplit(tmp[i],'eq')[[1]]
      suppressWarnings(attempt <- as.numeric(tmp2))
      if(any(!is.na(attempt))){
        value <- attempt[!is.na(attempt)]
        tmp3 <- tmp2[is.na(attempt)]
        tmp3 <- strsplit(tmp3,"@")								
        for(j in 1:length(tmp3)){					
          loc1 <- as.numeric(tmp3[[j]][1])
          loc2 <- as.numeric(tmp3[[j]][2])
          constvalues[zetaind2[[loc1]][loc2], ] <- c(1,value)
        }					
      } else {
        tmp3 <- strsplit(tmp2,"@")				
        equalconstr[[equalind]] <- rep(0,length(tmp3)) 
        for(j in 1:length(tmp3)){	
          loc1 <- as.numeric(tmp3[[j]][1])
          loc2 <- as.numeric(tmp3[[j]][2])
          equalconstr[[equalind]][j] <- zetaind2[[loc1]][loc2]					
        }
        if(any(equalconstr[[equalind]] == 0)) stop("Improper constraint specification.")
        equalind <- equalind + 1					
      }
    }
  }
  if(!all(sort(abs(colSums(estlam[ ,1:nfact,drop=FALSE]) - J)) >= 0:(nfact-1)) && 
    length(equalconstr) == 0 && !exploratory) stop('Slope parameters are not uniquely identified.')
  if(!all(sort(abs(colSums(estlam[ ,1:nfact,drop=FALSE]) - J)) >= 0:(nfact-1)) && !exploratory)
    warning('Slope parameters may not be uniquely identified.')	
  
  #PRIOR, 1 == norm, 2== beta
  parpriors <- list()
  parpriorscount <- 1
  if(sum(estGuess) > 0){
    for(i in 1:J){
      if(estGuess[i]){
        a <- guess[i] * guess.prior.n[i]
        b <- (1 - guess[i]) * guess.prior.n
        parpriors[[parpriorscount]] <- c(2,guessind[i],a,b)						
        parpriorscount <- parpriorscount + 1			
      }
    }
  }		
  if(any(model[,1] == 'PRIOR')){
    tmp <- model[model[,1] == 'PRIOR',2]
    if(any(regexpr(",",tmp)))
      tmp <- strsplit(tmp,",")[[1]]
    tmp <- gsub('\\s+','', tmp, perl = TRUE)	
    for(i in seq(1,length(tmp),by=2)){
      tmp2 <- strsplit(tmp[i],"\\(")[[1]]
      tmp3 <- as.numeric(strsplit(tmp[i+1],"\\)@")[[1]])			
      if(tmp2[1] == 'N')				
        parpriors[[parpriorscount]] <- c(1,tmp3[2],as.numeric(tmp2[2]),tmp3[1])
      if(tmp2[1] == 'B')
        parpriors[[parpriorscount]] <- c(2,tmp3[2],as.numeric(tmp2[2]),tmp3[1])
      parpriorscount <- parpriorscount + 1	
    }
  }  
  
  nconstvalues <- sum(constvalues[,1] == 1)
  if(length(equalconstr) > 0)    
      for(i in 1:length(equalconstr))
          nconstvalues <- nconstvalues + length(equalconstr[[i]]) - 1
  itemtype <- rep('2PL', J)
  itemtype[estGuess & estUpper] <- '4PL'
  itemtype[estGuess &! estUpper] <- '3PL'
  itemtype[estUpper &! estGuess] <- '3PLu'
  itemtype[estComp &! estGuess] <- 'N2PL'
  itemtype[estComp & estGuess] <- 'N3PL'
  itemtype[K > 2] <- 'ordinal' 
      
  val <- list(pars=pars, lambdas=lambdas, zetas=zetas, gmeans=gmeans, gcov=gcov, 
    guess=guess, upper=upper, constvalues=constvalues)
  est <- list(estlam=estlam, estComp=estComp, estzetas=estzetas, estzetas2=estzetas2, 
    estgcov=estgcov, estgmeans=estgmeans, estGuess=estGuess, estUpper=estUpper)
  ind <- list(equalind=equalind, equalconstr=equalconstr, parpriorscount=parpriorscount, 
    prodlist=prodlist, parpriors=parpriors, sind=sind, lamind=lamind, zetaind=zetaind, 
    zetaindlist=zetaind2, guessind=guessind, upperind=upperind, groupind=groupind, 
    meanind=meanind, covind=covind, parind=parind)
  ret <- list(val=val, est=est, ind=ind, parcount=parcount, npars=npars, nconstvalues=nconstvalues,
              itemtype=itemtype)
  ret
}

# Take long parameter form and return list of pars for polymirt (obsolete)
sortPars <- function(pars, indlist, nfact, estGuess)
{
	lambdas <- matrix(pars[indlist$lamind],ncol=nfact,byrow=TRUE)	
	J <- nrow(lambdas)		
	zetas <- list()
	for(i in 1:J)
		zetas[[i]] <- pars[indlist$zetaind[[i]]]
	guess <- upper <- rep(0,J)
	guess[estGuess] <- pars[indlist$gind]		
	
	return(list(lambdas=lambdas, zetas=zetas, guess=guess))
}

# Take long parameter form and return list of pars for confmirt
sortParsConfmirt <- function(pars, indlist, nfact, nfactNames)
{
	J <- length(indlist$guessind)
	lambdas <- matrix(pars[indlist$lamind],J,nfactNames,byrow=TRUE)
	zetas <- list()
	for(i in 1:J)
		zetas[[i]] <- pars[indlist$zetaindlist[[i]]]
	guess <- pars[indlist$guessind]
    upper <- pars[indlist$upperind]
	mu <- pars[indlist$meanind]
	sig <- matrix(0, nfact, nfact)
	sig[lower.tri(sig, diag=TRUE)] <- pars[indlist$covind]
	if(nfact > 1)
		sig <- sig + t(sig) - diag(diag(sig))							
	
	return(list(lambdas=lambdas, zetas=zetas, guess=guess, upper=upper, mu=mu, sig=sig))
}

# Ramsey rate acceleration adjustment for EM
rateChange <- function(pars, lastpars1, lastpars2)
{
	p <- unlist(pars)	
	lp1 <- unlist(lastpars1)
	lp2 <- unlist(lastpars2)
	rate <- rep(0, length(p))
	d1 <- lp1 - p
	d2 <- lp2 - p
	rate <- ifelse(abs(d1) > 0.001 & (d1*d2 > 0.0) & (d1/d2 < 1.0),
		(1 - (1 - rate) * (d1/d2)),
		0)	    		
	rate[p > 4] <- 0
	rate[p < -4] <- 0    
	p <- lp1*rate*(-2) + (1 - rate*(-2))*p
	parsret <- rebuildPars(p, pars)	
	parsret
}

# Rebuild parameters given a list
rebuildPars <- function(p, pars)
{
	names(p) <- NULL
	pars2 <- pars
	pars2$lambdas <- matrix(p[1:length(pars$lambdas)], ncol=ncol(pars$lambdas), 
		nrow=nrow(pars$lambdas)) 
	ind1 <- length(pars$lambdas) + 1
	for(i in 1:length(pars$zetas) ){
		ind2 <- ind1 + length(pars$zetas[[i]]) - 1
		pars2$zetas[[i]] <- p[ind1:ind2]
		ind1 <- ind1 + length(pars$zetas[[i]])
	}			
	return(pars2)
}

# Rotate lambda coefficients
rotateLambdas <- function(so){
    F <- 
    F <- so$rotF %*% t(chol(so$fcor))
    h2 <- so$h2
    h <- matrix(rep(sqrt(1 - h2), ncol(F)), ncol = ncol(F))
    a <- F / h
    a    
}

d2r <-function(d) pi*d/180
    
