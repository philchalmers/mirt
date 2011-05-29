# theta combinations
thetaComb <- function(theta, nfact)
{
  if (nfact == 1) Theta <- matrix(theta)
    else if (nfact == 2) Theta <- expand.grid(theta,theta)   
    else if (nfact == 3) Theta <- expand.grid(theta,theta,theta)  
    else if (nfact == 4) Theta <- expand.grid(theta,theta,theta,theta)
    else if (nfact == 5) Theta <- expand.grid(theta,theta,theta,theta,theta)        
  return(Theta)     
}

# start values
start.values <- function(fulldata, guess, Rpoly, nfact=2, bfactor = FALSE, nowarn = TRUE)
{
  if (length(guess) == 1) guess <- rep(guess,ncol(fulldata))
    else if (length(guess) > ncol(fulldata) || length(guess) < ncol(fulldata)) 
      stop("The number of guessing parameters is incorrect.")  
  if(nowarn) options(warn = -1)	  
  if (bfactor)
  { 
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
  if(nowarn) options(warn = 0)
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
P.poly <- function(lambda, zetas, Thetas, itemexp = FALSE){	
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
P.mirt <- function(a, d, Theta, g){ 
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
  Estep.mirt <- function(pars, tabdata, Theta, prior, guess) {
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
  
 P.bfactor <- function(a, d, Theta, g, patload){ 
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
		prior.t.var = diag(ncol(theta0))) { 		
		
		N <- nrow(fulldata)
		J <- length(K)
		nfact <- 1:ncol(theta0)		
		locz <- 1		
		P0 <- P1 <- matrix(0,N,J)		
		unif <- runif(N)
        if(length(nfact) > 1)		
		  theta1 <- theta0 + rmvnorm(N,rep(0,ncol(theta0)), 
			diag(rep(sqrt(cand.t.var),ncol(theta0)))) 
        else
          theta1 <- theta0 + rnorm(N,0,sqrt(cand.t.var))							
		den0 <- dmvnorm(theta0,rep(0,length(nfact)),prior.t.var)
		den1 <- dmvnorm(theta1,rep(0,length(nfact)),prior.t.var)						
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
						as.integer(ncol(lambdas)))
		log.lik <- accept[N+1]			
		accept <- as.logical(accept[-(N+1)])				
		theta1[!accept,] <- theta0[!accept,]			
		attr(theta1, "Proportion Accepted") <- sum(accept)/N 				
		attr(theta1, "log.lik") <- log.lik		
		return(theta1) 
	}	
	
	d.group <- function(grouplist,theta){		
		tr <- function(x) sum(diag(x))
		x <- theta
		u <- grouplist$u	
		sig <- grouplist$sig
		N <- nrow(x)
		nfact <- length(u)
		selcov <- matrix(FALSE,nfact,nfact)
		npars <- length(sig) + nfact
		for(i in 1:nfact)
			for(j in 1:nfact)
				if(i <= j) selcov[j,i] <- TRUE
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
  
  

 
 