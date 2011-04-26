# theta combinations
thetaComb <- function(theta, nfact)
{
  if (nfact == 1) Theta <- theta
    else if (nfact == 2) Theta <- expand.grid(theta,theta)   
    else if (nfact == 3) Theta <- expand.grid(theta,theta,theta)  
    else if (nfact == 4) Theta <- expand.grid(theta,theta,theta,theta)
    else if (nfact == 5) Theta <- expand.grid(theta,theta,theta,theta,theta)
    else if (nfact == 6) Theta <- expand.grid(theta,theta,theta,theta,theta,theta)
    else if (nfact == 7) Theta <- expand.grid(theta,theta,theta,theta,theta,theta,theta)
    else if (nfact == 8) Theta <- expand.grid(theta,theta,theta,theta,theta,theta,theta,theta)
    else if (nfact == 9) Theta <- expand.grid(theta,theta,theta,theta,theta,theta,theta,theta,theta)
    else if (nfact == 10) Theta <- expand.grid(theta,theta,theta,theta,theta,theta,theta,theta,theta,theta)    
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
    u[u < 0 ] <- .25
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
    if(patdata[j] == 1) L <- log(Pmirt(a[j, ],d[j],Theta,guess[j])) + L
	  else L <- log(1 - Pmirt(a[j, ],d[j],Theta,guess[j])) + L	
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
    if(patdata[j] == 1) L <- log(Pbfactor(a[j, ],d[j],Theta,guess[j],logicalfact[j, ])) + L
	  else L <- log(1 - Pbfactor(a[j, ],d[j],Theta,guess[j],logicalfact[j, ])) + L	
  }
  mu <- 0
  sigma <- 1
  L <- (-1)*(L + sum(log(exp(-0.5*((Theta - mu)/sigma)^2))))
  L  
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
 
mirt.MHRM <- function(fulldata, nfact, pars, guess, Rpoly,
	mk = 3, SEM.cycles = 20, max.cycles = 1000, tol = .001){

	all.P.mirt <- function(Thetas, pars, guess, N, J) {
		nfact <- ncol(pars) - 1
		P <- matrix(0,N, J)
		for(i in 1:J)	
		   P[ ,i] <- P.mirt(pars[i,1:nfact], pars[i,nfact + 1], Thetas, guess[i])
		return(P)	
	}
	draw.thetas <- function(theta0,pars,guess,fulldata,cand.t.var) {     
		N <- nrow(fulldata)
		J <- ncol(fulldata)
		nfact <- ncol(pars) - 1
		prior.t.var <- diag(nfact)	
		if(nfact > 1) 
			theta1 <- theta0 + rmvnorm(N,rep(0,nfact), diag(rep(sqrt(cand.t.var),nfact))) 
		else 
			theta1 <- theta0 + rnorm(N,0,sqrt(cand.t.var))		
		tmp <- all.P.mirt(theta0,pars,guess,N,J) 
		tmp <- ifelse(fulldata,tmp,1-tmp)	
		irt0 <- rowSums(log(tmp)) + dmvnorm(theta0,rep(0,nfact),prior.t.var,log=TRUE)
		tmp <- all.P.mirt(theta1,pars,guess,N,J)
		tmp <- ifelse(fulldata,tmp,1-tmp)	
		irt1 <- rowSums(log(tmp)) + dmvnorm(theta1,rep(0,nfact),prior.t.var,log=TRUE)		
		accept <- irt1 - irt0
		accept <- ifelse(accept>0,0,accept)
		accept <- ifelse(runif(N) < exp(accept),TRUE,FALSE) 		
		theta1[!accept,] <- theta0[!accept,]	
		attr(theta1, "Proportion Accepted") <- sum(accept)/N 		
		return(theta1) 
	}
	dpars <- function(pars,guess,item,Thetas) {
		nfact <- length(pars) - 1
		P <- P.mirt(pars[1:nfact], pars[nfact + 1], Thetas, guess)								
		PQ <- P*(1-P)
		L1 <- L11 <- L12 <- c()
		for(i in 1:nfact){
			L1[i] <- sum((item-P)*Thetas[,i])
			L11[i] <- (-1)*sum(PQ*Thetas[,i]^2)		
		}	
		L2 <- sum(item-P)	
		L22<- (-1)*sum(PQ)
		dL <- c(L1,L2)
		d2L<- diag(c(L11,L22))
		for(i in 1:nfact)
			for(j in 1:nfact)
				if(j > i) d2L[i,j] <- d2L[j,i] <- (-1)*sum(PQ * Thetas[,i] * Thetas[,j])
		for(i in 1:nfact)
			d2L[nfact+1,i] <- d2L[i,nfact+1] <- (-1)*sum(PQ * Thetas[,i])	
		list(grad = dL, hess = d2L)
	} 
	
    #preamble	
	N <- nrow(fulldata)
	J <- ncol(fulldata)
	nfact <- ncol(pars) - 1
	FA <- factor.minres(Rpoly,nfact)
	theta0 <- factor.scores(fulldata,FA$loadings)
    npars <- nfact + 1
	cand.t.var <- 1
	for(i in 1:20){
		theta0 <- draw.thetas(theta0,pars,guess,fulldata,cand.t.var)
		if(attr(theta0,"Proportion Accepted") > .5 && nfact < 5) cand.t.var <- cand.t.var + .1 
		else if(attr(theta0,"Proportion Accepted") > .3) cand.t.var <- cand.t.var + .1 
	}	
	m.thetas <- SEM.stores <- list()	
	phi <- g <- rep(0,J*npars)
	Tau <- info <- h <- matrix(0,J*npars,J*npars)
    m.list <- list()	  
	conv <- 0
    k <- 1 	
	
	for(cycles in 1:max.cycles)
	{
		if(cycles == (SEM.cycles + 1)){
		    pars <- matrix(0, ncol=npars, nrow = J)
			for(i in 1:SEM.cycles) pars <- pars + SEM.stores[[i]]
			pars <- pars/SEM.cycles			
		}
		if(cycles > SEM.cycles){
			gamma <- 1/(cycles - SEM.cycles)        		
			k <- mk
		}		
		
		#Step 1. Generate m_k datasets of theta 
		for(i in 1:k)
			m.thetas[[i]] <- draw.thetas(theta0,pars,guess,fulldata,cand.t.var)
		theta0 <- m.thetas[[1]]
		
		#Step 2. Find average of simulated data gradients and hessian 
		g.m <- h.m <- list()				
		for(j in 1:k){
			for(i in 0:(J - 1)){
			    temp <- dpars(pars[i+1,], guess[i+1], fulldata[,i+1], m.thetas[[j]])
				g[1:npars + i*npars] <- temp$grad
				h[1:npars + i*npars,1:npars + i*npars] <- temp$hess
			} 
			g.m[[j]] <- g
			h.m[[j]] <- h
		}
		ave.g <- rep(0,ncol(fulldata)*npars)
		ave.h <- matrix(0,ncol(fulldata)*npars,ncol(fulldata)*npars)		
		for(i in 1:k){
		  ave.g <- ave.g + g.m[[i]]
		  ave.h <- ave.h + h.m[[i]]
		}
		grad <- ave.g/k
		ave.h <- (-1)*ave.h/k
		if(cycles <= SEM.cycles){
			SEM.stores[[cycles]] <- pars <- pars + 
				matrix(solve(ave.h) %*% grad, ncol=npars, byrow=TRUE)
			next
		}	
		
		#Step 3. Update R-M step
		Tau <- Tau + gamma*(ave.h - Tau)		
		correction <- gamma*(solve(Tau) %*% grad)		
		if(all(correction < tol)) conv <- conv + 1
			else conv <- 0		
		if(conv == 3) break		
		pars <- pars + matrix(correction, ncol = npars, byrow=TRUE)	
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE 	
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	}
	SE <- matrix(sqrt(diag(solve(info))),ncol=npars,byrow=TRUE)
	tmp <- all.P.mirt(theta0, pars, guess, N, J)
	tmp <- ifelse(fulldata,tmp,1-tmp)
	Pl <- apply(tmp,1,prod)
	
	list(pars = pars, SE = SE, cycles = cycles - SEM.cycles, Theta = theta0, Pl = Pl)    
}
 
 
 
 