# quadrature weights
AXk <- function(mu, sigma, Theta) {   
  A <- matrix(0, ncol = ncol(Theta), nrow = nrow(Theta))
  for (i in 1:ncol(Theta)) {  
    A[ ,i] <- exp(-0.5*((Theta[ ,i] - mu)/sigma)^2)    
  }  
  A <- apply(A,1,prod)  
  return(A/sum(A))  
}

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
 