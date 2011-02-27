##########################################

residuals.bfactor <- function(object, digits = 3, res.cor=FALSE, ...)
{
  if(res.cor) {
    cormat <- object$cormat
    F <- object$F
    Rrep <- F %*% t(F)
	residual <- cormat - Rrep	
	RMR <- 0
	for(i in 1:ncol(cormat))
	  for(j in 1:ncol(cormat))
	    if(i < j) RMR <- RMR + (cormat[i,j] - Rrep[i,j])^2
	RMR <- sqrt(RMR/(ncol(cormat)*(ncol(cormat) - 1) / 2))	
    cat("Residual correlations: \n")
	print(residual,digits)  
	cat("\nRMR : ", round(RMR,3),"\n")
  
  } else {     
    r <- object$tabdata[ ,ncol(object$tabdata)]
    res <- (r - object$Pl * nrow(object$fulldata)) / 
      sqrt(object$Pl * nrow(object$fulldata))
    print(res,digits)
    invisible(res)
  }
}

fitted.bfactor <- function(object, digits = 3, ...)
{  
  expected <- round(object$sampsize * object$Pl,digits)  
  tabdata <- cbind(object$tabdata,expected)
  colnames(tabdata) <- c(object$itemnames, "freq", "exp")	
  print(tabdata)
  invisible(expected)
}

fscores.bfactor <- function(object, full.scores = FALSE, 
  method = "EAP", ...)
{  
  g <- object$guess
  d <- object$pars[ ,ncol(object$pars)]
  a <- as.matrix(object$pars[ ,1:(ncol(object$pars) - 1)])
  nfact <- 2
  theta <- as.matrix(seq(-4,4,length.out = 15))
  Theta <- thetaComb(theta,nfact)
  fulldata <- object$fulldata  
  tabdata <- object$tabdata[ ,1:ncol(fulldata)]
  colnames(tabdata) <- colnames(fulldata) <- object$itemnames
  scores <- matrix(0,ncol=ncol(Theta),nrow=nrow(tabdata))
  thetas <- rep(0,nfact)
  logicalfact <- object$logicalfact
  W <- AXk(0,1,Theta)  
  for (i in 1:nrow(scores)){
    L <- 0  
    for (j in 1:nrow(a)){
      if(tabdata[i,j] == 1) 
	    L <- log(P.bfactor(a[j, ],d[j],Theta,g[j],logicalfact[j, ])) + L
	  else 
        L <- log(1 - P.bfactor(a[j, ],d[j],Theta,g[j],logicalfact[j, ])) + L	
    }	
	for (k in 1:ncol(Theta))
	  thetas[k] <- sum(Theta[ ,k] * exp(L) * W / sum(exp(L) * W))
    scores[i, ] <- thetas
  }
  if(method == "MAP"){
    for (i in 1:nrow(scores)) {       
      Theta <- scores[i, ]	  
      thetas <- nlm(MAP.bfactor,Theta,a=a,d=d,guess=g,
	    patdata=tabdata[i, ],logicalfact=logicalfact)$estimate 
      scores[i, ] <- thetas
    }  
  }  
  scores <- as.matrix(scores[ ,1])
  colnames(scores) <- "g"
  if (full.scores){
    TFvec <- rep(FALSE,nrow(fulldata))  
    scoremat <- matrix(0,nrow=nrow(fulldata),ncol=1) 
    for (j in 1:nrow(tabdata)){
      for (i in 1:nrow(fulldata)){
        TFvec <- rep(FALSE,nrow(fulldata))
        TFvec[i] <- all(fulldata[i, ] == tabdata[j, ])
	    scoremat[TFvec, ] <- scores[j]
      }  
    }    
    return(cbind(fulldata,scoremat))
  } else {  
	cat("Method: ", method,"\n")
    return(cbind(tabdata,scores))
  }   
}  

coef.bfactor <- function(object, digits = 3, ...)
{
  a <- as.matrix(object$pars[ ,1:(ncol(object$pars)-1)])
  d <- object$pars[ ,ncol(object$pars)]
  A <- sqrt(apply(a^2,1,sum))
  B <- -d/A 
  fac <- object$facility  
  parameters <- round(cbind(object$pars,object$guess,fac,A,B),digits)
  colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),"d", "guess", 
    "facility","mvdisc", "mvint")  
  cat("Parameters with multivariate discrimination and intercept: \n")
  print(parameters)	    
  invisible(parameters)
}

summary.bfactor <- function(object, digits = 3, ...)
{
  F <- round(object$F,digits)
  h2 <- round(object$h2,digits)
  SS <- colSums(F^2)
  fac <- as.matrix(object$facility)
  colnames(F) <- c('g',paste("F_", 1:(ncol(F)-1),sep=""))
  names(h2) <- "h2"
  colnames(fac) <- "facility"
  loads <- round(cbind(F,h2,fac),digits)
  rownames(loads) <- object$itemnames  
  cat("Factor loadings: \n")
  print(loads)
  cat("\nSS loadings: ",round(SS,digits), "\n")
  cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
  if(any(h2 > 1)) 
	warning("Solution has heywood cases. Interpret with caution.")
 }

anova.bfactor <- function(object, object2, ...) 
{
  df <- object$df - object2$df
  X2 <- 2*object2$log.lik - 2*object$log.lik 
  AICdiff <- object$AIC - object2$AIC  
  cat("\tChi-squared difference \n\nX2 = ", round(X2,3), ", df = ",
    df, ", p = ", round(1 - pchisq(X2,df),4), "\n")
  cat("AIC difference = ", round(AICdiff,3), "\n")  	
}

print.bfactor <- function(x, ...) 
{
  cat("Call: ") 
  print(x$Call)
  cat("\nFull-information bifactor analysis with ", 
    length(unique(x$specific)), " specific factors \n", sep='')
  if(x$converge == 1)	
    cat("Converged in ", x$EMiter, " iterations.\n", sep="")
  else 	
    cat("Estimation stopped after ", x$EMiter, " iterations.\n", sep="")
  cat("Log-likelihood = ", x$log.lik, "\n")
  cat("AIC = ", x$AIC, "\n")
  cat("Chi-squared = ", round(x$X2,2), ", df = ", 
    x$df, ", p = ", round(x$p,4), "\n")
}

########################################## 

bfactor <- function(fulldata, specific, guess = 0, prev.cor=NULL, par.prior = FALSE,
  startvalues = NULL, quadpts = NULL, ncycles = 50, EMtol=.005, nowarn = TRUE, debug = FALSE, ...)
{ 
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
  if (length(guess) == 1) guess <- rep(guess,nitems)
  else if (length(guess) > nitems || length(guess) < nitems) 
    stop("The number of guessing parameters is incorrect.")
  pats <- apply(fulldata, 1, paste, collapse = "/") 
  freqs <- table(pats)
  nfreqs <- length(freqs)
  r <- as.vector(freqs)
  sampsize <- nrow(fulldata) 
  tabdata <- unlist(strsplit(cbind(names(freqs)), "/"))
  tabdata <- matrix(as.numeric(tabdata), nfreqs, nitems, TRUE)
  tabdata <- cbind(tabdata,r)  
  logicalfact <- matrix(FALSE,nitems,nfact - 1)
  is.na(specific) <- FALSE
  for (i in 1:nitems) logicalfact[i,specific[i]] <- TRUE
  logicalfact <- cbind(rep(TRUE,nitems),logicalfact)     
  if (is.null(quadpts)) quadpts <- 9
  theta <- as.matrix(seq(-4,4,length.out = quadpts))
  Theta <- expand.grid(theta,theta)
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
	} else Rpoly <- tetrachor(fulldata,guess)       
  pars <- matrix(0,nrow=nitems, ncol=nfact + 1)
  if (is.null(startvalues))
  {  
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
  prior <- AXk(0,1,as.matrix(theta))
  Prior <- AXk(0,1,Theta)  
  startvalues <- pars  
  converge <- 1
  index <- 1:nitems
  sitems <- matrix(0,ncol=nitems,nrow=(nfact-1))
  for(i in 1:32) sitems[specific[i],i] <- 1 
  if(debug) {
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
	maxim <- .Call("Mstep",                
                as.double(rlist[[1]]),
                as.double(rlist[[2]]),
				as.double(Prior),
			    as.double(mpars),
			    as.double(guess),    
			    as.double(as.matrix(Theta)),
			    as.double(par.prior))    
    sload <- maxim[(nitems+1):(2*nitems)]		
    for(i in 1:nitems){
	  temp <- selvec[pars[i,2:nfact] > 0]
      pars[i,temp] <- sload[i]
	}  
	pars[ ,1] <- maxim[1:nitems]
	pars[ ,nfact+1] <- maxim[(2*nitems+1):(3*nitems)]	
	if(any(is.na(pars))) converge <- 0
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
    if (cycles %% 3 == 0 & cycles > 6) 
	{
      d1 <- lastpars1 - pars
      d2 <- lastpars2 - pars      
      for (i in 1:nitems) {
        for(j in 1:ncol(pars)){      
          if((abs(d1[i,j]) > 0.001) & (d1[i,j]*d2[i,j] > 0.0) & 
            (d1[i,j]/d2[i,j] < 1.0)) rate[i,j] <- (1 - (1 - rate[i,j]) * (d1[i,j]/d2[i,j]))
		  else rate[i,j] <- 0
        }        
      }      
    }  
    rate[pars > 4] <- 0
	rate[pars < -4] <- 0
	pars <- lastpars1*rate*(-2) + (1 - rate*(-2))*pars          
  }
  
  if(any(par.prior[,1] != 1)) cat("Slope prior for item(s):",as.character(index[par.prior[,1] > 1]), 
    "\n")
  if(any(par.prior[,3] != 0)) cat("Intercept prior for item(s):",as.character(index[par.prior[,3] > 0]), 
    "\n")  
  if(converge == 0) 
    warning("Parameter estimation reached unacceptable values. Model probably did not 
	converged.") 	  	    
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
  norm <- as.matrix(sqrt(1 + pars[ ,1]^2))
  gam <- (-1)*pars[ ,nfact + 1]/norm  
  F <- matrix(0,ncol = nfact, nrow = nitems)
  for (i in 1:nitems) F[i,1:nfact] <- pars[i,1:nfact]/norm[i]  
  h2 <- rowSums(F^2)  
  
  mod <- list(EMiter=cycles, pars=pars, guess=guess, AIC=AIC, X2=X2, df=df, 
    log.lik=log.lik, p=p, F=F, h2=h2, itemnames=itemnames, 
    tabdata=tabdata, sampsize=sampsize, Pl=Pl, Theta=Theta, fulldata=fulldata, 
    logicalfact=logicalfact, facility=facility, specific=specific,
	cormat=Rpoly, converge=converge, par.prior=par.prior,Call=Call) 
    
  class(mod) <- "bfactor"
  mod  
}
