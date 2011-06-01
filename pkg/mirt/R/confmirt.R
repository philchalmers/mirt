coef.confmirt <- function(object, SE = TRUE, print.gmeans = FALSE, digits = 3, ...)
{  
	nfact <- ncol(object$Theta)	
	a <- matrix(object$pars[ ,1:nfact],ncol=nfact)
	d <- matrix(object$pars[,(nfact+1):ncol(object$pars)],
		ncol = ncol(object$pars)-nfact)    	

	parameters <- cbind(object$pars,object$guess)
	SEs <- cbind(object$SEpars,object$SEg)	
	colnames(SEs) <- colnames(parameters) <- c(paste("a_",1:nfact,sep=""),
		paste("d_",1:(ncol(object$pars)-nfact),sep=""),"guess")					
	cat("Item Parameters: \n")
	print(parameters, digits)
	if(SE){
		cat("\nStd. Errors: \n")	
		print(SEs, digits)
	}	
	u <- object$gpars$u	
	sig <- object$gpars$sig	
	cat("\nGroup Parameters: \n")
	if(print.gmeans){
		cat("Means: \n")
		print(u,digits)
		cat("\nStd. Errors: \n")	
		print(object$SEgpars$SEu, digits)	
	}
	cat("Covariance: \n")
	print(sig,digits)
	if(SE){
		cat("\nStd. Errors: \n")	
		print(object$SEgpars$SEsig, digits)	
	}		
}

print.confmirt <- function(x, ...){
	cat("Call: ")
	print(x$Call)
	cat("\nFull-information item factor analysis with ", ncol(x$Theta), " factors \n", sep="")
	if(x$converge == 1)	
		cat("Converged in ", x$cycles, " iterations.\n", sep="")
	else 	
		cat("Estimation stopped after ", x$cycles, " iterations.\n", sep="")	
} 

confmirt <- function(data, sem.mod, guess = 0, gmeans = 0, ncycles = 2000, 
	burnin = 200, SEM.cycles = 100, kdraws = 1, tol = .001, printcycles = TRUE, 
	debug = FALSE, ...){
		
	Call <- match.call()   
	itemnames <- colnames(data)
	data <- as.matrix(data)		
	colnames(data) <- itemnames	
	J <- ncol(data)
	N <- nrow(data)	
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
	Rpoly <- cormod(na.omit(data),K,guess)
	sem.mod <- unclass(sem.mod)
	itemnames <- colnames(data)
	for(i in 1:J){
		tmp <- paste(itemnames[i], "<->", itemnames[i])
		sem.mod <- rbind(sem.mod,c(tmp,paste("th",i,sep=""),NA))		
	}	
	suppressWarnings(SEM <- sem:::sem.mod(sem.mod,Rpoly,N)) 		
	ram <- SEM$ram
	coefs <- SEM$coef
	ramloads <- ram[ram[,1]==1,]
	ramloads[,3] <- ramloads[,3] - J
	groups <- ram[ram[,2] > J,]	
	groups[,2:3] <- groups[,2:3] - J
	nfact <- sum(groups[,2] == groups[,3])	
	if(length(gmeans) == 1)	gmeans <- rep(gmeans,nfact)					
	if(length(gmeans) > J || length(guess) < J) 
		stop("The number of gmeans parameters is incorrect.")
	estgmeans <- gmeans != 0	
	gcov <- selgcov <- estgcov <- matrix(FALSE,nfact,nfact)
	for(i in 1:nfact)
		for(j in 1:nfact)
			if(i <= j) selgcov[j,i] <- TRUE
	est <- is.na(groups[,5])
	for(i in 1:nrow(groups)){ 
		i1 <- groups[i,2]
		i2 <- groups[i,3]
		if(groups[i,4] != 0)
			gcov[i1,i2] <- gcov[i2,i1] <- coefs[groups[i,4]] 
		else 
			gcov[i1,i2] <- gcov[i2,i1] <- groups[i,5] 	
		if(est[i]) estgcov[i2,i1] <- TRUE					
	}	
	loads <- estlam <- matrix(FALSE,J,nfact)
	for(i in 1:nrow(ramloads)){
		item <- ramloads[i,2]
		estlam[item,ramloads[i,3]] <- TRUE
		loads[item,ramloads[i,3]] <- coefs[ramloads[i,4]]
	}		
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
	fulldata[is.na(fulldata)] <- 0
	cs <- sqrt(abs(1-rowSums(loads^2)))
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
	npars <- length(lambdas) + length(zetas) + sum(estGuess) 
	parind <- 1:npars
	pars <- rep(NA,npars)
	Ksum <- cumsum(K + nfact - 1 + estGuess) - (nfact-1)
	sind <- lamind	<- gind <- c()	 
	for(i in 1:J){
		pars[Ksum[i]:(Ksum[i] + nfact - 1)] <- lambdas[i,]
		lamind <- c(lamind,Ksum[i]:(Ksum[i] + nfact - 1))
		sind <- c(sind, rep(TRUE,K[i]-1), estlam[i,])
		if(estGuess[i]){
			pars[Ksum[i] - 2] <- guess[i]
			gind <- c(gind,Ksum[i] - 2)
			sind <- c(sind,TRUE)
		}	
	}
	sind <- c(sind, estgmeans, estgcov[selgcov])
	zetaind <- parind[is.na(pars)]		
	pars[is.na(pars)] <- zetas
	gmeansind <- (length(pars) + 1):(length(pars) + nfact)
	gcovind <- (gmeansind[length(gmeansind)] + 1):(gmeansind[length(gmeansind)] 
		+ nfact*(nfact+1)/2)	
	pars <- c(pars, gmeans, gcov[selgcov]) 
	npars <- length(pars)
	ngpars <- nfact + nfact*(nfact + 1)/2
	converge <- 1    	
	if(debug){
		print(lambdas)
		print(zetas)
		print(guess)
	}	
	
	#preamble for MRHM algorithm			
	theta0 <- matrix(0,N,nfact)	    
	cand.t.var <- 1	
	for(i in 1:20){
		theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var,gcov)
		if(attr(theta0,"Proportion Accepted") > .5) cand.t.var <- cand.t.var + .05 
		else if(attr(theta0,"Proportion Accepted") < .35) cand.t.var <- cand.t.var - .05
		if (cand.t.var < 0)	cand.t.var <- 0.1
	}	
	m.thetas <- grouplist <- list()		
	SEM.stores <- matrix(0,SEM.cycles,npars)
	phi <- rep(0,sum(sind))
	h <- matrix(0,npars,npars)		
	Tau <- info <- matrix(0,sum(sind),sum(sind))
	parind <- 1:npars	
	m.list <- list()	  
	conv <- 0
	k <- 1	
	gamma <- .1
	startvalues <- pars	
	stagecycle <- 1
	
	for(cycles in 1:(ncycles + burnin + SEM.cycles))
	{ 
		if(cycles == burnin + 1) stagecycle <- 2			
		if(stagecycle == 3)
			gamma <- sqrt(0.1/(2*(cycles - SEM.cycles - burnin - 1)))
		if(cycles == (burnin + SEM.cycles + 1)){ 
			stagecycle <- 3		
		    pars <- rep(0,npars)
			for(i in 1:SEM.cycles) pars <- pars + SEM.stores[i,]
			pars <- pars/SEM.cycles				
			k <- kdraws	
			gamma <- 1
		}	
				
		lambdas <- matrix(pars[lamind],J,nfact,byrow=TRUE)
		zetas <- pars[zetaind]
		guess <- rep(0,J)
		guess[estGuess] <- pars[gind]		
		grouplist$u <- pars[gmeansind]
		sig <- matrix(0,nfact,nfact)	
		tmp <- pars[gcovind]
		loc <- 1
		for(i in 1:nfact){
			for(j in 1:nfact){
				if(i <= j) {
					sig[i,j] <- tmp[loc]
					loc <- loc + 1
				}
			}
		}		
		sig <- sig + t(sig) - diag(diag(sig))		
		grouplist$sig <- sig		
		for(j in 1:5) theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)		
		
		#Step 1. Generate m_k datasets of theta 
		for(i in 1:k)
			m.thetas[[i]] <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var,sig)
		theta0 <- m.thetas[[1]]
		
		#Step 2. Find average of simulated data gradients and hessian 		
		g.m <- h.m <- group.m <- list()				
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
			tmp <- d.group(grouplist,m.thetas[[j]])
			g[is.na(g)] <- tmp$g
			h[(npars - ngpars + 1):npars,(npars - ngpars + 1):npars] <- tmp$h
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
		grad <- grad[parind[sind]]		
		ave.h <- ave.h[parind[sind],parind[sind]] 
		if(printcycles){
			if((cycles + 1) %% 10 == 0){
				if(cycles < burnin)
					cat("Stage 1: Cycle = ", cycles + 1, ", Log-Lik = ", 
						attr(theta0,"log.lik"), sep="")
				if(cycles > burnin && cycles < burnin + SEM.cycles)
					cat("Stage 2: Cycle = ", cycles-burnin+1, ", Log-Lik = ",
						attr(theta0,"log.lik"), sep="")
				if(cycles > burnin + SEM.cycles)
					cat("Stage 3: Cycle = ", cycles-burnin-SEM.cycles+1, 
						", Log-Lik = ", attr(theta0,"log.lik"), sep="")
			}
		}			
		if(stagecycle < 3){
		    correction <- solve(ave.h) %*% grad								
			parsold <- pars
			correct <- rep(0,npars)
			correct[sind] <- correction
			pars <- pars + gamma*correct
			if(printcycles && (cycles + 1) %% 10 == 0) 
				cat(", Max Change =", round(max(abs(gamma*correction)),5), "\n")
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
		if(all(gamma*correction < tol)) conv <- conv + 1
			else conv <- 0		
		if(conv == 3) break		
		parsold <- pars
		correct <- rep(0,npars)
		correct[sind] <- correction
		pars <- pars + gamma*correct
		if(printcycles && (cycles + 1) %% 10 == 0) 
			cat(", Max Change =", round(max(abs(gamma*correction)),5), "\n")
		pars[pars[gind] < 0] <- parsold[pars[gind] < 0]
		pars[pars[gind] > .4] <- parsold[pars[gind] > .4]
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE 			
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	}
	cat("\n")
	SEtmp <- diag(solve(info))
	if(any(SEtmp < 0)){
		warning("Information matrix is not positive definite, negative SEs set to 'NA'.\n")
		SEtmp[SEtmp < 0] <- NA
	}
	if(any(guess < 0)) warning("Negative lower asymptote parameter(s). \n")		
	SEtmp <- sqrt(SEtmp)
	SE <- rep(NA,npars) 
	SE[parind[sind]] <- SEtmp
	estpars <- pars[sind]
	lambdas <- matrix(pars[lamind],J,nfact,byrow=TRUE)	
	lambdas[!estlam] <- NA	
	guess <- rep(NA,J)
	guess[estGuess] <- pars[gind]
	zetas <- pars[zetaind]
	u <- pars[gmeansind]	
	sig <- matrix(0,nfact,nfact)
	SElam <- matrix(SE[lamind],J,nfact,byrow=TRUE)
	SEzetas <- SE[zetaind]	
	SEg <- rep(NA,J)	
	SEg[estGuess] <- SE[gind]	
	SEu <- SE[gmeansind]	
	SEsig <- matrix(0,nfact,nfact)	
	tmp <- pars[gcovind]
	tmp2 <- SE[gcovind]
	loc <- 1
	for(i in 1:nfact){
		for(j in 1:nfact){
			if(i <= j) {
				sig[i,j] <- tmp[loc]
				SEsig[i,j] <- tmp2[loc]
				loc <- loc + 1
			}
		}
	}		
	sig <- sig + t(sig) - diag(diag(sig))		
	SEsig <- SEsig + t(SEsig) - diag(diag(SEsig))	
	tmp1 <- tmp2 <- matrix(NA,J,(max(K)-1))
	loc <- 1
	for(i in 1:J){
		for(j in 1:(K[i]-1)){
			tmp1[i,j] <- zetas[loc] 
			tmp2[i,j] <- SEzetas[loc]
			loc <- loc + 1
		}
	}	 
	zetas <- tmp1
	SEzetas <- tmp2	
	pars <- cbind(lambdas,zetas)
	SEpars <- cbind(SElam,SEzetas)
	gpars <- list(u = u, sig = sig)	
	SEgpars <- list(SEu = SEu, SEsig = SEsig)
	estpars <- list(estlam=estlam,estGuess=estGuess,estgcov=estgcov,
		estgmeans=estgmeans)
	mod <- list(pars=pars, guess=guess, SEpars=SEpars, SEg = SEg, gpars=gpars, 
		SEgpars=SEgpars, estpars=estpars,cycles=cycles - SEM.cycles - burnin,
		Theta=theta0, fulldata=fulldata, K=K, itemloc=itemloc, 
		converge = converge, Call=Call)	 
	class(mod) <- 'confmirt'
	mod
}	
