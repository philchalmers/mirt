setMethod(
	f = "print",
	signature = signature(x = 'confmirtClass'),
	definition = function(x, ...){
		cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information item factor analysis with ", ncol(x@Theta), " factors \n", sep="")
		if(length(x@logLik) > 0){
			cat("Log-likelihood = ", x@logLik,", SE = ",round(x@SElogLik,3), "\n",sep='')			
			cat("df =", x@df, "\nAIC = ", x@AIC, "\n")
		}
		if(x@converge == 1)	
			cat("Converged in ", x@cycles, " iterations.\n", sep="")
		else 	
			cat("Estimation stopped after ", x@cycles, " iterations.\n", sep="")	
	} 
)

setMethod(
	f = "show",
	signature = signature(object = 'confmirtClass'),
	definition = function(object){
		cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information item factor analysis with ", ncol(object@Theta), " factors \n", sep="")
		if(length(object@logLik) > 0){
			cat("Log-likelihood = ", object@logLik,", SE = ",round(object@SElogLik,3), "\n",sep='')			
			cat("df =", object@df, "\nAIC = ", object@AIC, "\n")
		}
		if(object@converge == 1)	
			cat("Converged in ", object@cycles, " iterations.\n", sep="")
		else 	
			cat("Estimation stopped after ", object@cycles, " iterations.\n", sep="")	
	} 
)

setMethod(
	f = "summary",
	signature = 'confmirtClass',
	definition = function(object, digits = 3, ...)
	{
		nfact <- ncol(object@F)		
		F <- object@F
		h2 <- as.matrix(object@h2)		
		colnames(F) <- paste("F_", 1:ncol(F),sep="")
		colnames(h2) <- "h2"				
		SS <- apply(F^2,2,sum)			
		cat("\nFactor loadings metric: \n")
		print(cbind(F,h2),digits)		
		cat("\nSS loadings: ",round(SS,digits), "\n")		
		cat("\nFactor correlations: \n")
		Phi <- cov2cor(object@gpars$sig)	  
		Phi <- round(Phi, digits)
		colnames(Phi) <- rownames(Phi) <- colnames(F)
		print(Phi)		
		if(any(h2 > 1)) 
			warning("Solution has heywood cases. Interpret with caution.") 
		invisible(F)  	  
	}
)

setMethod(
	f = "coef",
	signature = 'confmirtClass',
	definition = function(object, SE = TRUE, print.gmeans = FALSE, digits = 3, ...)
	{  
		nfact <- ncol(object@Theta)	
		a <- matrix(object@pars[ ,1:nfact],ncol=nfact)
		d <- matrix(object@pars[,(nfact+1):ncol(object@pars)],
			ncol = ncol(object@pars)-nfact)    	

		parameters <- cbind(object@pars,object@guess)
		SEs <- cbind(object@SEpars,object@SEg)	
		colnames(SEs) <- colnames(parameters) <- c(paste("a_",1:nfact,sep=""),
			paste("d_",1:(ncol(object@pars)-nfact),sep=""),"guess")					
		cat("ITEM PARAMTERS: \n")
		print(parameters, digits)
		if(SE){
			cat("\nStd. Errors: \n")	
			print(SEs, digits)
		}	
		u <- object@gpars$u	
		sig <- object@gpars$sig	
		cat("\nGROUP PARAMETERS: \n")
		if(print.gmeans){
			cat("Means: \n")
			print(u,digits)
			cat("\nStd. Errors: \n")	
			print(object@SEgpars$SEu, digits)	
		}
		cat("Covariance: \n")
		print(sig,digits)
		if(SE){
			cat("\nStd. Errors: \n")	
			print(object@SEgpars$SEsig, digits)	
		}		
	}
)

setMethod(
	f = "residuals",
	signature = signature(object = 'confmirtClass'),
	definition = function(object, digits = 3, ...)
	{ 
		fulldata <- object@fulldata	
		data <- object@data
		data[data==99] <- NA
		N <- nrow(fulldata)
		K <- object@K
		J <- length(K)
		sig <- object@gpars$sig	
		nfact <- ncol(object@F)
		theta <- seq(-4,4, length.out = round(20/nfact))
		Theta <- thetaComb(theta,nfact)		
		lambdas <- matrix(object@pars[,1:nfact], J)
		lambdas[is.na(lambdas)] <- 0
		zetas <- as.vector(t(object@pars[,(nfact+1):ncol(object@pars)]))
		zetas <- na.omit(zetas)
		guess <- object@guess
		guess[is.na(guess)] <- 0	
		Ksums <- cumsum(K) - 1	
		itemloc <- object@itemloc
		res <- matrix(0,J,J)
		diag(res) <- NA
		colnames(res) <- rownames(res) <- colnames(data)
		prior <- dmvnorm(Theta,rep(0,nfact),sig)
		prior <- prior/sum(prior)
		loc <- loc2 <- 1	
		for(i in 1:J){
			if(i > 1) loc <- loc + K[i-1] - 1	
			loc2 <- 1
			for(j in 1:J){			
				if(i < j){
					if(K[i] > 2) P1 <- P.poly(lambdas[i,],zetas[loc:(loc+K[i]-2)],Theta,itemexp=TRUE)
					else { 
						P1 <- P.mirt(lambdas[i,],zetas[loc], Theta, guess[i])
						P1 <- cbind(1 - P1, P1)
					}	
					if(K[j] > 2) P2 <- P.poly(lambdas[j,],zetas[loc2:(loc2+K[j]-2)],Theta,itemexp=TRUE)
					else {
						P2 <- P.mirt(lambdas[j,],zetas[loc2], Theta, guess[j])	
						P2 <- cbind(1 - P2, P2)
					}
					tab <- table(data[,i],data[,j])		
					Etab <- matrix(0,K[i],K[j])
					for(k in 1:K[i])
						for(m in 1:K[j])						
							Etab[k,m] <- N * sum(P1[,k] * P2[,m] * prior)	
					s <- gamma.cor(tab) - gamma.cor(Etab)
					if(s == 0) s <- 1				
					res[j,i] <- sum(((tab - Etab)^2)/Etab) * sign(s)
					res[i,j] <- sqrt( abs(res[j,i]) / (N * min(c(K[i],K[j]) - 1))) 					
				}
			loc2 <- loc2 + K[j] - 1 	
			}
		}		
		cat("LD matrix:\n\n")	
		print(res,digits)    	
	}
)

setMethod(
	f = "logLik",
	signature = signature(object = 'confmirtClass'),
	definition = function(object, draws = 2000){	
		nfact <- ncol(object@Theta)
		N <- nrow(object@Theta)
		J <- length(object@K)
		pars <- object@pars
		lambdas <- pars[,1:nfact]
		lambdas[is.na(lambdas)] <- 0
		zetas <- pars[,(nfact+1):ncol(pars)]
		zetas <- t(zetas)[!is.na(t(zetas))]		
		mu <- object@gpars$u
		sigma <- object@gpars$sig		
		LL <- matrix(0,N,draws)
		theta <- matrix(0,N,nfact*draws)
		guess <- object@guess
		guess[is.na(guess)] <- 0
		K <- object@K
		df <- object@df	
		for(i in 1:draws){
			theta <- rmvnorm(N,mu,sigma)				
			LL[,i] <- .Call('logLik', 					
						as.numeric(lambdas),
						as.numeric(zetas),
						as.numeric(guess),
						as.numeric(theta),
						as.integer(object@fulldata),
						as.integer(object@itemloc-1),
						as.integer(object@K),
						as.integer(J),
						as.integer(N),
						as.integer(nfact))		
		}
		logLik <- sum(log(rowMeans(LL)))
		SElogLik <- sqrt(var(log(rowMeans(LL))) / draws)
		AIC <- (-2) * logLik + 2 * (prod(K) - df) 
		object@logLik <- logLik
		object@SElogLik <- SElogLik		
		object@AIC <- AIC		
		return(object)
	} 	
)

setMethod(
	f = "anova",
	signature = signature(object = 'confmirtClass'),
	definition = function(object, object2, ...){
		dots <- list(...)				
		nitems <- length(object@K)
		if(length(object@df) == 0 || length(object2@df) == 0) 
			stop('Use \'logLik\' to obtain likelihood values') 	
		df <- object@df - object2@df 
		if(df < 0){
			df <- abs(df)
			tmp <- object
			object <- object2
			object2 <- tmp
		}
		X2 <- 2*object2@logLik - 2*object@logLik 
		AICdiff <- object@AIC - object2@AIC  
		se <- round(object@SElogLik + object2@SElogLik,3)		
		cat("\nChi-squared difference: \n\nX2 = ", round(X2,3), 
			" (SE = ",se,"), df = ", df, ", p = ", round(1 - pchisq(X2,df),4), "\n", sep="")
		cat("AIC difference = ", round(AICdiff,3)," (SE = ", se,")\n", sep='')
	}		
)

####################
#Main Function

confmirt <- function(data, sem.model, guess = 0, gmeans = 0, ncycles = 2000, 
	burnin = 100, SEM.cycles = 50, kdraws = 1, tol = .001, printcycles = TRUE, 
	calcLL = TRUE, draws = 2000, debug = FALSE, ...){
		
	Call <- match.call()   
	itemnames <- colnames(data)
	data <- as.matrix(data)		
	colnames(data) <- itemnames	
	J <- ncol(data)
	N <- nrow(data)	
	if(length(guess) == 1) guess <- rep(guess,J)
	if(length(guess) > J || length(guess) < J) 
		stop("The number of guessing parameters is incorrect.")					
	uniques <- list()
	for(i in 1:J)
		uniques[[i]] <- sort(unique(data[,i]))
	K <- rep(0,J)
	for(i in 1:J) K[i] <- length(uniques[[i]])	
	guess[K > 2] <- 0
	estGuess <- guess > 0	
	Rpoly <- cormod(na.omit(data),K,guess)
	sem.model <- unclass(sem.model)
	itemnames <- colnames(data)
	for(i in 1:J){
		tmp <- paste(itemnames[i], "<->", itemnames[i])
		sem.model <- rbind(sem.model,c(tmp,paste("th",i,sep=""),NA))		
	}	
	SEM <- sem.mod(sem.model,Rpoly,N)
	ram <- SEM$ram	
	coefs <- rep(.5,nrow(ram))
	ramloads <- ram[ram[,1]==1,]
	ramloads[,3] <- ramloads[,3] - J
	constvalues <- unique(ramloads[,4])[table(ramloads[,4]) > 1]
	nconstvalues <- 0
	if(length(constvalues) > 0)
		for(i in 1:length(constvalues))
			nconstvalues <- nconstvalues + sum(ramloads[,4] == constvalues[i]) - 1
	groups <- matrix(ram[ram[,2] > J,],ncol=5)	
	groups[,2:3] <- groups[,2:3] - J	
	nfact <- sum(groups[,2] == groups[,3])		
	if(length(gmeans) == 1)	gmeans <- rep(gmeans,nfact)					
	if(length(gmeans) > J || length(guess) < J) 
		stop("The number of gmeans parameters is incorrect.")
	estgmeans <- gmeans != 0	
	gcov <- selgcov <- estgcov <- matrix(FALSE,nfact,nfact)
	selgcov <- lower.tri(selgcov, diag = TRUE)	
	est <- is.na(groups[,5])
	for(i in 1:nrow(groups)){ 
		i1 <- groups[i,2]
		i2 <- groups[i,3]
		if(groups[i,4] != 0)
			gcov[i1,i2] <- gcov[i2,i1] <- .1
		else 
			gcov[i1,i2] <- gcov[i2,i1] <- groups[i,5] 	
		if(est[i]) estgcov[i1,i2] <- estgcov[i2,i1] <- TRUE					
	}	
	estgcov <- (estgcov + selgcov) == 2
	loads <- tmplambdas <- matrix(0,J,nfact)
	estlam <- matrix(FALSE,J,nfact)
	constlam <- matrix(0,J,nfact)
	for(i in 1:nrow(ramloads)){
		item <- ramloads[i,2]
		if(!is.na(ramloads[i,5])){
			tmplambdas[item,ramloads[i,3]] <- ramloads[i,5]
		} else {	
			estlam[item,ramloads[i,3]] <- TRUE			
			loads[item,ramloads[i,3]] <- coefs[ramloads[i,4]]
			if(any(ramloads[i,4] == constvalues)) 
				constlam[item,ramloads[i,3]] <- constvalues[ramloads[i,4] == constvalues]
		}
	}
	constlam <- as.vector(t(constlam))
	itemloc <- cumsum(c(1,K))
	index <- 1:J	
	fulldata <- fulldata2 <- matrix(0,N,sum(K))
	Names <- NULL
	for(i in 1:J)
        Names <- c(Names, paste("Item.",i,"_",1:K[i],sep=""))				
	colnames(fulldata) <- Names			
	for(i in 1:J){
		ind <- index[i]
		if(setequal(uniques[[i]], c(0,1))){
			fulldata[ ,itemloc[ind]:(itemloc[ind]+1)] <- cbind(data[,ind],abs(1-data[,ind]))
			fulldata2[ ,itemloc[ind]:(itemloc[ind]+1)] <- cbind(abs(1-data[,ind]),data[,ind])
			next
		}
		dummy <- matrix(0,N,K[ind])
		for (j in 0:(K[ind]-1))  
			dummy[,j+1] <- as.integer(data[,ind] == uniques[[ind]][j+1])  		
		fulldata[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- dummy
		fulldata2[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- dummy	
	}	
	fulldata[is.na(fulldata)] <- fulldata2[is.na(fulldata2)] <- 0
	cs <- sqrt(abs(1-rowSums(loads^2)))	
	lambdas <- loads + tmplambdas
	zetas <- rep(0,ncol(fulldata) - J)	
	loc <- 1	
	for(i in 1:J){
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
	npars <- length(lambdas) + length(zetas) + sum(estGuess) 
	parind <- 1:npars
	pars <- constrained <- rep(NA,npars)
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
	tmp <- .05
	for(i in 1:30){			
		theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var,gcov)
		if(i > 5){		
			if(attr(theta0,"Proportion Accepted") > .35) cand.t.var <- cand.t.var + tmp 
			else if(attr(theta0,"Proportion Accepted") > .25 && nfact > 3) cand.t.var <- cand.t.var + tmp	
			else if(attr(theta0,"Proportion Accepted") < .2 && nfact < 4) cand.t.var <- cand.t.var - tmp
			else if(attr(theta0,"Proportion Accepted") < .1) cand.t.var <- cand.t.var - tmp
			if (cand.t.var < 0){
				cand.t.var <- tmp		
				tmp <- tmp / 2
			}		
		}
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
	gamma <- .25
	startvalues <- pars	
	stagecycle <- 1		
	
	for(cycles in 1:(ncycles + burnin + SEM.cycles))		
	{ 
		if(cycles == burnin + 1) stagecycle <- 2			
		if(stagecycle == 3)
			gamma <- (0.05/(cycles - SEM.cycles - burnin - 1))^(0.5) - .004
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
		
		#Step 1. Generate m_k datasets of theta 
		for(j in 1:4) theta0 <- draw.thetas(theta0,lambdas,zetas,guess,
			fulldata,K,itemloc,cand.t.var,sig)	
		for(i in 1:k) m.thetas[[i]] <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,
			K,itemloc,cand.t.var,sig)
		theta0 <- m.thetas[[1]]
		
		#Step 2. Find average of simulated data gradients and hessian 		
		g.m <- h.m <- group.m <- list()				
		for(j in 1:k){
			g <- rep(NA,npars)
			loc <- 1
			for(i in 0:(J - 1)){
				if(estGuess[i+1]){
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
						fulldata2[,itemloc[i+1]:(itemloc[i+2]-1)],m.thetas[[j]])
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
						sprintf("%.1f",attr(theta0,"log.lik")), sep="")
				if(cycles > burnin && cycles < burnin + SEM.cycles)
					cat("Stage 2: Cycle = ", cycles-burnin+1, ", Log-Lik = ",
						sprintf("%.1f",attr(theta0,"log.lik")), sep="")
				if(cycles > burnin + SEM.cycles)
					cat("Stage 3: Cycle = ", cycles-burnin-SEM.cycles+1, 
						", Log-Lik = ", sprintf("%.1f",attr(theta0,"log.lik")), sep="")					
			}
		}			
		if(stagecycle < 3){			
			correction <- SparseM::solve(ave.h) %*% grad
			correction[correction > 1] <- 1
			correction[correction < -1] <- -1								
			parsold <- pars
			correct <- rep(0,npars)
			correct[sind] <- correction	
			if(any(estGuess))
				correct[gind] <- 0
			for(i in 1:length(constvalues)){
				tmp <- correct[lamind]
				tmp[constlam == constvalues[i]] <- 
					mean(tmp[constlam == constvalues[i]])
				correct[lamind] <- tmp
			}						
			pars <- pars + gamma*correct
			if(printcycles && (cycles + 1) %% 10 == 0){ 
				cat(", Max Change =", sprintf("%.4f",max(abs(gamma*correction))), "\n")
				flush.console()
			}				
			pars[pars[gcovind] > 1] <- parsold[pars[gcovind] > 1]
			pars[pars[gcovind] < -1] <- parsold[pars[gcovind] < -1]			
			pars[pars[gcovind] > 1] <- parsold[pars[gcovind] > 1]
			pars[pars[gcovind] < -1] <- parsold[pars[gcovind] < -1]		
			if(stagecycle == 2) SEM.stores[cycles - burnin,] <- pars
			next
		}	 
		
		#Step 3. Update R-M step		
		Tau <- Tau + gamma*(ave.h - Tau)			
		correction <- SparseM::solve(Tau) %*% grad												
		parsold <- pars
		correct <- rep(0,npars)
		correct[sind] <- correction
		if(any(estGuess))
			correct[gind] <- 0
		for(i in 1:length(constvalues)){
			tmp <- correct[lamind]
			tmp[constlam == constvalues[i]] <- 
				mean(tmp[constlam == constvalues[i]])
			correct[lamind] <- tmp
		}
		if(all(gamma*correct < tol)) conv <- conv + 1
			else conv <- 0		
		if(conv == 3) break	
		pars <- pars + gamma*correct
		if(printcycles && (cycles + 1) %% 10 == 0){ 
			cat(", gam = ",sprintf("%.3f",gamma),", Max Change = ", 
				sprintf("%.4f",max(abs(gamma*correction))), "\n", sep = '')
			flush.console()		
		}	
		pars[pars[gind] < 0] <- parsold[pars[gind] < 0]
		pars[pars[gind] > .4] <- parsold[pars[gind] > .4]
		pars[pars[gcovind] > 1] <- parsold[pars[gcovind] > 1]
		pars[pars[gcovind] < -1] <- parsold[pars[gcovind] < -1]
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE 			
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	}
	cat("\n\n")
	SEtmp <- diag(solve(info))		
	if(any(SEtmp < 0)){
		warning("Information matrix is not positive definite, negative SEs set to 'NA'.\n")
		SEtmp[SEtmp < 0] <- NA
	}	
	SEtmp <- sqrt(SEtmp)	
	SE <- rep(NA,npars) 
	SE[parind[sind]] <- SEtmp
	for(i in 1:length(constvalues)){
		tmp <- SE[lamind]
		tmp[constlam == constvalues[i]] <- 
			mean(tmp[constlam == constvalues[i]])
		SE[lamind] <- tmp
	}
	estpars <- pars[sind]
	lambdas <- matrix(pars[lamind],J,nfact,byrow=TRUE)	
	lambdas[!estlam] <- NA	
	guess <- rep(NA,J)
	guess[estGuess] <- pars[gind]
	guess[K == 2 & !estGuess] <- 0
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
	if(nfact > 1) SEsig <- SEsig + t(SEsig) - diag(diag(SEsig))	
		else SEsig <- NA
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
	df <- as.integer(prod(K) - sum(estlam) - sum(estgcov) - 
		sum(estgmeans) - sum(!is.na(zetas)) - nconstvalues - 1)		
		
	if (nfact > 1) norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2,na.rm = TRUE))
		else norm <- as.matrix(sqrt(1 + pars[ ,1]^2))  
	F <- as.matrix(pars[ ,1:nfact]/norm)
	F[is.na(F)] <- 0	
	h2 <- rowSums(F^2)	

	mod <- new('confmirtClass', pars=pars, guess=guess, SEpars=SEpars, SEg = SEg, 
		gpars=gpars, SEgpars=SEgpars, estpars=estpars,cycles=cycles - SEM.cycles 
		- burnin, Theta=theta0, fulldata=fulldata, data=data, K=K, itemloc=itemloc, 
		h2=h2,F=F,converge = converge, df = df, Call=Call)
	if(calcLL){
		cat("Calculating log-likelihood...\n")
		flush.console()
		mod <- logLik(mod,draws)		
	}	
	return(mod)
}	
