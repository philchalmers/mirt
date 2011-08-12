setMethod(
	f = "print",
	signature = signature(x = 'confmirtClass'),
	definition = function(x, ...)
	{
		cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information item factor analysis with ", ncol(x@Theta), " factors \n", sep="")
		if(length(x@logLik) > 0){
			cat("Log-likelihood = ", x@logLik,", SE = ",round(x@SElogLik,3), "\n",sep='')			
			cat("AIC =", x@AIC, "\n")			
			cat("BIC =", x@BIC, "\n")
			if(x@p < 1)
				cat("G^2 = ", round(x@G2,2), ", df = ", 
					x@df, ", p = ", round(x@p,4), "\n", sep="")
			else 
				cat("G^2 = ", NA, ", df = ", 
					x@df, ", p = ", NA, "\n", sep="")		
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
	definition = function(object)
	{
		cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information item factor analysis with ", ncol(object@Theta), " factors \n", sep="")
		if(length(object@logLik) > 0){
			cat("Log-likelihood = ", object@logLik,", SE = ",round(object@SElogLik,3), "\n",sep='')
			cat("AIC =", object@AIC, "\n")			
			cat("BIC =", object@BIC, "\n")
			if(object@p < 1)	
				cat("G^2 = ", round(object@G2,2), ", df = ", 
					object@df, ", p = ", round(object@p,4), "\n", sep="")
			else 
				cat("G^2 = ", NA, ", df = ", 
					object@df, ", p = ", NA, "\n", sep="")
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
		if(any(object@estComp)) stop('No factor metric for non-compensatory models')
		nfact <- ncol(object@F)		
		F <- object@F		
		colnames(F) <- paste("F_", 1:ncol(F),sep="")						
		SS <- apply(F^2,2,sum)			
		cat("\nFactor loadings metric: \n")
		print(cbind(F),digits)		
		cat("\nSS loadings: ",round(SS,digits), "\n")		
		cat("\nFactor correlations: \n")
		Phi <- cov2cor(object@gpars$sig)	  
		Phi <- round(Phi, digits)
		colnames(Phi) <- rownames(Phi) <- colnames(F)
		print(Phi)				
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
		cat("\nITEM PARAMTERS: \n")
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
	definition = function(object, restype = 'LD', digits = 3, ...)
	{ 
		fulldata <- object@fulldata	
		data <- object@data
		data[data==99] <- NA		
		N <- nrow(data)
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
		zetalist <- list()
		loc <- 1
		for(i in 1:J){
			zetalist[[i]] <- zetas[loc:(loc+K[i]-2)]
			loc <- loc + K[i] - 1		
		}
		guess <- object@guess
		guess[is.na(guess)] <- 0	
		Ksums <- cumsum(K) - 1	
		itemloc <- object@itemloc
		res <- matrix(0,J,J)
		diag(res) <- NA
		colnames(res) <- rownames(res) <- colnames(data)
		prior <- dmvnorm(Theta,rep(0,nfact),sig)
		prior <- prior/sum(prior)		
		if(restype == 'LD'){	
			for(i in 1:J){				
				for(j in 1:J){			
					if(i < j){
						if(K[i] > 2) P1 <- P.poly(lambdas[i,],zetalist[[i]],Theta,itemexp=TRUE)
						else { 
							P1 <- P.mirt(lambdas[i,],zetalist[[i]], Theta, guess[i])
							P1 <- cbind(1 - P1, P1)
						}	
						if(K[j] > 2) P2 <- P.poly(lambdas[j,],zetalist[[j]],Theta,itemexp=TRUE)
						else {
							P2 <- P.mirt(lambdas[j,],zetalist[[j]], Theta, guess[j])	
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
						res[i,j] <- sqrt(abs(res[j,i]) / (N * min(c(K[i],K[j]) - 1))) 					
					}					
				}
			}		
			cat("LD matrix:\n\n")	
			res <- round(res,digits)    	
			return(res)
		}
		if(restype == 'exp'){
			if(length(object@tabdata) == 0) stop('Expected response vectors cannot be computed because logLik() 
				has not been run or the data contains missing responses.')
			tabdata <- object@tabdata
			res <- (tabdata[,J+1] - tabdata[,J+2]) / sqrt(tabdata[,J+2])
			tabdata <- round(cbind(tabdata,res),digits)
			colnames(tabdata) <- c(colnames(object@data), 'freq', 'exp', 'std_res')
			return(tabdata)
		}
	}
)

setMethod(
	f = "logLik",
	signature = signature(object = 'confmirtClass'),
	definition = function(object, draws = 2000, G2 = TRUE)
	{	
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
		guess <- object@guess
		guess[is.na(guess)] <- 0
		K <- object@K	
		fulldata <- object@fulldata	
		for(i in 1:draws){
			theta <- rmvnorm(N,mu,sigma)				
			LL[,i] <- .Call('logLik', 					
						as.numeric(lambdas),
						as.numeric(zetas),
						as.numeric(guess),
						as.numeric(theta),
						as.integer(fulldata),
						as.integer(object@itemloc-1),
						as.integer(object@K),
						as.integer(J),
						as.integer(N),
						as.integer(nfact),
						as.integer(object@estComp))		
		}		
		rwmeans <- rowMeans(LL)
		logLik <- sum(log(rwmeans))				
		pats <- apply(fulldata,1,paste,collapse = "/")
		freqs <- table(pats)
		nfreqs <- length(freqs)		
		r <- as.vector(freqs)
		ncolfull <- ncol(fulldata)
		tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
		tabdata <- matrix(as.numeric(tabdata),nfreqs,ncolfull,TRUE)
		tabdata <- cbind(tabdata,r)		 		 				
		pats <- apply(fulldata,1,paste,collapse = "/")
		freqs <- table(pats)			
		r <- as.vector(freqs)
		logN <- 0
		logr <- rep(0,length(r))
		for (i in 1:N) logN <- logN + log(i)
		for (i in 1:length(r)) 
			for (j in 1:r[i]) 
				logr[i] <- logr[i] + log(j)    		
		if(sum(logr) != 0)		
			logLik <- logLik + logN/sum(logr)			
		SElogLik <- sqrt(var(log(rowMeans(LL))) / draws)
		x <- object@estpars	
		df <- as.integer(length(r) - sum(x$estlam) - sum(x$estgcov) - 
			sum(x$estgmeans) - length(zetas) + object@nconstvalues + 
			nfact*(nfact - 1)/2 - 1)			
		AIC <- (-2) * logLik + 2 * (length(r) - df - 1)
		BIC <- (-2) * logLik + (length(r) - df - 1)*log(N)
		if(G2){			
			data <- object@data
			if(any(is.na(data))){
				object@G2 <- 0	
				object@p <- 1					
			} else {			
				pats <- apply(data,1,paste,collapse = "/")			
				freqs <- table(pats)
				nfreqs <- length(freqs)		
				r <- as.vector(freqs)
				ncolfull <- ncol(data)
				tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
				tabdata <- matrix(as.numeric(tabdata),nfreqs,ncolfull,TRUE)
				tabdata <- cbind(tabdata,r)	
				expected <- rep(0,nrow(tabdata))	
				for (j in 1:nrow(tabdata)){          
					TFvec <- colSums(ifelse(t(data) == tabdata[j,1:ncolfull],1,0)) == ncolfull        
					expected[j] <- mean(rwmeans[TFvec])
					rwmeans[TFvec] <- rwmeans[TFvec]/r[j]
				}
				tabdata <- cbind(tabdata,expected*N)
				G2 <- 2 * sum(log(1/(N*rwmeans)))
				p <- 1 - pchisq(G2,df) 
				object@G2 <- G2	
				object@p <- p
				object@tabdata <- tabdata
			}	
		}	
		object@logLik <- logLik
		object@SElogLik <- SElogLik		
		object@AIC <- AIC
		object@BIC <- BIC
		object@df <- df		
		return(object)
	} 	
)

setMethod(
	f = "anova",
	signature = signature(object = 'confmirtClass'),
	definition = function(object, object2, ...)
	{
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
		BICdiff <- object@BIC - object2@BIC  
		se <- round(object@SElogLik + object2@SElogLik,3)		
		cat("\nChi-squared difference: \n\nX2 = ", round(X2,3), 
			" (SE = ",se,"), df = ", df, ", p = ", round(1 - pchisq(X2,df),4), "\n", sep="")
		cat("AIC difference = ", round(AICdiff,3)," (SE = ", se,")\n", sep='')
		cat("BIC difference = ", round(BICdiff,3)," (SE = ", se,")\n", sep='')
	}		
)

####################
#Main Function

confmirt <- function(data, model, guess = 0, ncycles = 2000, 
	burnin = 150, SEM.cycles = 50, kdraws = 1, tol = .001, printcycles = TRUE, 
	calcLL = TRUE, draws = 2000, returnindex = FALSE, debug = FALSE, ...)
{		
	Call <- match.call()   
	itemnames <- colnames(data)
	keywords <- c('SLOPE','INT','COV','MEAN','COMP')
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
	itemloc <- cumsum(c(1,K))	
	model <- matrix(model$x,ncol=2)
	factorNames <- setdiff(model[,1],keywords) 
	nfact <- length(factorNames)	
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
		
	#slopes specification
	estlam <- matrix(FALSE, ncol = nfact, nrow = J)	
	for(i in 1:nfact){
		tmp <- model[model[,1] == factorNames[i],2]
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

	#COMP
	estComp <- rep(FALSE,J)
	if(any(model[,1] == 'COMP')){
		tmp <- model[model[,1] == 'COMP',2]		
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
	lamind <- zetaind <- guessind <- sind <- c()	
	for(i in 1:J){
		if(estComp[i])
			zetaind <- c(zetaind, loc1:(loc1+(length(estzetas[[i]])-1)))		
		else zetaind <- c(zetaind,loc1:(loc1 + K[i] - 2))
		lamind <- c(lamind,max(zetaind + 1):(max(zetaind)+nfact))		
		guessind <- c(guessind,max(lamind + 1):max(lamind + 1 ))
		sind <- c(sind, estzetas[[i]], estlam[i,], estGuess[i])
		loc1 <- loc1 + nfact + sum(estzetas[[i]]) + 1	
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
	pars[groupind] <- c(gmeans,gcov[lower.tri(gcov,diag=TRUE)])
	parcount <- list(lam = estlam, zeta = estzetas2, guess = estGuess, cov = estgcov, mean = estgmeans)
	parind <- 1:npars
	loc1 <- 1
	parcount$lam <- matrix(lamind,J,byrow=TRUE)
	parcount$zeta <- zetaind
	parcount$guess <- guessind
	parcount$cov <- covind
	parcount$mean <- meanind
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
	if(any(model[,1] == 'SLOPE')){
		tmp <- model[model[,1] == "SLOPE",2]
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
	if(returnindex) return(parcount)	
		
	#Preamble for MRHM algorithm
	pars[constvalues[,1] == 1] <- constvalues[constvalues[,1] == 1,2]
	theta0 <- matrix(0,N,nfact)	    
	cand.t.var <- 1			
	tmp <- .1
	for(i in 1:30){			
		theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var,gcov,gmeans,estComp)
		if(i > 5){		
			if(attr(theta0,"Proportion Accepted") > .35) cand.t.var <- cand.t.var + 2*tmp 
			else if(attr(theta0,"Proportion Accepted") > .25 && nfact > 3) cand.t.var <- cand.t.var + tmp	
			else if(attr(theta0,"Proportion Accepted") < .2 && nfact < 4) cand.t.var <- cand.t.var - tmp
			else if(attr(theta0,"Proportion Accepted") < .1) cand.t.var <- cand.t.var - 2*tmp
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
	m.list <- list()	  
	conv <- 0
	k <- 1	
	gamma <- .25
	startvalues <- pars	
	stagecycle <- 1	
	converge <- 1
	nconstvalues <- sum(constvalues[,1] == 1)		
	if(length(equalconstr) > 0)	
		for(i in 1:length(equalconstr))
			nconstvalues <- nconstvalues + length(equalconstr[[i]]) - 1
	if(debug){
		print(lambdas)
		print(zetas)
		print(guess)
		print(gmeans)
		print(gcov)		
	}		
	
	####Big MHRM loop
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
		zetas <- list()
		ind1 <- 1
		for(i in 1:J){ 
			zetas[[i]] <- pars[zetaind][ind1:(ind1+sum(estzetas[[i]])-1)]			
			ind1 <- ind1 + sum(estzetas[[i]])
		}		
		guess <- pars[guessind]		
		mu <- grouplist$u <- pars[meanind]
		sig <- matrix(0,nfact,nfact)
		sig[lower.tri(sig,diag=TRUE)] <- pars[covind]
		if(nfact > 1)
			sig <- sig + t(sig) - diag(diag(sig))				
		grouplist$sig <- sig			
		
		#Step 1. Generate m_k datasets of theta 
		for(j in 1:4) theta0 <- draw.thetas(theta0,lambdas,pars[zetaind],guess,
			fulldata,K,itemloc,cand.t.var,sig,mu,estComp)	
		for(i in 1:k) m.thetas[[i]] <- draw.thetas(theta0,lambdas,pars[zetaind],guess,fulldata,
			K,itemloc,cand.t.var,sig,mu,estComp)
		theta0 <- m.thetas[[1]]
		
		#Step 2. Find average of simulated data gradients and hessian 		
		g.m <- h.m <- group.m <- list()
		g <- rep(0,npars)
		h <- matrix(0,npars,npars)	
		for (j in 1:k) {
            g <- rep(NA, npars)
            loc <- 1
            for (i in 0:(J - 1)) {
				if(estComp[i+1]){
					if (estGuess[i + 1]) {
						
					} else {
						temp <- dpars.comp(lambdas[i + 1,][estlam[i+1,]], zetas[[i+1]], 
							guess[i+1], fulldata[, itemloc[i + 1]], m.thetas[[j]])
						ind <- parind[is.na(g)][1]	
						if(i > 0){
							g[is.na(g)][1] <- 0
							ind <- ind + 1
						}
						ind2 <- ind + length(zetas[[i+1]])*2 - 1
						g[ind:ind2] <- temp$grad
						h[ind:ind2, ind:ind2] <- temp$hess
						loc <- loc + length(zetas[[i+1]])*2 - 1
					}				
					next
				}
                if (estGuess[i + 1]) {
					temp <- dpars.dich(lambdas[i + 1, ], zetas[loc], 
						guess[i + 1], fulldata[, itemloc[i + 1]], 
						m.thetas[[j]], estGuess[i + 1])
					ind <- parind[is.na(g)][1]
					ind2 <- ind + nfact + estGuess[i + 1]
					g[ind:ind2] <- temp$grad
					h[ind:ind2, ind:ind2] <- temp$hess
					loc <- loc + 1
				} else {					
					temp <- dpars.poly(lambdas[i + 1, ], zetas[[i+1]], 
						fulldata2[, itemloc[i + 1]:(itemloc[i + 2] - 1)], m.thetas[[j]])
					ind <- parind[is.na(g)][1]
					if(i > 0){
						g[is.na(g)][1] <- 0
						ind <- ind + 1
					}						
					ind2 <- ind + nfact + K[i + 1] - 2
					g[ind:ind2] <- temp$grad
					h[ind:ind2, ind:ind2] <- temp$hess
					loc <- loc + K[i + 1] - 1
                }
            }
			g[is.na(g)] <- 0
			tmp <- d.group(grouplist,m.thetas[[j]])
			g[groupind] <- tmp$g
			h[groupind,groupind] <- tmp$h
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
		if(is.na(attr(theta0,"log.lik"))) stop('Estimation halted. Model did not converge.')		
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
			correct[constvalues[,1] == 1] <- 0
			if(length(equalconstr) > 0)	
				for(i in 1:length(equalconstr))
					correct[equalconstr[[i]]] <- mean(correct[equalconstr[[i]]])			
			pars <- pars + gamma*correct
			if(printcycles && (cycles + 1) %% 10 == 0){ 
				cat(", Max Change =", sprintf("%.4f",max(abs(gamma*correction))), "\n")
				flush.console()
			}			
			pars[covind][pars[covind] > .95] <- parsold[covind][pars[covind] > .95]
			pars[covind][pars[covind] < -.95] <- parsold[covind][pars[covind] < -.95]						
			if(stagecycle == 2) SEM.stores[cycles - burnin,] <- pars
			next
		}	 
		
		#Step 3. Update R-M step		
		Tau <- Tau + gamma*(ave.h - Tau)			
		correction <- SparseM::solve(Tau) %*% grad																
		parsold <- pars
		correct <- rep(0,npars)
		correct[sind] <- correction
		correct[constvalues[,1] == 1] <- 0
		if(length(equalconstr) > 0)		
			for(i in 1:length(equalconstr))
				correct[equalconstr[[i]]] <- mean(correct[equalconstr[[i]]])	
		if(printcycles && (cycles + 1) %% 10 == 0){ 
			cat(", gam = ",sprintf("%.3f",gamma),", Max Change = ", 
				sprintf("%.4f",max(abs(gamma*correction))), "\n", sep = '')
			flush.console()		
		}	
		if(all(gamma*correct < tol)) conv <- conv + 1
			else conv <- 0		
		if(conv == 3) break	
		pars <- pars + gamma*correct	
		pars[covind][pars[covind] > .95] <- parsold[covind][pars[covind] > .95]
		pars[covind][pars[covind] < -.95] <- parsold[covind][pars[covind] < -.95]
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE 			
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	} ###END BIG LOOP
	
	cat("\n\n")
	SEtmp <- diag(solve(info))		
	if(any(SEtmp < 0)){
		warning("Information matrix is not positive definite, negative SEs set to 'NA'.\n")
		SEtmp[SEtmp < 0] <- NA
	}	
	SEtmp <- sqrt(SEtmp)	
	SE <- rep(NA,npars) 
	SE[parind[sind]] <- SEtmp
	SE[constvalues[,1]==1] <- NA
	if(length(equalconstr) > 0)
		for(i in 1:length(equalconstr))
			SE[equalconstr[[i]]] <- mean(SE[equalconstr[[i]]])
	estpars <- pars[sind]
	lambdas <- matrix(pars[lamind],J,nfact,byrow=TRUE)	
	lambdas[!estlam & !lambdas != 0] <- NA	
	guess <- rep(NA,J)
	guess[estGuess] <- pars[guessind]
	guess[K == 2 & !estGuess] <- 0
	zetas <- pars[zetaind]
	u <- pars[meanind]	
	sig <- matrix(0,nfact,nfact)
	SElam <- matrix(SE[lamind],J,nfact,byrow=TRUE)
	SEzetas <- SE[zetaind]	
	SEg <- rep(NA,J)	
	SEg[estGuess] <- SE[guessind]	
	SEu <- SE[meanind]	
	SEsig <- matrix(0,nfact,nfact)	
	tmp <- pars[covind]
	tmp2 <- SE[covind]
	loc <- 1
	for(i in 1:nfact){
		for(j in 1:nfact){
			if(i <= j){
				sig[i,j] <- tmp[loc]
				SEsig[i,j] <- tmp2[loc]
				loc <- loc + 1
			}
		}
	}
	if(nfact > 1) {	
		sig <- sig + t(sig) - diag(diag(sig))
		SEsig <- SEsig + t(SEsig) - diag(diag(SEsig))	
	} else SEsig <- NA
	if(any(estComp)){
		if((max(K)-1) > nfact) tmp1 <- tmp2 <- matrix(NA,J,(max(K)-1))
		else tmp1 <- tmp2 <- matrix(NA,J,nfact)
	} else tmp1 <- tmp2 <- matrix(NA,J,(max(K)-1))
	
	loc <- 1
	for(i in 1:J){
		if(!estComp[i]){
			for(j in 1:(K[i]-1)){
				tmp1[i,j] <- zetas[loc] 
				tmp2[i,j] <- SEzetas[loc]
				loc <- loc + 1
			}
		} else {
			for(j in 1:nfact){
				tmp1[i,j] <- zetas[loc]
				tmp2[i,j] <- SEzetas[loc]
				loc <- loc + 1
			}	
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
		
	if (nfact > 1) norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2,na.rm = TRUE))
		else norm <- as.matrix(sqrt(1 + pars[ ,1]^2))  
	F <- as.matrix(pars[ ,1:nfact]/norm)
	F[is.na(F)] <- 0		
	h2 <- rowSums(F^2)

	mod <- new('confmirtClass', pars=pars, guess=guess, SEpars=SEpars, SEg = SEg, 
		gpars=gpars, SEgpars=SEgpars, estpars=estpars,cycles=cycles - SEM.cycles 
		- burnin, Theta=theta0, fulldata=fulldata, data=data, K=K, itemloc=itemloc, 
		h2=h2,F=F,converge = converge, nconstvalues = as.integer(nconstvalues), 
		estComp=estComp, Call=Call)
	if(calcLL){
		cat("Calculating log-likelihood...\n")
		flush.console()
		mod <- logLik(mod,draws)		
	}	
	return(mod)
}	

