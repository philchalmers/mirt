# Class "polymirtClass"
# 
# Defines the object returned from \code{\link{polymirt}}.
# 
# 
# @name polymirtClass-class
# @aliases polymirtClass-class coef,polymirtClass-method
# plot,polymirtClass,missing-method print,polymirtClass-method
# residuals,polymirtClass-method show,polymirtClass-method
# summary,polymirtClass-method anova,polymirtClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("polymirtClass", ...)}.
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass polymirtClass
# @keywords classes
setClass(
	Class = 'polymirtClass',
	representation = representation(pars = 'matrix', guess = 'numeric', SEpars = 'matrix', 
		cycles = 'numeric', Theta = 'matrix', fulldata = 'matrix', data = 'matrix', 
		K = 'numeric', F = 'matrix', h2 = 'numeric', itemloc = 'numeric', AIC = 'numeric',
		converge = 'numeric', logLik = 'numeric', SElogLik = 'numeric', df = 'integer', 
		G2 = 'numeric', p = 'numeric', tabdata = 'matrix', BIC = 'numeric', estGuess = 'logical', 
		Call = 'call'),	
	validity = function(object) return(TRUE)
)	

#' Full-Information Item Factor Analysis for Mixed Data Formats
#' 
#' \code{polymirt} fits an unconditional (exploratory) full-information
#' maximum-likelihood factor analysis model to dichotomous and polychotomous
#' data under the item response theory paradigm using Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm.
#'
#' 
#' \code{polymirt} follows the item factor analysis strategy by a stochastic
#' version of maximum likelihood estimation described by Cai (2010). The
#' general equation used for multidimensional item response theory in this
#' package is in the logistic form with a scaling correction of 1.702. This
#' correction is applied to allow comparison to mainstream programs such as
#' TESTFACT (2003) and POLYFACT. Missing data are treated as 'missing at
#' random' so that each response vector is included in the estimation (i.e.,
#' full-information). Residuals are computed using the LD statistic (Chen &
#' Thissen, 1997) in the lower diagonal of the matrix returned by
#' \code{residuals}, and Cramer's V above the diagonal. For computing the
#' log-likelihood more accurately see \code{\link{logLik}}.
#' 
#' Use of \code{plot} will display the test information function for 1 and 2
#' dimensional solutions. To examine individuals item plots use
#' \code{\link{itemplot}} (although the \code{\link[plink]{plink}} package is
#' much more suitable for IRT graphics) which will also plot information and
#' surface functions.
#' 
#' \code{coef} displays the item parameters with their associated standard
#' errors, while use of \code{summary} transforms the slopes into a factor
#' loadings metric. Also, factor loading values below a specified constant can
#' be also be suppressed in \code{summary} to allow better visual clarity.
#' Models may be compared by using the \code{anova} function, where a
#' Chi-squared difference test and AIC difference values are displayed.
#' 
#' @aliases polymirt summary,polymirt-method coef,polymirt-method
#' plot,polymirt-method residuals,polymirt-method anova,polymirt-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data
#' @param nfact number of factors to be extracted
#' @param guess fixed values for the pseudo-guessing parameter. Can be entered
#' as a single value to assign a global guessing parameter or may be entered as
#' a numeric vector for each item
#' @param estGuess a logical vector indicating which lower-asymptote parameters
#' to be estimated (default is null, and therefore is contigent on the values
#' in \code{guess}). By default, if any value in \code{guess} is greater than 0
#' then its respective \code{estGuess} value is set to \code{TRUE}.
#' Additionally, beta priors are automatically imposed for estimated parameters
#' that correspond to the input guessing value.
#' @param prev.cor use a previously computed correlation matrix to be used to
#' estimate starting values the estimation. The input could be any correlation
#' matrix, but it is advised to use a matrix of polychoric correlations.
#' @param rotate type of rotation to perform after the initial orthogonal
#' parameters have been extracted. See \code{\link{mirt}} for a list of
#' possible rotations
#' @param ncycles the maximum number of iterations to be performed
#' @param burnin number of burn-in cycles to perform before beginning the SEM
#' stage
#' @param SEM.cycles number of stochastic EM cycles to perform before beginning
#' the MH-RM algorithm
#' @param kdraws number of Metropolis-Hastings imputations of the factor scores
#' at each iteration. Default is 1
#' @param tol tolerance that will terminate the model estimation; must occur in
#' 3 consecutive iterations
#' @param SE logical; display the standard errors?
#' @param x an object of class \code{polymirt} to be plotted or printed
#' @param object a model estimated from \code{polymirt} of class
#' \code{polymirt}
#' @param object2 a model estimated from \code{polymirt} of class
#' \code{polymirt}
#' @param suppress a numeric value indicating which (possibly rotated) factor
#' loadings should be suppressed. Typical values are around .3 in most
#' statistical software
#' @param digits the number of significant digits to be rounded
#' @param npts number of quadrature points to be used for plotting features.
#' Larger values make plots look smoother
#' @param rot allows rotation of the 3D graphics
#' @param printcycles logical; display iteration history during estimation?
#' @param calcLL logical; calculate the log-likelihood?
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param draws the number of Monte Carlo draws to estimate the log-likelihood
#' @param type either \code{'info'} or \code{'infocontour'} to plot test
#' information plots
#' @param debug logical; turn on debugging features?
#' @param technical list specifying subtle parameters that can be adjusted
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{polymirt}},
#' \code{\link{itemplot}}
#' @references
#' 
#' Cai, L. (2010). High-Dimensional exploratory item factor analysis by a
#' Metropolis-Hastings Robbins-Monro algorithm. \emph{Psychometrika, 75},
#' 33-57.
#' 
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., &
#' Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics,
#' and Full-information Item Factor Analysis [Computer software]. Lincolnwood,
#' IL: Scientific Software International.
#' @keywords models
#' @usage 
#' polymirt(data, nfact, guess = 0, estGuess = NULL, prev.cor = NULL, ncycles = 2000, 
#'   burnin = 100, SEM.cycles = 50, kdraws = 1, tol = .001, printcycles = TRUE, calcLL = TRUE, 
#'   draws = 2000, debug = FALSE, technical = list(), ...)
#' 
#' \S4method{summary}{polymirt}(object, rotate='varimax', suppress = 0, digits = 3, ...)
#' 
#' \S4method{coef}{polymirt}(object, SE = TRUE, digits = 3, ...)
#' 
#' \S4method{plot}{polymirt}(x, npts = 50, type = 'info', rot = list(x = -70, y = 30, z = 10), ...)
#' 
#' \S4method{residuals}{polymirt}(object, restype = 'LD', digits = 3, ...)
#' 
#' \S4method{anova}{polymirt}(object, object2, ...)
#'
#' @export polymirt
#' @examples
#' 
#' \dontrun{
#' #load LSAT section 7 data and compute 1 and 2 factor models
#' data(LSAT7)
#' fulldata <- expand.table(LSAT7)
#' 
#' (mod1 <- polymirt(fulldata, 1))
#' summary(mod1)
#' residuals(mod1)
#' 
#' (mod2 <- polymirt(fulldata, 2))
#' summary(mod2)
#' coef(mod2)
#' anova(mod1,mod2)
#' 
#' ###########
#' #data from the 'ltm' package in numeric format
#' data(Science)
#' (mod1 <- polymirt(Science, 1))
#' summary(mod1)
#' residuals(mod1)
#' coef(mod1)
#' 
#' (mod2 <- polymirt(Science, 2, calcLL = FALSE)) #don't calculate log-likelihood
#' mod2 <- logLik(mod2,5000) #calc log-likelihood here with more draws
#' summary(mod2, 'promax', suppress = .3)
#' coef(mod2)
#' anova(mod1,mod2)
#' 
#' 
#'      }
#' 
polymirt <- function(data, nfact, guess = 0, estGuess = NULL, prev.cor = NULL, ncycles = 2000, 
	burnin = 100, SEM.cycles = 50, kdraws = 1, tol = .001, printcycles = TRUE,
	calcLL = TRUE, draws = 2000, debug = FALSE, technical = list(), ...)
{		
	Call <- match.call()
	set.seed(12345)
	if(!is.null(technical$set.seed)) set.seed(technical$set.seed)
	ifelse(!is.null(technical$guess.prior.n), guess.prior.n <- technical$guess.prior.n,
		guess.prior.n <- 20)
	itemnames <- colnames(data)
	data <- as.matrix(data)		
	J <- ncol(data)
	N <- nrow(data)	
	if(length(guess) == 1) guess <- rep(guess,J)
	colnames(data) <- itemnames
	if(length(guess) > J || length(guess) < J) 
		stop("The number of guessing parameters is incorrect.")
	if(is.null(estGuess))
		estGuess <- guess > 0					
	uniques <- list()
	for(i in 1:J)
		uniques[[i]] <- sort(unique(data[,i]))
	K <- rep(0,J)
	for(i in 1:J) K[i] <- length(uniques[[i]])	
	guess[K > 2] <- 0
	estGuess[K > 2] <- FALSE	
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
	if(!is.null(prev.cor)){
		if (ncol(prev.cor) == nrow(prev.cor)) Rpoly <- prev.cor
			else stop("Correlation matrix is not square.\n")
	} 	else Rpoly <- cormod(na.omit(data),K,guess)
	FA <- fa(Rpoly,nfact,rotate = 'none', warnings= FALSE, fm="minres")	
	loads <- unclass(loadings(FA))
	u <- FA$unique
	u[u < .001 ] <- .2
	cs <- sqrt(u)
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
	npars <- length(c(lambdas,zetas)) + sum(estGuess) 
	parind <- 1:npars
	pars <- rep(NA,npars)
	Ksum <- cumsum(K-1 + nfact + estGuess)
	Ksum <- Ksum - min(Ksum) + K[1]
	lamind	<- gind <- c()	 
	for(i in 1:J){
		pars[Ksum[i]:(Ksum[i] + nfact - 1)] <- lambdas[i,]
		lamind <- c(lamind,Ksum[i]:(Ksum[i] + nfact - 1))
		if(estGuess[i]){
			pars[Ksum[i] + nfact] <- guess[i]
			gind <- c(gind,Ksum[i] + nfact)
		}	
	}	
	zetaind <- parind[is.na(pars)]			
	pars[is.na(pars)] <- zetas
	diag(Rpoly) <- 1	
	converge <- 1
	guessPrior <- list()
	guessPriorCount <- 1
	if(sum(estGuess) > 0){
		for(i in 1:J){
			if(estGuess[i]){
				guessPrior[[guessPriorCount]] <- c(gind[i],guess[i]*guess.prior.n,
					(1-guess[i])*guess.prior.n)
				guessPriorCount <- guessPriorCount + 1			
			}
		}	
	}
	if(debug){
		print(lambdas)
		print(zetas)
	}	
	
    #preamble for MRHM algorithm		
	theta0 <- matrix(0,N,nfact)	
	cand.t.var <- 1	
	tmp <- .1
	for(i in 1:30){			
		theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)
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
	m.thetas <- list()		
	SEM.stores <- matrix(0,SEM.cycles,npars)
	phi <- rep(0,npars)
	Tau <- info <- h <- matrix(0,npars,npars)
	m.list <- list()	  
	conv <- noninvcount <- 0
	k <- 1	
	gamma <- 0.25
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
		lambdas <- matrix(pars[lamind],ncol=nfact,byrow=TRUE)
		zetas <- pars[zetaind]
		guess <- rep(0,J)
		guess[estGuess] <- pars[gind]		
		
		#Step 1. Generate m_k datasets of theta 
		for(j in 1:4) theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)
		for(i in 1:k)			
			m.thetas[[i]] <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)		
		theta0 <- m.thetas[[1]]
		
		#Step 2. Find average of simulated data gradients and hessian 
		g.m <- h.m <- list()					
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
		if(length(guessPrior) > 0){
			for(i in 1:length(guessPrior)){
				tmp <- guessPrior[[i]]				
				tmp2 <- betaprior(tmp[2],tmp[3],pars[tmp[1]])				
				grad[tmp[1]] <- grad[tmp[1]] + tmp2$g
				ave.h[tmp[1],tmp[1]] <- ave.h[tmp[1],tmp[1]] + tmp2$h
			}		
		}
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
			ave.h <- as(ave.h,'sparseMatrix')
			inv.ave.h <- try(solve(ave.h))		    
			if(class(inv.ave.h) == 'try-error'){
				inv.ave.h <- try(solve(ave.h + 2*diag(ncol(ave.h))))
				noninvcount <- noninvcount + 1
				if(noninvcount == 3) 
					stop('\nEstimation halted during burn in stages, solution is unstable')
			}
			correction <- inv.ave.h %*% grad
			parsold <- pars
			correction[correction > .5] <- .5
			correction[correction < -0.5] <- -0.5				
			pars <- pars + gamma*as.vector(correction)
			if(printcycles && (cycles + 1) %% 10 == 0){ 
				cat(", Max Change =", sprintf("%.4f",max(abs(gamma*correction))), "\n")
				flush.console()			
			}	
			pars[gind][pars[gind] < 0] <- parsold[gind][pars[gind] < 0]			
			if(stagecycle == 2) SEM.stores[cycles - burnin,] <- pars
			next
		}	
		
		#Step 3. Update R-M step		
		Tau <- Tau + gamma*(ave.h - Tau)
		Tau <- as(Tau,'sparseMatrix')	
		inv.Tau <- try(solve(Tau))		
		if(class(inv.Tau) == 'try-error'){
			inv.Tau <- try(solve(Tau + 2 * diag(ncol(Tau))))
			noninvcount <- noninvcount + 1
			if(noninvcount == 3) 
				stop('\nEstimation halted during burn stage 3, solution is unstable')
		}
		correction <- inv.Tau %*% grad	
		correction[correction > .5] <- .5
		correction[correction < -0.5] <- -0.5										
		if(printcycles && (cycles + 1) %% 10 == 0){ 
			cat(", gam = ",sprintf("%.3f",gamma),", Max Change = ", 
				sprintf("%.4f",max(abs(gamma*correction))), "\n", sep='')
			flush.console()			
		}	
		if(all(abs(parsold - pars) < tol)) conv <- conv + 1
			else conv <- 0	
		if(conv == 3) break		
		parsold <- pars
		pars <- pars + gamma*as.vector(correction)
		pars[gind][pars[gind] < 0] <- parsold[gind][pars[gind] < 0]	
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE		
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	}
	cat("\n\n")	
	SE <- diag(solve(info))
	if(any(SE < 0)){
		warning("Information matrix is not positive definite.\n")
		SE <- rep(0,npars)
	}
	if(any(guess < 0)) warning("Negative lower asymptote parameter(s). \n")					
	SE <- sqrt(SE)	
	lambdas <- matrix(pars[lamind],ncol=nfact,byrow=TRUE)
	SElam <- matrix(SE[lamind],ncol=nfact,byrow=TRUE)
	SEg <- guess <- rep(NA,J)
	guess[estGuess] <- pars[gind]
	SEg[estGuess] <- SE[gind]
	zetas <- SEzeta <- matrix(NA,J,(max(K)-1))
	temp <- pars[zetaind]
	temp1 <- SE[zetaind]	
	k <- 1
	for(i in 1:J){
		for(j in 1:(K[i]-1)){
			zetas[i,j] <- temp[k] 
			SEzeta[i,j] <- temp1[k]
			k <- k + 1
		}
	}	 
	guess[K == 2 & !estGuess] <- 0
	pars <- cbind(lambdas,zetas)
	SEpars <- cbind(SElam,SEzeta,SEg)
	
	if (nfact > 1) norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2))
		else norm <- as.matrix(sqrt(1 + pars[ ,1]^2))  
	alp <- as.matrix(pars[ ,1:nfact]/norm)
	FF <- alp %*% t(alp)
	V <- eigen(FF)$vector[ ,1:nfact]
	L <- eigen(FF)$values[1:nfact]
	if (nfact == 1) F <- as.matrix(V * sqrt(L))
		else F <- V %*% sqrt(diag(L))  
	if (sum(F[ ,1] < 0)) F <- (-1) * F  
	h2 <- rowSums(F^2) 	
		
	mod <- new('polymirtClass',pars=pars, guess=guess, SEpars=SEpars, 
		cycles=cycles-SEM.cycles-burnin, Theta=theta0, fulldata=fulldata, 
		data=data, K=K, F=F, h2=h2, itemloc=itemloc, converge = converge,
		estGuess=estGuess, Call=Call)
	if(calcLL){
		cat("Calculating log-likelihood...\n")
		flush.console()
		mod <- logLik(mod,draws,...)		
	}	
	return(mod)	
}

# Methods

setMethod(
	f = "print",
	signature = signature(x = 'polymirtClass'),
	definition = function(x, ...)
	{
		cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information factor analysis with ", ncol(x@F), " factor",
			if(ncol(x@F)>1) "s", "\n", sep="")
		if(x@converge == 1)	
			cat("Converged in ", x@cycles, " iterations.\n", sep="")
		else 	
			cat("Estimation stopped after ", x@cycles, " iterations.\n", sep="")		
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
	} 
)

setMethod(
	f = "show",
	signature = signature(object = 'polymirtClass'),
	definition = function(object)
	{
		cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information factor analysis with ", ncol(object@F), " factor",
			if(ncol(object@F)>1) "s", "\n", sep="")
		if(object@converge == 1)	
			cat("Converged in ", object@cycles, " iterations.\n", sep="")
		else 	
			cat("Estimation stopped after ", object@cycles, " iterations.\n", sep="")	
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
	} 
)

setMethod(
	f = "summary",
	signature = 'polymirtClass',
	definition = function(object, rotate = 'varimax', suppress = 0, digits = 3, ...)
	{
		nfact <- ncol(object@F)
		if (rotate == 'none' || nfact == 1) {
			F <- object@F
			F[abs(F) < suppress] <- NA
			h2 <- as.matrix(object@h2)    	
			SS <- apply(F^2,2,sum)
			colnames(h2) <- "h2"	
			colnames(F) <- names(SS) <- paste("F_", 1:ncol(F),sep="")
			cat("\nUnrotated factor loadings: \n\n")
			loads <- round(cbind(F,h2),digits)
			print(loads)	    	 
			cat("\nSS loadings: ",round(SS,digits), "\n")
			cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
			invisible(list(F,h2))
		} else {	
			F <- object@F
			h2 <- as.matrix(object@h2)		
			colnames(F) <- paste("F_", 1:ncol(F),sep="")
			colnames(h2) <- "h2"		
			cat("\nRotation: ", rotate, "\n")
			rotF <- Rotate(F,rotate)
			SS <- apply(rotF$loadings^2,2,sum)
			L <- rotF$loadings
			L[abs(L) < suppress] <- NA	
			loads <- round(cbind(L,h2),digits)		
			cat("\nRotated factor loadings: \n\n")
			print(loads,digits)		
			if(attr(rotF, "oblique")){
				cat("\nFactor correlations: \n\n")
				Phi <- rotF$Phi	  
				Phi <- round(Phi, digits)
				colnames(Phi) <- rownames(Phi) <- colnames(F)
				print(Phi)			    
			}		
			cat("\nSS loadings: ",round(SS,digits), "\n")		
			if(!attr(rotF, "oblique")) 
				cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
			if(any(h2 > 1)) 
				warning("Solution has heywood cases. Interpret with caution.") 
			invisible(list(rotF$loadings,h2))  
		}  
	}
)

setMethod(
	f = "coef",
	signature = 'polymirtClass',
	definition = function(object, SE = TRUE, digits = 3, ...)
	{  
		nfact <- ncol(object@Theta)	
		a <- matrix(object@pars[ ,1:nfact],ncol=nfact)
		d <- matrix(object@pars[,(nfact+1):ncol(object@pars)],
			ncol = ncol(object@pars)-nfact)    
		A <- sqrt(apply(a^2,1,sum))
		B <- -d/A  
		if (nfact > 1){  
			parameters <- cbind(object@pars,object@guess,A,B)
			SEs <- object@SEpars	
			colnames(parameters) <- c(paste("a_",1:nfact,sep=""),
				paste("d_",1:(ncol(object@pars)-nfact),sep=""),"guess","mvdisc",
				paste("mvint_",1:(ncol(object@pars)-nfact),sep=""))	
			colnames(SEs) <- c(paste("a_",1:nfact,sep=""),
				paste("d_",1:(ncol(object@pars)-nfact),sep=""),"guess")		
			cat("\nUnrotated parameters, multivariate discrimination and intercept: \n\n")
			print(round(parameters, digits))
			if(SE){
				cat("\nStd. Errors: \n\n")	
				print(round(SEs, digits))
			}				
		} else {
			parameters <- cbind(object@pars,object@guess)
			SEs <- object@SEpars	
			colnames(parameters) <- colnames(SEs) <- c(paste("a_",1:nfact,sep=""),
				paste("d_",1:(ncol(object@pars)-nfact),sep=""),"guess")			
			cat("\nParameter slopes and intercepts: \n\n")	
			print(round(parameters, digits))
			if(SE){
				cat("\nStd. Errors: \n\n")	
				print(round(SEs, digits))
			}
		}
		invisible(parameters)
	}
)

setMethod(
	f = "plot",
	signature = signature(x = 'polymirtClass', y = "missing"),
	definition = function(x, y, type = 'info', npts = 50, 
		rot = list(xaxis = -70, yaxis = 30, zaxis = 10))
	{  		
		if (!type %in% c('info','infocontour')) stop(type, " is not a valid plot type.")
		rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
		K <- x@K		
		nfact <- ncol(x@Theta)
		if(nfact >2) stop("Can't plot high dimensional solutions.")
		a <- as.matrix(x@pars[ ,1:nfact])
		d <- as.matrix(x@pars[ ,(nfact+1):ncol(x@pars)])	
		guess <- x@guess
		guess[is.na(guess)] <- 0
		A <- as.matrix(sqrt(apply(a^2,1,sum)))	
		theta <- seq(-4,4,length.out=npts)
		Theta <- thetaComb(theta, nfact)
		info <- rep(0,nrow(Theta))
		for(j in 1:length(K)){
			if(K[j] > 2){
				P <- P.poly(a[j,], d[j,],Theta, itemexp = FALSE)		
				for(i in 1:K[j]){
					w1 <- P[,i]*(1-P[,i])*A[j]
					w2 <- P[,i+1]*(1-P[,i+1])*A[j]
					I <- ((w1 - w2)^2) / (P[,i] - P[,i+1]) * P[,i]
					info <- info + I
				}
			} else {
				P <- P.mirt(a[j,], d[j,],Theta, guess[j])
				Pstar <- P.mirt(a[j,], d[j,],Theta, 0)
				info <- info + A[j]^2 * P * (1-P) * Pstar/P
			}			
		}		
		plt <- data.frame(cbind(info,Theta))
		if(nfact == 2){						
			colnames(plt) <- c("info", "Theta1", "Theta2")			
			if(type == 'infocontour')												
				contour(theta, theta, matrix(info,length(theta),length(theta)), 
					main = paste("Test Information Contour"), xlab = "Theta 1", ylab = "Theta 2")
			if(type == 'info')
				return(wireframe(info ~ Theta1 + Theta2, data = plt, main = "Test Information", 
					zlab = "I", xlab = "Theta 1", ylab = "Theta 2", scales = list(arrows = FALSE),
					screen = rot))
		} else {
			if(type == 'info')
				plot(Theta, info, type='l',main = 'Test Information', xlab = 'Theta', ylab='Information')
			if(type == 'infocontour') 
				cat('No \'contour\' plots for 1-dimensional models\n')
		}		
	}	  
)	

setMethod(
	f = "residuals",
	signature = signature(object = 'polymirtClass'),
	definition = function(object, restype = 'LD', digits = 3, ...)
	{ 	
		fulldata <- object@fulldata	
		data <- object@data
		data[data==99] <- NA
		N <- nrow(fulldata)
		K <- object@K
		J <- length(K)
		nfact <- ncol(object@F)
		theta <- seq(-4,4, length.out = round(20/nfact))
		Theta <- thetaComb(theta,nfact)
		lambdas <- matrix(object@pars[,1:nfact], J)
		zetas <- as.vector(t(object@pars[,(nfact+1):ncol(object@pars)]))
		zetas <- na.omit(zetas)
		guess <- object@guess
		guess[is.na(guess)] <- 0	
		Ksums <- cumsum(K) - 1	
		itemloc <- object@itemloc
		res <- matrix(0,J,J)
		diag(res) <- NA
		colnames(res) <- rownames(res) <- colnames(data)
		prior <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
		prior <- prior/sum(prior)
		loc <- loc2 <- 1
		if(restype == 'LD'){	
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
						res[j,i] <- sum(((tab - Etab)^2)/Etab) /
							((K[i] - 1) * (K[j] - 1)) * sign(s)
						res[i,j] <- sqrt( abs(res[j,i]) / (N - min(c(K[i],K[j]) - 1)))	
					}
				loc2 <- loc2 + K[j] - 1 	
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
	f = "anova",
	signature = signature(object = 'polymirtClass'),
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
			" (SE = ", se,"), df = ", df, ", p = ", round(1 - pchisq(X2,df),4), "\n", sep="")
		cat("AIC difference = ", round(AICdiff,3)," (SE = ", se,")\n", sep='')  
		cat("BIC difference = ", round(BICdiff,3)," (SE = ", se,")\n", sep='') 
	}		
) 

# @rdname itemplot-methods  
setMethod(
	f = "itemplot",
	signature = signature(object = 'polymirtClass', item = 'numeric'),
	definition = function(object, item, type = 'info', npts = 50,
		rot = list(), ...)
	{		
		if (!type %in% c('info','infocontour')) stop(type, " is not a valid plot type.")
		if(object@K[item] > 2){
			K <- object@K		
			nfact <- ncol(object@Theta)
			a <- as.matrix(object@pars[ ,1:nfact])
			d <- as.matrix(object@pars[ ,(nfact+1):ncol(object@pars)])			
			A <- as.matrix(sqrt(apply(a^2,1,sum)))[item,]
			nzeta <- K[item] - 1
			theta <- seq(-4,4,length.out=npts)
			Theta <- thetaComb(theta, nfact)		
			P <- P.poly(a[item,], d[item,], Theta, itemexp = FALSE)
			info <- rep(0,nrow(P))
			for(i in 1:K[item]){
				w1 <- P[,i]*(1-P[,i])*A
				w2 <- P[,i+1]*(1-P[,i+1])*A
				I <- ((w1 - w2)^2) / (P[,i] - P[,i+1]) * P[,i]
				info <- info + I
			}	
			plt <- data.frame(cbind(info,Theta))		
			if(nfact == 1)	
				plot(Theta, info, type='l',main = paste('Item', item,'Information'), 
					xlab = 'Theta', ylab='Information')
			else {					
				colnames(plt) <- c('info','Theta1','Theta2')
				if(type == 'info')
					return(wireframe(info ~ Theta1 + Theta2, data = plt, main = paste("Item",item,"Information"), 
						zlab = "I", xlab = "Theta 1", ylab = "Theta 2", scales = list(arrows = FALSE),
						screen = rot))				
				if(type == 'infocontour'){										
					contour(theta, theta, matrix(info,length(theta),length(theta)), 
						main = paste("Item", item,"Information Contour"), xlab = "Theta 1", ylab = "Theta 2")					
				}
			}	
		} else {
			class(object) <- 'mirtClass'
			itemplot(object,item,type,npts,rot)		 
		}	
	}
)

#' @rdname fscores-methods  
setMethod(
	f = "fscores",
	signature = 'polymirtClass',
	definition = function(object, full.scores = FALSE, ndraws = 3000, thin = 5, ...)
	{ 	
		cand.t.var <- 1
		theta0 <- object@Theta
		K <- object@K
		nfact <- ncol(theta0)
		lambdas <- matrix(object@pars[,1:nfact],ncol=nfact)
		zetas <- na.omit(as.numeric(t(object@pars[,(nfact+1):ncol(object@pars)])))
		guess <- object@guess
		guess[is.na(guess)] <- 0
		data <- cbind(object@data,object@fulldata)
		Names <- c(colnames(object@data[,1:length(K)]),paste("F",1:nfact,sep=''),paste("SE_F",1:nfact,sep=''))
		tabdata <- unique(data)[,-c(1:length(K))]			
		itemloc <- object@itemloc
		Theta <- list()
		for(i in 1:nfact)
			Theta[[i]] <- matrix(0,ncol=ndraws/thin,nrow=nrow(tabdata))		
		theta0 <- matrix(0,nrow(tabdata),nfact)
		for(i in 1:30){			
			theta0 <- draw.thetas(theta0,lambdas,zetas,guess,tabdata,K,itemloc,cand.t.var)
			if(attr(theta0,'Proportion Accepted') > .4) cand.t.var <- cand.t.var + .2
			if(attr(theta0,'Proportion Accepted') < .3) cand.t.var <- cand.t.var - .2
		}
		ind <- 1
		for(i in 1:ndraws){			
			theta0 <- draw.thetas(theta0,lambdas,zetas,guess,tabdata,K,itemloc,cand.t.var)
			if(i %% thin == 0){
				for(j in 1:nfact)
					Theta[[j]][,ind] <- theta0[,j]									
				ind <- ind + 1
			}			
		}

		expscores <- matrix(0,ncol=nfact,nrow=nrow(tabdata))
		sdscores <- matrix(0,ncol=nfact,nrow=nrow(tabdata))
		for(i in 1:nfact){
			expscores[,i] <- rowMeans(Theta[[i]])
			sdscores[,i] <- apply(Theta[[i]],1,sd)
		}
				
		ret <- cbind(unique(data)[,1:length(K)],expscores,sdscores)
		colnames(ret) <- Names
		
		if(!full.scores){ 
			ret <- ret[order(expscores[,1]),]
			rownames(ret) <- NULL
			return(ret)
		} else {
			fulldata <- object@data
			scoremat <- matrix(0,nrow=nrow(fulldata),ncol=nfact)
			colnames(scoremat) <- paste("F",1:nfact,sep='')
			tmp <- unique(data)[,1:length(K)]
			for (j in 1:nrow(tabdata)){          
				TFvec <- colSums(ifelse(t(fulldata) == tmp[j, ],1,0)) == ncol(fulldata)        
				scoremat[TFvec, ] <- expscores[j, ]
			}              
			return(cbind(object@data,scoremat))
		}	
	}	
)
