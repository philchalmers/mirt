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
	representation = representation(pars = 'list', guess = 'numeric', 
		SEpars = 'matrix', cycles = 'numeric', Theta = 'matrix', fulldata = 'matrix', data = 'matrix', 
		K = 'numeric', F = 'matrix', h2 = 'numeric', itemloc = 'numeric', AIC = 'numeric',
		converge = 'numeric', logLik = 'numeric', SElogLik = 'numeric', df = 'integer', 
		G2 = 'numeric', p = 'numeric', tabdata = 'matrix', BIC = 'numeric', estGuess = 'logical', 
		RMSEA = 'numeric', rotate='character', Call = 'call'),	
	validity = function(object) return(TRUE)
)	

#' Full-Information Item Factor Analysis for Mixed Data Formats
#' 
#' \code{polymirt} fits an unconditional (exploratory) full-information
#' maximum-likelihood factor analysis model to dichotomous and polychotomous
#' data under the item response theory paradigm using Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm. If requested, lower asymptote
#' parameters are estimated with a beta prior included automatically.
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
#' fitted,polymirt-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data
#' @param nfact number of factors to be extracted
#' @param guess starting (or fixed) values for the pseudo-guessing parameter. Can be 
#' entered as a single value to assign a global guessing parameter or may be entered as
#' a numeric vector for each item
#' @param estGuess a logical vector indicating which lower-asymptote parameters
#' to be estimated (default is null, and therefore is contingent on the values
#' in \code{guess}). By default, if any value in \code{guess} is greater than 0
#' then its respective \code{estGuess} value is set to \code{TRUE}.
#' Additionally, beta priors are automatically imposed for estimated parameters
#' that correspond to the input guessing value.
#' @param prev.cor use a previously computed correlation matrix to be used to
#' estimate starting values the estimation. The input could be any correlation
#' matrix, but it is advised to use a matrix of polychoric correlations.
#' @param rotate type of rotation to perform after the initial orthogonal
#' parameters have been extracted by using \code{summary}; default is \code{'varimax'}. 
#' See \code{\link{mirt}} for list of possible rotations. If \code{rotate != ''} in the 
#' \code{summary} input then the default from the object is ignored and the new rotation 
#' from the list is used instead
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
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#' option. Only prints patterns that have standardized residuals greater than 
#' \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param print logical; print output to console?
#' @param x an object of class \code{polymirtClass} to be plotted or printed
#' @param object a model estimated from \code{polymirtClass} of class
#' \code{polymirt}
#' @param object2 a model estimated from \code{polymirt} of class
#' \code{polymirtClass}
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
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6),
#' 1-29.
#' 
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., &
#' Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics,
#' and Full-information Item Factor Analysis [Computer software]. Lincolnwood,
#' IL: Scientific Software International.
#' @keywords models
#' @usage 
#' polymirt(data, nfact, guess = 0, estGuess = NULL, prev.cor = NULL, rotate = 'varimax', 
#'    ncycles = 2000, burnin = 100, SEM.cycles = 50, kdraws = 1, tol = .001, 
#'    printcycles = TRUE,	calcLL = TRUE, draws = 2000, debug = FALSE, technical = list(), ...)
#'  
#' 
#' 
#' \S4method{summary}{polymirt}(object, rotate='', suppress = 0, digits = 3, print = FALSE, ...)
#' 
#' \S4method{coef}{polymirt}(object, rotate = '', SE = TRUE, digits = 3, ...)
#' 
#' \S4method{plot}{polymirt}(x, npts = 50, type = 'info', rot = list(x = -70, y = 30, z = 10), ...)
#' 
#' \S4method{residuals}{polymirt}(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
#' 
#' \S4method{anova}{polymirt}(object, object2, ...)
#'
#' \S4method{fitted}{polymirt}(object, digits = 3, ...)
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
polymirt <- function(data, nfact, guess = 0, estGuess = NULL, prev.cor = NULL, rotate = 'varimax', 
    ncycles = 2000, burnin = 100, SEM.cycles = 50, kdraws = 1, tol = .001, 
    printcycles = TRUE,	calcLL = TRUE, draws = 2000, debug = FALSE, technical = list(), ...)
{		
	Call <- match.call()
	set.seed(12345)
	if(!is.null(technical$set.seed)) set.seed(technical$set.seed)
	guess.prior.n <- ifelse(!is.null(technical$guess.prior.n),  
                            technical$guess.prior.n, 20)
	itemnames <- colnames(data)
	data <- as.matrix(data)	
	if(!any(data %in% c(0:20,NA))) 
		stop("Data must contain only numeric values (including NA).")	
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
	if(!is.null(prev.cor)){
		if (ncol(prev.cor) == nrow(prev.cor)) Rpoly <- prev.cor
			else stop("Correlation matrix is not square.\n")
	} 	else Rpoly <- cormod(na.omit(data),K,guess)
	FA <- psych::fa(Rpoly,nfact,rotate = 'none', warnings= FALSE, fm="minres")	
	loads <- unclass(loadings(FA))
	u <- FA$unique
	u[u < .001 ] <- .2
	cs <- sqrt(u)
	lambdas <- loads/cs	
    zetas <- zetaindlist <- list()
	zetalong <- c()
    for(i in 1:J){
        if(K[i] == 2){
            zetas[[i]] <- qnorm(mean(fulldata[,itemloc[i]]))/cs[i]            
			zetalong <- c(zetalong, zetas[[i]])
        } else {
            temp <- table(data[,i])[1:(K[i]-1)]/N
            temp <- cumsum(temp)			
            zetas[[i]] <- qnorm(1 - temp)/cs[i]        
			zetalong <- c(zetalong, zetas[[i]])
        }       
    }
    nzetas <- 0
    for(i in 1:J) nzetas <- nzetas + length(zetas[[i]])
	npars <- length(lambdas) + nzetas + sum(estGuess) 
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
	tmp <- 1
	for(i in 1:J){
		zetaindlist[[i]] <- zetaind[tmp:(tmp + length(zetas[[i]]) - 1)]
		tmp <- tmp + length(zetas[[i]])
	}
	pars[is.na(pars)] <- zetalong
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
	indlist <- list(lamind=lamind,zetaind=zetaindlist,gind=gind)
	if(debug){
		print(indlist)
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
		
		normpars <- sortPars(pars, indlist, nfact, estGuess)
		lambdas <- normpars$lambdas
		zetas <- normpars$zetas		 
		guess <- normpars$guess		
		
		#Step 1. Generate m_k datasets of theta 
		for(j in 1:4) theta0 <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)
		for(i in 1:k)			
			m.thetas[[i]] <- draw.thetas(theta0,lambdas,zetas,guess,fulldata,K,itemloc,cand.t.var)		
		theta0 <- m.thetas[[1]]
		
		#Step 2. Find average of simulated data gradients and hessian 
		g.m <- h.m <- list()					
		for(j in 1:k){
			g <- rep(NA,npars)			
			for(i in 1:J){
				if(K[i] == 2){
					temp <- dpars.dich(lambdas[i, ], zetas[[i]],guess[i],
						fulldata[ ,itemloc[i]],m.thetas[[j]],estGuess[i])
					ind <- parind[is.na(g)][1]
					ind2 <- ind + length(temp$g) - 1		
					g[ind:ind2] <- temp$grad
					h[ind:ind2,ind:ind2] <- temp$hess					
				} else {						
					temp <- dpars.poly(lambdas[i, ],zetas[[i]],
						fulldata[ ,itemloc[i]:(itemloc[i+1]-1)],m.thetas[[j]])
					ind <- parind[is.na(g)][1]	
					ind2 <- ind + length(temp$g) - 1		
					g[ind:ind2] <- temp$grad
					h[ind:ind2,ind:ind2] <- temp$hess					
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
	SEpars <- sortPars(SE, indlist, nfact, estGuess)
	normpars <- sortPars(pars, indlist, nfact, estGuess)
	lambdas <- normpars$lambdas
	zetas <- normpars$zetas		 
	guess <- normpars$guess		
	SElam <- SEpars$lambdas
	SEzetas <- SEpars$zetas		 
	SEg <- SEpars$guess
	SEg[!estGuess] <- NA	
		
	zetatable <- SEzetatable <- matrix(NA,J,(max(K)-1))		
	for(i in 1:J){
		for(j in 1:(K[i]-1)){
			zetatable[i,j] <- zetas[[i]][j]
			SEzetatable[i,j] <- SEzetas[[i]][j]
			
		}
	}	 
	guess[K == 2 & !estGuess] <- 0
	pars <- cbind(lambdas,zetatable)
	SEpars <- cbind(SElam,SEzetatable,SEg)
	
	if (nfact > 1) norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2))
		else norm <- as.matrix(sqrt(1 + pars[ ,1]^2))  
	alp <- as.matrix(pars[ ,1:nfact]/norm)
	FF <- alp %*% t(alp)
	V <- eigen(FF)$vector[ ,1:nfact]
	L <- eigen(FF)$values[1:nfact]
	if (nfact == 1) F <- as.matrix(V * sqrt(L))
		else F <- V %*% sqrt(diag(L))  
	if (sum(F[ ,1] < 0)) F <- (-1) * F 
	colnames(F) <- paste("F_", 1:ncol(F),sep="")	
	h2 <- rowSums(F^2) 	
	names(h2) <- itemnames
		
	mod <- new('polymirtClass',pars=normpars, guess=guess, SEpars=SEpars, 
		cycles=cycles-SEM.cycles-burnin, Theta=theta0, fulldata=fulldata, 
		data=data, K=K, F=F, h2=h2, itemloc=itemloc, converge = converge,
		estGuess=estGuess, rotate=rotate, Call=Call)
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
			if(!is.nan(x@p))
				cat("G^2 = ", round(x@G2,2), ", df = ", 
					x@df, ", p = ", round(x@p,4), ", RMSEA = ", round(x@RMSEA,3), "\n", sep="")
			else 
				cat("G^2 = ", NA, ", df = ", 
					x@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="")	
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
			if(!is.nan(object@p))
				cat("G^2 = ", round(object@G2,2), ", df = ", 
					object@df, ", p = ", round(object@p,4), ", RMSEA = ", round(object@RMSEA,3), 
                    "\n", sep="")
			else 
				cat("G^2 = ", NA, ", df = ", 
					object@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="")
		}			
	} 
)

setMethod(
	f = "summary",
	signature = 'polymirtClass',
	definition = function(object, rotate = '', suppress = 0, digits = 3, print = TRUE, ...)
	{
		nfact <- ncol(object@F)
		itemnames <- colnames(object@data)
		if (rotate == 'none' || nfact == 1) {
			F <- object@F
			F[abs(F) < suppress] <- NA
			h2 <- as.matrix(object@h2)    	
			SS <- apply(F^2,2,sum)
			colnames(h2) <- "h2"	
			names(SS) <- colnames(F) 
			loads <- round(cbind(F,h2),digits)
			rownames(loads) <- itemnames
			if(print){
			    cat("\nUnrotated factor loadings: \n\n")
			    print(loads)	    	 
			    cat("\nSS loadings: ",round(SS,digits), "\n")
			    cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
            }
			invisible(list(F,h2))
		} else {	
			F <- object@F
			h2 <- as.matrix(object@h2)		
			colnames(h2) <- "h2"
			if(rotate == '') rotate <- object@rotate
			rotF <- Rotate(F,rotate)
			SS <- apply(rotF$loadings^2,2,sum)
			L <- rotF$loadings
			L[abs(L) < suppress] <- NA	
			loads <- round(cbind(L,h2),digits)
			rownames(loads) <- itemnames			
			Phi <- diag(nfact)
			if(attr(rotF, "oblique")){
				Phi <- rotF$Phi	  
				Phi <- round(Phi, digits)
				colnames(Phi) <- rownames(Phi) <- colnames(F)
				if(print){
				    cat("\nFactor correlations: \n\n")
				    print(Phi)
                }
			}		
			if(print){
			    cat("\nRotation: ", rotate, "\n")
			    cat("\nRotated factor loadings: \n\n")
			    print(loads,digits)	
			    cat("\nSS loadings: ",round(SS,digits), "\n")		
			    if(!attr(rotF, "oblique")) 
				    cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
            }
			if(any(h2 > 1)) 
				warning("Solution has heywood cases. Interpret with caution.") 
			invisible(list(rotF=rotF$loadings,h2=h2,fcor=Phi))  
		}  
	}
)

setMethod(
	f = "coef",
	signature = 'polymirtClass',
	definition = function(object, rotate = '', SE = TRUE, digits = 3, ...)
	{  
		K <- object@K
		a <- object@pars$lambdas		
		d <- matrix(NA, nrow(a), max(K-1))
		zetas <- object@pars$zetas
		for(i in 1:length(K)){
			d[i, 1:(K[i] - 1)] <- zetas[[i]]
		}
		A <- sqrt(apply(a^2,1,sum))
		B <- -d/A  
		if (ncol(a) > 1){  
		    rotname <- ifelse(rotate == '', object@rotate, rotate)
		    so <- summary(object, rotate = rotate, print = FALSE)             
		    a <- rotateLambdas(so)
			parameters <- cbind(a,d,object@guess,A,B)    
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),paste("d_",1:max(K-1),sep=""),"guess", 
				"mvdisc",paste("mvint_",1:max(K-1),sep=""))	  
			rownames(parameters) <- colnames(object@data)
		    cat("\nParameters with", rotname, "rotation, multivariate discrimination and intercept: \n\n")
			print(round(parameters, digits))  	
		} else {
			parameters <- cbind(a,d,object@guess)
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),paste("d_",1:max(K-1),sep=""),"guess")   
			rownames(parameters) <- colnames(object@data)
			cat("\nParameter slopes and intercepts: \n\n")	
			print(round(parameters, digits))	  
		}
		ret <- list(parameters)
		if(length(object@SEpars) > 1){
			if(SE){
				cat("\nStd. Errors: \n\n")	
				SEs <- object@SEpars
				colnames(SEs) <- c(paste("a_",1:ncol(a),sep=""),paste("d_",1:max(K-1),sep=""),"guess") 	
				rownames(SEs) <- rownames(parameters)
				print(SEs, digits)
				ret <- list(parameters,SEs)
			}
		}
		invisible(ret)
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
		if(nfact > 2) stop("Can't plot high dimensional solutions.")
		a <- x@pars$lambdas
		d <- x@pars$zetas	
		guess <- x@guess
		guess[is.na(guess)] <- 0
		A <- as.matrix(sqrt(apply(a^2,1,sum)))	
		theta <- seq(-4,4,length.out=npts)
		Theta <- thetaComb(theta, nfact)
		info <- rep(0,nrow(Theta))
		for(j in 1:length(K)){
			if(K[j] > 2){
				P <- P.poly(a[j,], d[[j]],Theta, itemexp = FALSE)		
				for(i in 1:K[j]){
					w1 <- P[,i]*(1-P[,i])*A[j]
					w2 <- P[,i+1]*(1-P[,i+1])*A[j]
					I <- ((w1 - w2)^2) / (P[,i] - P[,i+1]) * P[,i]
					info <- info + I
				}
			} else {
				P <- P.mirt(a[j,], d[[j]],Theta, guess[j])
				Pstar <- P.mirt(a[j,], d[[j]],Theta, 0)
				info <- info + A[j]^2 * P * (1-P) * Pstar/P
			}			
		}		
		plt <- data.frame(cbind(info,Theta))
		if(nfact == 2){						
			colnames(plt) <- c("info", "Theta1", "Theta2")			
			if(type == 'infocontour')												
				return(contourplot(info ~ Theta1 * Theta2, data = plt, 
					main = paste("Test Information Contour"), xlab = expression(theta[1]), 
					ylab = expression(theta[2])))
			if(type == 'info')
				return(wireframe(info ~ Theta1 + Theta2, data = plt, main = "Test Information", 
					zlab = expression(I(theta)), xlab = expression(theta[1]), ylab = expression(theta[2]), 
					scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE))
		} else {
			colnames(plt) <- c("info", "Theta")
			if(type == 'info')
				return(xyplot(info~Theta, plt, type='l',main = 'Test Information', xlab = expression(theta), 
					ylab = expression(I(theta))))
			if(type == 'infocontour') 
				cat('No \'contour\' plots for 1-dimensional models\n')
		}		
	}	  
)	

setMethod(
	f = "residuals",
	signature = signature(object = 'polymirtClass'),
	definition = function(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
	{ 	
		K <- object@K		
		data <- object@data	
		N <- nrow(data)	
		J <- ncol(data)
		nfact <- ncol(object@F)
		lambdas <- object@pars$lambdas
		zetas <- object@pars$zetas
		guess <- object@guess		
		itemloc <- object@itemloc
		theta <- seq(-4,4, length.out = round(20/nfact))
		Theta <- thetaComb(theta,nfact)
		res <- matrix(0,J,J)
		diag(res) <- NA
		colnames(res) <- rownames(res) <- colnames(data)
		prior <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact))
		prior <- prior/sum(prior)	
		if(restype == 'LD'){	
			for(i in 1:J){								
				for(j in 1:J){			
					if(i < j){
						if(K[i] > 2) P1 <- P.poly(lambdas[i,],zetas[[i]],Theta,itemexp=TRUE)
						else { 
							P1 <- P.mirt(lambdas[i,],zetas[[i]], Theta, guess[i])
							P1 <- cbind(1 - P1, P1)
						}	
						if(K[j] > 2) P2 <- P.poly(lambdas[j,],zetas[[j]],Theta,itemexp=TRUE)
						else {
							P2 <- P.mirt(lambdas[j,],zetas[[j]], Theta, guess[j])	
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
				}
			}	
			cat("LD matrix:\n\n")	
			res <- round(res,digits)
			return(res)
		} 
		if(restype == 'exp'){	
			r <- object@tabdata[ ,ncol(object@tabdata)-1]
			rexp <- object@tabdata[ ,ncol(object@tabdata)]
			res <- round((r - rexp) / sqrt(rexp),digits)
			tabdata <- object@tabdata
			freq <- tabdata[ ,ncol(tabdata)]			
			tabdata[tabdata[ ,1:ncol(object@data)] == 99] <- NA
			tabdata[ ,ncol(tabdata)] <- freq
			tabdata <- cbind(tabdata,res)
			colnames(tabdata) <- c(colnames(data),'freq', 'exp', 'std_res')	
			if(!is.null(printvalue)){
				if(!is.numeric(printvalue)) stop('printvalue is not a number.')
				tabdata <- tabdata[abs(tabdata[ ,ncol(tabdata)]) > printvalue, ]
			}			
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
			tmp <- object
			object <- object2
			object2 <- tmp
		}
		X2 <- 2*object2@logLik - 2*object@logLik 
		AICdiff <- object@AIC - object2@AIC
		BICdiff <- object@BIC - object2@BIC
		se <- round(object@SElogLik + object2@SElogLik,3)	
		cat("\nChi-squared difference: \n\nX2 = ", round(X2,3), 
			" (SE = ", se,"), df = ", df, ", p = ", round(1 - pchisq(X2,abs(df)),4),
			"\n", sep="")
		cat("AIC difference = ", round(AICdiff,3)," (SE = ", se,")\n", sep='')  
		cat("BIC difference = ", round(BICdiff,3)," (SE = ", se,")\n", sep='') 
	}		
) 

setMethod(
	f = "fitted",
	signature = signature(object = 'polymirtClass'),
	definition = function(object, digits = 3, ...){  		  
		tabdata <- object@tabdata		
		colnames(tabdata) <- c(colnames(object@data),"freq","exp")
		r <- round(tabdata[,ncol(tabdata)], digits)	
		print(cbind(tabdata[,-ncol(tabdata)],r))
		invisible(tabdata)
	}
)

