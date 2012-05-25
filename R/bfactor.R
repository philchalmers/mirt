# Class "bfactorClass"
# 
# Defines the object returned from \code{\link{bfactor}}.
# 
# 
# @name bfactorClass-class
# @aliases bfactorClass-class coef,bfactorClass-method
# fitted,bfactorClass-method print,bfactorClass-method
# residuals,bfactorClass-method show,bfactorClass-method
# summary,bfactorClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("bfactorClass", ...)}.
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass bfactorClass
# @keywords classes
setClass(
	Class = 'bfactorClass',
	representation = representation(EMiter = 'numeric', pars = 'list', 
		guess = 'numeric', parsSE='list', AIC = 'numeric', X2 = 'numeric', df = 'numeric', 
		logLik = 'numeric', p = 'numeric', F = 'matrix', h2 = 'numeric', 
		itemnames = 'character', tabdata = 'matrix', N = 'numeric', K='numeric',
		Pl = 'numeric', Theta = 'matrix', data = 'matrix', itemloc = 'numeric',
		logicalfact = 'matrix', facility = 'numeric', specific = 'numeric', tabdatalong='matrix',
		BIC = 'numeric', cormat = 'matrix', converge = 'numeric', RMSEA = 'numeric',
		par.prior = 'matrix', quadpts = 'numeric', vcov = 'matrix', Call = 'call'),	
	validity = function(object) return(TRUE)
)	

#' Full-Information Item Bifactor Analysis
#' 
#' \code{bfactor} fits a confirmatory maximum likelihood bifactor model to
#' dichotomous and polychotomous data under the item response theory paradigm. 
#' Pseudo-guessing parameters may be included but must be declared as constant.
#' 
#' 
#' 
#' \code{bfactor} follows the item factor analysis strategy explicated by
#' Gibbons and Hedeker (1992) and Gibbons et al. (2007). 
#' Nested models may be compared via an approximate
#' chi-squared difference test or by a reduction in AIC or BIC (accessible via
#' \code{\link{anova}}); note that this only makes sense when comparing class
#' \code{bfactorClass} models to class \code{mirtClass} or \code{polymirtClass}. 
#' The general equation used for item bifactor analysis in this package is in the logistic 
#' form with a scaling correction of 1.702. This correction is applied to allow
#' comparison to mainstream programs such as TESTFACT 4 (2003) and POLYFACT.
#' 
#' Unlike TESTFACT 4 (2003) initial start values are computed by using
#' information from a quasi-tetrachoric correlation matrix, potentially
#' with Carroll's (1945) adjustment for chance responses. To begin, a MINRES
#' factor analysis with one factor is extracted, and the transformed loadings
#' and intercepts (see \link{mirt} for more details) are used as starting
#' values for the general factor loadings and item intercepts. Values for the
#' specific factor loadings are taken to be half the magnitude of the extracted
#' general factor loadings. Note that while the sign of the loading may be
#' incorrect for specific factors (and possibly for some of the general factor
#' loadings) the intercepts and general factor loadings will be relatively
#' close to the final solution. These initial values should be an improvement
#' over the TESTFACT initial starting values of 1.414 for all the general
#' factor slopes, 1 for all the specific factor slopes, and 0 for all the
#' intercepts.
#' 
#' Factor scores are estimated assuming a normal prior distribution and can be
#' appended to the input data matrix (\code{full.scores = TRUE}) or displayed
#' in a summary table for all the unique response patterns. Fitted and residual
#' values can be observed by using the \code{fitted} and \code{residuals}
#' functions. To examine individuals item plots use \code{\link{itemplot}}
#' which will also plot information and surface functions.
#' Residuals are computed using the LD statistic (Chen & Thissen, 1997) in the
#' lower diagonal of the matrix returned by \code{residuals}, and Cramer's V
#' above the diagonal.
#' 
#' @aliases bfactor summary,bfactor-method coef,bfactor-method
#' fitted,bfactor-method residuals,bfactor-method
#' @param data a complete \code{matrix} or \code{data.frame} of item
#' responses that consists of only 0, 1, and \code{NA} values to be factor
#' analyzed. If scores have been recorded by the response pattern then they can
#' be recoded to dichotomous format using the \code{\link{key2binary}}
#' function.
#' @param specific a numeric vector specifying which factor loads on which
#' item. For example, if for a 4 item test with two specific factors, the first
#' specific factor loads on the first two items and the second specific factor
#' on the last two, then the vector is \code{c(1,1,2,2)}.
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#' option. Only prints patterns that have standardized residuals greater than 
#' \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param guess fixed pseudo-guessing parameter. Can be entered as a single
#' value to assign a global value or may be entered as a numeric vector for
#' each item of length \code{ncol(data)}.
#' @param SE logical; estimate parameter standard errors?
#' @param prev.cor uses a previously computed correlation matrix to be used to
#' estimate starting values for the EM estimation
#' @param par.prior a list declaring which items should have assumed priors
#' distributions, and what these prior weights are. Elements are \code{slope}
#' and \code{int} to specify the coefficients beta prior for the slopes and
#' normal prior for the intercepts, and \code{slope.items} and \code{int.items}
#' to specify which items to constrain. The value in \code{slope} is the
#' \emph{p} meta-parameter for the beta distribution (where \emph{p} > 1
#' constrains the slopes), and the two values in \code{int} are the normal
#' distribution intercept and variance. Larger values of the variance have less
#' impact on the solution. For example, if items 2 and 3 were Heywood cases
#' with no extreme item facilities, and item 4 had a very large item facility
#' (say, greater than .95) then a possible constraint might be \code{par.prior
#' = list(int = c(0,2), slope = 1.2, int.items = 4, slope.items = c(2,3))}
#' @param startvalues user declared start values for parameters
#' @param quadpts number of quadrature points per dimension. 
#' @param ncycles the number of EM iterations to be performed
#' @param tol if the largest change in the EM cycle is less than this value
#' then the EM iterations are stopped
#' @param object a model estimated from \code{bfactor} of class \code{bfactorClass}
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param digits number of significant digits to be rounded
#' @param nowarn logical; suppress warnings from dependent packages?
#' @param debug logical; turn on debugging features?
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{confmirt}},
#' \code{\link{fscores}}
#' @references
#' 
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6),
#' 1-29.
#'
#' Gibbons, R. D., & Hedeker, D. R. (1992). Full-information Item Bi-Factor
#' Analysis. \emph{Psychometrika, 57}, 423-436.
#'
#' Gibbons, R. D., Darrell, R. B., Hedeker, D., Weiss, D. J., Segawa, E., Bhaumik, D. K., 
#' Kupfer, D. J., Frank, E., Grochocinski, V. J., & Stover, A. (2007).
#' Full-Information item bifactor analysis of graded response data. 
#' \emph{Applied Psychological Measurement, 31}, 4-19
#' 
#' Carroll, J. B. (1945). The effect of difficulty and chance success on
#' correlations between items and between tests. \emph{Psychometrika, 26},
#' 347-372.
#' 
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., &
#' Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics,
#' and Full-information Item Factor Analysis [Computer software]. Lincolnwood,
#' IL: Scientific Software International.
#' @keywords models
#' @usage
#' bfactor(data, specific, guess = 0, SE = FALSE, prev.cor = NULL, par.prior = FALSE, 
#'   startvalues = NULL,  quadpts = 15, ncycles = 300, tol = .001, nowarn = TRUE, 
#'   debug = FALSE, ...)
#' 
#' \S4method{summary}{bfactor}(object, digits = 3, ...)
#' 
#' \S4method{coef}{bfactor}(object, digits = 3, ...)
#' 
#' \S4method{fitted}{bfactor}(object, digits = 3, ...)
#' 
#' \S4method{residuals}{bfactor}(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
#'
#'
#' @export bfactor
#' @examples
#' 
#' \dontrun{
#' 
#' ###load SAT12 and compute bifactor model with 3 specific factors
#' data(SAT12)
#' data <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
#' mod1 <- bfactor(data, specific)
#' coef(mod1)
#' 
#' ###Try with guessing parameters added
#' guess <- rep(.1,32)
#' mod2 <- bfactor(data, specific, guess = guess)
#' coef(mod2) 
#'
#' #########
#' #simulate data
#' a <- matrix(c(
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,0.5,NA,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5,
#' 1,NA,0.5),ncol=3,byrow=TRUE)
#' 
#' d <- matrix(c(
#' -1.0,NA,NA,
#' -1.5,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA,
#' 2.5,1.0,-1,
#' 3.0,2.0,-0.5,
#' 3.0,2.0,-0.5,
#' 3.0,2.0,-0.5,
#' 2.5,1.0,-1,
#' 2.0,0.0,NA,
#' -1.0,NA,NA,
#' -1.5,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA,
#' 1.0,NA,NA),ncol=3,byrow=TRUE)
#' 
#' sigma <- diag(3)
#' dataset <- simdata(a,d,2000,sigma)
#'
#' specific <- c(rep(1,7),rep(2,7))
#' simmod <- bfactor(dataset, specific)
#' coef(simmod)
#'
#'     }
#' 
bfactor <- function(data, specific, guess = 0, SE = FALSE, prev.cor = NULL, 
	par.prior = FALSE, startvalues = NULL, quadpts = 15, ncycles = 300, 
	tol = .001, nowarn = TRUE, debug = FALSE, ...)
{ 
	#local functions	
	fn <- function(par, rs, gues, Theta, Prior, parprior, nzeta){		
		a <- par[1:(length(par)-nzeta)]
		d <- par[(length(a)+1):length(par)]	
		rs <- rs * Prior
		if(ncol(rs) == 2){
			itemtrace <- P.mirt(a, d, Theta, gues) 
			itemtrace <- cbind(itemtrace, 1.0 - itemtrace)
		} else {
			itemtrace <- P.poly(a, d, Theta, TRUE)	
		}
		result <- (-1) * sum(rs * log(itemtrace))		
		if(parprior[1] > 1){
			sigma <- 1
			d <- sqrt(a %*% a)
			anew <- a/d
			sigma <- sigma - sum(anew)
			l <- log(sigma^(parprior[1] - 1.0) / beta(parprior[1],1.0))
			result <- result - l
		}
		if(parprior[3] > 0 && nzeta == 1){
			l <- log(dnorm(d,parprior[2],parprior[3]))
			result <- result - l
		}
		result
	}   
	
	#Main
	Call <- match.call()	
	itemnames <- colnames(data)
	data <- as.matrix(data)
	data.original <- data		
	if(!any(data %in% c(0:20,NA))) 
		stop("Data must contain only numeric values (including NA).")	
	J <- ncol(data)
	N <- nrow(data)	
	if(length(guess) == 1) guess <- rep(guess,J)
	colnames(data) <- itemnames
	if(length(guess) > J || length(guess) < J) 
		stop("The number of guessing parameters is incorrect.")
	facility <- colMeans(na.omit(data))		
	uniques <- list()
	for(i in 1:J)
		uniques[[i]] <- sort(unique(data[,i]))
	K <- rep(0,J)
	for(i in 1:J) K[i] <- length(uniques[[i]])	
	guess[K > 2] <- 0	
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
	pats <- apply(fulldata, 1, paste, collapse = "/") 
	freqs <- table(pats)
	nfreqs <- length(freqs)
	r <- as.vector(freqs)	
	tabdata <- unlist(strsplit(cbind(names(freqs)), "/"))
	tabdata <- matrix(as.numeric(tabdata), nfreqs, sum(K), TRUE)	
	tabdata <- cbind(tabdata,r) 
	colnames(tabdata) <- c(Names,'Freq')	
	#for return
	pats <- apply(data, 1, paste, collapse = "/") 
	freqs <- table(pats)		
	tabdata2 <- unlist(strsplit(cbind(names(freqs)), "/"))
	tabdata2 <- matrix(as.numeric(tabdata2), nfreqs, J, TRUE)	
	tabdata2 <- cbind(tabdata2,r) 
	colnames(tabdata2) <- c(itemnames,'Freq')
	if(is.logical(par.prior)) 
	    if(par.prior) suppressAutoPrior <- FALSE  
	        temp <- matrix(c(1,0,0),ncol = 3, nrow=J, byrow=TRUE)
	if(!is.logical(par.prior)){
		if(!is.null(par.prior$slope.items))
			for(i in 1:length(par.prior$slope.items))
				temp[par.prior$slope.items[i],1] <- par.prior$slope		
		if(!is.null(par.prior$int.items))
			for(i in 1:length(par.prior$int.items))
				temp[par.prior$int.items[i],2:3] <- par.prior$int		 
	}  
	par.prior <- temp 
	if(!is.null(prev.cor)){
		if (ncol(prev.cor) == nrow(prev.cor)) Rpoly <- prev.cor
			else stop("Correlation matrix is not square.\n")
	} else Rpoly <- cormod(na.omit(data.original),K,guess)	
	if(det(Rpoly) < 1e-15) Rpoly <- cor(na.omit(data.original))
	FA <- suppressWarnings(psych::fa(Rpoly,1,rotate = 'none', warnings= FALSE, fm="minres"))	
	loads <- unclass(loadings(FA))
	u <- FA$unique
	u[u < .1 ] <- .25	
	cs <- sqrt(u)
	lambdas <- loads/cs	
	slambdas <- matrix(0, nrow = J, ncol = length(unique(specific)))
	logicalfact <- matrix(FALSE, nrow = J, ncol = ncol(slambdas) + 1)
	for(i in 1:J){
		temp <- rep(FALSE, ncol(slambdas))
		slambdas[i, specific[i]] <- lambdas[i] / 2
		temp[specific[i]] <- TRUE
		logicalfact[i, ] <- c(TRUE,temp)
	}	
	lambdas <- cbind(lambdas, slambdas)	
    zetas <- list()	
    for(i in 1:J){
        if(K[i] == 2){
            zetas[[i]] <- qnorm(mean(fulldata[,itemloc[i]]))/cs[i]            			
        } else {
            temp <- table(data[,i])[1:(K[i]-1)]/N
            temp <- cumsum(temp)			
            zetas[[i]] <- qnorm(1 - temp)/cs[i]        			
        }       
    }    		
	pars <- list(lambdas=lambdas, zetas=zetas)
	attr(pars, 'lamsel') <- logicalfact
	npars <- length(unlist(pars))	

	lastpars2 <- lastpars1 <- pars 	
	theta <- as.matrix(seq(-4, 4, length.out = quadpts))
	Theta <- thetaComb(theta, 2)
	prior <- dnorm(theta) 
	Prior <- mvtnorm::dmvnorm(Theta,rep(0,2),diag(2))  
	startvalues <- pars  
	converge <- 1
	problemitems <- c()
	index <- 1:J
	nfact <- ncol(lambdas)
	temp <- matrix(0,nrow=J,ncol=(nfact-1))
	sitems <- matrix(0, nrow=sum(K), ncol=(nfact-1))
	for(i in 1:J) temp[i,specific[i]] <- 1
	ind <- 1
	for(i in 1:J){
		for(j in 1:K[i]){
			sitems[ind, ] <- temp[i, ]
			ind <- ind + 1
		}		
	}		
	if(debug) print(startvalues)			 		
	
	#EM  loop  
	for (cycles in 1:ncycles) 
	{    
		rlist <- Estep.bfactor(pars, tabdata, Theta, prior, guess, 
			specific, sitems, itemloc)
		if(debug) print(sum(r * log(rlist$expected)))			
		lastpars2 <- lastpars1
		lastpars1 <- pars			
		for(i in 1:J){ 
			par <- c(pars$lambdas[i, logicalfact[i, ]], pars$zetas[[i]])
			itemsel <- c(itemloc[i]:(itemloc[i+1] - 1))							
			maxim <- try(optim(par, fn=fn, rs=rlist$r1[, itemsel], gues=guess[i], Theta=Theta, 
				Prior=Prior, parprior=par.prior[i, ], nzeta=K[i]-1, control=list(maxit=25)))			
			if(class(maxim) == "try-error") {
				problemitems <- c(problemitems, i)	  
				converge <- 0
				next
			}	
			pars$lambdas[i, logicalfact[i, ]] <- maxim$par[1:2]
			pars$zetas[[i]] <- maxim$par[3:length(par)]	  
		}
		maxdif <- max(abs(unlist(lastpars1) - unlist(pars)))	
		if (maxdif < tol && cycles > 5) break 	
		# apply rate acceleration every third cycle    
		if (cycles %% 3 == 0 & cycles > 6)		 
			pars <- rateChange(pars, lastpars1, lastpars2)       
	}
	
	if(any(par.prior[,1] != 1)) cat("Slope prior for item(s):",
		as.character(index[par.prior[,1] > 1]), "\n")
	if(any(par.prior[,3] != 0)) cat("Intercept prior for item(s):",
		as.character(index[par.prior[,3] > 0]), "\n")  
	if(converge == 0) 
		warning("Parameter estimation reached unacceptable values. 
		Model probably did not converge.")
	if(length(problemitems) > 0) warning("Problem with the M-step for item(s): ", 
		paste(unique(problemitems), " "))	
	lastchange <- unlist(lastpars1) - unlist(pars)
	if (cycles == ncycles){ 
		converge <- 0
		message("Estimation terminated after ", cycles, " EM loops. Maximum changes: 
			\n slopes = ", round(max(abs(lastchange[,1:nfact])),4), ", intercepts = ", 
			round(max(abs(lastchange[,ncol(pars)])),4) ,"\n")
	}	
	rlist <- Estep.bfactor(pars, tabdata, Theta, prior, guess, 
			specific, sitems, itemloc)
	Pl <- rlist$expected
	logLik <- sum(r * log(Pl))
	vcovpar <- matrix(999)
	parsSE <- list()
	if(SE){		
		LLfun <- function(p, pars, tabdata, Theta, prior, guess, specific, sitems, itemloc){
			pars2 <- rebuildPars(p, pars)		
			rlist <- Estep.bfactor(pars2, tabdata, Theta, prior, guess, 
				specific, sitems, itemloc)    	  
			Pl <- rlist$expected
			logLik <- sum(r*log(Pl))
			-1*logLik		
		}
		fmin <- nlm(LLfun, unlist(pars), pars=pars,tabdata=tabdata,Theta=Theta,prior=prior,
			guess=guess, specific=specific, sitems=sitems, itemloc=itemloc, hessian=TRUE, gradtol=.1)		
		vcovpar <- solve(fmin$hessian)
		parsSE <- rebuildPars(sqrt(diag(vcovpar)), pars)	
	}
	logN <- 0
	logr <- rep(0,length(r))
	for (i in 1:N) logN <- logN + log(i)
	for (i in 1:length(r)) 
		for (j in 1:r[i]) 
			logr[i] <- logr[i] + log(j)	
	logLik <- logLik + logN/sum(logr)
	AIC <- (-2) * logLik + 2 * npars
	BIC <- (-2) * logLik + npars*log(N)
	X2 <- 2 * sum(r * log(r / (N*Pl)))  
	df <- length(r) - 1 + nfact*(nfact - 1)/2 - npars 
	p <- 1 - pchisq(X2,df)
	if(any(is.na(data.original))) p <- RMSEA <- X2 <- NaN
	RMSEA <- ifelse((X2 - df) > 0, 
	    sqrt(X2 - df) / sqrt(df * (N-1)), 0)

	#from last EM cycle pars to FA
	norm <- sqrt(1 + rowSums(pars$lambdas[ ,1:nfact]^2))	 
	F <- matrix(0,ncol = nfact, nrow = J)
	for (i in 1:J) 
		F[i,1:nfact] <- pars$lambdas[i,1:nfact]/norm[i]  
	colnames(F) <- c('G',paste("F_", 1:(ncol(F)-1),sep=""))
	h2 <- rowSums(F^2)  

	mod <- new('bfactorClass',EMiter=cycles, pars=pars, guess=guess, AIC=AIC, X2=X2, 
		parsSE=parsSE, df=df, logLik=logLik, p=p, F=F, h2=h2, itemnames=itemnames, BIC=BIC,
		tabdata=tabdata2, N=N, Pl=Pl, Theta=Theta, data=data.original, tabdatalong=tabdata, 
		logicalfact=logicalfact, facility=facility, specific=specific, itemloc=itemloc,
		cormat=Rpoly, converge=converge, par.prior=par.prior, quadpts=quadpts,
		vcov=vcovpar, RMSEA=RMSEA, K=K, Call=Call)  
	return(mod)  
} 

# Methods

setMethod(
	f = "print",
	signature = signature(x = 'bfactorClass'),
	definition = function(x, ...){  
		cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")		
		cat("Full-information bifactor analysis with ", 
			length(unique(x@specific)), " specific factors \n", sep='')
		if(x@converge == 1)	
			cat("Converged in ", x@EMiter, " iterations using ",x@quadpts,
			" quadrature points. \n", sep="")
		else 	
			cat("Estimation stopped after ", x@EMiter, " iterations using ",x@quadpts,
			" quadrature points. \n", sep="")
		cat("Log-likelihood = ", x@logLik, "\n")
		cat("AIC = ", x@AIC, "\n")		
		cat("BIC = ", x@BIC, "\n")
		if(!is.nan(x@p))
			cat("G^2 = ", round(x@X2,2), ", df = ", 
				x@df, ", p = ", round(x@p,4), ", RMSEA = ", round(x@RMSEA,3), "\n", sep="")
		else 
			cat("G^2 = ", NA, ", df = ", 
				x@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="")		
	}
)

setMethod(
	f = "show",
	signature = signature(object = 'bfactorClass'),
	definition = function(object){  
		cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")		
		cat("Full-information bifactor analysis with ", 
			length(unique(object@specific)), " specific factors \n", sep='')
		if(object@converge == 1)	
			cat("Converged in ", object@EMiter, " iterations using ", object@quadpts,
				" quadrature points.\n", sep="")
		else 	
			cat("Estimation stopped after ", object@EMiter, " iterations using ", 
				object@quadpts,	" quadrature points.\n", sep="")
		cat("Log-likelihood = ", object@logLik, "\n")
		cat("AIC = ", object@AIC, "\n")		
		cat("BIC = ", object@BIC, "\n")
		if(!is.nan(object@p))
			cat("G^2 = ", round(object@X2,2), ", df = ", 
				object@df, ", p = ", round(object@p,4), ", RMSEA = ", round(object@RMSEA,3),
                "\n", sep="")
		else 
			cat("G^2 = ", NA, ", df = ", 
				object@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="")			
	}
)

setMethod(
	f = "summary",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, digits = 3, ...){
		F <- round(object@F,digits)
		SS <- colSums(F^2)	
		F[!object@logicalfact] <- NA
		h2 <- round(object@h2,digits)					
		names(h2) <- "h2"		
		loads <- round(cbind(F,h2),digits)
		rownames(loads) <- colnames(object@data)	 
		cat("\nFactor loadings: \n\n")
		print(loads)
		cat("\nSS loadings: ",round(SS,digits), "\n")
		cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
		if(any(h2 > 1)) 
			warning("Solution has heywood cases. Interpret with caution.")
	}
)

setMethod(
	f = "coef",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, digits = 3, ...){
		K <- object@K
		a <- object@pars$lambdas			
		d <- matrix(NA, nrow(a), max(K-1))
		zetas <- object@pars$zetas
		for(i in 1:length(K)){
			d[i, 1:(K[i] - 1)] <- zetas[[i]]
		}
		A <- sqrt(apply(a^2,1,sum))
		B <- -d/A 
		a[!attr(object@pars,'lamsel')] <- NA	
		parameters <- round(cbind(a,d,object@guess,A,B),digits)
		colnames(parameters) <- c('a_G',paste("a_", 1:(ncol(object@F)-1),sep=""),
			paste("d_", 1:(max(K)-1),sep=""), "guess", "mvdisc", paste("mvint_", 1:(max(K)-1),sep=""))  
		rownames(parameters) <- colnames(object@data)	
		cat("\nParameters with multivariate discrimination and intercept: \n\n")		
		print(parameters)
		ret <- list(parameters)
		if(length(object@parsSE) > 1){
			cat("\nStd. Errors: \n\n")	
			a <- object@parsSE$lambdas			
			d <- matrix(NA, nrow(a), max(K-1))
			zetas <- object@parsSE$zetas
			for(i in 1:length(K)){
				d[i, 1:(K[i] - 1)] <- zetas[[i]]
			}
			SEs <- cbind(a,d)
			colnames(SEs) <- c('a_G',paste("a_", 1:(ncol(object@F)-1),sep=""),
				paste("d_", 1:(max(K)-1),sep=""))
			rownames(SEs) <- colnames(object@data)		
			print(SEs, digits)
			ret <- list(parameters, SEs)
		}	
		invisible(ret)
	}
)

setMethod(
	f = "residuals",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
	{       
		K <- object@K
		lf <- attr(object@pars, 'lamsel')
		Theta <- object@Theta
		data <- object@data	
		N <- nrow(data)	
		J <- ncol(data)		
		lambdas <- object@pars$lambdas
		zetas <- object@pars$zetas
		guess <- object@guess
		guess[is.na(guess)] <- 0
		itemloc <- object@itemloc
		res <- matrix(0,J,J)
		diag(res) <- NA
		colnames(res) <- rownames(res) <- colnames(data)
		prior <- mvtnorm::dmvnorm(Theta,rep(0,2),diag(2))
		prior <- prior/sum(prior)	
		if(restype == 'LD'){	
			for(i in 1:J){								
				for(j in 1:J){			
					if(i < j){
						P1 <- P.bfactor(lambdas[i, ], zetas[[i]], Theta, guess[i], lf[i, ])
						P2 <- P.bfactor(lambdas[j, ], zetas[[j]], Theta, guess[j], lf[j, ])
						if(K[i] == 2) P1 <- cbind(1-P1, P1)
						if(K[j] == 2) P2 <- cbind(1-P2, P2)						
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
			r <- object@tabdata[ ,ncol(object@tabdata)]
			res <- round((r - object@Pl * nrow(object@data)) / 
				sqrt(object@Pl * nrow(object@data)),digits)
			expected <- round(N * object@Pl/sum(object@Pl),digits)  
			tabdata <- object@tabdata
			freq <- tabdata[ ,ncol(tabdata)]			
			tabdata[tabdata[ ,1:ncol(object@data)] == 99] <- NA
			tabdata[ ,ncol(tabdata)] <- freq
			tabdata <- cbind(tabdata,expected,res)
			colnames(tabdata) <- c(colnames(object@tabdata),"exp","res")	
			if(!is.null(printvalue)){
				if(!is.numeric(printvalue)) stop('printvalue is not a number.')
				tabdata <- tabdata[abs(tabdata[ ,ncol(tabdata)]) > printvalue, ]
			}			
			return(tabdata)				
		}
	}
)

setMethod(
	f = "fitted",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, digits = 3, ...){  
		expected <- round(object@N * object@Pl/sum(object@Pl),digits)  
		tabdata <- object@tabdata
		freq <- tabdata[ ,ncol(tabdata)]
		tabdata[tabdata[ ,1:ncol(object@data)] == 9] <- NA
		tabdata[ ,ncol(tabdata)] <- freq
		tabdata <- cbind(tabdata,expected)
		colnames(tabdata) <- c(colnames(object@tabdata)[1:(ncol(tabdata)-2)],"freq","exp")		
		tabdata
	}
)



