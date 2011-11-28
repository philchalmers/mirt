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
	representation = representation(EMiter = 'numeric', pars = 'matrix', guess = 'numeric', 
		AIC = 'numeric', X2 = 'numeric', df = 'numeric', log.lik = 'numeric', p = 'numeric', 
		F = 'matrix', h2 = 'numeric', itemnames = 'character', tabdata = 'matrix', 
		N = 'numeric', Pl = 'numeric', Theta = 'matrix', fulldata = 'matrix', 
		logicalfact = 'matrix', facility = 'numeric', specific = 'numeric', BIC = 'numeric',
		cormat = 'matrix', converge = 'numeric', par.prior = 'matrix', quadpts = 'numeric', 
		Call = 'call'),	
	validity = function(object) return(TRUE)
)	

#' Full-Information Item Bifactor Analysis
#' 
#' \code{bfactor} fits a confirmatory maximum likelihood bifactor model to
#' dichotomous data under the item response theory paradigm. Pseudo-guessing
#' parameters may be included but must be declared as constant, since the
#' estimation of these parameters often leads to unacceptable solutions.
#' Missing values are automatically assumed to be 0.
#' 
#' 
#' \code{bfactor} follows the item factor analysis strategy explicated by
#' Gibbons and Hedeker (1992). Nested models may be compared via an approximate
#' chi-squared difference test or by a reduction in AIC or BIC (accessible via
#' \code{\link{anova}}); note that this only makes sense when comparing class
#' \code{bfactorClass} models to class \code{mirtClass}. The general equation
#' used for item bifactor analysis in this package is in the logistic form with
#' a scaling correction of 1.702. This correction is applied to allow
#' comparison to mainstream programs such as TESTFACT 4 (2003).
#' 
#' Unlike TESTFACT 4 (2003) initial start values are computed by using
#' information from the matrix of quasi-tetrachoric correlations, potentially
#' with Carroll's (1945) adjustment for chance responses. To begin, a MINRES
#' factor analysis with one factor is extracted, and the transformed loadings
#' and intercepts (see \link{mirt} for more details) are used as starting
#' values for the general factor loadings and item intercepts. Values for the
#' specific factor loadings are taken to be half the magnitude of the extracted
#' general factor loadings. Note that while the sign of the loading may be
#' incorrect for specific factors (and possibly for some of the general factor
#' loadings) the intercepts and general factor loadings will be relatively
#' close to the final solution. These initial values should be an improvement
#' over the TESTFACT 4 initial starting values of 1.414 for all the general
#' factor slopes, 1 for all the specific factor slopes, and 0 for all the
#' intercepts.
#' 
#' Factor scores are estimated assuming a normal prior distribution and can be
#' appended to the input data matrix (\code{full.scores = TRUE}) or displayed
#' in a summary table for all the unique response patterns. Fitted and residual
#' values can be observed by using the \code{fitted} and \code{residuals}
#' functions. To examine individuals item plots use \code{\link{itemplot}}
#' which will also plot information and surface functions (although the
#' \code{\link[plink]{plink}} package may be more suitable for IRT graphics).
#' Residuals are computed using the LD statistic (Chen & Thissen, 1997) in the
#' lower diagonal of the matrix returned by \code{residuals}, and Cramer's V
#' above the diagonal.
#' 
#' @aliases bfactor summary,bfactor-method coef,bfactor-method
#' fitted,bfactor-method residuals,bfactor-method
#' @param fulldata a complete \code{matrix} or \code{data.frame} of item
#' responses that consists of only 0, 1, and \code{NA} values to be factor
#' analyzed. If scores have been recorded by the response pattern then they can
#' be recoded to dichotomous format using the \code{\link{key2binary}}
#' function.
#' @param specific a numeric vector specifying which factor loads on which
#' item. For example, if for a 4 item test with two specific factors, the first
#' specific factor loads on the first two items and the second specific factor
#' on the last two, then the vector is \code{c(1,1,2,2)}.
#' @param guess fixed pseudo-guessing parameter. Can be entered as a single
#' value to assign a global value or may be entered as a numeric vector for
#' each item of length \code{ncol(fulldata)}.
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
#' @param quadpts number of quadrature points per dimension. If \code{NULL}
#' then the number of quadrature points is set to 15
#' @param ncycles the number of EM iterations to be performed
#' @param EMtol if the largest change in the EM cycle is less than this value
#' then the EM iteration are stopped early
#' @param object a model estimated from \code{bfactor} of class \code{bfactor}
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
#' Gibbons, R. D., & Hedeker, D. R. (1992). Full-information Item Bi-Factor
#' Analysis. \emph{Psychometrika, 57}, 423-436.
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
#' bfactor(fulldata, specific, guess = 0, prev.cor = NULL, par.prior = FALSE, 
#'   startvalues = NULL,  quadpts = NULL, ncycles = 300, EMtol = .001, nowarn = TRUE, 
#'   debug = FALSE, ...)
#' 
#' \S4method{summary}{bfactor}(object, digits = 3, ...)
#' 
#' \S4method{coef}{bfactor}(object, digits = 3, ...)
#' 
#' \S4method{fitted}{bfactor}(object, digits = 3, ...)
#' 
#' \S4method{residuals}{bfactor}(object, restype = 'LD', digits = 3, ...)
#'
#'
#' @export bfactor
#' @examples
#' 
#' \dontrun{
#' 
#' ###load SAT12 and compute bifactor model with 3 specific factors
#' data(SAT12)
#' fulldata <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
#' mod1 <- bfactor(fulldata, specific)
#' coef(mod1)
#' 
#' ###Try with guessing parameters added
#' guess <- rep(.1,32)
#' mod2 <- bfactor(fulldata, specific, guess = guess)
#' coef(mod2) #item 32 too difficult to include guessing par
#' 
#' #fix by imposing a weak intercept prior
#' mod3a <- bfactor(fulldata, specific, guess = guess, par.prior =
#'     list(int = c(0,4), int.items = 32))
#' coef(mod3a)
#' 
#' #...or by removing guessing parameter
#' guess[32] <- 0
#' mod3b <- bfactor(fulldata, specific, guess = guess)
#' coef(mod3b)
#'     }
#' 
bfactor <- function(fulldata, specific, guess = 0, prev.cor = NULL, 
  par.prior = FALSE, startvalues = NULL, quadpts = NULL, ncycles = 300, 
  EMtol = .001, nowarn = TRUE, debug = FALSE, ...)
{ 
	#local functions
	fn <- function(pars, r1, N, guess, Theta, prior, parprior){
		a <- pars[1:(length(pars)-1)]
		d <- pars[length(pars)]
		r1 <- r1 * prior	
		N <- N * prior
		result <- .Call("loglik", 	                
						as.double(a),				
						as.double(d),
						as.double(r1),
						as.double(N),
						as.double(guess),
						as.double(as.matrix(Theta)),
						as.integer(parprior))					
	}  
	gr <- function(pars, r1, N, guess, Theta, prior, parprior){
		a <- pars[1:(length(pars)-1)]
		d <- pars[length(pars)]			
		result <- .Call("grad", 	                
						as.double(a),				
						as.double(d),
						as.double(r1),
						as.double(N),
						as.double(guess),
						as.double(as.matrix(Theta)),
						as.double(prior),
						as.integer(parprior))	    				
	} 
	
	#Main
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
	if(2*nfact >= nitems) stop('Model is not identified.')
	if (length(guess) == 1) guess <- rep(guess,nitems)
		else if (length(guess) > nitems || length(guess) < nitems) 
			stop("The number of guessing parameters is incorrect.")
	pats <- apply(fulldata, 1, paste, collapse = "/") 
	freqs <- table(pats)
	nfreqs <- length(freqs)
	r <- as.vector(freqs)
	N <- nrow(fulldata)
	K <- rep(2,nitems)  
	tabdata <- unlist(strsplit(cbind(names(freqs)), "/"))
	tabdata <- matrix(as.numeric(tabdata), nfreqs, nitems, TRUE)
	tabdata <- cbind(tabdata,r)  
	logicalfact <- matrix(FALSE,nitems,nfact - 1)
	is.na(specific) <- FALSE
	for (i in 1:nitems) logicalfact[i,specific[i]] <- TRUE
	logicalfact <- cbind(rep(TRUE,nitems),logicalfact)     
	if (is.null(quadpts)) quadpts <- 15
	theta <- as.matrix(seq(-4,4,length.out = quadpts))
	Theta <- as.matrix(expand.grid(theta,theta))
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
		} else Rpoly <- cormod(fulldata,K,guess)       
	pars <- matrix(0,nrow=nitems, ncol=nfact + 1)
	if (is.null(startvalues)){  
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
	prior <- dnorm(theta) 
	Prior <- dmvnorm(Theta,rep(0,2),diag(2))  
	startvalues <- pars  
	converge <- 1
	problemitems <- c()
	index <- 1:nitems
	sitems <- matrix(0,ncol=nitems,nrow=(nfact-1))
	for(i in 1:nitems) sitems[specific[i],i] <- 1      
	if(debug){
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
		for(i in 1:nitems){ 
			if(guess[i] == 0)	
				maxim <- try(optim(mpars[i, ],fn=fn,gr=gr,r1=rlist[[1]][i, ],N=rlist[[2]][i, ],
					guess=guess[i],Theta=Theta,prior=Prior,parprior=par.prior[i, ],method="BFGS"))
			else	  
				maxim <- try(optim(mpars[i, ],fn=fn,r1=rlist[[1]][i, ],N=rlist[[2]][i, ],
					guess=guess[i],Theta=Theta,prior=Prior,parprior=par.prior[i, ],method="BFGS"))
			if(class(maxim) == "try-error") {
				problemitems <- c(problemitems, i)	  
				converge <- 0
				next
			}	
			mpars[i, ] <- maxim$par	  
		}			
		sload <- mpars[(nitems+1):(2*nitems)]		
		for(i in 1:nitems){
			temp <- selvec[pars[i,2:nfact] > 0]
			pars[i,temp] <- sload[i]
		}  
		pars[ ,1] <- mpars[,1]
		pars[ ,nfact+1] <- mpars[ ,3]		
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
		if (cycles %% 3 == 0 & cycles > 6){
			d1 <- lastpars1 - pars
			d2 <- lastpars2 - pars      
			for (i in 1:nitems) {
				for(j in 1:ncol(pars)){      
					if((abs(d1[i,j]) > 0.001) & (d1[i,j]*d2[i,j] > 0.0) & (d1[i,j]/d2[i,j] < 1.0))
						rate[i,j] <- (1 - (1 - rate[i,j]) * (d1[i,j]/d2[i,j]))
						else rate[i,j] <- 0
				}        
			}      
		}  
		rate[pars > 4] <- 0
		rate[pars < -4] <- 0
		pars <- lastpars1*rate*(-2) + (1 - rate*(-2))*pars          
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
	for (i in 1:N) logN <- logN + log(i)
	for (i in 1:length(r)) 
		for (j in 1:r[i]) 
			logr[i] <- logr[i] + log(j)	
	log.lik <- log.lik + logN/sum(logr)
	AIC <- (-2) * log.lik + 6 * length(specific)
	BIC <- (-2) * log.lik + 3 * length(specific)*log(N)
	X2 <- 2 * sum(r * log(r / (N*Pl)))  
	df <- length(r) - 1 - 2*nitems - length(specific)
	p <- 1 - pchisq(X2,df)

	#from last EM cycle pars to FA
	norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2))
	gam <- (-1)*pars[ ,nfact + 1]/norm  
	F <- matrix(0,ncol = nfact, nrow = nitems)
	for (i in 1:nitems) F[i,1:nfact] <- pars[i,1:nfact]/norm[i]  
	h2 <- rowSums(F^2)  

	mod <- new('bfactorClass',EMiter=cycles, pars=pars, guess=guess, AIC=AIC, X2=X2, 
		df=df, log.lik=log.lik, p=p, F=F, h2=h2, itemnames=itemnames, BIC=BIC,
		tabdata=tabdata, N=N, Pl=Pl, Theta=Theta, fulldata=fulldata, 
		logicalfact=logicalfact, facility=facility, specific=specific,
		cormat=Rpoly, converge=converge, par.prior=par.prior, quadpts=quadpts,Call=Call)  
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
		cat("Log-likelihood = ", x@log.lik, "\n")
		cat("AIC = ", x@AIC, "\n")		
		cat("BIC = ", x@BIC, "\n")
		cat("G^2 = ", round(x@X2,2), ", df = ", 
		x@df, ", p = ", round(x@p,4), "\n")
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
		cat("Log-likelihood = ", object@log.lik, "\n")
		cat("AIC = ", object@AIC, "\n")		
		cat("BIC = ", object@BIC, "\n")
		cat("G^2 = ", round(object@X2,2), ", df = ", 
		object@df, ", p = ", round(object@p,4), "\n")
	}
)

setMethod(
	f = "summary",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, digits = 3, ...){
		F <- round(object@F,digits)
		h2 <- round(object@h2,digits)
		SS <- colSums(F^2)		
		colnames(F) <- c('G',paste("F_", 1:(ncol(F)-1),sep=""))
		names(h2) <- "h2"		
		loads <- round(cbind(F,h2),digits)
		rownames(loads) <- object@itemnames  
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
		a <- as.matrix(object@pars[ ,1:(ncol(object@pars)-1)])
		d <- object@pars[ ,ncol(object@pars)]
		A <- sqrt(apply(a^2,1,sum))
		B <- -d/A 
		fac <- object@facility  
		parameters <- round(cbind(object@pars,object@guess,fac,A,B),digits)
		colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),"d", "guess", 
			"facility","mvdisc", "mvint")  
		cat("\nParameters with multivariate discrimination and intercept: \n\n")
		print(parameters)	    
		invisible(parameters)
	}
)

setMethod(
	f = "residuals",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, restype = 'LD', digits = 3, p = FALSE, ...)
	{       
		Theta <- object@Theta
		fulldata <- object@fulldata	
		N <- nrow(fulldata)	
		J <- ncol(fulldata)
		nfact <- ncol(object@F)
		lambdas <- matrix(object@pars[,1:nfact], J)
		zetas <- object@pars[,(nfact+1)]
		guess <- object@guess
		guess[is.na(guess)] <- 0
		logicalfact <- object@logicalfact
		if(restype == 'LD'){
			res <- matrix(0,J,J)
			diag(res) <- NA
			colnames(res) <- rownames(res) <- colnames(fulldata)
			prior <- dmvnorm(Theta,rep(0,2),diag(2))
			prior <- prior/sum(prior)
			for(i in 1:J){			
				for(j in 1:J){
					if(i < j){
						P1 <- P.bfactor(lambdas[i,],zetas[i], Theta, guess[i],logicalfact[i,])
						P2 <- P.bfactor(lambdas[j,],zetas[j], Theta, guess[j],logicalfact[j,])
						E22 <- N * sum(P1 * P2 * prior)
						E12 <- N * sum(P1 * (1-P2) * prior)
						E21 <- N * sum((1-P1) * P2 * prior)
						E11 <- N * sum((1-P1) * (1-P2) * prior)
						tab <- table(fulldata[,i],fulldata[,j])
						Etab <- matrix(c(E11,E12,E21,E22),2)
						s <- phi(tab) - phi(Etab)
						if(s == 0) s <- 1
						res[j,i] <- sum(((tab - Etab)^2)/Etab) * sign(s)
						res[i,j] <- sqrt( abs(res[j,i]) / N ) 
					}
				}
			}
			cat("\nLD matrix:\n\n")			
			res <- round(res,digits)	
			return(res)	
		}
		if(restype == 'exp'){
			r <- object@tabdata[ ,ncol(object@tabdata)]
			res <- round((r - object@Pl * nrow(object@fulldata)) / 
				sqrt(object@Pl * nrow(object@fulldata)),digits)
			expected <- round(object@N * object@Pl/sum(object@Pl),digits)  
			tabdata <- cbind(object@tabdata,expected,res)
			colnames(tabdata) <- c(object@itemnames, "freq", "exp", "std_res")								
			return(tabdata)
		}				
	}
)

setMethod(
	f = "fitted",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, digits = 3, ...){  
		expected <- round(object@N * object@Pl/sum(object@Pl),digits)  
		tabdata <- cbind(object@tabdata,expected)
		colnames(tabdata) <- c(object@itemnames, "freq", "exp")	
		print(tabdata)
		invisible(tabdata)
	}
)



