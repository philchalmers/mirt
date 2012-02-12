# Class "mirtClass"
# 
# Defines the object returned from \code{\link{mirt}}.
# 
# 
# @name mirtClass-class
# @aliases mirtClass-class anova,mirtClass-method coef,mirtClass-method
# fitted,mirtClass-method plot,mirtClass,missing-method print,mirtClass-method
# residuals,mirtClass-method show,mirtClass-method summary,mirtClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("mirtClass", ...).}.
# @method Emiter number of EM iterations
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass mirtClass
# @keywords classes
setClass(
	Class = 'mirtClass',
	representation = representation(EMiter = 'numeric', pars = 'matrix', guess = 'numeric', 
		X2 = 'numeric', df = 'numeric', p = 'numeric', AIC = 'numeric', logLik = 'numeric',
		F = 'matrix', h2 = 'numeric', tabdata = 'matrix', Theta = 'matrix', Pl = 'numeric',
		fulldata = 'matrix', cormat = 'matrix', facility = 'numeric', converge = 'numeric', 
		quadpts = 'numeric', BIC = 'numeric', vcov = 'matrix', Call = 'call'),	
	validity = function(object) return(TRUE)
)	

#' Full-Information Item Factor Analysis (Multidimensional Item Response
#' Theory)
#' 
#' \code{mirt} fits an unconditional maximum likelihood factor analysis model
#' to dichotomous data under the item response theory paradigm. Pseudo-guessing
#' parameters may be included but must be declared as constant, since the
#' estimation of these parameters often leads to unacceptable solutions.
#' 
#' 
#' 
#' \code{mirt} follows the item factor analysis strategy by marginal maximum
#' likelihood estimation (MML) outlined in Bock and Aiken (1981) and Bock,
#' Gibbons and Muraki (1988). Nested models may be compared via the approximate
#' chi-squared difference test or by a reduction in AIC/BIC values (comparison
#' via \code{\link{anova}}). The general equation used for dichotomous
#' multidimensional item response theory item is a logistic form with a scaling
#' correction of 1.702. This correction is applied to allow comparison to
#' mainstream programs such as TESTFACT (2003). The general IRT equation is
#' 
#' \deqn{P(X | \theta; \bold{a}_i; d_i; g_i) = g_j + (1 - g_j) / (1 +
#' exp(-1.702(\bold{a}_j' \theta + d_j)))}
#' 
#' where \emph{j} is the item index, \eqn{\bold{a}_j} is the vector of
#' discrimination parameters (i.e., slopes), \deqn{\theta} is the vector of
#' factor scores, \eqn{d_j} is the intercept, and \eqn{g_j} is the
#' pseudo-guessing parameter. To avoid estimation difficulties the \eqn{g_j}'s
#' must be specified by the user.
#' 
#' Estimation begins by computing a matrix of quasi-tetrachoric correlations,
#' potentially with Carroll's (1945) adjustment for chance responds. A MINRES
#' factor analysis with \code{nfact} is then extracted and item parameters are
#' estimated by \eqn{a_{ij} = f_{ij}/u_j}, where \eqn{f_{ij}} is the factor
#' loading for the \emph{j}th item on the \emph{i}th factor, and \eqn{u_j} is
#' the square root of the factor uniqueness, \eqn{\sqrt{1 - h_j^2}}. The
#' initial intercept parameters are determined by calculating the inverse
#' normal of the item facility (i.e., item easiness), \eqn{q_j}, to obtain
#' \eqn{d_j = q_j / u_j}. Following these initial estimates the model is
#' iterated using the EM estimation strategy with fixed quadrature points.
#' Implicit equation accelerations described by Ramsey (1975) are also added to
#' facilitate parameter convergence speed, and these are adjusted every third
#' cycle.
#' 
#' Factor scores are estimated assuming a normal prior distribution and can be
#' appended to the input data matrix (\code{full.data = TRUE}) or displayed in
#' a summary table for all the unique response patterns. \code{summary} allows
#' for various rotations available from the \code{GPArotation} package. These
#' are:
#' 
#' \describe{ \item{orthogonal: }{\code{"varimax", "quartimax", "tandemI",
#' "tandemII", "entropy", "mccammon"}} \item{oblique: }{\code{"promax",
#' "oblimin", "quartimin", "oblimax", "simplimax"}} }
#' 
#' Using \code{plot} will plot the either the test surface function or the test
#' information function for 1 and 2 dimensional solutions. To examine
#' individual item plots use \code{\link{itemplot}} (although the
#' \code{\link[plink]{plink}} package may be more suitable for IRT graphics)
#' which will also plot information and surface functions. Residuals are
#' computed using the LD statistic (Chen \& Thissen, 1997) in the lower
#' diagonal of the matrix returned by \code{residuals}, and Cramer's V above
#' the diagonal.
#' 
#' @aliases mirt summary,mirt-method coef,mirt-method anova,mirt-method
#' fitted,mirt-method plot,mirt-method residuals,mirt-method
#' @param fulldata a \code{matrix} or \code{data.frame} that consists of only
#' 0, 1, and \code{NA} values to be factor analyzed. If scores have been
#' recorded by the response pattern then they can be recoded to dichotomous
#' format using the \code{\link{key2binary}} function
#' @param nfact number of factors to be extracted
#' @param SE logical, estimate the standard errors?
#' @param guess fixed pseudo-guessing parameters. Can be entered as a single
#' value to assign a global guessing parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param prev.cor use a previously computed correlation matrix to be used to
#' estimate starting values for the EM estimation? Default in \code{NULL}
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
#' @param rotate type of rotation to perform after the initial orthogonal
#' parameters have been extracted. See below for list of possible rotations
#' @param startvalues user declared start values for parameters
#' @param quadpts number of quadrature points per dimension
#' @param ncycles the number of EM iterations to be performed
#' @param tol if the largest change in the EM cycle is less than this value
#' then the EM iteration are stopped early
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#' option. Only prints patterns that have standardized residuals greater than 
#' \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param x an object of class \code{mirt} to be plotted or printed
#' @param object a model estimated from \code{mirt} of class \code{mirt}
#' @param object2 a second model estimated from \code{mirt} of class
#' \code{mirt} with more estimated parameters than \code{object}
#' @param suppress a numeric value indicating which (possibly rotated) factor
#' loadings should be suppressed. Typical values are around .3 in most
#' statistical software. Default is 0 for no suppression
#' @param digits number of significant digits to be rounded
#' @param type type of plot to view; can be \code{'curve'} for the total test
#' score as a function of two dimensions, or \code{'info'} to show the test
#' information function for two dimensions
#' @param npts number of quadrature points to be used for plotting features.
#' Larger values make plots look smoother
#' @param rot allows rotation of the 3D graphics
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param nowarn logical; suppress warnings from dependent packages?
#' @param debug logical; turn on debugging features?
#' @param ... additional arguments to be passed
#' @section Convergence:
#' 
#' Unrestricted full-information factor analysis is known to have problems with
#' convergence, and some items may need to be constrained or removed entirely
#' to allow for an acceptable solution. Be mindful of the item facility values
#' that are printed with \code{coef} since these will be helpful in determining
#' whether a guessing parameter should be removed (item facility value is too
#' close to the guessing parameter) or if an item should be constrained or
#' removed entirely (values too close to 0 or 1). As a general rule, items with
#' facilities greater than .95, or items that are only .05 greater than the
#' guessing parameter, should be considered for removal from the analysis or
#' treated with prior distributions. Also, increasing the number of quadrature
#' points per dimension may help to stabilize the estimation process.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{polymirt}},
#' \code{\link{confmirt}}, \code{\link{polymirt}}, \code{\link{itemplot}}
#' @references
#' 
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of
#' item parameters: Application of an EM algorithm. \emph{Psychometrika,
#' 46}(4), 443-459.
#' 
#' Bock, R. D., Gibbons, R., & Muraki, E. (1988). Full-Information Item Factor
#' Analysis. \emph{Applied Psychological Measurement, 12}(3), 261-280.
#' 
#' Carroll, J. B. (1945). The effect of difficulty and chance success on
#' correlations between items and between tests. \emph{Psychometrika, 26},
#' 347-372.
#' 
#' Ramsay, J. O. (1975). Solving implicit equations in psychometric data
#' analysis. \emph{Psychometrika, 40}(3), 337-360.
#' 
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., &
#' Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics,
#' and Full-information Item Factor Analysis [Computer software]. Lincolnwood,
#' IL: Scientific Software International.
#' @keywords models
#' @usage 
#' mirt(fulldata, nfact, guess = 0, SE = FALSE, prev.cor = NULL, par.prior = FALSE,
#'   startvalues = NULL, quadpts = NULL, ncycles = 300, tol = .001, nowarn = TRUE, 
#'   debug = FALSE, ...)
#' 
#' \S4method{summary}{mirt}(object, rotate='varimax', suppress = 0, digits = 3, ...)
#' 
#' \S4method{coef}{mirt}(object, digits = 3, ...)
#' 
#' \S4method{anova}{mirt}(object, object2, ...)
#' 
#' \S4method{fitted}{mirt}(object, digits = 3, ...)
#' 
#' \S4method{plot}{mirt}(x, type = 'info', npts = 50, rot = list(x = -70, y = 30, z = 10), ...)
#' 
#' \S4method{residuals}{mirt}(object, restype = 'LD', digits=3, printvalue = NULL, ...)
#' @export mirt
#' @examples
#' 
#' \dontrun{
#' #load LSAT section 7 data and compute 1 and 2 factor models
#' data(LSAT7)
#' fulldata <- expand.table(LSAT7)
#' 
#' (mod1 <- mirt(fulldata, 1))
#' summary(mod1)
#' residuals(mod1)
#' plot(mod1) #test information function
#' 
#' (mod2 <- mirt(fulldata, 2))
#' summary(mod2)
#' coef(mod2)
#' residuals(mod2)
#' 
#' anova(mod1, mod2) #compare the two models
#' scores <- fscores(mod2) #save factor score table
#' 
#' ###########
#' data(SAT12)
#' fulldata <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' 
#' #without guessing scree(tmat) #looks like a 2 factor solution
#' mod1 <- mirt(fulldata, 1)
#' mod2 <- mirt(fulldata, 2)
#' mod3 <- mirt(fulldata, 3)
#' anova(mod1,mod2)
#' anova(mod2, mod3) #negative AIC, 2 factors probably best
#' 
#' #with guessing
#' mod1g <- mirt(fulldata, 1, guess = .1)
#' coef(mod1g)
#' mod2g <- mirt(fulldata, 2, guess = .1)
#' coef(mod2g)
#' anova(mod1g, mod2g)
#' summary(mod2g, rotate='promax')
#'      }
#' 
mirt <- function(fulldata, nfact, guess = 0, SE = FALSE, prev.cor = NULL, par.prior = FALSE, 
	startvalues = NULL, quadpts = NULL, ncycles = 300, tol = .001, nowarn = TRUE, 
	debug = FALSE, ...)
{ 
	fn <- function(pars, r1, N, guess, Theta, prior, parprior){
		a <- pars[1:(length(pars)-1)]
		d <- pars[length(pars)]		
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
  
	Call <- match.call()    
	itemnames <- colnames(fulldata)
	fulldata <- as.matrix(fulldata)	
	fulldata.original <- fulldata 
	fulldata[is.na(fulldata)] <- 9	
	if (!any(fulldata.original %in% c(0,1,NA))) stop("Data must contain only 0, 1, or NA.")	
	nitems <- ncol(fulldata)  
	colnames(fulldata) <- itemnames
	if (length(guess) == 1) guess <- rep(guess,nitems)
		else if (length(guess) > nitems || length(guess) < nitems) 
			stop("The number of guessing parameters is incorrect.")	
	pats <- apply(fulldata,1,paste,collapse = "/")
	freqs <- table(pats)
	nfreqs <- length(freqs)
	K <- rep(2,nitems)
	r <- as.vector(freqs)
	N <- nrow(fulldata) 
	tabdata <- unlist(strsplit(cbind(names(freqs)),"/"))
	tabdata <- matrix(as.numeric(tabdata),nfreqs,nitems,TRUE)
	tabdata <- cbind(tabdata,r)    
	if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))  
	theta <- as.matrix(seq(-4,4,length.out = quadpts))
	if(quadpts^nfact <= 10000){
		Theta <- thetaComb(theta,nfact)
		prior <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
		prior <- prior/sum(prior)
	} else stop('Greater than 10000 quadrature points, reduce number.')
	facility <- colMeans(na.omit(fulldata.original))
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
	if (any(class(prev.cor) == c('mirt','bmirt'))) Rpoly <- prev.cor$cormat
		else if(!is.null(prev.cor)) {
			if (ncol(prev.cor) == nrow(prev.cor)) Rpoly <- prev.cor
				else stop("Correlation matrix is not square.\n")
	} else 
		Rpoly <- cormod(na.omit(fulldata.original),K,guess)   
	if (is.null(startvalues)){ 
		suppressMessages(pars <- start.values(na.omit(fulldata.original),guess,Rpoly,
			nfact=nfact,nowarn=nowarn))
		pars[pars > 3] <- 3
		pars[pars < -3] <- -3	
	} else {
		if ((ncol(startvalues) != (nfact + 1)) || (nrow(startvalues) != nitems))
			stop("Startvalues are declared incorrectly.")  
		pars <- startvalues  
	} 
	diag(Rpoly) <- 1
	item <- 1
	lastpars2 <- lastpars1 <- rate <- matrix(0,nrow=nitems,ncol=ncol(pars))    
	startvalues <- pars
	converge <- 1  
	problemitems <- c()
	index <- 1:nitems  
	if(debug) print(startvalues)      

	# EM loop
	for (cycles in 1:ncycles)
	{       
		rlist <- Estep.mirt(pars,tabdata,Theta,prior,guess)		
		if (debug) print(sum(r*log(rlist[[3]])))
		lastpars2 <- lastpars1
		lastpars1 <- pars	
		for(i in 1:nitems){
			if(guess[i] == 0)	
				maxim <- try(optim(pars[i, ],fn=fn,gr=gr,r1=rlist[[1]][i, ],N=rlist[[2]][i, ],
				guess=guess[i],Theta=Theta,prior=prior,parprior=par.prior[i, ],method="BFGS"))
			else 
				maxim <- try(optim(pars[i, ],fn=fn,r1=rlist[[1]][i, ],N=rlist[[2]][i, ],
				guess=guess[i],Theta=Theta,prior=prior,parprior=par.prior[i, ],method="BFGS"))
			if(class(maxim) == "try-error"){
				problemitems <- c(problemitems, i)
				converge <- 0
				next
			}		  
			pars[i, ] <- maxim$par	  
		}	
		if(!suppressAutoPrior){
			if(any(abs(pars[ ,nfact+1]) > 4)){
				ints <- index[abs(pars[ ,nfact+1]) > 4] 	
				par.prior[ints,3] <- 2
				if(any(abs(pars[ ,nfact+1]) > 5.5)){
					ints <- index[abs(pars[ ,nfact+1]) > 5.5] 	
					par.prior[ints,3] <- 1
				} 
			}
			if(nfact > 1){ 
				norm <- sqrt(1 + rowSums(pars[ ,1:nfact]^2))
				alp <- as.matrix(pars[ ,1:nfact]/norm)
				FF <- alp %*% t(alp)
				V <- eigen(FF)$vector[ ,1:nfact]
				L <- eigen(FF)$values[1:nfact]      
				F <- V %*% sqrt(diag(L))
				h2 <- rowSums(F^2)
				if(any(h2 > .95)){
					ind <- index[h2 > .95]
					par.prior[ind,1] <- 1.2
					if(any(h2 > .98)){
						ind <- index[h2 > .98]
						par.prior[ind,1] <- 1.5
					} 
				}
			}
		}  
		maxdif <- max(abs(lastpars1 - pars))	
		if (maxdif < tol) break    
		# rate acceleration adjusted every third cycle
		if (cycles %% 3 == 0 & cycles > 6) 
		{
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
		Model probably did not converged.")  
	if(length(problemitems) > 0) warning("Problem with the M-step for item(s): ", 
		paste(unique(problemitems), " "))	
	lastchange <- lastpars1 - pars
	if (cycles == ncycles){
		converge <- 0  
		message("Estimation terminated after ", cycles, " EM loops. Maximum changes:") 
		message("\n slopes = ", round(max(abs(lastchange[ ,1:nfact])),4), ", intercepts = ", 
			round(max(abs(lastchange[ ,ncol(pars)])),4) ,"\n", sep="")
	}	    
	prior <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
	prior <- prior/sum(prior)  
	rlist <- Estep.mirt(pars,tabdata,Theta,prior,guess)      	  
	Pl <- rlist[[3]]  
	logLik <- sum(r*log(Pl))
	vcovpar <- matrix(999)
	if(SE){
		LLfun <- function(pars,tabdata,Theta,prior,guess){
			nfact <- ncol(Theta)
			pars2 <- matrix(pars, ncol=nfact+1)
			rlist <- Estep.mirt(pars2,tabdata,Theta,prior,guess)      	  
			Pl <- rlist[[3]]  
			logLik <- sum(r*log(Pl))
			-1*logLik		
		}
		fmin <- nlm(LLfun, as.numeric(pars), tabdata=tabdata,Theta=Theta,prior=prior,
			guess=guess, hessian=TRUE, gradtol=1)
		vcovpar <- solve(fmin$hessian)		
	}	
	logN <- 0
	logr <- rep(0,length(r))
	for (i in 1:N) logN <- logN + log(i)
	for (i in 1:length(r)) 
		for (j in 1:r[i]) 
			logr[i] <- logr[i] + log(j)    
	df <- (length(r) - 1) - nitems*(nfact + 1) + nfact*(nfact - 1)/2 
	X2 <- 2 * sum(r * log(r/(N*Pl)))	
	logLik <- logLik + logN/sum(logr)	
	p <- 1 - pchisq(X2,df)  
	AIC <- (-2) * logLik + 2 * length(pars)
	BIC <- (-2) * logLik + length(pars)*log(N)
	if(any(is.na(fulldata.original))) p <- 2	

	# pars to FA loadings
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

	mod <- new('mirtClass', EMiter=cycles, pars=pars, guess=guess, X2=X2, df=df, p=p,
		AIC=AIC, BIC=BIC, logLik=logLik, F=F, h2=h2, tabdata=tabdata, Theta=Theta, Pl=Pl, 
		fulldata=fulldata.original, cormat=Rpoly, facility=facility, converge=converge, 
		quadpts=quadpts, vcov=vcovpar, Call=Call)	  
	return(mod)    
}

#Methods 

setMethod(
	f = "print",
	signature = signature(x = 'mirtClass'),
	definition = function(x, ...){  
		cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information factor analysis with ", ncol(x@F), " factor",
			if(ncol(x@F)>1) "s", "\n", sep="")
			if(x@converge == 1)	
				cat("Converged in ", x@EMiter, " iterations using ", x@quadpts,
				" quadrature points.\n", sep="")
		else 	
			cat("Estimation stopped after ", x@EMiter, " iterations using ", 
				x@quadpts, " quadrature points.\n", sep="")
		cat("Log-likelihood =", x@logLik, "\n")
		cat("AIC =", x@AIC, "\n")		
		cat("BIC =", x@BIC, "\n")
		if(x@p < 1)
			cat("G^2 = ", round(x@X2,2), ", df = ", 
				x@df, ", p = ", round(x@p,4), "\n", sep="")
		else 
			cat("G^2 = ", NA, ", df = ", 
				x@df, ", p = ", NA, "\n", sep="")		
	}
)

setMethod(
	f = "show",
	signature = signature(object = 'mirtClass'),
	definition = function(object){  
		cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		cat("Full-information factor analysis with ", ncol(object@F), " factor",
			if(ncol(object@F)>1) "s", "\n", sep="")
			if(object@converge == 1)	
				cat("Converged in ", object@EMiter, " iterations using ", object@quadpts,
				" quadrature points.\n", sep="")
		else 	
			cat("Estimation stopped after ", object@EMiter, " iterations using ", 
				object@quadpts,	" quadrature points.\n", sep="")
		cat("Log-likelihood =", object@logLik, "\n")
		cat("AIC =", object@AIC, "\n")		
		cat("BIC =", object@BIC, "\n")
		if(object@p < 1)
			cat("G^2 = ", round(object@X2,2), ", df = ", 
				object@df, ", p = ", round(object@p,4), "\n", sep="")
		else 
			cat("G^2 = ", NA, ", df = ", 
				object@df, ", p = ", NA, "\n", sep="")			
	}
)

setMethod(
	f = "summary",
	signature = 'mirtClass',
	definition = function(object, rotate = 'varimax', suppress = 0, digits = 3, ...){
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
			rownames(loads) <- rownames(object@pars)
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
			rownames(loads) <- rownames(object@pars)			
			cat("\nRotated factor loadings: \n\n")
			print(loads,digits)		
			if(attr(rotF, "oblique")){
				cat("\nFactor correlations: \n\n")
				Phi <- rotF$Phi	  
				Phi <- round(Phi, digits)
				colnames(Phi) <- rownames(Phi) <- colnames(F)
				print(Phi)            
			}	
			cat("\nRotated SS loadings: ",round(SS,digits), "\n")		
			if(any(h2 > 1)) 
				warning("Solution has heywood cases. Interpret with caution.") 
			invisible(list(rotF$loadings,h2))  
		}  
	}
)

setMethod(
	f = "coef",
	signature = 'mirtClass',
	definition = function(object, SE = TRUE, digits = 3, ...){  
		a <- as.matrix(object@pars[ ,1:(ncol(object@pars)-1)])
		d <- object@pars[ ,ncol(object@pars)]
		A <- sqrt(apply(a^2,1,sum))
		B <- -d/A  
		if (ncol(a) > 1){  
			parameters <- cbind(object@pars,object@guess,object@facility,A,B)    
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),"d","guess", 
			"facility","mvdisc","mvint")	  
			cat("\nUnrotated parameters, multivariate discrimination and intercept: \n\n")
			print(round(parameters, digits))  	
		} else {
			parameters <- cbind(object@pars,object@guess,object@facility) 
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),"d","guess","facility")    
			cat("\nParameter slopes and intercepts: \n\n")	
			print(round(parameters, digits))	  
		}
		ret <- list(parameters)
		if(ncol(object@vcov) != 1){
			cat("\nStd. Errors: \n\n")	
			SEs <- matrix(sqrt(diag(object@vcov)), ncol = ncol(a) + 1)
			colnames(SEs) <- colnames(parameters)[1:(ncol(a) + 1)]
			rownames(SEs) <- rownames(parameters)
			print(SEs, digits)
			ret <- list(parameters,SEs)
		}
		invisible(ret)
	}
)

setMethod(
	f = "anova",
	signature = signature(object = 'mirtClass'),
	definition = function(object, object2, ...){
		dots <- list(...)		
		df <- object@df - object2@df  
		if(df < 0){
			temp <- object
			object <- object2
			object2 <- temp
		}
		X2 <- 2*object2@logLik - 2*object@logLik 		
		AICdiff <- object@AIC - object2@AIC    
		BICdiff <- object@BIC - object2@BIC
		cat("\nChi-squared difference: \n\nX2 = ", round(X2,3), ", df = ",
			df, ", p = ", round(1 - pchisq(X2,abs(df)),4), "\n", sep="")
		cat("AIC difference = ", round(AICdiff,3), "\n")  
		cat("BIC difference = ", round(BICdiff,3), "\n")
	}
)

setMethod(
	f = "residuals",
	signature = signature(object = 'mirtClass'),
	definition = function(object, restype = 'LD', digits = 3, printvalue = NULL, ...){   	
		Theta <- object@Theta
		fulldata <- object@fulldata	
		N <- nrow(fulldata)	
		J <- ncol(fulldata)
		nfact <- ncol(object@F)
		lambdas <- matrix(object@pars[,1:nfact], J)
		zetas <- object@pars[,(nfact+1)]
		guess <- object@guess
		guess[is.na(guess)] <- 0		
		if(restype == 'LD'){
			res <- matrix(0,J,J)
			diag(res) <- NA
			colnames(res) <- rownames(res) <- colnames(fulldata)
			prior <- dmvnorm(Theta,rep(0,nfact),diag(nfact))
			prior <- prior/sum(prior)
			for(i in 1:J){			
				for(j in 1:J){
					if(i < j){
						P1 <- P.mirt(lambdas[i,],zetas[i], Theta, guess[i])
						P2 <- P.mirt(lambdas[j,],zetas[j], Theta, guess[j])
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
			expected <- round(N * object@Pl/sum(object@Pl),digits)  
			tabdata <- object@tabdata
			freq <- tabdata[ ,ncol(tabdata)]
			tabdata[tabdata[ ,1:ncol(object@fulldata)] == 9] <- NA
			tabdata[ ,ncol(tabdata)] <- freq
			tabdata <- cbind(tabdata,expected,res)
			colnames(tabdata) <- c(colnames(fulldata), "freq", "exp", "std_res")
			if(!is.null(printvalue)){
				if(!is.numeric(printvalue)) stop('printvalue is not a number.')
				tabdata <- tabdata[abs(tabdata[ ,ncol(tabdata)]) > printvalue, ]
			}			
			return(tabdata)				
		}					
	}
)

setMethod(
	f = "plot",
	signature = signature(x = 'mirtClass', y = 'missing'),
	definition = function(x, y, type = 'info', npts = 50, 
		rot = list(xaxis = -70, yaxis = 30, zaxis = 10))
	{  
		if (!type %in% c('curve','info','contour','infocontour')) stop(type, " is not a valid plot type.")
		rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
		a <- as.matrix(x@pars[ ,1:(ncol(x@pars) - 1)])
		d <- x@pars[ ,ncol(x@pars)]
		g <- x@guess
		A <- as.matrix(sqrt(apply(a^2,1,sum)))
		B <- -d/A
		if(ncol(a) > 2 ) stop("Can't plot high dimentional solutions.\n")
		theta <- seq(-4,4,length.out=npts)
		Theta <- thetaComb(theta, ncol(a))
		P <- matrix(0, ncol=length(g), nrow = nrow(as.matrix(Theta)))
		for(i in 1:nrow(a)) P[ ,i] <- P.mirt(a[i, ],d[i],as.matrix(Theta),g[i])  
		Ptot <- rowSums(P) 
		if(ncol(a) == 2) {			
			I <- (P * (1 - P)) %*% A^2	
			if(type == 'curve'){
				plt <- data.frame(cbind(Ptot,Theta))
				colnames(plt) <- c("Ptot", "Theta1", "Theta2")
				return(wireframe(Ptot ~ Theta1 + Theta2, data = plt, main = "Test Score Surface", 
					zlab = "Test \nScore", xlab = "Theta 1", ylab = "Theta 2", 
					scales = list(arrows = FALSE), screen = rot))
			}	
			if(type == 'contour'){
				plt <- data.frame(cbind(Ptot,Theta))
				colnames(plt) <- c('Ptot','Theta1','Theta2')
				contour(theta, theta, matrix(Ptot,length(theta),length(theta)), 
					main = "Test Scores Contour", xlab = "Theta 1", ylab = "Theta 2")			
			}			
			if(type == 'infocontour'){
				plt <- data.frame(cbind(I,Theta))
				colnames(plt) <- c('I','Theta1','Theta2')
				contour(theta, theta, matrix(I,length(theta),length(theta)), 
					main = "Test Information Contour", xlab = "Theta 1", ylab = "Theta 2")		
			}				
			if(type == 'info'){				 
				plt <- data.frame(cbind(I,Theta))
				colnames(plt) <- c("I", "Theta1", "Theta2")
				return(wireframe(I ~ Theta1 + Theta2, data = plt, main = "Test Information", 
					zlab = "I", xlab = "Theta 1", ylab = "Theta 2", 
					scales = list(arrows = FALSE), screen = rot))
			}
		} else {
			if(type == 'curve')  
				plot(Theta, Ptot, type='l', main = 'Test score plot', xlab = 'Theta', ylab='Test Score')
			if(type == 'info'){
				I <- (P * (1 - P)) %*% a^2 
				plot(Theta, I, type='l', main = 'Test Information', xlab = 'Theta', ylab='Information')			
			}	
			if(type == 'contour' || type == 'infocontour') 
				cat('No \'contour\' plots for 1-dimensional models\n')					
		} 			
	}		
)	

setMethod(
	f = "fitted",
	signature = signature(object = 'mirtClass'),
	definition = function(object, digits = 3, ...){  
		expected <- round(nrow(object@fulldata) * object@Pl,digits)  
		tabdata <- object@tabdata
		freq <- tabdata[ ,ncol(tabdata)]
		tabdata[tabdata[ ,1:ncol(object@fulldata)] == 9] <- NA
		tabdata[ ,ncol(tabdata)] <- freq
		tabdata <- cbind(tabdata,expected)
		colnames(tabdata) <- c(colnames(object@fulldata),"freq","exp")	
		print(tabdata)
		invisible(tabdata)
	}
)
    



