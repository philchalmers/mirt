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
	representation = representation(EMiter='numeric', pars='list', guess='numeric', 
		K='numeric', parsSE='list', X2='numeric', df='numeric', p='numeric', AIC='numeric', logLik='numeric',
		F='matrix', h2='numeric', tabdata='matrix', tabdatalong='matrix', Theta='matrix', Pl='numeric',
		data='matrix', cormat='matrix', facility='numeric', converge='numeric', itemloc = 'numeric',
		quadpts='numeric', BIC='numeric', vcov='matrix', RMSEA='numeric', Call='call'),	
	validity = function(object) return(TRUE)
)	

#' Full-Information Item Factor Analysis (Multidimensional Item Response
#' Theory)
#' 
#' \code{mirt} fits an unconditional maximum likelihood factor analysis model
#' to dichotomous data under the item response theory paradigm. Pseudo-guessing
#' parameters may be included but must be declared as constant.
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
#' @param data a \code{matrix} or \code{data.frame} that consists of only
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
#' @param object a model estimated from \code{mirt} of class \code{mirtClass}
#' @param object2 a second model estimated from \code{mirt} of class
#' \code{mirtClass} with more estimated parameters than \code{object}
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
#' whether a guessing parameter is causing problems (item facility value is too
#' close to the guessing parameter) or if an item should be constrained or
#' removed entirely (values too close to 0 or 1). As a general rule, items with
#' facilities greater than .95, or items that are only .05 greater than the
#' guessing parameter, should be considered for removal from the analysis or
#' treated with prior distributions. Also, increasing the number of quadrature
#' points per dimension may help to stabilize the estimation process.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{polymirt}},
#' \code{\link{confmirt}}, \code{\link{bfactor}}, \code{\link{itemplot}}
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
#' mirt(data, nfact, guess = 0, SE = FALSE, prev.cor = NULL, par.prior = FALSE,
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
#' data <- expand.table(LSAT7)
#' 
#' (mod1 <- mirt(data, 1))
#' summary(mod1)
#' residuals(mod1)
#' plot(mod1) #test information function
#' 
#' (mod2 <- mirt(data, 2))
#' summary(mod2)
#' coef(mod2)
#' residuals(mod2)
#' plot(mod2)
#' 
#' anova(mod1, mod2) #compare the two models
#' scores <- fscores(mod2) #save factor score table
#' 
#' ###########
#' data(SAT12)
#' data <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' 
#' mod1 <- mirt(data, 1)
#' mod2 <- mirt(data, 2)
#' mod3 <- mirt(data, 3)
#' anova(mod1,mod2)
#' anova(mod2, mod3) #negative AIC, 2 factors probably best
#' 
#' #with guessing
#' mod1g <- mirt(data, 1, guess = .1)
#' coef(mod1g)
#' mod2g <- mirt(data, 2, guess = .1)
#' coef(mod2g)
#' anova(mod1g, mod2g)
#' summary(mod2g, rotate='promax')
#'      }
#' 
mirt <- function(data, nfact, guess = 0, SE = FALSE, prev.cor = NULL, par.prior = FALSE, 
	startvalues = NULL, quadpts = NULL, ncycles = 300, tol = .001, nowarn = TRUE, 
	debug = FALSE, ...)
{ 
	fn <- function(par, rs, gues, Theta, prior, parprior){
		nzeta <- ncol(rs) - 1
		a <- par[1:(length(par)-nzeta)]
		d <- par[(length(a)+1):length(par)]				
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
	} else Rpoly <- cormod(na.omit(data),K,guess)
	FA <- psych::fa(Rpoly,nfact,rotate = 'none', warnings= FALSE, fm="minres")	
	loads <- unclass(loadings(FA))
	u <- FA$unique
	u[u < .1 ] <- .25	
	cs <- sqrt(u)
	lambdas <- loads/cs		
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
	npars <- length(unlist(pars))
	if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))  
	theta <- as.matrix(seq(-4,4,length.out = quadpts))
	if(quadpts^nfact <= 10000){
		Theta <- thetaComb(theta,nfact)
		prior <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact))
		prior <- prior/sum(prior)
	} else stop('Greater than 10000 quadrature points, reduce number.')	  
	lastpars2 <- lastpars1 <- pars	    
	startvalues <- pars
	converge <- 1  
	problemitems <- c()
	index <- 1:J  
	if(debug) print(startvalues)
	
	# EM loop 
	for (cycles in 1:ncycles)
	{       
		rlist <- Estep.mirt(pars, tabdata, Theta, prior, guess, itemloc)						
		lastpars2 <- lastpars1
		lastpars1 <- pars				
		for(i in 1:J){
			par <- c(pars$lambdas[i, ], pars$zetas[[i]])
			itemsel <- c(itemloc[i]:(itemloc[i+1] - 1))					
			maxim <- try(optim(par, fn=fn, rs=rlist$r1[, itemsel], gues=guess[i], Theta=Theta, prior=prior, 
				parprior=par.prior[i, ], control=list(maxit=25)))
			if(class(maxim) == "try-error"){
				problemitems <- c(problemitems, i)
				converge <- 0
				next
			}		  
			pars$lambdas[i, ] <- maxim$par[1:nfact]
			pars$zetas[[i]] <- maxim$par[(nfact+1):length(par)]	  
		}				
		maxdif <- max(abs(unlist(lastpars1) - unlist(pars)))	
		if (maxdif < tol && cycles > 5) break    
		# rate acceleration adjusted every third cycle
		if (cycles %% 3 == 0 & cycles > 6)		 
			pars <- rateChange(pars, lastpars1, lastpars2)			     	
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
	lastchange <- unlist(lastpars1) - unlist(pars)
	if (cycles == ncycles){
		converge <- 0  
		message("Estimation terminated after ", cycles, " EM loops. Maximum changes:") 
		message("\n slopes = ", round(max(abs(lastchange[ ,1:nfact])),4), ", intercepts = ", 
			round(max(abs(lastchange[ ,ncol(pars)])),4) ,"\n", sep="")
	}	    	 
	rlist <- Estep.mirt(pars, tabdata, Theta, prior, guess, itemloc)      	  
	Pl <- rlist$expected  
	logLik <- sum(r*log(Pl))
	vcovpar <- matrix(999)
	parsSE <- list()
	if(SE && nfact == 1){
		LLfun <- function(p, pars, tabdata, Theta, prior, guess, itemloc){
			pars2 <- rebuildPars(p, pars)		
			rlist <- Estep.mirt(pars2, tabdata, Theta, prior, guess, itemloc)     	  
			Pl <- rlist$expected
			logLik <- sum(r*log(Pl))
			-1*logLik		
		}
		fmin <- nlm(LLfun, unlist(pars), pars=pars,tabdata=tabdata,Theta=Theta,prior=prior,
			guess=guess, itemloc=itemloc, hessian=TRUE, gradtol=.1)
		vcovpar <- solve(fmin$hessian)
		parsSE <- rebuildPars(sqrt(diag(vcovpar)), pars)	
	}	
	logN <- 0
	logr <- rep(0,length(r))	
	for (i in 1:N) logN <- logN + log(i)
	for (i in 1:length(r)) 
		for (j in 1:r[i]) 
			logr[i] <- logr[i] + log(j)    	
	df <- (length(r) - 1) - npars + nfact*(nfact - 1)/2 
	X2 <- 2 * sum(r * log(r/(N*Pl)))	
	logLik <- logLik + logN/sum(logr)	
	p <- 1 - pchisq(X2,df)  
	AIC <- (-2) * logLik + 2 * npars
	BIC <- (-2) * logLik + npars*log(N)
	RMSEA <- ifelse((X2 - df) > 0, 
	    sqrt(X2 - df) / sqrt(df * (N-1)), 0)
	if(any(is.na(data.original))) p <- 2
	guess[K > 2] <- NA	

	# pars to FA loadings
	if (nfact > 1) norm <- sqrt(1 + rowSums(pars$lambdas[ ,1:nfact]^2))
		else norm <- as.matrix(sqrt(1 + pars$lambdas[ ,1]^2))  
	alp <- as.matrix(pars$lambdas[ ,1:nfact]/norm)
	FF <- alp %*% t(alp)
	V <- eigen(FF)$vector[ ,1:nfact]
	L <- eigen(FF)$values[1:nfact]
	if (nfact == 1) F <- as.matrix(V * sqrt(L))
		else F <- V %*% sqrt(diag(L))  
	if (sum(F[ ,1] < 0)) F <- (-1) * F 
	colnames(F) <- paste("F_", 1:ncol(F),sep="")	
	h2 <- rowSums(F^2) 

	mod <- new('mirtClass', EMiter=cycles, pars=pars, guess=guess, parsSE=parsSE, X2=X2, df=df, 
		p=p, itemloc=itemloc, AIC=AIC, BIC=BIC, logLik=logLik, F=F, h2=h2, tabdata=tabdata2, 
		Theta=Theta, Pl=Pl, data=data.original, cormat=Rpoly, facility=facility, converge=converge, 
		quadpts=quadpts, vcov=vcovpar, RMSEA=RMSEA, K=K, tabdatalong=tabdata, Call=Call)	  
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
				x@df, ", p = ", round(x@p,4),", RMSEA = ", round(x@RMSEA,3), "\n", sep="")
		else 
			cat("G^2 = ", NA, ", df = ", 
				x@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="" )		
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
				object@df, ", p = ", round(object@p,4),", RMSEA = ", round(object@RMSEA,3),
                "\n", sep="")
		else 
			cat("G^2 = ", NA, ", df = ", 
				object@df, ", p = ", NA, "RMSEA = ", NA, "\n", sep="")			
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
			names(SS) <- colnames(F)
			cat("\nUnrotated factor loadings: \n\n")
			loads <- round(cbind(F,h2),digits)
			rownames(loads) <- rownames(object@pars$lambdas)
			print(loads)	    	 
			cat("\nSS loadings: ",round(SS,digits), "\n")
			cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
			invisible(list(F,h2))
		} else {	
			F <- object@F
			h2 <- as.matrix(object@h2)						
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
			parameters <- cbind(a,d,object@guess,A,B)    
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),paste("d_",1:max(K-1),sep=""),"guess", 
				"mvdisc",paste("mvint_",1:max(K-1),sep=""))	  
			cat("\nUnrotated parameters, multivariate discrimination and intercept: \n\n")
			print(round(parameters, digits))  	
		} else {
			parameters <- cbind(a,d,object@guess)
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),paste("d_",1:max(K-1),sep=""),"guess")   
			cat("\nParameter slopes and intercepts: \n\n")	
			print(round(parameters, digits))	  
		}
		ret <- list(parameters)
		if(length(object@parsSE) > 1){
			if(SE){
				cat("\nStd. Errors: \n\n")	
				SEs <- matrix(sqrt(diag(object@vcov)), ncol = ncol(a) + 1)
				colnames(SEs) <- colnames(parameters)[1:(ncol(a) + 1)]
				rownames(SEs) <- rownames(parameters)
				print(SEs, digits)
				ret <- list(parameters,SEs)
			}
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
		K <- object@K
		Theta <- object@Theta
		data <- object@data	
		N <- nrow(data)	
		J <- ncol(data)
		nfact <- ncol(object@F)
		lambdas <- object@pars$lambdas
		zetas <- object@pars$zetas
		guess <- object@guess
		guess[is.na(guess)] <- 0
		itemloc <- object@itemloc
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
			r <- object@tabdata[ ,ncol(object@tabdata)]
			res <- round((r - object@Pl * nrow(object@data)) / 
				sqrt(object@Pl * nrow(object@data)),digits)
			expected <- round(N * object@Pl/sum(object@Pl),digits)  
			tabdata <- object@tabdata
			freq <- tabdata[ ,ncol(tabdata)]			
			tabdata[tabdata[ ,1:ncol(object@data)] == 99] <- NA
			tabdata[ ,ncol(tabdata)] <- freq
			tabdata <- cbind(tabdata,expected,res)
			colnames(tabdata) <- c(colnames(object@tabdata),"freq","exp")	
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
				P <- P.poly(a[j,], d[[j]], Theta, itemexp = FALSE)		
				for(i in 1:K[j]){
					w1 <- P[,i]*(1-P[,i])*A[j]
					w2 <- P[,i+1]*(1-P[,i+1])*A[j]
					I <- ((w1 - w2)^2) / (P[,i] - P[,i+1]) * P[,i]
					info <- info + I
				}
			} else {
				P <- P.mirt(a[j,], d[[j]], Theta, guess[j])
				Pstar <- P.mirt(a[j,], d[[j]], Theta, 0)
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
				return(lattice::wireframe(info ~ Theta1 + Theta2, data = plt, main = "Test Information", 
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
	f = "fitted",
	signature = signature(object = 'mirtClass'),
	definition = function(object, digits = 3, ...){  
		expected <- round(nrow(object@data) * object@Pl,digits)  
		tabdata <- object@tabdata
		freq <- tabdata[ ,ncol(tabdata)]
		tabdata[tabdata[ ,1:ncol(object@data)] == 9] <- NA
		tabdata[ ,ncol(tabdata)] <- freq
		tabdata <- cbind(tabdata,expected)
		colnames(tabdata) <- c(colnames(object@tabdata),"freq","exp")	
		print(tabdata)
		invisible(tabdata)
	}
)
    



