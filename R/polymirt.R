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
	representation = representation(pars = 'list', guess = 'numeric', SEg = 'numeric',
        upper = 'numeric', SEup = 'numeric',
		SEpars = 'matrix', cycles = 'numeric', Theta = 'matrix', fulldata = 'matrix', data = 'matrix', 
		K = 'numeric', F = 'matrix', h2 = 'numeric', itemloc = 'numeric', AIC = 'numeric',
		converge = 'numeric', logLik = 'numeric', SElogLik = 'numeric', df = 'integer', 
		G2 = 'numeric', p = 'numeric', tabdata = 'matrix', BIC = 'numeric', estGuess = 'logical', 
		RMSEA = 'numeric', rotate='character', null.mod='S4', TLI = 'numeric', Call = 'call'),	
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
#' @param upper initial (or fixed) upper bound parameters for 4-PL model. Can be 
#' entered as a single value to assign a global upper bound parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param estGuess a logical vector indicating which lower-asymptote parameters
#' to be estimated (default is null, and therefore is contingent on the values
#' in \code{guess}). By default, if any value in \code{guess} is greater than 0
#' then its respective \code{estGuess} value is set to \code{TRUE}.
#' Additionally, beta priors are automatically imposed for estimated parameters
#' that correspond to the input guessing value.
#' @param estUpper same function as \code{estGuess}, but for upper bound parameters
#' @param prev.cor use a previously computed correlation matrix to be used to
#' estimate starting values the estimation. The input could be any correlation
#' matrix, but it is advised to use a matrix of polychoric correlations.
#' @param rotate type of rotation to perform after the initial orthogonal
#' parameters have been extracted by using \code{summary}; default is \code{'varimax'}. 
#' See \code{\link{mirt}} for list of possible rotations. If \code{rotate != ''} in the 
#' \code{summary} input then the default from the object is ignored and the new rotation 
#' from the list is used instead
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
#' @param verbose logical; display iteration history during estimation?
#' @param calcLL logical; calculate the log-likelihood?
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param draws the number of Monte Carlo draws to estimate the log-likelihood
#' @param type either \code{'info'} or \code{'infocontour'} to plot test
#' information plots
#' @param debug logical; turn on debugging features?
#' @param technical list specifying subtle parameters that can be adjusted. Can be
#' \describe{
#' \item{NCYCLES}{max number of MH-RM cycles; default 2000}
#' \item{BURNIN}{number of burn in cycles (stage 1); default 150}
#' \item{SEMCYCLES}{number of SEM cycles (stage 2); default 50}
#' \item{KDRAWS}{number of paralell MH sets to be drawn; default 1}
#' \item{TOL}{minimum threshold tolerance for convergence of MH-RM, must occur on three consecutive
#' occations; default .001}
#' }
#' @param ... additional arguments to be passed to the \code{\link{confmirt}} estimation engine
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
#' polymirt(data, nfact, guess = 0, upper = 1, estGuess = NULL, estUpper = NULL,
#' prev.cor = NULL, rotate = 'varimax', verbose = TRUE, calcLL = TRUE, draws = 2000, debug = FALSE, 
#' technical = list(), ...)
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
polymirt <- function(data, nfact, guess = 0, upper = 1, estGuess = NULL, estUpper = NULL,
    prev.cor = NULL, rotate = 'varimax', verbose = TRUE, calcLL = TRUE, draws = 2000, debug = FALSE, 
    technical = list(), ...)
{     
    mod <- confmirt(data, nfact, guess=guess, rotate=rotate, verbose=verbose,
                    calcLL=calcLL, draws=draws, debug=debug, technical=technical, ...)    
    mod@Call <- match.call()
    mod
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
					x@df, ", p = ", round(x@p,4), "\nTLI = ", round(x@TLI,3),
                    ", RMSEA = ", round(x@RMSEA,3), "\n", sep="")
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
					object@df, ", p = ", round(object@p,4), "\nTLI = ", round(object@TLI,3),
                    ", RMSEA = ", round(object@RMSEA,3), 
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
        browser()
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
			parameters <- cbind(a,d,object@guess, object@upper)
			colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""), 
                                      paste("d_",1:max(K-1),sep=""),"guess","upper")   
			rownames(parameters) <- colnames(object@data)
			cat("\nParameter slopes and intercepts: \n\n")	
			print(round(parameters, digits))	  
		}
		ret <- list(parameters)
		if(length(object@SEpars) > 1){
			if(SE){
				cat("\nStd. Errors: \n\n")	
				SEs <- cbind(object@SEpars,object@SEg,object@SEup)
				colnames(SEs) <- c(paste("a_",1:ncol(a),sep=""),paste("d_",1:max(K-1),sep=""),
                                   "SEguess", "SEupper") 	
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
		upper <- x@upper
		upper[is.na(upper)] <- 1
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
				P <- P.mirt(a[j,], d[[j]],Theta, guess[j], upper[j])
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

