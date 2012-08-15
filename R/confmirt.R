#' Confirmatory Full-Information Item Factor Analysis for Mixed Data Formats
#' 
#' \code{confmirt} fits a conditional (i.e., confirmatory) full-information
#' maximum-likelihood factor analysis model to dichotomous and polychotomous
#' data under the item response theory paradigm using Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm. If requested, lower and upper asymptote
#' parameters are estimated with a beta priors included automatically.
#'  
#' \code{confmirt} follows a confirmatory and exploratory item factor analysis strategy that
#' uses a stochastic version of maximum likelihood estimation described by Cai
#' (2010a, 2010b). The general equation used for multidimensional item response theory
#' in this function is in the logistic form with a scaling correction of 1.702.
#' This correction is applied to allow comparison to mainstream programs such
#' as TESTFACT (2003) and POLYFACT. Missing data are treated as 'missing at
#' random' so that each response vector is included in the estimation (i.e.,
#' full-information). Residuals are computed using the LD statistic (Chen &
#' Thissen, 1997) in the lower diagonal of the matrix returned by
#' \code{residuals}, and Cramer's V above the diagonal. For computing the
#' log-likelihood more accurately see \code{\link{logLik}}.
#' 
#' #' \code{coef} displays the item parameters with their associated standard
#' errors, while use of \code{summary} transforms the slopes into a factor
#' loadings metric and if the model is exploratory allows for rotating the parameters. 
#' Also, nested models may be compared by using the
#' \code{anova} function, where a Chi-squared difference test and AIC/BIC
#' difference values are displayed.
#' 
#' @section Confirmatory IRT
#' 
#' Specification of the confirmatory item factor analysis model follows many of
#' the rules in the SEM framework for confirmatory factor analysis. The
#' variances of the latent factors are automatically fixed to 1 to help
#' facilitate model identification. All parameters may be fixed to constant
#' values or set equal to other parameters using the appropriate declarations.
#' Guessing parameters may be specified for dichotomous items and are estimated
#' with beta priors automatically, and if a guessing parameter is declared for
#' a polychotomous item it is ignored.
#' 
#' @section Exploratory IRT
#' 
#' Specifying a number as the second input to confmirt an exploratory IRT model is estimatated and 
#' can be viewed as a stocastic analogue of \code{mirt}, with much of the same behaviour and 
#' specifications. Rotatation and target matrix options will be used in this subroutine and will be
#' passed to the returned object for use in generic functions such as \code{summary()} and 
#' \code{fscores}. Again, factor means and variances are fixed to ensure proper identification. See
#' \code{\link{mirt}} for more details.
#' 
#' 
#' @aliases confmirt coef,confmirt-method summary,confmirt-method
#' residuals,confmirt-method anova,confmirt-method fitted,confmirt-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data
#' @param model an object returned from \code{specifyModel()} declarating how
#' the factor model is to be estimated, or a single numeric value indicating the number 
#' of exploratory factors to estimate. See \code{\link{specifyModel}} for
#' more details
#' @param guess initial (or fixed) values for the pseudo-guessing parameter. Can be 
#' entered as a single value to assign a global guessing parameter or may be entered as
#' a numeric vector for each item
#' @param upper initial (or fixed) upper bound parameters for 4-PL model. Can be 
#' entered as a single value to assign a global upper bound parameter or may be entered as a 
#' numeric vector corresponding to each item
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#' option. Only prints patterns that have standardized residuals greater than 
#' \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param verbose logical; display iteration history during estimation?
#' @param calcLL logical; calculate the log-likelihood via Monte Carlo
#' integration?
#' @param draws the number of Monte Carlo draws to estimate the log-likelihood
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param returnindex logical; return the list containing the item paramter
#' locations? To be used when specifying prior parameter distributions
#' @param debug logical; turn on debugging features?
#' @param object an object of class \code{confmirtClass}
#' @param object2 an object of class \code{confmirtClass}
#' @param SE logical; print standard errors?
#' @param print.gmeans logical; print latent factor means?
#' @param digits the number of significant digits to be rounded
#' @param rotate if \code{model} is numeric (indicating an exploratory item FA) then this 
#' rotation is used. Default is \code{'varimax'}
#' @param Target a dummy variable matrix indicing a target rotation pattern
#' @param technical list specifying subtle parameters that can be adjusted. These 
#' values are 
#' \describe{
#' \item{NCYCLES}{max number of MH-RM cycles; default 2000}
#' \item{BURNIN}{number of burn in cycles (stage 1); default 150}
#' \item{SEMCYCLES}{number of SEM cycles (stage 2); default 50}
#' \item{KDRAWS}{number of paralell MH sets to be drawn; default 1}
#' \item{TOL}{minimum threshold tolerance for convergence of MH-RM, must occur on three consecutive
#' occations; default .001} 
#'   \item{set.seed}{seed number used during estimation. Default is 12345}
#' 	 \item{guess.prior.n}{a scalar or vector for the weighting of the beta priors for 
#'		guessing parameters (default is 50, typical ranges are from 2 to 500). If a 
#'      scalar is specified this is used globally, otherwise a numeric vector of size
#' 	    \code{ncol(data)} can be used to correspond to particualr items (NA values use 
#'      the default)} 
#'   \item{gain}{a vector of three values specifying the numerator, exponent, and subtracted
#'      values for the RM gain value. Default is \code{c(0.05,0.5,0.004)}}   	
#' }
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{simdata}},
#' \code{\link{fscores}}, \code{\link{confmirt.model}}
#' @references
#' 
#' Cai, L. (2010a). High-Dimensional exploratory item factor analysis by a
#' Metropolis-Hastings Robbins-Monro algorithm. \emph{Psychometrika, 75},
#' 33-57.
#' 
#' Cai, L. (2010b). Metropolis-Hastings Robbins-Monro algorithm for confirmatory
#' item factor analysis. \emph{Journal of Educational and Behavioral
#' Statistics, 35}, 307-335.
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
#' confmirt(data, model, guess = 0, upper = 1, estGuess = NULL, estUpper = NULL, 
#' verbose = TRUE, calcLL = TRUE, draws = 2000, returnindex = FALSE, debug = FALSE, 
#' rotate = 'varimax', Target = NULL, technical = list(),  ...)
#' 
#' \S4method{coef}{confmirt}(object, SE = TRUE, print.gmeans = FALSE, digits = 3, ...)
#' 
#' \S4method{summary}{confmirt}(object, digits = 3, ...)
#' 
#' \S4method{residuals}{confmirt}(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
#' 
#' \S4method{anova}{confmirt}(object, object2, ...)
#'
#' \S4method{fitted}{confmirt}(object, digits = 3, ...)
#'
#' @export confmirt
#' @examples
#'  
#' \dontrun{
#' #Exploratory model estimation, similar to mirt()
#' data(LSAT7)
#' fulldata <- expand.table(LSAT7) 
#' (mod1 <- confmirt(fulldata, 1))
#' 
#' #Confirmatory models
#' 
#' #simulate data
#' a <- matrix(c(
#' 1.5,NA,
#' 0.5,NA,
#' 1.0,NA,
#' 1.0,0.5,
#'  NA,1.5,
#'  NA,0.5,
#'  NA,1.0,
#'  NA,1.0),ncol=2,byrow=TRUE)
#' 
#' d <- matrix(c(
#' -1.0,NA,NA,
#' -1.5,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA,
#' 3.0,2.0,-0.5,
#' 2.5,1.0,-1,
#' 2.0,0.0,NA,
#' 1.0,NA,NA),ncol=3,byrow=TRUE)
#' 
#' sigma <- diag(2)
#' sigma[1,2] <- sigma[2,1] <- .4
#' dataset <- simdata(a,d,2000,sigma)
#' 
#' #analyses
#' #CIFA for 2 factor crossed structure
#' 
#' model.1 <- specifyModel()
#'   F1 = 1-4
#'   F2 = 4-8
#'   COV = F1*F2
#' 
#' 
#' mod1 <- confmirt(dataset,model.1)
#' coef(mod1)
#' summary(mod1)
#' residuals(mod1)
#' 
#' #####
#' #bifactor 	
#' model.3 <- specifyModel()
#'   G = 1-8
#'   F1 = 1-4
#'   F2 = 5-8
#' 
#' 
#' mod3 <- confmirt(dataset,model.3)
#' coef(mod3)
#' summary(mod3)
#' residuals(mod3)
#' anova(mod1,mod3)
#'
#' #####
#' #polynomial and combinations
#' model.linear <- specifyModel()
#'       F = 1-8
#'
#' 
#' model.quad <- specifyModel()
#'       F = 1-8
#'   (F*F) = 1-8
#'
#' 
#' model.cube <- specifyModel()
#'         F = 1-8
#'     (F*F) = 1-8
#'   (F*F*F) = 1-8
#'
#' 
#' model.combo <- specifyModel()
#'        F1 = 1-4
#'        F2 = 5-8
#'   (F1*F2) = 1-8
#'
#' 
#' mod.linear <- confmirt(dataset, model.linear)
#' mod.quad <- confmirt(dataset, model.quad)
#' mod.cube <- confmirt(dataset, model.cube)
#' mod.combo <- confmirt(dataset, model.combo)
#'
#' anova(mod.linear,mod.quad)
#' anova(mod.linear,mod.cube)
#' anova(mod.linear,mod.combo)
#' anova(mod.cube,mod.combo)
#' }
#' 
confmirt <- function(data, model, itemtype = NULL, guess = 0, upper = 1, startvalues = NULL, 
                     constrain = NULL, freepars = NULL, parprior = NULL, verbose = TRUE, calcLL = TRUE, 
                     draws = 2000, debug = FALSE, rotate = 'varimax', Target = NULL, 
                     technical = list(),  ...)
{		
	Call <- match.call()           
    ##########
    if(any(upper < 1)) stop('Upper bound estimation is not currently available.')
    ##########
	set.seed(12345)	
	itemnames <- colnames(data)
	keywords <- c('COV')
	data <- as.matrix(data)		
	colnames(data) <- itemnames	
	J <- ncol(data)
	N <- nrow(data)
	exploratory <- FALSE
    if(is(model, 'numeric')){
        tmp <- tempfile('tempfile')
        cat(paste('F',1:model,' = 1-', J, "\n", sep=''), file=tmp)
        model <- specifyModel(tmp, quiet = TRUE)
        exploratory <- TRUE
        unlink(tmp)
    }    
	##technical
	NCYCLES <- ifelse(is.null(technical$NCYCLES), 2000, technical$NCYCLES)
    BURNIN <- ifelse(is.null(technical$BURNIN), 150, technical$BURNIN)
    SEMCYCLES <- ifelse(is.null(technical$SEMCYCLES), 50, technical$SEMCYCLES)
    KDRAWS  <- ifelse(is.null(technical$KDRAWS), 1, technical$KDRAWS)
    TOL <- ifelse(is.null(technical$TOL), .001, technical$TOL)        
	if(!is.null(technical$set.seed)) set.seed(technical$set.seed)	
	gain <- c(0.05,0.5,0.004)
	if(!is.null(technical$gain)) {
		if(length(technical$gain) == 3 && is.numeric(technical$gain))
			gain <- technical$gain
	}
	##
	Target <- ifelse(is.null(Target), NaN, Target)
	if(length(guess) == 1) guess <- rep(guess,J)
	if(length(guess) > J || length(guess) < J) 
		stop("The number of guessing parameters is incorrect.")					
	if(length(upper) == 1) upper <- rep(upper,J)
	if(length(upper) > J || length(upper) < J) 
	    stop("The number of upper bound parameters is incorrect.")
	uniques <- list()
	for(i in 1:J)
		uniques[[i]] <- sort(unique(data[,i]))
	K <- rep(0,J)
	for(i in 1:J) K[i] <- length(uniques[[i]])	
	guess[K > 2] <- 0
	upper[K > 2] <- 1		
	if(is.null(itemtype)) {
	    itemtype <- rep('', J)
	    for(i in 1:J){
	        if(K[i] > 2) itemtype[i] <- 'graded'
	        if(K[i] == 2) itemtype[i] <- '2PL'                            
	    }        
	} 
	if(length(itemtype) != J) stop('itemtype specification is not the correct length')
	if(length(itemtype) == 1) itemtype <- rep(itemtype, J)
	itemloc <- cumsum(c(1,K))	
	model <- matrix(model$x,ncol=2)
	factorNames <- setdiff(model[,1],keywords)
	nfactNames <- length(factorNames)
	nfact <- sum(!grepl('\\(',factorNames))
	index <- 1:J	
    fulldata <- matrix(0,N,sum(K))
	Names <- NULL
	for(i in 1:J)
	    Names <- c(Names, paste("Item.",i,"_",1:K[i],sep=""))				
	colnames(fulldata) <- Names			
	for(i in 1:J){
	    ind <- index[i]		
	    dummy <- matrix(0,N,K[ind])
	    for (j in 0:(K[ind]-1))  
	        dummy[,j+1] <- as.integer(data[,ind] == uniques[[ind]][j+1])  		
	    fulldata[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- dummy		
	}	
	fulldata[is.na(fulldata)] <- 0    
    parnumber <- 1 #to be used later when looping over more than 1 group    
	pars <- model.elements(model=model, itemtype=itemtype, factorNames=factorNames, 
                           nfactNames=nfactNames, nfact=nfact, J=J, K=K, fulldata=fulldata, 
                           itemloc=itemloc, data=data, N=N, guess=guess, upper=upper,  
                           itemnames=itemnames, exploratory=exploratory, constrain=constrain,
                           startvalues=startvalues, freepars=freepars, parprior=parprior, 
                           parnumber=parnumber)    
	#Contraints, startvalues, and estimation
	if(!is.null(constrain)){
	    if(constrain == 'index'){
	        returnedlist <- list()                        
	        for(i in 1:J)
	            returnedlist[[i]] <- pars[[i]]@parnum 
	        names(returnedlist) <- itemnames
	        return(returnedlist)
	    }
	}    
	if(!is.null(startvalues)){
	    if(startvalues == 'index'){
	        returnedlist <- list()                        
	        for(i in 1:J){
	            par <- pars[[i]]@par
	            names(par) <- names(pars[[i]]@parnum)
	            returnedlist[[i]] <- par
	        }
	        names(returnedlist) <- itemnames
	        return(returnedlist)
	    }
	}
	if(!is.null(freepars)){
	    if(freepars == 'index'){
	        returnedlist <- list()                        
	        for(i in 1:J){
	            est <- pars[[i]]@est
	            names(est) <- names(pars[[i]]@parnum)
	            returnedlist[[i]] <- est
	        }
	        names(returnedlist) <- itemnames
	        return(returnedlist)
	    }
	}    
    npars <- 0
    for(i in 1:length(pars))
        npars <- npars + sum(pars[[i]]@est)			        
	if(exploratory){        
	    Rpoly <- cormod(na.omit(data),K,guess)
	    FA <- psych::fa(Rpoly, nfact, rotate = 'none', warnings= FALSE, fm="minres")    
	    loads <- unclass(loadings(FA))
	    u <- FA$unique
	    u[u < .001 ] <- .2
	    cs <- sqrt(u)
	    lambdas <- loads/cs        
        for(i in 1:J)
            pars[[i]]@par[1:nfact] <- lambdas[i, ]        
	}        
	if(debug) browser()       
	ESTIMATE <- MHRM(pars=pars, NCYCLES=NCYCLES, BURNIN=BURNIN, SEMCYCLES=SEMCYCLES, KDRAWS=KDRAWS,
                     TOL=TOL, gain=gain, nfactNames=nfactNames, itemloc=itemloc, fulldata=fulldata, 
                     nfact=nfact, N=N,  K=K, J=J, npars=npars, constrain=constrain, verbose=verbose)
    pars <- ESTIMATE$pars
	if(verbose) cat("\n\n")    
	lambdas <- Lambdas(pars)
	if (nfactNames > 1){
        norm <- sqrt(1 + rowSums(lambdas[ ,1:nfactNames]^2,na.rm = TRUE))
	} else norm <- as.matrix(sqrt(1 + lambdas[ ,1]^2))  
	F <- as.matrix(lambdas[ ,1:nfactNames]/norm)
	F[is.na(F)] <- 0		
	h2 <- rowSums(F^2)
	colnames(F) <- factorNames
	names(h2) <- itemnames  
	null.mod <- unclass(new('mirtClass'))
    if(!any(is.na(data))) null.mod <- unclass(mirt(data, 0, itemtype = 'NullModel'))
    
    ret <- new('confmirtClass', pars=pars, K=K, itemloc=itemloc, cycles=ESTIMATE$cycles,                
               fulldata=fulldata, data=data, h2=h2, F=F, converge=ESTIMATE$converge,                 
               null.mod=null.mod, Call=Call)    
	if(calcLL){
		if(verbose) cat("Calculating log-likelihood...\n")
		flush.console()
		ret <- logLik(ret, draws)		        
	}	
	return(ret)
}
