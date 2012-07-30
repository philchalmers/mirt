#' Confirmatory Full-Information Item Factor Analysis for Mixed Data Formats
#' 
#' \code{confmirt} fits a conditional (i.e., confirmatory) full-information
#' maximum-likelihood factor analysis model to dichotomous and polychotomous
#' data under the item response theory paradigm using Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm. If requested, lower and upper asymptote
#' parameters are estimated with a beta priors included automatically.
#' 
#' 
#' \code{confmirt} follows a confirmatory item factor analysis strategy that
#' uses a stochastic version of maximum likelihood estimation described by Cai
#' (2010). The general equation used for multidimensional item response theory
#' in this function is in the logistic form with a scaling correction of 1.702.
#' This correction is applied to allow comparison to mainstream programs such
#' as TESTFACT (2003) and POLYFACT. Missing data are treated as 'missing at
#' random' so that each response vector is included in the estimation (i.e.,
#' full-information). Residuals are computed using the LD statistic (Chen &
#' Thissen, 1997) in the lower diagonal of the matrix returned by
#' \code{residuals}, and Cramer's V above the diagonal. For computing the
#' log-likelihood more accurately see \code{\link{logLik}}.
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
#' \code{coef} displays the item parameters with their associated standard
#' errors, while use of \code{summary} transforms the slopes into a factor
#' loadings metric. Also, nested models may be compared by using the
#' \code{anova} function, where a Chi-squared difference test and AIC/BIC
#' difference values are displayed.
#' 
#' @aliases confmirt coef,confmirt-method summary,confmirt-method
#' residuals,confmirt-method anova,confmirt-method fitted,confmirt-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data
#' @param model an object returned from \code{confmirt.model()} declarating how
#' the factor model is to be estimated, or a single numeric value indicating the number 
#' of exploratory factors to estimate. See \code{\link{confmirt.model}} for
#' more details
#' @param guess initial (or fixed) values for the pseudo-guessing parameter. Can be 
#' entered as a single value to assign a global guessing parameter or may be entered as
#' a numeric vector for each item
#' @param upper initial (or fixed) upper bound parameters for 4-PL model. Can be 
#' entered as a single value to assign a global upper bound parameter or may be entered as a 
#' numeric vector corresponding to each item
#' @param estGuess a logical vector indicating which lower-asymptote parameters
#' to be estimated (default is null, and therefore is contingent on the values
#' in \code{guess}). By default, if any value in \code{guess} is greater than 0
#' then its respective \code{estGuess} value is set to \code{TRUE}.
#' Additionally, beta priors are automatically imposed for estimated parameters
#' which correspond to the input guessing values.
#' @param estUpper same function as \code{estGuess}, but for upper bound parameters
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
#' Cai, L. (2010). Metropolis-Hastings Robbins-Monro algorithm for confirmatory
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
#' model.1 <- confmirt.model()
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
#' #fix first slope at 1.5, and set slopes 7 & 8 to be equal
#' model.2 <- confmirt.model()
#'   F1 = 1-4
#'   F2 = 4-8
#'   COV = F1*F2
#'   SLOPE = F1@@1 eq 1.5, F2@@7 eq F2@@8
#'   
#' 
#' mod2 <- confmirt(dataset,model.2)
#' anova(mod2,mod1)
#' 
#' #####
#' #bifactor 	
#' model.3 <- confmirt.model()
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
#' model.linear <- confmirt.model()
#'       F = 1-8
#'
#' 
#' model.quad <- confmirt.model()
#'       F = 1-8
#'   (F*F) = 1-8
#'
#' 
#' model.cube <- confmirt.model()
#'         F = 1-8
#'     (F*F) = 1-8
#'   (F*F*F) = 1-8
#'
#' 
#' model.combo <- confmirt.model()
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
confmirt <- function(data, model, guess = 0, upper = 1, estGuess = NULL, estUpper = NULL, 
    verbose = TRUE, calcLL = TRUE, draws = 2000, returnindex = FALSE, debug = FALSE, 
    rotate = 'varimax', Target = NULL, technical = list(),  ...)
{		
	Call <- match.call()       
    ##########
    if(any(upper < 1)) stop('Upper bound estimation is not currently available.')
    ##########
	set.seed(12345)	
	itemnames <- colnames(data)
	keywords <- c('SLOPE','INT','COV','MEAN','PARTCOMP','PRIOR')
	data <- as.matrix(data)		
	colnames(data) <- itemnames	
	J <- ncol(data)
	N <- nrow(data)
	exploratory <- FALSE
    if(is(model, 'numeric')){
        tmp <- tempfile('tempfile')
        cat(paste('F',1:model,' = 1-', J, "\n", sep=''), file=tmp)
        model <- confmirt.model(tmp, quiet = TRUE)
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
	guess.prior.n <- ifelse(!is.null(technical$guess.prior.n), technical$guess.prior.n,
		rep(50,J))
	if(length(guess.prior.n) == 1) guess.prior.n <- rep(guess.prior.n,J)	
	if(length(guess.prior.n) != J) 
		stop('technical$guess.prior.n does not have the same number of values as items')
	guess.prior.n[is.na(guess.prior.n)] <- 50
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
	if(is.null(estGuess))
		estGuess <- guess > 0
	if(is.null(estUpper))
	    estUpper <- upper < 1
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
		if(setequal(uniques[[i]], c(0,1))){
			fulldata[ ,itemloc[ind]:(itemloc[ind]+1)] <- cbind(data[ ,ind],abs(1-data[ ,ind]))			
			next
		}
		dummy <- matrix(0,N,K[ind])
		for (j in 0:(K[ind]-1))  
			dummy[,j+1] <- as.integer(data[,ind] == uniques[[ind]][j+1])  		
		fulldata[ ,itemloc[ind]:(itemloc[ind+1]-1)] <- dummy			
	}	
	fulldata[is.na(fulldata)] <- 0	
	mod <- model.elements(model=model, factorNames=factorNames, nfactNames=nfactNames, nfact=nfact, 
                          J=J, K=K, fulldata=fulldata, itemloc=itemloc, data=data, N=N, 
                          estGuess=estGuess, guess=guess, upper=upper, estUpper=estUpper, 
                          guess.prior.n=guess.prior.n, itemnames=itemnames, exploratory=exploratory)
	parcount <- mod$parcount
	npars <- mod$npars
	if(returnindex) return(parcount)
	if(debug) print(mod)      
    if(any(rowSums(mod$est$estlam) == 0)){
		tmp <- 1:J
		tmp <- tmp[rowSums(mod$est$estlam) == 0]
		stop('Item(s) ', paste(tmp,''), 'have no factor loadings specified.')
	}	    
	if(exploratory){
        lamind <- mod$ind$lamind
	    Rpoly <- cormod(na.omit(data),K,guess)
	    FA <- psych::fa(Rpoly, nfact, rotate = 'none', warnings= FALSE, fm="minres")    
	    loads <- unclass(loadings(FA))
	    u <- FA$unique
	    u[u < .001 ] <- .2
	    cs <- sqrt(u)
	    lambdas <- mod$val$lambdas <- loads/cs        
        mod$val$pars[lamind] <- t(lambdas)
	}        
	ESTIMATE <- MHRM(mod=mod, NCYCLES=NCYCLES, BURNIN=BURNIN, KDRAWS=KDRAWS, SEMCYCLES=SEMCYCLES,  
                     TOL=TOL, nfactNames=nfactNames, itemloc=itemloc, fulldata=fulldata, nfact=nfact,
                     N=N, gain=gain, K=K, J=J, verbose=verbose)	    	
	if(verbose) cat("\n\n")
	SEtmp <- diag(solve(ESTIMATE$info))		
	if(any(SEtmp < 0)){
		warning("Information matrix is not positive definite, negative SEs set to 'NA'.\n")
		SEtmp[SEtmp < 0] <- NA
	}	    
    pars <- ESTIMATE$pars
	SEtmp <- sqrt(SEtmp)	    
	SE <- rep(NA,npars) 
	SE[mod$ind$parind[mod$ind$sind]] <- SEtmp
	SE[mod$val$constvalues[ ,1]==1] <- NA
	if(length(mod$ind$equalconstr) > 0)
		for(i in 1:length(mod$ind$equalconstr))
			SE[mod$ind$equalconstr[[i]]] <- mean(SE[mod$ind$equalconstr[[i]]])
	estpars <- pars[mod$ind$sind]
	lambdas <- ESTIMATE$normpars$lambdas
	lambdas[!mod$est$estlam & !lambdas != 0] <- NA	
	guess <- rep(NA,J)
	guess <- pars[mod$ind$guessind]
	guess[!mod$est$estGuess] <- NA
	guess[K == 2 & !mod$est$estGuess] <- 0
	upper <- rep(NA,J)
	upper <- pars[mod$ind$upperind]
	upper[!mod$est$estUpper] <- NA
	upper[K == 2 & !mod$est$estUpper] <- 1
	zetas <- pars[mod$ind$zetaind]
	u <- pars[mod$ind$meanind]	
	sig <- matrix(0,nfact,nfact)
	SElam <- matrix(SE[mod$ind$lamind],J,nfactNames,byrow=TRUE)
	SEzetas <- SE[mod$ind$zetaind]	
	SEg <- rep(NA,J)	
	SEg <- SE[mod$ind$guessind]	
	SEup <- SE[mod$ind$upperind]
	SEg[!mod$est$estGuess] <- NA
	SEup[!mod$est$estUpper] <- NA
	SEu <- SE[mod$ind$meanind]	
	SEsig <- matrix(0,nfact,nfact)	
	tmp <- pars[mod$ind$covind]
	tmp2 <- SE[mod$ind$covind]
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
	if(any(mod$est$estComp)){
		if((max(K)-1) > nfactNames) tmp1 <- tmp2 <- matrix(NA,J,(max(K)-1))
		else tmp1 <- tmp2 <- matrix(NA,J,nfactNames)
	} else tmp1 <- tmp2 <- matrix(NA,J,(max(K)-1))	
	
	#reload zetas to matrix
	zetas <- pars[mod$ind$zetaind]    
	loc <- 1
	for(i in 1:J){
		if(!mod$est$estComp[i]){
			for(j in 1:(K[i]-1)){
				tmp1[i,j] <- zetas[loc] 
				tmp2[i,j] <- SEzetas[loc]
				loc <- loc + 1
			}
		} else {
			for(j in 1:nfactNames){
				tmp1[i,j] <- zetas[loc]
				tmp2[i,j] <- SEzetas[loc]
				loc <- loc + 1
			}	
		}	
	}	 
	zetas <- tmp1
	SEzetas <- tmp2	
	parsprint <- cbind(ESTIMATE$normpars$lambdas, zetas)
	SEpars <- cbind(SElam, SEzetas)
	gpars <- list(u = ESTIMATE$normpars$mu, sig = ESTIMATE$normpars$sig)	
	SEgpars <- list(SEu = SEu, SEsig = SEsig)
	estpars <- mod$est
		
	if (nfactNames > 1){
        norm <- sqrt(1 + rowSums(lambdas[ ,1:nfactNames]^2,na.rm = TRUE))
	} else norm <- as.matrix(sqrt(1 + lambdas[ ,1]^2))  
	F <- as.matrix(lambdas[ ,1:nfactNames]/norm)
	F[is.na(F)] <- 0		
	h2 <- rowSums(F^2)
	colnames(F) <- factorNames
	names(h2) <- itemnames  
	null.mod <- unclass(new('mirtClass'))
    if(!any(is.na(data))) null.mod <- unclass(mirt(data, 0))
    
    if(exploratory){
        ret <- new('polymirtClass',pars=ESTIMATE$normpars, guess=guess, SEpars=SEpars, SEg=SEg,
                   upper=upper, SEup=SEup, cycles=ESTIMATE$cycles, Theta=ESTIMATE$theta0, 
                   fulldata=fulldata, data=data, K=K, F=F, h2=h2, itemloc=itemloc, 
                   converge=ESTIMATE$converge, estGuess=estGuess, rotate=rotate, null.mod=null.mod, 
                   Target=Target, Call=Call)
    } else {
    	ret <- new('confmirtClass', pars=ESTIMATE$normpars, parsprint=parsprint, guess=guess, upper=upper, 
                   SEg=SEg, SEup=SEup, gpars=gpars, SEgpars=SEgpars, estpars=estpars, K=K, 
                   itemloc=itemloc, cycles=ESTIMATE$cycles, Theta=ESTIMATE$theta0, 
                   fulldata=fulldata, data=data, h2=h2, F=F, converge=ESTIMATE$converge, 
                   nconstvalues=as.integer(mod$nconstvalues), SEpars=SEpars, estComp=mod$est$estComp, 
                   prodlist=as.list(mod$ind$prodlist), null.mod=null.mod, Call=Call)
    }
	if(calcLL){
		if(verbose) cat("Calculating log-likelihood...\n")
		flush.console()
		ret <- logLik(ret, draws)		        
	}	
	return(ret)
}
