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
    
	#pars
	pars <- mod$val$pars
	lambdas <- mod$val$lambdas
	zetas <- mod$val$zetas
	gmeans <- mod$val$gmeans
	gcov <- mod$val$gcov
	constvalues <- mod$val$constvalues      

	#est
	estlam <- mod$est$estlam
	estComp <- mod$est$estComp		
	estgcov <- mod$est$estgcov
	estgmeans <- mod$est$estgmeans	
  
	#ind
	parind <- mod$ind$parind	
	equalconstr <- mod$ind$equalconstr	
	prodlist <- mod$ind$prodlist
	parpriors <- mod$ind$parpriors
	sind <- mod$ind$sind
	lamind <- mod$ind$lamind
	zetaind <- mod$ind$zetaindlist
	guessind <- mod$ind$guessind
	upperind <- mod$ind$upperind
	groupind <- mod$ind$groupind
	meanind <- mod$ind$meanind
	covind <- mod$ind$covind    
	if(any(rowSums(estlam) == 0)){
		tmp <- 1:J
		tmp <- tmp[rowSums(estlam) == 0]
		stop('Item(s) ', paste(tmp,''), 'have no factor loadings specified.')
	}
	indlist <- mod$ind 
    
	if(exploratory){
	    Rpoly <- cormod(na.omit(data),K,guess)
	    FA <- psych::fa(Rpoly,nfact,rotate = 'none', warnings= FALSE, fm="minres")    
	    loads <- unclass(loadings(FA))
	    u <- FA$unique
	    u[u < .001 ] <- .2
	    cs <- sqrt(u)
	    lambdas <- mod$val$lambdas <- loads/cs        
        pars[lamind] <- t(lambdas)
	}
	
	#Preamble for MRHM algorithm
	pars[constvalues[,1] == 1] <- constvalues[constvalues[,1] == 1,2]
	theta0 <- matrix(0,N,nfact)	    
	cand.t.var <- 1			
	tmp <- .1
	for(i in 1:30){			
		theta0 <- draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, upper=upper,
                              fulldata=fulldata, K=K, itemloc=itemloc, cand.t.var=cand.t.var, 
		                      prior.t.var=gcov, prior.mu=gmeans, estComp=estComp, prodlist=prodlist)
		if(i > 5){		
			if(attr(theta0,"Proportion Accepted") > .35) cand.t.var <- cand.t.var + 2*tmp 
			else if(attr(theta0,"Proportion Accepted") > .25 && nfact > 3) 
                cand.t.var <- cand.t.var + tmp
			else if(attr(theta0,"Proportion Accepted") < .2 && nfact < 4) 
                cand.t.var <- cand.t.var - tmp
			else if(attr(theta0,"Proportion Accepted") < .1) 
                cand.t.var <- cand.t.var - 2*tmp
			if (cand.t.var < 0){
				cand.t.var <- tmp		
				tmp <- tmp / 2
			}		
		}
	}	
	m.thetas <- grouplist <- list()		
	SEM.stores <- matrix(0,SEMCYCLES,npars)
	SEM.stores2 <- list()
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
	noninvcount <- 0		
	if(length(equalconstr) > 0)	
		for(i in 1:length(equalconstr))
			nconstvalues <- nconstvalues + length(equalconstr[[i]]) - 1			
			
	####Big MHRM loop 
	for(cycles in 1:(NCYCLES + BURNIN + SEMCYCLES))								
	{ 
		if(cycles == BURNIN + 1) stagecycle <- 2			
		if(stagecycle == 3)
			gamma <- (gain[1]/(cycles - SEMCYCLES - BURNIN - 1))^(gain[2]) - gain[3]
		if(cycles == (BURNIN + SEMCYCLES + 1)){ 
			stagecycle <- 3		
		    pars <- rep(0,npars)
			for(i in 1:SEMCYCLES){
				pars <- pars + SEM.stores[i,]
				Tau <- Tau + SEM.stores2[[i]]
			}	
			pars <- pars/SEMCYCLES	
			Tau <- Tau/SEMCYCLES	
			k <- KDRAWS	
			gamma <- .25
		}	
				
		normpars <- sortParsConfmirt(pars, indlist, nfact, estGuess, nfactNames)		
		lambdas <- normpars$lambdas
		zetas <- normpars$zetas
		guess <- normpars$guess	
        upper <- normpars$upper
		grouplist$u <- mu <- normpars$mu					
		grouplist$sig <- sig <- normpars$sig			
		
		#Step 1. Generate m_k datasets of theta 
		for(j in 1:4) 
            theta0 <- draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, 
                            upper=upper, fulldata=fulldata, K=K, itemloc=itemloc, 
                            cand.t.var=cand.t.var, prior.t.var=sig, prior.mu=mu, estComp=estComp, 
                            prodlist=prodlist)
		for(i in 1:k) m.thetas[[i]] <- 
            draw.thetas(theta0=theta0, lambdas=lambdas, zetas=zetas, guess=guess, 
                                upper=upper, fulldata=fulldata, K=K, itemloc=itemloc, 
                                cand.t.var=cand.t.var, prior.t.var=sig, prior.mu=mu, 
                                estComp=estComp, prodlist=prodlist)
		theta0 <- m.thetas[[1]]
		        
		#Step 2. Find average of simulated data gradients and hessian 		
		g.m <- h.m <- group.m <- list()
		g <- rep(0, npars)
		h <- matrix(0, npars, npars)	
		for (j in 1:k) {
            g <- rep(NA, npars)            
			thetatemp <- m.thetas[[j]]
			if(!is.null(prodlist)) thetatemp <- prodterms(thetatemp,prodlist)	
            for (i in 1:J){			
				if(estComp[i]){
					if (estGuess[i]) {
						temp <- dpars.comp(lambdas[i,][estlam[i,]], zetas[[i]], 
							guess[i], fulldata[, itemloc[i]], thetatemp, estGuess[i])
						ind <- parind[is.na(g)][1]
						ind2 <- ind + length(temp$grad) - 1
						g[ind:ind2] <- temp$grad
						h[ind:ind2, ind:ind2] <- temp$hess						
					} else {
						temp <- dpars.comp(lambdas[i,][estlam[i,]], zetas[[i]], 
							guess[i], fulldata[, itemloc[i]], thetatemp)
						ind <- parind[is.na(g)][1]							
						ind2 <- ind + length(zetas[[i]])*2 - 1
						g[ind:ind2] <- temp$grad
						h[ind:ind2, ind:ind2] <- temp$hess
						g[ind2 + 1] <- 0	
					}				
					next
				}
                if(K[i] == 2){
					temp <- dpars.dich(lambda=lambdas[i, ], zeta=zetas[[i]], g=guess[i], u=upper[i],
						dat=fulldata[ ,itemloc[i]], Thetas=thetatemp, estGuess=estGuess[i])
					ind <- parind[is.na(g)][1]					
					ind2 <- ind + length(temp$g) - 1		
					if(!estGuess[i]) g[ind2 + 1] <- 0
					if(!estUpper[i]) g[ind2 + 2] <- 0
					g[ind:ind2] <- temp$grad
					h[ind:ind2,ind:ind2] <- temp$hess						
				} else {						
					temp <- dpars.poly(lambdas[i, ],zetas[[i]],
						fulldata[ ,itemloc[i]:(itemloc[i+1]-1)],thetatemp)
					ind <- parind[is.na(g)][1]	
					ind2 <- ind + length(temp$g) - 1		
					g[ind:ind2] <- temp$grad
					h[ind:ind2,ind:ind2] <- temp$hess
					g[ind2 + 1] <- g[ind2 + 2] <- 0	#zeros for guess + upper
				}
            }
			g[is.na(g)] <- 0
			tmp <- d.group(grouplist,as.matrix(thetatemp[ ,1:nfact]))
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
		if(length(parpriors) > 0){
			for(i in 1:length(parpriors)){
				tmp <- parpriors[[i]]
				if(tmp[1] == 1){
					grad[tmp[2]] <- grad[tmp[2]] - (pars[tmp[2]] - tmp[3])/ tmp[4]^2
					ave.h[tmp[2],tmp[2]] <- ave.h[tmp[2],tmp[2]] +  1/tmp[4]^2
				}				
				else if(tmp[1] == 2){		
					tmp2 <- betaprior(tmp[3],tmp[4],pars[tmp[2]])					
					grad[tmp[2]] <- grad[tmp[2]] + tmp2$g
					ave.h[tmp[2],tmp[2]] <- ave.h[tmp[2],tmp[2]] + tmp2$h
				}				
			}
		}		
		grad <- grad[parind[sind]]		
		ave.h <- ave.h[parind[sind],parind[sind]] 
		if(is.na(attr(theta0,"log.lik"))) stop('Estimation halted. Model did not converge.')		
		if(verbose){
			if((cycles + 1) %% 10 == 0){
				if(cycles < BURNIN)
					cat("Stage 1: Cycle = ", cycles + 1, ", Log-Lik = ", 
						sprintf("%.1f",attr(theta0,"log.lik")), sep="")
				if(cycles > BURNIN && cycles < BURNIN + SEMCYCLES)
					cat("Stage 2: Cycle = ", cycles-BURNIN+1, ", Log-Lik = ",
						sprintf("%.1f",attr(theta0,"log.lik")), sep="")
				if(cycles > BURNIN + SEMCYCLES)
					cat("Stage 3: Cycle = ", cycles-BURNIN-SEMCYCLES+1, 
						", Log-Lik = ", sprintf("%.1f",attr(theta0,"log.lik")), sep="")					
			}
		}			
		if(stagecycle < 3){	
			ave.h <- as(ave.h,'sparseMatrix')
			inv.ave.h <- try(solve(ave.h))			
			if(class(inv.ave.h) == 'try-error'){
				inv.ave.h <- try(qr.solve(ave.h + 2*diag(ncol(ave.h))))
				noninvcount <- noninvcount + 1
				if(noninvcount == 3) 
					stop('\nEstimation halted during burn in stages, solution is unstable')
			}
			correction <-  inv.ave.h %*% grad	
			correction[correction > 1] <- 1
			correction[correction < -1] <- -1			
			parsold <- pars
			correct <- rep(0,npars)
			correct[sind] <- as.vector(correction)
			correct[constvalues[,1] == 1] <- 0
			if(length(equalconstr) > 0)	
				for(i in 1:length(equalconstr))
					correct[equalconstr[[i]]] <- mean(correct[equalconstr[[i]]])			
			correct[correct[guessind] > .05] <- .05		
			correct[correct[guessind] < -.05] <- -.05
			correct[correct[upperind] > .05] <- .05    	
			correct[correct[upperind] < -.05] <- -.05
			pars <- pars + gamma*correct
			if(verbose && (cycles + 1) %% 10 == 0){ 
				cat(", Max Change =", sprintf("%.4f",max(abs(gamma*correction))), "\n")
				flush.console()
			}			
			pars[covind][pars[covind] > .95] <- parsold[covind][pars[covind] > .95]
			pars[covind][pars[covind] < -.95] <- parsold[covind][pars[covind] < -.95]
			pars[guessind][pars[guessind] < 0] <- parsold[guessind][pars[guessind] < 0]
			pars[upperind][pars[upperind] > 1] <- parsold[upperind][pars[upperind] > 1]
			if(stagecycle == 2){
				SEM.stores[cycles - BURNIN,] <- pars
				SEM.stores2[[cycles - BURNIN]] <- ave.h
			}	
			next
		}	 
		
		#Step 3. Update R-M step		
		Tau <- Tau + gamma*(ave.h - Tau)
		Tau <- as(Tau,'sparseMatrix')	
		inv.Tau <- solve(Tau)
		if(class(inv.Tau) == 'try-error'){
			inv.Tau <- try(qr.solve(Tau + 2 * diag(ncol(Tau))))
			noninvcount <- noninvcount + 1
			if(noninvcount == 3) 
				stop('\nEstimation halted during stage 3, solution is unstable')
		}		
		correction <-  inv.Tau %*% grad
		parsold <- pars
		correct <- rep(0,npars)
		correct[sind] <- as.vector(correction)
		correct[constvalues[,1] == 1] <- 0
		if(length(equalconstr) > 0)		
			for(i in 1:length(equalconstr))
				correct[equalconstr[[i]]] <- mean(correct[equalconstr[[i]]])	
		if(verbose && (cycles + 1) %% 10 == 0){ 
			cat(", gam = ",sprintf("%.3f",gamma),", Max Change = ", 
				sprintf("%.4f",max(abs(gamma*correction))), "\n", sep = '')
			flush.console()		
		}	
		if(all(gamma*correct < TOL)) conv <- conv + 1
			else conv <- 0		
		if(conv == 3) break
		correct[correct[guessind] > .025] <- .025		
		correct[correct[guessind] < -.025] <- -.025	
		correct[correct[upperind] > .025] <- .025    	
		correct[correct[upperind] < -.025] <- -.025
		pars <- pars + gamma*correct	
		pars[covind][pars[covind] > .95] <- parsold[covind][pars[covind] > .95]
		pars[covind][pars[covind] < -.95] <- parsold[covind][pars[covind] < -.95]
		pars[guessind][pars[guessind] < 0] <- parsold[guessind][pars[guessind] < 0]
		pars[upperind][pars[upperind] > 1] <- parsold[upperind][pars[upperind] > 1]
		
		#Extra: Approximate information matrix.	sqrt(diag(solve(info))) == SE
		if(gamma == .25) gamma <- 1	
		phi <- phi + gamma*(grad - phi)
		info <- info + gamma*(Tau - phi %*% t(phi) - info)		
	} ###END BIG LOOP	
	    
	normpars <- sortParsConfmirt(pars, indlist, nfact, estGuess, nfactNames)
	if(verbose) cat("\n\n")
	SEtmp <- diag(solve(info))		
	if(any(SEtmp < 0)){
		warning("Information matrix is not positive definite, negative SEs set to 'NA'.\n")
		SEtmp[SEtmp < 0] <- NA
	}	
	SEtmp <- sqrt(SEtmp)	
	SE <- rep(NA,npars) 
	SE[parind[sind]] <- SEtmp
	SE[constvalues[ ,1]==1] <- NA
	if(length(equalconstr) > 0)
		for(i in 1:length(equalconstr))
			SE[equalconstr[[i]]] <- mean(SE[equalconstr[[i]]])
	estpars <- pars[sind]
	lambdas <- matrix(pars[lamind],J,nfactNames,byrow=TRUE)	
	lambdas[!estlam & !lambdas != 0] <- NA	
	guess <- rep(NA,J)
	guess <- pars[guessind]
	guess[!estGuess] <- NA
	guess[K == 2 & !estGuess] <- 0
	upper <- rep(NA,J)
	upper <- pars[upperind]
	upper[!estUpper] <- NA
	upper[K == 2 & !estUpper] <- 1
	zetas <- pars[indlist$zetaind]
	u <- pars[meanind]	
	sig <- matrix(0,nfact,nfact)
	SElam <- matrix(SE[lamind],J,nfactNames,byrow=TRUE)
	SEzetas <- SE[indlist$zetaind]	
	SEg <- rep(NA,J)	
	SEg <- SE[guessind]	
	SEup <- SE[upperind]
	SEg[!estGuess] <- NA
	SEup[!estUpper] <- NA
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
		if((max(K)-1) > nfactNames) tmp1 <- tmp2 <- matrix(NA,J,(max(K)-1))
		else tmp1 <- tmp2 <- matrix(NA,J,nfactNames)
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
			for(j in 1:nfactNames){
				tmp1[i,j] <- zetas[loc]
				tmp2[i,j] <- SEzetas[loc]
				loc <- loc + 1
			}	
		}	
	}	 
	zetas <- tmp1
	SEzetas <- tmp2	
	parsprint <- cbind(lambdas,zetas)
	SEpars <- cbind(SElam,SEzetas)
	gpars <- list(u = u, sig = sig)	
	SEgpars <- list(SEu = SEu, SEsig = SEsig)
	estpars <- list(estlam=estlam,estGuess=estGuess,estUpper=estUpper, estgcov=estgcov,
		estgmeans=estgmeans,estComp=estComp)		
		
	if (nfactNames > 1) norm <- sqrt(1 + rowSums(lambdas[ ,1:nfactNames]^2,na.rm = TRUE))
		else norm <- as.matrix(sqrt(1 + lambdas[ ,1]^2))  
	F <- as.matrix(lambdas[ ,1:nfactNames]/norm)
	F[is.na(F)] <- 0		
	h2 <- rowSums(F^2)
	colnames(F) <- factorNames
	names(h2) <- itemnames    
	null.mod <- unclass(mirt(data, 0))
    
    if(exploratory){
        mod <- new('polymirtClass',pars=normpars, guess=guess, SEpars=SEpars, SEg=SEg,
                   upper=upper, SEup=SEup, cycles=cycles-SEMCYCLES-BURNIN, Theta=theta0, 
                   fulldata=fulldata, data=data, K=K, F=F, h2=h2, itemloc=itemloc, 
                   converge=converge, estGuess=estGuess, rotate=rotate, null.mod=null.mod, 
                   Target=Target, Call=Call)
    } else {
    	mod <- new('confmirtClass', pars=normpars, parsprint=parsprint, guess=guess, upper=upper, 
                   SEg=SEg, SEup=SEup, gpars=gpars, SEgpars=SEgpars, estpars=estpars,K=K, 
                   itemloc=itemloc, cycles=cycles - SEMCYCLES - BURNIN, Theta=theta0, 
                   fulldata=fulldata, data=data, h2=h2,F=F,converge=converge, 
                   nconstvalues=as.integer(nconstvalues), SEpars=SEpars, estComp=estComp, 
                   prodlist=as.list(prodlist), null.mod=null.mod, Call=Call)
    }
	if(calcLL){
		if(verbose) cat("Calculating log-likelihood...\n")
		flush.console()
		mod <- logLik(mod,draws)		        
	}	
	return(mod)
}
