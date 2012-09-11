#' Confirmatory Full-Information Item Factor Analysis
#' 
#' \code{confmirt} fits a conditional (i.e., confirmatory) full-information
#' maximum-likelihood factor analysis model to dichotomous and polytomous
#' data under the item response theory paradigm using Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm. Fits univariate and multivariate Rasch, 
#' 1-4PL, graded, (generalized) partial credit, nominal, multiple choice, and partially-compensatory models, 
#' potentially with polynomial and product constructed latent traits.
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
#' log-likelihood more accurately see \code{\link{calcLogLik}}.
#' 
#' \code{coef} displays the item parameters with their associated standard
#' errors, while use of \code{summary} transforms the slopes into a factor
#' loadings metric and if the model is exploratory allows for rotating the parameters. 
#' Also, nested models may be compared by using the
#' \code{anova} function, where a Chi-squared difference test and AIC/BIC
#' difference values are displayed.
#' 
#' @section Confirmatory IRT:
#' 
#' Specification of the confirmatory item factor analysis model follows many of
#' the rules in the SEM framework for confirmatory factor analysis. The
#' variances of the latent factors are automatically fixed to 1 to help
#' facilitate model identification. All parameters may be fixed to constant
#' values or set equal to other parameters using the appropriate declarations.
#' 
#' @section Exploratory IRT:
#' 
#' Specifying a number as the second input to confmirt an exploratory IRT model is estimated and 
#' can be viewed as a stochastic analogue of \code{mirt}, with much of the same behaviour and 
#' specifications. 
#' Rotation and target matrix options will be used in this subroutine and will be
#' passed to the returned object for use in generic functions such as \code{summary()} and 
#' \code{fscores}. Again, factor means and variances are fixed to ensure proper identification. See
#' \code{\link{mirt}} for more details.
#' 
#' 
#' @aliases confmirt coef,ConfirmatoryClass-method summary,ConfirmatoryClass-method
#' residuals,ConfirmatoryClass-method anova,ConfirmatoryClass-method fitted,ConfirmatoryClass-method
#' plot,ConfirmatoryClass-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data, with missing data coded as \code{NA}
#' @param model an object returned from \code{confmirt.model()} declaring how
#' the factor model is to be estimated, or a single numeric value indicating the number 
#' of exploratory factors to estimate. See \code{\link{confmirt.model}} for
#' more details
#' @param guess initial (or fixed) values for the pseudo-guessing parameter. Can be 
#' entered as a single value to assign a global guessing parameter or may be entered as
#' a numeric vector for each item
#' @param upper initial (or fixed) upper bound parameters for 4-PL model. Can be 
#' entered as a single value to assign a global upper bound parameter or may be entered as a 
#' numeric vector corresponding to each item
#' @param free.start a list containing the start value and logical indicating whether a given parameter 
#' is to be freely estimated. Each element of the list consists of three components, the parameter
#' number, the starting (or fixed) value, and a logical to indicate whether the parameter is free. For
#' example, \code{free.start = list(c(20,0,FALSE), c(10,1.5,TRUE))} would fix parameter 20 to 0 while
#' parameter 10 would be freely estimated with a starting value of 1.5. Note that this will override 
#' the values specified by a user defined \code{startvalues} or \code{freepars} input for the specified
#' parameters
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#' option. Only prints patterns that have standardized residuals greater than 
#' \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param verbose logical; display iteration history during estimation?
#' @param calcLL logical; calculate the log-likelihood via Monte Carlo
#' integration?
#' @param draws the number of Monte Carlo draws to estimate the log-likelihood
#' @param allpars logical; print all the item parameters instead of just the slopes?
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param itemtype type of items to be modeled, declared as a vector for each item or a single value
#' which will be repeated globally. The NULL default assumes that the items are ordinal or 2PL,
#' however they may be changed to the following: 'Rasch', '1PL', '2PL', '3PL', '3PLu', 
#' '4PL', 'graded', 'gpcm', 'nominal', 'mcm', and 'partcomp', for the Rasch/partial credit, 1 and 2 parameter logistic, 
#' 3 parameter logistic (lower asymptote and upper), 4 parameter logistic, graded response model, 
#' generalized partial credit model, nominal model, multiple choice model, and partially compensatory model,
#' respectively. The default assumes that items follow a '2PL' or 'graded' format
#' If \code{NULL} the default assumes that the data follow a '2PL' or 'graded' format
#' @param constrain a list of user declared equality constraints. To see how to define the
#' parameters correctly use \code{constrain = 'index'} initially to see how the parameters are labeled.
#' To constrain parameters to be equal create a list with separate concatenated vectors signifying which
#' parameters to constrain. For example, to set parameters 1 and 5 equal, and also set parameters 2, 6, and 10 equal
#' use \code{constrain = list(c(1,5), c(2,6,10))}
#' @param parprior a list of user declared prior item probabilities. To see how to define the
#' parameters correctly use \code{parprior = 'index'} initially to see how the parameters are labeled.
#' Can define either normal (normally for slopes and intercepts) or beta (for guessing and upper bounds) prior
#' probabilities. Note that for upper bounds the value used in the prior is 1 - u so that the lower and upper 
#' bounds can function the same. To specify a prior the form is c('priortype', ...), where normal priors 
#' are \code{parprior = list(c(parnumber, 'norm', mean, sd))} and betas are 
#' \code{parprior = list(c(parnumber, 'beta', alpha, beta))}. 
#' @param freepars a list of user declared logical values indicating which parameters to estimate. 
#' To see how to define the parameters correctly use \code{freepars = 'index'} initially to see how the parameters
#' are labeled. These values may be modified and input back into the function by using 
#' \code{freepars=newfreepars}. Note that user input values must match what the default structure 
#' would have been
#' @param startvalues a list of user declared start values for parameters. To see how to define the
#' parameters correctly use \code{startvalues = 'index'} initially to see what the defaults would 
#' noramlly be. These values may be modified and input back into the function by using 
#' \code{startavlues=newstartvalues}. Note that user input values must match what the default structure 
#' would have been
#' @param debug logical; turn on debugging features?
#' @param object an object of class \code{ConfirmatoryClass}
#' @param object2 an object of class \code{ConfirmatoryClass}
#' @param digits the number of significant digits to be rounded
#' @param rotate if \code{model} is numeric (indicating an exploratory item FA) then this 
#' rotation is used. Default is \code{'varimax'}
#' @param Target a dummy variable matrix indicting a target rotation pattern
#' @param suppress a numeric value indicating which factor
#' loadings should be suppressed. Typical values are around .3 in most
#' statistical software. Default is 0 for no suppression
#' @param technical list specifying subtle parameters that can be adjusted. These 
#' values are 
#' @param x an object of class \code{mirt} to be plotted or printed
#' @param y an unused variable to be ignored
#' @param type type of plot to view; can be \code{'curve'} for the total test
#' score as a function of two dimensions, or \code{'info'} to show the test
#' information function for two dimensions
#' @param theta_angle numeric values ranging from 0 to 90 used in \code{plot}. If a vector is 
#' used then a bubble plot is created with the summed information across the angles specified 
#' (e.g., \code{theta_angle = seq(0, 90, by=10)})
#' @param npts number of quadrature points to be used for plotting features.
#' Larger values make plots look smoother
#' @param rot allows rotation of the 3D graphics
#' \describe{
#' \item{NCYCLES}{max number of MH-RM cycles; default 2000}
#' \item{BURNIN}{number of burn in cycles (stage 1); default 150}
#' \item{SEMCYCLES}{number of SEM cycles (stage 2); default 50}
#' \item{KDRAWS}{number of parallel MH sets to be drawn; default 1}
#' \item{TOL}{minimum threshold tolerance for convergence of MH-RM, must occur on three consecutive
#' occations; default .001} 
#'   \item{set.seed}{seed number used during estimation. Default is 12345} 	 
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
#' confmirt(data, model, itemtype = NULL, guess = 0, upper = 1, free.start = NULL, startvalues = NULL, 
#' constrain = NULL, freepars = NULL, parprior = NULL, verbose = TRUE, calcLL = TRUE, 
#' draws = 2000, debug = FALSE, rotate = 'varimax', Target = NULL, 
#' technical = list(),  ...)
#' 
#' \S4method{summary}{ConfirmatoryClass}(object, suppress = 0, digits = 3, verbose = TRUE, ...)
#' 
#' \S4method{coef}{ConfirmatoryClass}(object, allpars = FALSE, digits = 3, verbose = TRUE, ...)
#' 
#' \S4method{anova}{ConfirmatoryClass}(object, object2)
#' 
#' \S4method{fitted}{ConfirmatoryClass}(object, digits = 3, ...)
#' 
#' \S4method{plot}{ConfirmatoryClass}(x, y, type = 'info', npts = 50, theta_angle = 45, 
#' rot = list(xaxis = -70, yaxis = 30, zaxis = 10), ...)
#' 
#' \S4method{residuals}{ConfirmatoryClass}(object, restype = 'LD', digits = 3, printvalue = NULL, 
#' verbose = TRUE, ...)
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
#' items <- c(rep('dich',4), rep('graded',3), 'dich')
#' dataset <- simdata(a,d,2000,items,sigma)
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
#' #polynomial/combinations
#' data(SAT12)
#' data <- key2binary(SAT12,
#'                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' 
#' model.quad <- confmirt.model()
#'        F1 = 1-32
#'   (F1*F1) = 1-32       
#'   
#' 
#' model.combo <- confmirt.model()
#'        F1 = 1-16
#'        F2 = 17-32
#'   (F1*F2) = 1-8
#'
#' 
#' (mod.quad <- confmirt(data, model.quad))
#' (mod.combo <- confmirt(data, model.combo))
#' anova(mod.quad, mod.combo)
#' 
#' }
#' 
confmirt <- function(data, model, itemtype = NULL, guess = 0, upper = 1, free.start = NULL,
                     startvalues = NULL, 
                     constrain = NULL, freepars = NULL, parprior = NULL, verbose = TRUE, calcLL = TRUE, 
                     draws = 2000, debug = FALSE, rotate = 'varimax', Target = NULL, 
                     technical = list(),  ...)
{    
    if(debug == 'Main') browser()
    ##technical
	Call <- match.call()               
	set.seed(12345)	    
    RETURN <- ifelse(any('index' == c(startvalues, freepars, parprior, constrain)), TRUE, FALSE)
    NCYCLES <- ifelse(is.null(technical$NCYCLES), 2000, technical$NCYCLES)
    BURNIN <- ifelse(is.null(technical$BURNIN), 150, technical$BURNIN)
    SEMCYCLES <- ifelse(is.null(technical$SEMCYCLES), 50, technical$SEMCYCLES)
    KDRAWS  <- ifelse(is.null(technical$KDRAWS), 1, technical$KDRAWS)
    TOL <- ifelse(is.null(technical$TOL), .001, technical$TOL)  
    EMSE <- ifelse(is.null(technical$EMSE), FALSE, technical$EMSE)
    if(!is.null(technical$set.seed)) set.seed(technical$set.seed)	
    gain <- c(0.05,0.5,0.004)
    if(!is.null(technical$gain)) {
        if(length(technical$gain) == 3 && is.numeric(technical$gain))
            gain <- technical$gain
    }	
    ##	
    Target <- ifelse(is.null(Target), NaN, Target)
    data <- as.matrix(data)
    parnumber <- 1
	PrepList <- PrepData(data=data, model=model, itemtype=itemtype, guess=guess, upper=upper, 
                         startvalues=startvalues, constrain=constrain, freepars=freepars, 
	                     parprior=parprior, verbose=verbose, debug=debug, free.start=free.start,
                         technical=technical, parnumber=parnumber)           
    if(RETURN) return(PrepList)
 	ESTIMATE <- MHRM(pars=PrepList$pars, 
                      list=list(NCYCLES=NCYCLES, BURNIN=BURNIN, SEMCYCLES=SEMCYCLES, 
                                KDRAWS=KDRAWS, TOL=TOL, gain=gain, nfactNames=PrepList$nfactNames, 
                                itemloc=PrepList$itemloc, fulldata=PrepList$fulldata, 
                                nfact=PrepList$nfact, npars=PrepList$npars, 
                                constrain=PrepList$constrain, verbose=verbose), 
                      debug=debug, startvalues=startvalues, EMSE=EMSE) 
    if(EMSE) return(ESTIMATE)
    null.mod <- unclass(mirt(data,1,itemtype='NullModel', SE = FALSE))
    # pars to FA loadings    
    pars <- ESTIMATE$pars    
    nfact <- pars[[1]]@nfact
    lambdas <- Lambdas(pars)
    if (nfact > 1) norm <- sqrt(1 + rowSums(lambdas[ ,1:nfact]^2))
    else norm <- as.matrix(sqrt(1 + lambdas[ ,1]^2))  
    alp <- as.matrix(lambdas[ ,1:nfact]/norm)
    if(PrepList$exploratory){
        FF <- alp %*% t(alp)
        V <- eigen(FF)$vector[ ,1:nfact]
        L <- eigen(FF)$values[1:nfact]
        if (nfact == 1) F <- as.matrix(V * sqrt(L))
        else F <- V %*% sqrt(diag(L))  
        if (sum(F[ ,1] < 0)) F <- (-1) * F 
        colnames(F) <- paste("F_", 1:ncol(F),sep="")    
        h2 <- rowSums(F^2)
        mod <- new('ExploratoryClass', iter=ESTIMATE$cycles, pars=ESTIMATE$pars, itemloc=PrepList$itemloc, 
                   F=F, h2=h2, tabdata=PrepList$tabdata2, data=data, converge=ESTIMATE$converge, esttype='MHRM',                
                   K=PrepList$K, tabdatalong=PrepList$tabdata, nfact=nfact, constrain=PrepList$constrain,
                   rotate=rotate, null.mod=null.mod, Target=Target, factorNames=PrepList$factorNames,
                   fulldata=PrepList$fulldata, information=ESTIMATE$info, longpars=ESTIMATE$longpars, 
                   Call=Call)
    } else {
        F <- alp
        colnames(F) <- PrepList$factorNames
        h2 <- rowSums(F^2)       
        mod <- new('ConfirmatoryClass', iter=ESTIMATE$cycles, pars=ESTIMATE$pars, itemloc=PrepList$itemloc, 
                   F=F, h2=h2, tabdata=PrepList$tabdata2, data=data, converge=ESTIMATE$converge, esttype='MHRM',                
                   K=PrepList$K, tabdatalong=PrepList$tabdata, nfact=nfact, constrain=PrepList$constrain,
                   fulldata=PrepList$fulldata, null.mod=null.mod, factorNames=PrepList$factorNames, 
                   information=ESTIMATE$info, longpars=ESTIMATE$longpars, Call=Call)
    }        
	if(calcLL){
		if(verbose) cat("\nCalculating log-likelihood...\n")
		flush.console()
		mod <- calcLogLik(mod, draws)
	}	
	return(mod)
}
