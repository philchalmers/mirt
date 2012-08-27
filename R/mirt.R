#' Full-Information Item Factor Analysis (Multidimensional Item Response
#' Theory)
#' 
#' \code{mirt} fits an unconditional maximum likelihood factor analysis model
#' to dichotomous and polytomous data under the item response theory paradigm. 
#' Fits univariate and multivariate Rasch, 1-4PL, graded, (generalized) partial credit, 
#' nominal, multiple choice, and partially compenatory models using the EM algorithm.
#' 
#' \code{mirt} follows the item factor analysis strategy by marginal maximum
#' likelihood estimation (MML) outlined in Bock and Aiken (1981), Bock,
#' Gibbons and Muraki (1988), and Muraki and Carlson (1995). 
#' Nested models may be compared via the approximate
#' chi-squared difference test or by a reduction in AIC/BIC values (comparison
#' via \code{\link{anova}}). The general equation used for 
#' multidimensional item response theory is a logistic form with a scaling
#' correction of 1.702. This correction is applied to allow comparison to
#' mainstream programs such as TESTFACT (2003) and POLYFACT. 
#' 
#' Estimation often begins by computing a matrix of quasi-tetrachoric correlations,
#' potentially with Carroll's (1945) adjustment for chance responds. A MINRES
#' factor analysis with \code{nfact} is then extracted and item parameters are
#' estimated by \eqn{a_{ij} = f_{ij}/u_j}, where \eqn{f_{ij}} is the factor
#' loading for the \emph{j}th item on the \emph{i}th factor, and \eqn{u_j} is
#' the square root of the factor uniqueness, \eqn{\sqrt{1 - h_j^2}}. The
#' initial intercept parameters are determined by calculating the inverse
#' normal of the item facility (i.e., item easiness), \eqn{q_j}, to obtain
#' \eqn{d_j = q_j / u_j}. A similar implementation is also used for obtaining 
#' initial values for polychotomous items. Following these initial estimates the model is
#' iterated using the EM estimation strategy with fixed quadrature points.
#' Implicit equation accelerations described by Ramsey (1975) are also added to
#' facilitate parameter convergence speed, and these are adjusted every third
#' cycle.
#' 
#' Factor scores are estimated assuming a normal prior distribution and can be
#' appended to the input data matrix (\code{full.data = TRUE}) or displayed in
#' a summary table for all the unique response patterns. \code{summary} and \code{coef} allow
#' for all the rotations available from the \code{GPArotation} package (e.g., \code{rotate = 'oblimin'})
#' as well as a \code{'promax'} rotation. 
#' 
#' Using \code{plot} will plot the either the test surface function or the test
#' information function for 1 and 2 dimensional solutions. To examine
#' individual item plots use \code{\link{itemplot}}. Residuals are
#' computed using the LD statistic (Chen & Thissen, 1997) in the lower
#' diagonal of the matrix returned by \code{residuals}, and Cramer's V above
#' the diagonal.
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
#' specifications. Rotation and target matrix options will be used in this subroutine and will be
#' passed to the returned object for use in generic functions such as \code{summary()} and 
#' \code{fscores}. Again, factor means and variances are fixed to ensure proper identification. See
#' \code{\link{mirt}} for more details.
#' 
#'
#' 
#' @aliases mirt summary,mirt-method coef,mirt-method anova,mirt-method fitted,mirt-method
#' plot,mirt-method residuals,mirt-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data, with missing data coded as \code{NA}
#' @param model an object returned from \code{confmirt.model()} declaring how
#' the factor model is to be estimated, or a single numeric value indicating the number 
#' of exploratory factors to estimate. See \code{\link{confmirt.model}} for
#' more details
#' @param itemtype type of items to be modeled, declared as a vector for each item or a single value
#' which will be repeated globally. The NULL default assumes that the items are ordinal or 2PL,
#' however they may be changed to the following: 'Rasch', '1PL', '2PL', '3PL', '3PLu', 
#' '4PL', 'graded', 'gpcm', 'nominal', 'mcm', and 'partcomp', for the Rasch/partial credit, 1 and 2 parameter logistic, 
#' 3 parameter logistic (lower asymptote and upper), 4 parameter logistic, graded response model, 
#' generalized partial credit model, nominal model, multiple choice model, and partially compensatory model,
#' respectively. If \code{NULL} the default assumes that the data follow a '2PL' or 'graded' format
#' @param SE logical, estimate the standard errors?
#' @param guess fixed pseudo-guessing parameters. Can be entered as a single
#' value to assign a global guessing parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param upper fixed upper bound parameters for 4-PL model. Can be entered as a single
#' value to assign a global guessing parameter or may be entered as a numeric
#' vector corresponding to each item
#' @param prev.cor use a previously computed correlation matrix to be used to
#' estimate starting values for the EM estimation? Default in \code{NULL} 
#' @param rotate type of rotation to perform after the initial orthogonal
#' parameters have been extracted by using \code{summary}; default is \code{'varimax'}. 
#' See below for list of possible rotations. If \code{rotate != ''} in the \code{summary} 
#' input then the default from the object is ignored and the new rotation from the list 
#' is used instead
#' @param allpars logical; print all the item parameters instead of just the slopes?
#' @param Target a dummy variable matrix indicting a target rotation pattern
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
#' @param quadpts number of quadrature points per dimension
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#' option. Only prints patterns that have standardized residuals greater than 
#' \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param print logical; print output to console?
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
#' @param theta_angle numeric values ranging from 0 to 90 used in \code{plot}. If a vector is 
#' used then a bubble plot is created with the summed information across the angles specified 
#' (e.g., \code{theta_angle = seq(0, 90, by=10)})
#' @param npts number of quadrature points to be used for plotting features.
#' Larger values make plots look smoother
#' @param rot allows rotation of the 3D graphics
#' @param restype type of residuals to be displayed. Can be either \code{'LD'}
#' for a local dependence matrix (Chen & Thissen, 1997) or \code{'exp'} for the
#' expected values for the frequencies of every response pattern
#' @param verbose logical; print observed log-likelihood value at each iteration?
#' @param debug logical; turn on debugging features?
#' @param technical a list containing lower level technical parameters for estimation. May be:
#' \describe{ 
#' \item{MAXQUAD}{maximum number of quadratures; default 10000}
#' \item{MSTEPMAXIT}{number of M-step iterations; default 25}
#' \item{TOL}{EM convergence threshold; default .001}
#' \item{NCYCLES}{maximum number of EM cycles; default 300}
#' \item{NOWARN}{a logical indicating whether dependent packages warnings should be printed; 
#' default \code{TRUE}}
#' }
#' @param ... additional arguments to be passed
#' @section Convergence:
#' 
#' Unrestricted full-information factor analysis is known to have problems with
#' convergence, and some items may need to be constrained or removed entirely
#' to allow for an acceptable solution. As a general rule dichotomous items with
#' means greater than .95, or items that are only .05 greater than the
#' guessing parameter, should be considered for removal from the analysis or
#' treated with prior distributions. The same type of reasoning is applicable when including 
#' upper bound parameters as well. Also, increasing the number of quadrature
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
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6),
#' 1-29.
#'
#' Muraki, E. & Carlson, E. B. (1995). Full-information factor analysis for polytomous 
#' item responses. \emph{Applied Psychological Measurement, 19}, 73-90.
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
#' mirt(data, model, itemtype = NULL, guess = 0, upper = 1, SE = FALSE, startvalues = NULL,
#' constrain = NULL, freepars = NULL,  parprior = NULL, rotate = 'varimax', Target = NULL, 
#' prev.cor = NULL, quadpts = NULL, verbose = FALSE, debug = FALSE, 
#' technical = list(), ...)
#' 
#' \S4method{summary}{mirt}(object, rotate = '', Target = NULL, suppress = 0, digits = 3, 
#' print = TRUE, ...)
#' 
#' \S4method{coef}{mirt}(object, rotate = '', Target = NULL, allpars = FALSE, digits = 3, ...)
#' 
#' \S4method{anova}{mirt}(object, object2, ...)
#' 
#' \S4method{fitted}{mirt}(object, digits = 3, ...)
#' 
#' \S4method{plot}{mirt}(x, type = 'info', npts = 50, theta_angle = 45, 
#' rot = list(x = -70, y = 30, z = 10), ...)
#' 
#' \S4method{residuals}{mirt}(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
#' @export mirt
#' @examples
#' 
#' \dontrun{
#' #load LSAT section 7 data and computmire 1 and 2 factor models
#' data(LSAT7)
#' data <- expand.table(LSAT7)
#' 
#' (mod1 <- mirt(data, 1))
#' summary(mod1)
#' residuals(mod1)
#' plot(mod1) #test information function
#' 
#' #estimated 3PL model for item 5 only
#' (mod1.3PL <- mirt(data, 1, itemtype = c('2PL', '2PL', '2PL', '2PL', '3PL')))
#' coef(mod1.3PL, allpars = TRUE)
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
#' #data from the 'ltm' package in numeric format
#' pmod1 <- mirt(Science, 1)
#' plot(pmod1)
#' summary(pmod1)
#'
#' #Constrain all slopes to be equal
#' #first obtain parameter index
#' mirt(Science,1, constrain = 'index') #note that slopes are numbered 1,5,9,13
#' (pmod1_equalslopes <- mirt(Science, 1, constrain = list(c(1,5,9,13))))
#' coef(pmod1_equalslopes)
#' 
#' pmod2 <- mirt(Science, 2)
#' coef(pmod2)
#' residuals(pmod2)
#' plot(pmod2, theta_angle = seq(0,90, by = 5)) #sum across angles of theta 1
#' itemplot(pmod2, 1)
#' anova(pmod1, pmod2)
#' 
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
#' #with fixed guessing parameters
#' mod1g <- mirt(data, 1, guess = .1)
#' coef(mod1g)
#' mod2g <- mirt(data, 2, guess = .1)
#' coef(mod2g)
#' anova(mod1g, mod2g)
#' summary(mod2g, rotate='promax')
#' }
#' 
mirt <- function(data, model, itemtype = NULL, guess = 0, upper = 1, SE = FALSE, startvalues = NULL,
                 constrain = NULL, freepars = NULL,  parprior = NULL, rotate = 'varimax', Target = NULL, 
                 prev.cor = NULL, quadpts = NULL, verbose = FALSE, debug = FALSE, 
                 technical = list(), ...)
{   
    if(debug == 'Main') browser()
	Call <- match.call()    
    ##technical
    MAXQUAD <- ifelse(is.null(technical$MAXQUAD), 10000, technical$MAXQUAD)
    MSTEPMAXIT <- ifelse(is.null(technical$MSTEPMAXIT), 25, technical$MSTEPMAXIT)
	TOL <- ifelse(is.null(technical$TOL), .001, technical$TOL)
	NCYCLES <- ifelse(is.null(technical$NCYCLES), 300, technical$NCYCLES)
    NOWARN <- ifelse(is.null(technical$NOWARN), TRUE, technical$NOWARN)        	       
    RETURN <- ifelse(any('index' == c(startvalues, freepars, parprior, constrain)), TRUE, FALSE)
    ##              
    Target <- ifelse(is.null(Target), NaN, Target)   
    data <- as.matrix(data)
    PrepList <- PrepData(data=data, model=model, itemtype=itemtype, guess=guess, upper=upper, 
                         startvalues=startvalues, constrain=constrain, freepars=freepars, 
                         parprior=parprior, verbose=verbose, debug=debug, 
                         technical=technical)
    if(RETURN) return(PrepList)
    NULL.MODEL <- ifelse(PrepList$itemtype[1] == 'NullModel' && !any(is.na(data)), TRUE, FALSE)    
    nfact <- PrepList$nfact
	if (is.null(quadpts)) quadpts <- ceiling(40/(PrepList$nfact^1.5))  
	Theta <- theta <- as.matrix(seq(-4,4,length.out = quadpts))
	if(quadpts^nfact <= MAXQUAD){
		Theta <- thetaComb(theta,nfact)		
	} else stop('Greater than ', MAXQUAD, ' quadrature points.')           
    ESTIMATE <- EM(pars=PrepList$pars, NCYCLES=NCYCLES, MSTEPMAXIT=MSTEPMAXIT, TOL=TOL,                    
                   tabdata=PrepList$tabdata, tabdata2=PrepList$tabdata2, npars=PrepList$npars,
                   Theta=Theta, theta=theta, itemloc=PrepList$itemloc, debug=debug, verbose=verbose, 
                   constrain=PrepList$constrain, data=data, NULL.MODEL=NULL.MODEL)	
	# pars to FA loadings
    pars <- ESTIMATE$pars
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
    	mod <- new('ExploratoryClass', iter=ESTIMATE$cycles, pars=ESTIMATE$pars, G2=ESTIMATE$G2, 
    	           df=ESTIMATE$df, p=ESTIMATE$p, itemloc=PrepList$itemloc, AIC=ESTIMATE$AIC, 
    	           BIC=ESTIMATE$BIC, logLik=ESTIMATE$logLik, F=F, h2=h2, tabdata=PrepList$tabdata2, 
    	           Theta=Theta, Pl=ESTIMATE$Pl, data=data, converge=ESTIMATE$converge, nfact=nfact,               
    	           quadpts=quadpts, RMSEA=ESTIMATE$RMSEA, K=PrepList$K, tabdatalong=PrepList$tabdata, 
    	           rotate=rotate, null.mod=ESTIMATE$null.mod, TLI=ESTIMATE$TLI, Target=Target, 
    	           factorNames=PrepList$factorNames, Call=Call)
    } else {
        F <- alp
        colnames(F) <- paste("F_", 1:ncol(F),sep="")    
        h2 <- rowSums(F^2)       
        mod <- new('ConfirmatoryClass', iter=ESTIMATE$cycles, pars=ESTIMATE$pars, G2=ESTIMATE$G2, 
                   df=ESTIMATE$df, p=ESTIMATE$p, itemloc=PrepList$itemloc, AIC=ESTIMATE$AIC, 
                   BIC=ESTIMATE$BIC, logLik=ESTIMATE$logLik, F=F, h2=h2, tabdata=PrepList$tabdata2, 
                   Theta=Theta, Pl=ESTIMATE$Pl, data=data, converge=ESTIMATE$converge, nfact=nfact,               
                   quadpts=quadpts, RMSEA=ESTIMATE$RMSEA, K=PrepList$K, tabdatalong=PrepList$tabdata, 
                   null.mod=ESTIMATE$null.mod, TLI=ESTIMATE$TLI, factorNames=PrepList$factorNames, 
                   Call=Call)
    }
		  
	return(mod)    
}
