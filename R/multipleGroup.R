#' Multiple Group Estimation
#' 
#' \code{multipleGroup} performes a full-information
#' maximum-likelihood multiple group analysis for dichotomous and polytomous
#' data under the item response theory paradigm using either Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm or with an EM approach. 
#'  
#' By default the estimation in \code{multipleGroup} assumes that the models are maximimally 
#' independent, and therefore could initially be performed by subsetting the data and running identical
#' models with \code{confmirt} or \code{mirt} and aggregating the results (e.g., log-likelihood). 
#' However, constrains may be imposed accross groups by invoking various \code{invariance} keywords
#' or by inputing user defined \code{freepars}, \code{constrain}, and \code{startvalues} lists.  
#' 
#' @aliases multipleGroup coef,MultipleGroupClass-method 
#' anova,MultipleGroupClass-method 
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data, with missing data coded as \code{NA}
#' @param model an object or named list of objects returned from \code{confmirt.model()} declaring how
#' the factor model is to be estimated. The names of the list input must correspond to the unique values 
#' in the \code{group} variable. See \code{\link{confmirt.model}} for more details
#' @param group a character vector indicating group membership
#' @param invariance a character vector containing the following possible options: 
#' \describe{ 
#' \item{\code{'free_means'}}{for freely estimating all latent means (reference group constrained to 0)}
#' \item{\code{'free_varcov'}}{for freely estimating the variance-covariance matrix accross groups 
#' (reference group has variances equal to 1, but
#' freely estimated covariance terms if specified in the model)}
#' \item{\code{'covariances'}}{to constrain all the covariance parameters to be equal, note that this only 
#' makes sense if the factor variances are the same (i.e., unity)}
#' \item{\code{'slopes'}}{to constrain all the slopes to be equal across all groups} 
#' \item{\code{'intercepts'}}{to constrain all the intercepts to be equal across all groups, note for 
#' nominal models this also includes the category specific slope parameters}
#'}
#' @param guess initial (or fixed) values for the pseudo-guessing parameter. Can be 
#' entered as a single value to assign a global guessing parameter or may be entered as
#' a numeric vector for each item
#' @param upper initial (or fixed) upper bound parameters for 4-PL model. Can be 
#' entered as a single value to assign a global upper bound parameter or may be entered as a 
#' numeric vector corresponding to each item
#' @param free.start a list containing the start value and logical indicating whether a given parameter 
#' is to be freely estimated. Each element of the list consists of three components, the parameter
#' number, the starting (or fixed) value, and a logical to indicate whether the parameter is free. For
#' example, \code{free.start = list(c(20,0.0,FALSE), c(10,1.5,TRUE))} would fix parameter 20 to 0.0 while
#' parameter 10 would be freely estimated from a starting value of 1.5. Note that this will override 
#' the values specified by a user defined \code{startvalues} or \code{freepars} input for the specified
#' parameters, and this may conflict with the \code{invariance} input (e.g., freeing slopes manually
#' while specifying \code{invariance = 'slopes'} is ambiguous and should be avoided). 
#' @param verbose logical; display iteration history during estimation?
#' @param draws the number of Monte Carlo draws to estimate the log-likelihood
#' @param quadpts the number of quadratures to be used per dimensions when \code{method = 'EM'}
#' @param method a character indicating whether to use the EM (\code{'EM'}) or the MH-RM 
#' (\code{'MHRM'}) algorithm
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
#' @param object an object of class \code{confmirtClass}
#' @param object2 an object of class \code{confmirtClass}
#' @param digits the number of significant digits to be rounded
#' @param technical list specifying subtle parameters that can be adjusted. These 
#' values are 
#' \describe{
#' \item{NCYCLES}{max number of cycles; default 2000 for MHRM and 300 for EM}
#' \item{MAXQUAD}{maximum number of quadratures; default 10000}
#' \item{MSTEPMAXIT}{number of M-step iterations; default 25}
#' \item{BURNIN}{number of burn in cycles (stage 1); default 150}
#' \item{SEMCYCLES}{number of SEM cycles (stage 2); default 50}
#' \item{KDRAWS}{number of parallel MH sets to be drawn; default 1}
#' \item{TOL}{minimum threshold tolerance for convergence of MH-RM, must occur on three consecutive
#' occations; default .001} 
#'   \item{set.seed}{seed number used during estimation. Default is 12345}       
#'   \item{gain}{a vector of three values specifying the numerator, exponent, and subtracted
#'      values for the RM gain value. Default is \code{c(0.05,0.5,0.004)}}   	
#'  \item{return_newconstrain}{if \code{TRUE} returns a list consisting of the constraints to be used
#'  just before estimation begins} 
#' }
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{simdata}},
#' \code{\link{confmirt.model}}
#' @keywords models
#' @usage 
#' multipleGroup(data, model, group, itemtype = NULL, guess = 0, upper = 1, free.start = NULL, 
#' invariance = '', method = 'MHRM', constrain = NULL, startvalues = NULL, 
#' parprior = NULL, freepars = NULL, draws = 2000, quadpts = NULL,
#' technical = list(), debug = FALSE, verbose = TRUE)
#' 
#' \S4method{coef}{MultipleGroupClass}(object, digits = 3, verbose = TRUE, ...)
#' 
#' \S4method{anova}{MultipleGroupClass}(object, object2)
#'
#' @export multipleGroup
#' @examples
#' \dontrun{
#' #single factor
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)    
#' itemtype <- rep('dich', nrow(a))
#' N <- 1000    
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))    
#' models <- confmirt.model()
#'    F1 = 1-15
#' 
#' 
#' mod_configural <- multipleGroup(dat, models, group = group, method = 'EM') #completely seperate analyses
#' mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), method = 'EM') #equal slopes
#' mod_scalar2 <- multipleGroup(dat, models, group = group, method = 'EM',  #equal intercepts, free variance and means
#'                              invariance=c('slopes', 'intercepts', 'free_varcov','free_means'))
#' mod_scalar1 <- multipleGroup(dat, models, group = group, method = 'EM', #fixed means
#'                              invariance=c('slopes', 'intercepts', 'free_varcov'))    
#' mod_fullconstrain <- mirt(data, models) #fix variance (equivelent to full constrain)
#' 
#' anova(mod_metric, mod_configural) #equal slopes only
#' anova(mod_scalar2, mod_metric) #equal intercepts, free variance and mean
#' anova(mod_scalar1, mod_scalar2) #fix mean 
#' anova(mod_fullconstrain, mod_scalar1) #fix variance
#'
#' 
#' #Wald test can be useful here too
#' #compare whether intercepts should be equal
#' index <- multipleGroup(dat, models, group = group, constrain = 'index') 
#' index
#' nitems <- ncol(dat)
#' L <- matrix(0, nitems, 124)
#' for(i in 1:nitems){
#'      L[i, index[[1]][[i]][2]] <- 1    
#'      L[i, index[[2]][[i]][2]] <- -1
#'  } 
#' wald(L, mod_configural)
#' 
#' 
#' #############
#' #multiple factors 
#' 
#' a <- matrix(c(abs(rnorm(5,1,.3)), rep(0,15),abs(rnorm(5,1,.3)),
#' rep(0,15),abs(rnorm(5,1,.3))), 15, 3)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' mu <- c(-.4, -.7, .1)
#' sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)   
#' itemtype <- rep('dich', nrow(a))
#' N <- 1000    
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = mu, sigma = sigma)
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N)) 
#'    
#' #group models
#' model1 <- confmirt.model()
#'    F1 = 1-5
#'    F2 = 6-10
#'    F3 = 11-15    
#' 
#' 
#' model2 <- confmirt.model()
#'    F1 = 1-5
#'    F2 = 6-10
#'    F3 = 11-15
#'    COV = F1*F2, F1*F3, F2*F3
#' 
#' 
#' models <- list(D1=model1, D2=model2) #note the names match the groups
#' 
#' mod_configural <- multipleGroup(dat, models, group = group) #completely seperate analyses
#' mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes')) #equal slopes
#' mod_scalar <- multipleGroup(dat, models, group = group, #equal means, slopes, intercepts
#'                              invariance=c('slopes', 'intercepts', 'free_varcov'))    
#' mod_fullconstrain <- confmirt(data, models)
#' 
#' anova(mod_metric, mod_configural)
#' anova(mod_scalar, mod_metric)
#' anova(mod_fullconstrain, mod_scalar)

#' }
multipleGroup <- function(data, model, group, itemtype = NULL, guess = 0, upper = 1, free.start = NULL, 
                          invariance = '', method = 'MHRM', constrain = NULL, 
                          startvalues = NULL, parprior = NULL, freepars = NULL, draws = 2000, 
                          quadpts = NULL,
                          technical = list(), debug = FALSE, verbose = TRUE)
{
    if(debug == 'Main') browser()
    ##technical
    Call <- match.call()               
    set.seed(12345)	   
    MAXQUAD <- ifelse(is.null(technical$MAXQUAD), 10000, technical$MAXQUAD)
    MSTEPMAXIT <- ifelse(is.null(technical$MSTEPMAXIT), 25, technical$MSTEPMAXIT)        
    RETURN <- ifelse(any('index' == c(startvalues, freepars, parprior, constrain)), TRUE, FALSE)
    NCYCLES <- ifelse(is.null(technical$NCYCLES), 2000, technical$NCYCLES)
    if(method == 'EM')
        NCYCLES <- ifelse(is.null(technical$NCYCLES), 300, technical$NCYCLES)
    BURNIN <- ifelse(is.null(technical$BURNIN), 150, technical$BURNIN)
    SEMCYCLES <- ifelse(is.null(technical$SEMCYCLES), 50, technical$SEMCYCLES)
    KDRAWS  <- ifelse(is.null(technical$KDRAWS), 1, technical$KDRAWS)
    TOL <- ifelse(is.null(technical$TOL), .001, technical$TOL)      
    if(!is.null(technical$set.seed)) set.seed(technical$set.seed)	
    gain <- c(0.05,0.5,0.004)
    if(!is.null(technical$gain)){
        if(length(technical$gain) == 3 && is.numeric(technical$gain))
            gain <- technical$gain
    }	 
    RETURNFREEPARS <- RETURNSTARTVALUES <- RETURNPARINDEX <- FALSE
    NULL.MODEL <- ifelse(!is.null(itemtype) && itemtype[1] == 'NullModel', TRUE, FALSE)
    USEEM <- ifelse(method == 'EM', TRUE, FALSE)
    ##	            
    data <- as.matrix(data)
    rownames(data) <- 1:nrow(data)
    group <- factor(group)
    groupNames <- unique(group)
    ngroups <- length(groupNames)
    oldmodel <- model
    if(length(model) == 1){
        model <- list()
        for(g in 1:ngroups)
            model[[g]] <- oldmodel
        names(model) <- groupNames
    } 
    parnumber <- 1
    PrepList <- vector('list', ngroups)
    if(!is.null(freepars) && freepars == 'index'){
        freepars <- rep('index', ngroups)    
        RETURNFREEPARS <- TRUE
    }
    if(!is.null(startvalues) && startvalues == 'index'){
        startvalues <- rep('index', ngroups)
        RETURNSTARTVALUES <- TRUE
    }
    if((!is.null(constrain) && constrain == 'index') || (!is.null(parprior) && parprior =='index'))
        RETURNPARINDEX <- TRUE 
    PrepListFull <- PrepData(data=data, model=model[[1]], itemtype=itemtype, guess=guess, upper=upper, 
                             startvalues=NULL, constrain=NULL, freepars=NULL, 
                             parprior=NULL, verbose=verbose, debug=debug, free.start=NULL,
                             technical=technical) #just a dummy model to collect fulldata stuff
    for(g in 1:ngroups){    
        select <- group == groupNames[g]        
        tmp <- 1:ngroups
        selectmod <- model[[tmp[names(model) == groupNames[g]]]]
        PrepList[[g]] <- PrepData(data=data[select,], model=selectmod, itemtype=itemtype, guess=guess, upper=upper, 
                             startvalues=startvalues[g], constrain=constrain, freepars=freepars[g], 
                             parprior=parprior, verbose=verbose, debug=debug, free.start=free.start,
                             technical=technical, parnumber=parnumber)
        if(RETURNFREEPARS || RETURNSTARTVALUES) next
        if(RETURNPARINDEX){
            tmp <- PrepList[[g]][[length(PrepList[[g]])]]
            parnumber <- tmp[length(tmp)] + 1
            next            
        }
        tmp <- PrepList[[g]]$pars[[length(PrepList[[g]]$pars)]]
        parnumber <- tmp@parnum[length(tmp@parnum)] + 1
    }    
    if(RETURN){
        names(PrepList) <- groupNames
        return(PrepList)
    }      
    pars <- vector('list', ngroups)
    for(g in 1:ngroups)
        pars[[g]] <- PrepList[[g]]$pars
    J <- length(PrepList[[1]]$itemtype)
    nfact <- PrepList[[1]]$pars[[J+1]]@nfact
    nLambdas <- PrepList[[1]]$pars[[1]]@nfact
    if(is.null(constrain)) constrain <- list()         
    #default MG uses configural model (independent groups but each identified)        
    if('free_means' %in% invariance ){ #Free factor means (means 0 for ref)
        for(g in 2:ngroups)
            pars[[g]][[J + 1]]@est[1:nfact] <- TRUE             
    }    
    if('free_varcov' %in% invariance){ #Free factor vars and covs (vars 1 for ref)        
        for(g in 2:ngroups)
            pars[[g]][[J + 1]]@est[(nfact+1):length(pars[[g]][[J + 1]]@est)] <- TRUE            
    } 
    constrain <- UpdateConstrain(pars=pars, constrain=constrain, invariance=invariance, nfact=nfact, 
                                 nLambdas=nLambdas, J=J, ngroups=ngroups)    
    if(!is.null(technical$return_newconstrain)) return(constrain)    
    startlongpars <- c()
    if(method == 'EM'){
        esttype <- 'EM'
        if(method == 'EM' && nLambdas > nfact) 
            stop('Polynominals and product terms not supported for EM method')
        if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))
        Theta <- theta <- as.matrix(seq(-4,4,length.out = quadpts))
        if(quadpts^nfact <= MAXQUAD){
            Theta <- thetaComb(theta,nfact)    	
        } else stop('Greater than ', MAXQUAD, ' quadrature points.')
        ESTIMATE <- EM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=NCYCLES, TOL=TOL, MSTEPMAXIT=MSTEPMAXIT,
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc,
                                           nfact=nfact, constrain=constrain, verbose=verbose), 
                               Theta=Theta, debug=debug)
        startlongpars <- ESTIMATE$longpars
        logLik <- ESTIMATE$logLik
        SElogLik <- 0                
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                      list = list(NCYCLES=NCYCLES, BURNIN=1, SEMCYCLES=5,
                                  KDRAWS=KDRAWS, TOL=.01, USEEM=USEEM, gain=gain, 
                                  nfactNames=PrepList[[1]]$nfactNames, 
                                  itemloc=PrepList[[1]]$itemloc,  
                                  nfact=nfact, constrain=constrain, verbose=FALSE,
                                  startlongpars=startlongpars), 
                      debug=debug)                 
    } else if(method == 'MHRM'){
        esttype <- 'MHRM'
        ESTIMATE <- MHRM.group(pars=pars, constrain=constrain, PrepList=PrepList,
                               list = list(NCYCLES=NCYCLES, BURNIN=BURNIN, SEMCYCLES=SEMCYCLES,
                                           KDRAWS=KDRAWS, TOL=TOL, USEEM=FALSE, gain=gain, 
                                           nfactNames=PrepList[[1]]$nfactNames, 
                                           itemloc=PrepList[[1]]$itemloc,  
                                           nfact=nfact, constrain=constrain, verbose=verbose,
                                           startlongpars=startlongpars), 
                               debug=debug)        
    }   
    cmods <- list()
    for(g in 1:ngroups){
        cmods[[g]] <- new('ConfirmatoryClass', pars=ESTIMATE$pars[[g]], itemloc=PrepList[[g]]$itemloc, 
                          tabdata=PrepList[[g]]$tabdata2, data=data[group == groupNames[[g]], ], 
                          converge=ESTIMATE$converge, esttype='MHRM',                
                          K=PrepList[[g]]$K, tabdatalong=PrepList[[g]]$tabdata, nfact=nfact, 
                          constrain=constrain,
                          fulldata=PrepList[[g]]$fulldata, factorNames=PrepList[[g]]$factorNames)        
    }            
    if(method =='MHRM'){
        if(verbose) cat("\nCalculating log-likelihood...\n")
        flush.console()      
        LLs <- matrix(0, ngroups, 2)
        for(g in 1:ngroups)
            LLs[g, ] <- calcLogLik(cmods[[g]], draws, G2 = 'return')[1:2]       
        LL <- colSums(LLs)  
        logLik <- LL[1]
        SElogLik <- LL[2]
    } 
    r <- PrepListFull$tabdata
    r <- r[, ncol(r)]
    N <- sum(r)
    logN <- 0
    logr <- rep(0,length(r))
    for (i in 1:N) logN <- logN + log(i)
    for (i in 1:length(r)) 
        for (j in 1:r[i]) 
            logr[i] <- logr[i] + log(j)    		
    if(sum(logr) != 0)		
        logLik <- logLik + logN/sum(logr)							
    nestpars <- nconstr <- 0
    for(g in 1:ngroups)
        for(i in 1:(J+1))
            nestpars <- nestpars + sum(pars[[g]][[i]]@est)
    if(length(constrain) > 0)
        for(i in 1:length(constrain))
            nconstr <- nconstr + length(constrain[[i]]) - 1     
    nmissingtabdata <- sum(is.na(rowSums(PrepListFull$tabdata2)))
    df <- length(r) - nestpars + nconstr + nfact*(nfact - 1)/2 - 1 - nmissingtabdata	
    AIC <- (-2) * logLik + 2 * (length(r) - df - 1)
    BIC <- (-2) * logLik + (length(r) - df - 1)*log(N)    			    	    
    mod <- new('MultipleGroupClass', iter=ESTIMATE$cycles, cmods=cmods, itemloc=PrepListFull$itemloc, 
               tabdata=PrepListFull$tabdata2, data=data, converge=ESTIMATE$converge, esttype=esttype,                
               K=PrepListFull$K, tabdatalong=PrepListFull$tabdata, constrain=constrain,               
               group=group, groupNames=groupNames, invariance=invariance, df=as.integer(df),
               logLik=logLik, SElogLik=SElogLik, AIC=AIC, BIC=BIC, information=ESTIMATE$info, 
               Call=Call)  
    return(mod)        
}
