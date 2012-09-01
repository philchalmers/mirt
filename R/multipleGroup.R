#' @param invariance a character vector containing the following possible options: 
#' \describe{ 
#' \item{'free_means'}{for freely estimating all latent means (reference group constrained to 0)}
#' \item{'free_varcov'}{for freely estimating the variance-covariance matrix accross groups 
#' (reference group has variances equal to 1, but
#' freely estimated covariance terms if specified in the model)}
#' \item{'covariances'}{to constrain all the covariance parameters to be equal, note that this only 
#' makes sense if the factor variances are the same (i.e., unity)}
#' \item{'slopes'}{to constrain all the slopes to be equal across all groups} 
#' \item{'intercepts'}{to constrain all the intercepts to be equal across all groups, note for 
#' nominal models this also includes the category specific slope parameters}
#'}
multipleGroup <- function(data, model, itemtype = NULL, guess = 0, upper = 1, group = NULL, 
                          invariance = '', constrain = NULL, startvalues = NULL, 
                          parprior = NULL, freepars = NULL, calcLL = TRUE, draws = 2000,
                          technical = list(), debug = FALSE, verbose = TRUE)
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
    if(!is.null(technical$set.seed)) set.seed(technical$set.seed)	
    gain <- c(0.05,0.5,0.004)
    if(!is.null(technical$gain)) {
        if(length(technical$gain) == 3 && is.numeric(technical$gain))
            gain <- technical$gain
    }	 
    RETURNFREEPARS <- RETURNSTARTVALUES <- RETURNPARINDEX <- FALSE
    ##	    
    data <- as.matrix(data)
    group <- factor(group)
    groupNames <- unique(group)
    ngroups <- length(groupNames)    
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
    PrepListFull <- PrepData(data=data, model=model, itemtype=itemtype, guess=guess, upper=upper, 
                             startvalues=NULL, constrain=NULL, freepars=NULL, 
                             parprior=NULL, verbose=verbose, debug=debug, 
                             technical=technical)
    for(g in 1:ngroups){    
        select <- group == groupNames[g]
        PrepList[[g]] <- PrepData(data=data[select,], model=model, itemtype=itemtype, guess=guess, upper=upper, 
                             startvalues=startvalues[g], constrain=constrain, freepars=freepars[g], 
                             parprior=parprior, verbose=verbose, debug=debug, 
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
    #ESTIMATE <- MHRM.group()
    
    
    
    
#     mod <- new('MultipleGroupClass', ...)
#     if(calcLL){
#         if(verbose) cat("\nCalculating log-likelihood...\n")
#         flush.console()
#         mod <- calcLogLik(mod, draws)
#     }
    
}