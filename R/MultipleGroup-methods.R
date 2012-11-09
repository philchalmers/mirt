# Methods
setMethod(
    f = "print",
    signature = signature(x = 'MultipleGroupClass'),
    definition = function(x)
    {
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
        cat("Full-information item factor analysis with ", x@nfact, " factors \n", sep="")
        EMquad <- ''
        if(x@method == 'EM') EMquad <- c(' with ', x@quadpts, ' quadrature')
        if(x@converge == 1)    
            cat("Converged in ", x@iter, " iterations", EMquad, ". \n", sep = "")
        else     
            cat("Estimation stopped after ", x@iter, " iterations", EMquad, ". \n", sep="")        
        if(length(x@logLik) > 0){
            cat("Log-likelihood = ", x@logLik, ifelse(length(x@SElogLik) > 0, 
                                                               paste(', SE = ', round(x@SElogLik,3)),
                                                               ''), "\n",sep='')			            
            cat("AIC =", x@AIC, "\n")			
            cat("BIC =", x@BIC, "\n")            	
        }		
        if(!is.nan(x@p)){            
            for(g in 1:length(x@cmods))
                cat('G^2 (', as.character(x@groupNames[g]), ') = ', 
                    round(x@cmods[[g]]@G2,2), '\n', sep = '')  
            cat("G^2 (Total) = ", round(x@G2,2), ", df = ", 
                x@df, ", p = ", round(x@p,4), "\nTLI = ", round(x@TLI,3),
                ", RMSEA = ", round(x@RMSEA,3), "\n", sep="")                                       
        } else { 
            cat("G^2 = ", NA, ", df = ", 
                x@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="")
        }
    } 
)

setMethod(
    f = "show",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object) {
        print(object)
    } 
)

setMethod(
    f = "coef",
    signature = 'MultipleGroupClass',
    definition = function(object, digits = 3, verbose = TRUE, ...)
    {                           
        ngroups <- length(object@cmods)
        allPars <- vector('list', ngroups)
        names(allPars) <- object@groupNames
        itemnames <- colnames(object@data)
        J <- length(itemnames)
        for(g in 1:ngroups){              
            allPars[[g]] <- list()
            if(length(object@cmods[[1]]@pars[[1]]@SEpar) > 0){
                for(i in 1:(J+1)){
                    allPars[[g]][[i]] <- round(matrix(c(object@cmods[[g]]@pars[[i]]@par, 
                                                   object@cmods[[g]]@pars[[i]]@SEpar), 
                                             2, byrow = TRUE), digits)
                    rownames(allPars[[g]][[i]]) <- c('pars', 'SE')
                    colnames(allPars[[g]][[i]]) <- names(object@cmods[[g]]@pars[[i]]@parnum)
                }
            } else {
                for(i in 1:(J+1)){
                    allPars[[g]][[i]] <- round(object@cmods[[g]]@pars[[i]]@par, digits)
                    names(allPars[[g]][[i]]) <- names(object@cmods[[g]]@pars[[i]]@parnum)            
                }
            }                      
            names(allPars[[g]]) <- c(itemnames, 'GroupPars')                
        }
        return(allPars)                  	
    }
)

setMethod(
    f = "summary",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, digits = 3, verbose = TRUE, ...) {        
        ngroups <- length(object@cmods)       
        ret <- list()
        for(g in 1:ngroups){
            if(verbose) cat('\n----------\nGROUP:', as.character(object@groupNames[g]), '\n')
            ret[[g]] <- summary(object@cmods[[g]], digits=digits, verbose=verbose, ...)
        }
        invisible(ret)
    } 
)

setMethod(
    f = "anova",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, object2)
    {        
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
            ", df = ", df, ", p = ", round(1 - pchisq(X2,abs(df)),4), 
            "\n", sep="")
        cat("AIC difference = ", round(AICdiff,3),"\n", sep='')
        cat("BIC difference = ", round(BICdiff,3),"\n", sep='')
    }		
)
