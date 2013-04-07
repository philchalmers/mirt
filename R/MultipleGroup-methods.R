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
            cat("AIC = ", x@AIC, "; AICc = ", x@AICc, "\n", sep='')
            cat("BIC = ", x@BIC, "; SABIC = ", x@SABIC, "\n", sep='')            
            if(!is.nan(x@p)){                
                cat("G2 (", x@df,") = ", round(x@G2,2), ", p = ", round(x@p,4), 
                    "\nX2 (", x@df,") = ", round(x@X2,2), ", p = ", round(x@p.X2,4), sep='')                     
                cat("\nRMSEA (G2) = ", round(x@RMSEA,3), "; RMSEA (X2) = ", round(x@RMSEA.X2,3), sep='')
                cat("\nCFI (G2) = ", round(x@CFI,3), "; CFI (X2) = ", round(x@CFI.X2,3), sep='')                    
                cat("\nTLI (G2) = ", round(x@TLI,3), "; TLI (X2) = ", round(x@TLI.X2,3), '\n\n', sep='') 
                for(g in 1:length(x@cmods))
                    cat(as.character(x@groupNames[g]), " group: G2 = ", round(x@cmods[[g]]@G2,2), 
                        ", X2 = ", round(x@cmods[[g]]@X2,2), '\n', sep='') 
            }
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
        groupind <- length(object@cmods[[1]]@pars)
        nfact <- object@nfact
        ret <- list()
        coeflist <- coef(object)        
        for(g in 1:ngroups){
            if(verbose) cat('\n----------\nGROUP:', as.character(object@groupNames[g]), '\n')
            ret[[g]] <- summary(object@cmods[[g]], digits=digits, verbose=verbose, ...)
            if(is(coeflist[[g]][[groupind]], 'matrix'))
                ret[[g]]$mean <- coeflist[[g]][[groupind]][1, 1:nfact]
            else ret[[g]]$mean <- coeflist[[g]][[groupind]][1:nfact]
            names(ret[[g]]$mean) <- colnames(ret[[g]]$fcor)
            if(verbose){
                cat('\nFactor means:\n') 
                print(round(ret[[g]]$mean, digits))
            }
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
        X2 <- round(2*object2@logLik - 2*object@logLik, 3)
        cat('\nModel 1: ')
        print(object@Call)
        cat('Model 2: ')
        print(object2@Call)
        cat('\n')
        ret <- cbind(Df = c(object@df, object2@df),
                          AIC = c(object@AIC, object2@AIC),
                          AICc = c(object@AICc, object2@AICc),
                          BIC = c(object@BIC, object2@BIC), 
                          SABIC = c(object@SABIC, object2@SABIC),
                          logLik = c(object@logLik, object2@logLik),
                          X2 = c(NA, X2),
                          df = c(NA, abs(df)),
                          p = c(NA, round(1 - pchisq(X2,abs(df)),3)))         
        ret
    }		
)

    
