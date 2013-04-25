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

setMethod(
    f = "plot",
    signature = signature(x = 'MultipleGroupClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45, 
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10), 
                          auto.key = TRUE, ...)
    {           
        if (!type %in% c('info','infocontour', 'SE', 'RE', 'score')) 
            stop(type, " is not a valid plot type.")
        if (any(theta_angle > 90 | theta_angle < 0)) 
            stop('Improper angle specifed. Must be between 0 and 90.')
        if(length(theta_angle) > 1) stop('No info-angle plot is available')
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])       
        ngroups <- length(x@cmods)          
        J <- length(x@cmods[[1]]@pars) - 1
        nfact <- x@nfact
        if(nfact > 2) stop("Can't plot high dimensional solutions.")
        if(nfact == 1) theta_angle <- 0        
        pars <- x@cmods
        theta <- seq(-4,4,length.out=npts)             
        ThetaFull <- Theta <- thetaComb(theta, nfact)        
        prodlist <- attr(x@pars, 'prodlist')
        if(length(prodlist) > 0)        
            ThetaFull <- prodterms(Theta,prodlist)        
        infolist <- vector('list', ngroups)            
        for(g in 1:ngroups){
            info <- 0
            for(i in 1:J){
                tmp <- extract.item(x, i, g)
                info <- info + iteminfo(tmp, Theta=ThetaFull, degrees=theta_angle)
            }
            infolist[[g]] <- info            
        }
        if(type == 'RE') infolist <- lapply(infolist, function(x) x / infolist[[1]])
        info <- do.call(rbind, infolist)
        Theta <- ThetaFull
        for(g in 2:ngroups) Theta <- rbind(Theta, ThetaFull)
        groups <- gl(ngroups, nrow(ThetaFull), labels=x@groupNames)        
        adj <- apply(x@data, 2, min)
        if(any(adj > 0) && type == 'score')
            message('Adjusted so that the lowest category score for every item is 0')
        gscore <- c()
        for(g in 1:ngroups){
            itemtrace <- computeItemtrace(x@cmods[[g]]@pars, ThetaFull, x@itemloc)        
            score <- c()
            for(i in 1:J)
                score <- c(score, 0:(x@K[i]-1))
            score <- matrix(score, nrow(itemtrace), ncol(itemtrace), byrow = TRUE)
            gscore <- c(gscore, rowSums(score * itemtrace))
        }
        plt <- data.frame(info=info, score=gscore, Theta, group=groups)        
        if(nfact == 2){    					
            colnames(plt) <- c("info", "score", "Theta1", "Theta2", "group")			
            plt$SE <- 1 / sqrt(plt$info)            
            if(type == 'infocontour')    											
                return(contourplot(info ~ Theta1 * Theta2|group, data = plt, 
                                   main = paste("Test Information Contour"), xlab = expression(theta[1]), 
                                   ylab = expression(theta[2]), ...))
            if(type == 'info')
                return(wireframe(info ~ Theta1 + Theta2|group, data = plt, main = "Test Information", 
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = TRUE, ...)) 
            if(type == 'RE')
                return(wireframe(info ~ Theta1 + Theta2|group, data = plt, main = "Relative Efficiency", 
                                 zlab=expression(RE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = TRUE, ...))
            if(type == 'SE')                
                return(wireframe(SE ~ Theta1 + Theta2|group, data = plt, main = "Test Standard Errors", 
                                 zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE, 
                                 auto.key = TRUE, ...))    
            if(type == 'score')
                return(wireframe(score ~ Theta1 + Theta2|group, data = plt, main = "Expected Total Score", 
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = TRUE, ...)) 
        } else {            
            colnames(plt) <- c("info", "score", "Theta", "group")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'info')
                return(xyplot(info~Theta, plt, type='l', group=group, main = 'Test Information', 
                              xlab = expression(theta), ylab=expression(I(theta)), auto.key = TRUE, ...))				
            if(type == 'RE')
                return(xyplot(info~Theta, plt, type='l', group=group, main = 'Relative Efficiency', 
                              xlab = expression(theta), ylab=expression(RE(theta)), auto.key = TRUE, ...))    			
            if(type == 'infocontour') 
                cat('No \'contour\' plots for 1-dimensional models\n')
            if(type == 'SE')                
                return(xyplot(SE~Theta, plt, type='l', group=group, main = 'Test Standard Errors', 
                              xlab = expression(theta), ylab=expression(SE(theta)), auto.key = TRUE, ...))
            if(type == 'score')
                return(xyplot(score~Theta, plt, type='l', group=group, main = 'Expected Total Score', 
                              xlab = expression(theta), ylab=expression(Total(theta)), auto.key = TRUE, ...))
        }
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, ...)
    {         
        ret <- vector('list', length(object@groupNames))
        names(ret) <- object@groupNames        
        for(g in 1:length(ret))
            ret[[g]] <- residuals(object@cmods[[g]], verbose = FALSE, ...)
        ret
    }
)

setMethod(
    f = "fitted",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, ...)
    {        
        ret <- vector('list', length(object@groupNames))
        names(ret) <- object@groupNames
        for(g in 1:length(ret))
            ret[[g]] <- fitted(object@cmods[[g]], ...)
        ret        
    }    
)
