# Methods
setMethod(
    f = "print",
    signature = signature(x = 'ConfirmatoryClass'),
    definition = function(x)
    {
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
        cat("Full-information item factor analysis with ", x@nfact, " factors \n", sep="")
        if(x@converge == 1)    
            cat("Converged in ", x@iter, " iterations.\n", sep="")
        else 	
            cat("Estimation stopped after ", x@iter, " iterations.\n", sep="")
        
        if(length(x@logLik) > 0){
            cat("Log-likelihood = ", x@logLik, ifelse(length(x@SElogLik) > 0, 
                                                               paste('SE = ', round(x@SElogLik,3)),
                                                               ''), "\n",sep='')			
            cat("AIC =", x@AIC, "\n")			
            cat("BIC =", x@BIC, "\n")
            if(!is.nan(x@p))
                cat("G^2 = ", round(x@G2,2), ", df = ", 
                    x@df, ", p = ", round(x@p,4), "\nTLI = ", round(x@TLI,3),
                    ", RMSEA = ", round(x@RMSEA,3), "\n", sep="")
            else 
                cat("G^2 = ", NA, ", df = ", 
                    x@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="")		
        }		
    } 
)

setMethod(
    f = "show",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object) {
        print(object)
    } 
)

setMethod(
    f = "summary",
    signature = 'ConfirmatoryClass',
    definition = function(object, rotate = '', Target = NULL, suppress = 0, digits = 3, 
                          print = TRUE, ...)
    {           
        nfact <- ncol(object@F)
        itemnames <- names(object@h2)	
        F <- object@F
        rownames(F) <- itemnames								
        SS <- apply(F^2,2,sum)			
        cat("\nFactor loadings metric: \n")
        print(cbind(F),digits)		
        cat("\nSS loadings: ",round(SS,digits), "\n")		
        cat("\nFactor correlations: \n")
        gpars <- ExtractGroupPars(object@pars[[length(object@pars)]])
        Phi <- cov2cor(gpars$gcov)	  
        Phi <- round(Phi, digits)
        colnames(Phi) <- rownames(Phi) <- paste('F',1:ncol(Phi), sep='')
        print(Phi)				
        invisible(F)  	          
    }
)

setMethod(
    f = "coef",
    signature = 'ConfirmatoryClass',
    definition = function(object, rotate = '', Target = NULL, allpars = FALSE, digits = 3, ...)
    {                             
        K <- object@K
        J <- length(K)
        nLambdas <- ncol(object@F)
        a <- matrix(0, J, nLambdas)        
        for(i in 1:J)
            a[i, ] <- ExtractLambdas(object@pars[[i]])        
        A <- sqrt(apply(a^2,1,sum)) 
        rownames(a) <- colnames(object@data)
        a <- round(a, digits)
        colnames(a) <- paste('a_', object@factorNames, sep='')            
        allPars <- list()
        if(allpars){
            if(length(object@pars[[1]]@SEpar) > 0){
                for(i in 1:(J+1)){
                    allPars[[i]] <- round(matrix(c(object@pars[[i]]@par, object@pars[[i]]@SEpar), 
                                             2, byrow = TRUE), digits)
                    rownames(allPars[[i]]) <- c('pars', 'SE')
                    colnames(allPars[[i]]) <- names(object@pars[[i]]@parnum)
                }
            } else {
                for(i in 1:(J+1))
                    allPars[[i]] <- round(object@pars[[i]]@par, digits)
            }                  
            names(allPars) <- c(rownames(a), 'GroupPars')                
        }        
        if(allpars) return(allPars)
            cat('\nItem parameters:\n')
        print(a)
        gpars <- ExtractGroupPars(object@pars[[J+1]])
        cat('\nGroup parameters:\n')
        cat('\nMeans:\n')            
        gmeans <- gpars$gmeans
        gcov <- gpars$gcov
        fnames <- object@factorNames
        fnames <- fnames[!grepl(pattern='\\(', fnames)]
        names(gmeans) <- colnames(gcov) <- rownames(gcov) <- fnames
        print(round(gmeans, digits))
        cat('\nCovariance:\n')
        print(round(gcov, digits))
        invisible(list(a,gpars))               	
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
    { 
        K <- object@K        
        data <- object@data    
        N <- nrow(data)	
        J <- ncol(data)
        nfact <- ncol(object@F) - length(attr(object@pars, 'prodlist'))
        if(object@pars[[1]]@bfactor) nfact <- 2
        itemloc <- object@itemloc
        res <- matrix(0,J,J)
        diag(res) <- NA
        colnames(res) <- rownames(res) <- colnames(data)
        theta <- seq(-4,4, length.out = round(20/nfact))
        Theta <- thetaComb(theta,nfact)    	
        ThetaShort <- Theta
        if(length(object@prodlist) > 0) Theta <- prodterms(Theta, object@prodlist)
        gpars <- ExtractGroupPars(object@pars[[length(object@pars)]])
        if(object@pars[[1]]@bfactor) prior <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact))
            else prior <- mvtnorm::dmvnorm(ThetaShort,gpars$gmeans,gpars$gcov)        
        prior <- prior/sum(prior)       	               
        if(restype == 'LD'){	
            for(i in 1:J){								
                for(j in 1:J){			
                    if(i < j){
                        P1 <- ProbTrace(x=object@pars[[i]], Theta=Theta)
                        P2 <- ProbTrace(x=object@pars[[j]], Theta=Theta)                        
                        tab <- table(data[,i],data[,j])		
                        Etab <- matrix(0,K[i],K[j])
                        for(k in 1:K[i])
                            for(m in 1:K[j])						
                                Etab[k,m] <- N * sum(P1[,k] * P2[,m] * prior)	
                        s <- gamma.cor(tab) - gamma.cor(Etab)
                        if(s == 0) s <- 1				
                        res[j,i] <- sum(((tab - Etab)^2)/Etab) /
                            ((K[i] - 1) * (K[j] - 1)) * sign(s)
                        res[i,j] <- sqrt( abs(res[j,i]) / (N - min(c(K[i],K[j]) - 1)))	
                    }
                }
            }	
            cat("LD matrix:\n\n")	
            res <- round(res,digits)
            return(res)
        } 
        if(restype == 'exp'){	
            r <- object@tabdata[ ,ncol(object@tabdata)]
            res <- round((r - object@Pl * nrow(object@data)) / 
                sqrt(object@Pl * nrow(object@data)),digits)
            expected <- round(N * object@Pl/sum(object@Pl),digits)  
            tabdata <- object@tabdata
            ISNA <- is.na(rowSums(tabdata))
            expected[ISNA] <- res[ISNA] <- NA
            tabdata <- data.frame(tabdata,expected,res)
            colnames(tabdata) <- c(colnames(object@tabdata),"exp","res")	
            if(!is.null(printvalue)){
                if(!is.numeric(printvalue)) stop('printvalue is not a number.')
                tabdata <- tabdata[abs(tabdata[ ,ncol(tabdata)]) > printvalue, ]
            }			
            return(tabdata)				
        }					
    }       
)

setMethod(
    f = "anova",
    signature = signature(object = 'ConfirmatoryClass'),
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

setMethod(
    f = "fitted",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, digits = 3, ...){  
        tabdata <- object@tabdata
        N <- nrow(object@data)
        expected <- round(N * object@Pl/sum(object@Pl),digits)
        print(cbind(tabdata,expected))
        invisible(tabdata)
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'ConfirmatoryClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45, 
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10), ...)
    {           
        class(x) <- 'ExploratoryClass'
        if(length(attr(x@pars, 'prodlist')) > 0 ) stop('No plots for models with polynomial and 
                                                       product terms')
        plot(x, type=type, npts=npts, theta_angle=theta_angle, rot=rot, ...)
        
    }		
)

