# Methods
setMethod(
    f = "print",
    signature = signature(x = 'confmirtClass'),
    definition = function(x)
    {
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
        cat("Full-information item factor analysis with ", x@nfact, " factors \n", sep="")
        if(x@converge == 1)    
            cat("Converged in ", x@cycles, " iterations.\n", sep="")
        else 	
            cat("Estimation stopped after ", x@cycles, " iterations.\n", sep="")
        
        if(length(x@logLik) > 0){
            cat("Log-likelihood = ", x@logLik,", SE = ",round(x@SElogLik,3), "\n",sep='')			
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
    signature = signature(object = 'confmirtClass'),
    definition = function(object) {
        print(object)
    } 
)

setMethod(
    f = "summary",
    signature = 'confmirtClass',
    definition = function(object, rotate = '', Target = NULL, suppress = 0, digits = 3, 
                          print = TRUE, ...)
    {        
        if(object@exploratory){
            nfact <- ncol(object@F)
            if (rotate == 'none' || nfact == 1) {
                F <- object@F
                F[abs(F) < suppress] <- NA
                h2 <- as.matrix(object@h2)    			
                SS <- apply(F^2,2,sum)
                colnames(h2) <- "h2"			
                names(SS) <- colnames(F)
                loads <- round(cbind(F,h2),digits)
                rownames(loads) <- colnames(object@data)
                if(print){
                    cat("\nUnrotated factor loadings: \n\n")
                    print(loads)	    	 
                    cat("\nSS loadings: ",round(SS,digits), "\n")
                    cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
                }
                invisible(list(rotF=F,h2=h2,fcor=matrix(1)))
            } else {	
                F <- object@F
                h2 <- as.matrix(object@h2)
                colnames(h2) <- "h2"
                if(rotate == ''){
                    rotate <- object@rotate
                    Target <- object@Target
                }
                rotF <- Rotate(F, rotate, Target = Target, ...)            
                SS <- apply(rotF$loadings^2,2,sum)
                L <- rotF$loadings
                L[abs(L) < suppress] <- NA	
                loads <- round(cbind(L,h2),digits)
                rownames(loads) <- colnames(object@data)			
                Phi <- diag(ncol(F))			
                if(!rotF$orthogonal){
                    Phi <- rotF$Phi	  
                    Phi <- round(Phi, digits)
                    colnames(Phi) <- rownames(Phi) <- colnames(F)
                    if(print){
                        cat("\nFactor correlations: \n\n")
                        print(Phi)            
                    }
                }			
                if(print){
                    cat("\nRotation: ", rotate, "\n")
                    cat("\nRotated factor loadings: \n\n")
                    print(loads,digits)
                    cat("\nRotated SS loadings: ",round(SS,digits), "\n")		
                }
                if(any(h2 > 1)) 
                    warning("Solution has heywood cases. Interpret with caution.") 
                invisible(list(rotF=rotF$loadings,h2=h2,fcor=Phi))  
            }
        } else {
            if(length(object@prodlist) > 0) stop('No factor metric for models with product terms')
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
            colnames(Phi) <- rownames(Phi) <- colnames(F)
            print(Phi)				
            invisible(F)  	  
        }
    }
)

setMethod(
    f = "coef",
    signature = 'confmirtClass',
    definition = function(object, rotate = '', Target = NULL, allpars = FALSE, digits = 3, ...)
    {          
        if(object@exploratory){
            K <- object@K
            J <- length(K)
            nfact <- ncol(object@F)
            a <- matrix(0, J, nfact)
            for(i in 1:J)
                a[i, ] <- ExtractLambdas(object@pars[[i]])        
            A <- sqrt(apply(a^2,1,sum))                        
            if (ncol(a) > 1){ 
                rotname <- ifelse(rotate == '', object@rotate, rotate)
                so <- summary(object, rotate = rotate, Target = Target, print = FALSE, ...)             
                a <- rotateLambdas(so)
            }   
            rownames(a) <- colnames(object@data)
            if(nfact > 1){
                a <- round(cbind(a, A), digits)
                colnames(a) <- c(paste('a', 1:nfact, sep=''), 'MV_disc')
            } else {
                a <- round(a, digits)
                colnames(a) <- paste('a', 1:nfact, sep='')
            }
            allPars <- list()
            if(allpars){
                if(length(object@pars[[1]]@SEpar) > 0){
                    for(i in 1:(J+1))
                        allPars[[i]] <- round(matrix(c(object@pars[[i]]@par, object@pars[[i]]@SEpar), 
                                                     2, byrow = TRUE), digits)
                } else {
                    for(i in 1:(J+1))
                        allPars[[i]] <- round(object@pars[[i]]@par, digits)
                }       
                names(allPars) <- rownames(a)
            }        
            ret <- if(allpars) allPars else a
            if(nfact > 1) cat('\nRotation:', rotname, '\n\n')
            print(ret)
            return(invisible(ret))
        } else {               
            K <- object@K
            J <- length(K)
            nfact <- ncol(object@F)
            a <- matrix(0, J, nfact)
            for(i in 1:J)
                a[i, ] <- ExtractLambdas(object@pars[[i]])        
            A <- sqrt(apply(a^2,1,sum))                                    
            rownames(a) <- colnames(object@data)
            a <- round(a, digits)
            colnames(a) <- paste('a', 1:nfact, sep='')            
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
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'confmirtClass'),
    definition = function(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
    { 
        K <- object@K        
        data <- object@data    
        N <- nrow(data)	
        J <- ncol(data)
        nfact <- ncol(object@F)        
        itemloc <- object@itemloc
        res <- matrix(0,J,J)
        diag(res) <- NA
        colnames(res) <- rownames(res) <- colnames(data)
        theta <- seq(-4,4, length.out = round(20/nfact))
        Theta <- thetaComb(theta,nfact)    	
        if(length(object@prodlist) > 0) Theta <- prodterms(Theta, object@prodlist)
        gpars <- ExtractGroupPars(object@pars[[length(object@pars)]])
        prior <- mvtnorm::dmvnorm(Theta,gpars$gmeans,gpars$gcov)
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
    signature = signature(object = 'confmirtClass'),
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
            " (SE = ",se,"), df = ", df, ", p = ", round(1 - pchisq(X2,abs(df)),4), 
            "\n", sep="")
        cat("AIC difference = ", round(AICdiff,3)," (SE = ", se,")\n", sep='')
        cat("BIC difference = ", round(BICdiff,3)," (SE = ", se,")\n", sep='')
    }		
)

setMethod(
    f = "fitted",
    signature = signature(object = 'confmirtClass'),
    definition = function(object, digits = 3, ...){  
        tabdata <- object@tabdata		
        colnames(tabdata) <- c(colnames(object@data),"freq","exp")
        r <- round(tabdata[,ncol(tabdata)], digits)	
        print(cbind(tabdata[,-ncol(tabdata)],r))
        invisible(tabdata)
    }
)
