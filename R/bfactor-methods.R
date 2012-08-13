# Methods
setMethod(
    f = "print",
    signature = signature(x = 'bfactorClass'),
    definition = function(x){  
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")		
        cat("Full-information bifactor analysis with ", 
            length(unique(x@specific)), " specific factors \n", sep='')
        if(x@converge == 1)	
            cat("Converged in ", x@EMiter, " iterations using ",x@quadpts,
                " quadrature points. \n", sep="")
        else 	
            cat("Estimation stopped after ", x@EMiter, " iterations using ",x@quadpts,
                " quadrature points. \n", sep="")
        cat("Log-likelihood = ", x@logLik, "\n")
        cat("AIC = ", x@AIC, "\n")		
        cat("BIC = ", x@BIC, "\n")
        if(!is.nan(x@p))
            cat("G^2 = ", round(x@X2,2), ", df = ", x@df, ", p = ", round(x@p,4), 
                "\nTLI = ", round(x@TLI,3), ", RMSEA = ", round(x@RMSEA,3), "\n", sep="")
        else 
            cat("G^2 = ", NA, ", df = ", 
                x@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="")		
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'bfactorClass'),
    definition = function(object){  
        print(object)			
    }
)

setMethod(
    f = "summary",
    signature = signature(object = 'bfactorClass'),
    definition = function(object, digits = 3, ...){
        F <- round(object@F,digits)
        SS <- colSums(F^2)	
        F[!object@logicalfact] <- NA
        h2 <- round(object@h2,digits)					
        names(h2) <- "h2"		
        loads <- round(cbind(F,h2),digits)
        rownames(loads) <- colnames(object@data)	 
        cat("\nFactor loadings: \n\n")
        print(loads)
        cat("\nSS loadings: ",round(SS,digits), "\n")
        cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
        if(any(h2 > 1)) 
            warning("Solution has heywood cases. Interpret with caution.")
    }
)

setMethod(
    f = "coef",
    signature = signature(object = 'bfactorClass'),
    definition = function(object, digits = 3, ...){
        K <- object@K
        a <- object@pars$lambdas			
        d <- matrix(NA, nrow(a), max(K-1))
        zetas <- object@pars$zetas
        for(i in 1:length(K)){
            d[i, 1:(K[i] - 1)] <- zetas[[i]]
        }
        A <- sqrt(apply(a^2,1,sum))
        B <- -d/A 
        a[!attr(object@pars,'lamsel')] <- NA	
        parameters <- round(cbind(a,d,object@guess,object@upper,A,B),digits)
        colnames(parameters) <- c('a_G',paste("a_", 1:(ncol(object@F)-1),sep=""),
                                  paste("d_", 1:(max(K)-1),sep=""), "guess", "upper", "mvdisc", 
                                  paste("mvint_", 1:(max(K)-1),sep=""))  
        rownames(parameters) <- colnames(object@data)	
        cat("\nParameters with multivariate discrimination and intercept: \n\n")		
        print(parameters)
        ret <- list(parameters)
        if(length(object@parsSE) > 1){
            cat("\nStd. Errors: \n\n")	
            a <- object@parsSE$lambdas			
            d <- matrix(NA, nrow(a), max(K-1))
            zetas <- object@parsSE$zetas
            for(i in 1:length(K)){
                d[i, 1:(K[i] - 1)] <- zetas[[i]]
            }
            SEs <- cbind(a,d)
            colnames(SEs) <- c('a_G',paste("a_", 1:(ncol(object@F)-1),sep=""),
                               paste("d_", 1:(max(K)-1),sep=""))
            rownames(SEs) <- colnames(object@data)		
            print(SEs, digits)
            ret <- list(parameters, SEs)
        }	
        invisible(ret)
    }
)

setMethod( 
    f = "residuals",
    signature = signature(object = 'bfactorClass'),
    definition = function(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
    {       
        K <- object@K
        lf <- attr(object@pars, 'lamsel')
        Theta <- object@Theta
        data <- object@data	
        N <- nrow(data)	
        J <- ncol(data)		
        lambdas <- object@pars$lambdas
        zetas <- object@pars$zetas
        guess <- object@guess
        upper <- object@upper
        guess[is.na(guess)] <- 0
        upper[is.na(upper)] <- 1
        itemloc <- object@itemloc
        res <- matrix(0,J,J)
        diag(res) <- NA
        colnames(res) <- rownames(res) <- colnames(data)
        prior <- mvtnorm::dmvnorm(Theta,rep(0,2),diag(2))
        prior <- prior/sum(prior)	
        if(restype == 'LD'){	
            for(i in 1:J){								
                for(j in 1:J){			
                    if(i < j){
                        P1 <- P.bfactor(lambdas[i, ], zetas[[i]], Theta, guess[i], upper[i], lf[i, ])
                        P2 <- P.bfactor(lambdas[j, ], zetas[[j]], Theta, guess[j], upper[j], lf[j, ])
                        if(K[i] == 2) P1 <- cbind(1-P1, P1)
                        if(K[j] == 2) P2 <- cbind(1-P2, P2)						
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
    f = "fitted",
    signature = signature(object = 'bfactorClass'),
    definition = function(object, digits = 3, ...){  
        Exp <- round(object@N * object@Pl/sum(object@Pl),digits)  
        tabdata <- object@tabdata
        Exp[is.na(rowSums(tabdata))] <- NA				
        tabdata <- cbind(tabdata,Exp)		
        tabdata
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'bfactorClass'),
    definition = function(object, object2){           	
        df <- object@df - object2@df  
        if(df < 0){
            temp <- object
            object <- object2
            object2 <- temp
        }
        X2 <- 2*object2@logLik - 2*object@logLik 		
        AICdiff <- object@AIC - object2@AIC    
        BICdiff <- object@BIC - object2@BIC
        cat("\nChi-squared difference: \n\nX2 = ", round(X2,3), ", df = ",
            df, ", p = ", round(1 - pchisq(X2,abs(df)),4), "\n", sep="")
        cat("AIC difference = ", round(AICdiff,3), "\n")  
        cat("BIC difference = ", round(BICdiff,3), "\n")
    }
)
