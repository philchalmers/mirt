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
    definition = function(object, digits = 3, ...)
    {
        if(any(object@estComp)) stop('No factor metric for noncompensatory models')
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
        Phi <- cov2cor(object@gpars$sig)	  
        Phi <- round(Phi, digits)
        colnames(Phi) <- rownames(Phi) <- colnames(F)
        print(Phi)				
        invisible(F)  	  
    }
)

setMethod(
    f = "coef",
    signature = 'confmirtClass',
    definition = function(object, SE = TRUE, print.gmeans = FALSE, digits = 3, ...)
    {  
        nfact <- ncol(object@nfact)
        nfactNames <- ifelse(length(object@prodlist) > 0, 
                             length(object@prodlist) + nfact, nfact)
        factorNames <- colnames(object@F)
        itemnames <- names(object@h2)
        a <- matrix(object@parsprint[ ,1:nfactNames], ncol=nfactNames)
        d <- matrix(object@parsprint[ ,(nfactNames+1):ncol(object@parsprint)],
                    ncol = ncol(object@parsprint)-nfactNames)    	
        parameters <- cbind(object@parsprint,object@guess,object@upper)
        SEs <- cbind(object@SEpars,object@SEg,object@SEup)
        rownames(parameters) <- itemnames
        rownames(SEs) <- itemnames
        colnames(SEs) <- colnames(parameters) <- c(paste("a_",factorNames[1:nfactNames],sep=""),
                                                   paste("d_",1:(ncol(object@parsprint)-nfactNames),sep=""),"guess",'upper')
        factorNames2 <- factorNames	
        if(nfact < nfactNames)
            factorNames2 <- factorNames[!grepl("\\(",factorNames)]			
        cat("\nITEM PARAMETERS: \n")
        print(parameters, digits)
        if(SE){
            cat("\nStd. Errors: \n")	
            print(SEs, digits)
        }	
        u <- object@gpars$u	
        SEu <- object@SEgpars$SEu
        sig <- object@gpars$sig
        SEsig <- as.matrix(object@SEgpars$SEsig)
        names(u) <- colnames(sig) <- rownames(sig) <- factorNames2	
        cat("\nGROUP PARAMETERS: \n")
        if(print.gmeans){
            cat("Means: \n")
            print(u,digits)
            cat("\nStd. Errors: \n")			
            names(SEu) <- names(u) 	
            print(SEu, digits)	
        }
        cat("Covariance: \n")
        print(sig,digits)
        if(SE){
            cat("\nStd. Errors: \n")			
            colnames(SEsig) <- rownames(SEsig) <- factorNames2	
            print(SEsig, digits)	
        }
        invisible(list(parsprint = parameters,mu = u,sigma = sig, sigmaSE = SEsig,
                       muSE = SEu))	
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'confmirtClass'),
    definition = function(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
    { 
        fulldata <- object@fulldata	
        data <- object@data
        data[data==99] <- NA		
        N <- nrow(data)
        K <- object@K
        J <- length(K)
        sig <- object@gpars$sig	
        nfact <- ncol(sig)
        nfactNames <- ncol(object@F)
        theta <- seq(-4,4, length.out = round(20/nfact))
        Theta <- thetaComb(theta,nfact)		
        if(length(object@prodlist) > 0) Theta <- prodterms(Theta, object@prodlist)
        lambdas <- object@pars$lambdas
        zetas <- object@pars$zetas
        guess <- object@guess
        guess[is.na(guess)] <- 0	
        upper <- object@upper
        upper[is.na(upper)] <- 1
        Ksums <- cumsum(K) - 1	
        itemloc <- object@itemloc
        res <- matrix(0,J,J)
        diag(res) <- NA
        colnames(res) <- rownames(res) <- colnames(data)
        prior <- mvtnorm::dmvnorm(Theta[,1:nfact,drop=FALSE],rep(0,nfact),sig)
        prior <- prior/sum(prior)		
        if(restype == 'LD'){	
            for(i in 1:J){				
                for(j in 1:J){			
                    if(i < j){
                        if(K[i] > 2) P1 <- P.poly(lambdas[i,],zetas[[i]],Theta,itemexp=TRUE)
                        else { 
                            P1 <- P.mirt(lambdas[i,],zetas[[i]], Theta, guess[i], upper[i])
                            P1 <- cbind(1 - P1, P1)
                        }	
                        if(K[j] > 2) P2 <- P.poly(lambdas[j,],zetas[[j]],Theta,itemexp=TRUE)
                        else {
                            P2 <- P.mirt(lambdas[j,],zetas[[j]], Theta, guess[j], upper[j])	
                            P2 <- cbind(1 - P2, P2)
                        }
                        tab <- table(data[,i],data[,j])		
                        Etab <- matrix(0,K[i],K[j])
                        for(k in 1:K[i])
                            for(m in 1:K[j])						
                                Etab[k,m] <- N * sum(P1[,k] * P2[,m] * prior)	
                        s <- gamma.cor(tab) - gamma.cor(Etab)
                        if(s == 0) s <- 1				
                        res[j,i] <- sum(((tab - Etab)^2)/Etab) * sign(s)
                        res[i,j] <- sqrt(abs(res[j,i]) / (N * min(c(K[i],K[j]) - 1))) 					
                    }					
                }
            }		
            if(is.null(printvalue)) cat("LD matrix:\n\n")	
            res <- round(res,digits)    	
            return(res)
        }
        if(restype == 'exp'){
            if(length(object@tabdata) == 0) stop('Expected response vectors cannot be computed 
                because logLik() has not been run or the data contains missing responses.')
			tabdata <- object@tabdata
            res <- (tabdata[,J+1] - tabdata[,J+2]) / sqrt(tabdata[,J+2])
            tabdata <- round(cbind(tabdata,res),digits)
            colnames(tabdata) <- c(colnames(object@data), 'freq', 'exp', 'std_res')
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
