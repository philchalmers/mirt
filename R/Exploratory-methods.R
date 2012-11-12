#Methods 
setMethod(
    f = "print",
    signature = signature(x = 'ExploratoryClass'),
    definition = function(x){  
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
        cat("Full-information factor analysis with ", ncol(x@F), " factor",
            if(ncol(x@F)>1) "s", "\n", sep="")
        EMquad <- ''
        if(x@method == 'EM') EMquad <- c(' with ', x@quadpts, ' quadrature')
        if(x@converge == 1)	
            cat("Converged in ", x@iter, " iterations", EMquad, ". \n", sep = "")
        else 	
            cat("Estimation stopped after ", x@iter, " iterations", EMquad, ". \n", sep="")
        cat("Log-likelihood =", x@logLik, "\n")
        cat("AIC =", x@AIC, "\n")		
        cat("BIC =", x@BIC, "\n")
        if(!is.nan(x@p))            		    
            cat("G^2 = ", round(x@G2,2), ", df = ", x@df, ", p = ", round(x@p,4),
                "\nTLI = ", round(x@TLI,3), ", RMSEA = ", round(x@RMSEA,3), "\n", sep="")
        else             
            cat("G^2 = ", NA, ", df = ", 
                x@df, ", p = ", NA, ", RMSEA = ", NA, "\n", sep="" )		
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'ExploratoryClass'),
    definition = function(object){  
        print(object)
    }
)

setMethod(
    f = "summary",
    signature = 'ExploratoryClass',
    definition = function(object, rotate = '', Target = NULL, suppress = 0, digits = 3, 
                          verbose = TRUE, ...){        
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
            if(verbose){
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
            }			
            if(verbose){
                cat("\nRotation: ", rotate, "\n")                            
                cat("\nRotated factor loadings: \n\n")
                print(loads,digits)
                cat("\nRotated SS loadings: ",round(SS,digits), "\n")		
                cat("\nFactor correlations: \n\n")
                print(Phi)
            }
            if(any(h2 > 1)) 
                warning("Solution has heywood cases. Interpret with caution.") 
            invisible(list(rotF=rotF$loadings,h2=h2,fcor=Phi))  
        }  
    }
)

setMethod(
    f = "coef",
    signature = 'ExploratoryClass',
    definition = function(object, rotate = '', Target = NULL, digits = 3, ...){         
        K <- object@K
        J <- length(K)
        nfact <- ncol(object@F)
        a <- matrix(0, J, nfact)
        for(i in 1:J)
            a[i, ] <- ExtractLambdas(object@pars[[i]])        
        if (ncol(a) > 1){ 
            rotname <- ifelse(rotate == '', object@rotate, rotate)
            so <- summary(object, rotate=rotate, Target=Target, verbose=FALSE, digits=digits, ...)             
            a <- rotateLambdas(so) * 1.702/object@pars[[1]]@D
            for(i in 1:J)
                object@pars[[i]]@par[1:nfact] <- a[i, ]            
            object@pars[[J + 1]]@par[-c(1:nfact)] <- so$fcor[lower.tri(so$fcor, TRUE)]            
        }   
        allPars <- list()        
        if(length(object@pars[[1]]@SEpar) > 0){
            for(i in 1:(J+1)){
                allPars[[i]] <- round(matrix(c(object@pars[[i]]@par, object@pars[[i]]@SEpar), 
                                             2, byrow = TRUE), digits)
                rownames(allPars[[i]]) <- c('pars', 'SE')
                colnames(allPars[[i]]) <- names(object@pars[[i]]@parnum)
            } 
        } else {
            for(i in 1:(J+1)){
                allPars[[i]] <- round(object@pars[[i]]@par, digits)
                names(allPars[[i]]) <- names(object@pars[[i]]@parnum)
            }
        }                  
        names(allPars) <- c(colnames(object@data), 'GroupPars')
        return(allPars) 
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'ExploratoryClass'),
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

setMethod(
    f = "residuals",
    signature = signature(object = 'ExploratoryClass'),
    definition = function(object, restype = 'LD', digits = 3, df.p = FALSE, printvalue = NULL, 
                          verbose = TRUE, ...)
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
        Theta <- object@Theta
        prior <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact))
        prior <- prior/sum(prior) 
        df <- (object@K - 1) %o% (object@K - 1)
        diag(df) <- NA
        colnames(df) <- rownames(df) <- colnames(res)
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
                        df[i,j] <- pchisq(abs(res[j,i]), df=df[j,i], lower.tail=FALSE)
                    }
                }
            }	
            if(df.p){
                cat("Degrees of freedom (lower triangle) and p-values:\n\n")
                print(round(df, digits))
                cat("\n")
            }
            if(verbose) cat("LD matrix (lower triangle) and standardized values:\n\n")    
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
    f = "plot",
    signature = signature(x = 'ExploratoryClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45, 
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10), ...)
    {           
        if (!type %in% c('info','infocontour', 'SE')) stop(type, " is not a valid plot type.")
        if (any(theta_angle > 90 | theta_angle < 0)) 
            stop('Improper angle specifed. Must be between 0 and 90.')
        if(length(theta_angle) > 1) type = 'infoangle'
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])       
        nfact <- x@pars[[1]]@nfact
        if(nfact > 2) stop("Can't plot high dimensional solutions.")
        if(nfact == 1) theta_angle <- 0        
        J <- length(x@pars) - 1        
        theta <- seq(-4,4,length.out=npts)             
        ThetaFull <- Theta <- thetaComb(theta, nfact)        
        prodlist <- attr(x@pars, 'prodlist')
        if(length(prodlist) > 0)        
            ThetaFull <- prodterms(Theta,prodlist)        
        info <- 0            
        for(l in 1:length(theta_angle)){
            ta <- theta_angle[l]
            if(nfact == 2) ta <- c(theta_angle[l], 90 - theta_angle[l])
            for(i in 1:J)
                info <- info + iteminfo(x=x@pars[[i]], Theta=ThetaFull, degrees=ta)            
        }
        plt <- data.frame(cbind(info,Theta))         
        if(nfact == 2){						
            colnames(plt) <- c("info", "Theta1", "Theta2")			
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'infocontour')												
                return(contourplot(info ~ Theta1 * Theta2, data = plt, 
                                   main = paste("Test Information Contour"), xlab = expression(theta[1]), 
                                   ylab = expression(theta[2])))
            if(type == 'info')
                return(wireframe(info ~ Theta1 + Theta2, data = plt, main = "Test Information", 
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE))
            if(type == 'infoangle')
                symbols(plt[,2], plt[,3], circles = sqrt(plt[,1]/pi), inches = .35, fg='white', bg='blue', 
                        xlab = expression(theta[1]), ylab = expression(theta[2]), 
                        main = 'Information across different angles')
            if(type == 'SE')                
                return(wireframe(SE ~ Theta1 + Theta2, data = plt, main = "Test Standard Errors", 
                                 zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE))            
        } else {
            colnames(plt) <- c("info", "Theta")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'info')
                return(xyplot(info~Theta, plt, type='l',main = 'Test Information', 
                              xlab = expression(theta), ylab=expression(I(theta))))				
            if(type == 'infocontour') 
                cat('No \'contour\' plots for 1-dimensional models\n')
            if(type == 'SE')                
                xyplot(SE~Theta, plt, type='l',main = 'Test Standard Errors', 
                       xlab = expression(theta), ylab=expression(SE(theta)))
        }		
    }		
)	

setMethod(
    f = "fitted",
    signature = signature(object = 'ExploratoryClass'),
    definition = function(object, digits = 3, ...){  
        Exp <- round(nrow(object@data) * object@Pl,digits)  
        tabdata <- object@tabdata
        Exp[is.na(rowSums(tabdata))] <- NA				
        tabdata <- cbind(tabdata,Exp)			
        return(tabdata)        
    }
)
