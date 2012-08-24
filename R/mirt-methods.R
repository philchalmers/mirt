#Methods 
setMethod(
    f = "print",
    signature = signature(x = 'mirtClass'),
    definition = function(x){  
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
        cat("Full-information factor analysis with ", ncol(x@F), " factor",
            if(ncol(x@F)>1) "s", "\n", sep="")
        if(x@converge == 1)	
            cat("Converged in ", x@EMiter, " iterations using ", x@quadpts,
                " quadrature points.\n", sep="")
        else 	
            cat("Estimation stopped after ", x@EMiter, " iterations using ", 
                x@quadpts, " quadrature points.\n", sep="")
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
    signature = signature(object = 'mirtClass'),
    definition = function(object){  
        print(object)
    }
)

setMethod(
    f = "summary",
    signature = 'mirtClass',
    definition = function(object, rotate = '', Target = NULL, suppress = 0, digits = 3, 
                          print = TRUE, ...){        
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
    }
)

setMethod(
    f = "coef",
    signature = 'mirtClass',
    definition = function(object, rotate = '', Target = NULL, allpars = FALSE, digits = 3, ...){  
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
                for(i in 1:J){
                    allPars[[i]] <- round(matrix(c(object@pars[[i]]@par, object@pars[[i]]@SEpar), 
                                                 2, byrow = TRUE), digits)
                    rownames(allPars[[i]]) <- c('pars', 'SE')
                    colnames(allPars[[i]]) <- names(object@pars[[i]]@parnum)
                } 
            } else {
                for(i in 1:J)
                    allPars[[i]] <- round(object@pars[[i]]@par, digits)                
            }       
            names(allPars) <- rownames(a)
        }        
        ret <- if(allpars) allPars else a
        if(nfact > 1) cat('\nRotation:', rotname, '\n\n')
        ret
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'mirtClass'),
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
    signature = signature(object = 'mirtClass'),
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
        Theta <- object@Theta
        prior <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact))
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
    f = "plot",
    signature = signature(x = 'mirtClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45, 
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10), ...)
    {           
        if (!type %in% c('info','infocontour')) stop(type, " is not a valid plot type.")
        if (any(theta_angle > 90 | theta_angle < 0)) 
            stop('Improper angle specifed. Must be between 0 and 90.')
        if(length(theta_angle) > 1) type = 'infoangle'
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
        K <- x@K		
        nfact <- ncol(x@Theta)
        if(nfact > 2) stop("Can't plot high dimensional solutions.")
        J <- length(x@pars)
        a <- coef(x, ...)  
        A <- list()
        if(nfact == 2){
            a <- a[,1:2]
            theta_angle2 <- c(90 - theta_angle)
            angles <- rbind(theta_angle, theta_angle2)
            cosalpha <- cos(d2r(angles))            
            if(length(theta_angle) == 1)
                A[[1]] <- as.matrix(sqrt(rowSums((a * matrix(cosalpha[ ,1], nrow(a), 2, TRUE))^2)))
            else                 
                for(i in 1:ncol(cosalpha))
                    A[[i]] <- as.matrix(sqrt(rowSums((a * matrix(cosalpha[ ,i], nrow(a), 2, TRUE))^2)))                                
        } else A <- list(a)  
        theta <- if(length(theta_angle) == 1) seq(-4,4,length.out=npts) 
            else seq(-4,4,length.out=9)
        Theta <- thetaComb(theta, nfact)        
        info <- test_info(pars=x@pars, Theta=Theta, Alist=A, K=K)         
        plt <- data.frame(cbind(info,Theta))
        if(nfact == 2){						
            colnames(plt) <- c("info", "Theta1", "Theta2")			
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
        } else {
            colnames(plt) <- c("info", "Theta")
            if(type == 'info')
                return(xyplot(info~Theta, plt, type='l',main = 'Test Information', 
                              xlab = expression(theta), ylab=expression(I(theta))))				
            if(type == 'infocontour') 
                cat('No \'contour\' plots for 1-dimensional models\n')
        }		
    }		
)	

setMethod(
    f = "fitted",
    signature = signature(object = 'mirtClass'),
    definition = function(object, digits = 3, ...){  
        Exp <- round(nrow(object@data) * object@Pl,digits)  
        tabdata <- object@tabdata
        Exp[is.na(rowSums(tabdata))] <- NA				
        tabdata <- cbind(tabdata,Exp)			
        print(tabdata)
        invisible(tabdata)
    }
)
