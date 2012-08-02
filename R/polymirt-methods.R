# Methods
setMethod(
    f = "print",
    signature = signature(x = 'polymirtClass'),
    definition = function(x)
    {
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
        cat("Full-information factor analysis with ", ncol(x@F), " factor",
            if(ncol(x@F)>1) "s", "\n", sep="")
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
    signature = signature(object = 'polymirtClass'),
    definition = function(object) {
        print(object)			
    } 
)

setMethod(
    f = "summary",
    signature = 'polymirtClass',
    definition = function(object, rotate = '', Target = NULL, suppress = 0, digits = 3, 
                          print = TRUE, ...)
    {
        nfact <- ncol(object@F)
        itemnames <- colnames(object@data)
        if (rotate == 'none' || nfact == 1) {
            F <- object@F
            F[abs(F) < suppress] <- NA
            h2 <- as.matrix(object@h2)    	
            SS <- apply(F^2,2,sum)
            colnames(h2) <- "h2"	
            names(SS) <- colnames(F) 
            loads <- round(cbind(F,h2),digits)
            rownames(loads) <- itemnames
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
            rownames(loads) <- itemnames			
            Phi <- diag(nfact)
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
                cat("\nSS loadings: ",round(SS,digits), "\n")		
                if(!rotF$orthogonal) 
                    cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
            }
            if(any(h2 > 1)) 
                warning("Solution has heywood cases. Interpret with caution.") 
            invisible(list(rotF=rotF$loadings,h2=h2,fcor=Phi))  
        }  
    }
)

setMethod(
    f = "coef",
    signature = 'polymirtClass',
    definition = function(object, rotate = '', Target = NULL, SE = TRUE, digits = 3, ...)
    {          
        K <- object@K
        a <- object@pars$lambdas		
        d <- matrix(NA, nrow(a), max(K-1))
        zetas <- object@pars$zetas
        for(i in 1:length(K)){
            d[i, 1:(K[i] - 1)] <- zetas[[i]]
        }
        A <- sqrt(apply(a^2,1,sum))
        B <- -d/A  
        if (ncol(a) > 1){  
            rotname <- ifelse(rotate == '', object@rotate, rotate)
            so <- summary(object, rotate = rotate, Target = Target, print = FALSE, ...)             
            a <- rotateLambdas(so)
            parameters <- cbind(a,d,object@guess,A,B)    
            colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""),paste("d_",1:max(K-1),sep=""),"guess", 
                                      "mvdisc",paste("mvint_",1:max(K-1),sep=""))	  
            rownames(parameters) <- colnames(object@data)
            cat("\nParameters with", rotname, "rotation, multivariate discrimination and intercept: \n\n")
            print(round(parameters, digits))  	
        } else {
            parameters <- cbind(a,d,object@guess, object@upper)
            colnames(parameters) <- c(paste("a_",1:ncol(a),sep=""), 
                                      paste("d_",1:max(K-1),sep=""),"guess","upper")   
            rownames(parameters) <- colnames(object@data)
            cat("\nParameter slopes and intercepts: \n\n")	
            print(round(parameters, digits))	  
        }
        ret <- list(parameters)
        if(length(object@SEpars) > 1){
            if(SE){
                cat("\nStd. Errors: \n\n")	
                SEs <- cbind(object@SEpars,object@SEg,object@SEup)
                colnames(SEs) <- c(paste("a_",1:ncol(a),sep=""),paste("d_",1:max(K-1),sep=""),
                                   "SEguess", "SEupper") 	
                rownames(SEs) <- rownames(parameters)
                print(SEs, digits)
                ret <- list(parameters,SEs)
            }
        }
        invisible(ret)
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'polymirtClass', y = "missing"),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45, 
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10))
    {  		
        if (!type %in% c('info','infocontour')) stop(type, " is not a valid plot type.")
        if (any(theta_angle > 90 | theta_angle < 0)) 
            stop('Improper angle specifed. Must be between 0 and 90.')
        if(length(theta_angle) > 1) type = 'infoangle'        
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
        K <- x@K		
        nfact <- ncol(x@Theta)
        if(nfact > 2) stop("Can't plot high dimensional solutions.")
        a <- x@pars$lambdas
        d <- x@pars$zetas	
        guess <- x@guess
        guess[is.na(guess)] <- 0
        upper <- x@upper
        upper[is.na(upper)] <- 1
        if(nfact == 2){
            theta_angle2 <- c(90 - theta_angle)
            angles <- rbind(theta_angle, theta_angle2)
            cosalpha <- cos(d2r(angles))
            A <- list()
            if(length(theta_angle) == 1)
                A[[1]] <- as.matrix(sqrt(rowSums((a * matrix(cosalpha[ ,1], nrow(a), 2, TRUE))^2)))
            else                 
                for(i in 1:ncol(cosalpha))
                    A[[i]] <- as.matrix(sqrt(rowSums((a * matrix(cosalpha[ ,i], nrow(a), 2, TRUE))^2)))                                
        }   
        theta <- if(length(theta_angle) == 1) seq(-4,4,length.out=npts) 
        else seq(-4,4,length.out=9)
        Theta <- thetaComb(theta, nfact)        
        info <- test_info(a=a, d=d, Theta=Theta, Alist=A, guess=guess, upper=upper, K=K)        
        plt <- data.frame(cbind(info,Theta))
        if(nfact == 2){						
            colnames(plt) <- c("info", "Theta1", "Theta2")			
            if(type == 'infocontour')												
                return(contourplot(info ~ Theta1 * Theta2, data = plt, 
                                   main = paste("Test Information Contour"), xlab = expression(theta[1]), 
                                   ylab = expression(theta[2])))
            if(type == 'info')
                return(wireframe(info ~ Theta1 + Theta2, data = plt, main = "Test Information", 
                                 zlab = expression(I(theta)), xlab = expression(theta[1]), ylab = expression(theta[2]), 
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE))
            if(type == 'infoangle')
                symbols(plt[,2], plt[,3], circles = sqrt(plt[,1]/pi), inches = .35, fg='white', bg='blue', 
                        xlab = expression(theta[1]), ylab = expression(theta[2]), 
                        main = 'Information across different angles')            
        } else {
            colnames(plt) <- c("info", "Theta")
            if(type == 'info')
                return(xyplot(info~Theta, plt, type='l',main = 'Test Information', xlab = expression(theta), 
                              ylab = expression(I(theta))))
            if(type == 'infocontour') 
                cat('No \'contour\' plots for 1-dimensional models\n')
        }		
    }	  
)	

setMethod(
    f = "residuals",
    signature = signature(object = 'polymirtClass'),
    definition = function(object, restype = 'LD', digits = 3, printvalue = NULL, ...)
    { 	
        K <- object@K		
        data <- object@data	
        N <- nrow(data)	
        J <- ncol(data)
        nfact <- ncol(object@F)
        lambdas <- object@pars$lambdas
        zetas <- object@pars$zetas
        guess <- object@guess		
        itemloc <- object@itemloc
        theta <- seq(-4,4, length.out = round(20/nfact))
        Theta <- thetaComb(theta,nfact)
        res <- matrix(0,J,J)
        diag(res) <- NA
        colnames(res) <- rownames(res) <- colnames(data)
        prior <- mvtnorm::dmvnorm(Theta,rep(0,nfact),diag(nfact))
        prior <- prior/sum(prior)	
        if(restype == 'LD'){	
            for(i in 1:J){								
                for(j in 1:J){			
                    if(i < j){
                        if(K[i] > 2) P1 <- P.poly(lambdas[i,],zetas[[i]],Theta,itemexp=TRUE)
                        else { 
                            P1 <- P.mirt(lambdas[i,],zetas[[i]], Theta, guess[i])
                            P1 <- cbind(1 - P1, P1)
                        }	
                        if(K[j] > 2) P2 <- P.poly(lambdas[j,],zetas[[j]],Theta,itemexp=TRUE)
                        else {
                            P2 <- P.mirt(lambdas[j,],zetas[[j]], Theta, guess[j])	
                            P2 <- cbind(1 - P2, P2)
                        }
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
            r <- object@tabdata[ ,ncol(object@tabdata)-1]
            rexp <- object@tabdata[ ,ncol(object@tabdata)]
            res <- round((r - rexp) / sqrt(rexp),digits)
            tabdata <- object@tabdata
            freq <- tabdata[ ,ncol(tabdata)]			
            tabdata[tabdata[ ,1:ncol(object@data)] == 99] <- NA
            tabdata[ ,ncol(tabdata)] <- freq
            tabdata <- cbind(tabdata,res)
            colnames(tabdata) <- c(colnames(data),'freq', 'exp', 'std_res')	
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
    signature = signature(object = 'polymirtClass'),
    definition = function(object, object2, ...)
    {
        dots <- list(...)				
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
            " (SE = ", se,"), df = ", df, ", p = ", round(1 - pchisq(X2,abs(df)),4),
            "\n", sep="")
        cat("AIC difference = ", round(AICdiff,3)," (SE = ", se,")\n", sep='')  
        cat("BIC difference = ", round(BICdiff,3)," (SE = ", se,")\n", sep='') 
    }		
) 

setMethod(
    f = "fitted",
    signature = signature(object = 'polymirtClass'),
    definition = function(object, digits = 3, ...){  		  
        tabdata <- object@tabdata		
        colnames(tabdata) <- c(colnames(object@data),"freq","exp")
        r <- round(tabdata[,ncol(tabdata)], digits)	
        print(cbind(tabdata[,-ncol(tabdata)],r))
        invisible(tabdata)
    }
)