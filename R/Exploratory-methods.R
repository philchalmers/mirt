#Methods
setMethod(
    f = "print",
    signature = signature(x = 'ExploratoryClass'),
    definition = function(x){
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
                cat("\nTLI (G2) = ", round(x@TLI,3), "; TLI (X2) = ", round(x@TLI.X2,3), '\n', sep='')
            }
        }
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
            Phi <- diag(ncol(F))
            colnames(h2) <- "h2"
            rownames(Phi) <- colnames(Phi) <- names(SS) <- colnames(F)
            loads <- round(cbind(F,h2),digits)            
            rownames(loads) <- colnames(object@data)            
            if(verbose){
                cat("\nUnrotated factor loadings: \n\n")
                print(loads)
                cat("\nSS loadings: ",round(SS,digits), "\n")
                cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
                cat("\nFactor correlations: \n\n")
                print(Phi)                
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
    definition = function(object, CI = .95, rotate = '', Target = NULL, digits = 3, 
                          verbose = TRUE, ...){
        if(CI >= 1 || CI <= 0)
            stop('CI must be between 0 and 1')
        z <- abs(qnorm((1 - CI)/2))
        SEnames <- paste0('CI_', c((1 - CI)/2*100, ((1 - CI)/2 + CI)*100))
        K <- object@K
        J <- length(K)
        nfact <- ncol(object@F)
        a <- matrix(0, J, nfact)
        for(i in 1:J)
            a[i, ] <- ExtractLambdas(object@pars[[i]])
        if (ncol(a) > 1){
            rotname <- ifelse(rotate == '', object@rotate, rotate)
            if(verbose) cat("\nRotation: ", rotname, "\n\n")
            so <- summary(object, rotate=rotate, Target=Target, verbose=FALSE, digits=digits, ...)
            a <- rotateLambdas(so) * 1.702
            for(i in 1:J)
                object@pars[[i]]@par[1:nfact] <- a[i, ]            
            if(rotname != 'none')
                object@pars[[J + 1]]@par[-c(1:nfact)] <- so$fcor[lower.tri(so$fcor, TRUE)]
        }
        allPars <- list()
        if(length(object@pars[[1]]@SEpar) > 0){
            for(i in 1:(J+1)){
                allPars[[i]] <- round(matrix(c(object@pars[[i]]@par, 
                                               object@pars[[i]]@par - z*object@pars[[i]]@SEpar,
                                               object@pars[[i]]@par + z*object@pars[[i]]@SEpar),
                                             3, byrow = TRUE), digits)
                rownames(allPars[[i]]) <- c('par', SEnames)
                colnames(allPars[[i]]) <- names(object@pars[[i]]@est)
            }
        } else {
            for(i in 1:(J+1)){
                allPars[[i]] <- matrix(round(object@pars[[i]]@par, digits), 1L)
                colnames(allPars[[i]]) <- names(object@pars[[i]]@est)
                rownames(allPars[[i]]) <- 'par'
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
        X2 <- round(2*object2@logLik - 2*object@logLik, 3)
        cat('\nModel 1: ')
        print(object@Call)
        cat('Model 2: ')
        print(object2@Call)
        cat('\n')
        ret <- data.frame(Df = c(object@df, object2@df),
                          AIC = c(object@AIC, object2@AIC),
                          AICc = c(object@AICc, object2@AICc),
                          BIC = c(object@BIC, object2@BIC),
                          SABIC = c(object@SABIC, object2@SABIC),
                          logLik = c(object@logLik, object2@logLik),
                          X2 = c('', X2),
                          df = c('', abs(df)),
                          p = c('', round(1 - pchisq(X2,abs(df)),3)))
        return(ret)
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'ExploratoryClass'),
    definition = function(object, restype = 'LD', digits = 3, df.p = FALSE, full.scores = FALSE, 
                          printvalue = NULL, verbose = TRUE, ...)
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
                        res[j,i] <- sum(((tab - Etab)^2)/Etab) * sign(s)
                        res[i,j] <- sign(res[j,i]) * sqrt( abs(res[j,i]) / (N*min(c(K[i],K[j]) - 1L)))
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
            expected <- round(N * object@Pl,digits)
            tabdata <- object@tabdata
            rownames(tabdata) <- NULL
            ISNA <- is.na(rowSums(tabdata))
            expected[ISNA] <- res[ISNA] <- NA
            tabdata <- data.frame(tabdata,expected,res)
            colnames(tabdata) <- c(colnames(object@tabdata),"exp","res")            
            if(full.scores){
                tabdata[, 'exp'] <- object@Pl / r * N
                tabdata2 <- object@tabdatalong
                tabdata2 <- tabdata2[,-ncol(tabdata2)]
                stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
                sfulldata <- apply(object@fulldata, 1, paste, sep='', collapse = '/')
                scoremat <- tabdata[match(sfulldata, stabdata2), 'exp', drop = FALSE]                
                res <- (1-scoremat) / sqrt(scoremat)
                colnames(res) <- 'res'
                ret <- cbind(object@data, scoremat, res)
                ret[is.na(rowSums(ret)), c('exp', 'res')] <- NA
                rownames(ret) <- NULL
                return(ret)
            } else {
                tabdata <- tabdata[do.call(order, as.data.frame(tabdata[,1:J])),]
                if(!is.null(printvalue)){
                    if(!is.numeric(printvalue)) stop('printvalue is not a number.')
                    tabdata <- tabdata[abs(tabdata[ ,ncol(tabdata)]) > printvalue, ]
                }                
                return(tabdata)
            }
        }
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'ExploratoryClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45, 
                          which.items = 1:ncol(x@data),
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                          auto.key = TRUE, ...)
    {
        if (!type %in% c('info','infocontour', 'SE', 'trace', 'infotrace', 'infoSE', 'score'))
            stop(type, " is not a valid plot type.")
        if (any(theta_angle > 90 | theta_angle < 0))
            stop('Improper angle specifed. Must be between 0 and 90.')
        if(length(theta_angle) > 1) type = 'infoangle'
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
        nfact <- x@nfact
        if(nfact > 2) stop("Can't plot high dimensional solutions.")
        if(nfact == 1) theta_angle <- 0
        J <- length(x@pars) - 1
        theta <- seq(-4,4,length.out=npts)
        ThetaFull <- Theta <- thetaComb(theta, nfact)
        prodlist <- attr(x@pars, 'prodlist')
        if(length(prodlist) > 0)
            ThetaFull <- prodterms(Theta,prodlist)
        info <- 0
        if(any(sapply(x@pars, is , 'custom')) && type != 'trace')
            stop('Information function for custom classes not available')
        if(all(!sapply(x@pars, is , 'custom'))){
            for(l in 1:length(theta_angle)){
                ta <- theta_angle[l]
                if(nfact == 2) ta <- c(theta_angle[l], 90 - theta_angle[l])
                for(i in 1:J)
                    info <- info + iteminfo(x=x@pars[[i]], Theta=ThetaFull, degrees=ta)
            }
        }
        adj <- apply(x@data, 2, min)
        if(any(adj > 0) && type == 'score')
            message('Adjusted so that the lowest category score for every item is 0')
        tmp <- try(x@rotate, silent = TRUE)
        if (x@nfact > 1 && !is(tmp,'try-error')){
            rotname <- x@rotate
            so <- summary(x, rotate=x@rotate, Target=NULL, verbose=FALSE, digits=5, ...)
            a <- rotateLambdas(so) * 1.702
            for(i in 1:J)
                x@pars[[i]]@par[1:nfact] <- a[i, ]
        }
        itemtrace <- computeItemtrace(x@pars, ThetaFull, x@itemloc)
        score <- c()
        for(i in 1:J)
            score <- c(score, 0:(x@K[i]-1))
        score <- matrix(score, nrow(itemtrace), ncol(itemtrace), byrow = TRUE)
        plt <- data.frame(cbind(info,score=rowSums(score*itemtrace),Theta=Theta))
        if(nfact == 2){
            colnames(plt) <- c("info", "score", "Theta1", "Theta2")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'infocontour')
                return(contourplot(info ~ Theta1 * Theta2, data = plt,
                                   main = paste("Test Information Contour"), xlab = expression(theta[1]),
                                   ylab = expression(theta[2])))
            if(type == 'info')
                return(wireframe(info ~ Theta1 + Theta2, data = plt, main = "Test Information",
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE))
            if(type == 'score')
                return(wireframe(score ~ Theta1 + Theta2, data = plt, main = "Expected Total Score",
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
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
            colnames(plt) <- c("info", "score", "Theta")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'info')
                return(xyplot(info~Theta, plt, type='l',main = 'Test Information',
                              xlab = expression(theta), ylab=expression(I(theta))))
            if(type == 'infocontour')
                cat('No \'contour\' plots for 1-dimensional models\n')
            if(type == 'SE')
                return(xyplot(SE~Theta, plt, type='l',main = 'Test Standard Errors',
                       xlab = expression(theta), ylab=expression(SE(theta))))
            if(type == 'infoSE'){
                obj1 <- xyplot(info~Theta, plt, type='l',main = 'Test Information and Standard Errors',
                               xlab = expression(theta), ylab=expression(I(theta)))
                obj2 <- xyplot(SE~Theta, plt, type='l', ylab=expression(SE(theta)))
                if(!require(latticeExtra)) require(latticeExtra)
                return(doubleYScale(obj1, obj2, add.ylab2 = TRUE))
            }
            if(type == 'trace'){
                P <- vector('list', length(which.items))
                names(P) <- colnames(x@data)[which.items]
                for(i in which.items){
                    tmp <- probtrace(extract.item(x, i), ThetaFull)
                    if(ncol(tmp) == 2L) tmp <- tmp[,2, drop=FALSE]
                    tmp2 <- data.frame(P=as.numeric(tmp), cat=gl(ncol(tmp), k=nrow(Theta), 
                                                           labels=paste0('cat', 1L:ncol(tmp))))
                    P[[i]] <- tmp2
                }
                nrs <- sapply(P, nrow)                
                Pstack <- do.call(rbind, P)
                names <- c()
                for(i in 1L:length(nrs))
                    names <- c(names, rep(names(P)[i], nrs[i]))
                plotobj <- data.frame(Pstack, item=names, Theta=Theta)
                return(xyplot(P ~ Theta, plotobj, group = item, ylim = c(-0.1,1.1),,
                       xlab = expression(theta), ylab = expression(P(theta)),
                       auto.key = auto.key, type = 'l', main = 'Item trace lines', ...))
            }
            if(type == 'infotrace'){                
                I <- matrix(NA, nrow(Theta), J)
                for(i in which.items)
                    I[,i] <- iteminfo(extract.item(x, i), ThetaFull)
                I <- t(na.omit(t(I)))
                items <- gl(n=length(unique(which.items)), k=nrow(Theta), 
                            labels = paste('Item', which.items))
                plotobj <- data.frame(I = as.numeric(I), Theta=Theta, item=items)
                return(xyplot(I ~ Theta, plotobj, group = item, 
                              xlab = expression(theta), ylab = expression(I(theta)),
                              auto.key = auto.key, type = 'l', main = 'Item information trace lines', ...))
            }
            if(type == 'score')
                return(xyplot(score ~ Theta, plt,
                              xlab = expression(theta), ylab = expression(Total(theta)),
                              type = 'l', main = 'Expected Total Score', ...))
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
