# Methods
setMethod(
    f = "print",
    signature = signature(x = 'ConfirmatoryClass'),
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
                cat("G2 (", x@df,") = ", round(x@G2,2), ", p = ", round(x@p,4), sep='')
                cat("\nRMSEA = ", round(x@RMSEA,3), ", CFI = ", round(x@CFI,3), 
                    ", TLI = ", round(x@TLI,3), sep='')
            }
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
    definition = function(object, suppress = 0, digits = 3, verbose = TRUE, ...)
    {
        nfact <- ncol(object@F)
        itemnames <- colnames(object@data)
        F <- object@F
        h2 <- as.matrix(object@h2)
        colnames(h2) <- 'h2'
        rownames(F) <- itemnames
        SS <- apply(F^2,2,sum)
        gpars <- ExtractGroupPars(object@pars[[length(object@pars)]])
        Phi <- gpars$gcov
        Phi <- round(Phi, digits)
        colnames(Phi) <- rownames(Phi) <- paste('F',1:ncol(Phi), sep='')
        if(verbose){
            cat("\nFactor loadings metric: \n")
            print(cbind(F, h2),digits)
            cat("\nSS loadings: ",round(SS,digits), "\n")
            cat("\nFactor covariance: \n")
            print(Phi)
        }
        invisible(list(F=F, fcor=Phi))
    }
)

setMethod(
    f = "coef",
    signature = 'ConfirmatoryClass',
    definition = function(object, CI = .95, digits = 3, IRTpars = FALSE, rawug = FALSE, ...)
    {
        if(CI >= 1 || CI <= 0)
            stop('CI must be between 0 and 1')
        z <- abs(qnorm((1 - CI)/2))
        SEnames <- paste0('CI_', c((1 - CI)/2*100, ((1 - CI)/2 + CI)*100))
        K <- object@K
        J <- length(K)
        nLambdas <- ncol(object@F)
        allPars <- list()
        if(IRTpars){
            if(object@nfact > 1L) 
                stop('traditional parameterization is only available for unidimensional models')
            for(i in 1:(J+1))
                allPars[[i]] <- round(mirt2traditional(object@pars[[i]]), digits)
        } else {
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
        }
        if(!rawug){
            allPars <- lapply(allPars, function(x, digits){
                x[ , colnames(x) %in% c('g', 'u')] <- round(antilogit(x[ , colnames(x) %in% c('g', 'u')]), digits)
                x
            },  digits=digits)
        }
        names(allPars) <- c(colnames(object@data), 'GroupPars')
        return(allPars)
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, restype = 'LD', digits = 3, df.p = FALSE, full.scores = FALSE, 
                          printvalue = NULL, verbose = TRUE, ...)
    {
        K <- object@K
        data <- object@data
        N <- nrow(data)
        J <- ncol(data)
        nfact <- ncol(object@F) - length(attr(object@pars, 'prodlist'))
        itemloc <- object@itemloc
        res <- matrix(0,J,J)
        diag(res) <- NA
        colnames(res) <- rownames(res) <- colnames(data)
        theta <- seq(-4,4, length.out = round(20/nfact))
        Theta <- thetaComb(theta,nfact)
        ThetaShort <- Theta
        if(length(object@prodlist) > 0) Theta <- prodterms(Theta, object@prodlist)
        gpars <- ExtractGroupPars(object@pars[[length(object@pars)]])
        prior <- mvtnorm::dmvnorm(ThetaShort,gpars$gmeans,gpars$gcov)
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
    f = "fitted",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, digits = 3, ...){
        tabdata <- object@tabdata
        N <- nrow(object@data)
        expected <- round(N * object@Pl,digits)
        return(cbind(tabdata,expected))
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'ConfirmatoryClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45,
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10), ...)
    {
        class(x) <- 'ExploratoryClass'
        plot(x, type=type, npts=npts, theta_angle=theta_angle, rot=rot, ...)

    }
)

mirt2traditional <- function(x){
    cls <- class(x)
    par <- x@par
    if(cls != 'GroupPars')
        ncat <- x@ncat
    if(cls == 'dich'){
        par[2] <- -par[2]/par[1]
        names(par) <- c('a', 'b', 'g', 'u')
    } else if(cls == 'graded'){
        for(i in 2:ncat)
            par[i] <- -par[i]/par[1]
        names(par) <- c('a', paste0('b', 1:(length(par)-1)))
    } else if(cls == 'gpcm'){
        ds <- par[-1]/par[1]        
        newd <- numeric(length(ds)-1)
        for(i in 1:length(newd))
            newd[i] <- -(ds[i+1] - ds[i])
        par <- c(par[1], newd)
        names(par) <- c('a', paste0('b', 1:length(newd)))
    } else if(cls == 'nominal'){
        as <- par[2:(ncat+1)] * par[1] 
        as <- as - mean(as)
        ds <- par[(ncat+2):length(par)]
        ds <- ds - mean(ds)
        par <- c(as, ds)
        names(par) <- c(paste0('a', 1:ncat), paste0('c', 1:ncat))
    } else if(cls == 'nestlogit'){
        par1 <- par[1:4]
        par1[2] <- -par1[2]/par1[1]
        names(par1) <- c('a', 'b', 'g', 'u')
        par2 <- par[5:length(par)]        
        as <- par2[1:(ncat-1)]
        as <- as - mean(as)
        ds <- par2[-c(1:(ncat-1))]
        ds <- ds - mean(ds)
        names(as) <- paste0('a', 1:(ncat-1))
        names(ds) <- paste0('c', 1:(ncat-1))        
        par <- c(par1, as, ds)
    } else {
        names(par) <- names(x@est)
    }    
    ret <- matrix(par, 1L, dimnames=list('par', names(par)))
    ret
}
