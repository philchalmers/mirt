# Methods
setMethod(
    f = "print",
    signature = signature(x = 'MixedClass'),
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
            cat("AIC =", x@AIC, "\n")
            cat("AICc =", x@AICc, "\n")
            cat("BIC =", x@BIC, "\n")
            cat("SABIC =", x@SABIC, "\n")
        }
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'MixedClass'),
    definition = function(object) {
        print(object)
    }
)

setMethod(
    f = "summary",
    signature = 'MixedClass',
    definition = function(object, digits = 3, ...)
    {
        cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")        
        nbetas <- ncol(object@pars[[1L]]@fixed.design)
        out <- data.frame()
        if(nbetas > 0L){
            out <- data.frame(Estimate=object@pars[[1L]]@par[1L:nbetas], 
                                    'Std.Error'=object@pars[[1L]]@SEpar[1L:nbetas],
                                    row.names=names(object@pars[[1L]]@est[1L:nbetas]))
            out$'z.value' <- out$Estimate / out$'Std.Error'
        }        
        if(all(dim(out) != 0L)){
            cat('--------------\nFIXED EFFECTS:\n')        
            print(round(out, digits))
        }        
        cat('\n--------------\nRANDOM EFFECT COVARIANCE(S):\n')
        cat('Correlations on upper diagonal\n')
        par <- object@pars[[length(object@pars)]]@par[-c(1L:object@nfact)]
        sigma <- matrix(0, object@nfact, object@nfact)
        sigma[lower.tri(sigma, TRUE)] <- par
        if(object@nfact > 1L){           
            sigma <- sigma + t(sigma) - diag(diag(sigma))
            csigma <- cov2cor(sigma)
            sigma[upper.tri(sigma, diag=FALSE)] <- csigma[upper.tri(sigma, diag=FALSE)]
        }
        colnames(sigma) <- rownames(sigma) <- object@factorNames
        rand <- list(sigma)
        listnames <- 'Theta'
        if(length(object@random) > 0L){
            for(i in 1L:length(object@random)){
                par <- object@random[[i]]@par                
                sigma <- matrix(0, object@random[[i]]@ndim, object@random[[i]]@ndim)
                sigma[lower.tri(sigma, TRUE)] <- par                
                if(ncol(sigma) > 1L){
                    sigma <- sigma + t(sigma) - diag(diag(sigma))
                    csigma <- cov2cor(sigma)
                    sigma[upper.tri(sigma, diag=FALSE)] <- csigma[upper.tri(sigma, diag=FALSE)]
                }
                colnames(sigma) <- rownames(sigma) <- 
                    paste0('COV_', colnames(object@random[[i]]@gdesign))
                rand[[length(rand) + 1L]] <- sigma
                listnames <- c(listnames, colnames(object@random[[i]]@gframe)[1L])
            }
        }
        names(rand) <- listnames
        cat('\n')
        print(rand, digits)
        return(invisible(list(random=rand, fixed=out)))
    }
)

setMethod(
    f = "coef",
    signature = 'MixedClass',
    definition = function(object, CI = .95, digits = 3, rawug = FALSE, ...)
    {
        if(CI >= 1 || CI <= 0)
            stop('CI must be between 0 and 1')
        z <- abs(qnorm((1 - CI)/2))
        SEnames <- paste0('CI_', c((1 - CI)/2*100, ((1 - CI)/2 + CI)*100))
        K <- object@K
        J <- length(K)
        nLambdas <- ncol(object@F)
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
        if(!rawug){
            allPars <- lapply(allPars, function(x, digits){
                x[ , colnames(x) %in% c('g', 'u')] <- round(antilogit(x[ , colnames(x) %in% c('g', 'u')]), digits)
                x
            },  digits=digits)
        }
        listnames <- c(colnames(object@data), 'GroupPars')
        if(length(object@random) > 0L){            
            for(i in 1L:length(object@random)){
                allPars[[length(allPars) + 1L]] <- 
                    round(matrix(c(object@random[[i]]@par, 
                                   object@random[[i]]@par - z*object@random[[i]]@SEpar,
                                   object@random[[i]]@par + z*object@random[[i]]@SEpar),
                                 3, byrow = TRUE), digits)
                rownames(allPars[[length(allPars)]]) <- c('par', SEnames)
                colnames(allPars[[length(allPars)]]) <- names(object@random[[i]]@est)
                listnames <- c(listnames, colnames(object@random[[i]]@gframe)[1L])
            }            
        }
        names(allPars) <- listnames        
        return(allPars)
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'MixedClass'),
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
        ret <- data.frame(AIC = c(object@AIC, object2@AIC),
                          AICc = c(object@AICc, object2@AICc),
                          SABIC = c(object@SABIC, object2@SABIC),
                          BIC = c(object@BIC, object2@BIC),
                          logLik = c(object@logLik, object2@logLik),
                          X2 = c('', X2),
                          df = c('', abs(df)),
                          p = c('', round(1 - pchisq(X2,abs(df)),3)))
        return(ret)
    }
)

