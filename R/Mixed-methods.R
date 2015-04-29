# Methods
setMethod(
    f = "print",
    signature = signature(x = 'MixedClass'),
    definition = function(x)
    {
        class(x) <- 'SingleGroupClass'
        print(x)
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
    definition = function(object, digits = 3, verbose = TRUE, ...)
    {
        if(verbose)
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
        if(verbose){
            if(all(dim(out) != 0L)){
                cat('--------------\nFIXED EFFECTS:\n')
                print(round(out, digits))
            }
            cat('\n--------------\nRANDOM EFFECT COVARIANCE(S):\n')
            cat('Correlations on upper diagonal\n')
        }
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
                listnames <- c(listnames, colnames(object@random[[i]]@gdesign)[1L])
            }
        }
        names(rand) <- listnames
        if(verbose){
            cat('\n')
            print(rand, digits)
        }
        if(length(object@lrPars)){
            betas <- object@lrPars@beta
            SE.betas <- matrix(object@lrPars@SEpar, nrow(betas), ncol(betas),
                               dimnames = list(rownames(betas), paste0('Std.Error_', colnames(betas))))
            z <- betas/SE.betas
            colnames(z) <- paste0('z_', colnames(betas))
            keep <- colSums(betas != 0) > 0
            lr.out <- data.frame(betas, SE.betas, z)
            if(verbose){
                cat('--------------\nLATENT REGRESSION FIXED EFFECTS:\n\n')
                print(round(betas[, keep, drop=FALSE], digits))
                cat("\n")
                print(round(data.frame(SE.betas[, keep, drop=FALSE], z[, keep, drop=FALSE]), digits))
            }
        } else lr.out <- NULL
        return(invisible(list(random=rand, fixed=out, lr.out=lr.out)))
    }
)

setMethod(
    f = "coef",
    signature = 'MixedClass',
    definition = function(object, CI = .95, printSE = FALSE, digits = 3, rawug = FALSE, ...)
    {
        if(printSE) rawug <- TRUE
        if(CI >= 1 || CI <= 0)
            stop('CI must be between 0 and 1', call.=FALSE)
        z <- abs(qnorm((1 - CI)/2))
        SEnames <- paste0('CI_', c((1 - CI)/2*100, ((1 - CI)/2 + CI)*100))
        K <- object@K
        J <- length(K)
        nLambdas <- ncol(object@F)
        allPars <- list()
        if(length(object@pars[[1]]@SEpar) > 0){
            if(printSE){
                for(i in 1:(J+1)){
                    allPars[[i]] <- round(matrix(c(object@pars[[i]]@par,
                                                   object@pars[[i]]@SEpar),
                                                 2, byrow = TRUE), digits)
                    rownames(allPars[[i]]) <- c('par', 'SE')
                    colnames(allPars[[i]]) <- names(object@pars[[i]]@est)
                }
            } else {
                for(i in 1:(J+1)){
                    allPars[[i]] <- round(matrix(c(object@pars[[i]]@par,
                                                   object@pars[[i]]@par - z*object@pars[[i]]@SEpar,
                                                   object@pars[[i]]@par + z*object@pars[[i]]@SEpar),
                                                 3, byrow = TRUE), digits)
                    rownames(allPars[[i]]) <- c('par', SEnames)
                    colnames(allPars[[i]]) <- names(object@pars[[i]]@est)
                }
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
        listnames <- c(colnames(object@Data$data), 'GroupPars')
        if(length(object@random) > 0L){
            if(printSE){
                for(i in 1L:length(object@random)){
                    allPars[[length(allPars) + 1L]] <-
                        round(matrix(c(object@random[[i]]@par,
                                       object@random[[i]]@SEpar),
                                     2, byrow = TRUE), digits)
                    rownames(allPars[[length(allPars)]]) <- c('par', 'SE')
                    colnames(allPars[[length(allPars)]]) <- names(object@random[[i]]@est)
                    listnames <- c(listnames, colnames(object@random[[i]]@gdesign)[1L])
                }
            } else {
                for(i in 1L:length(object@random)){
                    allPars[[length(allPars) + 1L]] <-
                        round(matrix(c(object@random[[i]]@par,
                                       object@random[[i]]@par - z*object@random[[i]]@SEpar,
                                       object@random[[i]]@par + z*object@random[[i]]@SEpar),
                                     3, byrow = TRUE), digits)
                    rownames(allPars[[length(allPars)]]) <- c('par', SEnames)
                    colnames(allPars[[length(allPars)]]) <- names(object@random[[i]]@est)
                    listnames <- c(listnames, colnames(object@random[[i]]@gdesign)[1L])
                }
            }
        }
        if(length(object@lrPars)){
            listnames <- c(listnames, 'lr.betas')
            if(!printSE){
                allPars[[length(allPars)+1L]] <- round(matrix(c(object@lrPars@par,
                                                                object@lrPars@par - z*object@lrPars@SEpar,
                                                                object@lrPars@par + z*object@lrPars@SEpar),
                                                              3, byrow = TRUE), digits)
                rownames(allPars[[length(allPars)]]) <- c('par', SEnames)
            } else {
                allPars[[length(allPars)+1L]] <- round(matrix(c(object@lrPars@par,
                                                                object@lrPars@SEpar),
                                                              2, byrow = TRUE), digits)
                rownames(allPars[[length(allPars)]]) <- c('par', 'SE')
            }
            colnames(allPars[[length(allPars)]]) <- names(object@lrPars@est)
        }
        names(allPars) <- listnames
        return(allPars)
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'MixedClass'),
    definition = function(object, object2, verbose = TRUE)
    {
        class(object) <- 'SingleGroupClass'
        anova(object, object2, verbose=verbose)
    }
)

