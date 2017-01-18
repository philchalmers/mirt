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
    definition = function(object, verbose = TRUE, ...)
    {
        if(verbose)
            cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"),
                "\n\n", sep = "")
        nbetas <- ncol(object@ParObjects$pars[[1L]]@fixed.design)
        out <- data.frame()
        if(nbetas > 0L){
            out <- data.frame(Estimate=object@ParObjects$pars[[1L]]@par[1L:nbetas],
                                    'Std.Error'=object@ParObjects$pars[[1L]]@SEpar[1L:nbetas],
                                    row.names=names(object@ParObjects$pars[[1L]]@est[1L:nbetas]))
            out$'z.value' <- out$Estimate / out$'Std.Error'
        }
        if(verbose){
            if(all(dim(out) != 0L)){
                cat('--------------\nFIXED EFFECTS:\n')
                print(round(out, 3))
            }
            cat('\n--------------\nRANDOM EFFECT COVARIANCE(S):\n')
            cat('Correlations on upper diagonal\n')
        }
        par <- object@ParObjects$pars[[length(object@ParObjects$pars)]]@par[-c(1L:object@Model$nfact)]
        sigma <- matrix(0, object@Model$nfact, object@Model$nfact)
        sigma[lower.tri(sigma, TRUE)] <- par
        if(object@Model$nfact > 1L){
            sigma <- sigma + t(sigma) - diag(diag(sigma))
            csigma <- cov2cor(sigma)
            sigma[upper.tri(sigma, diag=FALSE)] <- csigma[upper.tri(sigma, diag=FALSE)]
        }
        nfact <- extract.mirt(object, 'nfact')
        colnames(sigma) <- rownames(sigma) <- object@Model$factorNames[1L:nfact]
        rand <- list(sigma)
        listnames <- 'Theta'
        if(length(object@ParObjects$random) > 0L){
            for(i in 1L:length(object@ParObjects$random)){
                par <- object@ParObjects$random[[i]]@par
                sigma <- matrix(0, object@ParObjects$random[[i]]@ndim, object@ParObjects$random[[i]]@ndim)
                sigma[lower.tri(sigma, TRUE)] <- par
                if(ncol(sigma) > 1L){
                    sigma <- sigma + t(sigma) - diag(diag(sigma))
                    csigma <- cov2cor(sigma)
                    sigma[upper.tri(sigma, diag=FALSE)] <- csigma[upper.tri(sigma, diag=FALSE)]
                }
                colnames(sigma) <- rownames(sigma) <-
                    paste0('COV_', colnames(object@ParObjects$random[[i]]@gdesign))
                rand[[length(rand) + 1L]] <- sigma
                listnames <- c(listnames, colnames(object@ParObjects$random[[i]]@gdesign)[1L])
            }
        }
        names(rand) <- listnames
        if(verbose){
            cat('\n')
            print(rand, 3)
        }
        if(length(object@Model$lrPars) && any(object@Model$lrPars@est)){
            betas <- object@Model$lrPars@beta
            SE.betas <- matrix(object@Model$lrPars@SEpar, nrow(betas), ncol(betas),
                               dimnames = list(rownames(betas), paste0('Std.Error_', colnames(betas))))
            z <- betas/SE.betas
            colnames(z) <- paste0('z_', colnames(betas))
            keep <- colSums(betas != 0) > 0
            lr.out <- data.frame(betas, SE.betas, z)
            if(verbose){
                cat('--------------\nLATENT REGRESSION FIXED EFFECTS:\n\n')
                print(round(betas[, keep, drop=FALSE], 3))
                cat("\n")
                print(round(data.frame(SE.betas[, keep, drop=FALSE], z[, keep, drop=FALSE]), 3))
            }
        } else lr.out <- NULL
        if(length(object@ParObjects$lr.random) > 0L){
            lr.rand <- list()
            listnames <- NULL
            for(i in 1L:length(object@ParObjects$lr.random)){
                par <- object@ParObjects$lr.random[[i]]@par
                sigma <- matrix(0, object@ParObjects$lr.random[[i]]@ndim,
                                object@ParObjects$lr.random[[i]]@ndim)
                sigma[lower.tri(sigma, TRUE)] <- par
                if(ncol(sigma) > 1L){
                    sigma <- sigma + t(sigma) - diag(diag(sigma))
                    csigma <- cov2cor(sigma)
                    sigma[upper.tri(sigma, diag=FALSE)] <- csigma[upper.tri(sigma, diag=FALSE)]
                }
                colnames(sigma) <- rownames(sigma) <-
                    paste0('COV_', colnames(object@ParObjects$lr.random[[i]]@gdesign))
                lr.rand[[length(lr.rand) + 1L]] <- sigma
                listnames <- c(listnames, colnames(object@ParObjects$lr.random[[i]]@gdesign)[1L])
            }
            names(lr.rand) <- listnames
        } else lr.rand <- NULL
        if(verbose && !is.null(lr.rand)){
            cat('\n--------------\nRANDOM EFFECT COVARIANCE(S):\n')
            cat('Correlations on upper diagonal\n')
            cat('\n')
            print(lr.rand, 3)
        }
        return(invisible(list(random=rand, fixed=out, lr.out=lr.out, lr.rand=lr.rand)))
    }
)

setMethod(
    f = "coef",
    signature = 'MixedClass',
    definition = function(object, CI = .95, printSE = FALSE, rawug = FALSE, ...)
    {
        if(printSE) rawug <- TRUE
        if(CI >= 1 || CI <= 0)
            stop('CI must be between 0 and 1', call.=FALSE)
        z <- abs(qnorm((1 - CI)/2))
        SEnames <- paste0('CI_', c((1 - CI)/2*100, ((1 - CI)/2 + CI)*100))
        K <- object@Data$K
        J <- length(K)
        allPars <- list()
        if(length(object@ParObjects$pars[[1]]@SEpar) > 0){
            if(printSE){
                for(i in 1L:(J+1L)){
                    allPars[[i]] <- matrix(c(object@ParObjects$pars[[i]]@par,
                                                   object@ParObjects$pars[[i]]@SEpar),
                                                 2, byrow = TRUE)
                    rownames(allPars[[i]]) <- c('par', 'SE')
                    colnames(allPars[[i]]) <- names(object@ParObjects$pars[[i]]@est)
                }
            } else {
                for(i in 1L:(J+1L)){
                    allPars[[i]] <- matrix(c(object@ParObjects$pars[[i]]@par,
                                                   object@ParObjects$pars[[i]]@par - z*object@ParObjects$pars[[i]]@SEpar,
                                                   object@ParObjects$pars[[i]]@par + z*object@ParObjects$pars[[i]]@SEpar),
                                                 3, byrow = TRUE)
                    rownames(allPars[[i]]) <- c('par', SEnames)
                    colnames(allPars[[i]]) <- names(object@ParObjects$pars[[i]]@est)
                }
            }
        } else {
            for(i in 1L:(J+1L)){
                allPars[[i]] <- matrix(object@ParObjects$pars[[i]]@par, 1L)
                colnames(allPars[[i]]) <- names(object@ParObjects$pars[[i]]@est)
                rownames(allPars[[i]]) <- 'par'
            }
        }
        if(!rawug){
            allPars <- lapply(allPars, function(x){
                x[ , colnames(x) %in% c('g', 'u')] <- antilogit(x[ , colnames(x) %in% c('g', 'u')])
                x})
        }
        listnames <- c(colnames(object@Data$data), 'GroupPars')
        if(length(object@ParObjects$random) > 0L){
            if(printSE){
                for(i in 1L:length(object@ParObjects$random)){
                    allPars[[length(allPars) + 1L]] <-
                        matrix(c(object@ParObjects$random[[i]]@par,
                                       object@ParObjects$random[[i]]@SEpar),
                                     2, byrow = TRUE)
                    rownames(allPars[[length(allPars)]]) <- c('par', 'SE')
                    colnames(allPars[[length(allPars)]]) <- names(object@ParObjects$random[[i]]@est)
                    listnames <- c(listnames, colnames(object@ParObjects$random[[i]]@gdesign)[1L])
                }
            } else {
                for(i in 1L:length(object@ParObjects$random)){
                    allPars[[length(allPars) + 1L]] <-
                        matrix(c(object@ParObjects$random[[i]]@par,
                                       object@ParObjects$random[[i]]@par - z*object@ParObjects$random[[i]]@SEpar,
                                       object@ParObjects$random[[i]]@par + z*object@ParObjects$random[[i]]@SEpar),
                                     3, byrow = TRUE)
                    rownames(allPars[[length(allPars)]]) <- c('par', SEnames)
                    colnames(allPars[[length(allPars)]]) <- names(object@ParObjects$random[[i]]@est)
                    listnames <- c(listnames, colnames(object@ParObjects$random[[i]]@gdesign)[1L])
                }
            }
        }
        if(length(object@Model$lrPars) && any(object@Model$lrPars@est)){
            listnames <- c(listnames, 'lr.betas')
            if(!printSE){
                allPars[[length(allPars)+1L]] <- matrix(c(object@Model$lrPars@par,
                                                                object@Model$lrPars@par - z*object@Model$lrPars@SEpar,
                                                                object@Model$lrPars@par + z*object@Model$lrPars@SEpar),
                                                              3, byrow = TRUE)
                rownames(allPars[[length(allPars)]]) <- c('par', SEnames)
            } else {
                allPars[[length(allPars)+1L]] <- matrix(c(object@Model$lrPars@par,
                                                                object@Model$lrPars@SEpar),
                                                              2, byrow = TRUE)
                rownames(allPars[[length(allPars)]]) <- c('par', 'SE')
            }
            colnames(allPars[[length(allPars)]]) <- names(object@Model$lrPars@est)
        }
        if(length(object@ParObjects$lr.random) > 0L){
            if(printSE){
                for(i in 1L:length(object@ParObjects$lr.random)){
                    allPars[[length(allPars) + 1L]] <-
                        matrix(c(object@ParObjects$lr.random[[i]]@par,
                                       object@ParObjects$lr.random[[i]]@SEpar),
                                     2, byrow = TRUE)
                    rownames(allPars[[length(allPars)]]) <- c('par', 'SE')
                    colnames(allPars[[length(allPars)]]) <- names(object@ParObjects$lr.random[[i]]@est)
                    listnames <- c(listnames, colnames(object@ParObjects$lr.random[[i]]@gdesign)[1L])
                }
            } else {
                for(i in 1L:length(object@ParObjects$lr.random)){
                    allPars[[length(allPars) + 1L]] <-
                        matrix(c(object@ParObjects$lr.random[[i]]@par,
                                       object@ParObjects$lr.random[[i]]@par - z*object@ParObjects$lr.random[[i]]@SEpar,
                                       object@ParObjects$lr.random[[i]]@par + z*object@ParObjects$lr.random[[i]]@SEpar),
                                     3, byrow = TRUE)
                    rownames(allPars[[length(allPars)]]) <- c('par', SEnames)
                    colnames(allPars[[length(allPars)]]) <- names(object@ParObjects$lr.random[[i]]@est)
                    listnames <- c(listnames, colnames(object@ParObjects$lr.random[[i]]@gdesign)[1L])
                }
            }
        }
        names(allPars) <- listnames
        class(allPars) <- c('mirt_list', 'list')
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

# Methods
setMethod(
    f = "vcov",
    signature = signature(object = 'MixedClass'),
    definition = function(object)
    {
        class(object) <- 'SingleGroupClass'
        vcov(object)
    }
)

setMethod(
    f = "logLik",
    signature = signature(object = 'MixedClass'),
    definition = function(object){
        extract.mirt(object, 'logLik')
    }
)
