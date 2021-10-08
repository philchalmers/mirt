#' Print the model objects
#'
#' Print model object summaries to the console.
#'
#' @param x an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#'
#' @name print-method
#' @aliases print,SingleGroupClass-method
#'   print,MultipleGroupClass-method print,MixedClass-method print,DiscreteClass-method
#'   print,MixtureClass-method
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @docType methods
#' @rdname print-method
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 1)
#' print(x)
#' }
setMethod(
    f = "print",
    signature = signature(x = 'SingleGroupClass'),
    definition = function(x){
        if(!length(x@time)){
            cat('An object of class \"SingleGroupClass\"\n')
            return(invisible(NULL))
        }
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
        cat("Full-information item factor analysis with ", x@Model$nfact, " factor(s).\n", sep="")
        EMquad <- ''
        if(x@Options$method == 'EM') EMquad <- c('\n     using ', x@Options$quadpts, ' quadrature')
        method <- x@Options$method
        if(method == 'MIXED') method <- 'MHRM'
        if(x@OptimInfo$converged)
            cat("Converged within ", x@Options$TOL, ' tolerance after ', x@OptimInfo$iter, ' ',
                method, " iterations.\n", sep = "")
        else
            cat("FAILED TO CONVERGE within ", x@Options$TOL, ' tolerance after ',
                x@OptimInfo$iter, ' ', method, " iterations.\n", sep="")
        cat('mirt version:', as.character(utils::packageVersion('mirt')), '\n')
        cat('M-step optimizer:', x@Options$Moptim, '\n')
        if(method %in% c('EM', 'QMCEM', 'BL', 'MCEM')){
            if(method == 'EM' || method == 'QMCEM')
                cat('EM acceleration:', x@Options$accelerate, '\n')
            if(method == 'EM' || method == 'BL')
                cat('Number of rectangular quadrature:', x@Options$quadpts)
            else if(method == 'QMCEM')
                cat('Number of quasi-Monte Carlo points:', x@Options$quadpts)
            else if(method == 'MCEM')
                cat('Number of Monte Carlo points:', x@Options$quadpts)
            cat('\n')
        }
        dentype <- switch(x@Options$dentype,
               "EH" = "Empirical histogram",
               "EHW" = 'Empirical histogram (scaled)',
               x@Options$dentype)
        cat('Latent density type:', dentype, '\n')
        if(method == 'MHRM')
            cat("Average MH acceptance ratio(s):", paste0(round(x@OptimInfo$aveAR,3), collapse=', '), "\n")
        if(!is.na(x@OptimInfo$secondordertest)){
            cat("\nInformation matrix estimated with method:", x@Options$SE.type)
            cat('\nSecond-order test: model ', if(!x@OptimInfo$secondordertest)
                'is not a maximum or the information matrix is too inaccurate' else
                    'is a possible local maximum', sep = "")
            if(x@OptimInfo$secondordertest)
                cat("\nCondition number of information matrix = ", x@OptimInfo$condnum)
            cat('\n')
        }
        if(length(x@Fit$logLik) > 0){
            if(x@Fit$logPrior != 0){
                cat("\nLog-posterior = ", x@Fit$logLik + x@Fit$logPrior, if(method == 'MHRM')
                    paste(', SE =', round(x@Fit$SElogLik,3)), "\n",sep='')
                cat('Estimated parameters:', length(extract.mirt(x, 'parvec')), '\n')
            } else {
                cat("\nLog-likelihood = ", x@Fit$logLik, if(method == 'MHRM')
                    paste(', SE =', round(x@Fit$SElogLik,3)), "\n",sep='')
                cat('Estimated parameters:', extract.mirt(x, 'nestpars'), '\n')
                cat("AIC = ", x@Fit$AIC, "\n", sep='')
                cat("BIC = ", x@Fit$BIC, "; SABIC = ", x@Fit$SABIC, "\n", sep='')
            }
            if(!is.nan(x@Fit$p)){
                cat("G2 (", x@Fit$df,") = ", round(x@Fit$G2,2), ", p = ", round(x@Fit$p,4), sep='')
                cat("\nRMSEA = ", round(x@Fit$RMSEA,3), ", CFI = ", round(x@Fit$CFI,3),
                    ", TLI = ", round(x@Fit$TLI,3), sep='')
            }
        }
    }
)

#' Show model object
#'
#' Print model object summaries to the console.
#'
#' @param object an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#'
#' @name show-method
#' @aliases show,SingleGroupClass-method
#'   show,MultipleGroupClass-method show,MixedClass-method show,DiscreteClass-method
#'   show,MixtureClass-method
#' @docType methods
#' @rdname show-method
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 1)
#' show(x)
#' }
setMethod(
    f = "show",
    signature = signature(object = 'SingleGroupClass'),
    definition = function(object){
        print(object)
    }
)

#' Summary of model object
#'
#' Transforms coefficients into a standardized factor loading's metric. For \code{MixedClass} objects,
#' the fixed and random coefficients are printed. Note that while the output to the console is rounded
#' to three digits, the returned list of objects is not. For simulations, use
#' \code{output <- summary(mod, verbose = FALSE)} to suppress the console messages.
#'
#' @param object an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#' @param rotate a string indicating which rotation to use for exploratory models, primarily
#'   from the \code{GPArotation} package (see documentation therein).
#'
#'   Rotations currently supported are: \code{'promax'}, \code{'oblimin'}, \code{'varimax'},
#'   \code{'quartimin'}, \code{'targetT'}, \code{'targetQ'}, \code{'pstT'}, \code{'pstQ'},
#'   \code{'oblimax'}, \code{'entropy'}, \code{'quartimax'}, \code{'simplimax'},
#'   \code{'bentlerT'}, \code{'bentlerQ'}, \code{'tandemI'}, \code{'tandemII'},
#'   \code{'geominT'}, \code{'geominQ'}, \code{'cfT'}, \code{'cfQ'}, \code{'infomaxT'},
#'   \code{'infomaxQ'}, \code{'mccammon'}, \code{'bifactorT'}, \code{'bifactorQ'}.
#'
#'   For models that are not exploratory this input will automatically be set to \code{'none'}
#' @param Target a dummy variable matrix indicting a target rotation pattern. This is required for
#'   rotations such as \code{'targetT'}, \code{'targetQ'}, \code{'pstT'}, and \code{'pstQ'}
#' @param suppress a numeric value indicating which (possibly rotated) factor
#'   loadings should be suppressed. Typical values are around .3 in most
#'   statistical software. Default is 0 for no suppression
#' @param verbose logical; allow information to be printed to the console?
#' @param ... additional arguments to be passed
#'
#' @name summary-method
#' @aliases summary,SingleGroupClass-method
#'   summary,MultipleGroupClass-method summary,MixedClass-method summary,DiscreteClass-method
#'   summary,MixtureClass-method
#' @docType methods
#' @export
#' @rdname summary-method
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @seealso \code{\link{coef-method}}
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 2)
#' summary(x)
#' summary(x, rotate = 'varimax')
#'
#' }
setMethod(
    f = "summary",
    signature = 'SingleGroupClass',
    definition = function(object, rotate = 'oblimin', Target = NULL, suppress = 0,
                          verbose = TRUE, ...){
        if (!object@Options$exploratory || rotate == 'none') {
            F <- object@Fit$F
            F[abs(F) < suppress] <- NA
            h2 <- as.matrix(object@Fit$h2)
            SS <- apply(F^2,2,sum)
            gp <- ExtractGroupPars(object@ParObjects$pars[[object@Data$nitems + 1L]])
            Phi <- cov2cor(gp$gcov)
            colnames(h2) <- "h2"
            rownames(Phi) <- colnames(Phi) <- names(SS) <- colnames(F)[seq_len(object@Model$nfact)]
            loads <- cbind(F,h2)
            if(verbose){
                if(object@Options$exploratory)
                    cat("\nUnrotated factor loadings: \n\n")
                print(loads, 3)
                cat("\nSS loadings: ", round(SS, 3), "\n")
                cat("Proportion Var: ",round(SS/nrow(F), 3), "\n")
                cat("\nFactor correlations: \n\n")
                print(round(Phi, 3))
            }
            invisible(list(rotF=F,h2=h2,fcor=Phi))
        } else {
            F <- object@Fit$F
            h2 <- as.matrix(object@Fit$h2)
            colnames(h2) <- "h2"
            rotF <- Rotate(F, rotate, Target = Target, ...)
            SS <- apply(rotF$loadings^2,2,sum)
            L <- rotF$loadings
            L[abs(L) < suppress] <- NA
            loads <- cbind(L,h2)
            Phi <- diag(ncol(F))
            if(!rotF$orthogonal){
                Phi <- rotF$Phi
            }
            colnames(Phi) <- rownames(Phi) <- colnames(F)
            if(verbose){
                cat("\nRotation: ", rotate, "\n")
                cat("\nRotated factor loadings: \n\n")
                print(loads, 3)
                cat("\nRotated SS loadings: ",round(SS,3), "\n")
                cat("\nFactor correlations: \n\n")
                print(round(Phi, 3))
            }
            if(any(h2 > 1))
                warning("Solution has Heywood cases. Interpret with caution.", call.=FALSE)
            invisible(list(rotF=rotF$loadings,h2=h2,fcor=Phi))
        }
    }
)

#' Extract raw coefs from model object
#'
#' Return a list (or data.frame) of raw item and group level coefficients. Note that while
#' the output to the console is rounded to three digits, the returned list of objects is not.
#' Hence, elements from \code{cfs <- coef(mod); cfs[[1]]} will contain the unrounded results (useful
#' for simulations).
#'
#' @param object an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#' @param CI the amount of converged used to compute confidence intervals; default is
#'   95 percent confidence intervals
#' @param IRTpars logical; convert slope intercept parameters into traditional IRT parameters?
#'   Only applicable to unidimensional models. If a suitable ACOV estimate was computed in the fitted
#'   model, and \code{printSE = FALSE}, then suitable CIs will be included based on the delta
#'   method (where applicable)
#' @param rotate see \code{summary} method for details. The default rotation is \code{'none'}
#' @param Target a dummy variable matrix indicting a target rotation pattern
#' @param printSE logical; print the standard errors instead of the confidence intervals? When
#'   \code{IRTpars = TRUE} then the delta method will be used to compute the associated standard errors
#'   from mirt's default slope-intercept form
#' @param as.data.frame logical; convert list output to a data.frame instead?
#' @param simplify logical; if all items have the same parameter names (indicating they are
#'   of the same class) then they are collapsed to a matrix, and a list of length 2 is returned
#'   containing a matrix of item parameters and group-level estimates
#' @param unique return the vector of uniquely estimated parameters
#' @param verbose logical; allow information to be printed to the console?
#' @param rawug logical; return the untransformed internal g and u parameters?
#'   If \code{FALSE}, g and u's are converted with the original format along with delta standard errors
#' @param ... additional arguments to be passed
#'
#' @name coef-method
#' @aliases coef,SingleGroupClass-method
#'   coef,MultipleGroupClass-method coef,MixedClass-method coef,DiscreteClass-method
#'   coef,MixtureClass-method
#' @docType methods
#' @rdname coef-method
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @export
#' @seealso \code{\link{summary-method}}
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' x <- mirt(dat, 1)
#' coef(x)
#' coef(x, IRTpars = TRUE)
#' coef(x, simplify = TRUE)
#'
#' #with computed information matrix
#' x <- mirt(dat, 1, SE = TRUE)
#' coef(x)
#' coef(x, printSE = TRUE)
#' coef(x, as.data.frame = TRUE)
#'
#' #two factors
#' x2 <- mirt(Science, 2)
#' coef(x2)
#' coef(x2, rotate = 'varimax')
#'
#' }
setMethod(
    f = "coef",
    signature = 'SingleGroupClass',
    definition = function(object, CI = .95, printSE = FALSE, rotate = 'none', Target = NULL,
                          IRTpars = FALSE, rawug = FALSE, as.data.frame = FALSE,
                          simplify = FALSE, unique = FALSE, verbose = TRUE, ...){
        dots <- list(...)
        discrete <- ifelse(is.null(dots$discrete), FALSE, TRUE)
        if(unique) return(extract.mirt(object, 'parvec'))
        if(printSE && length(object@ParObjects$pars[[1L]]@SEpar)) rawug <- TRUE
        if(CI >= 1 || CI <= 0)
            stop('CI must be between 0 and 1', call.=FALSE)
        z <- abs(qnorm((1 - CI)/2))
        SEnames <- paste0('CI_', c((1 - CI)/2*100, ((1 - CI)/2 + CI)*100))
        J <- object@Data$nitems
        nfact <- object@Model$nfact + length(object@Model$prodlist)
        a <- matrix(0, J, nfact)
        for(i in 1:J)
            a[i, ] <- ExtractLambdas(object@ParObjects$pars[[i]])

        if (object@Options$exploratory && rotate != 'none'){
            if(verbose) cat("\nRotation: ", rotate, "\n\n")
            so <- summary(object, rotate=rotate, Target=Target, verbose=FALSE, ...)
            a <- rotateLambdas(so) * 1.702
            for(i in 1:J){
                object@ParObjects$pars[[i]]@par[1:nfact] <- a[i, ]
                object@ParObjects$pars[[i]]@SEpar <- numeric(0L)
            }
            object@ParObjects$pars[[J + 1]]@par[-c(1:nfact)] <- so$fcor[lower.tri(so$fcor, TRUE)]
        }
        if(IRTpars){
            if(object@Model$nfact > 1L)
                stop('traditional parameterization is only available for unidimensional models',
                     call.=FALSE)
            vcov <- vcov(object)
            for(i in 1L:J){
                if(class(object@ParObjects$pars[[i]]) %in% c('gpcmIRT')) next
                object@ParObjects$pars[[i]] <- mirt2traditional(object@ParObjects$pars[[i]], vcov=vcov)

            }
        }
        allPars <- list()
        if(length(object@ParObjects$pars[[1L]]@SEpar) && !simplify){
            if(printSE){
                for(i in seq_len(J+1L)){
                    allPars[[i]] <- matrix(c(object@ParObjects$pars[[i]]@par,
                                                   object@ParObjects$pars[[i]]@SEpar),
                                                 2, byrow = TRUE)
                    rownames(allPars[[i]]) <- c('par', 'SE')
                    nms <- names(object@ParObjects$pars[[i]]@est)
                    if(i <= J && object@Model$itemtype[i] != 'custom' && !IRTpars){
                        nms[nms == 'g'] <- 'logit(g)'
                        nms[nms == 'u'] <- 'logit(u)'
                    }
                    colnames(allPars[[i]]) <- nms
                }
            } else {
                for(i in seq_len(J+1L)){
                    allPars[[i]] <- matrix(c(object@ParObjects$pars[[i]]@par,
                                                   object@ParObjects$pars[[i]]@par - z*object@ParObjects$pars[[i]]@SEpar,
                                                   object@ParObjects$pars[[i]]@par + z*object@ParObjects$pars[[i]]@SEpar),
                                                 3, byrow = TRUE)
                    rownames(allPars[[i]]) <- c('par', SEnames)
                    colnames(allPars[[i]]) <- object@ParObjects$pars[[i]]@parnames
                }
            }
        } else {
            for(i in seq_len(J+1L)){
                allPars[[i]] <- matrix(object@ParObjects$pars[[i]]@par, 1L)
                colnames(allPars[[i]]) <- object@ParObjects$pars[[i]]@parnames
                rownames(allPars[[i]]) <- 'par'
            }
        }
        if(!rawug && !IRTpars){
            allPars <- lapply(allPars, function(x, digits){
                x[ , colnames(x) %in% c('g', 'u')] <- antilogit(x[ , colnames(x) %in% c('g', 'u')])
                x
            })
        }
        names(allPars) <- c(colnames(object@Data$data), 'GroupPars')
        if(as.data.frame)
            allPars <- t(as.data.frame(allPars))
        if(simplify && !as.data.frame){
            allPars <- lapply(allPars, function(x) x[1L, , drop=FALSE])
            items.old <- allPars[seq_len(length(allPars)-1L)]
            nms <- lapply(items.old, colnames)
            unms <- unique(do.call(c, nms))
            items <- matrix(NA, length(items.old), length(unms))
            rownames(items) <- names(items.old)
            colnames(items) <- unms
            for(i in seq_len(nrow(items)))
                items[i, nms[[i]]] <- items.old[[i]]
            nfact <- object@Model$nfact
            means <- allPars$GroupPars[seq_len(nfact)]
            if(discrete){
                allPars <- list(items=items, group.intercepts=allPars$GroupPars)
            } else {
                if(object@ParObjects$pars[[J+1L]]@dentype == "Davidian"){
                    covs <- matrix(NA, nfact, nfact)
                    covs[lower.tri(covs, TRUE)] <- allPars$GroupPars[2L]
                    covs[upper.tri(covs, FALSE)] <- covs[lower.tri(covs, FALSE)]
                    colnames(covs) <- rownames(covs) <- names(means) <- object@Model$factorNames[seq_len(nfact)]
                    allPars <- list(items=items, means=means, cov=covs,
                                    Davidian_phis=allPars$GroupPars[-c(1:2)])
                } else {
                    covs <- matrix(NA, nfact, nfact)
                    if(object@ParObjects$pars[[J+1L]]@dentype == "mixture")
                        covs[lower.tri(covs, TRUE)] <- allPars$GroupPars[-c(seq_len(nfact), length(allPars$GroupPars))]
                    else covs[lower.tri(covs, TRUE)] <- allPars$GroupPars[-seq_len(nfact)]
                    covs <- makeSymMat(covs)
                    colnames(covs) <- rownames(covs) <- names(means) <- object@Model$factorNames[seq_len(nfact)]
                    allPars <- list(items=items, means=means, cov=covs)
                }
            }
        }
        if(.hasSlot(object@Model$lrPars, 'beta')){
            allPars$lr.betas <- object@Model$lrPars@beta
            if(!all(is.nan(object@Model$lrPars@SEpar))){
                tmp <- allPars$lr.betas
                if(printSE){
                    tmp[] <- object@Model$lrPars@SEpar
                    allPars$lr.betas <- list(betas=allPars$lr.betas, SE=tmp)
                } else {
                    low <- tmp - z*object@Model$lrPars@SEpar
                    high <- tmp + z*object@Model$lrPars@SEpar
                    allPars$lr.betas <- list(betas=allPars$lr.betas, low, high)
                    names(allPars$lr.betas) <- c('betas', SEnames)
                }
            }
        }
        if(!as.data.frame)
            class(allPars) <- c('mirt_list', 'list')
        return(allPars)
    }
)

#' Compare nested models with likelihood-based statistics
#'
#' Compare nested models using likelihood ratio test (X2), Akaike Information Criterion (AIC),
#' Bayesian Information Criterion (BIC),
#' Sample-Size Adjusted BIC (SABIC), and Hannan-Quinn (HQ) Criterion.
#' When given a sequence of objects, \code{anova} tests the models against one another
#' in the order specified.
#'
#' @param object an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#' @param object2 a second model estimated from any of the mirt package estimation methods
#' @param ... additional model objects to be sequentially compared
#' @param bounded logical; are the two models comparing a bounded parameter (e.g., comparing a single
#'   2PL and 3PL model with 1 df)? If \code{TRUE} then a 50:50 mix of chi-squared distributions
#'   is used to obtain the p-value
#' @param mix proportion of chi-squared mixtures. Default is 0.5
#' @param verbose logical; print additional information to console?
#'
#' @return a \code{data.frame}/\code{mirt_df} object
#'
#' @name anova-method
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @aliases anova,SingleGroupClass-method
#'   anova,MultipleGroupClass-method anova,MixedClass-method anova,DiscreteClass-method
#'   anova,MixtureClass-method
#' @export
#' @docType methods
#' @rdname anova-method
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 1)
#' x2 <- mirt(Science, 2)
#' anova(x, x2)
#'
#' # compare three models sequentially
#' x2 <- mirt(Science, 1, 'gpcm')
#' x3 <- mirt(Science, 1, 'nominal')
#' anova(x, x2, x3)
#'
#' # in isolation
#' anova(x)
#'
#' # with priors on first model
#' model <- "Theta = 1-4
#'           PRIOR = (1-4, a1, lnorm, 0, 10)"
#' xp <- mirt(Science, model)
#' anova(xp, x2)
#' anova(xp)
#'
#' # bounded parameter
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1)
#' mod2 <- mirt(dat, 1, itemtype = c(rep('2PL', 4), '3PL'))
#' anova(mod, mod2) #unbounded test
#' anova(mod, mod2, bounded = TRUE) #bounded
#'
#' # priors
#' model <- 'F = 1-5
#'           PRIOR = (5, g, norm, -1, 1)'
#' mod1b <- mirt(dat, model, itemtype = c(rep('2PL', 4), '3PL'))
#' anova(mod1b)
#'
#' model2 <- 'F = 1-5
#'           PRIOR = (1-5, g, norm, -1, 1)'
#' mod2b <- mirt(dat, model2, itemtype = '3PL')
#' anova(mod1b, mod2b)
#'
#' }
setMethod(
    f = "anova",
    signature = signature(object = 'SingleGroupClass'),
    definition = function(object, object2, ...,
                          bounded = FALSE, mix = 0.5, verbose = TRUE){
        dots <- list(...)
        if(length(dots)){
            dots <- c(object, object2, dots)
            ret <- vector('list', length(dots)-1L)
            for(i in 1L:length(ret)){
                ret[[i]] <- anova(dots[[i]], dots[[i+1L]], bounded=bounded,
                                  mix=mix, verbose=FALSE)
                if(i > 1L)
                    ret[[i]] <- ret[[i]][2L, ]
            }
            ret <- do.call(rbind, ret)
            rownames(ret) <- 1L:nrow(ret)
            return(ret)
        }
        if(missing(object2)){
            hasPriors <- object@Fit$logPrior != 0
            ret <- data.frame(AIC = object@Fit$AIC,
                              SABIC = object@Fit$SABIC,
                              HQ = object@Fit$HQ,
                              BIC = object@Fit$BIC,
                              logLik = object@Fit$logLik)
            if(hasPriors){
                ret <- ret[!(colnames(ret) %in% c('AIC'))]
                ret$logPost = object@Fit$logPrior + object@Fit$logLik
            }
            class(ret) <- c('mirt_df', 'data.frame')
            return(ret)
        }
        df <- object@Fit$df - object2@Fit$df
        if(df < 0){
            temp <- object
            object <- object2
            object2 <- temp
        } else if(df == 0 && !any(object2@Fit$logPrior != 0 || object@Fit$logPrior != 0)){
            if((2*object2@Fit$logLik - 2*object@Fit$logLik) < 0){
                temp <- object
                object <- object2
                object2 <- temp
            }
        }
        if(verbose){
            cat('\nModel 1: ')
            print(object@Call)
            cat('Model 2: ')
            print(object2@Call)
            cat('\n')
        }
        if(any(object2@Fit$logPrior != 0 || object@Fit$logPrior != 0)){
            ret <- data.frame(SABIC = c(object@Fit$SABIC, object2@Fit$SABIC),
                              HQ = c(object@Fit$HQ, object2@Fit$HQ),
                              BIC = c(object@Fit$BIC, object2@Fit$BIC),
                              df = c(NaN, abs(df)),
                              logLik = c(object@Fit$logLik, object2@Fit$logLik),
                              logPost = c(object@Fit$logLik + object@Fit$logPrior,
                                          object2@Fit$logLik + object2@Fit$logPrior))
        } else {
            X2 <- 2*object2@Fit$logLik - 2*object@Fit$logLik
            ret <- data.frame(AIC = c(object@Fit$AIC, object2@Fit$AIC),
                              SABIC = c(object@Fit$SABIC, object2@Fit$SABIC),
                              HQ = c(object@Fit$HQ, object2@Fit$HQ),
                              BIC = c(object@Fit$BIC, object2@Fit$BIC),
                              logLik = c(object@Fit$logLik, object2@Fit$logLik),
                              X2 = c(NaN, X2),
                              df = c(NaN, abs(df)),
                              p = c(NaN, 1 - pchisq(X2,abs(df))))
            if(bounded)
                ret$p[2L] <- 1 - mixX2(X2, df=abs(df), mix=mix)
        }
        class(ret) <- c('mirt_df', 'data.frame')
        ret
    }
)

#' Compute model residuals
#'
#' Return model implied residuals for linear dependencies between items or at the person level.
#' If the latent trait density was approximated (e.g., Davidian curves, Empirical histograms, etc)
#' then passing \code{use_dentype_estimate = TRUE} will use the internally saved quadrature and
#' density components (where applicable).
#'
#' @param object an object of class \code{SingleGroupClass} or
#'   \code{MultipleGroupClass}. Bifactor models are automatically detected and utilized for
#'   better accuracy
#' @param type type of residuals to be displayed.
#'   Can be either \code{'LD'} or \code{'LDG2'} for a local dependence matrix based on the
#'   X2 or G2 statistics (Chen & Thissen, 1997), \code{'Q3'} for the statistic proposed by
#'   Yen (1984), \code{'JSI'} for the jack-knife statistic proposed Edwards et al. (2018),
#'   \code{'exp'} for the expected values for the frequencies of every response pattern,
#'   and \code{'expfull'} for the expected values for every theoretically observable response pattern.
#'   For the 'LD' and 'LDG2' types, the upper diagonal elements represent the standardized
#'   residuals in the form of signed Cramers V coefficients
#' @param tables logical; for LD type, return the observed, expected, and standardized residual
#'   tables for each item combination?
#' @param df.p logical; print the degrees of freedom and p-values?
#' @param full.scores logical; compute relevant statistics
#'  for each subject in the original data?
#' @param printvalue a numeric value to be specified when using the \code{res='exp'}
#'   option. Only prints patterns that have standardized residuals greater than
#'   \code{abs(printvalue)}. The default (NULL) prints all response patterns
#' @param verbose logical; allow information to be printed to the console?
#' @param Theta a matrix of factor scores used for statistics that require empirical estimates (i.e., Q3).
#'   If supplied, arguments typically passed to \code{fscores()} will be ignored and these values will
#'   be used instead
#' @param theta_lim range for the integration grid
#' @param fold logical; apply the sum 'folding' described by Edwards et al. (2018) for the JSI statistic?
#' @param quadpts number of quadrature nodes to use. The default is extracted from model (if available)
#'   or generated automatically if not available
#' @param QMC logical; use quasi-Monte Carlo integration? If \code{quadpts} is omitted the
#'   default number of nodes is 5000
#' @param suppress a numeric value indicating which parameter local dependency combinations
#'   to flag as being too high. Absolute values for the standardized estimates greater than
#'   this value will be returned, while all values less than this value will be set to NA
#' @param technical list of technical arguments when models are re-estimated (see \code{\link{mirt}}
#'   for details)
#' @param ... additional arguments to be passed to \code{fscores()}
#'
#' @name residuals-method
#' @aliases residuals,SingleGroupClass-method residuals,MixtureClass-method
#'   residuals,MultipleGroupClass-method residuals,DiscreteClass-method
#' @docType methods
#' @rdname residuals-method
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Chen, W. H. & Thissen, D. (1997). Local dependence indices for item pairs using item
#' response theory. \emph{Journal of Educational and Behavioral Statistics, 22}, 265-289.
#'
#' Edwards, M. C., Houts, C. R. & Cai, L. (2018). A Diagnostic Procedure to Detect Departures
#' From Local Independence in Item Response Theory Models.
#' \emph{Psychological Methods, 23}, 138-149.
#'
#' Yen, W. (1984). Effects of local item dependence on the fit and equating performance of the three
#' parameter logistic model. \emph{Applied Psychological Measurement, 8}, 125-145.
#' @examples
#'
#' \dontrun{
#'
#' x <- mirt(Science, 1)
#' residuals(x)
#' residuals(x, tables = TRUE)
#' residuals(x, type = 'exp')
#' residuals(x, suppress = .15)
#' residuals(x, df.p = TRUE)
#'
#' # Pearson's X2 estimate for goodness-of-fit
#' full_table <- residuals(x, type = 'expfull')
#' head(full_table)
#' X2 <- with(full_table, sum((freq - exp)^2 / exp))
#' df <- nrow(full_table) - extract.mirt(x, 'nest') - 1
#' p <- pchisq(X2, df = df, lower.tail=FALSE)
#' data.frame(X2, df, p, row.names='Pearson-X2')
#'
#' # above FOG test as a function
#' PearsonX2 <- function(x){
#'    full_table <- residuals(x, type = 'expfull')
#'    X2 <- with(full_table, sum((freq - exp)^2 / exp))
#'    df <- nrow(full_table) - extract.mirt(x, 'nest') - 1
#'    p <- pchisq(X2, df = df, lower.tail=FALSE)
#'    data.frame(X2, df, p, row.names='Pearson-X2')
#' }
#' PearsonX2(x)
#'
#'
#' # extract results manually
#' out <- residuals(x, df.p = TRUE, verbose=FALSE)
#' str(out)
#' out$df.p[1,2]
#'
#' # with and without supplied factor scores
#' Theta <- fscores(x)
#' residuals(x, type = 'Q3', Theta=Theta)
#' residuals(x, type = 'Q3', method = 'ML')
#'
#' # Edwards et al. (2018) JSI statistic
#' N <- 250
#' a <- rnorm(10, 1.7, 0.3)
#' d <- rnorm(10)
#' dat <- simdata(a, d, N=250, itemtype = '2PL')
#'
#' mod <- mirt(dat, 1)
#' residuals(mod, type = 'JSI')
#' residuals(mod, type = 'JSI', fold=FALSE) # unfolded
#'
#' # LD between items 1-2
#' aLD <- numeric(10)
#' aLD[1:2] <- rnorm(2, 2.55, 0.15)
#' a2 <- cbind(a, aLD)
#' dat <- simdata(a2, d, N=250, itemtype = '2PL')
#'
#' mod <- mirt(dat, 1)
#'
#' # JSI executed in parallel over multiple cores
#' mirtCluster()
#' residuals(mod, type = 'JSI')
#'
#' }
setMethod(
    f = "residuals",
    signature = signature(object = 'SingleGroupClass'),
    definition = function(object, type = 'LD', df.p = FALSE, full.scores = FALSE, QMC = FALSE,
                          printvalue = NULL, tables = FALSE, verbose = TRUE, Theta = NULL,
                          suppress = 1, theta_lim = c(-6, 6), quadpts = NULL, fold = TRUE,
                          technical = list(), ...)
    {
        dots <- list(...)
        if(.hasSlot(object@Model$lrPars, 'beta'))
            stop('Latent regression models not yet supported')
        discrete <- use_dentype_estimate <- FALSE
        if(!is.null(dots$use_dentype_estimate))
            use_dentype_estimate <- dots$use_dentype_estimate
        if(!is.null(dots$discrete) || use_dentype_estimate) discrete <- TRUE
        K <- object@Data$K
        data <- object@Data$data
        N <- nrow(data)
        J <- ncol(data)
        nfact <- ncol(object@Fit$F)
        res <- matrix(0,J,J)
        diag(res) <- NA
        colnames(res) <- rownames(res) <- colnames(data)
        if(!is.null(Theta)){
            if(!is.matrix(Theta)) Theta <- matrix(Theta)
            if(nrow(Theta) > nrow(data))
                Theta <- Theta[-extract.mirt(object, 'completely_missing'), , drop=FALSE]
        }
        if(!discrete){
            if(is.null(quadpts)){
                if(QMC) quadpts <- 5000L
                else quadpts <- object@Options$quadpts
            }
            if(is.nan(quadpts))
                quadpts <- select_quadpts(nfact)
            theta <- as.matrix(seq(theta_lim[1L], theta_lim[2L], length.out = quadpts))
            if(type != 'Q3'){
                Theta <- if(QMC)
                    QMC_quad(npts=quadpts, nfact=nfact, lim=theta_lim)
                else {
                    if(nfact > 3L)
                        warning('High-dimensional models should use QMC integration instead',
                                call.=FALSE)
                    thetaComb(theta, nfact)
                }
            } else if(is.null(Theta)){
                Theta <- fscores(object, verbose=FALSE, full.scores=TRUE, leave_missing=TRUE, ...)
            }
        } else {
            Theta <- object@Model$Theta
        }
        itemnames <- colnames(data)
        listtabs <- list()
        calcG2 <- ifelse(type == 'LDG2', TRUE, FALSE)
        if(type %in% c('LD', 'LDG2')){
            if(!discrete){
                groupPars <- ExtractGroupPars(object@ParObjects$pars[[object@Data$nitems + 1L]])
                if(QMC){
                    prior <- rep(1/nrow(Theta), nrow(Theta))
                } else {
                    prior <- mirt_dmvnorm(Theta,groupPars$gmeans, groupPars$gcov)
                    prior <- prior/sum(prior)
                }
            } else {
                prior <- object@Internals$Prior[[1L]]
            }
            df <- (object@Data$K - 1) %o% (object@Data$K - 1)
            diag(df) <- NA
            colnames(df) <- rownames(df) <- colnames(res)
            for(i in seq_len(J)){
                for(j in seq_len(J)){
                    if(i < j){
                        P1 <- ProbTrace(x=object@ParObjects$pars[[i]], Theta=Theta)
                        P2 <- ProbTrace(x=object@ParObjects$pars[[j]], Theta=Theta)
                        tab <- table(data[,i], data[,j], useNA = 'no')
                        Etab <- matrix(0,K[i],K[j])
                        NN <- sum(tab)
                        for(k in seq_len(K[i]))
                            for(m in seq_len(K[j]))
                                Etab[k,m] <- NN * sum(P1[,k] * P2[,m] * prior)
                        s <- try(gamma.cor(tab) - gamma.cor(Etab), TRUE)
                        if(is.nan(s) || is(s, 'try-error')){
                            res[i,j] <- res[j,i] <- NaN
                            next
                        }
                        if(s == 0) s <- 1
                        if(calcG2){
                            tmp <- tab
                            tmp[tab == 0] <- NA
                            res[j,i] <- 2 * sum(tmp * log(tmp/Etab), na.rm=TRUE)
                        } else {
                            res[j,i] <- sum(((tab - Etab)^2)/Etab)
                        }
                        res[i,j] <- sign(s) * sqrt( abs(res[j,i]) / (NN * min(c(K[i],K[j]) - 1L)))
                        df[i,j] <- pchisq(abs(res[j,i]), df=df[j,i], lower.tail=FALSE)
                        if(tables){
                            tmp <- paste0(itemnames[i], '_', itemnames[j])
                            listtabs[[tmp]] <- list(Obs=tab, Exp=Etab, std_res=(tab-Etab)/sqrt(Etab))
                        }
                    }
                }
            }
            if(tables) return(listtabs)
            if(df.p){
                class(df) <- c('mirt_matrix', 'matrix')
                if(verbose){
                    cat("Degrees of freedom (lower triangle) and p-values:\n\n")
                    print(df, ...)
                    cat("\n")
                }
            }
            if(verbose) cat("LD matrix (lower triangle) and standardized values:\n\n")
            class(res) <- c('mirt_matrix', 'matrix')
            if(suppress < 1){
                pick <- abs(res[upper.tri(res)]) < suppress
                res[lower.tri(res)] <- res[upper.tri(res)][pick] <- NA
            }
            if(verbose) print(res, ...)
            if(df.p){
                ret <- list(df, res)
                names(ret) <- c('df.p', type)
                return(invisible(ret))
            }
            return(invisible(res))
        } else if(type == 'exp'){
            r <- object@Data$Freq[[1L]]
            res <- (r - object@Internals$Pl * nrow(object@Data$data)) /
                             sqrt(object@Internals$Pl * nrow(object@Data$data))
            expected <- N * object@Internals$Pl
            tabdata <- object@Data$tabdata
            rownames(tabdata) <- NULL
            ISNA <- is.na(rowSums(tabdata))
            expected[ISNA] <- res[ISNA] <- NA
            tabdata <- data.frame(tabdata,object@Data$Freq[[1L]],expected,res)
            colnames(tabdata) <- c(colnames(object@Data$tabdata),"freq","exp","std.res")
            if(full.scores){
                tabdata[, 'exp'] <- object@Internals$Pl / r * N
                tabdata2 <- object@Data$tabdata
                stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
                sfulldata <- apply(object@Data$data, 1, paste, sep='', collapse = '/')
                scoremat <- tabdata[match(sfulldata, stabdata2), 'exp', drop = FALSE]
                res <- (1-scoremat) / sqrt(scoremat)
                colnames(res) <- 'std.res'
                ret <- cbind(object@Data$data, scoremat, res)
                ret[is.na(rowSums(ret)), c('exp', 'std.res')] <- NA
                rownames(ret) <- NULL
                class(ret) <- c('mirt_df', 'data.frame')
                return(ret)
            } else {
                tabdata <- tabdata[do.call(order, as.data.frame(tabdata[,1:J])),]
                if(!is.null(printvalue)){
                    if(!is.numeric(printvalue)) stop('printvalue is not a number.', call.=FALSE)
                    tabdata <- tabdata[abs(tabdata[ ,ncol(tabdata)]) > printvalue, ]
                }
                class(tabdata) <- c('mirt_df', 'data.frame')
                return(tabdata)
            }
        } else if(type == 'expfull'){
            K <- extract.mirt(object, 'K')
            nitems <- length(K)
            resp <- vector('list', nitems)
            for(i in seq_len(nitems))
                resp[[i]] <- 0L:(K[i]-1L)
            tabdata <- expand.grid(resp)
            rownames(tabdata) <- NULL
            tabdata <- t(t(tabdata) + extract.mirt(object, 'mins'))
            colnames(tabdata) <- colnames(extract.mirt(object, 'data'))
            sv <- mod2values(object)
            itemtype <- extract.mirt(object, 'itemtype')
            nfact <- extract.mirt(object, 'nfact')
            tmpdat <- matrix(0, nrow=2, ncol=nitems)
            colnames(tmpdat) <- colnames(tabdata)
            large <- mirt(tmpdat, nfact, itemtype=itemtype, pars=sv, TOL=NaN, large='return',
                                      technical = list(customK=K))
            large$tabdata <- poly2dich(tabdata)
            large$Freq$all <- rep(1L, nrow(tabdata))
            large$tabdata2 <- matrix(1L)
            full_object <- mirt(tabdata, nfact, itemtype=itemtype, pars=sv, TOL=NaN, large=large)
            Pl <- full_object@Internals$Pl
            r <- integer(length(Pl))
            ro <- object@Data$Freq[[1L]]
            N <- sum(ro)
            for(i in seq_len(length(ro))){
                pick <- colSums(t(tabdata) == object@Data$tabdata[i,]) == ncol(tabdata)
                r[pick] <- ro[i]
            }
            res <- (r - Pl * N) / sqrt(Pl * N)
            expected <- N * Pl
            tabdata <- data.frame(tabdata,r,expected,res)
            colnames(tabdata) <- c(colnames(object@Data$data),"freq","exp","res")
            tabdata <- tabdata[do.call(order, as.data.frame(tabdata[,1:J])),]
            rownames(tabdata) <- 1:length(r)
            if(!is.null(printvalue)){
                if(!is.numeric(printvalue)) stop('printvalue is not a number.', call.=FALSE)
                tabdata <- tabdata[abs(tabdata[ ,ncol(tabdata)]) > printvalue, ]
            }
            class(tabdata) <- c('mirt_df', 'data.frame')
            return(tabdata)
        } else if(type == 'Q3'){
            if(discrete && !use_dentype_estimate)
                stop('residual type not supported for discrete density forms', call.=FALSE)
            dat <- matrix(NA, N, 2L)
            diag(res) <- 1
            for(i in seq_len(J)){
                ei <- extract.item(object, item=i)
                EI <- expected.item(ei, Theta=Theta)
                dat[ ,1L] <- object@Data$data[ ,i] - EI
                for(j in seq_len(J)){
                    if(i < j){
                        ej <- extract.item(object, item=j)
                        EJ <- expected.item(ej, Theta=Theta)
                        dat[,2L] <- object@Data$data[ ,j] - EJ
                        tmpdat <- na.omit(dat)
                        res[i,j] <- res[j,i] <- cor(tmpdat)[1L,2L]
                    }
                }
            }
            if(verbose) cat("Q3 matrix:\n\n")
            if(suppress < 1){
                pick <- abs(res[upper.tri(res)]) < suppress
                res[lower.tri(res)] <- res[upper.tri(res)][pick] <- NA
            }
            class(res) <- c('mirt_matrix', 'matrix')
            if(verbose) print(res, ...)
            return(invisible(res))
        } else if(type == 'JSI'){
            nfact <- extract.mirt(object, 'nfact')
            stopifnot(nfact == 1L)
            nitems <- extract.mirt(object, 'nitems')
            technical$omp <- FALSE
            as_drop <- myLapply(seq_len(nitems), function(item, mod, technical, ...){
                itemtype <- extract.mirt(mod, 'itemtype')[-item]
                tmpdat <- extract.mirt(mod, 'data')[,-item]
                tmpmod <- mirt(tmpdat, 1L, itemtype=itemtype, SE=TRUE, verbose=FALSE,
                               technical=technical, ...)
                ret <- sapply(coef(tmpmod, printSE=TRUE)[1:ncol(tmpdat)], function(x) x[1L:2L, 'a1'])
                ret
            }, mod=object, technical=technical, ...)
            as <- sapply(coef(object)[1:nitems], function(x) x[1L, 'a1'])
            retmat <- matrix(NA, nitems, nitems)
            colnames(retmat) <- rownames(retmat) <- extract.mirt(object, 'itemnames')
            for(i in seq_len(nitems)){
                tmp <- as_drop[[i]]
                pick <- colnames(tmp)
                zs <- (as[pick] - tmp[1L, ]) / tmp[2L, ]
                retmat[i, pick] <- zs
            }
            if(fold) retmat <- retmat + t(retmat)
            class(retmat) <- c('mirt_matrix', 'matrix')
            retmat

        } else {
            stop('specified type does not exist', call.=FALSE)
        }
    }
)

#' Plot various test-implied functions from models
#'
#' Plot various test implied response functions from models estimated in the mirt package.
#'
#' @param x an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{DiscreteClass}
#' @param y an arbitrary missing argument required for \code{R CMD check}
#' @param type type of plot to view. Can be
#'   \describe{
#'     \item{\code{'info'}}{test information function}
#'     \item{\code{'rxx'}}{for the reliability function}
#'     \item{\code{'infocontour'}}{for the test information contours}
#'     \item{\code{'SE'}}{for the test standard error function}
#'     \item{\code{'infotrace'}}{item information traceline plots}
#'     \item{\code{'infoSE'}}{a combined test information and standard error plot}
#'     \item{\code{'trace'}}{item probability traceline plots}
#'     \item{\code{'itemscore'}}{item scoring traceline plots}
#'     \item{\code{'score'}}{expected total score surface}
#'     \item{\code{'scorecontour'}}{expected total score contour plot}
#'     \item{\code{'EAPsum'}}{compares sum-scores to the expected values based on the EAP for sum-scores method (see \code{\link{fscores}})}
#'   }
#'
#'   Note that if \code{dentype = 'empiricalhist'} was used in estimation then
#'   the type \code{'empiricalhist'}
#'   also will be available to generate the empirical histogram plot, and if
#'   \code{dentype = 'Davidian-#'} was used then the type \code{'Davidian'}
#'   will also be available to generate the curve estimates at the quadrature
#'   nodes used during estimation
#' @param drop2 logical; where appropriate, for dichotomous response items drop the lowest category
#'   and provide information pertaining only to the second response option?
#' @param degrees numeric value ranging from 0 to 90 used in \code{plot} to compute angle
#'   for information-based plots with respect to the first dimension.
#'   If a vector is used then a bubble plot is created with the summed information across the angles specified
#'   (e.g., \code{degrees = seq(0, 90, by=10)})
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is
#'   used in conjunction with \code{npts}
#' @param npts number of quadrature points to be used for plotting features.
#'   Larger values make plots look smoother
#' @param MI a single number indicating how many imputations to draw to form bootstrapped confidence
#'   intervals for the selected test statistic. If greater than 0 a plot will be drawn with a shaded
#'   region for the interval
#' @param CI a number from 0 to 1 indicating the confidence interval to select when MI input is
#'   used. Default uses the 95\% confidence (CI = .95)
#' @param rot allows rotation of the 3D graphics
#' @param which.items numeric vector indicating which items to be used when plotting. Default is
#'   to use all available items
#' @param facet_items logical; apply grid of plots across items? If \code{FALSE}, items will be
#'   placed in one plot for each group
#' @param profile logical; provide a profile plot of response probabilities (objects returned from
#'   \code{\link{mdirt}} only)
#' @param auto.key plotting argument passed to \code{\link{lattice}}
#' @param par.strip.text plotting argument passed to \code{\link{lattice}}
#' @param par.settings plotting argument passed to \code{\link{lattice}}
#' @param ehist.cut a probability value indicating a threshold for excluding cases in empirical
#'   histogram plots. Values larger than the default will include more points in the tails of the
#'   plot, potentially squishing the 'meat' of the plot to take up less area than visually desired
#' @param main argument passed to lattice. Default generated automatically
#' @param drape logical argument passed to lattice. Default generated automatically
#' @param colorkey logical argument passed to lattice. Default generated automatically
#' @param add.ylab2 logical argument passed to lattice. Default generated automatically
#' @param ... additional arguments to be passed to lattice
#'
#' @name plot-method
#' @export
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @aliases plot,SingleGroupClass-method
#'   plot,MultipleGroupClass-method plot,SingleGroupClass,missing-method
#'   plot,DiscreteClass,missing-method plot,MixtureClass,missing-method
#' @docType methods
#' @rdname plot-method
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 1, SE=TRUE)
#' plot(x)
#' plot(x, type = 'info')
#' plot(x, type = 'infotrace')
#' plot(x, type = 'infotrace', facet_items = FALSE)
#' plot(x, type = 'infoSE')
#' plot(x, type = 'rxx')
#'
#' # confidence interval plots when information matrix computed
#' plot(x)
#' plot(x, MI=100)
#' plot(x, type='info', MI=100)
#' plot(x, type='SE', MI=100)
#' plot(x, type='rxx', MI=100)
#'
#' # use the directlabels package to put labels on tracelines
#' library(directlabels)
#' plt <- plot(x, type = 'trace')
#' direct.label(plt, 'top.points')
#'
#' set.seed(1234)
#' group <- sample(c('g1','g2'), nrow(Science), TRUE)
#' x2 <- multipleGroup(Science, 1, group)
#' plot(x2)
#' plot(x2, type = 'trace')
#' plot(x2, type = 'trace', which.items = 1:2)
#' plot(x2, type = 'itemscore', which.items = 1:2)
#' plot(x2, type = 'trace', which.items = 1, facet_items = FALSE) #facet by group
#' plot(x2, type = 'info')
#'
#' x3 <- mirt(Science, 2)
#' plot(x3, type = 'info')
#' plot(x3, type = 'SE', theta_lim = c(-3,3))
#'
#' }
setMethod(
    f = "plot",
    signature = signature(x = 'SingleGroupClass', y = 'missing'),
    definition = function(x, y, type = 'score', npts = 200, drop2 = TRUE, degrees = 45,
                          theta_lim = c(-6,6), which.items = 1:extract.mirt(x, 'nitems'),
                          MI = 0, CI = .95, rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                          facet_items = TRUE, main = NULL,
                          drape = TRUE, colorkey = TRUE, ehist.cut = 1e-10, add.ylab2 = TRUE,
                          par.strip.text = list(cex = 0.7),
                          par.settings = list(strip.background = list(col = '#9ECAE1'),
                                              strip.border = list(col = "black")),
                          auto.key = list(space = 'right', points=FALSE, lines=TRUE),
                          profile = FALSE, ...)
    {
        dots <- list(...)
        if(!(type %in% c('info', 'SE', 'infoSE', 'rxx', 'trace', 'score', 'itemscore',
                       'infocontour', 'infotrace', 'scorecontour', 'empiricalhist', 'Davidian',
                       'EAPsum')))
            stop('type supplied is not supported')
        if (any(degrees > 90 | degrees < 0))
            stop('Improper angle specified. Must be between 0 and 90.', call.=FALSE)
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
        nfact <- x@Model$nfact
        if(nfact > 3) stop("Can't plot high dimensional solutions.", call.=FALSE)
        J <- x@Data$nitems
        theta <- seq(theta_lim[1L],theta_lim[2L],length.out=npts/(nfact^2))
        ThetaFull <- Theta <- thetaComb(theta, nfact)
        prodlist <- attr(x@ParObjects$pars, 'prodlist')
        if(all(x@Data$K[which.items] == 2L) && facet_items) auto.key <- FALSE
        if(length(prodlist))
            ThetaFull <- prodterms(Theta,prodlist)
        if(length(degrees) > ncol(ThetaFull)) type <- 'infoangle'
        if(length(degrees) == 1L) degrees <- rep(degrees, ncol(ThetaFull))
        info <- numeric(nrow(ThetaFull))
        if(type %in% c('info', 'infocontour', 'rxx', 'SE', 'infoSE'))
            info <- testinfo(x, ThetaFull, degrees = degrees, which.items=which.items)
        if(type == 'infoangle'){
            for(i in seq_len(length(degrees)))
                info <- info + testinfo(x, ThetaFull, degrees = rep(degrees[i], ncol(ThetaFull)),
                                        which.items=which.items)
        }
        mins <- x@Data$mins
        maxs <- extract.mirt(x, 'K') + mins - 1
        rotate <- if(is.null(dots$rotate)) 'none' else dots$rotate
        if (x@Options$exploratory){
            if(!is.null(dots$rotate)){
                so <- summary(x, verbose=FALSE, digits=5, ...)
                a <- rotateLambdas(so) * 1.702
                for(i in 1:J)
                    x@ParObjects$pars[[i]]@par[1:nfact] <- a[i, ]
            }
        }
        itemtrace <- computeItemtrace(x@ParObjects$pars, ThetaFull, x@Model$itemloc,
                                      CUSTOM.IND=x@Internals$CUSTOM.IND)
        score <- c()
        for(i in 1:J)
            score <- c(score, (0:(x@Data$K[i]-1) + mins[i]) * (i %in% which.items))
        score <- matrix(score, nrow(itemtrace), ncol(itemtrace), byrow = TRUE)
        plt <- data.frame(cbind(info,score=rowSums(score*itemtrace),Theta=Theta))
        bundle <- length(which.items) != J
        gp <- ExtractGroupPars(x@ParObjects$pars[[J+1]])
        if(MI > 0L && nfact == 1L){
            tmpx <- x
            if(!x@Options$SE)
                stop('Must compute an information matrix', call.=FALSE)
            covB <- x@vcov
            pre.ev <- eigen(covB)
            names <- colnames(covB)
            tmp <- lapply(names, function(x, split){
                as.numeric(strsplit(x, split=split)[[1L]][-1L])
            }, split='\\.')
            imputenums <- do.call(c, tmp)
            CIscore <- CIinfo <- CIrxx <- matrix(0, MI, length(plt$score))
            for(i in seq_len(MI)){
                while(TRUE){
                    tmp <- try(imputePars(pars=x@ParObjects$pars, pre.ev=pre.ev,
                                          imputenums=imputenums, constrain=x@Model$constrain),
                               silent=TRUE)
                    if(!is(tmp, 'try-error')) break
                }
                tmpx@ParObjects$pars <- tmp
                gp2 <- ExtractGroupPars(tmp[[J+1]])
                itemtrace <- computeItemtrace(tmpx@ParObjects$pars, ThetaFull, x@Model$itemloc,
                                              CUSTOM.IND=x@Internals$CUSTOM.IND)
                tmpscore <- rowSums(score * itemtrace)
                CIscore[i, ] <- tmpscore
                CIinfo[i, ] <- testinfo(tmpx, ThetaFull)
                CIrxx[i, ] <- CIinfo[i, ] / (CIinfo[i, ] + 1/gp2$gcov[1L,1L])
            }
        }
        mins <- mins[which.items]
        maxs <- maxs[which.items]
        ybump <- (max(maxs) - min(mins))/15
        ybump_full <- (sum(maxs) - sum(mins))/15
        if(type == 'EAPsum'){
            if(is.null(main))
                main <- "Expected vs Observed Sum-Scores"
            fs <- fscores(x, method = 'EAPsum', full.scores=FALSE, verbose=FALSE, ...)
            plt <- with(fs, data.frame(Scores=Sum.Scores, y=c(observed, expected),
                                       group = rep(c('observed', 'expected'), each=nrow(fs))))
            return(xyplot(y~Scores, plt, type='l', main = main, group=plt$group,
                   auto.key=auto.key, xlab = expression(Sum-Score), ylab=expression(n),
                   par.strip.text=par.strip.text, par.settings=par.settings, ...))
        }
        if(nfact == 3){
            colnames(plt) <- c("info", "score", "Theta1", "Theta2", "Theta3")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'infocontour'){
                if(is.null(main)){
                    main <- paste("Test Information Contour")
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(contourplot(info ~ Theta1 * Theta2 | Theta3, data = plt,
                                   main = main, xlab = expression(theta[1]),
                                   ylab = expression(theta[2]),
                                   par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'scorecontour'){
                if(is.null(main)){
                    main <- paste("Expected Score Contour")
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(contourplot(score ~ Theta1 * Theta2 | Theta3, data = plt,
                                   ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                                   main = main, xlab = expression(theta[1]),
                                   ylab = expression(theta[2]), ylim=c(sum(mins)-.1, sum(maxs)+.1),
                                   par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'info'){
                if(is.null(main)){
                    main <- "Test Information"
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(wireframe(info ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'SEcontour'){
                if(is.null(main)){
                    main <- "Test Standard Errors"
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(contourplot(score ~ Theta1 * Theta2 | Theta3, data = plt,
                                          main = main, xlab = expression(theta[1]),
                                          ylab = expression(theta[2]),
                                   par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'score'){
                if(is.null(main)){
                    main <- if(bundle) "Expected Bundle Score" else "Expected Total Score"
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(wireframe(score ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                                 ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings,
                                 ylim=c(sum(mins)-.1, sum(maxs)+.1), ...))
            } else if(type == 'SE'){
                if(is.null(main)){
                    main <- "Test Standard Errors"
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(wireframe(SE ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                                 zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                stop('plot type not supported for three dimensional model', call.=FALSE)
            }
        } else if(nfact == 2){
            colnames(plt) <- c("info", "score", "Theta1", "Theta2")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'infocontour'){
                if(is.null(main)){
                    main <- paste("Test Information Contour")
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(contourplot(info ~ Theta1 * Theta2, data = plt,
                                   main = main, xlab = expression(theta[1]),
                                   ylab = expression(theta[2]),
                                   par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'scorecontour'){
                if(is.null(main)){
                    main <- paste("Expected Score Contour")
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(contourplot(score ~ Theta1 * Theta2, data = plt,
                                   ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                                   main = main, xlab = expression(theta[1]),
                                   ylab = expression(theta[2]),
                                   par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'info'){
                if(is.null(main)){
                    main <- "Test Information"
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(wireframe(info ~ Theta1 + Theta2, data = plt, main = main,
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'trace'){
                if(is.null(main))
                    main <- 'Item Probability Functions'
                P <- vector('list', length(which.items))
                names(P) <- colnames(x@Data$data)[which.items]
                ind <- 1L
                for(i in which.items){
                    tmp <- probtrace(extract.item(x, i), ThetaFull)
                    if(ncol(tmp) == 2L && drop2) tmp <- tmp[,2, drop=FALSE]
                    tmp2 <- data.frame(P=as.numeric(tmp), cat=gl(ncol(tmp), k=nrow(ThetaFull),
                                                                 labels=paste0('P', seq_len(ncol(tmp)))))
                    P[[ind]] <- tmp2
                    ind <- ind + 1L
                }
                nrs <- sapply(P, nrow)
                Pstack <- do.call(rbind, P)
                names <- c()
                for(i in seq_len(length(nrs)))
                    names <- c(names, rep(names(P)[i], nrs[i]))
                plotobj <- data.frame(Pstack, item=names, Theta=ThetaFull)
                plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
                return(wireframe(P ~ Theta.1 * Theta.2 | item, plotobj, zlim = c(-0.1,1.1), group=plotobj$cat,
                                 zlab=expression(P(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'infotrace'){
                if(is.null(main))
                    main <- 'Item Information'
                I <- matrix(NA, nrow(Theta), J)
                for(i in which.items)
                    I[,i] <- iteminfo(extract.item(x, i), ThetaFull, degrees=degrees)
                I <- t(na.omit(t(I)))
                items <- rep(colnames(x@Data$data)[which.items], each=nrow(Theta))
                plotobj <- data.frame(I = as.numeric(I), Theta=ThetaFull, item=items)
                plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
                return(wireframe(I ~ Theta.1 * Theta.2 | item, plotobj, group=plotobj$cat,
                                 zlab=expression(P(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'itemscore'){
                if(is.null(main))
                    main <- 'Expected Item Score'
                S <- vector('list', length(which.items))
                names(S) <- colnames(x@Data$data)[which.items]
                ind <- 1L
                for(i in which.items){
                    S[[ind]] <- expected.item(extract.item(x, i), ThetaFull, mins[i])
                    ind <- ind + 1L
                }
                Sstack <- do.call(c, S)
                names <- rep(names(S), each = nrow(ThetaFull))
                plotobj <- data.frame(S=Sstack, item=names, Theta=ThetaFull)
                plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
                return(wireframe(S ~ Theta.1 * Theta.2 | item, plotobj, group=plotobj$cat,
                                 zlab=expression(P(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'SEcontour'){
                if(is.null(main)){
                    main <- "Test Standard Errors"
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(contourplot(score ~ Theta1 * Theta2, data = plt,
                                          main = main, xlab = expression(theta[1]),
                                          ylab = expression(theta[2]),
                                   par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'score'){
                if(is.null(main)){
                    main <- if(bundle) "Expected Bundle Score" else "Expected Total Score"
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(wireframe(score ~ Theta1 + Theta2, data = plt, main = main,
                                 zlim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'infoangle'){
                if(is.null(main)){
                    main <- 'Information across different angles'
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                graphics::symbols(plt$Theta1, plt$Theta2, circles = sqrt(plt$info/pi), inches = .35, fg='white', bg='blue',
                        xlab = expression(theta[1]), ylab = expression(theta[2]),
                        main = main)
            } else if(type == 'SE'){
                if(is.null(main)){
                    main <- "Test Standard Errors"
                    if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
                }
                return(wireframe(SE ~ Theta1 + Theta2, data = plt, main = main,
                                 zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                                 par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                stop('plot type not supported for two dimensional model', call.=FALSE)
            }
        } else {
            colnames(plt) <- c("info", "score", "Theta")
            plt$SE <- 1 / sqrt(plt$info)
            plt$rxx <- plt$info / (plt$info + 1/gp$gcov[1L,1L])
            if(MI > 0){
                bs_range <- function(x, CI){
                    ss <- sort(x)
                    N <- length(ss)
                    ret <- c(upper = ss[ceiling(N * (1 - (1-CI)/2))],
                             middle = median(x),
                             lower = ss[floor(N * (1-CI)/2)])
                    ret
                }
                tmp <- apply(CIscore, 2, bs_range, CI=CI)
                plt$CIscoreupper <- tmp['upper', ]
                plt$CIscorelower <- tmp['lower', ]
                tmp <- apply(CIinfo, 2, bs_range, CI=CI)
                plt$CIinfoupper <- tmp['upper', ]
                plt$CIinfolower <- tmp['lower', ]
                plt$CISElower <- 1/sqrt(tmp['upper', ])
                plt$CISEupper <- 1/sqrt(tmp['lower', ])
                tmp <- apply(CIrxx, 2, bs_range, CI=CI)
                plt$CIrxxupper <- tmp['upper', ]
                plt$CIrxxlower <- tmp['lower', ]
            }
            if(type == 'info'){
                if(is.null(main))
                    main <- 'Test Information'
                if(MI > 0){
                    return(xyplot(info ~ Theta, data=plt,
                                  upper=plt$CIinfoupper, lower=plt$CIinfolower,
                                  panel = function(x, y, lower, upper, ...){
                                      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                    col="#E6E6E6", border = FALSE, ...)
                                      panel.xyplot(x, y, type='l', lty=1,...)
                                  },
                                  main = main, ylim=c(min(plt$CIinfolower), max(plt$CIinfoupper)),
                                  ylab = expression(I(theta)), xlab = expression(theta),
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(info~Theta, plt, type='l', main = main,
                                  xlab = expression(theta), ylab=expression(I(theta)),
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            } else if(type == 'rxx'){
                if(is.null(main))
                    main <- 'Reliability'
                if(MI > 0){
                    return(xyplot(rxx ~ Theta, data=plt,
                                  upper=plt$CIrxxupper, lower=plt$CIrxxlower,
                                  panel = function(x, y, lower, upper, ...){
                                      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                    col="#E6E6E6", border = FALSE, ...)
                                      panel.xyplot(x, y, type='l', lty=1,...)
                                  },
                                  main = main, ylim=c(-0.1, 1.1),
                                  ylab = expression(r[xx](theta)), xlab = expression(theta),
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(rxx~Theta, plt, type='l', main = main, ylim=c(-0.1, 1.1),
                                  xlab = expression(theta), ylab=expression(r[xx](theta)),
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            } else if(type == 'SE'){
                if(is.null(main))
                    main <- 'Test Standard Errors'
                if(MI > 0){
                    return(xyplot(SE ~ Theta, data=plt,
                                  upper=plt$CISEupper, lower=plt$CISElower,
                                  panel = function(x, y, lower, upper, ...){
                                      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                    col="#E6E6E6", border = FALSE, ...)
                                      panel.xyplot(x, y, type='l', lty=1,...)
                                  },
                                  main = main, ylim=c(min(plt$CISElower), max(plt$CISEupper)),
                                  ylab = expression(I(theta)), xlab = expression(theta),
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(SE~Theta, plt, type='l', main = main,
                           xlab = expression(theta), ylab=expression(SE(theta)),
                           par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            } else if(type == 'infoSE'){
                if(is.null(main))
                    main <- 'Test Information and Standard Errors'
                par.settings <- c(par.settings,
                                  lattice::simpleTheme(col = c("#0080ff",'red'), lty = 1:2))
                obj1 <- xyplot(info~Theta, plt, type='l', main = main,
                               xlab = expression(theta), ylab=expression(I(theta)),
                               par.strip.text=par.strip.text, par.settings=par.settings)
                obj2 <- xyplot(SE~Theta, plt, type='l', ylab=expression(SE(theta)),
                               par.strip.text=par.strip.text, par.settings=par.settings)
                if(requireNamespace("latticeExtra", quietly = TRUE)){
                    return(latticeExtra::doubleYScale(obj1, obj2, add.ylab2 = add.ylab2, ...))
                } else {
                    stop('latticeExtra package is not available. Please install.', call.=FALSE)
                }
            } else if(type == 'trace'){
                if(is.null(main))
                    main <- 'Item Probability Functions'
                P <- vector('list', length(which.items))
                names(P) <- colnames(x@Data$data)[which.items]
                ind <- 1L
                alltwocats <- all(extract.mirt(x, 'K')[which.items] == 2L)
                for(i in which.items){
                    tmp <- probtrace(extract.item(x, i), ThetaFull)
                    if(ncol(tmp) == 2L && (facet_items || (!facet_items && alltwocats)) && drop2)
                        tmp <- tmp[,2, drop=FALSE]
                    tmp2 <- data.frame(P=as.numeric(tmp), cat=gl(ncol(tmp), k=nrow(Theta),
                                                           labels=paste0('P', seq_len(ncol(tmp)))))
                    P[[ind]] <- tmp2
                    ind <- ind + 1L
                }
                nrs <- sapply(P, nrow)
                Pstack <- do.call(rbind, P)
                names <- c()
                for(i in seq_len(length(nrs)))
                    names <- c(names, rep(names(P)[i], nrs[i]))
                plotobj <- data.frame(Pstack, item=names, Theta=Theta)
                plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
                if(facet_items){
                    return(xyplot(P ~ Theta|item, plotobj, ylim = c(-0.1,1.1), groups = cat,
                                  xlab = expression(theta), ylab = expression(P(theta)),
                                  auto.key = auto.key, type = 'l', main = main,
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(P ~ Theta|cat, plotobj, groups=plotobj$item, ylim = c(-0.1,1.1),
                                  xlab = expression(theta), ylab = expression(P(theta)),
                                  auto.key = auto.key, type = 'l', main = main,
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            } else if(type == 'itemscore'){
                if(is.null(main))
                    main <- 'Expected Item Score'
                S <- vector('list', length(which.items))
                names(S) <- colnames(x@Data$data)[which.items]
                ind <- 1L
                for(i in which.items){
                    S[[ind]] <- expected.item(extract.item(x, i), ThetaFull, mins[i])
                    ind <- ind + 1L
                }
                Sstack <- do.call(c, S)
                names <- rep(names(S), each = nrow(ThetaFull))
                plotobj <- data.frame(S=Sstack, item=names, Theta=Theta)
                plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
                if(facet_items){
                    return(xyplot(S ~ Theta|item, plotobj, ylim=c(min(mins)-ybump, max(maxs)+ybump),
                                  xlab = expression(theta), ylab = expression(S(theta)),
                                  auto.key = auto.key, type = 'l', main = main,
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(S ~ Theta, plotobj, groups=plotobj$item, ylim=c(min(mins)-.1, max(maxs)+.1),
                                  xlab = expression(theta), ylab = expression(S(theta)),
                                  auto.key = auto.key, type = 'l', main = main,
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            } else if(type == 'infotrace'){
                if(is.null(main))
                    main <- 'Item Information'
                I <- matrix(NA, nrow(Theta), J)
                for(i in which.items)
                    I[,i] <- iteminfo(extract.item(x, i), ThetaFull)
                I <- t(na.omit(t(I)))
                items <- rep(colnames(x@Data$data)[which.items], each=nrow(Theta))
                plotobj <- data.frame(I = as.numeric(I), Theta=Theta, item=items)
                plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
                if(facet_items){
                    return(xyplot(I ~ Theta|item, plotobj,
                                  xlab = expression(theta), ylab = expression(I(theta)),
                                  auto.key = auto.key, type = 'l', main = main,
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(I ~ Theta, plotobj, groups = plotobj$item,
                                  xlab = expression(theta), ylab = expression(I(theta)),
                                  auto.key = auto.key, type = 'l', main = main,
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            } else if(type == 'score'){
                if(is.null(main))
                    main <- if(bundle) "Expected Bundle Score" else "Expected Total Score"
                if(MI > 0){
                    return(xyplot(score ~ Theta, data=plt,
                                  ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                                  upper=plt$CIscoreupper, lower=plt$CIscorelower,
                                  panel = function(x, y, lower, upper, ...){
                                      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                    col="#E6E6E6", border = FALSE, ...)
                                      panel.xyplot(x, y, type='l', lty=1,...)
                                  },
                                  main = main,
                                  ylab = expression(T(theta)), xlab = expression(theta),
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(score ~ Theta, plt, ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                                  xlab = expression(theta), ylab = expression(T(theta)),
                                  type = 'l', main = main,
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            } else if(type == 'empiricalhist'){
                if(!(x@Options$dentype %in% c('EH', 'EHW')))
                    stop('Empirical histogram was not estimated for this object', call.=FALSE)
                if(is.null(main))
                    main <- 'Empirical Histogram'
                Prior <- x@Internals$Prior[[1L]]
                Theta <- x@Model$Theta
                # Prior <- Prior * nrow(x@Data$data)
                cuts <- cut(Theta, floor(npts/2))
                Prior <- do.call(c, lapply(split(Prior, cuts), mean))
                Theta <- do.call(c, lapply(split(Theta, cuts), mean))
                keep1 <- min(which(Prior > ehist.cut))
                keep2 <- max(which(Prior > ehist.cut))
                plt <- data.frame(Theta = Theta, Prior = Prior)
                plt <- plt[keep1:keep2, , drop=FALSE]
                return(xyplot(Prior ~ Theta, plt,
                              xlab = expression(theta), ylab = 'Density',
                              type = 'b', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else if(type == 'Davidian'){
                if(x@Options$dentype != 'Davidian')
                   stop('Davidian curve was not estimated for this object', call.=FALSE)
                if(is.null(main))
                    main <- 'Davidian Curve'
                Prior <- x@Internals$Prior[[1L]]
                Theta <- x@Model$Theta
                # Prior <- Prior * nrow(x@Data$data)
                plt <- data.frame(Theta = Theta, Prior = Prior)
                return(xyplot(Prior ~ Theta, plt,
                              xlab = expression(theta), ylab = 'Density',
                              type = 'b', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                stop('plot not supported for unidimensional models', call.=FALSE)
            }
        }
    }
)

mirt2traditional <- function(x, vcov){
    cls <- class(x)
    opar <- par <- x@par
    if(cls != 'GroupPars')
        ncat <- x@ncat
    if(cls == 'dich'){
        fns <- vector('list', 4L)
        fns[[2]] <- function(par, index, opar){
            if(index == 2L){
                opar[1L:2L] <- par
                ret <- -opar[2L]/opar[1L]
            }
            ret
        }
        fns[[3]] <- function(par, index, opar){
            if(index == 3L)
                ret <- plogis(par)
            ret
        }
        fns[[4]] <- function(par, index, opar){
            if(index == 4L)
                ret <- plogis(par)
            ret
        }
        delta_index <- list(NA, 1L:2L, 3L, 4L)
        par[2] <- -par[2]/par[1]
        par[3] <- plogis(par[3])
        par[4] <- plogis(par[4])
        names(par) <- c('a', 'b', 'g', 'u')
    } else if(cls == 'graded'){
        fns <- vector('list', ncat+1L)
        for(i in 2L:ncat){
            fns[[i]] <- function(par, index, opar){
                if(index > 1L){
                    opar[c(1L, index)] <- par
                    ret <- -opar[index]/opar[1L]
                }
                ret
            }
        }
        delta_index <- vector('list', ncat)
        delta_index[[1L]] <- NA
        for(i in 2:ncat){
            par[i] <- -par[i]/par[1]
            delta_index[[i]] <- c(1L, i)
        }
        names(par) <- c('a', paste0('b', 1:(length(par)-1)))
    } else if(cls == 'gpcm'){
        fns <- vector('list', ncat+1L)
        for(i in 2L:ncat){
            fns[[i]] <- function(par, index, opar){
                if(index > 1L){
                    if(index == 2L) opar[c(1, ncat + 3)] <- par
                    else opar[c(1, ncat + index, ncat + index + 1)] <- par
                    par <- opar
                    ds <- par[-1]/par[1]
                    ds <- ds[-seq_len(ncat)]
                    newd <- numeric(length(ds)-1L)
                    for(i in 2:length(ds))
                        newd[i-1L] <- -(ds[i] - ds[i-1L])
                    ret <- c(par[1], newd)
                    ret <- ret[index]
                }
                ret
            }
        }
        delta_index <- vector('list', ncat)
        delta_index[[1L]] <- NA
        ds <- par[-1]/par[1]
        ds <- ds[-seq_len(ncat)]
        newd <- numeric(length(ds)-1L)
        tmp <- rbind(1:(ncat-1), 2:ncat)
        for(i in 2:length(ds)){
            newd[i-1L] <- -(ds[i] - ds[i-1L])
            delta_index[[i]] <- c(1L, tmp[,i-1L] + ncat + 1L)
        }
        delta_index[[2L]] <- c(1L, ncat + 3L)
        par <- c(par[1], newd)
        names(par) <- c('a', paste0('b', 1:length(newd)))
        x@est <- x@est[c(1, (ncat+3L):length(x@est))]
    } else if(cls == 'nominal'){
        as <- par[2:(ncat+1)] * par[1]
        as <- as - mean(as)
        ds <- par[(ncat+2):length(par)]
        ds <- ds - mean(ds)
        par <- c(as, ds)
        names(par) <- c(paste0('a', 1:ncat), paste0('c', 1:ncat))
        x@est <- x@est[-2]
    } else if(cls == 'nestlogit'){
        par1 <- par[1:4]
        par1[2] <- -par1[2]/par1[1]
        par1[3] <- plogis(par1[3])
        par1[4] <- plogis(par1[4])
        names(par1) <- c('a', 'b', 'g', 'u')
        par2 <- par[5:length(par)]
        as <- par2[1:(ncat-1)]
        as <- as - mean(as)
        ds <- par2[-c(1:(ncat-1))]
        ds <- ds - mean(ds)
        names(as) <- paste0('a', 1:(ncat-1))
        names(ds) <- paste0('c', 1:(ncat-1))
        par <- c(par1, as, ds)
    }
    x@par <- par
    names(x@est) <- names(par)
    x@parnames <- names(x@par)
    if(length(vcov) == 0L || (is.na(vcov[1L,1L]) || !(cls %in% c('dich', 'graded', 'gpcm')))){
        x@SEpar <- numeric()
    } else {
        nms <- colnames(vcov)
        splt <- strsplit(nms, "\\.")
        splt <- lapply(splt, function(x) as.integer(x[-1L]))
        for(i in 1L:length(delta_index)){
            if(!x@est[i]) next
            if(is.na(delta_index[[i]][1L])) next
            grad <- numerical_deriv(opar[delta_index[[i]]], fns[[i]], opar=opar, index=i)
            parnum <- x@parnum[delta_index[[i]]]
            pick <- numeric(length(grad))
            for(j in 1L:length(parnum))
                pick[j] <- suppressWarnings(min(which(do.call(c, lapply(splt, function(x)
                    any(x %in% parnum[j]))))))
            grad <- grad[is.finite(pick)]
            pick <- pick[is.finite(pick)]
            x@SEpar[i] <- as.vector(sqrt(grad %*% vcov[pick, pick, drop=FALSE] %*% grad))
        }
        if(cls == 'gpcm') x@SEpar <- x@SEpar[1L:length(x@par)]
    }
    x
}

#' Convert traditional IRT metric into slope-intercept form used in mirt
#'
#' This is a helper function for users who have previously available traditional/classical
#' IRT parameters and want to know the equivalent slope-intercept translation used in \code{mirt}.
#' Note that this function assumes that the supplied models are unidimensional by definition (i.e.,
#' will have only one slope/discrimination). If there is no supported slope-interecept transformation
#' available then the original vector of parameters will be returned by default.
#'
#' Supported class transformations for the \code{cls} input are:
#'
#' \describe{
#'   \item{Rasch, 2PL, 3PL, 3PLu, 4PL}{
#'     Form must be: (discrimination, difficulty, lower-bound, upper-bound)
#'     }
#'   \item{graded}{
#'     Form must be: (discrimination, difficulty 1, difficulty 2, ..., difficulty k-1)
#'     }
#'   \item{gpcm}{
#'     Form must be: (discrimination, difficulty 1, difficulty 2, ..., difficulty k-1)
#'     }
#'   \item{nominal}{
#'     Form must be: (discrimination 1, discrimination 2, ..., discrimination k,
#'       difficulty 1, difficulty 2, ..., difficulty k)
#'     }
#' }
#'
#' @param x a vector of parameters to tranform
#' @param cls the class or itemtype of the supplied model
#' @param ncat the number of categories implied by the IRT model
#'
#' @return a named vector of slope-intercept parameters (if supported)
#' @export
#'
#' @examples
#'
#' # classical 3PL model
#' vec <- c(a=1.5, b=-1, g=.1, u=1)
#' slopeint <- traditional2mirt(vec, '3PL', ncat=2)
#' slopeint
#'
#' # classical graded model (four category)
#' vec <- c(a=1.5, b1=-1, b2=0, b3=1.5)
#' slopeint <- traditional2mirt(vec, 'graded', ncat=4)
#' slopeint
#'
#' # classical generalize partial credit model (four category)
#' vec <- c(a=1.5, b1=-1, b2=0, b3=1.5)
#' slopeint <- traditional2mirt(vec, 'gpcm', ncat=4)
#' slopeint
#'
#' # classical nominal model (4 category)
#' vec <- c(a1=.5, a2 = -1, a3=1, a4=-.5, d1=1, d2=-1, d3=-.5, d4=.5)
#' slopeint <- traditional2mirt(vec, 'nominal', ncat=4)
#' slopeint
#'
#'
traditional2mirt <- function(x, cls, ncat){
    cls <- toInternalItemtype(cls)
    if(cls == 'dich'){
        par <- x
        par[2L] <- -par[2L]*par[1L]
        names(par) <- c('a1', 'd', 'g', 'u')
    } else if(cls == 'graded'){
        par <- x
        for(i in 2L:ncat)
            par[i] <- -par[i]*par[1L]
        names(par) <- c('a1', paste0('d', 1:(length(par)-1)))
    } else if(cls %in% c('gpcm', 'gpcmIRT')){
        if(cls == 'gpcmIRT'){
            x[-c(1, length(x))] <- x[-c(1, length(x))] - x[length(x)]
            x <- x[-length(x)]
        }
        par <- c(x[1L], 0L:(ncat-1L), 0, x[-1L])
        ds <- -par[-c(1:(ncat+1))]*par[1]
        newd <- numeric(length(ds))
        for(i in length(ds):2L)
            newd[i] <- (ds[i] + ds[i-1L])
        for(i in length(newd):3L){
            pick <- if(i %% 2 == 0) seq(2, i, by = 2) else seq(1, i, by = 2)
            newd[i] <- sum(newd[pick])
        }
        par <- c(par[1:(ncat+1)], newd)
        names(par) <- c('a1', paste0('ak', 0:(ncat-1)), paste0('d', 0:(ncat-1)))
    } else if(cls == 'nominal'){
        as <- x[seq_len(length(x)/2)]
        ds <- x[-seq_len(length(x)/2)]
        a1 <- (as[ncat] - as[1L]) / (ncat-1L)
        ak <- 1:ncat - 1
        for(i in 2:(ncat-1))
            ak[i] <- -(as[1L] - as[i]) / a1
        dk <- ak
        for(i in 2:ncat)
            dk[i] <- ds[i] - ds[1L]
        par <- c(a1, ak, dk)
        names(par) <- c('a1', paste0('ak', 0:(ncat-1)), paste0('d', 0:(ncat-1)))
    } else {
        stop('traditional2mirt item class not supported', call.=FALSE)
    }
    par
}

#' Extract parameter variance covariance matrix
#'
#' Extract parameter variance covariance matrix
#'
#' @param object an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#'
#' @name vcov-method
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @export
#' @aliases vcov,SingleGroupClass-method vcov,MixtureClass-method
#'   vcov,MultipleGroupClass-method vcov,MixedClass-method vcov,DiscreteClass-method
#' @docType methods
#' @rdname vcov-method
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 1, SE=TRUE)
#' vcov(x)
#'
#' }
setMethod(
    f = "vcov",
    signature = signature(object = 'SingleGroupClass'),
    definition = function(object){
        extract.mirt(object, 'vcov')
    }
)

#' Extract log-likelihood
#'
#' Extract the observed-data log-likelihood.
#'
#' @param object an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#'
#' @name logLik-method
#' @aliases logLik,SingleGroupClass-method logLik,MixtureClass-method
#'   logLik,MultipleGroupClass-method logLik,MixedClass-method logLik,DiscreteClass-method
#' @docType methods
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @rdname logLik-method
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 1)
#' logLik(x)
#'
#' }
setMethod(
    f = "logLik",
    signature = signature(object = 'SingleGroupClass'),
    definition = function(object){
        extract.mirt(object, 'logLik')
    }
)
