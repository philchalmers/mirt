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
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
        cat("Full-information item factor analysis with ", x@nfact, " factor(s).\n", sep="")
        EMquad <- ''
        if(x@method == 'EM') EMquad <- c('\n     using ', x@quadpts, ' quadrature')
        method <- x@method
        if(method == 'MIXED') method <- 'MHRM'
        if(x@converge == 1)
            cat("Converged within ", x@TOL, ' tolerance after ', x@iter, ' ',
                method, " iterations.\n", sep = "")
        else
            cat("FAILED TO CONVERGE within ", x@TOL, ' tolerance after ',
                x@iter, ' ', method, " iterations.\n", sep="")
        cat('mirt version:', as.character(packageVersion('mirt')), '\n')
        cat('M-step optimizer:', x@Moptim, '\n')
        if(method == 'EM' || method == 'BL'){
            cat('EM acceleration:', x@accelerate)
            cat('\nNumber of rectangular quadrature:', x@quadpts)
            cat('\n')
        }
        if(!is.nan(x@condnum)){
            cat("\nInformation matrix estimated with method:", x@infomethod)
            cat("\nCondition number of information matrix = ", x@condnum,
                '\nSecond-order test: model ', if(!x@secondordertest)
                    'is not a maximum, or the information matrix is too inaccurate' else
                        'is a possible local maximum', '\n', sep = "")
        }
        if(length(x@logLik) > 0){
            cat("\nLog-likelihood = ", x@logLik, if(method == 'MHRM')
                paste(', SE =', round(x@SElogLik,3)), "\n",sep='')
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
#' @docType methods
#' @rdname show-method
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
#' the fixed and random coefficients are printed.
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
#' @param printCI print a confidence interval for standardized loadings
#'   (e.g., \code{printCI = .95} gives a 95\% confidence interval)
#' @param digits number of significant digits to be rounded
#' @param verbose logical; allow information to be printed to the console?
#' @param ... additional arguments to be passed
#'
#' @name summary-method
#' @aliases summary,SingleGroupClass-method
#'   summary,MultipleGroupClass-method summary,MixedClass-method summary,DiscreteClass-method
#' @docType methods
#' @rdname summary-method
#' @seealso \code{\link{coef-method}}
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 2)
#' summary(x)
#' summary(x, rotate = 'varimax')
#'
#' #print confidence interval (requires computed information matrix)
#' x2 <- mirt(Science, 1, SE=TRUE)
#' summary(x2, printCI=.95)
#' }
setMethod(
    f = "summary",
    signature = 'SingleGroupClass',
    definition = function(object, rotate = 'oblimin', Target = NULL, suppress = 0, digits = 3,
                          printCI = FALSE, verbose = TRUE, ...){
        nfact <- ncol(object@F)
        if (!object@exploratory || rotate == 'none') {
            F <- object@F
            F[abs(F) < suppress] <- NA
            h2 <- as.matrix(object@h2)
            SS <- apply(F^2,2,sum)
            gp <- ExtractGroupPars(object@pars[[length(object@pars)]])
            Phi <- cov2cor(gp$gcov)
            colnames(h2) <- "h2"
            rownames(Phi) <- colnames(Phi) <- names(SS) <- colnames(F)[1L:object@nfact]
            loads <- round(cbind(F,h2),digits)
            rownames(loads) <- colnames(object@Data$data)
            if(verbose){
                if(object@exploratory)
                    cat("\nUnrotated factor loadings: \n\n")
                print(loads)
                cat("\nSS loadings: ",round(SS,digits), "\n")
                cat("Proportion Var: ",round(SS/nrow(F),digits), "\n")
                cat("\nFactor correlations: \n\n")
                print(round(Phi, digits))
            }
            invisible(list(rotF=F,h2=h2,fcor=Phi))
        } else {
            F <- object@F
            h2 <- as.matrix(object@h2)
            colnames(h2) <- "h2"
            rotF <- Rotate(F, rotate, Target = Target, ...)
            SS <- apply(rotF$loadings^2,2,sum)
            L <- rotF$loadings
            L[abs(L) < suppress] <- NA
            loads <- round(cbind(L,h2),digits)
            rownames(loads) <- colnames(object@Data$data)
            Phi <- diag(ncol(F))
            if(!rotF$orthogonal){
                Phi <- rotF$Phi
            }
            colnames(Phi) <- rownames(Phi) <- colnames(F)
            if(verbose){
                cat("\nRotation: ", rotate, "\n")
                cat("\nRotated factor loadings: \n\n")
                print(loads,digits)
                cat("\nRotated SS loadings: ",round(SS,digits), "\n")
                cat("\nFactor correlations: \n\n")
                print(round(Phi, digits))
            }
            if(any(h2 > 1))
                warning("Solution has Heywood cases. Interpret with caution.", call.=FALSE)
            invisible(list(rotF=rotF$loadings,h2=h2,fcor=Phi))
        }
    }
)

#' Extract raw coefs from model object
#'
#' Return a list (or data.frame) of raw item and group level coefficients.
#'
#' @param object an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#' @param CI the amount of converged used to compute confidence intervals; default is
#'   95 percent confidence intervals
#' @param IRTpars logical; convert slope intercept parameters into traditional IRT parameters?
#'   Only applicable to unidimensional models
#' @param rotate see \code{\link{mirt}} for details. The default rotation is \code{'none'}
#' @param Target a dummy variable matrix indicting a target rotation pattern
#' @param printSE logical; print the standard errors instead of the confidence intervals?
#' @param digits number of significant digits to be rounded
#' @param as.data.frame logical; convert list output to a data.frame instead?
#' @param simplify logical; if all items have the same parameter names (indicating they are
#'   of the same class) then they are collapsed to a matrix, and a list of length 2 is returned
#'   containing a matrix of item parameters and group-level estimates
#' @param verbose logical; allow information to be printed to the console?
#' @param rawug logical; return the untransformed internal g and u parameters?
#'   If \code{FALSE}, g and u's are converted with the original format along with delta standard errors
#' @param ... additional arguments to be passed
#'
#' @name coef-method
#' @aliases coef,SingleGroupClass-method
#'   coef,MultipleGroupClass-method coef,MixedClass-method coef,DiscreteClass-method
#' @docType methods
#' @rdname coef-method
#' @seealso \code{\link{summary-method}}
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' x <- mirt(dat, 1)
#' coef(x)
#' coef(x, IRTpars = TRUE)
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
    definition = function(object, CI = .95, printSE = FALSE, rotate = 'none', Target = NULL, digits = 3,
                          IRTpars = FALSE, rawug = FALSE, as.data.frame = FALSE,
                          simplify=FALSE, verbose = TRUE, ...){
        if(printSE && length(object@pars[[1L]]@SEpar)) rawug <- TRUE
        if(CI >= 1 || CI <= 0)
            stop('CI must be between 0 and 1', call.=FALSE)
        z <- abs(qnorm((1 - CI)/2))
        SEnames <- paste0('CI_', c((1 - CI)/2*100, ((1 - CI)/2 + CI)*100))
        K <- object@K
        J <- length(K)
        nfact <- ncol(object@F)
        a <- matrix(0, J, nfact)
        for(i in 1:J)
            a[i, ] <- ExtractLambdas(object@pars[[i]])

        if (object@exploratory && rotate != 'none'){
            if(verbose) cat("\nRotation: ", rotate, "\n\n")
            so <- summary(object, rotate=rotate, Target=Target, verbose=FALSE, digits=digits, ...)
            a <- rotateLambdas(so) * 1.702
            for(i in 1:J)
                object@pars[[i]]@par[1:nfact] <- a[i, ]
            object@pars[[J + 1]]@par[-c(1:nfact)] <- so$fcor[lower.tri(so$fcor, TRUE)]
        }
        allPars <- list()
        if(IRTpars){
            if(object@nfact > 1L)
                stop('traditional parameterization is only available for unidimensional models',
                     call.=FALSE)
            for(i in 1:(J+1))
                allPars[[i]] <- round(mirt2traditional(object@pars[[i]]), digits)
        } else {
            if(length(object@pars[[1L]]@SEpar)){
                if(printSE){
                    for(i in 1L:(J+1L)){
                        allPars[[i]] <- round(matrix(c(object@pars[[i]]@par,
                                                       object@pars[[i]]@SEpar),
                                                     2, byrow = TRUE), digits)
                        rownames(allPars[[i]]) <- c('par', 'SE')
                        nms <- names(object@pars[[i]]@est)
                        if(i <= J && object@itemtype[i] != 'custom'){
                            nms[nms == 'g'] <- 'logit(g)'
                            nms[nms == 'u'] <- 'logit(u)'
                        }
                        colnames(allPars[[i]]) <- nms
                    }
                } else {
                    for(i in 1L:(J+1L)){
                        allPars[[i]] <- round(matrix(c(object@pars[[i]]@par,
                                                       object@pars[[i]]@par - z*object@pars[[i]]@SEpar,
                                                       object@pars[[i]]@par + z*object@pars[[i]]@SEpar),
                                                     3, byrow = TRUE), digits)
                        rownames(allPars[[i]]) <- c('par', SEnames)
                        colnames(allPars[[i]]) <- names(object@pars[[i]]@est)
                    }
                }
            } else {
                for(i in 1L:(J+1L)){
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
        names(allPars) <- c(colnames(object@Data$data), 'GroupPars')
        if(as.data.frame)
            allPars <- t(as.data.frame(allPars))
        if(simplify && !as.data.frame){
            allPars <- lapply(allPars, function(x) x[1L, , drop=FALSE])
            items.old <- allPars[1L:(length(allPars)-1L)]
            nms <- lapply(items.old, colnames)
            unms <- unique(do.call(c, nms))
            items <- matrix(NA, length(items.old), length(unms))
            rownames(items) <- names(items.old)
            colnames(items) <- unms
            for(i in 1L:nrow(items))
                items[i, nms[[i]]] <- items.old[[i]]
            means <- allPars$GroupPars[1L:object@nfact]
            covs <- matrix(NA, object@nfact, object@nfact)
            covs[lower.tri(covs, TRUE)] <- allPars$GroupPars[-c(1L:object@nfact)]
            colnames(covs) <- rownames(covs) <- names(means) <- object@factorNames[1L:object@nfact]
            allPars <- list(items=items, means=means, cov=covs)
        }
        if(.hasSlot(object@lrPars, 'beta'))
            allPars$lr.betas <- round(object@lrPars@beta, digits)
        return(allPars)
    }
)

#' Compare nested models with likelihood-based statistics
#'
#' Compare nested models using likelihood ratio, AIC, BIC, etc.
#'
#' @param object an object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#' @param object2 a second model estimated from any of the mirt package estimation methods
#' @param verbose logical; print additional information to console?
#'
#' @name anova-method
#' @aliases anova,SingleGroupClass-method
#'   anova,MultipleGroupClass-method anova,MixedClass-method anova,DiscreteClass-method
#' @docType methods
#' @rdname anova-method
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 1)
#' x2 <- mirt(Science, 2)
#' anova(x, x2)
#' }
setMethod(
    f = "anova",
    signature = signature(object = 'SingleGroupClass'),
    definition = function(object, object2, verbose = TRUE){
        df <- object@df - object2@df
        if(df < 0){
            temp <- object
            object <- object2
            object2 <- temp
        } else if(df == 0){
            if((2*object2@logLik - 2*object@logLik) < 0){
                temp <- object
                object <- object2
                object2 <- temp
            }
        }
        X2 <- round(2*object2@logLik - 2*object@logLik, 3)
        if(verbose){
            cat('\nModel 1: ')
            print(object@Call)
            cat('Model 2: ')
            print(object2@Call)
            cat('\n')
        }
        ret <- data.frame(AIC = c(object@AIC, object2@AIC),
                          AICc = c(object@AICc, object2@AICc),
                          SABIC = c(object@SABIC, object2@SABIC),
                          BIC = c(object@BIC, object2@BIC),
                          logLik = c(object@logLik, object2@logLik),
                          X2 = c(NaN, X2),
                          df = c(NaN, abs(df)),
                          p = c(NaN, round(1 - pchisq(X2,abs(df)),4)))
        return(ret)
    }
)

#' Compute model residuals
#'
#' Return model implied residuals for linear dependencies between items or at the person level.
#'
#' @param object an object of class \code{SingleGroupClass} or
#'   \code{MultipleGroupClass}. Bifactor models are automatically detected and utilized for
#'   better accuracy
#' @param type type of residuals to be displayed.
#'   Can be either \code{'LD'} or \code{'LDG2'} for a local dependence matrix based on the
#'   X2 or G2 statistics (Chen & Thissen, 1997), \code{'Q3'} for the statistic proposed by
#'   Yen (1984), or \code{'exp'} for the expected values for the frequencies of every response pattern
#' @param tables logical; for LD type, return the observed, expected, and standardized residual
#'   tables for each item combination?
#' @param digits number of significant digits to be rounded
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
#' @param quadpts number of quadrature nodes to use. The default is extracted from model (if available)
#'   or generated automatically if not available
#' @param suppress a numeric value indicating which parameter local dependency combinations
#'   to flag as being too high. Absolute values for the standardized estimates greater than
#'   this value will be returned, while all values less than this value will be set to NA
#' @param ... additional arguments to be passed to \code{fscores()}
#'
#' @name residuals-method
#' @aliases residuals,SingleGroupClass-method
#'   residuals,MultipleGroupClass-method residuals,DiscreteClass-method
#' @docType methods
#' @rdname residuals-method
#' @references
#'
#' Chen, W. H. & Thissen, D. (1997). Local dependence indices for item pairs using item
#' response theory. \emph{Journal of Educational and Behavioral Statistics, 22}, 265-289.
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
#'
#' # with and without supplied factor scores
#' Theta <- fscores(x, full.scores=TRUE, scores.only=TRUE)
#' residuals(x, type = 'Q3', Theta=Theta)
#' residuals(x, type = 'Q3', method = 'ML')
#'
#' }
setMethod(
    f = "residuals",
    signature = signature(object = 'SingleGroupClass'),
    definition = function(object, type = 'LD', digits = 3, df.p = FALSE, full.scores = FALSE,
                          printvalue = NULL, tables = FALSE, verbose = TRUE, Theta = NULL,
                          suppress = 1, theta_lim = c(-6, 6), quadpts = NULL, ...)
    {
        dots <- list(...)
        if(.hasSlot(object@lrPars, 'beta'))
            stop('Latent regression models not yet supported')
        discrete <- FALSE
        if(!is.null(dots$discrete)) discrete <- TRUE
        K <- object@K
        data <- object@Data$data
        N <- nrow(data)
        J <- ncol(data)
        nfact <- ncol(object@F)
        itemloc <- object@itemloc
        res <- matrix(0,J,J)
        diag(res) <- NA
        colnames(res) <- rownames(res) <- colnames(data)
        if(!discrete){
            if(is.null(quadpts))
                quadpts <- object@quadpts
            if(is.nan(quadpts))
                quadpts <- select_quadpts2(nfact)
            bfactorlist <- object@bfactor
            theta <- as.matrix(seq(theta_lim[1L], theta_lim[2L], length.out = quadpts))
            if(type != 'Q3'){
                if(is.null(bfactorlist$Priorbetween[[1L]])){
                    Theta <- thetaComb(theta, nfact)
                } else {
                    Theta <- object@Theta
                }
            } else if(is.null(Theta)){
                Theta <- fscores(object, verbose=FALSE, full.scores=TRUE, scores.only=TRUE, ...)
            }
        } else {
            Theta <- object@Theta
            if(!any(type %in% c('exp', 'LD', 'LDG2')))
                stop('residual type not supported for discrete latent variables', call.=FALSE)
        }
        itemnames <- colnames(data)
        listtabs <- list()
        calcG2 <- ifelse(type == 'LDG2', TRUE, FALSE)
        if(type %in% c('LD', 'LDG2')){
            if(!discrete){
                groupPars <- ExtractGroupPars(object@pars[[length(object@pars)]])
                prior <- mirt_dmvnorm(Theta,groupPars$gmeans, groupPars$gcov)
                prior <- prior/sum(prior)
            } else {
                prior <- object@Prior[[1L]]
            }
            df <- (object@K - 1) %o% (object@K - 1)
            diag(df) <- NA
            colnames(df) <- rownames(df) <- colnames(res)
            for(i in 1L:J){
                for(j in 1L:J){
                    if(i < j){
                        P1 <- ProbTrace(x=object@pars[[i]], Theta=Theta)
                        P2 <- ProbTrace(x=object@pars[[j]], Theta=Theta)
                        tab <- table(data[,i],data[,j])
                        Etab <- matrix(0,K[i],K[j])
                        for(k in 1L:K[i])
                            for(m in 1:K[j])
                                Etab[k,m] <- N * sum(P1[,k] * P2[,m] * prior)
                        s <- gamma.cor(tab) - gamma.cor(Etab)
                        if(s == 0) s <- 1
                        if(calcG2){
                            tmp <- tab
                            tmp[tab == 0] <- NA
                            res[j,i] <- 2 * sum(tmp * log(tmp/Etab), na.rm=TRUE) * sign(s)
                        } else {
                            res[j,i] <- sum(((tab - Etab)^2)/Etab) * sign(s)
                        }
                        res[i,j] <- sign(res[j,i]) * sqrt( abs(res[j,i]) / (N*min(c(K[i],K[j]) - 1L)))
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
                cat("Degrees of freedom (lower triangle) and p-values:\n\n")
                print(round(df, digits))
                cat("\n")
            }
            if(verbose) cat("LD matrix (lower triangle) and standardized values:\n\n")
            if(suppress < 1){
                pick <- abs(res[upper.tri(res)]) < suppress
                res[lower.tri(res)] <- res[upper.tri(res)][pick] <- NA
            }
            res <- round(res,digits)
            return(res)
        } else if(type == 'exp'){
            r <- object@Data$Freq[[1L]]
            res <- round((r - object@Pl * nrow(object@Data$data)) /
                             sqrt(object@Pl * nrow(object@Data$data)),digits)
            expected <- round(N * object@Pl,digits)
            tabdata <- object@Data$tabdata
            rownames(tabdata) <- NULL
            ISNA <- is.na(rowSums(tabdata))
            expected[ISNA] <- res[ISNA] <- NA
            tabdata <- data.frame(tabdata,object@Data$Freq[[1L]],expected,res)
            colnames(tabdata) <- c(colnames(object@Data$tabdata),"freq","exp","res")
            if(full.scores){
                tabdata[, 'exp'] <- object@Pl / r * N
                tabdata2 <- object@Data$tabdata
                stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
                sfulldata <- apply(object@Data$data, 1, paste, sep='', collapse = '/')
                scoremat <- tabdata[match(sfulldata, stabdata2), 'exp', drop = FALSE]
                res <- (1-scoremat) / sqrt(scoremat)
                colnames(res) <- 'res'
                ret <- cbind(object@Data$data, scoremat, res)
                ret[is.na(rowSums(ret)), c('exp', 'res')] <- NA
                rownames(ret) <- NULL
                return(ret)
            } else {
                tabdata <- tabdata[do.call(order, as.data.frame(tabdata[,1:J])),]
                if(!is.null(printvalue)){
                    if(!is.numeric(printvalue)) stop('printvalue is not a number.', call.=FALSE)
                    tabdata <- tabdata[abs(tabdata[ ,ncol(tabdata)]) > printvalue, ]
                }
                return(tabdata)
            }
        } else if(type == 'Q3'){
            dat <- matrix(NA, N, 2L)
            diag(res) <- 1
            for(i in 1L:J){
                ei <- extract.item(object, item=i)
                EI <- expected.item(ei, Theta=Theta)
                dat[ ,1L] <- object@Data$data[ ,i] - EI
                for(j in 1L:J){
                    if(i < j){
                        ej <- extract.item(object, item=i)
                        EJ <- expected.item(ej, Theta=Theta)
                        dat[,2L] <- object@Data$data[ ,j] - EJ
                        tmpdat <- na.omit(dat)
                        n <- nrow(tmpdat)
                        Sz <- sqrt(1 / (n-3))
                        res[i,j] <- res[j,i] <- cor(tmpdat)[1L,2L]
                    }
                }
            }
            if(verbose) cat("Q3 matrix:\n\n")
            if(suppress < 1){
                pick <- abs(res[upper.tri(res)]) < suppress
                res[lower.tri(res)] <- res[upper.tri(res)][pick] <- NA
            }
            res <- round(res,digits)
            return(res)
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
#' @param type type of plot to view; can be \code{'info'} to show the test
#'   information function, \code{'infocontour'} for the test information contours,
#'   \code{'SE'} for the test standard error function, \code{'trace'} and \code{'infotrace'}
#'   for all item probability information or trace lines (only available when all items are dichotomous),
#'   \code{'infoSE'} for a combined test information and standard error plot, and \code{'score'} and
#'   \code{'scorecontour'} for the expected total score surface and contour plots.
#'   If \code{empiricalhist = TRUE} was used in estimation then the type \code{'empiricalhist'}
#'   also will be available to generate the empirical histogram plot
#' @param theta_angle numeric values ranging from 0 to 90 used in \code{plot}.
#'   If a vector is used then a bubble plot is created with the summed information across the angles specified
#'   (e.g., \code{theta_angle = seq(0, 90, by=10)})
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
#' @param auto.key logical parameter passed to the \code{lattice} package
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
#' @aliases plot,SingleGroupClass-method
#'   plot,MultipleGroupClass-method plot,SingleGroupClass,missing-method
#'   plot,DiscreteClass,missing-method
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
#'
#' # confidence interval plots when information matrix computed
#' plot(x)
#' plot(x, MI=100)
#' plot(x, type='info', MI=100)
#' plot(x, type='SE', MI=100)
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
#' plot(x2, type = 'trace', which.items = 1, facet_items = FALSE) #facet by group
#' plot(x2, type = 'info')
#'
#' x3 <- mirt(Science, 2)
#' plot(x3, type = 'info')
#' plot(x3, type = 'SE')
#'
#' }
setMethod(
    f = "plot",
    signature = signature(x = 'SingleGroupClass', y = 'missing'),
    definition = function(x, y, type = 'score', npts = 50, theta_angle = 45,
                          theta_lim = c(-6,6), which.items = 1:ncol(x@Data$data),
                          MI = 0, CI = .95, rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                          facet_items = TRUE, auto.key = TRUE, main = NULL,
                          drape = TRUE, colorkey = TRUE, ehist.cut = 1e-10, add.ylab2 = TRUE, ...)
    {
        dots <- list(...)
        if (any(theta_angle > 90 | theta_angle < 0))
            stop('Improper angle specified. Must be between 0 and 90.', call.=FALSE)
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
        nfact <- x@nfact
        if(length(theta_angle) > nfact) type = 'infoangle'
        if(nfact > 3) stop("Can't plot high dimensional solutions.", call.=FALSE)
        if(nfact == 2 && length(theta_angle) == 1L)
            theta_angle <- c(theta_angle, 90 - theta_angle)
        if(nfact == 3 && length(theta_angle) == 1L) theta_angle <- rep(90/3, 3)
        if(nfact == 1) theta_angle <- 0
        J <- length(x@pars) - 1
        theta <- seq(theta_lim[1L],theta_lim[2L],length.out=npts)
        if(nfact == 3) theta <- seq(theta_lim[1L],theta_lim[2L], length.out=20)
        ThetaFull <- Theta <- thetaComb(theta, nfact)
        prodlist <- attr(x@pars, 'prodlist')
        if(length(prodlist) > 0)
            ThetaFull <- prodterms(Theta,prodlist)
        info <- 0
        if(any(sapply(x@pars, is , 'custom')) && type != 'trace' && type != 'score')
            stop('Information function for custom classes not available', call.=FALSE)
        if(any(sapply(x@pars, is , 'ideal')) && type != 'trace' && type != 'score') ##TODO
            warning('Information function for ideal point models are currently experimental',
                    call.=FALSE)
        if(all(!sapply(x@pars, is , 'custom'))){
            for(l in 1:length(theta_angle)){
                ta <- theta_angle[l]
                if(nfact == 2) ta <- c(theta_angle[l], 90 - theta_angle[l])
                if(nfact == 3) ta <- theta_angle
                for(i in 1:J)
                    info <- info + iteminfo(x=x@pars[[i]], Theta=ThetaFull, degrees=ta)
            }
        }
        adj <- x@Data$mins
        if (x@exploratory){
            if(!is.null(dots$rotate)){
                so <- summary(x, verbose=FALSE, digits=5, ...)
                a <- rotateLambdas(so) * 1.702
                for(i in 1:J)
                    x@pars[[i]]@par[1:nfact] <- a[i, ]
            }
        }
        itemtrace <- computeItemtrace(x@pars, ThetaFull, x@itemloc, CUSTOM.IND=x@CUSTOM.IND)
        score <- c()
        for(i in 1:J)
            score <- c(score, 0:(x@K[i]-1) + adj[i])
        score <- matrix(score, nrow(itemtrace), ncol(itemtrace), byrow = TRUE)
        plt <- data.frame(cbind(info,score=rowSums(score*itemtrace),Theta=Theta))
        if(MI > 0L && nfact == 1L){
            tmpx <- x
            if(is(try(chol(x@information), silent=TRUE), 'try-error'))
                stop('Proper information matrix must be precomputed in model', call.=FALSE)
            tmppars <- x@pars
            covB <- solve(x@information)
            names <- colnames(covB)
            tmp <- lapply(names, function(x, split){
                as.numeric(strsplit(x, split=split)[[1L]][-1L])
            }, split='\\.')
            imputenums <- do.call(c, tmp)
            CIscore <- CIinfo <- matrix(0, MI, length(plt$score))
            for(i in 1L:MI){
                while(TRUE){
                    tmp <- try(imputePars(pars=x@pars, covB=covB,
                                          imputenums=imputenums, constrain=x@constrain),
                               silent=TRUE)
                    if(!is(tmp, 'try-error')) break
                }
                tmpx@pars <- tmp
                itemtrace <- computeItemtrace(tmpx@pars, ThetaFull, x@itemloc, CUSTOM.IND=x@CUSTOM.IND)
                tmpscore <- rowSums(score * itemtrace)
                CIscore[i, ] <- tmpscore
                CIinfo[i, ] <- testinfo(tmpx, ThetaFull)[,1L]
            }
        }
        if(nfact == 3){
            colnames(plt) <- c("info", "score", "Theta1", "Theta2", "Theta3")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'infocontour'){
                if(is.null(main))
                    main <- paste("Test Information Contour")
                return(contourplot(info ~ Theta1 * Theta2 | Theta3, data = plt,
                                   main = main, xlab = expression(theta[1]),
                                   ylab = expression(theta[2]), ...))
            } else if(type == 'scorecontour'){
                if(is.null(main))
                    main <- paste("Expected Score Contour")
                return(contourplot(score ~ Theta1 * Theta2 | Theta3, data = plt,
                                   main = main, xlab = expression(theta[1]),
                                   ylab = expression(theta[2]), ...))
            } else if(type == 'info'){
                if(is.null(main))
                    main <- "Test Information"
                return(wireframe(info ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape, ...))
            } else if(type == 'SEcontour'){
                if(is.null(main))
                    main <- "Test Standard Errors"
                return(contourplot(score ~ Theta1 * Theta2 | Theta3, data = plt,
                                          main = main, xlab = expression(theta[1]),
                                          ylab = expression(theta[2]), ...))
            } else if(type == 'score'){
                if(is.null(main))
                    main <- "Expected Total Score"
                return(wireframe(score ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape, ...))
            } else if(type == 'SE'){
                if(is.null(main))
                    main <- "Test Standard Errors"
                return(wireframe(SE ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                                 zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape, ...))
            } else {
                stop('plot type not supported for three dimensional model', call.=FALSE)
            }
        } else if(nfact == 2){
            colnames(plt) <- c("info", "score", "Theta1", "Theta2")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'infocontour'){
                if(is.null(main))
                    main <- paste("Test Information Contour")
                return(contourplot(info ~ Theta1 * Theta2, data = plt,
                                   main = main, xlab = expression(theta[1]),
                                   ylab = expression(theta[2]), ...))
            } else if(type == 'scorecontour'){
                if(is.null(main))
                    main <- paste("Expected Score Contour")
                    return(contourplot(score ~ Theta1 * Theta2, data = plt,
                                       main = main, xlab = expression(theta[1]),
                                       ylab = expression(theta[2]), ...))
            } else if(type == 'info'){
                if(is.null(main))
                    main <- "Test Information"
                return(wireframe(info ~ Theta1 + Theta2, data = plt, main = main,
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape, ...))
            } else if(type == 'SEcontour'){
                if(is.null(main))
                    main <- "Test Standard Errors"
                return(contourplot(score ~ Theta1 * Theta2, data = plt,
                                          main = main, xlab = expression(theta[1]),
                                          ylab = expression(theta[2]), ...))
            } else if(type == 'score'){
                if(is.null(main))
                    main <- "Expected Total Score"
                return(wireframe(score ~ Theta1 + Theta2, data = plt, main = main,
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape, ...))
            } else if(type == 'infoangle'){
                if(is.null(main))
                    main <- 'Information across different angles'
                symbols(plt[,2], plt[,3], circles = sqrt(plt[,1]/pi), inches = .35, fg='white', bg='blue',
                        xlab = expression(theta[1]), ylab = expression(theta[2]),
                        main = main)
            } else if(type == 'SE'){
                if(is.null(main))
                    main <- "Test Standard Errors"
                return(wireframe(SE ~ Theta1 + Theta2, data = plt, main = main,
                                 zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape, ...))
            } else {
                stop('plot type not supported for two dimensional model', call.=FALSE)
            }
        } else {
            colnames(plt) <- c("info", "score", "Theta")
            plt$SE <- 1 / sqrt(plt$info)
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
            }
            if(type == 'info'){
                if(is.null(main))
                    main <- 'Test Information'
                if(MI > 0){
                    return(xyplot(info ~ Theta, data=plt,
                                  upper=plt$CIinfoupper, lower=plt$CIinfolower,
                                  panel = function(x, y, lower, upper, ...){
                                      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                    col=grey(.9), border = FALSE, ...)
                                      panel.xyplot(x, y, type='l', lty=1,...)
                                  },
                                  main = main, ylim=c(min(plt$CIinfolower), max(plt$CIinfoupper)),
                                  ylab = expression(I(theta)), xlab = expression(theta), ...))
                } else {
                    return(xyplot(info~Theta, plt, type='l', main = main,
                                  xlab = expression(theta), ylab=expression(I(theta)), ...))
                }
            } else if(type == 'SE'){
                if(is.null(main))
                    main <- 'Test Standard Errors'
                if(MI > 0){
                    return(xyplot(SE ~ Theta, data=plt,
                                  upper=plt$CISEupper, lower=plt$CISElower,
                                  panel = function(x, y, lower, upper, ...){
                                      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                    col=grey(.9), border = FALSE, ...)
                                      panel.xyplot(x, y, type='l', lty=1,...)
                                  },
                                  main = main, ylim=c(min(plt$CISElower), max(plt$CISEupper)),
                                  ylab = expression(I(theta)), xlab = expression(theta), ...))
                } else {
                    return(xyplot(SE~Theta, plt, type='l', main = main,
                           xlab = expression(theta), ylab=expression(SE(theta)), ...))
                }
            } else if(type == 'infoSE'){
                if(is.null(main))
                    main <- 'Test Information and Standard Errors'
                obj1 <- xyplot(info~Theta, plt, type='l', main = main,
                               xlab = expression(theta), ylab=expression(I(theta)))
                obj2 <- xyplot(SE~Theta, plt, type='l', ylab=expression(SE(theta)))
                if(!require(latticeExtra)) require(latticeExtra)
                return(doubleYScale(obj1, obj2, add.ylab2 = add.ylab2, ...))
            } else if(type == 'trace'){
                if(is.null(main))
                    main <- 'Item trace lines'
                P <- vector('list', length(which.items))
                names(P) <- colnames(x@Data$data)[which.items]
                ind <- 1L
                for(i in which.items){
                    tmp <- probtrace(extract.item(x, i), ThetaFull)
                    if(ncol(tmp) == 2L) tmp <- tmp[,2, drop=FALSE]
                    tmp2 <- data.frame(P=as.numeric(tmp), cat=gl(ncol(tmp), k=nrow(Theta),
                                                           labels=paste0('P', 1L:ncol(tmp))))
                    P[[ind]] <- tmp2
                    ind <- ind + 1L
                }
                nrs <- sapply(P, nrow)
                Pstack <- do.call(rbind, P)
                names <- c()
                for(i in 1L:length(nrs))
                    names <- c(names, rep(names(P)[i], nrs[i]))
                plotobj <- data.frame(Pstack, item=names, Theta=Theta)
                if(facet_items){
                    return(xyplot(P ~ Theta|item, plotobj, ylim = c(-0.1,1.1), group = cat,
                                  xlab = expression(theta), ylab = expression(P(theta)),
                                  auto.key = auto.key, type = 'l', main = main, ...))
                } else {
                    return(xyplot(P ~ Theta, plotobj, group=item, ylim = c(-0.1,1.1),
                                  xlab = expression(theta), ylab = expression(P(theta)),
                                  auto.key = auto.key, type = 'l', main = main, ...))
                }
            } else if(type == 'infotrace'){
                if(is.null(main))
                    main <- 'Item information trace lines'
                I <- matrix(NA, nrow(Theta), J)
                for(i in which.items)
                    I[,i] <- iteminfo(extract.item(x, i), ThetaFull)
                I <- t(na.omit(t(I)))
                items <- rep(colnames(x@Data$data)[which.items], each=nrow(Theta))
                plotobj <- data.frame(I = as.numeric(I), Theta=Theta, item=items)
                if(facet_items){
                    return(xyplot(I ~ Theta|item, plotobj,
                                  xlab = expression(theta), ylab = expression(I(theta)),
                                  auto.key = auto.key, type = 'l', main = main, ...))
                } else {
                    return(xyplot(I ~ Theta, plotobj, group = item,
                                  xlab = expression(theta), ylab = expression(I(theta)),
                                  auto.key = auto.key, type = 'l', main = main, ...))
                }
            } else if(type == 'score'){
                if(is.null(main))
                    main <- 'Expected Total Score'
                if(MI > 0){
                    return(xyplot(score ~ Theta, data=plt,
                                  upper=plt$CIscoreupper, lower=plt$CIscorelower,
                                  panel = function(x, y, lower, upper, ...){
                                      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                    col=grey(.9), border = FALSE, ...)
                                      panel.xyplot(x, y, type='l', lty=1,...)
                                  },
                                  main = main,
                                  ylab = expression(T(theta)), xlab = expression(theta), ...))
                } else {
                    return(xyplot(score ~ Theta, plt,
                                  xlab = expression(theta), ylab = expression(T(theta)),
                                  type = 'l', main = main, ...))
                }
            } else if(type == 'empiricalhist'){
                if(is.null(main))
                    main <- 'Empirical Histogram'
                Prior <- x@Prior[[1L]]
                if(!x@empiricalhist)
                    stop('Empirical histogram was not estimated for this object', call.=FALSE)
                Theta <- as.matrix(seq(-(.8 * sqrt(x@quadpts)), .8 * sqrt(x@quadpts),
                                    length.out = x@quadpts))
                Prior <- Prior * nrow(x@Data$data)
                cuts <- cut(Theta, floor(npts/2))
                Prior <- do.call(c, lapply(split(Prior, cuts), mean))
                Theta <- do.call(c, lapply(split(Theta, cuts), mean))
                keep1 <- min(which(Prior > ehist.cut))
                keep2 <- max(which(Prior > ehist.cut))
                plt <- data.frame(Theta = Theta, Prior = Prior)
                plt <- plt[keep1:keep2, , drop=FALSE]
                return(xyplot(Prior ~ Theta, plt,
                              xlab = expression(theta), ylab = 'Expected Frequency',
                              type = 'b', main = main, ...))
            } else {
                stop('plot not supported for unidimensional models', call.=FALSE)
            }
        }
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
        ds <- ds[-c(1L:ncat)]
        newd <- numeric(length(ds)-1L)
        for(i in 2:length(ds))
            newd[i-1L] <- -(ds[i] - ds[i-1L])
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

traditional2mirt <- function(x, cls, ncat, digits = 3){
    if(cls == 'dich'){
        par <- x
        par[2L] <- -par[2L]*par[1L]
        names(par) <- c('a1', 'd', 'g', 'u')
    } else if(cls == 'graded'){
        par <- x
        for(i in 2L:ncat)
            par[i] <- -par[i]*par[1L]
        names(par) <- c('a1', paste0('d', 1:(length(par)-1)))
    } else if(cls == 'gpcm'){
        par <- c(x[1L], 0L:(ncat-1L), 0, x[-1L])
        ds <- -par[-c(1:(ncat+1))]*par[1]
        newd <- numeric(length(ds))
        for(i in length(ds):2L)
            newd[i] <- (ds[i] + ds[i-1L])
        for(i in length(newd):3L)
            newd[i] <- newd[i] + newd[i-2L]
        par <- c(par[1:(ncat+1)], newd)
        names(par) <- c('a1', paste0('ak', 0:(ncat-1)), paste0('d', 0:(ncat-1)))
    } else if(cls == 'nominal'){
        as <- x[1L:(length(x)/2)]
        ds <- x[-c(1L:(length(x)/2))]
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
