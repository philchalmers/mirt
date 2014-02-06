#' Print the model objects
#'
#' \code{print(x)}
#'
#' Print model object summaries to the console.
#'
#' @param x an object of class \code{ExploratoryClass}, \code{ConfirmatoryClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#'
#' @name print-method
#' @aliases print,ExploratoryClass-method print,ConfirmatoryClass-method
#'   print,MultipleGroupClass-method print,MixedClass-method
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
        if(!is.nan(x@condnum))
            cat("Condition number of information matrix = ", x@condnum,
                '\nSecond-order test: model ', if(!x@secondordertest)
                    'is not a maximum, or the information matrix is too inaccurate' else
                        'is a possible local maximum', '\n', sep = "")
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

#' Show model object
#'
#' \code{show(object)}
#'
#' Print model object summaries to the console.
#'
#' @param x an object of class \code{ExploratoryClass}, \code{ConfirmatoryClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#'
#' @name show-method
#' @aliases show,ExploratoryClass-method show,ConfirmatoryClass-method
#'   show,MultipleGroupClass-method show,MixedClass-method
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
    signature = signature(object = 'ExploratoryClass'),
    definition = function(object){
        print(object)
    }
)

#' Summary of model object
#'
#' \code{summary(object, rotate = '', Target = NULL, suppress = 0, digits = 3, verbose = TRUE, ...)}
#'
#' Transforms coefficients into a standardized factor loading's metric. For \code{MixedClass} objects,
#' the fixed and random coefficients are printed.
#'
#' @param object an object of class \code{ExploratoryClass}, \code{ConfirmatoryClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#' @param rotate see \code{\link{mirt}} for details
#' @param Target a dummy variable matrix indicting a target rotation pattern
#' @param suppress a numeric value indicating which (possibly rotated) factor
#'   loadings should be suppressed. Typical values are around .3 in most
#'   statistical software. Default is 0 for no suppression
#' @param digits number of significant digits to be rounded
#' @param verbose logical; allow information to be printed to the console?
#' @param ... additional arguments to be passed
#'
#' @name summary-method
#' @aliases summary,ExploratoryClass-method summary,ConfirmatoryClass-method
#'   summary,MultipleGroupClass-method summary,MixedClass-method
#' @docType methods
#' @rdname summary-method
#' @seealso \code{\link{coef-method}}
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 2)
#' summary(x)
#' summary(x, rotate = 'varimax')
#' }
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
                warning("Solution has Heywood cases. Interpret with caution.")
            invisible(list(rotF=rotF$loadings,h2=h2,fcor=Phi))
        }
    }
)

#' Extract raw coefs from model object
#'
#' \code{coef(object, CI = .95, printSE = FALSE, rotate = '', Target = NULL, digits = 3,
#'    IRTpars = FALSE, rawug = FALSE, as.data.frame = FALSE, verbose = TRUE, ...)}
#'
#' @param object an object of class \code{ExploratoryClass}, \code{ConfirmatoryClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#' @param CI the amount of converged used to compute confidence intervals; default is
#'   95 percent confidence intervals
#' @param IRTpars logical; convert slope intercept parameters into traditional IRT parameters?
#'   Only applicable to unidimensional models
#' @param rotate see \code{\link{mirt}} for details
#' @param Target a dummy variable matrix indicting a target rotation pattern
#' @param suppress a numeric value indicating which (possibly rotated) factor
#'   loadings should be suppressed. Typical values are around .3 in most
#'   statistical software. Default is 0 for no suppression
#' @param printSE logical; print the standard errors instead of the confidence intervals?
#' @param digits number of significant digits to be rounded
#' @param as.data.frame logical; convert list output to a data.frame instead?
#' @param verbose logical; allow information to be printed to the console?
#' @param rawug logical; return the untransformed internal g and u parameters?
#'   If \code{FALSE}, g and u's are converted with the original format along with delta standard errors
#' @param ... additional arguments to be passed
#'
#' @name coef-method
#' @aliases coef,ExploratoryClass-method coef,ConfirmatoryClass-method
#'   coef,MultipleGroupClass-method coef,MixedClass-method
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
    signature = 'ExploratoryClass',
    definition = function(object, CI = .95, printSE = FALSE, rotate = '', Target = NULL, digits = 3,
                          IRTpars = FALSE, rawug = FALSE, as.data.frame = FALSE, verbose = TRUE, ...){
        if(printSE) rawug <- TRUE
        if(CI >= 1 || CI <= 0)
            stop('CI must be between 0 and 1')
        if(rotate == ''){
            rotate <- try(slot(object, 'rotate'), TRUE)
            if(is(rotate, 'try-error')) rotate <- 'none'
        }
        z <- abs(qnorm((1 - CI)/2))
        SEnames <- paste0('CI_', c((1 - CI)/2*100, ((1 - CI)/2 + CI)*100))
        K <- object@K
        J <- length(K)
        nfact <- ncol(object@F)
        a <- matrix(0, J, nfact)
        for(i in 1:J)
            a[i, ] <- ExtractLambdas(object@pars[[i]])

        if (ncol(a) > 1 && rotate != 'none'){
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
        if(IRTpars){
            if(object@nfact > 1L)
                stop('traditional parameterization is only available for unidimensional models')
            for(i in 1:(J+1))
                allPars[[i]] <- round(mirt2traditional(object@pars[[i]]), digits)
        } else {
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
        }
        if(!rawug){
            allPars <- lapply(allPars, function(x, digits){
                x[ , colnames(x) %in% c('g', 'u')] <- round(antilogit(x[ , colnames(x) %in% c('g', 'u')]), digits)
                x
            },  digits=digits)
        }
        names(allPars) <- c(colnames(object@data), 'GroupPars')
        if(as.data.frame)
            allPars <- t(as.data.frame(allPars))
        return(allPars)
    }
)

#' Compare nested models
#'
#' \code{anova(object, object2, verbose = TRUE)}
#'
#' @param object an object of class \code{ExploratoryClass}, \code{ConfirmatoryClass},
#'   \code{MultipleGroupClass}, or \code{MixedClass}
#' @param object2 a second model estimated from any of the mirt package estimation methods
#' @param verbose logical; print additional information to console?
#'
#' @name anova-method
#' @aliases anova,ExploratoryClass-method anova,ConfirmatoryClass-method
#'   anova,MultipleGroupClass-method anova,MixedClass-method
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
    signature = signature(object = 'ExploratoryClass'),
    definition = function(object, object2, verbose = TRUE){
        df <- object@df - object2@df
        if(df < 0){
            temp <- object
            object <- object2
            object2 <- temp
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
                          X2 = c('', X2),
                          df = c('', abs(df)),
                          p = c('', round(1 - pchisq(X2,abs(df)),3)))
        return(ret)
    }
)

#' Compute model residuals
#'
#' \code{residuals(object, type = 'LD', digits = 3, df.p = FALSE, full.scores = FALSE,
#'                          printvalue = NULL, tables = FALSE, verbose = TRUE, Theta = NULL, ...)}
#'
#' @param object an object of class \code{ExploratoryClass}, \code{ConfirmatoryClass} or
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
#' @param ... additional arguments to be passed to \code{fscores()}
#'
#' @name residuals-method
#' @aliases residuals,ExploratoryClass-method residuals,ConfirmatoryClass-method
#'   residuals,MultipleGroupClass-method
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
#'
#' # with and without supplied factor scores
#' Theta <- fscores(x, full.scores=TRUE, scores.only=TRUE)
#' residuals(x, type = 'Q3', Theta=Theta)
#' residuals(x, type = 'Q3', method = 'ML')
#'
#' }
setMethod(
    f = "residuals",
    signature = signature(object = 'ExploratoryClass'),
    definition = function(object, type = 'LD', digits = 3, df.p = FALSE, full.scores = FALSE,
                          printvalue = NULL, tables = FALSE, verbose = TRUE, Theta = NULL, ...)
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
        quadpts <- object@quadpts
        if(is.nan(quadpts)) 
            quadpts <- switch(as.character(nfact), '1'=41, '2'=21, '3'=11, '4'=7, '5'=5, 3)
        bfactorlist <- object@bfactor
        theta <- as.matrix(seq(-(.8 * sqrt(quadpts)), .8 * sqrt(quadpts), length.out = quadpts))
        if(type != 'Q3'){
            if(is.null(bfactorlist$Priorbetween[[1L]])){
                Theta <- thetaComb(theta, nfact)
            } else {
                Theta <- object@Theta
            }
        } else if(is.null(Theta)){
            Theta <- fscores(object, verbose=FALSE, full.scores=TRUE, scores.only=TRUE, ...)
        }
        itemnames <- colnames(data)
        listtabs <- list()
        calcG2 <- ifelse(type == 'LDG2', TRUE, FALSE)
        if(type %in% c('LD', 'LDG2')){
            groupPars <- ExtractGroupPars(object@pars[[length(object@pars)]])
            prior <- mvtnorm::dmvnorm(Theta,groupPars$gmeans, groupPars$gcov)
            prior <- prior/sum(prior)
            df <- (object@K - 1) %o% (object@K - 1)
            diag(df) <- NA
            colnames(df) <- rownames(df) <- colnames(res)
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
            res <- round(res,digits)
            return(res)
        } else if(type == 'exp'){
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
        } else if(type == 'Q3'){
            dat <- matrix(NA, N, 2L)
            diag(res) <- 1
            for(i in 1L:J){
                ei <- extract.item(object, item=i)
                EI <- expected.item(ei, Theta=Theta)
                dat[ ,1L] <- object@data[ ,i] - EI
                for(j in 1L:J){
                    if(i < j){
                        ej <- extract.item(object, item=i)
                        EJ <- expected.item(ej, Theta=Theta)
                        dat[,2L] <- object@data[ ,j] - EJ
                        tmpdat <- na.omit(dat)
                        n <- nrow(tmpdat)
                        Sz <- sqrt(1 / (n-3))
                        res[i,j] <- res[j,i] <- cor(tmpdat)[1L,2L]
                    }
                }
            }
            if(verbose) cat("Q3 matrix:\n\n")
            res <- round(res,digits)
            return(res)
        } else {
            stop('specified type does not exist')
        }

    }
)

#' Plot various test implied functions from models
#'
#' \code{plot(x, y, type = 'info', npts = 50, theta_angle = 45,
#'                          which.items = 1:ncol(x@@data),
#'                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
#'                          facet_items = TRUE, auto.key = TRUE, ehist.cut = 1e-10, ...)}
#'
#' @param x an object of class \code{ExploratoryClass}, \code{ConfirmatoryClass} or
#'   \code{MultipleGroupClass}
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
#' @param npts number of quadrature points to be used for plotting features.
#'   Larger values make plots look smoother
#' @param rot allows rotation of the 3D graphics
#' @param which.items numeric vector indicating which items to be used when plotting. Default is
#'   to use all available items
#' @param facet_items logical; apply grid of plots across items? If \code{FALSE}, items will be
#'   placed in one plot for each group
#' @param auto.key logical parameter passed to the \code{lattice} package
#' @param ehist.cut a probability value indicating a threshold for excliding cases in empirical 
#'   histogram plots. Values larger than the default will include more points in the tails of the 
#'   plot, potentially squishing the 'meat' of the plot to take up less area than visually desired
#' @param ... additional arguments to be passed to lattice
#'
#' @name plot-method
#' @aliases plot,ExploratoryClass-method plot,ConfirmatoryClass-method
#'   plot,MultipleGroupClass-method
#' @docType methods
#' @rdname plot-method
#' @examples
#'
#' \dontrun{
#' x <- mirt(Science, 1)
#' plot(x)
#' plot(x, type = 'trace')
#' plot(x, type = 'infotrace')
#' plot(x, type = 'infotrace', facet_items = FALSE)
#' plot(x, type = 'infoSE')
#'
#' set.seed(1234)
#' group <- sample(c('g1','g2'), nrow(Science), TRUE)
#' x2 <- multipleGroup(Science, 1, group)
#' plot(x2)
#' plot(x2, type = 'trace')
#' plot(x2, type = 'trace', which.items = 1:2)
#' plot(x2, type = 'trace', which.items = 1, facet_items = FALSE) #facet by group
#' plot(x2, type = 'score')
#'
#' x3 <- mirt(Science, 2)
#' plot(x3)
#' plot(x3, type = 'SE')
#'
#' }
setMethod(
    f = "plot",
    signature = signature(x = 'ExploratoryClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45,
                          which.items = 1:ncol(x@data),
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                          facet_items = TRUE, auto.key = TRUE, main = NULL,
                          drape = TRUE, colorkey = TRUE, ehist.cut = 1e-10, add.ylab2 = TRUE, ...)
    {
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
        itemtrace <- computeItemtrace(x@pars, ThetaFull, x@itemloc, CUSTOM.IND=x@CUSTOM.IND)
        score <- c()
        for(i in 1:J)
            score <- c(score, 0:(x@K[i]-1))
        score <- matrix(score, nrow(itemtrace), ncol(itemtrace), byrow = TRUE)
        plt <- data.frame(cbind(info,score=rowSums(score*itemtrace),Theta=Theta))
        if(nfact == 2){
            colnames(plt) <- c("info", "score", "Theta1", "Theta2")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'infocontour'){
                if(is.null(main))
                    main <- paste("Test Information Contour")
                return(contourplot(info ~ Theta1 * Theta2, data = plt,
                                   main = main, xlab = expression(theta[1]),
                                   ylab = expression(theta[2])))
            } else if(type == 'scorecontour'){
                if(is.null(main))
                    main <- paste("Expected Score Contour")
                    return(contourplot(score ~ Theta1 * Theta2, data = plt,
                                       main = main, xlab = expression(theta[1]),
                                       ylab = expression(theta[2])))
            } else if(type == 'info'){
                if(is.null(main))
                    main <- "Test Information"
                return(wireframe(info ~ Theta1 + Theta2, data = plt, main = main,
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape))
            } else if(type == 'score'){
                if(is.null(main))
                    main <- "Expected Total Score"
                return(wireframe(score ~ Theta1 + Theta2, data = plt, main = main,
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape))
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
                                 scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape))
            } else {
                stop('plot type not supported for two dimensional model')
            }
        } else {
            colnames(plt) <- c("info", "score", "Theta")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'info'){
                if(is.null(main))
                    main <- 'Test Information'
                return(xyplot(info~Theta, plt, type='l', main = main,
                              xlab = expression(theta), ylab=expression(I(theta))))
            } else if(type == 'SE'){
                if(is.null(main))
                    main <- 'Test Standard Errors'
                return(xyplot(SE~Theta, plt, type='l', main = main,
                       xlab = expression(theta), ylab=expression(SE(theta))))
            } else if(type == 'infoSE'){
                if(is.null(main))
                    main <- 'Test Information and Standard Errors'
                obj1 <- xyplot(info~Theta, plt, type='l', main = main,
                               xlab = expression(theta), ylab=expression(I(theta)))
                obj2 <- xyplot(SE~Theta, plt, type='l', ylab=expression(SE(theta)))
                if(!require(latticeExtra)) require(latticeExtra)
                return(doubleYScale(obj1, obj2, add.ylab2 = add.ylab2))
            } else if(type == 'trace'){
                if(is.null(main))
                    main <- 'Item trace lines'
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
                items <- gl(n=length(unique(which.items)), k=nrow(Theta),
                            labels = paste('Item', which.items))
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
                return(xyplot(score ~ Theta, plt,
                              xlab = expression(theta), ylab = expression(Total(theta)),
                              type = 'l', main = main, ...))
            } else if(type == 'empiricalhist'){
                if(is.null(main))
                    main <- 'Empirical Histogram'
                if(all(is.nan(x@Prior))) stop('Empirical histogram was not estimated for this object')
                Theta <- as.matrix(seq(-(.8 * sqrt(x@quadpts)), .8 * sqrt(x@quadpts),
                                    length.out = x@quadpts))
                Prior <- x@Prior * nrow(x@data)
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
                stop('plot not supported for unidimensional models')
            }
        }
    }
)
