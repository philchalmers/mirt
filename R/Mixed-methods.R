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

# setMethod(
#     f = "summary",
#     signature = 'MixedClass',
#     definition = function(object, digits = 3, ...)
#     {
#         cat("\nCall:\n", paste(deparse(object@Call), sep = "\n", collapse = "\n"),
#             "\n\n", sep = "")
#         cat('Fixed effects:\n')
#         if(object@mixedlist$fixed.constrain){
#             fixed <- data.frame(Estimate=object@mixedlist$betas, 'Std.Error'=object@mixedlist$SEbetas,
#                                 row.names=names(object@mixedlist$SEbetas))
#         } else {
#             nfixed <- length(object@mixedlist$betas)
#             nms <- names(object@mixedlist$betas)
#             itemnames <- colnames(data)
#             J <- ncol(object@data)
#             betas <- SEbetas <- names <- c()
#             for(i in 1:J){
#                 tmp1 <- object@pars[[i]]@par[1:nfixed]
#                 tmp2 <- object@pars[[i]]@SEpar[1:nfixed]
#                 names <- c(names, paste(itemnames[i], '.', nms, sep = ''))
#                 betas <- c(betas, tmp1)
#                 SEbetas <- c(SEbetas, tmp2)
#             }
#             fixed <- data.frame(Estimate=betas, 'Std.Error'=SEbetas,
#                                 row.names=names)
#         }
#         fixed$'t.value' <- fixed$Estimate / fixed$'Std.Error'
#         print(round(fixed, digits))
#         cat('\n')
#     }
# )

setMethod(
    f = "coef",
    signature = 'MixedClass',
    definition = function(object, digits = 3, ...)
    {
        K <- object@K
        J <- length(K)
        nLambdas <- ncol(object@F)
        allPars <- list()
        if(length(object@pars[[1]]@SEpar) > 0){
            for(i in 1:(J+1)){
                allPars[[i]] <- round(matrix(c(object@pars[[i]]@par, object@pars[[i]]@SEpar),
                                             2, byrow = TRUE), digits)
                rownames(allPars[[i]]) <- c('pars', 'SE')
                colnames(allPars[[i]]) <- names(object@pars[[i]]@est)
            }
        } else {
            for(i in 1:(J+1)){
                allPars[[i]] <- round(object@pars[[i]]@par, digits)
                names(allPars[[i]]) <- names(object@pars[[i]]@est)
            }
        }
        names(allPars) <- c(colnames(object@data), 'GroupPars')
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

