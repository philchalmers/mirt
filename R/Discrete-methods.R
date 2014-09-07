# Methods
setMethod(
    f = "print",
    signature = signature(x = 'DiscreteClass'),
    definition = function(x)
    {
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
        cat("Latent class model with ", x@nfact, " classes.\n", sep="")
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
        cat('M-step optimizer used:', x@Moptim, '\n')
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
                cat(", RMSEA = ", round(x@RMSEA,3), sep = '')
            }
        }
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'DiscreteClass'),
    definition = function(object) {
        print(object)
    }
)

setMethod(
    f = "summary",
    signature = 'DiscreteClass',
    definition = function(object, printSE=FALSE, ...)
    {
        ngroups <- length(object@pars)
        Theta <- object@Theta
        ret <- vector('list', ngroups)
        items <- vector('list', object@Data$nitems + 1L)
        names(items) <- c(colnames(object@Data$data), 'Class.Proportions')
        for(g in 1L:ngroups){
            ret[[g]] <- items
            pars <- object@pars[[g]]
            for(i in 1L:object@Data$nitems){
                item <- extract.item(pars, i)
                P <- probtrace(item, Theta)
                colnames(P) <- paste0('category_', 1L:ncol(P))
                rownames(P) <- paste0('Class_', 1L:nrow(P))
                ret[[g]][[i]] <- P
            }
            ret[[g]][[i+1L]] <- object@Prior[[g]]
        }
        if(length(ret) == 1L) ret <- ret[[1L]]
        ret
    }
)
setMethod(
    f = "coef",
    signature = 'DiscreteClass',
    definition = function(object, printSE=FALSE, drop = TRUE, ...){
        class(object) <- 'MultipleGroupClass'
        ret <- coef(object,  ...)
        for(g in 1L:length(ret))
            ret[[g]][[length(ret[[g]])]] <- NULL
        if(drop)
            if(length(ret) == 1L) ret <- ret[[1L]]
        ret
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'DiscreteClass'),
    definition = function(object, object2, ...)
    {
        class(object) <- 'ExploratoryClass'
        anova(object, object2, ...)
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'DiscreteClass'),
    definition = function(object, ...)
    {
        class(object) <- 'MultipleGroupClass'
        ret <- residuals(object, discrete = TRUE, ...)
        if(length(ret) == 1L) ret <- ret[[1L]]
        ret
    }
)