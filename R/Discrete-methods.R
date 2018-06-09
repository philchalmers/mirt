# Methods
setMethod(
    f = "print",
    signature = signature(x = 'DiscreteClass'),
    definition = function(x)
    {
        cat("\nCall:\n", paste(deparse(x@Call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
        cat("Latent class model with ", x@Model$nfact, " classes and ", nrow(x@Model$Theta),
            " profiles.\n", sep="")
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
        cat('EM acceleration:', x@Options$accelerate)
        cat('\nLatent density type:', 'discrete')
        # cat('\nNumber of classes quadrature:', x@Options$quadpts)
        cat('\n')
        if(!is.na(x@OptimInfo$condnum)){
            cat("\nInformation matrix estimated with method:", x@Options$SE.type)
            cat("\nCondition number of information matrix = ", x@OptimInfo$condnum,
                '\nSecond-order test: model ', if(!x@OptimInfo$secondordertest)
                    'is not a maximum, or the information matrix is too inaccurate' else
                        'is a possible local maximum', '\n', sep = "")
        }
        if(length(x@Fit$logLik) > 0){
            cat("\nLog-likelihood = ", x@Fit$logLik, if(method == 'MHRM')
                paste(', SE =', round(x@Fit$SElogLik,3)), "\n",sep='')
            cat('Estimated parameters:', extract.mirt(x, 'nestpars'), '\n')
            cat("AIC = ", x@Fit$AIC, "; AICc = ", x@Fit$AICc, "\n", sep='')
            cat("BIC = ", x@Fit$BIC, "; SABIC = ", x@Fit$SABIC, "\n", sep='')
            if(!is.nan(x@Fit$p)){
                cat("G2 (", x@Fit$df,") = ", round(x@Fit$G2,2), ", p = ", round(x@Fit$p,4), sep='')
                cat(", RMSEA = ", round(x@Fit$RMSEA,3), sep = '')
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
        ngroups <- object@Data$ngroups
        Theta <- object@Model$Theta
        colnames(Theta) <- extract.mirt(object, 'factorNames')[1:ncol(Theta)]
        ret <- vector('list', ngroups)
        items <- vector('list', object@Data$nitems + 1L)
        names(items) <- c(colnames(object@Data$data), 'Class.Probability')
        prof <- apply(Theta, 1L, function(x) sprintf("P[%s]", paste0(x, collapse=' ')))
        for(g in seq_len(ngroups)){
            ret[[g]] <- items
            pars <- object@ParObjects$pars[[g]]
            for(i in seq_len(object@Data$nitems)){
                item <- extract.item(pars, i)
                P <- probtrace(item, Theta)
                colnames(P) <- paste0('category_', 1L:ncol(P))
                rownames(P) <- prof
                ret[[g]][[i]] <- P
            }
            if(is.matrix(object@Internals$Prior[[g]])){
                ret[[g]][[i+1L]] <- data.frame(Theta, prob=colMeans(object@Internals$Prior[[g]]))
            } else {
                ret[[g]][[i+1L]] <- data.frame(Theta, prob=object@Internals$Prior[[g]])
            }
            rownames(ret[[g]][[i+1L]]) <- paste0('Profile_', 1L:nrow(ret[[g]][[i+1L]]))
        }
        names(ret) <- extract.mirt(object, 'groupNames')
        for(i in seq_len(length(ret))) class(ret[[i]]) <- c('mirt_list', 'list')
        if(length(ret) == 1L) ret <- ret[[1L]]
        ret
    }
)
setMethod(
    f = "coef",
    signature = 'DiscreteClass',
    definition = function(object, drop = TRUE, ...){
        class(object) <- 'MultipleGroupClass'
        ret <- coef(object, discrete = TRUE, ...)
        for(g in seq_len(length(ret)))
            if(!is.null(ret[[g]]$lr.betas))
                ret[[g]]$GroupPars <- NULL
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
        class(object) <- 'SingleGroupClass'
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

setMethod(
    f = "plot",
    signature = signature(x = 'DiscreteClass', y = 'missing'),
    definition = function(x, which.items = 1:extract.mirt(x, 'nitems'),
                          facet_items = TRUE, type = 'b', profile = FALSE,
                          par.strip.text = list(cex = 0.7),
                          par.settings = list(strip.background = list(col = '#9ECAE1'),
                                              strip.border = list(col = "black")),
                          auto.key = list(space = 'right', points=FALSE, lines=TRUE), ...)
    {
        if(extract.mirt(x, 'ngroups') > 1L)
            stop('plot methods do not support multiple group latent class models yet',
                 call.=FALSE)
        so <- summary(x)
        index <- which.items
        names <- colnames(x@Data$data)
        mlt <- lapply(index, function(x, so, names){
            pick <- so[[x]]
            colnames(pick) <- paste0('cat', 1:ncol(pick))
            ret <- data.frame(item=names[x],
                              class=rep(rownames(pick), ncol(pick)),
                              cat=rep(colnames(pick), each=nrow(pick)),
                              prob=as.numeric(pick))
            ret
        }, so=so, names=names)
        mlt <- do.call(rbind, mlt)
        mlt$item <- factor(mlt$item, levels = colnames(x@Data$data)[which.items])
        if(profile){
            if(all(x@Data$K == 2L)){
                mlt <- mlt[mlt$cat == 'cat2', ]
                return(xyplot(prob ~ item, data=mlt, groups = class, type = type,
                              auto.key = auto.key, ylab = 'Probability', ylim = c(-.1, 1.1),
                              par.settings=par.settings, par.strip.text=par.strip.text, ...))
            } else {
                return(xyplot(prob ~ item|cat, data=mlt, groups = class, type = type,
                              auto.key = auto.key, ylab = 'Probability', ylim = c(-.1, 1.1),
                              par.settings=par.settings, par.strip.text=par.strip.text, ...))
            }
        }
        if(facet_items){
            return(xyplot(prob ~ cat|item, data=mlt, groups = class, type = type,
                          auto.key = auto.key, ylab = 'Probability', ylim = c(-.1, 1.1),
                          par.settings=par.settings, par.strip.text=par.strip.text, ...))
        } else {
            return(xyplot(prob ~ cat|class, data=mlt, groups = mlt$item, type = type,
                          auto.key = auto.key, ylab = 'Probability', ylim = c(-.1, 1.1),
                          par.settings=par.settings, par.strip.text=par.strip.text, ...))
        }
    }
)

# Methods
setMethod(
    f = "vcov",
    signature = signature(object = 'DiscreteClass'),
    definition = function(object)
    {
        class(object) <- 'SingleGroupClass'
        vcov(object)
    }
)

setMethod(
    f = "logLik",
    signature = signature(object = 'DiscreteClass'),
    definition = function(object){
        extract.mirt(object, 'logLik')
    }
)
