#' Calculate bootstrapped standard errors for estimated models
#'
#' Given an internal mirt object estimate the bootstrapped standard errors. It may
#' be beneficial to run the computations using multi-core architecture (e.g., the \code{parallel}
#' package).
#'
#' @aliases boot.mirt
#' @param x an estimated model object
#' @param R number of draws to use (passed to the \code{boot()} function)
#' @param ... additional arguments to be passed on to \code{boot(...)}
#' @keywords bootstrapped standard errors
#' @export boot.mirt
#' @seealso
#' \code{\link{PLCI.mirt}}
#' @examples
#'
#' \dontrun{
#'
#' #standard
#' mod <- mirt(Science, 1)
#' booted <- boot.mirt(mod, R=20)
#' plot(booted)
#' booted
#'
#' #run in parallel using snow back-end using all available cores
#' mod <- mirt(Science, 1)
#' booted <- boot.mirt(mod, parallel = 'snow', ncpus = parallel::detectCores())
#' booted
#'
#'
#' }
boot.mirt <- function(x, R = 100, ...){
    boot.draws <- function(orgdat, ind, npars, constrain, parprior, model, itemtype, group,
                           class, LR, obj, ...) {
        ngroup <- length(unique(group))
        dat <- orgdat[ind, ]
        rownames(dat) <- NULL
        g <- group[ind]
        if(length(unique(g)) != ngroup) return(rep(NA, npars))
        if(class == 'MixedClass'){
            fm <- obj@formulas
            itemdesign <- if(length(obj@itemdesign)) obj@itemdesign else NULL
            covdata <- if(length(obj@covdata)) obj@covdata[ind, , drop=FALSE] else NULL
            mod <- try(mixedmirt(data=dat, model=model, covdata=covdata, itemtype=itemtype,
                                 itemdesign=itemdesign, fixed=fm$fixed, random=fm$random,
                                 lr.fixed=fm$lr.fixed, lr.random=fm$lr.random,
                                 constrain=constrain, parprior=parprior,
                                 verbose=FALSE, technical=list(parallel=FALSE),
                                 SE=FALSE, ...))
        } else if(class == 'DiscreteClass'){
            mod <- try(mdirt(data=dat, model=model, itemtype=itemtype, group=g,
                             constrain=constrain, parprior=parprior,
                             verbose=FALSE, technical=list(parallel=FALSE),
                             SE=FALSE, ...))
        } else {
            if(class == 'MultipleGroupClass'){
                mod <- try(multipleGroup(data=dat, model=model, itemtype=itemtype, group=g,
                                     constrain=constrain, parprior=parprior,
                                     calcNull=FALSE, verbose=FALSE, technical=list(parallel=FALSE),
                                     method=x@method, draws=1, SE=FALSE, ...))
            } else {
                if(.hasSlot(LR, 'beta')){
                    formula <- LR@formula
                    if(length(formula) == 1L) formula <- formula[[1]]
                    df <- LR@df[ind, ]
                } else {
                    formula = ~ 1
                    df <- NULL
                }
                mod <- try(mirt(data=dat, model=model, itemtype=itemtype, constrain=constrain,
                            parprior=parprior, calcNull=FALSE, verbose=FALSE,
                            technical=list(parallel=FALSE), formula=formula,
                            covdata=df, method=x@method, draws=1, SE=FALSE, ...))
            }
        }
        if(is(mod, 'try-error')) return(rep(NA, npars))
        structure <- mod2values(mod)
        longpars <- structure$value[structure$est]
        if(length(longpars) != npars) return(rep(NA, npars)) #in case intercepts dropped
        return(longpars)
    }

    return.boot <- TRUE
    dat <- x@Data$data
    method <- x@method
    itemtype <- x@itemtype
    class <- class(x)
    group <- if(class == 'MultipleGroupClass') x@Data$group else NULL
    model <- x@model[[1L]]
    parprior <- x@parprior
    constrain <- x@constrain
    LR <- x@lrPars
    if(length(parprior) == 0L) parprior <- NULL
    if(length(constrain) == 0L) constrain <- NULL
    prodlist <- x@prodlist
    ret <- x
    if(!require(boot)) require('boot')
    structure <- mod2values(x)
    longpars <- structure$value
    npars <- sum(structure$est)
    boots <- boot::boot(dat, boot.draws, R=R, npars=npars, constrain=constrain, class=class,
                  parprior=parprior, model=model, itemtype=itemtype, group=group, LR=LR,
                  obj=x, ...)
    names(boots$t0) <- paste(paste(structure$item[structure$est],
                             structure$name[structure$est], sep='.'),
                             structure$parnum[structure$est], sep='_')
    class(boots) <- c('boot.mirt', 'boot')
    if(class == 'ExploratoryClass')
        message('Note: bootstrapped standard errors for slope parameters for exploratory
                       models are not meaningful.')
    return(boots)
}

#' @rdname boot.mirt
#' @param digits number of decimal points to display
#' @param index which bootstrap draws to use. Default uses all draws
#' @method print boot.mirt
print.boot.mirt <- function(x, digits = getOption("digits"),
                           index = 1L:ncol(boot.out$t), ...)
{
    # adapted from the boot package version 1.3-13, December 10 2014
    if(requireNamespace("boot", quietly = TRUE)){

        #
        # Print the output of a bootstrap
        #
        boot.out <- x
        sim <- boot.out$sim
        cl <- boot.out$call
        t <- matrix(boot.out$t[, index], nrow = nrow(boot.out$t))
        allNA <- apply(t,2L,function(t) all(is.na(t)))
        ind1 <- index[allNA]
        index <- index[!allNA]
        t <- matrix(t[, !allNA], nrow = nrow(t))
        # the following two lines are the only edits
        #     rn <- paste("t",index,"*",sep="")
        rn <- names(boot.out$t0)
        if (length(index) == 0L)
            op <- NULL
        else if (is.null(t0 <- boot.out$t0)) {
            if (is.null(boot.out$call$weights))
                op <- cbind(apply(t,2L,mean,na.rm=TRUE),
                            sqrt(apply(t,2L,function(t.st) var(t.st[!is.na(t.st)]))))
            else {
                op <- NULL
                for (i in index)
                    op <- rbind(op, boot::imp.moments(boot.out,index=i)$rat)
                op[,2L] <- sqrt(op[,2])
            }
            dimnames(op) <- list(rn,c("mean", "std. error"))
        }
        else {
            t0 <- boot.out$t0[index]
            if (is.null(boot.out$call$weights)) {
                op <- cbind(t0,apply(t,2L,mean,na.rm=TRUE)-t0,
                            sqrt(apply(t,2L,function(t.st) var(t.st[!is.na(t.st)]))))
                dimnames(op) <- list(rn, c("original"," bias  "," std. error"))
            }
            else {
                op <- NULL
                for (i in index)
                    op <- rbind(op, boot::imp.moments(boot.out,index=i)$rat)
                op <- cbind(t0,op[,1L]-t0,sqrt(op[,2L]),
                            apply(t,2L,mean,na.rm=TRUE))
                dimnames(op) <- list(rn,c("original", " bias  ",
                                          " std. error", " mean(t*)"))
            }
        }
        if (cl[[1L]] == "boot") {
            if (sim == "parametric")
                cat("\nPARAMETRIC BOOTSTRAP\n\n")
            else if (sim == "antithetic") {
                if (is.null(cl$strata))
                    cat("\nANTITHETIC BOOTSTRAP\n\n")
                else    cat("\nSTRATIFIED ANTITHETIC BOOTSTRAP\n\n")
            }
            else if (sim == "permutation") {
                if (is.null(cl$strata))
                    cat("\nDATA PERMUTATION\n\n")
                else    cat("\nSTRATIFIED DATA PERMUTATION\n\n")
            }
            else if (sim == "balanced") {
                if (is.null(cl$strata) && is.null(cl$weights))
                    cat("\nBALANCED BOOTSTRAP\n\n")
                else if (is.null(cl$strata))
                    cat("\nBALANCED WEIGHTED BOOTSTRAP\n\n")
                else if (is.null(cl$weights))
                    cat("\nSTRATIFIED BALANCED BOOTSTRAP\n\n")
                else	cat("\nSTRATIFIED WEIGHTED BALANCED BOOTSTRAP\n\n")
            }
            else {
                if (is.null(cl$strata) && is.null(cl$weights))
                    cat("\nORDINARY NONPARAMETRIC BOOTSTRAP\n\n")
                else if (is.null(cl$strata))
                    cat("\nWEIGHTED BOOTSTRAP\n\n")
                else if (is.null(cl$weights))
                    cat("\nSTRATIFIED BOOTSTRAP\n\n")
                else 	cat("\nSTRATIFIED WEIGHTED BOOTSTRAP\n\n")
            }
        }
        else if (cl[[1L]] == "tilt.boot") {
            R <- boot.out$R
            th <- boot.out$theta
            if (sim == "balanced")
                cat("\nBALANCED TILTED BOOTSTRAP\n\n")
            else	cat("\nTILTED BOOTSTRAP\n\n")
            if ((R[1L] == 0) || is.null(cl$tilt) || eval(cl$tilt))
                cat("Exponential tilting used\n")
            else	cat("Frequency Smoothing used\n")
            i1 <- 1
            if (boot.out$R[1L]>0)
                cat(paste("First",R[1L],"replicates untilted,\n"))
            else {
                cat(paste("First ",R[2L]," replicates tilted to ",
                          signif(th[1L],4),",\n",sep=""))
                i1 <- 2
            }
            if (i1 <= length(th)) {
                for (j in i1:length(th))
                    cat(paste("Next ",R[j+1L]," replicates tilted to ",
                              signif(th[j],4L),
                              ifelse(j!=length(th),",\n",".\n"),sep=""))
            }
            op <- op[, 1L:3L]
        }
        else if (cl[[1L]] == "tsboot") {
            if (!is.null(cl$indices))
                cat("\nTIME SERIES BOOTSTRAP USING SUPPLIED INDICES\n\n")
            else if (sim == "model")
                cat("\nMODEL BASED BOOTSTRAP FOR TIME SERIES\n\n")
            else if (sim == "scramble") {
                cat("\nPHASE SCRAMBLED BOOTSTRAP FOR TIME SERIES\n\n")
                if (boot.out$norm)
                    cat("Normal margins used.\n")
                else	cat("Observed margins used.\n")
            }
            else if (sim == "geom") {
                if (is.null(cl$ran.gen))
                    cat("\nSTATIONARY BOOTSTRAP FOR TIME SERIES\n\n")
                else	cat(paste("\nPOST-BLACKENED STATIONARY",
                               "BOOTSTRAP FOR TIME SERIES\n\n"))
                cat(paste("Average Block Length of",boot.out$l,"\n"))
            }
            else {	if (is.null(cl$ran.gen))
                cat("\nBLOCK BOOTSTRAP FOR TIME SERIES\n\n")
                else	cat(paste("\nPOST-BLACKENED BLOCK",
                               "BOOTSTRAP FOR TIME SERIES\n\n"))
                cat(paste("Fixed Block Length of",boot.out$l,"\n"))
            }
        }
        else {
            cat("\n")
            if (sim == "weird") {
                if (!is.null(cl$strata)) cat("STRATIFIED ")
                cat("WEIRD BOOTSTRAP FOR CENSORED DATA\n\n")
            }
            else if ((sim == "ordinary") ||
                         ((sim == "model") && is.null(boot.out$cox))) {
                if (!is.null(cl$strata)) cat("STRATIFIED ")
                cat("CASE RESAMPLING BOOTSTRAP FOR CENSORED DATA\n\n")
            }
            else if (sim == "model") {
                if (!is.null(cl$strata)) cat("STRATIFIED ")
                cat("MODEL BASED BOOTSTRAP FOR COX REGRESSION MODEL\n\n")
            }
            else if (sim == "cond") {
                if (!is.null(cl$strata)) cat("STRATIFIED ")
                cat("CONDITIONAL BOOTSTRAP ")
                if (is.null(boot.out$cox))
                    cat("FOR CENSORED DATA\n\n")
                else	cat("FOR COX REGRESSION MODEL\n\n")
            }
        }
        cat("\nCall:\n")
        dput(cl, control=NULL)
        cat("\n\nBootstrap Statistics :\n")
        if (!is.null(op)) print(op,digits=digits)
        if (length(ind1) > 0L)
            for (j in ind1)
                cat(paste("WARNING: All values of t", j, "* are NA\n", sep=""))
    }
    invisible(boot.out)
}
