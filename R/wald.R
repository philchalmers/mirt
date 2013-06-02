#' Wald test for mirt models
#'
#' Compute a Wald test given an \code{L} vector or matrix of numeric contrasts. Requires that the
#' model information matrix be computed (including \code{SE = TRUE} when using the EM method). Use
#' \code{wald(model)} to observe how the information matrix columns are named, especially if
#' the estimated model contains constrained parameters (e.g., 1PL). The information matrix names
#' are labelled according to which parameter number(s) they correspond to (to check the original
#' numbering use the option \code{pars = 'values'} in the original estimation function).
#'
#'
#' @aliases wald
#' @param L a coefficient matrix with dimensions nconstrasts x npars, or a vector if only one
#' set of contrasts is being tested. Omitting this value will return the column names of the
#' information matrix used to identify the (potentially constrained) parameters
#' @param object estimated object from \code{mirt}, \code{bfactor}, \code{confmirt},
#' \code{multipleGroup}, or \code{mixedmirt}
#' @param C a constant vector/matrix to be compared along side L
#' @keywords wald
#' @export wald
#' @examples
#' \dontrun{
#' #View parnumber index
#' data(LSAT7)
#' data <- expand.table(LSAT7)
#' cmodel <- confmirt.model()
#'    F1 = 1,4,5
#'    F2 = 2,3
#'
#'
#' mod <- mirt(data, cmodel, SE = TRUE)
#' coef(mod)
#'
#' #see how the information matrix relates to estimated parameters, and how it lines up with the index
#' (infonames <- wald(mod))
#' index <- mirt(data, cmodel, pars = 'values')
#' index
#'
#' #second factor slope equal to 0?
#' L <- rep(0, 10)
#' names(L) <- infonames
#' L[3] <- 1
#' wald(mod, L)
#'
#' #simultaneously test equal factor slopes for item 2 and 3, and 4 and 5
#' L <- matrix(0, 2, 10)
#' colnames(L) <- infonames #colnames() not required
#' L[1,3] <- L[2, 7] <- 1
#' L[1,5] <- L[2, 9] <- -1
#' L
#' wald(mod, L)
#'
#' #logLiklihood tests (requires estimating a new model)
#' mod2 <- mirt(data, cmodel, constrain = list(c(7,12), c(16,21)))
#' anova(mod2, mod)
#' }
wald <- function(object, L, C = 0){
    Names <- colnames(object@information)
    if(missing(L)){
        names(Names) <- 1:length(Names)
        return(Names)
    }
    if(!is.matrix(L))
        L <- matrix(L, 1)
    pars <- object@pars
    covB <- solve(object@information)
    B <- parnum <- c()
    if(is(object, 'MultipleGroupClass')){
        for(g in 1L:length(pars)){
            for(i in 1L:length(pars[[g]]@pars)){
                B <- c(B, pars[[g]]@pars[[i]]@par)
                parnum <- c(parnum, pars[[g]]@pars[[i]]@parnum)
            }
        }
    } else {
        for(i in 1L:length(pars)){
            B <- c(B, pars[[i]]@par)
            parnum <- c(parnum, pars[[i]]@parnum)
        }
    }
    keep <- c()
    for(i in 1L:length(Names))
        keep <- c(keep, as.numeric(strsplit(Names[i], '.', fixed = TRUE)[[1]][2]))
    B <- B[keep]
    W <- t(L %*% B - C) %*% solve(L %*% covB %*% t(L)) %*% (L %*% B - C)
    ret <- list(W=W, df = nrow(L))
    class(ret) <- 'wald'
    ret
}

#' @S3method print wald
#' @rdname wald
#' @method print wald
#' @param x an object of class 'wald'
#' @param ... additional arguments to be passed
print.wald <- function(x, ...){
    W <- x$W
    df <- x$df
    p <- 1 - pchisq(x$W, x$df)
    cat('\nWald test: \nW = ', round(W, 3), ', df = ', df, ', p = ',
        round(p, 3), '\n', sep='')
}

