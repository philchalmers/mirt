#' Wald statistics for mirt models
#'
#' Compute a Wald test given an \code{L} vector or matrix of numeric contrasts. Requires that the
#' model information matrix be computed (including \code{SE = TRUE} when using the EM method). Use
#' \code{wald(model)} to observe how the information matrix columns are named, especially if
#' the estimated model contains constrained parameters (e.g., 1PL).
#'
#'
#' @aliases wald
#' @param L a coefficient matrix with dimensions nconstrasts x npars.
#'   Omitting this value will return the column names of the
#'   information matrix used to identify the (potentially constrained) parameters
#' @param object estimated object from \code{mirt}, \code{bfactor},
#'   \code{multipleGroup}, \code{mixedmirt}, or \code{mdirt}
#' @param C a constant vector of population parameters to be compared along side L, where
#'   \code{length(C) == ncol(L)}. By default a vector of 0's is constructed
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords wald
#' @export wald
#' @examples
#' \dontrun{
#' #View parnumber index
#' data(LSAT7)
#' data <- expand.table(LSAT7)
#' mod <- mirt(data, 1, SE = TRUE)
#' coef(mod)
#'
#' # see how the information matrix relates to estimated parameters, and how it lines up
#' #   with the parameter index
#' (infonames <- wald(mod))
#' index <- mod2values(mod)
#' index[index$est, ]
#'
#' #second item slope equal to 0?
#' L <- matrix(0, 1, 10)
#' L[1,3] <- 1
#' wald(mod, L)
#'
#' #simultaneously test equal factor slopes for item 1 and 2, and 4 and 5
#' L <- matrix(0, 2, 10)
#' L[1,1] <- L[2, 7] <- 1
#' L[1,3] <- L[2, 9] <- -1
#' L
#' wald(mod, L)
#'
#' #logLiklihood tests (requires estimating a new model)
#' cmodel <- 'theta = 1-5
#'            CONSTRAIN = (1,2, a1), (4,5, a1)'
#' mod2 <- mirt(data, cmodel)
#' #or, eqivalently
#' #mod2 <- mirt(data, 1, constrain = list(c(1,5), c(13,17)))
#' anova(mod2, mod)
#'
#' }
wald <- function(object, L, C = 0){
    if(missing(object)) missingMsg('object')
    covB <- extract.mirt(object, 'vcov')
    if(all(dim(covB) == c(1,1)))
        if(covB[1,1] == 0L)
            stop('No information matrix has been calculated for the model', call.=FALSE)
    Names <- colnames(covB)
    B <- extract.mirt(object, 'parvec')
    if(missing(L)){
        index <- 1L:length(Names)
        ret <- as.data.frame(t(data.frame(infoname=Names, par = round(B, 3))))
        colnames(ret) <- index
        return(ret)
    }
    if(!is.matrix(L)){
        stop('L must be a matrix', call.=FALSE)
    } else if(ncol(L) != length(Names)){
        stop('L does not have an appropriate number of columns', call.=FALSE)
    }
    if(!is.vector(C))
        stop('C must be a vector of constant population parameters', call.=FALSE)
    if(length(C) == 1L){
        C <- numeric(ncol(L))
    } else if(length(C) != ncol(L)){
        stop('length(C) must be the same as ncol(L)', call.=FALSE)
    }
    W <- t(L %*% (B - C)) %*% solve(L %*% covB %*% t(L)) %*% (L %*% (B - C))
    W <- ifelse(W < 0, 0, W)
    ret <- list(W=W, df = nrow(L))
    p <- 1 - pchisq(ret$W, ret$df)
    ret$p <- p
    as.data.frame(ret)
}
