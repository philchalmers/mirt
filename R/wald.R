#' Wald test for mirt models
#' 
#' Compute a Wald test given an \code{L} vector or matrix of contrasts.
#' 
#' @aliases wald
#' @param L a coefficient matrix with dimensions nconstrasts x npars. Use \code{constrain = 'index'}
#' on the initially estimated model to obtain the parameter indicators  
#' @param object estimated object from mirt, confmirt, or multipleGroup
#' @param C a constant vector to be compared along side L. Default is 0
#' @keywords wald
#' @export wald
#' @examples
#' \dontrun{
#' #View parnumber index
#' mirt(Science, 2, constrain = 'index')
#' mod <- mirt(Science, 2)
#' 
#' #all first factor slopes equal to 0
#' L <- rep(0, 25)
#' L[c(1,6,11,16)] <- 1
#' wald(L, mod)
#' 
#' #first two items same factor 1 slope, last two items same factor 1 slope
#' L <- matrix(0, 2, 25)
#' L[1,1] <- L[2, 11] <- 1
#' L[1,6] <- L[2, 16] <- -1
#' wald(L, mod)
#' }
wald <- function(L, object, C = 0){
    pars <- object@pars
    covB <- solve(object@information)
    estB <- B <- c()
    if(is(object, 'MultipleGroupClass')){
        for(g in 1:length(pars)){
            for(i in 1:length(pars[[g]])){
                B <- c(B, pars[[g]][[i]]@par)        
                estB <- c(estB, pars[[g]][[i]]@est)
            }
        }        
    } else {
        for(i in 1:length(pars)){
            B <- c(B, pars[[i]]@par)        
            estB <- c(estB, pars[[i]]@est)
        }
    }
    estB <- matrix(estB, 1)
    B <- B[estB[1,]]       
    if(!is.matrix(L))
        L <- matrix(L, 1)    
    if(ncol(L) == ncol(estB)) L <- L[, estB, drop = FALSE]    
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
    cat('Wald test: \nW = ', round(W, 3), ', df = ', df, ', p = ', 
        round(p, 3), sep='')       
}

