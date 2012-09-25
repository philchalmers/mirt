#' Wald test for mirt models
#' 
#' Compute a Wald test given an \code{L} vector or matrix of contrasts.
#' 
#' @aliases wald
#' @param L a coefficient matrix with dimensions nconstrasts x npars. Use \code{pars = 'values'}
#' on the initially estimated model to obtain the parameter indicators  
#' @param object estimated object from mirt, confmirt, or multipleGroup
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
#' mod <- mirt(data, cmodel)
#' coef(mod, allpars = TRUE)
#' index <- mirt(data, cmodel, pars = 'values')
#' index
#' 
#'        
#' #second factor slopes equal to 0?
#' L <- rep(0, 30)
#' L[c(7, 12)] <- 1
#' wald(L, mod)
#' 
#' #simultaniously test equal factor slopes for item 2 and 3, and 4 and 5
#' L <- matrix(0, 2, 30)
#' L[1,16] <- L[2, 7] <- 1
#' L[1,21] <- L[2, 12] <- -1
#' wald(L, mod)
#' 
#' #logLiklihood tests (requires estimating a new model)
#' mod2 <- mirt(data, cmodel, constrain = list(c(7,12), c(16,21)))
#' anova(mod2, mod)
#' }
wald <- function(L, object, C = 0){
    pars <- object@pars
    covB <- solve(object@information)
    estB <- B <- c()
    if(is(object, 'MultipleGroupClass')){
        pars <- object@cmods
        for(g in 1:length(pars)){
            for(i in 1:length(pars[[g]]@pars)){
                B <- c(B, pars[[g]]@pars[[i]]@par)        
                estB <- c(estB, pars[[g]]@pars[[i]]@est)
            }
        }        
    } else {
        for(i in 1:length(pars)){
            B <- c(B, pars[[i]]@par)        
            estB <- c(estB, pars[[i]]@est)
        }
    }
    if(length(object@constrain) > 0){        
        constr <- object@constrain
        for(i in 1:length(constr))
            for(j in 2:length(constr[[i]]))
                estB[constr[[i]][j]] <- FALSE        
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

