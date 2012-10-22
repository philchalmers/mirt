#' Person fit statistics
#' 
#' \code{personfit} calculates the Zh values from Drasgow, Levine and Williams (1985) for 
#' unidimensional and multidimensional models, as well as infit and outfit statistics. 
#' The returned object is a \code{data.frame}
#' consisting either of the tabulated data or full data with the statistics appended to the last 
#' columns. 
#' 
#' 
#' @aliases personfit
#' @param x a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or 
#' \code{MultipleGroupClass}
#' @param degrees the degrees angle to be passed to the \code{\link{iteminfo}} function
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords person fit
#' @export personfit
#' 
#' @seealso
#' \code{\link{itemfit}}
#' 
#' @references
#' 
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with 
#' polychotomous item response models and standardized indices. 
#' \emph{Journal of Mathematical and Statistical Psychology, 38}, 67-86.
#' 
#' Reise, S. P. (1990). A comparison of item- and person-fit methods of assessing model-data fit 
#' in IRT. \emph{Applied Psychological Measurement, 14}, 127-137.
#' 
#' @examples
#' 
#' \dontrun{
#'
#' #make some data
#' set.seed(1234)
#' a <- matrix(rnorm(20),ncol=1)
#' d <- matrix(rnorm(20),ncol=1)
#' items <- rep('dich', 20)
#' data <- simdata(a,d, 2000, items)
#'  
#' x <- mirt(data, 1)
#' fit <- personfit(x)
#' head(fit)
#' 
#'   }
#' 
personfit <- function(x, degrees = NULL){        
    if(is(x, 'MultipleGroupClass')){
        ret <- list()   
        for(g in 1:length(x@cmods))
            ret[[g]] <- personfit(x@cmods[[g]], degrees=degrees)
        names(ret) <- names(x@cmods)
        return(ret)
    }
    full.scores <- TRUE
    sc <- fscores(x, verbose = FALSE, full.scores=full.scores) 
    J <- ncol(x@data)
    itemloc <- x@itemloc
    pars <- x@pars
    prodlist <- attr(pars, 'prodlist')    
    nfact <- x@nfact + length(prodlist)   
    if(full.scores){
        fulldata <- x@fulldata    
        Theta <- sc[ ,ncol(sc):(ncol(sc) - nfact + 1), drop = FALSE]
    } else {
        fulldata <- x@tabdatalong
        fulldata <- fulldata[,-ncol(fulldata)]
        Theta <- sc[ ,ncol(sc):(ncol(sc) - nfact + 1) - nfact, drop = FALSE]        
    }
    N <- nrow(Theta)
    itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)        
    for (i in 1:J)
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)
    LL <- itemtrace * fulldata
    LL[LL < 1e-8] <- 1
    LL <- rowSums(log(LL)) 
    Zh <- rep(0, length(LL))
    for(n in 1:nrow(Theta)){    
        mu <- sigma2 <- 0
        for(item in 1:J){
            P <- itemtrace[n ,itemloc[item]:(itemloc[item+1]-1)]
            mu <- mu + sum(P * log(P))            
            for(i in 1:length(P)){
                for(j in 1:length(P)){
                    if(i != j)
                        sigma2 <- sigma2 + P[i] * P[j] * log(P[i]) * log(P[i]/P[j])
                }
            }            
        }
        Zh[n] <- (LL[n] - mu) / sqrt(sigma2)
    }               
    V <- Z <- info <- matrix(0, ncol=J, nrow=N)        
    if(is.null(degrees))
        if(nfact > 1) degrees <- rep(90/nfact, nfact)            
    for (i in 1:J){
        P <- ProbTrace(x=pars[[i]], Theta=Theta)
        dat <- fulldata[ ,itemloc[i]:(itemloc[i+1] - 1)]            
        item <- extract.item(x, i)
        info[,i] <- iteminfo(item, Theta, degrees=degrees)
        Z[ ,i] <- rowSums(dat - dat * P) / sqrt(info[,i])                         
    }
    if(!is.null(attr(x, 'inoutfitreturn'))) return(list(Z=Z, info=info))
    outfit <- rowSums(Z^2) / J
    infit <- rowSums(Z^2 * info) / rowSums(info)    
    if(full.scores) ret <- data.frame(x@data, outfit=outfit, infit=infit, Zh=Zh)
    else ret <- data.frame(x@tabdata, outfit=outfit, infit=infit, Zh=Zh)        
    return(ret)
}
    
