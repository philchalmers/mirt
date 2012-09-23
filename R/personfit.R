#' Person fit statistics
#' 
#' \code{personfit} calculates the Zh values from Drasgow, Levine and Williams (1985) for 
#' unidimensional and multidimensional models. The returned values approximate a standard normal
#' distribution and therefore p-values are also returned. The returned object is a \code{data.frame}
#' consisting either of the tabulated data or full data with the Zh and p-values appended to the last 
#' columns. 
#' 
#' 
#' @aliases personfit
#' @param x a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or 
#' \code{MultipleGroupClass}
#' @param full.scores if \code{FALSE} (default) then a summary table with the Zh and p-values is returned,
#' otherwise the original data matrix is returned with values appended to the rightmost column 
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
#' x <- mirt(data, 1, SE = FALSE)
#' tabdatafit <- personfit(x)
#' head(tabdatafit)
#' 
#'   }
#' 
personfit <- function(x, full.scores = FALSE){
    if(is(x, 'MultipleGroupClass')){
        ret <- list()   
        for(g in 1:length(x@cmods))
            ret[[g]] <- personfit(x@cmods[[g]], full.scores=full.scores)
        names(ret) <- names(x@cmods)
        return(ret)
    }    
    sc <- fscores(x, verbose = FALSE, full.scores=full.scores) 
    J <- ncol(x@data)
    itemloc <- x@itemloc
    nfact <- x@nfact        
    pars <- x@pars
    if(full.scores){
        fulldata <- x@fulldata    
        Theta <- sc[ ,ncol(sc):(ncol(sc) - nfact + 1), drop = FALSE]
    } else {
        fulldata <- x@tabdatalong
        fulldata <- fulldata[,-ncol(fulldata)]
        Theta <- sc[ ,ncol(sc):(ncol(sc) - nfact + 1) - nfact, drop = FALSE]        
    }
    itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=nrow(Theta))        
    for (i in 1:J){
        if(x@pars[[1]]@bfactor) pars[[i]]@bfactor <- FALSE
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)
    }
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
    p <- round(pnorm(Zh), 3) 
    p[Zh > 0] <- round(pnorm(Zh[Zh>0], lower.tail = FALSE), 3)
    p <- 2*p
    if(full.scores) ret <- data.frame(x@data, Zh=Zh, p=p)
    else ret <- data.frame(x@tabdata, Zh=Zh, p=p)
    ret
}
    