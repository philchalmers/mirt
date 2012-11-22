#' Person fit statistics
#' 
#' \code{personfit} calculates the Zh values from Drasgow, Levine and Williams (1985) for 
#' unidimensional and multidimensional models. For Rasch models infit and outfit statistics are 
#' also produced.
#' The returned object is a \code{data.frame}
#' consisting either of the tabulated data or full data with the statistics appended to the last 
#' columns. 
#' 
#' 
#' @aliases personfit
#' @param x a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or 
#' \code{MultipleGroupClass}
#' @param method type of factor score estimation method. Can be expected
#' a-posteriori (\code{"EAP"}), Bayes modal (\code{"MAP"}), weighted likelihood estimation 
#' (\code{"WLE"}), or maximum likelihood (\code{"ML"}) 
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
#' a <- matrix(rlnorm(20),ncol=1)
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
personfit <- function(x, method = 'EAP'){        
    if(is(x, 'MultipleGroupClass')){
        ret <- list()   
        for(g in 1:length(x@cmods))
            ret[[g]] <- personfit(x@cmods[[g]])
        names(ret) <- names(x@cmods)
        return(ret)
    }
    full.scores <- TRUE
    sc <- fscores(x, verbose = FALSE, full.scores=full.scores, method=method) 
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
    if(all(x@itemtype %in% 'Rasch')){
        V <- resid2 <- info <- matrix(0, ncol=J, nrow=N)  
        K <- x@K
        for (i in 1:J){            
            P <- ProbTrace(x=pars[[i]], Theta=Theta)
            Emat <- matrix(0:(K[i]-1), nrow(P), ncol(P), byrow = TRUE)            
            dat <- fulldata[ ,itemloc[i]:(itemloc[i+1] - 1)]               
            item <- extract.item(x, i)
            info[,i] <- iteminfo(item, Theta)
            resid2[ ,i] <- rowSums((dat*Emat - Emat * P))^2 
        }
        if(!is.null(attr(x, 'inoutfitreturn'))) return(list(resid2=resid2, info=info))
        outfit <- rowSums(resid2/info) / J
        infit <- rowSums(resid2) / rowSums(info)    
        if(full.scores) ret <- data.frame(x@data, outfit=outfit, infit=infit, Zh=Zh)
        else ret <- data.frame(x@tabdata, outfit=outfit, infit=infit, Zh=Zh)
    } else {
        if(full.scores) ret <- data.frame(x@data, Zh=Zh)
        else ret <- data.frame(x@tabdata, Zh=Zh)        
    }
    return(ret)
}
    
