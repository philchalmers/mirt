#' Item fit statistics
#' 
#' \code{itemfit} calculates the Zh values from Drasgow, Levine and Williams (1985) for 
#' unidimensional and multidimensional models, or \eqn{\chi^2} values for unidimensional models.
#' 
#' 
#' @aliases itemfit
#' @param x a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or 
#' \code{MultipleGroupClass}
#' @param type a character specifying whether the Zh (\code{'Zh'}) or \eqn{\chi^2} (\code{'X2'}) statistic
#' should be computed. Not that \code{'X2'} can only be used for unidimensional models
#' @param ngroups the number of theta groupings to use when computing \code{'X2'}. Cells that have 
#' any expected values less than 5 are dropped and the degrees of freedom are adjusted accordingly
#' @param empirical.plot a single numeric value indicating which item to plot (via \code{itemplot}) and
#' overlay with the empirical \eqn{\theta} groupings. Only applicable when \code{type = 'X2'}. 
#' The default is \code{NULL}, therefore no plots are drawn 
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords item fit
#' @export itemfit
#' 
#' @seealso
#' \code{\link{personfit}}
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
#' fit <- itemfit(x)
#' fit
#' 
#' itemfit(x, type = 'X2', empirical.plot = 1) #empirical item plot
#' 
#'   }
#'
itemfit <- function(x, type = 'Zh', ngroups = 10, empirical.plot = NULL){
    if(is(x, 'MultipleGroupClass')){
        ret <- list()   
        for(g in 1:length(x@cmods))
            ret[[g]] <- itemfit(x@cmods[[g]], type=type, ngroups=ngroups)
        names(ret) <- names(x@cmods)
        return(ret)
    }    
    sc <- fscores(x, verbose = FALSE, full.scores = TRUE) 
    J <- ncol(x@data)
    itemloc <- x@itemloc
    nfact <- x@nfact
    pars <- x@pars
    fulldata <- x@fulldata    
    Theta <- sc[ ,ncol(sc):(ncol(sc) - nfact + 1), drop = FALSE]
    if(type == 'X2'){
        if(nfact > 1) stop('Item chi-squared values are only available for unidimensional models')
        ord <- order(Theta[,1])    
        fulldata <- fulldata[ord,]
        Theta <- Theta[ord, , drop = FALSE]
        den <- dnorm(Theta, 0, .5)
        den <- den / sum(den)
        cumTheta <- cumsum(den)
        Groups <- rep(20, length(ord))
        weight <- 1/ngroups        
        for(i in 1:20)
            Groups[round(cumTheta,2) >= weight*(i-1) & round(cumTheta,2) < weight*i] <- i        
        n.uniqueGroups <- length(unique(Groups))
        X2 <- df <- RMSEA <- rep(0, J)            
        if(!is.null(empirical.plot)) itemplot(x, empirical.plot)
        for (i in 1:J){    
            if(!is.null(empirical.plot) && i != empirical.plot) next            
            for(j in 1:n.uniqueGroups){                                    
                dat <- fulldata[Groups == j, itemloc[i]:(itemloc[i+1] - 1), drop = FALSE]            
                r <- colSums(dat)
                N <- nrow(dat)                  
                mtheta <- matrix(mean(Theta[Groups == j,]), nrow=1)
                if(!is.null(empirical.plot)){
                    tmp <- r/N
                    col <- 2:length(tmp)
                    points(rep(mtheta, length(tmp) - 1), tmp[-1], col = col)                  
                }
                P <- ProbTrace(x=pars[[i]], Theta=mtheta)            
                if(any(N * P < 2)){
                    df[i] <- df[i] - 1
                    next
                }                
                X2[i] <- X2[i] + sum((r - N*P)^2 / N*P)
            }
            df[i] <- df[i] + n.uniqueGroups*(length(r) - 1) - sum(pars[[i]]@est)            
        }
        if(!is.null(empirical.plot)) return(invisible(NULL))
        p <- pchisq(X2, df, lower.tail = FALSE)               
        ret <- data.frame(item=colnames(x@data), X2=X2, df=df, p=round(p,3)) 
        return(ret)
    }    
    if(type == 'Zh'){
        N <- nrow(Theta)
        itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)        
        for (i in 1:J){
            if(x@pars[[1]]@bfactor) pars[[i]]@bfactor <- FALSE
            itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)
        }
        LL <- itemtrace * fulldata        
        LL[LL < 1e-8] <- 1
        Lmatrix <- matrix(log(LL[as.logical(fulldata)]), N, J)              
        mu <- sigma2 <- rep(0, J)
        for(item in 1:J){              
            for(n in 1:N){                
                P <- itemtrace[n ,itemloc[item]:(itemloc[item+1]-1)]
                mu[item] <- mu[item] + sum(P * log(P))            
                for(i in 1:length(P)){
                    for(j in 1:length(P)){
                        if(i != j)
                            sigma2[item] <- sigma2[item] + P[i] * P[j] * log(P[i]) * log(P[i]/P[j])
                    }
                }                               
            }               
        }        
        Zh <- (colSums(Lmatrix) - mu) / sqrt(sigma2)
        p <- round(pnorm(Zh), 3)
        p[Zh > 0] <- 1 - p[Zh > 0]
        p <- 2*p
        ret <- data.frame(item=colnames(x@data), Zh=Zh, p=p)                
        return(ret)
    }    
}
