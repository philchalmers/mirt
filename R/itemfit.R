#' Item fit statistics
#' 
#' \code{itemfit} calculates the Zh values from Drasgow, Levine and Williams (1985) 
#' and \eqn{\chi^2} values for unidimensional models. For Rasch, partial credit, and rating scale models 
#' infit and outfit statistics are also produced.
#' 
#' @aliases itemfit
#' @param x a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or 
#' \code{MultipleGroupClass}
#' @param X2 logical; calculate the X2 statistic for unidimensional models?
#' @param group.size approximate size of each group to be used in calculating the \eqn{\chi^2} statistic
#' @param empirical.plot a single numeric value or character of the item name  indicating which item to plot
#'  (via \code{itemplot}) and
#' overlay with the empirical \eqn{\theta} groupings. Only applicable when \code{type = 'X2'}. 
#' The default is \code{NULL}, therefore no plots are drawn 
#' @param method type of factor score estimation method. Can be expected
#' a-posteriori (\code{"EAP"}), Bayes modal (\code{"MAP"}), weighted likelihood estimation 
#' (\code{"WLE"}), or maximum likelihood (\code{"ML"}) 
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
#' Reise, S. P. (1990). A comparison of item- and person-fit methods of assessing model-data fit 
#' in IRT. \emph{Applied Psychological Measurement, 14}, 127-137.
#' 
#' @examples
#' 
#' \dontrun{
#'
#' #make some data
#' set.seed(1234)
#' a <- matrix(rlnorm(20, meanlog=0, sdlog = .1),ncol=1)
#' d <- matrix(rnorm(20),ncol=1)
#' items <- rep('dich', 20)
#' data <- simdata(a,d, 2000, items)
#'  
#' x <- mirt(data, 1)
#' raschfit <- mirt(data, 1, itemtype='Rasch')
#' fit <- itemfit(x)
#' fit
#' 
#' itemfit(x, empirical.plot = 1) #empirical item plot
#' itemfit(raschfit, method = 'ML') #infit and outfit stats (method='ML' agrees better with eRm package)
#' 
#'   }
#'
itemfit <- function(x, X2 = FALSE, group.size = 150, empirical.plot = NULL, method = 'EAP'){    
    if(any(is.na(x@data)))
        stop('Fit statistics cannot be computed when there are missing data.')
    if(is(x, 'MultipleGroupClass')){
        ret <- list()   
        for(g in 1:length(x@cmods))
            ret[[g]] <- itemfit(x@cmods[[g]], group.size=group.size)
        names(ret) <- names(x@cmods)
        return(ret)
    }        
    sc <- fscores(x, verbose = FALSE, full.scores = TRUE) 
    J <- ncol(x@data)
    itemloc <- x@itemloc    
    pars <- x@pars
    prodlist <- attr(pars, 'prodlist')    
    nfact <- x@nfact + length(prodlist)
    fulldata <- x@fulldata    
    Theta <- sc[ ,ncol(sc):(ncol(sc) - nfact + 1), drop = FALSE]        
    N <- nrow(Theta)
    itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)        
    for (i in 1:J)            
        itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)        
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
    #if all Rasch models, infit and outfit        
    if(all(x@itemtype %in% c('Rasch', 'rsm', 'gpcm'))){
        for(i in 1:length(x@itemtype))
            if((x@pars[[i]]@par[1] * x@pars[[1]]@D) != 1) break
        attr(x, 'inoutfitreturn') <- TRUE
        pf <- personfit(x, method=method)        
        outfit <- colSums(pf$resid2 / pf$info) / N
        infit <- colSums(pf$resid2) / colSums(pf$info)        
        ret <- data.frame(item=colnames(x@data), outfit=outfit, infit=infit, Zh=Zh)        
    } else {
        ret <- data.frame(item=colnames(x@data), Zh=Zh)    
    }    
    if((X2 || !is.null(empirical.plot)) && nfact == 1){                
        ord <- order(Theta[,1])    
        fulldata <- fulldata[ord,]
        Theta <- Theta[ord, , drop = FALSE]
        den <- dnorm(Theta, 0, .5)
        den <- den / sum(den)
        cumTheta <- cumsum(den)
        Groups <- rep(20, length(ord))
        ngroups <- ceiling(nrow(fulldata) / group.size)        
        weight <- 1/ngroups        
        for(i in 1:20)
            Groups[round(cumTheta,2) >= weight*(i-1) & round(cumTheta,2) < weight*i] <- i        
        n.uniqueGroups <- length(unique(Groups))
        X2 <- df <- RMSEA <- rep(0, J)                    
        if(!is.null(empirical.plot)){
            if(nfact > 1) stop('Cannot make empirical plot for multidimensional models')
            theta <- seq(-4,4, length.out=40)
            ThetaFull <- thetaComb(theta, nfact)
            if(!is.numeric(empirical.plot)){
                inames <- colnames(x@data)
                ind <- 1:length(inames)
                empirical.plot <- ind[inames == empirical.plot]              
            }
            P <- ProbTrace(pars[[empirical.plot]], ThetaFull)            
            plot(ThetaFull, P[,1], type = 'l', ylim = c(0,1), las = 1, 
                 main =paste('Item', empirical.plot), ylab = expression(P(theta)), xlab = expression(theta))
            for(i in 2:ncol(P))
                lines(ThetaFull, P[,i], col = i)
        }
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
                    points(rep(mtheta, length(tmp) - 1), tmp[-1], col = col, pch = col+2)                  
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
        X2[X2 == 0] <- NA
        if(!is.null(empirical.plot)) return(invisible(NULL))
        ret$df <- df
        ret$X2 <- X2
    }                    
    return(ret)
}
