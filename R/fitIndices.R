#' Compute Extra Model Fit Indices
#' 
#' Compute additional model fit indecies that do not come as direct results following parameter
#' convergence. Will always compute the M2 (Maydeu-Olivares & Joe, 2006) statistic by default.
#' 
#' 
#'
#' @aliases fitIndices
#' @param obj an estimated model object from the mirt package
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Maydeu-Olivares, A. & Joe, H. (2006). Limited information goodness-of-fit testing in 
#' multidimensional contingency tables Psychometrika, 71, 713-732. 
#' @keywords model fit
#' @export fitIndices
#' @examples 
#' \dontrun{
#' #LSAT6 example
#' dat <- expand.table(LSAT6)
#' (mod1 <- mirt(dat, 1, itemtype = '1PL'))
#' fitIndices(mod1)
#' 
#' #Science data, much more sparce so M2 would be more informative
#' (mod2 <- mirt(Science, 1))
#' fitIndices(mod2)
#' }
fitIndices <- function(obj){
    #if MG loop    
    if(is(obj, 'MixedClass'))
        stop('mixedmirt objects not yet supported')       
    if(is(obj, 'MultipleGroupClass')){        
        cmods <- obj@cmods
        r <- obj@tabdata[, ncol(obj@tabdata)]
        ngroups <- length(cmods)
        ret <- vector('list', length(cmods))
        for(g in 1:ngroups){
            attr(cmods[[g]], 'MG') <- TRUE
            ret[[g]] <- fitIndices(cmods[[g]])
        }        
        newret <- list()
        newret$M2 <- numeric(ngroups)   
        names(newret$M2) <- obj@groupNames
        for(g in 1:ngroups)            
            newret$M2[g] <- ret[[g]]$M2
        newret$M2Total <- sum(newret$M2)
        Tsum <- 0
        for(g in 1:ngroups) Tsum <- Tsum + ret[[g]]$nrowT
        newret$dfM2 <- obj@df - (nrow(obj@tabdata) - Tsum)
        newret$p.M2 <- 1 - pchisq(newret$M2Total, newret$dfM2)
        newret$RMSEA.M2 <- ifelse((newret$M2Total - newret$dfM2) > 0, 
                           sqrt(newret$M2Total - newret$dfM2) / sqrt(newret$dfM2 * (sum(r)-1)), 0) 
        return(newret)
    }
    ret <- list()        
    tabdata <- obj@tabdatalong
    K <- obj@K
    nitems <- length(K)    
    r <- tabdata[, ncol(tabdata)]
    N <- sum(r)
    p <- r/N
    p_theta <- obj@Pl
    tabdata <- tabdata[, -ncol(tabdata)]
    itemloc <- obj@itemloc
    T <- matrix(NA, nrow(tabdata), nrow(tabdata))
    Gamma <- diag(p_theta) - outer(p_theta, p_theta)    
    ind <- 1
    
    ## M2 stat
    #find univariate marginals
    for(i in 1:nitems){
        for(j in 1:(K[i]-1)){
            loc <- itemloc[i] + j    
            T[ind, ] <- tabdata[, loc]            
            ind <- ind + 1
        }       
    }    
    #find bivariate marginals
    for(i in 1:nitems){
        for(j in 1:nitems){
            if(i < j){
                for(k1 in 1:(K[i]-1)){
                    for(k2 in 1:(K[j]-1)){
                        loc1 <- itemloc[i] + k1
                        loc2 <- itemloc[j] + k2
                        T[ind, ] <- as.numeric(tabdata[, loc1] & tabdata[, loc2])
                        ind <- ind + 1                        
                    }
                }
            }
        }
    }    
    T <- na.omit(T)
    Eta <- T %*% Gamma %*% t(T)
    T.p <- T %*% p
    T.p_theta <- T %*% p_theta       
    inv.Eta <- try(solve(Eta), silent = TRUE)
    if(is(inv.Eta, 'try-error')){        
        diag(Eta) <- diag(Eta) + .01 * diag(Eta)
        inv.Eta <- try(solve(Eta), silent = TRUE)
        if(is(inv.Eta, 'try-error'))
            stop('M2 cannot be computed')
    }
    pars <- obj@pars
    theta <- seq(-4, 4, length.out = 40)
    Theta <- thetaComb(theta, obj@nfact)
    gstructgrouppars <- ExtractGroupPars(pars[[nitems+1]])
    Prior <- mvtnorm::dmvnorm(Theta,gstructgrouppars$gmeans,
                                   gstructgrouppars$gcov)    
    itemloc <- obj@itemloc
    Prior <- Prior/sum(Prior)
    delta <- NULL       
    for(pat in 1:nrow(tabdata)){
        DX <- c()
        rlist <- Estep.mirt(pars=pars, tabdata=matrix(c(tabdata[pat, ], r[pat]), 1), 
                            Theta=Theta, prior=Prior, itemloc=itemloc, 
                            debug='NULL')      
        for(i in 1:nitems){                         
            tmp <- c(itemloc[i]:(itemloc[i+1] - 1))
            pars[[i]]@rs <- rlist$r1[, tmp]                       
            dx <- Deriv(pars[[i]], Theta=Theta, EM = TRUE, prior=Prior)$grad
            DX <- c(DX, dx[pars[[i]]@est])
        } 
        if(is.null(delta)) delta <- matrix(NA, nrow(tabdata), length(DX), byrow = TRUE)
        delta[pat, ] <- DX
    }
    delta2 <- T %*% delta    
    C2 <- inv.Eta - inv.Eta %*% delta2 %*% solve(t(delta2) %*% inv.Eta %*% delta2) %*% 
        t(delta2) %*% inv.Eta
    M2 <- N * t(T.p - T.p_theta) %*% C2 %*% (T.p - T.p_theta)
    ret$M2 <- M2  
    if(is.null(attr(obj, 'MG'))){
        ret$dfM2 <- obj@df - (nrow(tabdata) -  nrow(T))    
        ret$p.M2 <- 1 - pchisq(M2, ret$dfM2)
        ret$RMSEA.M2 <- ifelse((M2 - ret$dfM2) > 0, 
                        sqrt(M2 - ret$dfM2) / sqrt(ret$dfM2 * (sum(r)-1)), 0)                  
    } else {
        ret$nrowT <- nrow(T)        
    }
    ret    
}
