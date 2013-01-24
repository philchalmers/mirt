#' Compute Extra Model Fit Indices
#' 
#' Compute additional model fit indecies that do not come as direct results following parameter
#' convergence. Will only compute the M2 (Maydeu-Olivares & Joe, 2006) statistic by default, and 
#' returns a list containing the requested statistics.
#' 
#' 
#'
#' @aliases fitIndices
#' @param obj an estimated model object from the mirt package
#' @param prompt logical; prompt user for input if the internal matricies are too large?
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
fitIndices <- function(obj, prompt = TRUE){
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
        newret$df.M2 <- obj@df - (nrow(obj@tabdata) - Tsum) + 1
        newret$p.M2 <- 1 - pchisq(newret$M2Total, newret$df.M2)
        newret$RMSEA.M2 <- ifelse((newret$M2Total - newret$df.M2) > 0, 
                           sqrt(newret$M2Total - newret$df.M2) / sqrt(newret$df.M2 * (sum(r)-1)), 0) 
        return(newret)
    }        
    ret <- list()        
    tabdata <- obj@tabdatalong
    if(any(is.na(obj@tabdata)))
        stop('M2 can not be calulated for data with missing values.')
    NOROWNA <- rowSums(is.na(obj@tabdata)) == 0
    tabdata <- tabdata[NOROWNA, ]
    K <- obj@K
    nitems <- length(K)    
    r <- tabdata[, ncol(tabdata)]
    N <- sum(r)
    p <- r/N
    p_theta <- obj@Pl[NOROWNA]
    p_theta <- p_theta/sum(p_theta)
    tabdata <- tabdata[, -ncol(tabdata)]
    itemloc <- obj@itemloc
    T <- matrix(NA, sum(K) + sum(K*(sum(K))), nrow(tabdata))
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
    T <- T[1:(ind-1), ]
    if(nrow(T) > 4000 && prompt){ 
        cat('Internal matricies are very large and computations will therefore take an extended 
            amount of time and require large amounts of RAM. The largest matrix has', nrow(T), 'columns. 
            Do you wish to continue anyways?')
        input <- readline("(yes/no): ")
        if(input == 'no' || input == 'n') stop('Execution halted.')
        if(input != 'yes' || input != 'y') stop('Illegal user input')        
    }        
    Eta <- T %*% Gamma %*% t(T)
    T.p <- T %*% p
    T.p_theta <- T %*% p_theta       
    Etarank <- qr(Eta)$rank
    while(Etarank < ncol(Eta)){
        diag(Eta) <- diag(Eta) + .001 * diag(Eta)
        Etarank <- qr(Eta)$rank
    }
    inv.Eta <- solve(Eta)
    pars <- obj@pars
    quadpts <- ceiling(40/(obj@nfact^1.5))
    theta <- seq(-4, 4, length.out = quadpts)
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
    deltarank <- qr(delta)$rank
    while(deltarank < ncol(delta)){
        diag(delta) <- diag(delta) + .001 * diag(delta)
        deltarank <- qr(delta)$rank      
    }    
    delta2 <- T %*% delta    
    C2 <- inv.Eta - inv.Eta %*% delta2 %*% solve(t(delta2) %*% inv.Eta %*% delta2) %*% 
        t(delta2) %*% inv.Eta
    M2 <- N * t(T.p - T.p_theta) %*% C2 %*% (T.p - T.p_theta)
    ret$M2 <- M2      
    if(is.null(attr(obj, 'MG'))){        
        ret$df.M2 <- obj@df - (nrow(tabdata) -  nrow(T)) + 1  
        ret$p.M2 <- 1 - pchisq(M2, ret$df.M2)
        ret$RMSEA.M2 <- ifelse((M2 - ret$df.M2) > 0, 
                        sqrt(M2 - ret$df.M2) / sqrt(ret$df.M2 * (sum(r)-1)), 0)                  
    } else {
        ret$nrowT <- nrow(T)        
    }
    ret    
}
