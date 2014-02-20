#' Compute M2 statistic
#'
#' Computes the M2 (Maydeu-Olivares & Joe, 2006) statistic and associated fit indicies.
#' For now, only dichotomous models are supported.
#'
#'
#' @aliases M2
#' @param obj an estimated model object from the mirt package
#' @param calcNull logical; calculate statistics for the null model as well?
#'   Allows for statistics such as the limited information TLI and CFI
#' @param collapse_poly logical; collapse across polytomous item categories to reduce 
#'   sparceness? Will also helo to reduce the internal matrix sizes. THIS FEATURE IS 
#'   CURRENTLY EXPERIMENTAL AND SHOULD NOT BE TRUSTED
#' @param prompt logical; prompt user for input if the internal matrices are too large?
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Maydeu-Olivares, A. & Joe, H. (2006). Limited information goodness-of-fit testing in
#' multidimensional contingency tables Psychometrika, 71, 713-732.
#' @keywords model fit
#' @export M2
#' @examples
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' (mod1 <- mirt(dat, 1))
#' M2(mod1)
#'
# #Science data with computing the null model M2 stat
# (mod2 <- mirt(Science, 1))
# M2(mod2, calcNull = TRUE)
#' }
M2 <- function(obj, calcNull = FALSE){
    
    #if MG loop
    collapse_poly <- TRUE
    if(is(obj, 'MixedClass'))
        stop('mixedmirt objects not yet supported')
    if(is(obj, 'MultipleGroupClass')){
        cmods <- obj@cmods
        r <- obj@tabdata[, ncol(obj@tabdata)]
        ngroups <- length(cmods)
        ret <- vector('list', length(cmods))
        for(g in 1L:ngroups){
            attr(cmods[[g]], 'MG') <- g
            cmods[[g]]@bfactor <- obj@bfactor
            cmods[[g]]@quadpts <- obj@quadpts
            ret[[g]] <- M2(cmods[[g]], calcNull=FALSE)
        }
        newret <- list()
        newret$M2 <- numeric(ngroups)
        names(newret$M2) <- obj@groupNames
        for(g in 1L:ngroups)
            newret$M2[g] <- ret[[g]]$M2
        newret$Total.M2 <- sum(newret$M2)
        Tsum <- 0
        for(g in 1L:ngroups) Tsum <- Tsum + ret[[g]]$nrowT
        newret$df.M2 <- Tsum - obj@nest
        newret$p.M2 <- 1 - pchisq(newret$Total.M2, newret$df.M2)
        newret$RMSEA.M2 <- ifelse((newret$Total.M2 - newret$df.M2) > 0,
                                  sqrt(newret$Total.M2 - newret$df.M2) / sqrt(newret$df.M2 * (sum(r)-1)), 0)
        if(calcNull){
            null.mod <- try(multipleGroup(obj@data, 1, group=obj@group, TOL=1e-3, technical=list(NULL.MODEL=TRUE),
                                          verbose=FALSE))
            null.fit <- M2(null.mod)
            newret$TLI.M2 <- (null.fit$Total.M2 / null.fit$df.M2 - newret$Total.M2/newret$df.M2) /
                (null.fit$Total.M2 / null.fit$df.M2 - 1)
            newret$CFI.M2 <- 1 - (newret$Total.M2 - newret$df.M2) / (null.fit$Total.M2 - null.fit$df.M2)
            if(newret$CFI.M2 > 1) newret$CFI.M2 <- 1
            if(newret$CFI.M2 < 0 ) newret$CFI.M2 <- 0
        }
        M2s <- as.numeric(newret$M2)
        names(M2s) <- paste0(obj@groupNames, '.M2')
        newret$M2 <- NULL
        return(data.frame(as.list(M2s), newret))
    }
    
    if(!all(sapply(obj@pars, class) %in% c('dich', 'GroupPars')))
       stop('M2 currently only supported for dichotomous objects')
    ret <- list()
    group <- if(is.null(attr(obj, 'MG'))) 1 else attr(obj, 'MG')
    tabdata <- obj@tabdatalong
    nitems <- ncol(obj@data)
    if(any(is.na(obj@tabdata)))
        stop('M2 can not be calulated for data with missing values.')
    adj <- apply(obj@data, 2, min)
    dat <- t(t(obj@data) - adj)
    N <- nrow(dat)
    if(!collapse_poly){
        dat <- expand.table(tabdata)
        dat <- dat[,-obj@itemloc]
    }
    p  <- colMeans(dat)
    cross <- crossprod(dat, dat)
    p <- c(p, cross[lower.tri(cross)]/N)
    p <- p[p != 0]
    K <- obj@K
    pars <- obj@pars
    quadpts <- obj@quadpts    
    if(is.nan(quadpts)) 
        quadpts <- switch(as.character(obj@nfact), '1'=41, '2'=21, '3'=11, '4'=7, '5'=5, 3)
    estpars <- c()
    for(i in 1L:(nitems+1L))
        estpars <- c(estpars, pars[[i]]@est)
    itemloc <- obj@itemloc
    bfactorlist <- obj@bfactor
    theta <- as.matrix(seq(-(.8 * sqrt(quadpts)), .8 * sqrt(quadpts), length.out = quadpts))
    if(is.null(bfactorlist$Priorbetween[[1L]])){
        Theta <- thetaComb(theta, obj@nfact)
        prior <- Priorbetween <- sitems <- specific <- NULL
        gstructgrouppars <- ExtractGroupPars(pars[[nitems+1L]])
        Prior <- Prior <- mvtnorm::dmvnorm(Theta,gstructgrouppars$gmeans,
                                           gstructgrouppars$gcov)
        Prior <- Prior/sum(Prior)
    } else {
        Theta <- obj@Theta        
        prior <- bfactorlist$prior[[group]]; Priorbetween <- bfactorlist$Priorbetween[[group]]
        sitems <- bfactorlist$sitems; specific <- bfactorlist$specific; Prior <- obj@Prior
    }
    if(collapse_poly){
        E1 <- numeric(nitems)
        E2 <- matrix(NA, nitems, nitems)
        names <- integer(nitems * (nitems -1)/2)
        names <- rbind(names, names)
        ind <- 1
        for(i in 1L:nitems){
            x <- extract.item(obj, i)
            Ex <- expected.item(x, Theta, min=0L)
            E1[i] <- sum(Ex * Prior)
            for(j in 1L:nitems){
                if(i > j){
                    y <- extract.item(obj, j)
                    Ey <- expected.item(y, Theta, min=0L)
                    E2[i,j] <- sum(Ex * Ey * Prior)
                    names[1,ind] <- i
                    names[2,ind] <- j
                    ind <- ind+1
                }
            }
        }
        e <- c(E1, E2[lower.tri(E2)])
        names(e) <- c(paste0('pi.', 1L:nitems), paste0('pi.', names[1L,], names[2L,]))
    } else {
        browser()
        #TODO direct method for polytomous items
    }
    
    delta1 <- matrix(0, nitems, length(estpars))
    delta2 <- matrix(0, length(p) - nitems, length(estpars))
    ind <- 1L
    offset <- pars[[1L]]@parnum[1L] - 1L
    for(i in 1L:nitems){
        x <- extract.item(obj, i)
        dp <- dP(x, Theta, Prior)
        delta1[i, pars[[i]]@parnum - offset] <- dp
        for(j in 1L:nitems){
            if(i < j){ 
                y <- extract.item(obj, j)
                P <- ProbTrace(y, Theta)[,2]
                dp <- dP(x, Theta, Prior, extra_term=P)
                delta2[ind, pars[[i]]@parnum - offset] <- dp
                
                P <- ProbTrace(x, Theta)[,2]
                dp <- dP(y, Theta, Prior, extra_term=P)
                delta2[ind, pars[[j]]@parnum - offset] <- dp
                ind <- ind + 1L
            }
        }
    }
    delta <- rbind(delta1, delta2)
    delta <- delta[, estpars, drop=FALSE]
    
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, 
                                  CUSTOM.IND=obj@CUSTOM.IND)
    itemtrace <- itemtrace[,-itemloc]
    Xi11 <- matrix(NA, nrow(delta1), nrow(delta1))
    Xi12 <- matrix(NA, nrow(delta1), nrow(delta2))
    Xi22 <- matrix(NA, nrow(delta2), nrow(delta2))    
    for(i in 1L:nrow(Xi11)){
        for(j in 1L:nrow(Xi11)){    
            if(i >= j){
                pab <- sum(itemtrace[,i] * itemtrace[,j] * Prior)
                pa <- sum(itemtrace[,i] * Prior)
                pb <- sum(itemtrace[,j] * Prior)
                if(i == j) pab <- pa
                Xi11[i,j] <- Xi11[j,i] <- pab - pa*pb
            }
        }
    }
    for(k in 1L:nitems){
        ind <- 1L
        for(i in 1L:nitems){
            for(j in 1L:nitems){
                if(i < j){
                    pabc <- sum(itemtrace[,i] * itemtrace[,j] * itemtrace[,k] * Prior)
                    pab <- sum(itemtrace[,i] * itemtrace[,j] * Prior)
                    pc <- sum(itemtrace[,k] * Prior)
                    if(i == k || j == k)
                        pabc <- sum(itemtrace[,i] * itemtrace[,j] * Prior)
                    Xi12[k, ind] <- pabc - pab*pc
                    ind <- ind + 1L
                }
            }
        }
    }   
    ind1 <- 1
    for(k in 1L:nitems){
        for(l in 1L:nitems){
            if(k < l){
                ind2 <- 1L
                for(i in 1L:nitems){
                    for(j in 1L:nitems){
                        if(i < j){
                            pabcd <- sum(itemtrace[,i] * itemtrace[,j] * 
                                             itemtrace[,k] * itemtrace[,l] * Prior)
                            pab <- sum(itemtrace[,i] * itemtrace[,j] * Prior)
                            pcd <- sum(itemtrace[,k] * itemtrace[,l] * Prior)
                            if(all(sort(c(i,j)) == sort(c(k,l)))){
                                pabcd <- pab
                            } else if(i == k || j == k){
                                pabcd <- sum(itemtrace[,i]*itemtrace[,j]*itemtrace[,l]*Prior)
                            } else if(i == l || j == l){
                                pabcd <- sum(itemtrace[,i]*itemtrace[,j]*itemtrace[,k]*Prior)
                            }
                            Xi22[ind1, ind2] <- pabcd - pab*pcd
                            ind2 <- ind2 + 1L
                        }
                    }
                }
                ind1 <- ind1 + 1L
            }
        }
    }
    Xi2 <- rbind(cbind(Xi11, Xi12), cbind(t(Xi12), Xi22))    
    tmp <- qr.Q(qr(delta), complete=TRUE)
    deltac <- tmp[,(ncol(delta) + 1L):ncol(tmp)]
    C2 <- deltac %*% solve(t(deltac) %*% Xi2 %*% deltac) %*% t(deltac)
    M2 <- N * t(p - e) %*% C2 %*% (p - e)
    ret$M2 <- M2
    if(is.null(attr(obj, 'MG'))){
        df <- length(p) - obj@nest
        ret$df.M2 <- df
        ret$p.M2 <- 1 - pchisq(M2, ret$df.M2)
        ret$RMSEA.M2 <- ifelse((M2 - ret$df.M2) > 0,
                               sqrt(M2 - ret$df.M2) / sqrt(ret$df.M2 * (N-1)), 0)
        if(calcNull){
            null.mod <- try(mirt(obj@data, 1, TOL=1e-3, technical=list(NULL.MODEL=TRUE),
                                 verbose=FALSE))
            null.fit <- M2(null.mod)
            ret$TLI.M2 <- (null.fit$M2 / null.fit$df.M2 - ret$M2/ret$df.M2) /
                (null.fit$M2 / null.fit$df.M2 - 1)
            ret$CFI.M2 <- 1 - (ret$M2 - ret$df.M2) / (null.fit$M2 - null.fit$df.M2)
            if(ret$CFI.M2 > 1) ret$CFI.M2 <- 1
            if(ret$CFI.M2 < 0) ret$CFI.M2 <- 0
        }
    } else {
        ret$nrowT <- length(p)
    }
    return(as.data.frame(ret))
}
