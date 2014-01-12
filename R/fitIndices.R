#' Compute Extra Model Fit Indices
#'
#' Compute additional model fit indices that do not come as direct results following parameter
#' convergence. Will compute the M2 (Maydeu-Olivares & Joe, 2006) statistic by default, and
#' returns a data.frame containing various model fit statistics.
#'
#'
#' @aliases fitIndices
#' @param obj an estimated model object from the mirt package
#' @param calcNull logical; calculate statistics for the null model as well?
#' Allows for statistics such as the limited information TLI and CFI
#' @param prompt logical; prompt user for input if the internal matrices are too large?
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
#' (mod1 <- mirt(dat, 1, itemtype = '2PL', constrain = list(c(1,5,9,13,17))))
#' fitIndices(mod1)
#'
#' #Science data with computing the null model M2 stat
#' (mod2 <- mirt(Science, 1))
#' fitIndices(mod2, calcNull = TRUE)
#' }
fitIndices <- function(obj, calcNull = FALSE, prompt = TRUE){
    #if MG loop
    if(is(obj, 'MixedClass'))
        stop('mixedmirt objects not yet supported')
    if(is(obj, 'MultipleGroupClass')){
        cmods <- obj@cmods
        r <- obj@tabdata[, ncol(obj@tabdata)]
        ngroups <- length(cmods)
        ret <- vector('list', length(cmods))
        for(g in 1L:ngroups){
            attr(cmods[[g]], 'MG') <- TRUE
            ret[[g]] <- fitIndices(cmods[[g]])
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
            null.fit <- fitIndices(null.mod, prompt=FALSE)
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
    p_theta <- p_theta
    tabdata <- tabdata[, -ncol(tabdata)]
    itemloc <- obj@itemloc
    T <- matrix(NA, sum(K) + sum(K*(sum(K))), nrow(tabdata))
    Gamma <- diag(p_theta) - outer(p_theta, p_theta)
    ind <- 1L

    ## M2 stat
    #find univariate marginals
    for(i in 1L:nitems){
        for(j in 1L:(K[i]-1L)){
            loc <- itemloc[i] + j
            T[ind, ] <- tabdata[, loc]
            ind <- ind + 1L
        }
    }
    #find bivariate marginals
    for(i in 1L:nitems){
        for(j in 1L:nitems){
            if(i < j){
                for(k1 in 1L:(K[i]-1L)){
                    for(k2 in 1L:(K[j]-1L)){
                        loc1 <- itemloc[i] + k1
                        loc2 <- itemloc[j] + k2
                        T[ind, ] <- as.integer(tabdata[, loc1] & tabdata[, loc2])
                        ind <- ind + 1L
                    }
                }
            }
        }
    }
    T <- T[1L:(ind-1L), ]
    if(nrow(T) > 4000L && prompt){
        cat('Internal matricies are very large and computations will therefore take an extended
            amount of time and require large amounts of RAM. The largest matrix has', nrow(T), 'columns.
            Do you wish to continue anyways?')
        input <- readline("(yes/no): ")
        if(input == 'no') stop('Execution halted.')
        if(input != 'yes') stop('Illegal user input')
    }
    Eta <- T %*% Gamma %*% t(T)
    T.p <- T %*% p
    T.p_theta <- T %*% p_theta
    inv.Eta <- ginv(Eta)
    pars <- obj@pars
    quadpts <- switch(as.character(obj@nfact), '1'=41, '2'=21, '3'=11, '4'=7, '5'=5, 3)
    theta <- as.matrix(seq(-(.8 * sqrt(quadpts)), .8 * sqrt(quadpts), length.out = quadpts))
    Theta <- thetaComb(theta, obj@nfact)
    gstructgrouppars <- ExtractGroupPars(pars[[nitems+1L]])
    Prior <- mvtnorm::dmvnorm(Theta,gstructgrouppars$gmeans,
                                   gstructgrouppars$gcov)
    itemloc <- obj@itemloc
    Prior <- Prior/sum(Prior)
    whichpar <- integer(nitems)
    for(i in 1L:nitems)
        whichpar[i] <- sum(pars[[i]]@est)
    npick <- sum(whichpar)
    whichpar <- c(0L, cumsum(whichpar)) + 1L
    delta <- matrix(NA, nrow(tabdata), npick, byrow = TRUE)
    DX <- rep(NA, npick)
    for(pat in 1L:nrow(tabdata)){
        rlist <- Estep.mirt(pars=pars, tabdata=matrix(c(tabdata[pat, ], r[pat]), 1),
                            Theta=Theta, prior=Prior, itemloc=itemloc, deriv=TRUE)
        for(i in 1L:nitems){
            tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
            pars[[i]]@dat <- rlist$r1[, tmp]
            pars[[i]]@itemtrace <- rlist$itemtrace[, tmp]
            dx <- Deriv(pars[[i]], Theta=Theta, EM = TRUE, estHess=FALSE)$grad
            DX[whichpar[i]:(whichpar[i+1L]-1L)] <- dx[pars[[i]]@est]
        }
        delta[pat, ] <- DX
    }
    delta2 <- T %*% delta
    delta2.invEta.delta2 <- t(delta2) %*% inv.Eta %*% delta2
    C2 <- inv.Eta - inv.Eta %*% delta2 %*% solve(delta2.invEta.delta2) %*%
        t(delta2) %*% inv.Eta
    M2 <- N * t(T.p - T.p_theta) %*% C2 %*% (T.p - T.p_theta)
    ret$M2 <- M2
    if(is.null(attr(obj, 'MG'))){
        ret$df.M2 <- nrow(T) - obj@nest
        ret$p.M2 <- 1 - pchisq(M2, ret$df.M2)
        ret$RMSEA.M2 <- ifelse((M2 - ret$df.M2) > 0,
                        sqrt(M2 - ret$df.M2) / sqrt(ret$df.M2 * (sum(r)-1)), 0)
        if(calcNull){
            null.mod <- try(mirt(obj@data, 1, TOL=1e-3, technical=list(NULL.MODEL=TRUE),
                                 verbose=FALSE))
            null.fit <- fitIndices(null.mod, prompt=FALSE)
            ret$TLI.M2 <- (null.fit$M2 / null.fit$df.M2 - ret$M2/ret$df.M2) /
                (null.fit$M2 / null.fit$df.M2 - 1)
            ret$CFI.M2 <- 1 - (ret$M2 - ret$df.M2) / (null.fit$M2 - null.fit$df.M2)
            if(ret$CFI.M2 > 1) ret$CFI.M2 <- 1
            if(ret$CFI.M2 < 0) ret$CFI.M2 <- 0
        }
    } else {
        ret$nrowT <- nrow(T)
    }
    return(as.data.frame(ret))
}
