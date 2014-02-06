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
fitIndices <- function(obj, calcNull = FALSE, collapse_poly = FALSE, prompt = TRUE){
    #if MG loop
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
            ret[[g]] <- fitIndices(cmods[[g]], prompt = if(g == 1L) prompt else FALSE)
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
    group <- if(is.null(attr(obj, 'MG'))) 1 else attr(obj, 'MG')
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
    p_theta <- p_theta / sum(p_theta)
    tabdata <- tabdata[, -ncol(tabdata)]
    itemloc <- obj@itemloc
    Tmat <- matrix(NA, sum(K-1L) + sum((K-1L)*(sum(K-1L))), nrow(tabdata))
    if(collapse_poly)
        Pmat <- matrix(0L, nitems*(nitems+1L)/2, sum(K-1L) + sum((K-1L)*(sum(K-1L))))
    Gamma <- diag(p_theta) - outer(p_theta, p_theta)
    ind <- ind2 <- 1L
    #find univariate marginals
    for(i in 1L:nitems){
        if(collapse_poly){
            Pmat[i, ind2:(ind2+K[i]-2L)] <- 2L:K[i] - 1L
            ind2 <- ind2 + K[i] - 1L
        }
        for(j in 1L:(K[i]-1L)){
            loc <- itemloc[i] + j
            Tmat[ind, ] <- as.integer(tabdata[, loc])
            ind <- ind + 1L
        }
    }
    #find bivariate marginals
    ind1 <- nitems + 1L
    for(i in 1L:nitems){
        for(j in 1L:nitems){
            if(i > j){
                if(collapse_poly){
                    tmp <- kronecker(2L:K[i] - 1L, 2L:K[j] - 1L)
                    Pmat[ind1, ind2:(ind2+length(tmp)-1L)] <- tmp
                    ind1 <- ind1 + 1L
                    ind2 <- ind2 + length(tmp)
                }
                for(k1 in 1L:(K[i]-1L)){
                    for(k2 in 1L:(K[j]-1L)){
                        loc1 <- itemloc[i] + k1
                        loc2 <- itemloc[j] + k2
                        Tmat[ind, ] <- as.integer(tabdata[, loc1] & tabdata[, loc2])
                        ind <- ind + 1L
                    }
                }
            }
        }
    }
    Tmat <- Tmat[1L:(ind-1L), ]
    if(collapse_poly){
        Pmat <- Pmat[ , which(colSums(Pmat) != 0L)]
        if(is.numeric(collapse_poly)) return(list(T=Tmat, P=Pmat))
        Tmat <- Pmat %*% Tmat
    }
    if(nrow(Tmat) > 4000L){
        if(prompt){
            cat('Internal matricies are very large and computations will therefore take an extended
                amount of time and require large amounts of RAM. The largest matrix has', nrow(Tmat), 'columns.
                Do you wish to continue anyways?')
            input <- readline("(yes/no): ")
            if(input == 'no') stop('Execution halted.')
            if(input != 'yes') stop('Illegal user input')
        }
    }
    Eta <- Tmat %*% Gamma %*% t(Tmat)
    T.p <- Tmat %*% p
    T.p_theta <- Tmat %*% p_theta
    inv.Eta <- ginv(Eta)
    pars <- obj@pars
    quadpts <- obj@quadpts
    if(is.nan(quadpts)) 
        quadpts <- switch(as.character(obj@nfact), '1'=41, '2'=21, '3'=11, '4'=7, '5'=5, 3)
    itemloc <- obj@itemloc
    bfactorlist <- obj@bfactor
    theta <- as.matrix(seq(-(.8 * sqrt(quadpts)), .8 * sqrt(quadpts), length.out = quadpts))
    if(is.null(bfactorlist$Priorbetween[[1L]])){
        Theta <- thetaComb(theta, obj@nfact)
        prior <- Priorbetween <- sitems <- specific <- NULL
    } else {
        Theta <- obj@Theta        
        prior <- bfactorlist$prior[[group]]; Priorbetween <- bfactorlist$Priorbetween[[group]]
        sitems <- bfactorlist$sitems; specific <- bfactorlist$specific
    }
    gstructgrouppars <- ExtractGroupPars(pars[[nitems+1L]])
    Prior <- Prior <- mvtnorm::dmvnorm(Theta,gstructgrouppars$gmeans,
                                       gstructgrouppars$gcov)
    Prior <- Prior/sum(Prior)
    whichpar <- integer(nitems)
    for(i in 1L:nitems)
        whichpar[i] <- sum(pars[[i]]@est)
    npick <- sum(whichpar)
    whichpar <- c(0L, cumsum(whichpar)) + 1L
    delta <- matrix(NA, nrow(tabdata), npick, byrow = TRUE)
    DX <- rep(NA, npick)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, 
                                  CUSTOM.IND=obj@CUSTOM.IND)
    for(pat in 1L:nrow(tabdata)){
        if(is.null(prior)){
            rlist <- Estep.mirt(pars=pars, tabdata=matrix(c(tabdata[pat, ], r[pat]), 1),
                                Theta=Theta, prior=Prior, itemloc=itemloc, deriv=TRUE,
                                CUSTOM.IND=obj@CUSTOM.IND, itemtrace=itemtrace)
        } else {
            rlist <- Estep.bfactor(pars=pars, tabdata=matrix(c(tabdata[pat, ], r[pat]), 1),
                                   Theta=Theta, prior=prior, Prior=Prior,
                                   Priorbetween=Priorbetween, specific=specific, 
                                   sitems=sitems, itemloc=itemloc, CUSTOM.IND=obj@CUSTOM.IND, 
                                   itemtrace=itemtrace)
        }
        for(i in 1L:nitems){
            tmp <- c(itemloc[i]:(itemloc[i+1L] - 1L))
            pars[[i]]@dat <- rlist$r1[, tmp]
            dx <- Deriv(pars[[i]], Theta=Theta, estHess=FALSE)$grad
            DX[whichpar[i]:(whichpar[i+1L]-1L)] <- dx[pars[[i]]@est]
        }
        delta[pat, ] <- DX
    }
    delta2 <- Tmat %*% delta
    delta2.invEta.delta2 <- t(delta2) %*% inv.Eta %*% delta2
    C2 <- inv.Eta - inv.Eta %*% delta2 %*% solve(delta2.invEta.delta2) %*%
        t(delta2) %*% inv.Eta
    M2 <- N * t(T.p - T.p_theta) %*% C2 %*% (T.p - T.p_theta)
    ret$M2 <- M2
    if(is.null(attr(obj, 'MG'))){
        ret$df.M2 <- nrow(Tmat) - obj@nest
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
        ret$nrowT <- nrow(Tmat)
    }
    return(as.data.frame(ret))
}
