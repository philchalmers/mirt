#' Compute M2 statistic
#'
#' Computes the M2 (Maydeu-Olivares & Joe, 2006) statistic for dichotomous data and the 
#' M2* statistic for polytomous data (collapsing over response categories for better stability;
#' see Cai and Hansen, 2013), as well as associated fit indicies that are based on 
#' fitting the null model.
#'
#'
#' @aliases M2
#' @param obj an estimated model object from the mirt package
#' @param quadpts number of quadrature points to use during estimation. If \code{NULL}, \code{quadpts}
#'   will be extracted from the \code{obj}; if not available, a suitable value will be chosen based
#'   on the rubric found in \code{\link{mirt}}
#' @param calcNull logical; calculate statistics for the null model as well?
#'   Allows for statistics such as the limited information TLI and CFI
#' @param Theta a matrix of factor scores for each person used for imputation
#' @param impute a number indicating how many imputations to perform (passed to \code{\link{imputeMissing}})
#'   when there are missing data present. This requires a precomputed \code{Theta} input. Will return
#'   a data.frame object with the mean estimates of the stats and their imputed standard deviations
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Cai, L. & Hansen, M. (2013). Limited-information goodness-of-fit testing of 
#' hierarchical item factor models. British Journal of Mathematical and Statistical 
#' Psychology, 66, 245-276.
#' 
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
#' #M2 imputed with missing data present (run in parallel)
#' dat[sample(1:prod(dim(dat)), 250)] <- NA
#' mod2 <- mirt(dat, 1)
#' mirtCluster()
#' Theta <- fscores(mod2, full.scores=TRUE)
#' M2(mod2, Theta=Theta, impute = 10)
#'
#' }
M2 <- function(obj, calcNull = TRUE, quadpts = NULL, Theta = NULL, impute = 0){
    
    fn <- function(collect, obj, Theta, ...){
        dat <- imputeMissing(obj, Theta)
        tmpobj <- obj
        tmpobj@data <- dat
        if(is(obj, 'MultipleGroupClass')){
            large <- multipleGroup(dat, 1, group=obj@group, large = TRUE)
            for(g in 1L:length(obj@groupNames)){
                tmpobj@cmods[[g]]@data <- dat[obj@groupNames[g] == obj@group, , drop=FALSE]
                tmpobj@cmods[[g]]@tabdata <- large$tabdata2[[g]]
                tmpobj@cmods[[g]]@tabdatalong <- large$tabdata[[g]]
            }
        } else {
            large <- mirt(dat, 1, large = TRUE)
            tmpobj@tabdata <- large$tabdata2[[1L]]
            tmpobj@tabdatalong <- large$tabdata[[1L]]
        }
        return(M2(tmpobj, ...))
    }
    
    #if MG loop
    if(is(obj, 'MixedClass'))
        stop('mixedmirt objects not yet supported')
    if(any(is.na(obj@data))){
        if(impute == 0 || is.null(Theta))
            stop('Fit statistics cannot be computed when there are missing data. Pass suitable
                 Theta and impute arguments to compute statistics following multiple data inputations')
        collect <- vector('list', impute)
        collect <- myLapply(collect, fn, obj=obj, Theta=Theta, calcNull=calcNull,
                            quadpts=quadpts)
        ave <- SD <- collect[[1L]]
        ave[ave!= 0] <- SD[SD!=0] <- 0
        for(i in 1L:impute)
            ave <- ave + collect[[i]]
        ave <- ave/impute
        for(i in 1L:impute)
            SD <- (ave - collect[[i]])^2
        SD <- sqrt(SD/impute)
        ret <- rbind(ave, SD)
        rownames(ret) <- c('stats', 'SD_stats')
        return(ret)
    }
    if(is(obj, 'MultipleGroupClass')){
        cmods <- obj@cmods
        r <- obj@tabdata[, ncol(obj@tabdata)]
        ngroups <- length(cmods)
        ret <- vector('list', length(cmods))
        for(g in 1L:ngroups){
            attr(cmods[[g]], 'MG') <- g
            cmods[[g]]@bfactor <- obj@bfactor
            cmods[[g]]@quadpts <- obj@quadpts
            ret[[g]] <- M2(cmods[[g]], calcNull=FALSE, quadpts=quadpts)
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
            null.fit <- M2(null.mod, calcNull=FALSE, quadpts=quadpts)
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
    
    if(!all(sapply(obj@pars, class) %in% c('dich', 'graded', 'gpcm', 'nominal', 'ideal', 'GroupPars')))
       stop('M2 currently only supported for \'dich\', \'ideal\', \'graded\', 
            \'gpcm\', and \'nominal\' objects')
    ret <- list()
    group <- if(is.null(attr(obj, 'MG'))) 1 else attr(obj, 'MG')
    tabdata <- obj@tabdatalong
    nitems <- ncol(obj@data)
    if(any(is.na(obj@tabdata)))
        stop('M2 can not be calulated for data with missing values.')
    adj <- apply(obj@data, 2, min)
    dat <- t(t(obj@data) - adj)
    N <- nrow(dat)
    p  <- colMeans(dat)
    cross <- crossprod(dat, dat)
    p <- c(p, cross[lower.tri(cross)]/N)
    K <- obj@K
    pars <- obj@pars
    quadpts2 <- obj@quadpts    
    if(is.nan(quadpts2)) 
        quadpts2 <- switch(as.character(obj@nfact), '1'=41, '2'=21, '3'=11, '4'=7, '5'=5, 3)
    if(is.null(quadpts)) quadpts <- quadpts2
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
        Prior <- Prior <- mirt_dmvnorm(Theta,gstructgrouppars$gmeans,
                                           gstructgrouppars$gcov)
        Prior <- Prior/sum(Prior)
    } else {
        Theta <- obj@Theta        
        prior <- bfactorlist$prior[[group]]; Priorbetween <- bfactorlist$Priorbetween[[group]]
        sitems <- bfactorlist$sitems; specific <- bfactorlist$specific; 
        Prior <- bfactorlist$Prior[[group]]
    }
    E1 <- numeric(nitems)
    E2 <- matrix(NA, nitems, nitems)
    EIs <- EIs2 <- matrix(0, nrow(Theta), nitems)
    DP <- matrix(0, nrow(Theta), length(estpars))
    wherepar <- c(1L, numeric(nitems))
    ind <- 1L
    for(i in 1L:nitems){
        x <- extract.item(obj, i)
        EIs[,i] <- expected.item(x, Theta, min=0L)
        tmp <- ProbTrace(x, Theta)
        for(j in ncol(tmp):2L)
            tmp[,j-1L] <- tmp[,j] + tmp[,j-1L]
        cfs <- c(0,1)
        if(K[i] > 2L) cfs <- c(cfs, 2:(ncol(tmp)-1L) * 2 - 1)
        EIs2[,i] <- t(cfs %*% t(tmp))
        tmp <- length(x@parnum)
        DP[ ,ind:(ind+tmp-1L)] <- dP(x, Theta)
        ind <- ind + tmp
        wherepar[i+1L] <- ind
    }
    ind <- 1L
    for(i in 1L:nitems){
        E1[i] <- sum(EIs[,i] * Prior)
        for(j in 1L:nitems){
            if(i > j){
                E2[i,j] <- sum(EIs[,i] * EIs[,j] * Prior)
                ind <- ind + 1L
            }
        }
    }
    e <- c(E1, E2[lower.tri(E2)])
    delta1 <- matrix(0, nitems, length(estpars))
    delta2 <- matrix(0, length(p) - nitems, length(estpars))
    ind <- 1L
    offset <- pars[[1L]]@parnum[1L] - 1L
    for(i in 1L:nitems){
        dp <- colSums(DP[ , wherepar[i]:(wherepar[i+1L]-1L), drop=FALSE] * Prior)
        delta1[i, pars[[i]]@parnum - offset] <- dp
        for(j in 1L:nitems){
            if(i < j){
                dp <- colSums(DP[ , wherepar[i]:(wherepar[i+1L]-1L), drop=FALSE] * EIs[,j] * Prior)
                delta2[ind, pars[[i]]@parnum - offset] <- dp
                dp <- colSums(DP[ , wherepar[j]:(wherepar[j+1L]-1L), drop=FALSE] * EIs[,i] * Prior)
                delta2[ind, pars[[j]]@parnum - offset] <- dp
                ind <- ind + 1L
            }
        }
    }
    delta <- rbind(delta1, delta2)
    delta <- delta[, estpars, drop=FALSE]
    Xi2els <- .Call('buildXi2els', nrow(delta1), nrow(delta2), nitems, EIs, EIs2, Prior)
    Xi2 <- rbind(cbind(Xi2els$Xi11, Xi2els$Xi12), cbind(t(Xi2els$Xi12), Xi2els$Xi22))    
    tmp <- qr.Q(qr(delta), complete=TRUE)
    if((ncol(delta) + 1L) > ncol(tmp))
        stop('M2 cannot be calulated since df is too low')
    deltac <- tmp[,(ncol(delta) + 1L):ncol(tmp), drop=FALSE]
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
            null.fit <- M2(null.mod, calcNull=FALSE, quadpts=quadpts)
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
