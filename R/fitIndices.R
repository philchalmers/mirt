# Compute Extra Model Fit Indices
#
# Compute additional model fit indices that do not come as direct results following parameter
# convergence. Will compute the M2 (Maydeu-Olivares & Joe, 2006) statistic by default, and
# returns a data.frame containing various model fit statistics.
#
#
# @aliases fitIndices
# @param obj an estimated model object from the mirt package
# @param calcNull logical; calculate statistics for the null model as well?
#   Allows for statistics such as the limited information TLI and CFI
# @param collapse_poly logical; collapse across polytomous item categories to reduce 
#   sparceness? Will also helo to reduce the internal matrix sizes. THIS FEATURE IS 
#   CURRENTLY EXPERIMENTAL AND SHOULD NOT BE TRUSTED
# @param prompt logical; prompt user for input if the internal matrices are too large?
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
# @references
# Maydeu-Olivares, A. & Joe, H. (2006). Limited information goodness-of-fit testing in
# multidimensional contingency tables Psychometrika, 71, 713-732.
# @keywords model fit
# @export fitIndices
# @examples
# \dontrun{
# #LSAT6 example
# dat <- expand.table(LSAT6)
# (mod1 <- mirt(dat, 1, itemtype = '2PL', constrain = list(c(1,5,9,13,17))))
# fitIndices(mod1)
#
# #Science data with computing the null model M2 stat
# (mod2 <- mirt(Science, 1))
# fitIndices(mod2, calcNull = TRUE)
# }
fitIndices <- function(obj, calcNull = FALSE, collapse_poly = TRUE, prompt = TRUE){
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
    nitems <- ncol(obj@data)
    if(any(is.na(obj@tabdata)))
        stop('M2 can not be calulated for data with missing values.')
    adj <- apply(obj@data, 2, min)
    dat <- t(t(obj@data) - adj)
    N <- nrow(dat)
#     if(!collapse_poly){
#         dat <- expand.table(tabdata)
#         browser()
#     }
    P  <- colMeans(dat)
    cross <- crossprod(dat, dat)
    P  <- c(P, cross[lower.tri(cross)]/N)    
    df <- length(P) - obj@nest
    if(df <= 0L)
        stop('Negative degrees of freedom')
    K <- obj@K
    pars <- obj@pars
    quadpts <- obj@quadpts    
    if(is.nan(quadpts)) 
        quadpts <- switch(as.character(obj@nfact), '1'=41, '2'=21, '3'=11, '4'=7, '5'=5, 3)
    npars <- 0
    for(i in 1L:(nitems+1L))
        npars <- npars + sum(pars[[i]]@est)
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
    E1 <- numeric(nitems)
    E2 <- numeric(nitems*(nitems - 1L)/2L)
    ind <- 1L
    for(i in 1L:nitems){
        x <- extract.item(obj, i)
        Ex <- expected.item(x, Theta, min=0L)
        E1[i] <- sum(Ex * Prior)
        for(j in 1L:nitems){
            if(i > j){
                y <- extract.item(obj, j)
                Ey <- expected.item(y, Theta, min=0L)
                E2[ind] <- sum(Ex * Ey * Prior)
                ind <- ind + 1L
            }
        }
    }
    E <- c(E1, E2)
    inv.Eta <- ginv(diag(E) - outer(E, E))
    #delta2.invEta.delta2 <- t(delta2) %*% inv.Eta %*% delta2
    #C2 <- inv.Eta - inv.Eta %*% delta2 %*% solve(delta2.invEta.delta2) %*% t(delta2) %*% inv.Eta
    C2 <- diag(length(P)) #this is just a placeholder
    M2 <- N * t(P - E) %*% C2 %*% (P - E) #the weight matrix C2 is missing...fack
    ret$M2 <- M2
    if(is.null(attr(obj, 'MG'))){
        ret$df.M2 <- df
        ret$p.M2 <- 1 - pchisq(M2, ret$df.M2)
        ret$RMSEA.M2 <- ifelse((M2 - ret$df.M2) > 0,
                        sqrt(M2 - ret$df.M2) / sqrt(ret$df.M2 * (N-1)), 0)
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
        ret$nrowT <- length(P)
    }
    return(as.data.frame(ret))
}
