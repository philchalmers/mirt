#' Compute the M2 model fit statistic
#'
#' Computes the M2 (Maydeu-Olivares & Joe, 2006) statistic when all data are dichotomous,
#' the collapsed M2* statistic (collapsing over univariate and bivariate response categories;
#' see Cai and Hansen, 2013), and the hybrid C2 statistic which only collapses only the bivariate
#' moments (Cai and Monro, 2014). The C2 variant is mainly useful when polytomous response models
#' do not have sufficient degrees of freedom to compute M2*. This function
#' also computes associated fit indices that are based on
#' fitting the null model. Supports single and multiple-group models.
#' If the latent trait density was approximated (e.g., Davidian curves, Empirical histograms, etc)
#' then passing \code{use_dentype_estimate = TRUE} will use the internally saved quadrature and
#' density components (where applicable).
#'
#' @return Returns a data.frame object with the M2-type statistic, along with the degrees of freedom,
#'   p-value, RMSEA (with 90\% confidence interval), SRMSR for each group,
#'   and optionally the TLI and CFI model fit statistics if \code{calcNull = TRUE}.
#'
#' @aliases M2
#' @param obj an estimated model object from the mirt package
#' @param type type of fit statistic to compute. Options are "M2", "M2*" for the univariate and
#'   bivariate collapsed version of the M2 statistic ("M2" currently limited to dichotomous
#'   response data only), and "C2" for a hybrid between
#'   M2 and M2* where only the bivariate moments are collapsed
#' @param quadpts number of quadrature points to use during estimation. If \code{NULL},
#'   a suitable value will be chosen based
#'   on the rubric found in \code{\link{fscores}}
#' @param calcNull logical; calculate statistics for the null model as well?
#'   Allows for statistics such as the limited information TLI and CFI. Only valid when items all
#'   have a suitable null model (e.g., those created via \code{\link{createItem}} will not)
#' @param theta_lim lower and upper range to evaluate latent trait integral for each dimension
# @param impute a number indicating how many imputations to perform
#   (passed to \code{\link{imputeMissing}}) when there are missing data present. This requires
#   a precomputed \code{Theta} input. Will return a data.frame object with the mean estimates
#   of the stats and their imputed standard deviations
#' @param CI numeric value from 0 to 1 indicating the range of the confidence interval for
#'   RMSEA. Default returns the 90\% interval
#' @param residmat logical; return the residual matrix used to compute the SRMSR statistic?
#'   Only the lower triangle of the residual correlation matrix will be returned
#'   (the upper triangle is filled with NA's)
#' @param QMC logical; use quasi-Monte Carlo integration? Useful for higher dimensional models.
#'   If \code{quadpts} not specified, 5000 nodes are used by default
#' @param suppress a numeric value indicating which parameter residual dependency combinations
#'   to flag as being too high. Absolute values for the standardized residuals greater than
#'   this value will be returned, while all values less than this value will be set to NA.
#'   Must be used in conjunction with the argument \code{residmat = TRUE}
#' @param ... additional arguments to pass
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Cai, L. & Hansen, M. (2013). Limited-information goodness-of-fit testing of
#' hierarchical item factor models. \emph{British Journal of Mathematical and Statistical
#' Psychology, 66}, 245-276.
#'
#' Cai, L. & Monro, S. (2014). \emph{A new statistic for evaluating item response theory
#' models for ordinal data}. National Center for Research on Evaluation, Standards,
#' & Student Testing. Technical Report.
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Maydeu-Olivares, A. & Joe, H. (2006). Limited information goodness-of-fit testing in
#' multidimensional contingency tables. \emph{Psychometrika, 71}, 713-732.
#' @keywords model fit
#' @export M2
#' @examples
#' \dontrun{
#' dat <- as.matrix(expand.table(LSAT7))
#' (mod1 <- mirt(dat, 1))
#' M2(mod1)
#' resids <- M2(mod1, residmat=TRUE) #lower triangle of residual correlation matrix
#' resids
#' summary(resids[lower.tri(resids)])
#'
#' # M2 with missing data present
#' dat[sample(1:prod(dim(dat)), 250)] <- NA
#' mod2 <- mirt(dat, 1)
#' M2(mod2)
#'
#' # C2 statistic (useful when polytomous IRT models have too few df)
#' pmod <- mirt(Science, 1)
#' # This fails with too few df:
#' # M2(pmod)
#' # This, however, works:
#' M2(pmod, type = 'C2')
#'
#' }
M2 <- function(obj, type="M2*", calcNull = TRUE, quadpts = NULL, theta_lim = c(-6, 6),
                CI = .9, residmat = FALSE, QMC = FALSE, suppress = 1, ...){
    impute <- 0
    if(is(obj, 'MixtureModel'))
        stop('Mixture IRT models not yet supported', call.=FALSE)
    fn <- function(Theta, obj, ...){
        dat <- imputeMissing(obj, Theta, warn=FALSE)
        tmpobj <- obj
        tmpobj@Data$data <- dat
        if(is(obj, 'MultipleGroupClass')){
            for(g in seq_len(length(obj@Data$groupNames)))
                tmpobj@ParObjects$pars[[g]]@Data$data <- dat[obj@Data$groupNames[g] == obj@Data$group,
                                                  , drop=FALSE]
        }
        return(M2(tmpobj, ...))
    }
    M2internal <- function(obj, calcNull, quadpts, theta_lim,
                           residmat, QMC, discrete, type,
                           use_dentype_estimate = FALSE, ...){
        ret <- list()
        group <- if(is.null(attr(obj, 'MG'))) 1 else attr(obj, 'MG')
        nitems <- ncol(obj@Data$data)
        # if(any(is.na(obj@Data$data)))
        #     stop('M2 can not be calculated for data with missing values.', call.=FALSE)
        adj <- obj@Data$mins
        dat <- t(t(obj@Data$data) - adj)
        N.full <- nrow(dat)
        N <- colSums(!is.na(dat))
        cN <- crossprod_miss(!is.na(dat), !is.na(dat))
        p  <- colMeans(dat, na.rm = TRUE)
        cross <- crossprod_miss(dat, dat)
        p <- c(p, cross[lower.tri(cross)]/cN[lower.tri(cross)])
        prodlist <- attr(obj@ParObjects$pars, 'prodlist')
        K <- obj@Data$K
        pars <- obj@ParObjects$pars
        if(is.null(quadpts))
            quadpts <- select_quadpts(obj@Model$nfact)
        if(obj@Model$nfact > 3L && !QMC && !discrete)
            warning('High-dimensional models should use quasi-Monte Carlo integration. Pass QMC=TRUE',
                    call.=FALSE)
        if(!discrete)
            gstructgrouppars <- ExtractGroupPars(pars[[nitems+1L]])
        estpars <- c()
        for(i in seq_len(nitems+1L)){
            if(i <= nitems)
                estpars <- c(estpars, pars[[i]]@est)
            else estpars <- c(estpars, rep(FALSE, length(pars[[i]]@est)))
        }
        bfactorlist <- obj@Internals$bfactor
        if(.hasSlot(obj@Model$lrPars, 'beta'))
            stop('Latent regression models not yet supported', call.=FALSE)
        # if(!discrete && obj@ParObjects$pars[[extract.mirt(obj, 'nitems')+1L]]@dentype == 'custom')
        #     stop('M2() does not currently support custom group densities', call.=FALSE)
        if(!discrete && !use_dentype_estimate){
            #         if(is.null(bfactorlist$Priorbetween[[1L]])){
            if(TRUE){ #TODO bifactor reduction possibility? Not as effective at computing marginals
                prior <- Priorbetween <- sitems <- specific <- NULL
                gstructgrouppars <- ExtractGroupPars(pars[[nitems+1L]])
                if(pars[[nitems+1L]]@dentype == 'custom'){
                    if(any(pars[[nitems+1L]]@est))
                        stop('custom group objects with estimated parameters not supported', call.=FALSE)
                    den_fun <- function(Theta, ...){
                        x <- obj@ParObjects$pars[[extract.mirt(obj, 'nitems')+1L]]
                        as.vector(x@den(x, Theta=Theta))
                    }
                    if(!is.null(obj@Internals$theta_lim))
                        theta_lim <- obj@Internals$theta_lim
                    theta <- as.matrix(seq(theta_lim[1L], theta_lim[2L], length.out = quadpts))
                    Theta <- thetaComb(theta, obj@Model$nfact)
                    Prior <- den_fun(Theta)
                    Prior <- Prior/sum(Prior)
                } else {
                    if(QMC){
                        Theta <- QMC_quad(npts=quadpts, nfact=obj@Model$nfact, lim=theta_lim)
                        Theta <- Theta_meanSigma_shift(Theta, gstructgrouppars$gmeans,
                                                       gstructgrouppars$gcov)
                        Prior <- rep(1/nrow(Theta), nrow(Theta))
                    } else {
                        theta <- as.matrix(seq(theta_lim[1L], theta_lim[2L], length.out = quadpts))
                        Theta <- thetaComb(theta, obj@Model$nfact)
                        Prior <- mirt_dmvnorm(Theta,gstructgrouppars$gmeans, gstructgrouppars$gcov)
                        Prior <- Prior/sum(Prior)
                    }
                }
                if(length(prodlist) > 0L)
                    Theta <- prodterms(Theta, prodlist)
            } else {
                Theta <- obj@Model$Theta
                prior <- bfactorlist$prior[[group]]; Priorbetween <- bfactorlist$Priorbetween[[group]]
                sitems <- bfactorlist$sitems; specific <- bfactorlist$specific;
                Prior <- bfactorlist$Prior[[group]]
            }
        } else {
            Theta <- obj@Model$Theta
            Prior <- obj@Internals$Prior[[1L]]
        }
        if(type == "M2*"){
            E1 <- E11 <- numeric(nitems)
            E2 <- matrix(NA, nitems, nitems)
            EIs <- EIs2 <- E11s <- matrix(0, nrow(Theta), nitems)
            DP <- matrix(0, nrow(Theta), length(estpars))
            wherepar <- c(1L, numeric(nitems))
            ind <- 1L
            for(i in seq_len(nitems)){
                x <- extract.item(obj, i)
                scs <- 1L:x@ncat - 1L
                Thetastar <- Theta
                if(x@nfixedeffects > 0)
                    Thetastar <- cbind(x@fixed.design[rep(1, nrow(Theta)), , drop=FALSE], Theta)
                EIs[,i] <- expected.item(x, Thetastar, min=0L)
                prob <- ProbTrace(x, Thetastar)
                E11s[,i] <- colSums((1L:ncol(prob)-1L)^2 * t(prob))
                cfs <- scs * scs
                EIs2[,i] <- t(cfs %*% t(prob))
                tmp <- length(x@parnum)
                DP[ ,ind:(ind+tmp-1L)] <- dP(x, Thetastar)
                ind <- ind + tmp
                wherepar[i+1L] <- ind
            }
            ind <- 1L
            for(i in seq_len(nitems)){
                E1[i] <- sum(EIs[,i] * Prior)
                E11[i] <- sum(E11s[,i] * Prior)
                for(j in seq_len(nitems)){
                    if(i >= j){
                        E2[i,j] <- sum(EIs[,i] * EIs[,j] * Prior)
                        ind <- ind + 1L
                    }
                }
            }
            e <- c(E1, E2[lower.tri(E2)])
            E2[is.na(E2)] <- 0
            E2 <- E2 + t(E2)
            diag(E2) <- E11
            R <- cov2cor(cross/cN - outer(colMeans(dat, na.rm=TRUE), colMeans(dat, na.rm=TRUE)))
            Kr <- cov2cor(E2 - outer(E1, E1))
            SRMSR <- sqrt( sum((R[lower.tri(R)] - Kr[lower.tri(Kr)])^2) / sum(lower.tri(R)))
            if(residmat){
                ret <- matrix(NA, nrow(R), nrow(R))
                ret[lower.tri(ret)] <- R[lower.tri(R)] - Kr[lower.tri(Kr)]
                colnames(ret) <- rownames(ret) <- colnames(obj@Data$dat)
                if(suppress < 1)
                    ret[lower.tri(ret)][abs(ret[lower.tri(ret)]) < suppress] <- NA
                return(ret)
            }
            delta1 <- matrix(0, nitems, length(estpars))
            delta2 <- matrix(0, length(p) - nitems, length(estpars))
            ind <- 1L
            offset <- pars[[1L]]@parnum[1L] - 1L
            for(i in seq_len(nitems)){
                dp <- colSums(DP[ , wherepar[i]:(wherepar[i+1L]-1L), drop=FALSE] * Prior)
                delta1[i, pars[[i]]@parnum - offset] <- dp
                for(j in seq_len(nitems)){
                    if(i < j){
                        dp <- colSums(DP[ , wherepar[i]:(wherepar[i+1L]-1L), drop=FALSE] * EIs[,j] * Prior)
                        delta2[ind, pars[[i]]@parnum - offset] <- dp
                        dp <- colSums(DP[ , wherepar[j]:(wherepar[j+1L]-1L), drop=FALSE] * EIs[,i] * Prior)
                        delta2[ind, pars[[j]]@parnum - offset] <- dp
                        ind <- ind + 1L
                    }
                }
            }
        } else if(type == "C2"){
            nK <- sum(K-1L)
            PIs <- matrix(0, nrow(Theta), nK)
            E1 <- E11 <- numeric(nitems)
            E2 <- matrix(NA, nitems, nitems)
            EIs <- EIs2 <- E11s <- matrix(0, nrow(Theta), nitems)
            DP <- matrix(0, nrow(Theta), length(estpars))
            DPIs <- vector('list', nitems)
            wherepar <- c(1L, numeric(nitems))
            ind <- pind <- 1L
            for(i in seq_len(nitems)){
                x <- extract.item(obj, i)
                scs <- 1L:x@ncat - 1L
                Thetastar <- Theta
                if(x@nfixedeffects > 0)
                    Thetastar <- cbind(x@fixed.design[rep(1, nrow(Theta)), , drop=FALSE], Theta)
                EIs[,i] <- expected.item(x, Thetastar, min=0L)
                prob <- ProbTrace(x, Thetastar)
                PIs[,pind:(pind+ncol(prob)-2L)] <- prob[,-1L]
                E11s[,i] <- colSums((1L:ncol(prob)-1L)^2 * t(prob))
                cfs <- scs * scs
                EIs2[,i] <- t(cfs %*% t(prob))
                tmp <- length(x@parnum)
                DP[ ,ind:(ind+tmp-1L)] <- dP(x, Thetastar)
                pind <- pind + K[i] - 1L
                ind <- ind + tmp
                wherepar[i+1L] <- ind
            }
            ind <- 1L
            for(i in seq_len(nitems)){
                E1[i] <- sum(EIs[,i] * Prior)
                E11[i] <- sum(E11s[,i] * Prior)
                for(j in seq_len(nitems)){
                    if(i >= j){
                        E2[i,j] <- sum(EIs[,i] * EIs[,j] * Prior)
                        ind <- ind + 1L
                    }
                }
            }
            P1 <- colSums(PIs * Prior)
            e <- c(P1, E2[lower.tri(E2)])
            if(all(sapply(obj@ParObjects$pars, class) %in% c(ordinal_itemtypes(), 'GroupPars'))){
                E2[is.na(E2)] <- 0
                E2 <- E2 + t(E2)
                diag(E2) <- E11
                R <- cov2cor(cross/cN - outer(colMeans(dat, na.rm=TRUE), colMeans(dat, na.rm=TRUE)))
                Kr <- cov2cor(E2 - outer(E1, E1))
                SRMSR <- sqrt( sum((R[lower.tri(R)] - Kr[lower.tri(Kr)])^2) / sum(lower.tri(R)))
            } else SRMSR <- NULL
            delta1 <- matrix(0, nK, length(estpars))
            delta2 <- matrix(0, length(p) - nitems, length(estpars))
            ind <- pind1 <- pind2 <- 1L
            offset <- pars[[1L]]@parnum[1L] - 1L
            for(i in seq_len(nitems)){
                x <- extract.item(obj, i)
                Thetastar <- Theta
                if(x@nfixedeffects > 0)
                    Thetastar <- cbind(x@fixed.design[rep(1, nrow(Theta)), , drop=FALSE], Theta)
                tmp <- lapply(numDeriv_dP2(x, Thetastar), function(x) colSums(x * Prior))
                dp <- if(length(tmp) == 1L) matrix(tmp[[1L]], nrow=1L) else do.call(rbind, tmp)
                delta1[pind1:(pind1+nrow(dp)-1L), pars[[i]]@parnum - offset] <- dp
                pind1 <- pind1 + nrow(dp)
                pind2 <- pind2 + ncol(dp)
                for(j in seq_len(nitems)){
                    if(i < j){
                        dp <- colSums(DP[ , wherepar[i]:(wherepar[i+1L]-1L), drop=FALSE] *
                                          EIs[,j] * Prior)
                        delta2[ind, pars[[i]]@parnum - offset] <- dp
                        dp <- colSums(DP[ , wherepar[j]:(wherepar[j+1L]-1L), drop=FALSE] *
                                          EIs[,i] * Prior)
                        delta2[ind, pars[[j]]@parnum - offset] <- dp
                        ind <- ind + 1L
                    }
                }
            }
            itemloc <- obj@Model$itemloc
            itemloc <- itemloc[-length(itemloc)]
            was_na <- is.na(extract.mirt(obj, 'data'))
            fulldata <- obj@Data$fulldata[[1L]]
            for(i in 1:(nitems)){
                pick <- if(i == nitems) itemloc[i]:ncol(fulldata)
                    else itemloc[i]:(itemloc[i+1]-1)
                fulldata[was_na[,i], pick] <- NA
            }
            N <- colSums(!is.na(fulldata))[-itemloc]
            p <- c(colMeans(fulldata[,-itemloc], na.rm=TRUE),
                   cross[lower.tri(cross)]/cN[lower.tri(cross)])
        } else {
            # M2 TODO
        }
        delta <- rbind(delta1, delta2)
        abcats <- do.call(c, lapply(1:nitems, function(x) rep(x, each=K[x]-1))) - 1L
        abcats2 <- do.call(c, lapply(K, function(x) 2L:x - 1L))
        Xi2els <- if(type == "M2*"){
            .Call('buildXi2els', nrow(delta1), nrow(delta2), nitems, EIs, EIs2, Prior)
        } else .Call('buildXi2els_C2', nrow(delta1), nrow(delta2), ncol(PIs), nitems,
                     PIs, EIs, EIs2, Prior, abcats, abcats2)
        Xi2 <- rbind(cbind(Xi2els$Xi11, Xi2els$Xi12), cbind(t(Xi2els$Xi12), Xi2els$Xi22))
        Nstar <- c(N, cN[lower.tri(cN)])
        ret <- list(Xi2=Xi2, delta=delta, estpars=estpars,
                    p=sqrt(Nstar)*p, e=sqrt(Nstar)*e, SRMSR=SRMSR, N=Nstar,
                    N.ratio=Nstar/N.full)
        ret
    }

    #main
    if(residmat) type <- "M2*"
    stopifnot(type %in% c('M2*', 'M2', 'C2'))
    if(all(extract.mirt(obj, 'K') == 2)) type <- 'M2*'
    if(type == "M2"){
        if(!all(extract.mirt(obj, 'K') == 2L))
            stop("M2 statistic currently not supported for polytomous data. Use M2* or C2 instead",
                 call.=FALSE)
    }
    if(missing(obj)) missingMsg('obj')
    if(is(obj, 'MixedClass'))
        stop('MixedClass objects are not yet supported', call.=FALSE)
    if(is(obj, 'MixtureClass'))
        stop('MixtureClass objects are not yet supported', call.=FALSE)
    if(QMC && is.null(quadpts)) quadpts <- 5000L
    if(nrow(extract.mirt(obj, 'data')) == 0L) stop('No data!', call.=FALSE)
    discrete <- FALSE
    if(is(obj, 'DiscreteClass')){
        discrete <- TRUE
        class(obj) <- 'MultipleGroupClass'
    }
    alpha <- (1 - CI)/2
    pars <- obj@ParObjects$pars
    ngroups <- extract.mirt(obj, 'ngroups')
    ret <- vector('list', ngroups)
    for(g in seq_len(ngroups)){
        if(ngroups > 1L || discrete){
            attr(pars[[g]], 'MG') <- g
            pars[[g]]@Internals$bfactor <- obj@Internals$bfactor
            if(discrete){
                pars[[g]]@Internals$Prior <- list(obj@Internals$Prior[[g]])
                pars[[g]]@Model$Theta <- obj@Model$Theta
            }
            pars[[g]]@Data <- list(data=obj@Data$data[obj@Data$group == obj@Data$groupNames[g], ],
                                   mins=obj@Data$mins, K=obj@Data$K,
                                   fulldata=list(obj@Data$fulldata[[g]]))
            if(is(obj, 'MixtureClass')) pars[[g]]@Data$data <- extract.mirt(obj, 'data')
            ret[[g]] <- M2internal(pars[[g]], calcNull=FALSE, quadpts=quadpts, theta_lim=theta_lim,
                                   residmat=residmat, QMC=QMC, discrete=discrete, type=type, ...)
        } else {
            ret[[g]] <- M2internal(obj, calcNull=FALSE, quadpts=quadpts, theta_lim=theta_lim,
                                   residmat=residmat, QMC=QMC, discrete=discrete, type=type, ...)
        }
        # build
        if(!residmat){
            if(g == 1L){
                estpars <- ret[[1L]]$estpars
                delta <- ret[[1L]]$delta
                Xi2 <- ret[[1L]]$Xi2
                dims_delta <- dim(delta)
                dims_Xi2 <- dim(Xi2)
                p <- ret[[1L]]$p
                e <- ret[[1L]]$e
                if(ngroups > 1L){
                    tmp <- matrix(0, dims_delta[1L]*ngroups, dims_delta[2L]*ngroups)
                    tmp[1L:dims_delta[1L], 1L:dims_delta[2L]] <- delta
                    delta <- tmp
                    tmp <- matrix(0, dims_Xi2[1L]*ngroups, dims_Xi2[2L]*ngroups)
                    tmp[1L:dims_Xi2[1L], 1L:dims_Xi2[2L]] <- Xi2
                    Xi2 <- tmp
                }
            } else {
                estpars <- c(estpars, ret[[g]]$estpars)
                p <- c(p, ret[[g]]$p)
                e <- c(e, ret[[g]]$e)
                delta[(1L + (g-1L)*dims_delta[1L]):(g*dims_delta[1L]),
                      (1L + (g-1L)*dims_delta[2L]):(g*dims_delta[2L])] <- ret[[g]]$delta
                Xi2[(1L + (g-1L)*dims_Xi2[1L]):(g*dims_Xi2[1L]),
                      (1L + (g-1L)*dims_Xi2[2L]):(g*dims_Xi2[2L])] <- ret[[g]]$Xi2
            }
        }
    }
    if(residmat){
        if(ngroups == 1L) return(ret[[1L]])
        names(ret) <- obj@Data$groupNames
        return(ret)
    }
    constrain <- extract.mirt(obj, 'constrain')
    if(length(constrain)){
        for(i in seq_len(length(constrain))){
            delta[ ,constrain[[i]][1L]] <- rowSums(delta[ ,constrain[[i]]])
            estpars[constrain[[i]][-1L]] <- FALSE
        }
    }
    delta <- delta[ ,estpars, drop=FALSE]
    tmp <- qr.Q(qr(delta), complete=TRUE)
    if((ncol(delta) + 1L) > ncol(tmp))
        stop('M2() statistic cannot be calculated due to too few degrees of freedom',
             call.=FALSE)
    deltac <- tmp[,(ncol(delta) + 1L):ncol(tmp), drop=FALSE]
    N <- nrow(extract.mirt(obj, 'data'))
    C2 <- try(deltac %*% solve(t(deltac) %*% Xi2 %*% deltac) %*% t(deltac), TRUE)
    if(is(C2, 'try-error'))
        stop('Could not invert orthogonal complement matrix', call.=FALSE)
    Ns.ratio <- do.call(c, lapply(ret, function(x) x$N.ratio))
    C2 <- outer(sqrt(Ns.ratio), sqrt(Ns.ratio)) * C2
    M2 <- abs(t(p - e) %*% C2 %*% (p - e))
    df <- length(p) - extract.mirt(obj, 'nest')
    # df <- qr(deltac)$rank
    newret <- list(M2=M2, df=df)
    newret$p <- 1 - pchisq(M2, df)
    newret$RMSEA <- rmsea(X2=M2, df=df, N=N)   # CHECKME this seems off with missing data
    RMSEA.90_CI <- RMSEA.CI(M2, df, N, ci.lower=alpha, ci.upper=1-alpha)
    newret[[paste0("RMSEA_", alpha*100)]]  <- RMSEA.90_CI[1L]
    newret[[paste0("RMSEA_", (1-alpha)*100)]] <- RMSEA.90_CI[2L]
    if(!is.null(ret[[1L]]$SRMSR)){
        SRMSR <- numeric(ngroups)
        for(g in seq_len(ngroups))
            SRMSR[g] <- ret[[g]]$SRMSR
        if(ngroups > 1){
            names(SRMSR) <- obj@Data$groupNames
            SRMSR <- as.list(SRMSR)
        }
        newret$SRMSR <- SRMSR
    }
    if(calcNull){
        null.mod <- try(computeNullModel(data=obj@Data$data, group=obj@Data$group,
                                         key=obj@Internals$key))
        if(is(null.mod, 'try-error'))
            stop('Null model did not converge or is not supported', call.=FALSE)
        null.fit <- M2(null.mod, calcNull=FALSE, type=type, quadpts=2)
        if(null.fit$M2 > newret$M2){
            newret$TLI <- tli(X2=newret$M2, X2.null=null.fit$M2, df=newret$df,
                              df.null=null.fit$df)
            newret$CFI <- cfi(X2=newret$M2, X2.null=null.fit$M2, df=newret$df,
                              df.null=null.fit$df)
        } else warning('Null model chi-squared value smaller than fitted model', call.=FALSE)
    }
    newret <- as.data.frame(newret)
    rownames(newret) <- 'stats'
    newret
}
