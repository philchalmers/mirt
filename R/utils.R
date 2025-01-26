#' Create all possible combinations of vector input
#'
#' This function constructs all possible k-way combinations of an input vector.
#' It is primarily useful when used in conjunction with the \code{\link{mdirt}} function,
#' though users may have other uses for it as well. See \code{\link{expand.grid}} for more
#' flexible combination formats.
#'
#' @param theta the vector from which all possible combinations should be obtained
#' @param nfact the number of observations (and therefore the number of columns to return in
#'   the matrix of combinations)
#' @param intercept logical; should a vector of 1's be appended to the first column of the
#'   result to include an intercept design component? Default is \code{FALSE}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @return a matrix with all possible combinations
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @export
#' @examples
#'
#' # all possible joint combinations for the vector -4 to 4
#' thetaComb(-4:4, 2)
#'
#' # all possible binary combinations for four observations
#' thetaComb(c(0,1), 4)
#'
#' # all possible binary combinations for four observations (with intercept)
#' thetaComb(c(0,1), 4, intercept=TRUE)
#'
thetaComb <- function(theta, nfact, intercept = FALSE)
{
	if (nfact == 1L){
        Theta <- matrix(theta)
	} else {
        thetalist <- vector('list', nfact)
        for(i in seq_len(nfact))
            thetalist[[i]] <- theta
        Theta <- as.matrix(expand.grid(thetalist))
	}
    if(intercept) Theta <- cbind(1, Theta)
    colnames(Theta) <- NULL
	return(Theta)
}

thetaStack <- function(theta, nclass){
    thetalist <- vector('list', nclass)
    for(i in seq_len(nclass))
        thetalist[[i]] <- theta
    as.matrix(do.call(rbind, thetalist))
}

#' Second-order test of convergence
#'
#' Test whether terminated estimation criteria for a given model passes
#' the second order test by checking the positive definiteness of the resulting
#' Hessian matrix. This function, which accepts the symmetric Hessian/information
#' matrix as the input, returns \code{TRUE} if the matrix is positive definite
#' and \code{FALSE} otherwise.
#'
#' @param mat symmetric matrix to test for positive definiteness (typically the Hessian at the
#'   highest point of model estimator, such as MLE or MAP)
#' @param ... arguments passed to either \code{\link{eigen}}, \code{\link{chol}}, or
#'   \code{'det'} for the positiveness of the eigen values, positiveness of leading minors
#'   via the Cholesky decomposition, or evaluation of whether the determinant
#'   is greater than 0
#' @param method method to use to test positive definiteness. Default is \code{'eigen'}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @return a matrix with all possible combinations
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @export
#' @examples
#'
#' \dontrun{
#'
#' # PD matrix
#' mod <- mirt(Science, 1, SE=TRUE)
#' info <- solve(vcov(mod))   ## observed information
#' secondOrderTest(info)
#' secondOrderTest(info, method = 'chol')
#' secondOrderTest(info, method = 'det')
#'
#' # non-PD matrix
#' mat <- matrix(c(1,0,0,0,1,1,0,1,1), ncol=3)
#' mat
#' secondOrderTest(mat)
#' secondOrderTest(mat, method = 'chol')
#' secondOrderTest(mat, method = 'det')
#'
#' }
secondOrderTest <- function(mat, ..., method = 'eigen'){
    if(method == 'eigen'){
        evs <- eigen(mat, ...)$values
        ret <- all(!sapply(evs, function(x) isTRUE(all.equal(x, 0))) & evs > 0)
    } else if(method == 'chol'){
        chl <- try(chol(mat, ...), silent = TRUE)
        ret <- if(is(chl, "try-error")) FALSE else TRUE
    } else if(method == 'det'){
        dt <- det(mat, ...)
        ret <- !isTRUE(all.equal(dt, 0)) && dt > 0
    }
    ret
}

# Product terms
prodterms <- function(theta0, prodlist)
{
    products <- matrix(1, ncol = length(prodlist), nrow = nrow(theta0))
    for(i in seq_len(length(prodlist))){
        tmp <- prodlist[[i]]
        for(j in 1L:length(tmp))
            products[ ,i] <- products[ ,i] * theta0[ ,tmp[j]]
    }
    ret <- cbind(theta0,products)
    ret
}

gain_fun <- function(gain, t) (gain[1L] / t)^gain[2L]

update_cand.var <- function(PAs, CTVs, target = .5){
    pick <- max(which(!is.na(PAs)))
    if(PAs[pick] < .4 & PAs[pick] > .3) return(CTVs[pick])
    if(PAs[pick] < .01) return(min(CTVs, na.rm = TRUE)/10)
    if(PAs[pick] > .99) return(max(CTVs, na.rm = TRUE)*5)
    if(length(unique(CTVs)) < 4L){
        return(CTVs[pick] * ifelse(PAs[pick] > .5, 1.1, .8))
    }
    df <- data.frame(PAs, CTVs)
    mod <- lm(plogis(PAs) ~ CTVs, data=df)
    fn <- function(x) (target - qlogis(predict(mod, data.frame(CTVs=x))))^2
    root <- nlm(f = fn, p = CTVs[pick])$estimate
    if(root < 0) root <- min(CTVs, na.rm = TRUE) / 2
    return(root)
}

# MH sampler for theta values
draw.thetas <- function(theta0, pars, fulldata, itemloc, cand.t.var, prior.t.var,
                        prior.mu, prodlist, OffTerm, CUSTOM.IND)
{
    makedraws <- try({
        N <- nrow(fulldata)
        unif <- runif(N)
        sigma <- if(ncol(theta0) == 1L) matrix(cand.t.var) else diag(rep(cand.t.var,ncol(theta0)))
        total_0 <- attr(theta0, 'log.lik_full')
        theta1prod <- theta1 <- theta0 + mirt_rmvnorm(N, sigma = sigma)
        if(is.null(total_0)) theta1prod <- theta1 <- theta0 #for intial draw
        if(length(prodlist))
            theta1prod <- prodterms(theta1, prodlist)
        total_1 <- complete.LL(theta=theta1prod, pars=pars, nfact=ncol(theta1), prior.mu=prior.mu,
                               prior.t.var=prior.t.var, OffTerm=OffTerm,
                               CUSTOM.IND=CUSTOM.IND, itemloc=itemloc, fulldata=fulldata)
        if(is.null(total_0)){ #for intial draw
            attr(theta1, 'log.lik_full') <- total_1
            return(theta1)
        }
        diff <- total_1 - total_0
        accept <- log(unif) < diff
        theta1[!accept, ] <- theta0[!accept, ]
        total_1[!accept] <- total_0[!accept]
        log.lik <- sum(total_1)
        attr(theta1, "Proportion Accepted") <- sum(accept)/N
        attr(theta1, "log.lik") <- log.lik
        attr(theta1, 'log.lik_full') <- total_1
    }, silent = TRUE)
    if(is(makedraws, 'try-error'))
        stop('MH sampler failed. Model is likely unstable or may need better starting values',
             .call=FALSE)
    return(theta1)
}

complete.LL <- function(theta, pars, nfact, prior.mu, prior.t.var,
                        OffTerm, CUSTOM.IND, itemloc, fulldata){
    log_den <- mirt_dmvnorm(theta[,seq_len(nfact), drop=FALSE], prior.mu, prior.t.var, log=TRUE)
    itemtrace <- computeItemtrace(pars=pars, Theta=theta, itemloc=itemloc,
                                   offterm=OffTerm, CUSTOM.IND=CUSTOM.IND)
    rowSums(fulldata * log(itemtrace)) + log_den
}

imputePars <- function(pars, imputenums, constrain, pre.ev){
    shift <- mirt_rmvnorm(1L, mean=numeric(length(pre.ev$values)), pre.ev=pre.ev)
    for(i in seq_len(length(pars))){
        pn <- pars[[i]]@parnum
        pick2 <- imputenums %in% pn
        pick1 <- pn %in% imputenums
        pars[[i]]@par[pick1] <- pars[[i]]@par[pick1] + shift[pick2]
        if(is(pars[[i]], 'graded')){
            where <- (length(pars[[i]]@par) - pars[[i]]@ncat + 2L):length(pars[[i]]@par)
            if(!(all(sort(pars[[i]]@par[where], decreasing=TRUE) == pars[[i]]@par[where])))
                stop('Drawn values out of order', call.=FALSE)
        } else if(is(pars[[i]], 'grsm')){
            where <- (length(pars[[i]]@par) - pars[[i]]@ncat + 1L):(length(pars[[i]]@par)-1L)
            if(!(all(sort(pars[[i]]@par[where], decreasing=TRUE) == pars[[i]]@par[where])))
                stop('Drawn values out of order', call.=FALSE)
        }
    }
    for(con in seq_len(length(constrain))){
        tmp <- shift[imputenums %in% constrain[[con]][1L]]
        if(length(tmp)){
            for(i in seq_len(length(pars))){
                pick <- pars[[i]]@parnum %in% constrain[[con]][-1L]
                pars[[i]]@par[pick] <- tmp + pars[[i]]@par[pick]
            }
        }
    }
    return(pars)
}

imputePars2 <- function(MGmod, shortpars, longpars, imputenums, pre.ev){
    while(TRUE){
        shift <- mirt_rmvnorm(1L, mean=shortpars, pre.ev=pre.ev)
        longpars[imputenums] <- shift[1L,]
        constrain <- MGmod@Model$constrain
        longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
        pars <- list(MGmod@ParObjects$pars[[1L]]@ParObjects$pars, MGmod@ParObjects$pars[[2L]]@ParObjects$pars)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=2L, J=length(pars[[1L]])-1L)
        if(any(MGmod@Model$itemtype %in% c('graded', 'grsm'))){
            pick <- c(MGmod@Model$itemtype %in% c('graded', 'grsm'), FALSE)
            if(!all(sapply(pars[[1L]][pick], CheckIntercepts) &
                    sapply(pars[[2L]][pick], CheckIntercepts))) next
        }
        break
    }
    MGmod@ParObjects$pars[[1L]]@ParObjects$pars <- pars[[1L]]
    MGmod@ParObjects$pars[[2L]]@ParObjects$pars <- pars[[2L]]
    MGmod
}

# Rotation function
Rotate <- function(F, rotate, Target = NULL, par.strip.text = NULL, par.settings = NULL, digits, ...)
{
    if(ncol(F) == 1L) rotF <- list()
    if(rotate == 'none') rotF <- list(loadings=F, Phi=diag(ncol(F)), orthogonal=TRUE)
	if(rotate == 'promax'){
        mypromax <- function (x, m = 4) {
                #borrowed and modified from stats::promax on Febuary 13, 2013
                if (ncol(x) < 2L)
                    return(x)
                dn <- dimnames(x)
                xx <- varimax(x)
                x <- xx$loadings
                Q <- x * abs(x)^(m - 1)
                U <- lm.fit(x, Q)$coefficients
                d <- diag(solve(t(U) %*% U))
                U <- U %*% diag(sqrt(d))
                dimnames(U) <- NULL
                z <- x %*% U
                U <- xx$rotmat %*% U
                ui <- solve(U)
                Phi <- ui %*% t(ui)
                dimnames(z) <- dn
                class(z) <- "loadings"
                list(loadings = z, rotmat = U, Phi = Phi, orthogonal = FALSE)
            }
        rotF <- mypromax(F, ...)
	}
    if(rotate == 'oblimin') rotF <- GPArotation::oblimin(F, ...)
	if(rotate == 'quartimin') rotF <- GPArotation::quartimin(F, ...)
	if(rotate == 'targetT') rotF <- GPArotation::targetT(F, Target = Target, ...)
	if(rotate == 'targetQ') rotF <- GPArotation::targetQ(F, Target = Target, ...)
	if(rotate == 'pstT') rotF <- GPArotation::pstT(F, Target = Target, ...)
	if(rotate == 'pstQ') rotF <- GPArotation::pstQ(F, Target = Target, ...)
	if(rotate == 'oblimax') rotF <- GPArotation::oblimax(F, ...)
	if(rotate == 'entropy') rotF <- GPArotation::entropy(F, ...)
	if(rotate == 'quartimax') rotF <- GPArotation::quartimax(F, ...)
	if(rotate == 'varimax') rotF <- GPArotation::Varimax(F, ...)
	if(rotate == 'simplimax') rotF <- GPArotation::simplimax(F, ...)
	if(rotate == 'bentlerT') rotF <- GPArotation::bentlerT(F, ...)
	if(rotate == 'bentlerQ') rotF <- GPArotation::bentlerQ(F, ...)
	if(rotate == 'tandemI') rotF <- GPArotation::tandemI(F, ...)
	if(rotate == 'tandemII') rotF <- GPArotation::tandemII(F, ...)
	if(rotate == 'geominT') rotF <- GPArotation::geominT(F, ...)
	if(rotate == 'geominQ') rotF <- GPArotation::geominQ(F, ...)
	if(rotate == 'cfT') rotF <- GPArotation::cfT(F, ...)
	if(rotate == 'cfQ') rotF <- GPArotation::cfQ(F, ...)
	if(rotate == 'infomaxT') rotF <- GPArotation::infomaxT(F, ...)
	if(rotate == 'infomaxQ') rotF <- GPArotation::infomaxQ(F, ...)
	if(rotate == 'mccammon') rotF <- GPArotation::mccammon(F, ...)
	if(rotate == 'bifactorT') rotF <- GPArotation::bifactorT(F, ...)
	if(rotate == 'bifactorQ') rotF <- GPArotation::bifactorQ(F, ...)
	s <- apply(rotF$loadings, 2L, function(x)
	    sign(x[which.max(abs(x))]))
	rotF$loadings <- t(s * t(rotF$loadings))
	if(is.null(rotF$Phi)) rotF$Phi <- diag(ncol(rotF$loadings))
	rotF$Phi <- diag(s) %*% rotF$Phi %*% diag(s)
	return(unclass(rotF))
}

# Gamma correlation, mainly for obtaining a sign
gamma_cor <- function(x)
{
	concordant <- function(x){
			mat.lr <- function(r, c, r.x, c.x){
				lr <- as.numeric(x[(r.x > r) & (c.x > c)])
				sum(lr)
			}
		r.x <- row(x)
		c.x <- col(x)
		sum(x * mapply(mat.lr, r = r.x, c = c.x, MoreArgs=list(r.x=r.x, c.x=c.x)))
	}
	discordant <- function(x){
		mat.ll <- function(r, c, r.x, c.x){
			ll <- as.numeric(x[(r.x > r) & (c.x < c)])
			sum(ll)
		}
		r.x <- row(x)
		c.x <- col(x)
		sum(x * mapply(mat.ll, r = r.x, c = c.x, MoreArgs=list(r.x=r.x, c.x=c.x)))
	}
	c <- concordant(x)
	d <- discordant(x)
	gamma <- (c - d) / (c + d)
	gamma
}

# Approximation to polychoric matrix for initial values
cormod <- function(fulldata, K, guess, smooth = TRUE, use = 'pairwise.complete.obs')
{
	fulldata <- as.matrix(fulldata)
	cormat <- suppressWarnings(cor(fulldata, use=use))
    diag(cormat) <- 1
    cormat[is.na(cormat)] <- 0
	cormat <- abs(cormat)^(1/1.15) * sign(cormat)
	if(smooth)
		cormat <- smooth.cor(cormat)
	cormat
}

# Rotate lambda coefficients
rotateLambdas <- function(so){
    F <- so$rotF
    h2 <- so$h2
    h <- matrix(rep(sqrt(1 - h2), ncol(F)), ncol = ncol(F))
    a <- F / h
    a
}

d2r <-function(d) pi*d/180

closeEnough <- function(x, low, up) all(x >= low & x <= up)

logit <- function(x){
    ret <- qlogis(x)
    ret <- ifelse(x == 0, -999, ret)
    ret <- ifelse(x == 1, 999, ret)
    ret
}

antilogit <- function(x) plogis(x)

MPinv <- function(mat){
    svd <- svd(mat)
    d <- svd$d; v <- svd$v; u <- svd$u
    P <- d > max(sqrt(.Machine$double.eps) * d[1L], 0)
    if(all(P)){
        mat <- v %*% (1/d * t(u))
    } else {
        mat <- v[ , P, drop=FALSE] %*% ((1/d[P]) * t(u[ , P, drop=FALSE]))
    }
    mat
}

test_info <- function(pars, Theta, Alist, K){
    infolist <- list()
    for(cut in seq_len(length(Alist))){
        A <- Alist[[cut]]
        info <- rep(0,nrow(Theta))
        for(j in seq_len(length(K))){
            info <- info + ItemInfo(pars[[j]], A[j,], Theta)
        }
        infolist[[cut]] <- info
    }
    tmp <- 0
    for(i in seq_len(length(infolist))){
        tmp <- tmp + infolist[[i]]
    }
    info <- tmp/length(infolist)
    info
}

Lambdas <- function(pars, Names){
    J <- length(pars) - 1L
    lambdas <- matrix(NA, J, length(ExtractLambdas(pars[[1L]], include_fixed = FALSE)))
    gcov <- ExtractGroupPars(pars[[J+1L]])$gcov
    if(ncol(gcov) < ncol(lambdas)){
        tmpcov <- diag(ncol(lambdas))
        tmpcov[seq_len(ncol(gcov)), seq_len(ncol(gcov))] <- gcov
        gcov <- tmpcov
    }
    rownames(lambdas) <- Names
    for(i in seq_len(J)){
        tmp <- pars[[i]]
        lambdas[i,] <- ExtractLambdas(tmp, include_fixed = FALSE) /1.702
        if(tmp@ncat > 2 && (is(tmp, 'gpcm') || is(tmp, 'nominal'))){
            aks <- tmp@par[ncol(lambdas) + 1:tmp@ncat]
            d <- nominal_rescale_d(tmp@par[1:(ncol(lambdas) + tmp@ncat)],
                              nfact=ncol(lambdas), ncat=tmp@ncat)
            lambdas[i,] <- lambdas[i,] * d
            lambdas[i,] <- lambdas[i,] * sqrt((tmp@ncat-1)* 1.702)*sign(lambdas[i, 1])
        }
    }
    dcov <- if(ncol(gcov) > 1L) diag(sqrt(diag(gcov))) else matrix(sqrt(diag(gcov)))
    lambdas <- lambdas %*% dcov
    norm <- sqrt(1 + rowSums(lambdas^2))
    F <- as.matrix(lambdas/norm)
    F
}

#change long pars for groups into mean in sigma
ExtractGroupPars <- function(x){
    if(x@itemclass < 0L) return(list(gmeans=0, gcov=matrix(1)))
    nfact <- x@nfact
    gmeans <- x@par[seq_len(nfact)]
    phi_matches <- grepl("PHI", x@parnames)
    if (x@dentype == "Davidian") {
        phi <- x@par[phi_matches]
        tmp <- x@par[-c(seq_len(nfact), which(phi_matches))]
        gcov <- matrix(0, nfact, nfact)
        gcov[lower.tri(gcov, diag=TRUE)] <- tmp
        gcov <- makeSymMat(gcov)
        return(list(gmeans=gmeans, gcov=gcov, phi=phi))
    } else {
        par <- x@par
        if(x@dentype == "mixture") par <- par[-length(par)] # drop pi
        tmp <- par[-seq_len(nfact)]
        gcov <- matrix(0, nfact, nfact)
        gcov[lower.tri(gcov, diag=TRUE)] <- tmp
        gcov <- makeSymMat(gcov)
        return(list(gmeans=gmeans, gcov=gcov))
    }
}

ExtractMixtures <- function(pars){
    pick <- length(pars[[1L]])
    logit_pi <- sapply(pars, function(x) x[[pick]]@par[length(x[[pick]]@par)])
    psumexp(logit_pi)
}

psumexp <- function(logit){
    max_logit <- max(logit)
    pi <- exp(logit - max_logit)
    pi / sum(pi)
}

reloadConstr <- function(par, constr, obj){
    par2 <- rep(NA, length(constr[[1L]]))
    notconstr <- rep(TRUE, length(par2))
    for(i in seq_len(length(constr))){
        par2[constr[[i]]] <- par[i]
        notconstr[constr[[i]]] <- FALSE
    }
    par2[notconstr] <- par[(length(constr)+1L):length(par)]
    ind <- 1L
    for(i in seq_len(length(obj))){
        obj[[i]]@par[obj[[i]]@est] <- par2[ind:(ind + sum(obj[[i]]@est) - 1L)]
        ind <- ind + sum(obj[[i]]@est)
    }
    return(obj)
}

bfactor2mod <- function(model, J){
    tmp <- tempfile('tempfile')
    unique <- sort(unique(model))
    index <- seq_len(J)
    tmp2 <- c()
    for(i in seq_len(length(unique))){
        ind <- na.omit(index[model == unique[i]])
        comma <- rep(',', 2*length(ind))
        TF <- rep(c(TRUE,FALSE), length(ind))
        comma[TF] <- ind
        comma[length(comma)] <- ""
        tmp2 <- c(tmp2, c(paste('\nS', i, ' =', sep=''), comma))
    }
    cat(tmp2, file=tmp)
    model <- mirt.model(file=tmp, quiet = TRUE)
    unlink(tmp)
    return(model)
}

Theta_meanSigma_shift <- function(Theta, mean, sigma){
    t(t(Theta) + mean) %*% chol(sigma)
}

updateTheta <- function(npts, nfact, pars, QMC = FALSE){
    ngroups <- length(pars)
    pick <- length(pars[[1L]])
    Theta <- vector('list', ngroups)
    for(g in seq_len(ngroups)){
        gp <- ExtractGroupPars(pars[[g]][[pick]])
        theta <- if(QMC){
            QMC_quad(npts=npts, nfact=nfact, lim = c(0,1))
        } else {
            MC_quad(npts=npts, nfact=nfact, lim = c(0,1))
        }
        Theta[[g]] <- Theta_meanSigma_shift(theta, gp$gmeans, gp$gcov)
    }
    Theta
}

updatePrior <- function(pars, gTheta, list, ngroups, nfact, J,
                        dentype, sitems, cycles, rlist, lrPars = list(), full=FALSE,
                        MC = FALSE){
    prior <- Prior <- Priorbetween <- vector('list', ngroups)
    if(dentype == 'mixture')
        pis <- ExtractMixtures(pars)
    if(dentype == 'custom'){
        for(g in seq_len(ngroups)){
            gp <- pars[[g]][[J+1L]]
            Prior[[g]] <- gp@den(gp, gTheta[[g]])
            Prior[[g]] <- Prior[[g]] / sum(Prior[[g]])
        }
    } else if(dentype == 'discrete'){
        for(g in seq_len(ngroups)){
            gp <- pars[[g]][[J+1L]]
            if(full){
                Prior[[g]] <- gp@den(gp, gTheta[[g]], mus=lrPars@mus)
                Prior[[g]] <- Prior[[g]]/rowSums(Prior[[g]])
            } else {
                Prior[[g]] <- gp@den(gp, gTheta[[g]])
                Prior[[g]] <- Prior[[g]] / sum(Prior[[g]])
            }
        }
    } else if(dentype == 'Davidian'){
        for(g in seq_len(ngroups)){
            gp <- pars[[g]][[J+1L]]
            Prior[[g]] <- gp@den(gp, gTheta[[g]])
            Prior[[g]] <- Prior[[g]] / sum(Prior[[g]])
        }
    } else {
        for(g in seq_len(ngroups)){
            gp <- ExtractGroupPars(pars[[g]][[J+1L]])
            if(dentype == 'bfactor'){
                theta <- pars[[g]][[J+1L]]@theta
                Thetabetween <- pars[[g]][[J+1L]]@Thetabetween
                p <- matrix(0, nrow(gTheta[[g]]), ncol(sitems))
                pp <- matrix(0, nrow(theta), ncol(sitems))
                for(i in seq_len(ncol(sitems))){
                    sel <- c(seq_len(nfact-ncol(sitems)), i + nfact - ncol(sitems))
                    p[,i] <- mirt_dmvnorm(gTheta[[g]][ ,sel], gp$gmeans[sel], gp$gcov[sel,sel,drop=FALSE])
                    pp[,i] <- dnorm(theta, gp$gmeans[sel[length(sel)]],
                                    sqrt(gp$gcov[sel[length(sel)],sel[length(sel)],drop=FALSE]))
                }
                pb <- mirt_dmvnorm(Thetabetween, gp$gmeans[seq_len(ncol(Thetabetween))],
                                   gp$gcov[seq_len(ncol(Thetabetween)), seq_len(ncol(Thetabetween)), drop=FALSE])
                Priorbetween[[g]] <- pb / sum(pb)
                Prior[[g]] <- t(t(p) / colSums(p))
                prior[[g]] <- t(t(pp) / colSums(pp))
                next
            }
            if(full){
                Prior[[g]] <- mirt_dmvnorm(gTheta[[g]][ ,seq_len(nfact),drop=FALSE], lrPars@mus, gp$gcov,
                                           quad=TRUE)
                Prior[[g]] <- Prior[[g]]/rowSums(Prior[[g]])
            } else {
                Prior[[g]] <- mirt_dmvnorm(gTheta[[g]][ ,seq_len(nfact),drop=FALSE], gp$gmeans, gp$gcov)
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
            }
        }
    }
    if(dentype %in% c('EH', 'EHW')){
        if(cycles > 1L){
            for(g in seq_len(ngroups)){
                rr <- if(dentype == 'EHW') standardizeQuadrature(gTheta[[g]], pars[[g]][[J+1L]]@rr,
                                                                 estmean=pars[[g]][[J+1L]]@est[1L],
                                                                 estsd=pars[[g]][[J+1L]]@est[2L])
                    else pars[[g]][[J+1L]]@rr
                Prior[[g]] <- rr / sum(rr)
                attr(Prior[[g]], 'mean_var') <- attr(rr, 'mean_var')
            }
        } else {
            for(g in seq_len(ngroups)){
                Prior[[g]] <- mirt_dmvnorm(gTheta[[g]], 0, matrix(1))
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
            }
        }
    } else if(!is.null(list$customPriorFun)){
        for(g in seq_len(ngroups)){
            Prior[[g]] <- list$customPriorFun(gTheta[[g]], Etable=rlist[[g]][[1L]])
            Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
        }
    }
    if(MC){
        if(full){
            for(g in seq_len(ngroups))
                Prior[[g]] <- matrix(rep(1 / length(gTheta[[g]])),
                                         nrow(lrPars@mus), nrow(gTheta[[g]]))
        } else {
            for(g in seq_len(ngroups))
                Prior[[g]] <- matrix(rep(1 / length(Prior[[g]]), length(Prior[[g]])))
        }
    }
    if(dentype == 'mixture'){
        for(g in seq_len(ngroups))
            Prior[[g]] <- pis[g] * Prior[[g]]
    }
    return(list(prior=prior, Prior=Prior, Priorbetween=Priorbetween))
}

fill_neg_groups_with_complement <- function(OptionalGroups, groupNames){
    has_neg <- grepl("^-", OptionalGroups)
    if(!any(has_neg)) return(OptionalGroups)
    for(i in seq_len(length(has_neg))){
        if(has_neg[i]){
            split <- strsplit(OptionalGroups[i], ',')[[1L]]
            if(!all(grepl("^-", split)))
                stop('Use of negation group syntax (-) cannot be mixed with non-negated syntax',
                     call.=FALSE)
            split <- gsub("^-", "", split)
            OptionalGroups[i] <- paste0(setdiff(groupNames, split), collapse=',')
        }
    }
    OptionalGroups
}

UpdateConstrain <- function(pars, constrain, invariance, nfact, nLambdas, J, ngroups, PrepList,
                            method, itemnames, model, groupNames, mixed.design)
{
    if(!is.numeric(model)){
        groupNames <- as.character(groupNames)
        names(pars) <- groupNames
        model$x[,'OptionalGroups'] <-
            fill_neg_groups_with_complement(model$x[,'OptionalGroups'], groupNames)
        for(row in 1L:nrow(model$x)){
            groupsPicked <- strsplit(model$x[row,'OptionalGroups'], split=',')[[1L]]
            groupsPicked <- which(groupNames %in% groupsPicked)
            input <- model$x[row,2L]
            if(model$x[row,1L] == 'CONSTRAIN'){
                elements <- strsplit(input, '\\),\\(')[[1L]]
                elements <- gsub('\\(', replacement='', x=elements)
                elements <- gsub('\\)', replacement='', x=elements)
                esplit <- strsplit(elements, ',')
                esplit <- lapply(esplit, function(x){
                                newx <- c()
                                if(length(x) < 2L)
                                    stop('CONTRAIN = ... has not been supplied enough arguments', call.=FALSE)
                                for(i in seq_len(length(x)-1L)){
                                    if(grepl('-', x[i])){
                                        tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                                        newx <- c(newx, tmp[1L]:tmp[2L])
                                    } else newx <- c(newx, x[i])
                                }
                                x <- c(newx, x[length(x)])
                                if(x[1L] == 'GROUP') x[1L] <- J + 1L
                                x
                            })
                for(i in seq_len(length(esplit))){
                    for(g in groupsPicked){
                        constr <- c()
                        p <- pars[[g]]
                        sel <- suppressWarnings(
                            as.numeric(esplit[[i]][seq_len(length(esplit[[i]]))]))
                        picknames <- esplit[[i]][is.na(sel)]
                        sel <- na.omit(sel)
                        if(length(picknames) > 1L){
                            constr <- numeric(length(sel))
                            if(length(sel) == 1L){
                                constr <- p[[sel]]@parnum[names(p[[sel]]@est) %in% picknames]
                            } else {
                                if(length(picknames) != length(sel))
                                    stop('Number of items selected not equal to number of parameter names', call.=FALSE)
                                for(j in seq_len(length(sel)))
                                    constr[j] <- p[[sel[j]]]@parnum[names(p[[sel[j]]]@est) == picknames[j]]
                            }
                        } else {
                            for(j in seq_len(length(sel))){
                                pick <- p[[sel[j]]]@parnum[names(p[[sel[j]]]@est) == picknames]
                                if(!length(pick))
                                    stop('CONSTRAIN = ... indexed a parameter that was not relevant for item ', sel[j],
                                         call.=FALSE)
                                constr <- c(constr, pick)
                            }
                        }
                        constrain[[length(constrain) + 1L]] <- constr
                    }
                }
            } else if(model$x[row,1L] == 'CONSTRAINB'){
                if(length(unique(groupNames)) == 1L)
                    stop('CONSTRAINB model argument not valid for single group models', call.=FALSE)
                elements <- strsplit(input, '\\),\\(')[[1L]]
                elements <- gsub('\\(', replacement='', x=elements)
                elements <- gsub('\\)', replacement='', x=elements)
                esplit <- strsplit(elements, ',')
                esplit <- lapply(esplit, function(x){
                    newx <- c()
                    if(length(x) < 2L)
                        stop('CONSTRAINB = ... has not been supplied enough arguments', call.=FALSE)
                    for(i in seq_len(length(x)-1L)){
                        if(grepl('-', x[i])){
                            tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                            newx <- c(newx, tmp[1L]:tmp[2L])
                        } else newx <- c(newx, x[i])
                    }
                    x <- c(newx, x[length(x)])
                    if(x[1L] == 'GROUP') x[1L] <- J + 1L
                    x
                })
                for(i in seq_len(length(esplit))){
                    pickname <- esplit[[i]][length(esplit[[i]])]
                    sel <- as.numeric(esplit[[i]][seq_len(length(esplit[[i]])-1L)])
                    for(j in sel){
                        constr <- c()
                        for(g in groupsPicked){
                            pick <- pars[[g]][[j]]@parnum[names(pars[[g]][[j]]@est) == pickname]
                            if(!length(pick))
                                stop('CONSTRAINB = ... indexed a parameter that was not relevant across groups',
                                     call.=FALSE)
                            constr <- c(constr, pick)
                        }
                        constrain[[length(constrain) + 1L]] <- constr
                    }
                }
            }
        }
    }

    #within group item constraints only
    for(g in seq_len(ngroups))
        for(i in seq_len(length(PrepList[[g]]$constrain)))
            constrain[[length(constrain) + 1L]] <- PrepList[[g]]$constrain[[i]]
    if('covariances' %in% invariance){ #Fix covariance across groups (only makes sense with vars = 1)
        tmp <- c()
        tmpmats <- tmpestmats <- matrix(NA, ngroups, nfact*(nfact+1L)/2)
        for(g in seq_len(ngroups)){
            tmpmats[g,] <- pars[[g]][[J + 1L]]@parnum[(nfact+1L):length(pars[[g]][[J + 1L]]@parnum)]
            tmpestmats[g,] <- pars[[g]][[J + 1L]]@est[(nfact+1L):length(pars[[g]][[J + 1L]]@est)]
        }
        select <- colSums(tmpestmats) == ngroups
        for(i in seq_len(length(select)))
            if(select[i])
                constrain[[length(constrain) + 1L]] <- tmpmats[seq_len(ngroups), i]

    }
    if(any(itemnames %in% invariance)){
        matched <- na.omit(match(invariance, itemnames))
        for(i in matched){
            jj <- sum(pars[[1L]][[i]]@est)
            if(!(jj > 0))
                stop('Equality constraints applied to items where no parameters were estimated. Please fix',
                     call.=FALSE)
            for(j in seq_len(jj)){
                tmp <- c()
                for(g in seq_len(ngroups))
                    tmp <- c(tmp, pars[[g]][[i]]@parnum[pars[[g]][[i]]@est][j])
                constrain[[length(constrain) + 1L]] <- tmp
            }
        }
    }
    if('slopes' %in% invariance){ #Equal factor loadings
        tmpmats <- tmpests <- list()
        for(g in seq_len(ngroups))
            tmpmats[[g]] <- tmpests[[g]] <- matrix(NA, J, nLambdas)
        for(g in seq_len(ngroups)){
            for(i in seq_len(J)){
                tmpmats[[g]][i,] <- pars[[g]][[i]]@parnum[1L:nLambdas]
                tmpests[[g]][i,] <- pars[[g]][[i]]@est[1L:nLambdas]
            }
        }
        for(i in seq_len(J)){
            for(j in seq_len(nLambdas)){
                tmp <- c()
                for(g in seq_len(ngroups)){
                    if(tmpests[[1L]][[i, j]])
                        tmp <- c(tmp, tmpmats[[g]][i,j])
                }
                constrain[[length(constrain) + 1L]] <- tmp
            }
        }
    }
    if('intercepts' %in% invariance){ #Equal item intercepts (and all other item pars)
        tmpmats <- tmpests <- list()
        for(g in seq_len(ngroups))
            tmpmats[[g]] <- tmpests[[g]] <- list()
        for(g in seq_len(ngroups)){
            for(i in seq_len(J)){
                ind <- (nLambdas+1L):length(pars[[g]][[i]]@parnum)
                if(is(pars[[g]][[i]], 'dich')) ind <- ind[seq_len(length(ind)-2L)]
                if(is(pars[[g]][[i]], 'partcomp')) ind <- ind[seq_len(length(ind)-1L)]
                tmpmats[[g]][[i]] <- pars[[g]][[i]]@parnum[ind]
                tmpests[[g]][[i]] <- pars[[g]][[i]]@est[ind]
            }
        }
        for(i in seq_len(J)){
            for(j in seq_len(length(tmpmats[[1L]][[i]]))){
                tmp <- c()
                for(g in seq_len(ngroups)){
                    if(tmpests[[1L]][[i]][j])
                        tmp <- c(tmp, tmpmats[[g]][[i]][j])
                }
                constrain[[length(constrain) + 1L]] <- tmp
            }
        }
    }
    if(!is.null(mixed.design) && mixed.design$from == 'mirt'){
        for(g in seq_len(ngroups)){
            are_partcomp <- sapply(pars[[g]], class) == 'partcomp'
            are_partcomp <- are_partcomp[-length(are_partcomp)] # no group
            cpow_unique <- NULL
            if(any(are_partcomp & mixed.design$has_idesign)){
                cpowmat <- lapply(1:length(are_partcomp), \(i){
                    if(are_partcomp[i]){
                        pars[[g]][[i]]@cpow
                    } else {
                        NA
                    }})
                cpowmat <- do.call(rbind, cpowmat)
                cpowgroup <- apply(cpowmat, 1, \(x) paste0(x, collapse=','))
                cpowgroup[rowSums(is.na(cpowmat)) == ncol(cpowmat)] <- NA
                cpow_unique <- unique(na.omit(cpowgroup))
            }
            if(!is.null(cpow_unique)){
                for(cp in cpow_unique){
                    for(p in 1:ncol(mixed.design$fixed)){
                        constr <- rep(NA, J)
                        for(i in which(mixed.design$has_idesign & are_partcomp)){
                            if(cp == cpowgroup[i] && pars[[g]][[i]]@est[p])
                                constr[i] <- pars[[g]][[i]]@parnum[p]
                        }
                        constr <- na.omit(constr)
                        if(length(constr))
                            constrain[[length(constrain) + 1L]] <- constr
                    }
                }
            }
            if(any(!are_partcomp & mixed.design$has_idesign)){
                for(p in 1:ncol(mixed.design$fixed)){
                    constr <- rep(NA, J)
                    for(i in which(mixed.design$has_idesign & !are_partcomp))
                        if(pars[[g]][[i]]@est[p])
                            constr[i] <- pars[[g]][[i]]@parnum[p]
                    constr <- na.omit(constr)
                    if(length(constr))
                        constrain[[length(constrain) + 1L]] <- constr
                }
            }

        }
    }
    #remove redundant constraints
    redun <- rep(FALSE, length(constrain))
    if(length(constrain)){
        for(i in seq_len(length(redun))){
            while(TRUE){
                lastredun <- redun
                for(j in seq_len(length(redun))){
                    if(i < j && !redun[j] && !redun[i]){
                        if(any(constrain[[i]] %in% constrain[[j]])){
                            constrain[[i]] <- unique(c(constrain[[i]], constrain[[j]]))
                            redun[j] <- TRUE
                        }
                    }
                }
                if(all(lastredun == redun)) break
            }
        }
    }
    constrain[redun] <- NULL
    return(constrain)
}

expbeta_sv <- function(val1, val2){
    ret <- qlogis((val1-1)/(val1 + val2-2))
    if(!is.finite(ret)) ret <- qlogis(val1/(val1 + val2))
    ret
}

UpdateParameters <- function(PrepList, model, groupNames){
    if(!is.numeric(model)){
        nitems <- length(PrepList[[1L]]$pars) - 1L
        model$x[,"Parameters"] <- gsub("\\(GROUP,",
                                       replacement = sprintf("(%i,", nitems + 1L),
                                       model$x[,"Parameters"])
        groupNames <- as.character(groupNames)
        pars <- vector('list', length(PrepList))
        for(g in seq_len(length(PrepList)))
            pars[[g]] <- PrepList[[g]]$pars
        names(pars) <- groupNames
        for(row in 1L:nrow(model$x)){
            groupsPicked <- strsplit(model$x[row,'OptionalGroups'], split=',')[[1L]]
            groupsPicked <- which(groupNames %in% groupsPicked)
            input <- model$x[row,2L]
            if(model$x[row,1L] == 'START'){
                elements <- strsplit(input, '\\),\\(')[[1L]]
                elements <- gsub('\\(', replacement='', x=elements)
                elements <- gsub('\\)', replacement='', x=elements)
                esplit <- strsplit(elements, ',')
                esplit <- lapply(esplit, function(x){
                    newx <- c()
                    if(length(x) < 3L)
                        stop('START = ... has not been supplied enough arguments', call.=FALSE)
                    for(i in seq_len(length(x)-2L)){
                        if(grepl('-', x[i])){
                            tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                            newx <- c(newx, tmp[1L]:tmp[2L])
                        } else newx <- c(newx, x[i])
                    }
                    x <- c(newx, x[length(x)-1L], x[length(x)])
                    x
                })
                picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-2)]))
                for(i in seq_len(length(picks))){
                    for(gpick in groupsPicked){
                        tmp <- pars[[gpick]][picks[[i]]]
                        len <- length(esplit[[i]])
                        tmp <- lapply(tmp, function(x, which, val){
                            if(which %in% c('g', 'u')) val <- qlogis(val)
                            x@par[names(x@est) == which] <- val
                            x
                        }, which=esplit[[i]][len-1L], val = as.numeric(esplit[[i]][len]))
                        pars[[gpick]][picks[[i]]] <- tmp
                    }
                }
            } else if(model$x[row,1L] == 'FIXED'){
                elements <- strsplit(input, '\\),\\(')[[1L]]
                elements <- gsub('\\(', replacement='', x=elements)
                elements <- gsub('\\)', replacement='', x=elements)
                esplit <- strsplit(elements, ',')
                esplit <- lapply(esplit, function(x){
                    newx <- c()
                    if(length(x) < 2L)
                        stop('FIXED = ... has not been supplied enough arguments', call.=FALSE)
                    for(i in seq_len(length(x)-1L)){
                        if(grepl('-', x[i])){
                            tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                            newx <- c(newx, tmp[1L]:tmp[2L])
                        } else newx <- c(newx, x[i])
                    }
                    x <- c(newx, x[length(x)])
                    x
                })
                picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-1L)]))
                for(i in seq_len(length(picks))){
                    for(gpick in groupsPicked){
                        tmp <- pars[[gpick]][picks[[i]]]
                        len <- length(esplit[[i]])
                        tmp <- lapply(tmp, function(x, which){
                            x@est[names(x@est) == which] <- FALSE
                            x
                        }, which=esplit[[i]][len])
                        pars[[gpick]][picks[[i]]] <- tmp
                    }
                }
            } else if(model$x[row,1L] == 'FREE'){
                elements <- strsplit(input, '\\),\\(')[[1L]]
                elements <- gsub('\\(', replacement='', x=elements)
                elements <- gsub('\\)', replacement='', x=elements)
                esplit <- strsplit(elements, ',')
                esplit <- lapply(esplit, function(x){
                    newx <- c()
                    if(length(x) < 2L)
                        stop('FREE = ... has not been supplied enough arguments', call.=FALSE)
                    for(i in seq_len(length(x)-1L)){
                        if(grepl('-', x[i])){
                            tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                            newx <- c(newx, tmp[1L]:tmp[2L])
                        } else newx <- c(newx, x[i])
                    }
                    x <- c(newx, x[length(x)])
                    x
                })
                picks <- lapply(esplit, function(x) as.integer(x[seq_len(length(x)-1L)]))
                for(i in seq_len(length(picks))){
                    for(gpick in groupsPicked){
                        tmp <- pars[[gpick]][picks[[i]]]
                        len <- length(esplit[[i]])
                        tmp <- lapply(tmp, function(x, which){
                            x@est[names(x@est) == which] <- TRUE
                            x
                        }, which=esplit[[i]][len])
                        pars[[gpick]][picks[[i]]] <- tmp
                    }
                }
            } else if(model$x[row,1L] == 'LBOUND'){
                elements <- strsplit(input, '\\),\\(')[[1L]]
                elements <- gsub('\\(', replacement='', x=elements)
                elements <- gsub('\\)', replacement='', x=elements)
                esplit <- strsplit(elements, ',')
                esplit <- lapply(esplit, function(x){
                    newx <- c()
                    if(length(x) < 3L)
                        stop('LBOUND = ... has not been supplied enough arguments', call.=FALSE)
                    for(i in seq_len(length(x)-2L)){
                        if(grepl('-', x[i])){
                            tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                            newx <- c(newx, tmp[1L]:tmp[2L])
                        } else newx <- c(newx, x[i])
                    }
                    x <- c(newx, x[length(x)-1L], x[length(x)])
                    x
                })
                picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-2)]))
                for(i in seq_len(length(picks))){
                    for(gpick in groupsPicked){
                        tmp <- pars[[gpick]][picks[[i]]]
                        len <- length(esplit[[i]])
                        tmp <- lapply(tmp, function(x, which, val){
                            if(which %in% c('g', 'u')) val <- qlogis(val)
                            x@lbound[names(x@est) == which] <- val
                            x
                        }, which=esplit[[i]][len-1L], val = as.numeric(esplit[[i]][len]))
                        pars[[gpick]][picks[[i]]] <- tmp
                    }
                }
            } else if(model$x[row,1L] == 'UBOUND'){
                elements <- strsplit(input, '\\),\\(')[[1L]]
                elements <- gsub('\\(', replacement='', x=elements)
                elements <- gsub('\\)', replacement='', x=elements)
                esplit <- strsplit(elements, ',')
                esplit <- lapply(esplit, function(x){
                    newx <- c()
                    if(length(x) < 3L)
                        stop('UBOUND = ... has not been supplied enough arguments', call.=FALSE)
                    for(i in seq_len(length(x)-2L)){
                        if(grepl('-', x[i])){
                            tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                            newx <- c(newx, tmp[1L]:tmp[2L])
                        } else newx <- c(newx, x[i])
                    }
                    x <- c(newx, x[length(x)-1L], x[length(x)])
                    x
                })
                picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-2)]))
                for(i in seq_len(length(picks))){
                    for(gpick in groupsPicked){
                        tmp <- pars[[gpick]][picks[[i]]]
                        len <- length(esplit[[i]])
                        tmp <- lapply(tmp, function(x, which, val){
                            if(which %in% c('g', 'u')) val <- qlogis(val)
                            x@ubound[names(x@est) == which] <- val
                            x
                        }, which=esplit[[i]][len-1L], val = as.numeric(esplit[[i]][len]))
                        pars[[gpick]][picks[[i]]] <- tmp
                    }
                }
            } else if(model$x[row,1L] == 'PRIOR'){
                input <- gsub(' ', replacement='', x=input)
                elements <- strsplit(input, '\\),\\(')[[1L]]
                elements <- gsub('\\(', replacement='', x=elements)
                elements <- gsub('\\)', replacement='', x=elements)
                esplit <- strsplit(elements, ',')
                esplit <- lapply(esplit, function(x){
                    newx <- c()
                    if(length(x) < 5L)
                        stop('PRIOR = ... has not been supplied enough arguments', call.=FALSE)
                    for(i in seq_len(length(x)-2L)){
                        if(grepl('-', x[i])){
                            tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                            newx <- c(newx, tmp[1L]:tmp[2L])
                        } else newx <- c(newx, x[i])
                    }
                    x <- c(newx, x[length(x)-1L], x[length(x)])
                    x
                })
                picks <- lapply(esplit, function(x) as.integer(x[1L:(length(x)-4L)]))
                for(i in seq_len(length(picks))){
                    for(gpick in groupsPicked){
                        tmp <- pars[[gpick]][picks[[i]]]
                        len <- length(esplit[[i]])
                        tmp <- lapply(tmp, function(x, name, type, val1, val2){
                            if(!(type %in% c('norm', 'beta', 'lnorm', 'expbeta')))
                                stop('Prior type specified in PRIOR = ... not available', call.=FALSE)
                            type <- switch(type, norm=1L, lnorm=2L, beta=3L, expbeta=4L, 0L)
                            which <- names(x@est) == name
                            if(!any(which)) stop('Parameter \'', name, '\' does not exist for respective item',
                                                 call.=FALSE)
                            x@any.prior <- TRUE
                            x@prior.type[which] <- type
                            x@prior_1[which] <- val1
                            x@prior_2[which] <- val2
                            x@par[which] <- switch(type,
                                                   '1'=val1,
                                                   '2'=exp(val1),
                                                   '3'=(val1-1)/(val1 + val2 - 2),
                                                   '4'=expbeta_sv(val1, val2))
                            if(type == '2')
                                x@lbound[which] <- 0
                            if(type == '3'){
                                x@lbound[which] <- 0
                                x@ubound[which] <- 1
                            }
                            x
                        }, name = esplit[[i]][len-3L],
                        type = esplit[[i]][len-2L],
                        val1 = as.numeric(esplit[[i]][len-1L]),
                        val2 = as.numeric(esplit[[i]][len]))
                        pars[[gpick]][picks[[i]]] <- tmp
                    }
                }
            }
        }
        for(g in seq_len(length(PrepList)))
            PrepList[[g]]$pars <- pars[[g]]
    }
    return(PrepList)
}

buildModelSyntax <- function(model, J, groupNames, itemtype){
    exploratory <- FALSE
    if(is(model, 'mirt.model') && any(model$x[,1L] == 'NEXPLORE')){
        oldmodel <- model
        model <- as.integer(model$x[model$x[,1L] == 'NEXPLORE', 2L])
        if(model != 1L) exploratory <- TRUE
        tmp <- tempfile('tempfile')
        for(i in 1L:model)
            cat(paste('F', i,' = 1-', (J-i+1L), "\n", sep=''), file=tmp, append = TRUE)
        model <- mirt.model(file=tmp, quiet = TRUE)
        model$x <- rbind(model$x, oldmodel$x[oldmodel$x[,1L] != 'NEXPLORE', ])
        unlink(tmp)
    } else if((is(model, 'numeric') && length(model) == 1L)){
        if(any(itemtype == 'lca')){
            tmp <- tempfile('tempfile')
            for(i in 1L:model)
                cat(paste('F', i,' = 1-', J, "\n", sep=''), file=tmp, append = TRUE)
            model <- mirt.model(file=tmp, quiet = TRUE)
            unlink(tmp)
        } else {
            if(model != 1L) exploratory <- TRUE
            tmp <- tempfile('tempfile')
            for(i in 1L:model)
                cat(paste('F', i,' = 1-', (J-i+1L), "\n", sep=''), file=tmp, append = TRUE)
            model <- mirt.model(file=tmp, quiet = TRUE)
            unlink(tmp)
        }
    }
    if(is(model, 'numeric') && length(model) > 1L)
        model <- bfactor2mod(model, J)
    if(!is.numeric(model)){
        model$x <- cbind(model$x, OptionalGroups="")
        for(i in 1L:nrow(model$x)){
            brackets <- grepl('\\[', model$x[i, 'Type'])
            if(!brackets){
                model$x[i,"OptionalGroups"] <- paste0(as.character(groupNames), collapse = ',')
            } else {
                tmp <- strsplit(model$x[i, 'Type'], '\\[')[[1L]]
                tmp[2L] <- gsub("\\]", "", tmp[2L])
                tmp[2L] <- gsub(" ", "", tmp[2L])
                model$x[i,"OptionalGroups"] <- tmp[2L]
                model$x[i,"Type"] <- tmp[1L]
            }
        }
        model$x[,'Type'] <- gsub(" ", "", model$x[,'Type'])
    }
    attr(model, 'exploratory') <- exploratory
    model
}

resetPriorConstrain <- function(pars, constrain, nconstrain){
    nitems <- length(pars[[1]])
    if(length(constrain)){
        for(g in seq_len(length(pars))){
            for(i in seq_len(nitems)){
                for(ci in seq_len(length(constrain))){
                    pick <- pars[[g]][[i]]@parnum %in% constrain[[ci]][-1L]
                    pars[[g]][[i]]@prior.type[pick] <- 0L
                    pars[[g]][[i]]@any.prior <- any(pars[[g]][[i]]@prior.type > 0L)
                }
            }
        }
    }
    pars
}

ReturnPars <- function(PrepList, itemnames, random, lrPars, clist, nclist,
                       itemtype, lr.random = NULL, MG = FALSE){
    parnum <- par <- est <- item <- parname <- gnames <- class <-
        lbound <- ubound <- prior.type <- prior_1 <- prior_2 <- c()
    if(!MG) PrepList <- list(full=PrepList)
    gnames_count <- integer(length(PrepList))
    for(g in seq_len(length(PrepList))){
        tmpgroup <- PrepList[[g]]$pars
        for(i in seq_len(length(tmpgroup))){
            if(i <= length(itemnames))
                item <- c(item, rep(itemnames[i], length(tmpgroup[[i]]@parnum)))
            class <- c(class, rep(class(tmpgroup[[i]]), length(tmpgroup[[i]]@parnum)))
            parname <- c(parname, tmpgroup[[i]]@parnames)
            parnum <- c(parnum, tmpgroup[[i]]@parnum)
            par <- c(par, tmpgroup[[i]]@par)
            est <- c(est, tmpgroup[[i]]@est)
            lbound <- c(lbound, tmpgroup[[i]]@lbound)
            ubound <- c(ubound, tmpgroup[[i]]@ubound)
            tmp <- sapply(as.character(tmpgroup[[i]]@prior.type),
                                 function(x) switch(x, '1'='norm', '2'='lnorm',
                                                    '3'='beta', '4'='expbeta', 'none'))
            gnames_count[g] <- gnames_count[g] + length(tmp)
            prior.type <- c(prior.type, tmp)
            prior_1 <- c(prior_1, tmpgroup[[i]]@prior_1)
            prior_2 <- c(prior_2, tmpgroup[[i]]@prior_2)
        }
        item <- c(item, rep('GROUP', length(tmpgroup[[i]]@parnum)))
    }
    for(i in seq_len(length(random))){
        parname <- c(parname, random[[i]]@parnames)
        parnum <- c(parnum, random[[i]]@parnum)
        par <- c(par, random[[i]]@par)
        est <- c(est, random[[i]]@est)
        lbound <- c(lbound, random[[i]]@lbound)
        ubound <- c(ubound, random[[i]]@ubound)
        tmp <- sapply(as.character(random[[i]]@prior.type),
                      function(x) switch(x, '1'='norm', '2'='lnorm',
                                         '3'='beta', '4'='expbeta', 'none'))
        prior.type <- c(prior.type, tmp)
        prior_1 <- c(prior_1, random[[i]]@prior_1)
        prior_2 <- c(prior_2, random[[i]]@prior_2)
        class <- c(class, rep('RandomPars', length(random[[i]]@parnum)))
        item <- c(item, rep('RANDOM', length(random[[i]]@parnum)))
        gnames_count[g] <- gnames_count[g] + length(tmp) # TODO this assumes one group
    }
    if(length(lrPars)){
        parname <- c(parname, lrPars@parnames)
        parnum <- c(parnum, lrPars@parnum)
        par <- c(par, lrPars@par)
        est <- c(est, lrPars@est)
        lbound <- c(lbound, lrPars@lbound)
        ubound <- c(ubound, lrPars@ubound)
        tmp <- sapply(as.character(lrPars@prior.type),
                      function(x) switch(x, '1'='norm', '2'='lnorm',
                                         '3'='beta', '4'='expbeta', 'none'))
        prior.type <- c(prior.type, tmp)
        prior_1 <- c(prior_1, lrPars@prior_1)
        prior_2 <- c(prior_2, lrPars@prior_2)
        class <- c(class, rep('lrPars', length(lrPars@parnum)))
        item <- c(item, rep('BETA', length(lrPars@parnum)))
        gnames_count[g] <- gnames_count[g] + length(tmp)
    }
    for(i in seq_len(length(lr.random))){
        parname <- c(parname, lr.random[[i]]@parnames)
        parnum <- c(parnum, lr.random[[i]]@parnum)
        par <- c(par, lr.random[[i]]@par)
        est <- c(est, lr.random[[i]]@est)
        lbound <- c(lbound, lr.random[[i]]@lbound)
        ubound <- c(ubound, lr.random[[i]]@ubound)
        tmp <- sapply(as.character(lr.random[[i]]@prior.type),
                      function(x) switch(x, '1'='norm', '2'='lnorm',
                                         '3'='beta', '4'='expbeta', 'none'))
        prior.type <- c(prior.type, tmp)
        prior_1 <- c(prior_1, lr.random[[i]]@prior_1)
        prior_2 <- c(prior_2, lr.random[[i]]@prior_2)
        class <- c(class, rep('LRRandomPars', length(lr.random[[i]]@parnum)))
        item <- c(item, rep('LRRANDOM', length(lr.random[[i]]@parnum)))
        gnames_count[g] <- gnames_count[g] + length(tmp)
    }
    gnames <- rep(names(PrepList), times = gnames_count)
    par[parname %in% c('g', 'u')] <- antilogit(par[parname %in% c('g', 'u')])
    lbound[parname %in% c('g', 'u')] <- antilogit(lbound[parname %in% c('g', 'u')])
    ubound[parname %in% c('g', 'u')] <- antilogit(ubound[parname %in% c('g', 'u')])
    constrain <- nconstrain <- rep("none", length(gnames))
    if(length(clist)){
        for(i in seq_len(length(clist)))
            constrain[clist[[i]]] <- i
    }
    if(length(nclist)){
        for(i in seq_len(length(nclist)))
            nconstrain[nclist[[i]]] <- i
    }
    ret <- data.frame(group=gnames, item=item, class=class, name=parname, parnum=parnum,
                      value=par, lbound=lbound, ubound=ubound, est=est, const=constrain, nconst=nconstrain,
                      prior.type=prior.type, prior_1=prior_1, prior_2=prior_2, stringsAsFactors = FALSE)
    ret <- as.mirt_df(ret)
    attr(ret, 'itemtype') <- itemtype
    ret
}

UpdatePrepList <- function(PrepList, pars, random, lr.random, clist, nclist,
                           itemtype, lrPars = list(), MG = FALSE){
    currentDesign <- ReturnPars(PrepList, PrepList[[1L]]$itemnames, random=random,
                                clist=clist, nclist=nclist, itemtype=itemtype,
                                lrPars=lrPars, lr.random=lr.random, MG = TRUE)
    if(nrow(currentDesign) != nrow(pars))
        stop('Rows in supplied and starting value data.frame objects do not match. Were the
             data or itemtype input arguments modified?', call.=FALSE)
    if(!all(as.matrix(currentDesign[,c('group', 'class', 'name', 'parnum')]) ==
                as.matrix(pars[,c('group', 'class', 'name', 'parnum')])))
        stop('Critical internal parameter labels do not match those returned from pars = \'values\'',
             call.=FALSE)
    if(!all(sapply(currentDesign, class) == sapply(pars, class)))
        stop('pars input does not contain the appropriate classes, which should match pars = \'values\'',
             call.=FALSE)
    if(!all(unique(pars$prior.type) %in% c('none', 'norm', 'beta', 'lnorm', 'expbeta')))
        stop('prior.type input in pars contains invalid prior types', call.=FALSE)
    if(!MG) PrepList <- list(PrepList)
    pars$value[pars$name %in% c('g', 'u')] <- logit(pars$value[pars$name %in% c('g', 'u')])
    pars$lbound[pars$name %in% c('g', 'u')] <- logit(pars$lbound[pars$name %in% c('g', 'u')])
    pars$ubound[pars$name %in% c('g', 'u')] <- logit(pars$ubound[pars$name %in% c('g', 'u')])
    if(PrepList[[1L]]$nfact > 1L){
        mat <- matrix(pars$est[pars$name %in% paste0('a', seq_len(PrepList[[1L]]$nfact))],
                      ncol = PrepList[[1L]]$nfact, byrow=TRUE)
        PrepList[[1L]]$exploratory <- all(sort(colSums(!mat)) == seq(0L, PrepList[[1L]]$nfact - 1L))
    }
    ind <- 1L
    for(g in seq_len(length(PrepList))){
        for(i in seq_len(length(PrepList[[g]]$pars))){
            for(j in seq_len(length(PrepList[[g]]$pars[[i]]@par))){
                PrepList[[g]]$pars[[i]]@par[j] <- pars[ind,'value']
                PrepList[[g]]$pars[[i]]@est[j] <- as.logical(pars[ind,'est'])
                PrepList[[g]]$pars[[i]]@lbound[j] <- pars[ind,'lbound']
                PrepList[[g]]$pars[[i]]@ubound[j] <- pars[ind,'ubound']
                tmp <- as.character(pars[ind,'prior.type'])
                PrepList[[g]]$pars[[i]]@prior.type[j] <-
                    switch(tmp, norm=1L, lnorm=2L, beta=3L, expbeta=4L, 0L)
                PrepList[[g]]$pars[[i]]@prior_1[j] <- pars[ind,'prior_1']
                PrepList[[g]]$pars[[i]]@prior_2[j] <- pars[ind,'prior_2']
                ind <- ind + 1L
            }
            if(is(PrepList[[g]]$pars[[i]], 'graded')){
                tmp <- ExtractZetas(PrepList[[g]]$pars[[i]])
                if(!all(tmp == sort(tmp, decreasing=TRUE)) || length(unique(tmp)) != length(tmp))
                    stop('Graded model intercepts for item ', i, ' in group ', g,
                         ' do not descend from highest to lowest. Please fix', call.=FALSE)
            }
            PrepList[[g]]$pars[[i]]@any.prior <- any(1L:3L %in%
                                                         PrepList[[g]]$pars[[i]]@prior.type)
        }
    }
    if(length(random)){
        for(i in seq_len(length(random))){
            for(j in seq_len(length(random[[i]]@par))){
                random[[i]]@par[j] <- pars[ind,'value']
                random[[i]]@est[j] <- as.logical(pars[ind,'est'])
                random[[i]]@lbound[j] <- pars[ind,'lbound']
                random[[i]]@ubound[j] <- pars[ind,'ubound']
                ind <- ind + 1L
            }
        }
        attr(PrepList, 'random') <- random
    }
    if(length(lrPars)){
        for(j in seq_len(length(lrPars@par))){
            lrPars@par[j] <- pars[ind,'value']
            lrPars@est[j] <- as.logical(pars[ind,'est'])
            lrPars@lbound[j] <- pars[ind,'lbound']
            lrPars@ubound[j] <- pars[ind,'ubound']
            ind <- ind + 1L
        }
    }
    if(length(lr.random)){
        for(i in seq_len(length(lr.random))){
            for(j in seq_len(length(lr.random[[i]]@par))){
                lr.random[[i]]@par[j] <- pars[ind,'value']
                lr.random[[i]]@est[j] <- as.logical(pars[ind,'est'])
                lr.random[[i]]@lbound[j] <- pars[ind,'lbound']
                lr.random[[i]]@ubound[j] <- pars[ind,'ubound']
                ind <- ind + 1L
            }
        }
        attr(PrepList, 'lr.random') <- lr.random
    }
    if(!MG) PrepList <- PrepList[[1L]]
    return(PrepList)
}

rebuild_clist <- function(parnum, cvec){
    ret <- list()
    if(!all(cvec == 'none')){
        uniq <- unique(cvec)
        uniq <- uniq[uniq != 'none']
        for(i in seq_len(length(uniq)))
            ret[[i]] <- parnum[uniq[i] == cvec]
    }
    ret
}

#new gradient and hessian with priors
DerivativePriors <- function(x, grad, hess){
    if(any(x@prior.type %in% 1L)){ #norm
        ind <- x@prior.type %in% 1L
        val <- x@par[ind]
        mu <- x@prior_1[ind]
        s <- x@prior_2[ind]
        g <- -(val - mu)/(s^2)
        h <- -1/(s^2)
        grad[ind] <- grad[ind] + g
        if(length(val) == 1L) hess[ind, ind] <- hess[ind, ind] + h
        else diag(hess[ind, ind]) <- diag(hess[ind, ind]) + h
    }
    if(any(x@prior.type %in% 2L)){ #lnorm
        ind <- x@prior.type %in% 2L
        val <- x@par[ind]
        val <- ifelse(val > 0, val, 1e-10)
        lval <- log(val)
        mu <- x@prior_1[ind]
        s <- x@prior_2[ind]
        g <- -(lval - mu)/(val * s^2) - 1/val
        h <- 1/(val^2) - 1/(val^2 * s^2) - (lval - mu)/(val^2 * s^2)
        grad[ind] <- grad[ind] + g
        if(length(val) == 1L) hess[ind, ind] <- hess[ind, ind] + h
        else diag(hess[ind, ind]) <- diag(hess[ind, ind]) + h
    }
    if(any(x@prior.type %in% c(3L, 4L))){ #beta
        tmp <- x@par
        ind <- x@prior.type == 4L
        tmp[ind] <- plogis(tmp[ind])
        ind <- x@prior.type %in% c(3L, 4L)
        val <- tmp[ind]
        val <- ifelse(val < 1e-10, 1e-10, val)
        val <- ifelse(val > 1-1e-10, 1-1e-10, val)
        a <- x@prior_1[ind]
        b <- x@prior_2[ind]
        g <- (a - 1)/val - (b-1)/(1-val)
        h <- -(a - 1)/(val^2) - (b-1) / (1-val)^2
        grad[ind] <- grad[ind] + g
        if(length(val) == 1L) hess[ind, ind] <- hess[ind, ind] + h
        else diag(hess[ind, ind]) <- diag(hess[ind, ind]) + h
    }
    return(list(grad=grad, hess=hess))
}

#new likelihood with priors
LL.Priors <- function(x, LL){
    if(any(x@prior.type %in% 1L)){
        ind <- x@prior.type %in% 1L
        val <- x@par[ind]
        u <- x@prior_1[ind]
        s <- x@prior_2[ind]
        for(i in seq_len(length(val))){
            tmp <- dnorm(val[i], u[i], s[i], log=TRUE)
            LL <- LL + ifelse(tmp == -Inf, log(1e-100), tmp)
        }
    }
    if(any(x@prior.type %in% 2L)){
        ind <- x@prior.type %in% 2L
        val <- x@par[ind]
        u <- x@prior_1[ind]
        s <- x@prior_2[ind]
        for(i in seq_len(length(val))){
            if(val[i] > 0)
                LL <- LL + dlnorm(val[i], u[i], s[i], log=TRUE)
            else LL <- LL + log(1e-100)
        }
    }
    if(any(x@prior.type %in% 3L)){
        ind <- x@prior.type %in% 3L
        val <- x@par[ind]
        a <- x@prior_1[ind]
        b <- x@prior_2[ind]
        for(i in seq_len(length(val))){
            if(val[i] > 0 && val[i] < 1)
                LL <- LL + dbeta(val[i], a[i], b[i], log=TRUE)
            else LL <- LL + log(1e-100)
        }
    }
    if(any(x@prior.type %in% 4L)){
        ind <- x@prior.type %in% 4L
        val <- plogis(x@par[ind])
        a <- x@prior_1[ind]
        b <- x@prior_2[ind]
        for(i in seq_len(length(val))){
            if(val[i] > 0 && val[i] < 1)
                LL <- LL + dbeta(val[i], a[i], b[i], log=TRUE)
            else LL <- LL + log(1e-100)
        }
    }
    return(LL)
}

ItemInfo <- function(x, Theta, cosangle, total.info = TRUE){
    P <- ProbTrace(x, Theta)
    dx <- DerivTheta(x, Theta)
    info <- matrix(0, nrow(Theta), ncol(P))
    cosanglefull <- matrix(cosangle, nrow(P), length(cosangle), byrow = TRUE)
    for(i in seq_len(x@ncat))
        dx$grad[[i]] <- matrix(rowSums(dx$grad[[i]] * cosanglefull))
    for(i in seq_len(x@ncat))
        info[,i] <- (dx$grad[[i]])^2 / P[ ,i]
    if(total.info) info <- rowSums(info)
    return(info)
}

ItemInfo2 <- function(x, Theta, total.info = TRUE, MD = FALSE, DERIV = NULL, P = NULL){
    if(is.null(P)) P <- ProbTrace(x, Theta)
    dx <- if(is.null(DERIV)) DerivTheta(x, Theta) else DERIV(x, Theta)
    if(MD){
        info <- matrix(0, length(Theta), length(Theta))
        for(i in seq_len(x@ncat))
            info <- info + outer(as.numeric(dx$grad[[i]]), as.numeric(dx$grad[[i]])) / P[ ,i]
    } else {
        grad <- do.call(cbind, dx$grad)
        hess <- do.call(cbind, dx$hess)
        info <- grad^2 / P - hess
        if(total.info) info <- rowSums(info)
    }
    return(info)
}

nameInfoMatrix <- function(info, correction, L, npars){
    #give info meaningful names for wald test
    parnames <- names(correction)
    tmp <- outer(seq_len(npars), rep(1L, npars))
    matind <- matrix(0, ncol(tmp), nrow(tmp))
    matind[lower.tri(matind, diag = TRUE)] <- tmp[lower.tri(tmp, diag = TRUE)]
    matind <- matind * as.matrix(L) # TODO as.matrix could be avoided
    matind[matind == 0 ] <- NA
    matind[!is.na(matind)] <- tmp[!is.na(matind)]
    shortnames <- c()
    for(i in seq_len(length(correction))){
        keep <- is.na(matind[ , 1L])
        while(all(keep)){
            matind <- matind[-1L, -1L, drop = FALSE]
            keep <- is.na(matind[ , 1L])
        }
        tmp <- paste0(parnames[i], paste0('.', matind[!keep, 1L], collapse=''))
        shortnames <- c(shortnames, tmp)
        matind <- matind[keep, keep, drop = FALSE]
    }
    colnames(info) <- rownames(info) <- shortnames
    return(info)
}

maketabData <- function(tmpdata, group, groupNames, nitem, K, itemloc,
                        Names, itemnames, survey.weights){
    tmpdata[is.na(tmpdata)] <- 99999L
    stringfulldata <- apply(tmpdata, 1L, paste, sep='', collapse = '/')
    stringtabdata <- unique(stringfulldata)
    tabdata2 <- lapply(strsplit(stringtabdata, split='/'), as.integer)
    tabdata2 <- do.call(rbind, tabdata2)
    tabdata2[tabdata2 == 99999L] <- NA
    tabdata <- matrix(0L, nrow(tabdata2), sum(K))
    for(i in seq_len(nitem)){
        uniq <- sort(na.omit(unique(tabdata2[,i])))
        if(length(uniq) < K[i]) uniq <- 0L:(K[i]-1L)
        for(j in seq_len(length(uniq)))
            tabdata[,itemloc[i] + j - 1L] <- as.integer(tabdata2[,i] == uniq[j])
    }
    tabdata[is.na(tabdata)] <- 0L
    colnames(tabdata) <- Names
    colnames(tabdata2) <- itemnames
    groupFreq <- vector('list', length(groupNames))
    names(groupFreq) <- groupNames
    for(g in seq_len(length(groupNames))){
        Freq <- integer(length(stringtabdata))
        tmpstringdata <- stringfulldata[group == groupNames[g]]
        if(!is.null(survey.weights)){
            Freq <- mySapply(seq_len(nrow(tabdata)), function(x, std, tstd, w)
                sum(w[stringtabdata[x] == tstd]), std=stringtabdata, tstd=tmpstringdata,
                w=survey.weights[group == groupNames[g]])
        } else {
            Freq[stringtabdata %in% tmpstringdata] <- as.integer(table(
                match(tmpstringdata, stringtabdata)))
        }
        groupFreq[[g]] <- Freq
    }
    ret <- list(tabdata=tabdata, tabdata2=tabdata2, Freq=groupFreq)
    ret
}

maketabDataLarge <- function(tmpdata, group, groupNames, nitem, K, itemloc,
                             Names, itemnames, survey.weights){
    tabdata2 <- tmpdata
    tabdata <- matrix(0L, nrow(tabdata2), sum(K))
    for(i in seq_len(nitem)){
        uniq <- sort(na.omit(unique(tabdata2[,i])))
        if(length(uniq) < K[i]) uniq <- 0L:(K[i]-1L)
        for(j in seq_len(length(uniq)))
            tabdata[,itemloc[i] + j - 1L] <- as.integer(tabdata2[,i] == uniq[j])
    }
    tabdata[is.na(tabdata)] <- 0L
    colnames(tabdata) <- Names
    colnames(tabdata2) <- itemnames
    groupFreq <- vector('list', length(groupNames))
    names(groupFreq) <- groupNames
    for(g in seq_len(length(groupNames))){
        Freq <- as.integer(group == groupNames[g])
        if(!is.null(survey.weights))
            Freq <- Freq * survey.weights
        groupFreq[[g]] <- Freq
    }
    ret <- list(tabdata=tabdata, tabdata2=tabdata2, Freq=groupFreq)
    ret
}

sparseLmat <- function(L, constrain, nconstrain){
    # no constrain
    L <- as.numeric(L)
    whc <- which(L == 1)
    vals <- L[whc]
    full_loc <- cbind(whc, whc)

    # constrain
    for(i in seq_len(length(constrain))){
        cexp <- expand.grid(constrain[[i]], constrain[[i]])
        cexp <- cexp[cexp[,1] != cexp[,2], ]
        full_loc <- rbind(full_loc, as.matrix(cexp))
        vals <- c(vals, rep(1, nrow(cexp)))
    }

    # nconstrain
    for(i in seq_len(length(nconstrain))){
        cexp <- expand.grid(nconstrain[[i]], nconstrain[[i]])
        cexp <- cexp[cexp[,1] != cexp[,2], ]
        full_loc <- rbind(full_loc, as.matrix(cexp))
        vals <- c(vals, c(-2,-2))
    }

    ret <- Matrix::sparseMatrix(i = full_loc[,1], j = full_loc[,2], x=vals,
                                dims = c(length(L), length(L)))
    attr(ret, 'diag') <- L
    ret
}

makeLmats <- function(pars, constrain, random = list(), lrPars = list(), lr.random = list(),
                      nconstrain = NULL){
    ngroups <- length(pars)
    J <- length(pars[[1L]]) - 1L
    LL <- c()
    for(g in seq_len(ngroups))
        for(i in seq_len(J+1L))
            LL <- c(LL, pars[[g]][[i]]@est)
    for(i in seq_len(length(random)))
        LL <- c(LL, random[[i]]@est)
    if(length(lrPars))
        LL <- c(LL, lrPars@est)
    for(i in seq_len(length(lr.random)))
        LL <- c(LL, lr.random[[i]]@est)
    redun_constr <- rep(FALSE, length(LL))
    for(i in seq_len(length(constrain)))
        for(j in 2L:length(constrain[[i]]))
            redun_constr[constrain[[i]][j]] <- TRUE
    if(!is.null(nconstrain)){
        for(i in seq_len(length(nconstrain))){
            stopifnot(length(nconstrain[[i]]) == 2L)
            redun_constr[nconstrain[[i]][2L]] <- TRUE
        }
    }
    L <- sparseLmat(LL, constrain=constrain, nconstrain=nconstrain)
    return(list(L=L, redun_constr=redun_constr))
}

updateGrad <- function(g, L) L %*% g
updateHess <- function(h, L) L %*% h %*% L

makeopts <- function(method = 'MHRM', draws = 2000L, calcLL = TRUE, quadpts = NULL,
                     SE = FALSE, verbose = TRUE, GenRandomPars, dentype = 'Gaussian',
                     SEtol = .001, grsm.block = NULL, D = 1, TOL = NULL,
                     rsm.block = NULL, calcNull = FALSE, BFACTOR = FALSE,
                     technical = list(), hasCustomGroup = FALSE,
                     SE.type = 'Oakes', large = NULL, accelerate = 'Ramsay',
                     optimizer = NULL, solnp_args = list(), nloptr_args = list(),
                     item.Q = NULL, ...)
{
    opts <- list()
    tnames <- names(technical)
    gnames <- c('MAXQUAD', 'NCYCLES', 'BURNIN', 'SEMCYCLES', 'set.seed', 'SEtol', 'symmetric',
                'gain', 'warn', 'message', 'customK', 'customPriorFun', 'customTheta', 'MHcand',
                'parallel', 'NULL.MODEL', 'theta_lim', 'RANDSTART', 'MHDRAWS', 'removeEmptyRows',
                'internal_constraints', 'SEM_window', 'delta', 'MHRM_SE_draws', 'Etable', 'infoAsVcov',
                'PLCI', 'plausible.draws', 'storeEtable', 'keep_vcov_PD', 'Norder', 'MCEM_draws',
                "zeroExtreme", 'mins', 'info_if_converged', 'logLik_if_converged', 'omp', 'nconstrain',
                'standardize_ref', "storeEMhistory", 'fixedEtable')
    if(!all(tnames %in% gnames))
        stop('The following inputs to technical are invalid: ',
             paste0(tnames[!(tnames %in% gnames)], ' '), call.=FALSE)
    if(any(tnames == 'removeEmptyRows'))
        warning(c('removeEmptyRows option has been deprecated. Complete NA response vectors now supported ',
                  'by using NA placeholders'), call.=FALSE)
    if((method %in% c('MHRM', 'MIXED', 'SEM')) && SE.type == 'Oakes') SE.type <- 'MHRM'
    if(method == 'MCEM' && SE && SE.type != 'complete')
        stop('SE.type not currently supported for MCEM method', call.=FALSE)
    if(method == 'QMCEM' && SE && SE.type == 'Fisher')
        stop('Fisher SE.type not supported for QMCEM method', call.=FALSE)
    if((method %in% c('MHRM', 'MIXED', 'SEM')) && !(SE.type %in% c('MHRM', 'FMHRM', 'none')))
        stop('SE.type not supported for stochastic method', call.=FALSE)
    if(!(method %in% c('MHRM', 'MIXED', 'BL', 'EM', 'QMCEM', 'SEM', 'MCEM')))
        stop('method argument not supported', call.=FALSE)
    if(!(SE.type %in% c('Richardson', 'forward', 'central', 'crossprod', 'Louis', 'sandwich',
                        'sandwich.Louis', 'Oakes', 'complete', 'SEM', 'Fisher', 'MHRM', 'FMHRM', 'numerical')))
        stop('SE.type argument not supported', call.=FALSE)
    if(!(method %in% c('EM', 'QMCEM', 'MCEM'))) accelerate <- 'none'
    if(grepl('Davidian', dentype)){
        tmp <- strsplit(dentype, '-')[[1]]
        dentype <- tmp[1L]
        opts$dcIRT_nphi <- as.integer(tmp[2L])
        stopifnot(opts$dcIRT_nphi > 1L)
    }
    if(grepl('mixture', dentype)){
        tmp <- strsplit(dentype, '-')[[1]]
        dentype <- tmp[1L]
        opts$ngroups <- as.integer(tmp[2L])
        stopifnot(opts$n_mixture > 1L)
    }
    if(!(dentype %in% c('Gaussian', 'empiricalhist', 'discrete', 'empiricalhist_Woods', "Davidian",
                        "EH", "EHW", 'mixture')))
        stop('dentype not supported', call.=FALSE)
    opts$method = method
    if(draws < 1) stop('draws must be greater than 0', call.=FALSE)
    opts$draws = draws
    opts$theta_lim = technical$theta_lim
    opts$calcLL = calcLL
    opts$quadpts = quadpts
    opts$SE = SE
    opts$SE.type = SE.type
    opts$verbose = verbose
    opts$SEtol = ifelse(is.null(technical$SEtol), .001, technical$SEtol)
    opts$storeEMhistory = ifelse(is.null(technical$storeEMhistory), FALSE, technical$storeEMhistory)
    opts$grsm.block = grsm.block
    opts$rsm.block = rsm.block
    opts$calcNull = calcNull
    opts$customPriorFun = technical$customPriorFun
    opts$item.Q = item.Q
    if(dentype == "empiricalhist") dentype <- 'EH'
    if(dentype == "empiricalhist_Woods") dentype <- 'EHW'
    opts$dentype <- opts$odentype <- dentype
    opts$zeroExtreme <- FALSE
    if(!is.null(technical$zeroExtreme)) opts$zeroExtreme <- technical$zeroExtreme
    if(BFACTOR) opts$dentype <- 'bfactor'
    if(hasCustomGroup) opts$dentype <- 'custom'
    opts$accelerate = accelerate
    opts$Norder <- ifelse(is.null(technical$Norder), 2L, technical$Norder)
    opts$delta <- ifelse(is.null(technical$delta), 1e-5, technical$delta)
    opts$Etable <- ifelse(is.null(technical$Etable), TRUE, technical$Etable)
    opts$plausible.draws <- ifelse(is.null(technical$plausible.draws), 0, technical$plausible.draws)
    opts$storeEtable <- ifelse(is.null(technical$storeEtable), FALSE, technical$storeEtable)
    if(!is.null(TOL))
        if(is.nan(TOL) || is.na(TOL)) opts$calcNull <- opts$verbose <- FALSE
    opts$TOL <- ifelse(is.null(TOL),
                       if(method %in% c('EM', 'QMCEM', 'MCEM')) 1e-4 else
                           if(method == 'BL') 1e-8 else 1e-3, TOL)
    if(SE.type == 'SEM' && SE){
        opts$accelerate <- 'none'
        if(is.null(TOL)) opts$TOL <- 1e-5
        if(is.null(technical$NCYCLES)) technical$NCYCLES <- 1000L
        if(method == 'QMCEM')
            stop('SEM information matrix not supported with QMCEM estimation', call.=FALSE)
    }
    if(is.null(technical$symmetric)) technical$symmetric <- TRUE
    opts$omp_threads <- ifelse(is.null(technical$omp), .mirtClusterEnv$omp_threads, 1L)
    opts$PLCI <- ifelse(is.null(technical$PLCI), FALSE, technical$PLCI)
    opts$warn <- if(is.null(technical$warn)) TRUE else technical$warn
    opts$message <- if(is.null(technical$message)) TRUE else technical$message
    opts$technical <- technical
    opts$technical$parallel <- ifelse(is.null(technical$parallel), TRUE, technical$parallel)
    opts$technical$omp <- ifelse(is.null(technical$omp), TRUE, technical$omp)
    opts$MAXQUAD <- ifelse(is.null(technical$MAXQUAD), 20000L, technical$MAXQUAD)
    opts$NCYCLES <- ifelse(is.null(technical$NCYCLES), 2000L, technical$NCYCLES)
    if(opts$method %in% c('EM', 'QMCEM', 'MCEM'))
        opts$NCYCLES <- ifelse(is.null(technical$NCYCLES), 500L, technical$NCYCLES)
    opts$BURNIN <- ifelse(is.null(technical$BURNIN), 150L, technical$BURNIN)
    opts$SEMCYCLES <- ifelse(is.null(technical$SEMCYCLES), 100L, technical$SEMCYCLES)
    opts$SEM_from <- ifelse(is.null(technical$SEM_window), 0, technical$SEM_window[1L])
    opts$SEM_to <- ifelse(is.null(technical$SEM_window), 1 - opts$SEtol, technical$SEM_window[2L])
    opts$KDRAWS  <- ifelse(is.null(technical$KDRAWS), 1L, technical$KDRAWS)
    opts$MHDRAWS  <- ifelse(is.null(technical$MHDRAWS), 5L, technical$MHDRAWS)
    opts$MHRM_SE_draws  <- ifelse(is.null(technical$MHRM_SE_draws), 2000L, technical$MHRM_SE_draws)
    opts$internal_constraints  <- ifelse(is.null(technical$internal_constraints),
                                         TRUE, technical$internal_constraints)
    opts$keep_vcov_PD  <- ifelse(is.null(technical$keep_vcov_PD), TRUE, technical$keep_vcov_PD)
    if(dentype == 'mixture'){
        if(opts$method != 'EM')
            stop('Mixture IRT densities only supported when method = \'EM\' ', call.=FALSE)
        if(SE && !(SE.type %in% c('complete', 'Oakes')))
            stop('Only Oakes and complete SE.types current supported for mixture models', call.=FALSE)
    }
    if(dentype %in% c("EH", 'EHW')){
        if(opts$method != 'EM')
            stop('empirical histogram method only applicable when method = \'EM\' ', call.=FALSE)
        if(is.null(opts$quadpts)) opts$quadpts <- 121L
        if(is.null(opts$technical$NCYCLES)) opts$NCYCLES <- 2000L
    }
    if(dentype == 'Davidian'){
        if(opts$method != 'EM')
            stop('Davidian curve method only applicable when method = \'EM\' ', call.=FALSE)
        if(is.null(opts$quadpts)) opts$quadpts <- 121L
    }
    if(is.null(opts$theta_lim)) opts$theta_lim <- c(-6,6)
    if(method == 'QMCEM' && is.null(opts$quadpts)) opts$quadpts <- 5000L
    if((opts$method %in% c('MHRM', 'MIXED') || SE.type == 'MHRM') && !GenRandomPars &&
       opts$plausible.draws == 0L)
        set.seed(12345L)
    if(!is.null(technical$set.seed)) set.seed(technical$set.seed)
    opts$gain <- c(0.1, 0.75)
    if(!is.null(technical$gain)){
        if(length(technical$gain) == 2L && is.numeric(technical$gain))
            opts$gain <- technical$gain
    }
    opts$NULL.MODEL <- ifelse(is.null(technical$NULL.MODEL), FALSE, TRUE)
    opts$USEEM <- ifelse(method == 'EM', TRUE, FALSE)
    opts$returnPrepList <- FALSE
    opts$PrepList <- NULL
    if(is.null(optimizer)){
        opts$Moptim <- if(method %in% c('EM','BL','QMCEM', 'MCEM')) 'BFGS' else 'NR1'
    } else {
        opts$Moptim <- optimizer
    }
    if(method == 'MIXED' && opts$Moptim != 'NR1')
        stop('optimizer currently cannot be changed for mixedmirt()', call.=FALSE)
    if(opts$Moptim == 'solnp'){
        if(is.null(solnp_args$control)) solnp_args$control <- list()
        if(is.null(solnp_args$control$trace)) solnp_args$control$trace <- 0
        if(!method %in% c('EM', 'QMCEM', 'MCEM'))
            stop('solnp only supported for optimization with EM estimation engine',
                                call.=FALSE)
        opts$solnp_args <- solnp_args
    } else if(opts$Moptim == 'nloptr'){
        if(!method %in% c('EM', 'QMCEM', 'MCEM'))
            stop('nloptr only supported for optimization with EM estimation engine',
                                call.=FALSE)
        opts$solnp_args <- nloptr_args
    }
    if(SE && opts$Moptim %in% c('solnp', 'nloptr')) #TODO
        stop('SE computations currently not supported for solnp or nloptr optimizers', call. = FALSE)
    if(!is.null(large)){
        if(is.logical(large))
            if(large) opts$returnPrepList <- TRUE
        if(is.list(large)) opts$PrepList <- large
    }
    if(!is.null(technical$customK)) opts$calcNull <- FALSE
    opts$logLik_if_converged <- ifelse(is.null(technical$logLik_if_converged), TRUE,
                                       technical$logLik_if_converged)
    opts$info_if_converged <- ifelse(is.null(technical$info_if_converged), TRUE,
                                       technical$info_if_converged)
    if(method == 'MCEM'){
        opts$accelerate <- 'none'
        opts$MCEM_draws <- if(is.null(technical$MCEM_draws))
            function(cycles) 500 + (cycles - 1)*2
        else technical$MCEM_draws

    }
    return(opts)
}

reloadPars <- function(longpars, pars, ngroups, J){
    nclasspars <- if(ngroups > 1L)
        do.call(rbind, lapply(pars, function(x) attr(x, 'nclasspars')))
    else attr(pars[[1L]], 'nclasspars')
    .Call('reloadPars', longpars, pars, ngroups, J, nclasspars)
}

computeItemtrace <- function(pars, Theta, itemloc, offterm = matrix(0L, 1L, length(itemloc)-1L),
                             CUSTOM.IND, pis = NULL){
    if(is.null(pis)){
        itemtrace <- .Call('computeItemTrace', pars, Theta, itemloc, offterm)
        if(length(CUSTOM.IND)){
            for(i in CUSTOM.IND){
                Thetas <- Theta
                if(pars[[i]]@nfixedeffects > 0 && nrow(pars[[i]]@fixed.design) == 1)
                    Thetas <- cbind(pars[[i]]@fixed.design[rep(1,nrow(Theta)),], Theta)
                itemtrace[,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(pars[[i]], Theta=Thetas)
            }
        }
    } else {
        tmp_itemtrace <- vector('list', length(pis))
        for(g in seq_len(length(pis))){
            tmp_itemtrace[[g]] <- .Call('computeItemTrace', pars[[g]]@ParObjects$pars, Theta, itemloc, offterm)
            if(length(CUSTOM.IND)){
                for(i in CUSTOM.IND){
                    Thetas <- Theta
                    if(pars[[i]]@nfixedeffects > 0 && nrow(pars[[i]]@fixed.design) == 1)
                        Thetas <- cbind(pars[[i]]@fixed.design[rep(1,nrow(Theta)),], Theta)
                    tmp_itemtrace[[g]][,itemloc[i]:(itemloc[i+1L] - 1L)] <-
                        ProbTrace(pars[[g]]@ParObjects$pars[[i]], Theta=Thetas)
                }
            }
        }
        itemtrace <- do.call(rbind, tmp_itemtrace)
    }
    return(itemtrace)
}

assignItemtrace <- function(pars, itemtrace, itemloc){
    for(i in seq_len(length(pars)-1L))
        pars[[i]]@itemtrace <- itemtrace[ ,itemloc[i]:(itemloc[i+1L] - 1L)]
    pars[[length(pars)]]@itemtrace <- itemtrace
    pars
}

loadESTIMATEinfo <- function(info, ESTIMATE, constrain, warn){
    longpars <- ESTIMATE$longpars
    pars <- ESTIMATE$pars
    ngroups <- length(pars)
    J <- length(pars[[1L]]) - 1L
    info <- nameInfoMatrix(info=info, correction=ESTIMATE$correction, L=ESTIMATE$L,
                           npars=length(longpars))
    ESTIMATE$info <- info
    isna <- is.na(diag(info))
    info <- info[!isna, !isna]
    acov <- try(solve(info), TRUE)
    if(is(acov, 'try-error')){
        if(warn)
            warning('Could not invert information matrix; model may not be (empirically) identified.',
                    call.=FALSE)
        ESTIMATE$fail_invert_info <- TRUE
        return(ESTIMATE)
    } else ESTIMATE$fail_invert_info <- FALSE
    SEtmp <- diag(solve(info))
    if(any(is.na(SEtmp) | is.nan(SEtmp)) || any(SEtmp < 0)){
        if(warn)
            warning('Could not invert information matrix; model may not be (empirically) identified.',
                    call.=FALSE)
        ESTIMATE$fail_invert_info <- TRUE
        return(ESTIMATE)
    }
    SEtmp <- sqrt(SEtmp)
    SE <- rep(NA, length(longpars))
    SE[ESTIMATE$estindex_unique[!isna]] <- SEtmp
    index <- seq_len(length(longpars))
    for(i in seq_len(length(constrain)))
        SE[index %in% constrain[[i]][-1L]] <- SE[constrain[[i]][1L]]
    ind1 <- 1L
    for(g in seq_len(ngroups)){
        for(i in seq_len(J+1L)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1L
            pars[[g]][[i]]@SEpar <- SE[ind1:ind2]
            ind1 <- ind2 + 1L
        }
    }
    ESTIMATE$pars <- pars
    if(length(ESTIMATE$lrPars))
        ESTIMATE$lrPars@SEpar <- SE[ESTIMATE$lrPars@parnum]
    return(ESTIMATE)
}

make.mixed.design <- function(item.formula, itemdesign, data){
    mixed.design <- NULL
    if(!is.null(itemdesign)){
        stopifnot('itemdesign only supported for dichotmous item tests' =
                      all(apply(data, 2, \(x) length(na.omit(unique(x)))) == 2))
        if(nrow(itemdesign) < ncol(data)){
            has_idesign <- colnames(data) %in% rownames(itemdesign)
            if(!any(has_idesign))
                stop('No rownames in itemdesign match colnames(data)', call.=FALSE)
            dummy <- as.data.frame(matrix(NA, sum(!has_idesign), ncol(itemdesign)))
            colnames(dummy) <- colnames(itemdesign)
            itemdesign <- rbind(dummy, itemdesign)
        } else {
            has_idesign <- rep(TRUE, nrow(itemdesign))
            rownames(itemdesign) <- colnames(data)
        }
        itemdesignold <- itemdesign
        if(is.list(item.formula)){
            mf <- lapply(item.formula, \(x){
                if(length(x) == 3){
                    ghost <- x[[2]]
                    itemdesign[[ghost]] <- 1
                }
                model.frame(x, itemdesign, na.action=NULL)
            })
            mm <- lapply(1:length(mf), \(i){
                ret <- model.matrix(item.formula[[i]], mf[[i]])
                ret[rowSums(is.na(ret)) > 0, ] <- NA
                ret
            })
            names(mm) <- do.call(c, lapply(item.formula,
                                           \(x) if(length(x) == 3) as.character(x[[2]]) else ""))
            for(i in 1:length(mm))
                if(names(mm)[i] != "")
                    colnames(mm[[i]]) <- paste0(names(mm)[i], '.', colnames(mm[[i]]))
            mm <- do.call(cbind, mm)
        } else {
            mf <- model.frame(item.formula, itemdesign, na.action = NULL)
            mm <- model.matrix(item.formula, mf)
        }

        mixed.design <- list(random=NULL, fixed=mm, from='mirt',
                             lr.random=NULL, lr.fixed=NULL, has_idesign=has_idesign)
        attr(mixed.design, 'itemdesign') <- itemdesignold
    }
    mixed.design
}

make.randomdesign <- function(random, longdata, covnames, itemdesign, N, LR=FALSE){
    ret <- vector('list', length(random))
    for(i in seq_len(length(random))){
        f <- gsub(" ", "", as.character(random[[i]])[2L])
        splt <- strsplit(f, '\\|')[[1L]]
        if(any(grepl('\\*', splt[2L]) | grepl('\\+', splt[2L])))
            stop('The + and * operators are not supported. Please specify
                 which effects you want to interact with the : operator, and specify
                 additional random effects in separate list elements', call.=FALSE)
        gframe <- model.frame(as.formula(paste0('~',splt[2L])), longdata)
        sframe <- model.frame(as.formula(paste0('~',splt[1L])), longdata)
        levels <- interaction(gframe)
        uniq_levels <- unique(levels)
        matpar <- diag(1L + ncol(sframe))
        estmat <- lower.tri(matpar, diag=TRUE)
        ndim <- ncol(matpar)
        if(strsplit(f, '+')[[1L]][[1L]] == '-')
            estmat[lower.tri(estmat)] <- FALSE
        fn <- paste0('COV_', c(splt[2L], colnames(sframe)))
        FNCOV <- outer(fn, c(splt[2L], colnames(sframe)), FUN=paste, sep='_')
        par <- matpar[lower.tri(matpar, diag=TRUE)]
        est <- estmat[lower.tri(estmat, diag=TRUE)]
        names(par) <- names(est) <- FNCOV[lower.tri(FNCOV, diag=TRUE)]
        drawvals <- matrix(0, length(uniq_levels), ndim,
                           dimnames=list(uniq_levels, NULL))
        mtch <- match(levels, rownames(drawvals))
        gdesign <- matrix(1, length(levels), 1L, dimnames = list(NULL, splt[2L]))
        if(ncol(sframe) != 0L){
            if(grepl('-1+', splt[1L])){
                splt[1L] <- strsplit(splt[1L], '-1\\+')[[1]][2]
            } else if(grepl('0+', splt[1L]))
                splt[1L] <- strsplit(splt[1L], '0\\+')[[1]][2]
            gdesign <- cbind(gdesign,
                             model.matrix(as.formula(paste0('~',splt[1L])), sframe)[,-1L,drop=FALSE])
        }
        tmp <- matrix(-Inf, ndim, ndim)
        diag(tmp) <- 1e-4
        lbound <- tmp[lower.tri(tmp, diag=TRUE)]
        ret[[i]] <- new('RandomPars',
                        par=par,
                        est=est,
                        parnames=names(est),
                        SEpar=rep(NaN,length(par)),
                        ndim=ndim,
                        lbound=lbound,
                        ubound=rep(Inf, length(par)),
                        gframe=gframe,
                        gdesign=gdesign,
                        cand.t.var=.5,
                        any.prior=FALSE,
                        prior.type=rep(0L, length(par)),
                        prior_1=rep(NaN,length(par)),
                        prior_2=rep(NaN,length(par)),
                        drawvals=drawvals,
                        mtch=mtch)
    }
    ret
}

make.lrdesign <- function(df, formula, factorNames, EM=FALSE, TOL){
    nfact <- length(factorNames)
    if(is.list(formula)){
        if(!all(names(formula) %in% factorNames))
            stop('List of fixed effect names do not match factor names', call.=FALSE)
        estnames <- X <- vector('list', length(formula))
        for(i in 1L:length(formula)){
            X[[i]] <- model.matrix(as.formula(formula[[i]]), df)
            estnames[[i]] <- colnames(X[[i]])
        }
        X <- do.call(cbind, X)
        X <- X[,unique(colnames(X))]
    } else {
        X <- model.matrix(as.formula(formula), df)
    }
    tXX <- t(X) %*% X
    qr_XX <- qr(0)
    if(!is.na(TOL)){
        qr_XX <- try(qr(tXX), silent = TRUE)
        if(!is.nan(TOL)){
            if(is(qr_XX, 'try-error'))
                stop('Latent regression design matrix contains problematic terms.', call. = FALSE)
        }
    }
    beta <- matrix(0, ncol(X), nfact)
    sigma <- matrix(0, nfact, nfact)
    diag(sigma) <- 1
    if(is.list(formula)){
        est <- matrix(FALSE, nrow(beta), ncol(beta))
        for(i in 1L:length(formula)){
            name <- names(formula)[[i]]
            pick <- which(name == factorNames)
            est[colnames(X) %in% estnames[[i]], pick] <- TRUE
        }
    } else est <- matrix(TRUE, nrow(beta), ncol(beta))
    est[1L, ] <- FALSE
    est <- as.logical(est)
    names(est) <- as.character(t(outer(factorNames, colnames(X),
                                     FUN = function(X, Y) paste(X,Y,sep="_"))))
    colnames(beta) <- factorNames
    rownames(beta) <- colnames(X)
    par <- as.numeric(beta)
    names(par) <- names(est)
    ret <- new('lrPars',
               par=par,
               SEpar=rep(NaN,length(par)),
               est=est,
               parnames=names(est),
               beta=beta,
               sigma=sigma,
               nfact=nfact,
               nfixed=ncol(X),
               df=df,
               X=as.matrix(X),
               tXX=tXX,
               qr_XX=list(qr = qr_XX),
               lbound=rep(-Inf,length(par)),
               ubound=rep(Inf,length(par)),
               any.prior=FALSE,
               prior.type=rep(0L, length(par)),
               prior_1=rep(NaN,length(par)),
               prior_2=rep(NaN,length(par)),
               formula=if(!is.list(formula)) list(formula) else formula,
               EM=EM)
    ret
}

update.lrPars <- function(df, lrPars){
    pick <- df$class == 'lrPars'
    df2 <- df[pick, , drop=FALSE]
    lrPars@est[] <- df2$est
    lrPars@par <- lrPars@beta[] <- df2$value
    if(!all(df2$lbound == -Inf))
        warning('latent regression parameters cannot be bounded. Ignoring constraint', call.=FALSE)
    if(!all(df2$ubound == Inf))
        warning('latent regression parameters cannot be bounded. Ignoring constraint', call.=FALSE)
    if(!all(df2$prior.type == 'none'))
        warning('latent regression parameters do not support prior distribution. Ignoring input.',
                call.=FALSE)
    lrPars
}


OffTerm <- function(random, J, N){
    ret <- numeric(N*J)
    for(i in seq_len(length(random))){
        tmp <- rowSums(random[[i]]@gdesign*random[[i]]@drawvals[random[[i]]@mtch, ,drop=FALSE])
        ret <- ret + tmp
    }
    return(matrix(ret, N, J))
}

reloadRandom <- function(random, longpars){
    for(i in seq_len(length(random))){
        parnum <- random[[i]]@parnum
        random[[i]]@par <- longpars[min(parnum):max(parnum)]
    }
    random
}

phi_bound <- function(phi, bound = c(-pi/2, pi/2)){
    for(i in 1L:length(phi)){
        while(TRUE){
            if(phi[i] > bound[1L] && phi[i] <= bound[2L]) break
            if(phi[i] < bound[1L]){
                phi[i] <- phi[i] + pi
            }
            if(phi[i] >= bound[2L]){
                phi[i] <- phi[i] - pi
            }
        }
    }
    phi
}

smooth.cor <- function(x){
    eig <- eigen(x)
    negvalues <- eig$values <= 0
    while (any(negvalues)) {
        eig2 <- ifelse(eig$values < 0, 100 * .Machine$double.eps, eig$values)
        x <- eig$vectors %*% diag(eig2) %*% t(eig$vectors)
        x <- x/sqrt(diag(x) %*% t(diag(x)))
        eig <- eigen(x)
        negvalues <- eig$values <= .Machine$double.eps
    }
    x
}

smooth.cov <- function(cov){
    ev <- eigen(cov)
    v <- ev$values
    off <- sum(v) * .01
    v[v < 0] <- off
    v <- sum(ev$values) * v/sum(v)
    v <- ifelse(v < 1e-4, 1e-4, v)
    return(ev$vectors %*% diag(v) %*% t(ev$vectors))
}

RMSEA.CI <- function(X2, df, N, ci.lower=.05, ci.upper=.95) {

    lower.lambda <- function(lambda) pchisq(X2, df=df, ncp=lambda) - ci.upper
    upper.lambda <- function(lambda) pchisq(X2, df=df, ncp=lambda) - ci.lower

    lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=X2)$root, silent=TRUE)
    lambda.u <- try(uniroot(f=upper.lambda, lower=0, upper=max(N, X2*5))$root, silent=TRUE)
    if(!is(lambda.l, 'try-error')){
        RMSEA.lower <- sqrt(lambda.l/(N*df))
    } else {
        RMSEA.lower <- 0
    }
    if(!is(lambda.u, 'try-error')){
        RMSEA.upper <- sqrt(lambda.u/(N*df))
    } else {
        RMSEA.upper <- 0
    }

    return(c(RMSEA.lower, RMSEA.upper))
}

longpars_constrain <- function(longpars, constrain, nconstrain = NULL){
    for(i in seq_len(length(constrain)))
        longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
    if(!is.null(nconstrain)){
        for(i in seq_len(length(nconstrain)))
            longpars[nconstrain[[i]][-1L]] <- -longpars[nconstrain[[i]][1L]]
    }
    longpars
}

BL.LL <- function(p, est, longpars, pars, ngroups, J, Theta, PrepList, specific, sitems, nconstrain,
               CUSTOM.IND, EHPrior, Data, dentype, itemloc, theta, constrain, lrPars, omp_threads){
    #TODO use nconstrain
    longpars[est] <- p
    longpars <- longpars_constrain(longpars=longpars, constrain=constrain)
    pars2 <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    gstructgrouppars <- prior <- Prior <- vector('list', ngroups)
    full <- length(lrPars) > 0L
    if(full){
        lrPars@par <- longpars[lrPars@parnum]
        lrPars@beta[] <- lrPars@par
        lrPars@mus <- lrPars@X %*% lrPars@beta
    }
    if(dentype %in% c('EH', 'EHW')){
        Prior[[1L]] <- EHPrior[[1L]]
    } else if(dentype == 'custom'){
        for(g in seq_len(ngroups)){
            gp <- pars[[g]][[J+1L]]
            Prior[[g]] <- gp@den(gp, Theta)
            Prior[[g]] <- Prior[[g]] / sum(Prior[[g]])
        }
    } else if(dentype == 'discrete'){
        for(g in seq_len(ngroups)){
            gp <- pars[[g]][[J+1L]]
            if(full){
                Prior[[g]] <- gp@den(gp, Theta, mus=lrPars@mus)
                Prior[[g]] <- Prior[[g]]/rowSums(Prior[[g]])
            } else {
                Prior[[g]] <- gp@den(gp, Theta)
                Prior[[g]] <- Prior[[g]] / sum(Prior[[g]])
            }
        }
    } else {
        for(g in seq_len(ngroups)){
            gstructgrouppars[[g]] <- ExtractGroupPars(pars2[[g]][[J+1L]])
            if(dentype == 'bfactor'){
                prior[[g]] <- dnorm(theta, 0, 1)
                prior[[g]] <- prior[[g]]/sum(prior[[g]])
                Prior[[g]] <- apply(expand.grid(prior[[g]], prior[[g]]), 1L, prod)
                next
            }
            if(full){
                Prior[[g]] <- mirt_dmvnorm(Theta[ ,1L:ncol(lrPars@mus),drop=FALSE],
                                           lrPars@mus, gstructgrouppars[[g]]$gcov, quad=TRUE)
                Prior[[g]] <- Prior[[g]]/rowSums(Prior[[g]])
            } else {
                Prior[[g]] <- mirt_dmvnorm(Theta,gstructgrouppars[[g]]$gmeans,
                                           gstructgrouppars[[g]]$gcov)
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
            }
        }
    }
    LL <- 0
    for(g in seq_len(ngroups)){
        expected <- Estep.mirt(pars=pars2[[g]],
                               tabdata=Data$tabdatalong, wmiss=Data$wmiss,
                               freq=if(full) rep(1L, nrow(Prior[[1L]])) else Data$Freq[[g]],
                               Theta=Theta, prior=Prior[[g]], itemloc=itemloc,
                               CUSTOM.IND=CUSTOM.IND, full=full, Etable=FALSE, omp_threads=omp_threads)$expected
        LL <- LL + sum(Data$Freq[[g]] * log(expected), na.rm = TRUE)
    }
    LL
}

select_quadpts <- function(nfact) switch(as.character(nfact),
                                         '1'=61, '2'=31, '3'=15, '4'=9, '5'=7, 3)

mirt_rmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                         check = FALSE, pre.ev=list())
{
    if(!length(pre.ev)){
        # Version modified from mvtnorm::rmvnorm, version 0.9-9996, 19-April, 2014.
        if(check){
            if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE))
                stop("sigma must be a symmetric matrix", call.=FALSE)
            if (length(mean) != nrow(sigma))
                stop("mean and sigma have non-conforming size", call.=FALSE)
        }
        ev <- eigen(sigma, symmetric = TRUE)
        NCOL <- ncol(sigma)
        if(check)
            if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1])))
                warning("sigma is numerically non-positive definite", call.=FALSE)
    } else {
        ev <- pre.ev
        NCOL <- length(ev$values)
    }
    retval <- ev$vectors %*% diag(sqrt(ev$values), NCOL) %*% t(ev$vectors)
    retval <- matrix(rnorm(n * NCOL), nrow = n) %*%  retval
    retval <- sweep(retval, 2L, mean, "+")
    colnames(retval) <- names(mean)
    retval
}

mirt_dmvnorm <- function(x, mean, sigma, log = FALSE, quad = FALSE, stable = TRUE, ...)
{
    if(quad && is.matrix(mean)){
        isigma <- solve(sigma)
        distval <- matrix(0, nrow(mean), nrow(x))
        for(i in seq_len(nrow(mean))){
            centered <- t(t(x) - mean[i,])
            distval[i, ] <- rowSums((centered %*% isigma) * centered)
        }
    } else {
        if(is.matrix(mean)){
            centered <- x - mean
            distval <- rowSums((centered %*% solve(sigma)) * centered)
        } else {
            distval <- mahalanobis(x, center = mean, cov = sigma)
        }
    }
    eigs <- eigen(sigma, symmetric=TRUE, only.values=TRUE)$values
    if(any(eigs < 0))
        stop('sigma matrix contains negative eigenvalues', call. = FALSE)
    logdet <- sum(log(eigs))
    logretval <- -(ncol(x)*log(2*pi) + logdet + distval)/2
    if(stable)
        logretval <- ifelse(logretval < -690.7755, -690.7755, logretval)
    if(log) return(logretval)
    exp(logretval)

}

# prior for latent class analysis
lca_prior <- function(Theta, Etable){
    TP <- nrow(Theta)
    if ( is.null(Etable) ){
        prior <- rep( 1/TP , TP )
    } else {
        prior <- rowSums(Etable)
    }
    prior <- prior / sum(prior)
    return(prior)
}

makeObstables <- function(dat, K, which.items){
    ret <- vector('list', ncol(dat))
    sumscore <- rowSums(dat)
    for(i in seq_len(length(ret))){
        if(!(i %in% which.items)) next
        ret[[i]] <- matrix(0, sum(K-1L)+1L, K[i])
        colnames(ret[[i]]) <- paste0(1L:K[i]-1L)
        rownames(ret[[i]]) <- paste0(1L:nrow(ret[[i]])-1L)
        split <- by(sumscore, dat[,i], table)
        for(j in seq_len(length(split))){
            m <- match(names(split[[j]]), rownames(ret[[i]]))
            ret[[i]][m,j] <- split[[j]]
        }
        ret[[i]] <- ret[[i]][-c(1L, nrow(ret[[i]])), ]
    }
    ret
}

collapseCells <- function(O, E, mincell = 1){
    for(i in seq_len(length(O))){
        On <- O[[i]]
        En <- E[[i]]
        if(is.null(En)) next
        drop <- which(rowSums(is.na(En)) > 0)
        En[is.na(En)] <- 0

        #collapse known upper and lower sparce cells
        if(length(drop)){
            up <- drop[1L]:drop[length(drop)/2]
            low <- drop[length(drop)/2 + 1L]:drop[length(drop)]
            En[max(up)+1, ] <- colSums(En[c(up, max(up)+1), , drop = FALSE])
            On[max(up)+1, ] <- colSums(On[c(up, max(up)+1), , drop = FALSE])
            En[min(low)-1, ] <- colSums(En[c(low, min(low)-1), , drop = FALSE])
            On[min(low)-1, ] <- colSums(On[c(low, min(low)-1), , drop = FALSE])
            En[c(up, low), ] <- On[c(up, low), ] <- NA
            En <- na.omit(En)
            On <- na.omit(On)
        }

        if(nrow(On) == 0){
            E[[i]] <- En
            O[[i]] <- On
            break
        }

        #drop 0's and 1's
        drop <- rowSums(On) == 0L
        On <- On[!drop,]
        En <- En[!drop,]
        L <- En < mincell
        drop <- c()
        for(j in seq_len(nrow(On)-1L)){
            ss <- sum(On[j,])
            if(ss == 1L){
                drop <- c(drop, j)
                On[j+1L, ] <- On[j+1L, ] + On[j, ]
                En[j+1L, ] <- En[j+1L, ] + En[j, ]
            }
        }
        if(length(drop)){
            On <- On[-drop,]
            En <- En[-drop,]
        }
        ss <- sum(On[nrow(On),])
        if(ss == 1L){
            On[nrow(On)-1L, ] <- On[nrow(On)-1L, ] + On[nrow(On), ]
            En[nrow(On)-1L, ] <- En[nrow(On)-1L, ] + En[nrow(On), ]
            On <- On[-nrow(On),]; En <- En[-nrow(En),]
        }

        #collapse accross as much as possible
        if(ncol(En) > 2L){
            for(j in seq_len(nrow(En))){
                if(!any(L[j,])) next
                tmp <- En[j, ]
                tmp2 <- On[j, ]
                while(length(tmp) > 2L){
                    m <- min(tmp)
                    whc <- max(which(m == tmp))
                    if(whc == 1L){
                        tmp[2L] <- tmp[2L] + tmp[1L]
                        tmp2[2L] <- tmp2[2L] + tmp2[1L]
                    } else if(whc == length(tmp)){
                        tmp[length(tmp)-1L] <- tmp[length(tmp)-1L] + tmp[length(tmp)]
                        tmp2[length(tmp2)-1L] <- tmp2[length(tmp2)-1L] + tmp2[length(tmp2)]
                    } else {
                        left <- min(tmp[whc-1L], tmp[whc+1L]) == c(tmp[whc-1L], tmp[whc+1L])[1L]
                        pick <- if(left) whc-1L else whc+1L
                        tmp[pick] <- tmp[pick] + tmp[whc]
                        tmp2[pick] <- tmp2[pick] + tmp2[whc]
                    }
                    tmp[whc] <- tmp2[whc] <- NA
                    tmp <- na.omit(tmp); tmp2 <- na.omit(tmp2)
                    if(all(tmp >= mincell)) break
                }
                tmp <- c(tmp, rep(NA, ncol(En)-length(tmp)))
                tmp2 <- c(tmp2, numeric(ncol(En)-length(tmp2)))
                En[j, ] <- tmp
                On[j, ] <- tmp2
            }
        }
        En[is.na(En)] <- 0
        # collapse right columns if they are too rare
        if(ncol(En) > 2L){
            while(TRUE){
                pick <- colSums(En) < mincell * ceiling(nrow(En) * .1)
                if(!pick[length(pick)] || ncol(En) == 2L) break
                if(pick[length(pick)]){
                    On[ ,length(pick)-1L] <- On[ ,length(pick)-1L] + On[ ,length(pick)]
                    En[ ,length(pick)-1L] <- En[ ,length(pick)-1L] + En[ ,length(pick)]
                    On <- On[ ,-length(pick)]; En <- En[ ,-length(pick)]
                }
            }
        }
        dropcol <- logical(ncol(En))
        #drop all other columns if they are very rare
        for(j in ncol(En):2L){
            tmp <- sum(En[,j] > 0) / nrow(En)
            if(tmp < .05){
                dropcol[j] <- TRUE
                En[,j-1L] <- En[,j-1L] + En[,j]
                On[,j-1L] <- On[,j-1L] + On[,j]
            }
        }
        En <- En[,!dropcol]; On <- On[,!dropcol]

        #merge across
        L <- En < mincell & En != 0
        while(any(L, na.rm = TRUE)){
            if(!is.matrix(L)) break
            whc <- min(which(rowSums(L) > 0L))
            if(whc == 1L){
                En[2L,] <- En[2L, ] + En[1L,]
                On[2L,] <- On[2L, ] + On[1L,]
                En <- En[-1L,]; On <- On[-1L,]
            } else if(whc == nrow(En)){
                En[nrow(En)-1L,] <- En[nrow(En)-1L, ] + En[nrow(En),]
                On[nrow(On)-1L,] <- On[nrow(On)-1L, ] + On[nrow(On),]
                En <- En[-nrow(En),]; On <- On[-nrow(On),]
            } else {
                ss <- c(sum(On[whc-1L,]), sum(On[whc+1L,]))
                up <- (min(ss) == ss)[1L]
                pick <- if(up) whc-1L else whc+1L
                En[pick,] <- En[pick, ] + En[whc,]
                On[pick,] <- On[pick, ] + On[whc,]
                En <- En[-whc,]; On <- On[-whc,]
            }
            L <- En < mincell & En != 0
        }
        En[En == 0] <- NA
        E[[i]] <- En
        O[[i]] <- On
    }
    return(list(O=O, E=E))
}

MGC2SC <- function(x, which){
    tmp <- x@ParObjects$pars[[which]]
    tmp@Model$lrPars <- x@ParObjects$lrPars
    ind <- 1L
    for(i in seq_len(x@Data$nitems) + 1L){
        tmp@ParObjects$pars[[i]]@parnum[] <- seq(ind, ind + length(tmp@ParObjects$pars[[i]]@parnum) - 1L)
        ind <- ind + length(tmp@ParObjects$pars[[i]]@parnum)
    }
    tmp@Data <- x@Data
    tmp@Data$completely_missing <- integer(0L)
    tmp@Data$data <- tmp@Data$data[tmp@Data$group == tmp@Data$groupNames[which], , drop=FALSE]
    tmp@Data$rowID <- 1L:nrow(tmp@Data$data)
    tmp@Data$Freq[[1L]] <- tmp@Data$Freq[[which]]
    tmp@Data$fulldata[[1L]] <- x@Data$fulldata[[which]]
    tmp@Data$ngroups <- 1L
    tmp@Model$model <- x@Model$model
    tmp
}

computeNullModel <- function(data, key, group=NULL){
    if(!is.null(group) && !all(group == 'all')){
        null.mod <- suppressMessages(multipleGroup(data, 1L, group=group, verbose=FALSE,
                                  key=key, quadpts=3, technical=list(NULL.MODEL=TRUE)))
    } else {
        null.mod <- suppressMessages(mirt(data, 1L, verbose=FALSE,
                                          key=key, quadpts=3, technical=list(NULL.MODEL=TRUE)))
    }
    null.mod
}

loadSplineParsItem <- function(x, Theta){
    sargs <- x@sargs
    Theta_prime <- if(x@stype == 'bs'){
        splines::bs(Theta, df=sargs$df, knots=sargs$knots,
                    degree=sargs$degree, intercept=sargs$intercept)
    } else if(x@stype == 'ns'){
        splines::ns(Theta, df=sargs$df, knots=sargs$knots,
                    intercept=sargs$intercept)
    }
    class(Theta_prime) <- 'matrix'
    x@Theta_prime <- Theta_prime
    x
}

loadSplinePars <- function(pars, Theta, MG = TRUE){
    fn <- function(pars, Theta){
        cls <- sapply(pars, class)
        pick <- which(cls == 'spline')
        if(length(pick)){
            for(i in pick)
                pars[[i]] <- loadSplineParsItem(pars[[i]], Theta)
        }
        return(pars)
    }
    if(MG){
        for(g in seq_len(length(pars))){
            pars[[g]] <- fn(pars[[g]], Theta)
        }
    } else {
        pars <- fn(pars, Theta)
    }
    return(pars)
}

latentRegression_obj <- function(data, covdata, formula, dentype, method){
    if(!is.null(covdata) && !is.null(formula)){
        if(!dentype %in% c("Gaussian", 'discrete'))
            stop('Only Guassian dentype currently supported for latent regression models',
                 call.=FALSE)
        if(!is.data.frame(covdata))
            stop('covdata must be a data.frame object', call.=FALSE)
        if(nrow(covdata) != nrow(data))
            stop('number of rows in covdata do not match number of rows in data', call.=FALSE)
        if(!(method %in% c('EM', 'QMCEM')))
            stop('method must be from the EM estimation family', call.=FALSE)
        tmp <- apply(covdata, 1, function(x) sum(is.na(x)) > 0)
        if(any(tmp)){
            message('removing rows with NAs in covdata')
            covdata <- covdata[!tmp, ]
            data <- data[!tmp, ]
        }
        completely_missing <- which(rowSums(is.na(data)) == ncol(data))
        if(length(completely_missing))
            covdata <- covdata[-completely_missing, , drop=FALSE]
        latent.regression <- list(df=covdata, formula=formula, data=data, EM=TRUE)
    } else latent.regression <- NULL
    latent.regression
}

bs_range <- function(x, CI){
    ss <- sort(x)
    N <- length(ss)
    ret <- c(lower = ss[floor(N * (1-CI)/2)],
             upper = ss[ceiling(N * (1 - (1-CI)/2))])
    ret
}

# borrowed and modified from emdbook package, March 1 2017
mixX2 <- function (p, df = 1, mix = 0.5, lower.tail = TRUE)
{
    df <- rep(df, length.out = length(p))
    mix <- rep(mix, length.out = length(p))
    c1 <- ifelse(df == 1, if (lower.tail)  1
                 else 0, pchisq(p, df - 1, lower.tail = lower.tail))
    c2 <- pchisq(p, df, lower.tail = lower.tail)
    r <- mix * c1 + (1 - mix) * c2
    r
}

# borrowed and modified from car package, July 31, 2020
makeHypothesis <- function (cnames, hypothesis, rhs = NULL)
{
    parseTerms <- function(terms) {
        component <- gsub("^[-\\ 0-9\\.]+", "", terms)
        component <- gsub(" ", "", component, fixed = TRUE)
        component
    }
    stripchars <- function(x) {
        x <- gsub("\\n", " ", x)
        x <- gsub("\\t", " ", x)
        x <- gsub(" ", "", x, fixed = TRUE)
        x <- gsub("*", "", x, fixed = TRUE)
        x <- gsub("-", "+-", x, fixed = TRUE)
        x <- strsplit(x, "+", fixed = TRUE)[[1]]
        x <- x[x != ""]
        x
    }
    char2num <- function(x) {
        x[x == ""] <- "1"
        x[x == "-"] <- "-1"
        as.numeric(x)
    }
    constants <- function(x, y) {
        with.coef <- unique(unlist(sapply(y, function(z) which(z ==
                                                                   parseTerms(x)))))
        if (length(with.coef) > 0)
            x <- x[-with.coef]
        x <- if (is.null(x))
            0
        else sum(as.numeric(x))
        if (any(is.na(x)))
            stop("The hypothesis \"", hypothesis,
                 "\" is not well formed: contains bad coefficient/variable names.", call.=FALSE)
        x
    }
    coefvector <- function(x, y) {
        rv <- gsub(" ", "", x, fixed = TRUE) == parseTerms(y)
        if (!any(rv))
            return(0)
        if (sum(rv) > 1)
            stop("The hypothesis \"", hypothesis, "\" is not well formed.")
        rv <- sum(char2num(unlist(strsplit(y[rv], x, fixed = TRUE))))
        if (is.na(rv))
            stop("The hypothesis \"", hypothesis,
                 "\" is not well formed: contains non-numeric coefficients.", call.=FALSE)
        rv
    }
    if (!is.null(rhs))
        rhs <- rep(rhs, length.out = length(hypothesis))
    if (length(hypothesis) > 1)
        return(rbind(Recall(cnames, hypothesis[1], rhs[1]), Recall(cnames,
                                                                   hypothesis[-1], rhs[-1])))
    cnames_symb <- sapply(c("@", "#", "~"),
                          function(x) length(grep(x, cnames)) < 1)
    if (any(cnames_symb)) {
        cnames_symb <- head(c("@", "#", "~")[cnames_symb],
                            1)
        cnames_symb <- paste(cnames_symb, seq_along(cnames),
                             cnames_symb, sep = "")
        hypothesis_symb <- hypothesis
        for (i in order(nchar(cnames), decreasing = TRUE)) hypothesis_symb <- gsub(cnames[i],
                                                                                   cnames_symb[i], hypothesis_symb, fixed = TRUE)
    }
    else {
        stop("The hypothesis \"", hypothesis,
             "\" is not well formed: contains non-standard coefficient names.", call.=FALSE)
    }
    lhs <- strsplit(hypothesis_symb, "=", fixed = TRUE)[[1]]
    if (is.null(rhs)) {
        if (length(lhs) < 2)
            rhs <- "0"
        else if (length(lhs) == 2) {
            rhs <- lhs[2]
            lhs <- lhs[1]
        }
        else stop("The hypothesis \"", hypothesis,
                  "\" is not well formed: contains more than one = sign.", call.=FALSE)
    }
    else {
        if (length(lhs) < 2)
            as.character(rhs)
        else stop("The hypothesis \"", hypothesis,
                  "\" is not well formed: contains a = sign although rhs was specified.", call.=FALSE)
    }
    lhs <- stripchars(lhs)
    rhs <- stripchars(rhs)
    rval <- sapply(cnames_symb, coefvector, y = lhs) - sapply(cnames_symb,
                                                              coefvector, y = rhs)
    rval <- c(rval, constants(rhs, cnames_symb) - constants(lhs,
                                                            cnames_symb))
    names(rval) <- c(cnames, "*rhs*")
    if (is.null(dim(rval)))
        rval <- t(rval)
    rval
}

get_deriv_coefs <- function(order, deriv = 1L){
    if(deriv == 1L){
        ret <- switch(as.character(order),
                      "1" = c(-1, 1),
                      "2" = c(-1/2, 1/2))
    }
    ret
}

cfi <- function(X2, X2.null, df, df.null){
    ret <- 1 - (X2 - df) / (X2.null - df.null)
    if(ret > 1) ret <- 1
    if(ret < 0) ret <- 0
    ret
}

tli <- function(X2, X2.null, df, df.null)
    (X2.null/df.null - X2/df) / (X2.null/df.null - 1)



rmsea <- function(X2, df, N){
    ret <- suppressWarnings(ifelse( (X2/df - 1) > 0 & df > 0,
                   sqrt((X2/df - 1) / (N-1)), 0))
    ret <- ifelse(df <= 0, NaN, ret)
    ret
}

removeMissing <- function(obj){
    dat <- extract.mirt(obj, 'data')
    obj@Data$data <- na.omit(dat)
    pick <- attr(obj@Data$data, 'na.action')
    if(is.null(pick)){
        message('Data does not contain missing values. Continuing normally')
        return(obj)
    }
    for(g in seq_len(length(obj@Data$groupNames))){
        whc <- obj@Data$group == obj@Data$groupNames[g]
        ind2 <- obj@Data$rowID[whc]
        pick2 <- ind2 %in% names(pick)
        obj@Data$fulldata[[g]] <- obj@Data$fulldata[[g]][!pick2, , drop=FALSE]
    }
    obj@Data$group <- obj@Data$group[-pick]
    if(is(obj, 'MultipleGroupClass')){
        for(g in seq_len(length(obj@Data$groupNames))){
            obj@ParObjects$pars[[g]]@Data$data <- dat[obj@Data$groupNames[g] == obj@Data$group,
                                                         , drop=FALSE]
        }
    }
    obj
}

controlCandVar <- function(PA, cand, min = .1, max = .6){
    if(PA > max) cand <- cand * 1.05
    else if(PA < min) cand <- cand * 0.9
    if(cand < .001) cand <- .001
    cand
}

toInternalItemtype <- function(itemtype){
    itemtype <- ifelse(itemtype %in% c('2PL', '3PL', '3PLu', '4PL'), 'dich', itemtype)
    itemtype <- ifelse(itemtype %in% c('PC2PL', 'PC3PL'), 'partcomp', itemtype)
    itemtype <- ifelse(itemtype %in% c("2PLNRM", "3PLNRM", "3PLuNRM", "4PLNRM"), 'nestlogit', itemtype)
    itemtype
}

# function borrowed and edited from the sfsmisc package, v1.1-1. Date: 18, Oct, 2017
QUnif <- function (n, min = 0, max = 1, n.min = 1, p, leap = 1, silent = FALSE)
{
    digitsBase <- function (x, base = 2, ndigits = 1 + floor(1e-09 + log(max(x), base)))
    {
        if (any(x < 0))
            stop("'x' must be non-negative integers", call.=FALSE)
        if (any(x != trunc(x)))
            stop("'x' must be integer-valued", call.=FALSE)
        r <- matrix(0, nrow = ndigits, ncol = length(x))
        if (ndigits >= 1)
            for (i in ndigits:1) {
                r[i, ] <- x%%base
                if (i > 1)
                    x <- x%/%base
            }
        class(r) <- "basedInt"
        attr(r, "base") <- base
        r
    }
    sHalton <- function (n.max, n.min = 1, base = 2, leap = 1)
    {
        stopifnot((leap <- as.integer(leap)) >= 1)
        nd <- as.integer(1 + log(n.max, base))
        dB <- digitsBase(if (leap == 1)
            n.min:n.max
            else seq(n.min, n.max, by = leap), base = base, ndigits = nd)
        colSums(dB/base^(nd:1))
    }
    stopifnot(1 <= (n <- as.integer(n)), length(n) == 1, 1 <=
                  (p <- as.integer(p)), length(p) == 1, length(min) ==
                  p || length(min) == 1, length(max) == p || length(max) ==
                  1, 1 <= (n.min <- as.integer(n.min)), 1 <= (leap <- as.integer(leap)),
              (n.max <- n.min + (n - 1:1) * leap) < .Machine$integer.max)
    pr. <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
             43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
             103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157,
             163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
             227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277,
             281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
             353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
             421, 431, 433, 439, 443, 449, 457)
    pr <- pr.[1:p]
    if (leap > 1 && any(leap == pr) && length(pr.) >= p + 1)
        pr <- c(pr[leap != pr], pr.[p + 1])
    max <- rep.int(max, p)
    min <- rep.int(min, p)
    dU <- max - min
    r <- matrix(0, n, p)
    for (j in 1:p) r[, j] <- min[j] + dU[j] * sHalton(n.max,
                                                      n.min, base = pr[j], leap = leap)
    r
}

nominal_rescale_d <- function(as, nfact, ncat){
    a <- sum(as[1:nfact])
    ask <- a * as[-(1:nfact)]
    d <- max(as[-(1:nfact)]) / (ncat - 1)
    d
}

CA <- function(dat, guess_correction = rep(0, ncol(dat))){
    n <- ncol(dat)
    C <- cov(dat)
    d <- diag(C)
    dich <- apply(dat, 2, function(x) length(unique(x)) == 2L)
    for(i in seq_len(ncol(dat))){
        if(dich[i] & guess_correction[i] > 0){
            U <- mean(dat[,i]) - min(dat[,i])
            v <- max((U - guess_correction[i]) / (1 - guess_correction[i]), 0)
            d[i] <- v * (1 - v)
        }
    }
    (n/(n - 1)) * (1 - sum(d)/sum(C))
}

add_completely.missing_back <- function(data, completely_missing){
    if(length(completely_missing)){
        tmp <- matrix(0L, nrow(data) + length(completely_missing), ncol(data))
        colnames(tmp) <- colnames(data)
        tmp[completely_missing, ] <- NA
        tmp[!(1:nrow(tmp) %in% completely_missing), ] <- data
        data <- tmp
    }
    data
}

replace_dash <- function(syntax){
    for(i in 1L:length(syntax)){
        tmp <- syntax[i]
        if(any(regexpr(",",tmp)))
            tmp <- strsplit(tmp,",")[[1L]]
        popout <- c()
        for(j in seq_len(length(tmp))){
            if(regexpr("-",tmp[j]) > 1L){
                popout <- c(popout,j)
                tmp2 <- as.numeric(strsplit(tmp[j],"-")[[1L]])
                tmp2 <- as.character(tmp2[1L]:tmp2[2L])
                tmp <- c(tmp[-j],tmp2)
            }
        }
        syntax[i] <- paste0(tmp, collapse=',')
    }
    syntax
}

hasConverged <- function(p0, p1, TOL){
    pick <- names(p0) %in% c('g', 'u')
    if(any(pick)){
        p0[pick] <- plogis(p0[pick])
        p1[pick] <- plogis(p1[pick])
    }
    pick <- names(p0) == c('PI')
    if(any(pick)){
        p0[pick] <- psumexp(p0[pick])
        p1[pick] <- psumexp(p1[pick])
    }
    all(abs(p0 - p1) < TOL)
}

QMC_quad <- function(npts, nfact, lim, leap=409, norm=FALSE){
    qnorm(QUnif(npts, min=0, max=1, p=nfact, leap=leap))
}

MC_quad <- function(npts, nfact, lim)
    qnorm(matrix(runif(n=npts * nfact, min = lim[1L], max = lim[2]), npts, nfact))

respSample <- function(P) .Call("respSample", P)

as.mirt_df <- function(df){
    class(df) <- c('mirt_df', class(df))
    df
}

as.mirt_matrix <- function(df){
    class(df) <- c('mirt_matrix', class(df))
    df
}

is.latent_regression <- function(mod){
    !is.null(mod@Data$covdata)
}

# subset of vcov for delta method
subset_vcov <- function(obj, vcov){
    nms <- colnames(vcov)
    splt <- strsplit(nms, "\\.")
    splt <- lapply(splt, function(x) as.integer(x[-1L]))
    pick <- sapply(splt, \(ind) ind %in% obj@parnum)
    vcov[pick, pick, drop=FALSE]
}

makeSymMat <- function(mat){
    if(ncol(mat) > 1L){
        mat[is.na(mat)] <- 0
        mat <- mat + t(mat) - diag(diag(mat))
    }
    mat
}

suppressMat <- function(res, suppress, upper = TRUE){
    if(!is.na(suppress)){
        if(upper){
            pick <- abs(res[upper.tri(res)]) < suppress
            res[upper.tri(res)][pick] <- NA
        } else {
            pick <- abs(res[lower.tri(res)]) < suppress
            res[lower.tri(res)][pick] <- NA
        }
        res[is.na(t(res) * res)] <- NA
    }
    res
}

missingMsg <- function(string)
    stop(paste0('\'', string, '\' argument is missing.'), call.=FALSE)

.mirtClusterEnv <- new.env(parent=emptyenv())
.mirtClusterEnv$ncores <- 1L
.mirtClusterEnv$omp_threads <- 1L

myApply <- function(X, MARGIN, FUN, progress = FALSE, ...){
    if(progress)
        return(t(pbapply::pbapply(X, MARGIN, FUN, ...,
                         cl=.mirtClusterEnv$MIRTCLUSTER)))
    if(.mirtClusterEnv$ncores > 1L){
        return(t(parallel::parApply(cl=.mirtClusterEnv$MIRTCLUSTER, X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    } else {
        return(t(apply(X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    }
}

myLapply <- function(X, FUN, progress = FALSE, ...){
    if(progress)
        return(t(pbapply::pblapply(X, FUN, ...,
                                  cl=.mirtClusterEnv$MIRTCLUSTER)))
    if(.mirtClusterEnv$ncores > 1L){
        return(parallel::parLapply(cl=.mirtClusterEnv$MIRTCLUSTER, X=X, fun=FUN, ...))
    } else {
        return(lapply(X=X, FUN=FUN, ...))
    }
}

mySapply <- function(X, FUN, progress = FALSE, ...){
    if(progress)
        return(t(pbapply::pbsapply(X, FUN, ...,
                                   cl=.mirtClusterEnv$MIRTCLUSTER)))
    if(.mirtClusterEnv$ncores > 1L){
        return(t(parallel::parSapply(cl=.mirtClusterEnv$MIRTCLUSTER, X=X, FUN=FUN, ...)))
    } else {
        return(t(sapply(X=X, FUN=FUN, ...)))
    }
}

crossprod_miss <- function(x, y){
    mat <- matrix(0, ncol(x), ncol(x))
    for(i in 1:ncol(x))
        for(j in 1:ncol(x))
            if(i <= j)
                mat[i,j] <- mat[j,i] <- sum(x[,i] * y[,j], na.rm=TRUE)
    mat
}

printf <- function(...) {
  catf(sprintf(...))
}

catf <- function(...) {
    cat(...)
    flush.console()
}
