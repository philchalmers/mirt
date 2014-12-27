# theta combinations
thetaComb <- function(theta, nfact)
{
	if (nfact == 1L){
        Theta <- matrix(theta)
	} else {
        thetalist <- vector('list', nfact)
        for(i in 1L:nfact)
            thetalist[[i]] <- theta
        Theta <- as.matrix(expand.grid(thetalist))
    }
	return(Theta)
}

# Product terms
prodterms <- function(theta0, prodlist)
{
    products <- matrix(1, ncol = length(prodlist), nrow = nrow(theta0))
    for(i in 1L:length(prodlist)){
        tmp <- prodlist[[i]]
        for(j in 1L:length(tmp))
            products[ ,i] <- products[ ,i] * theta0[ ,tmp[j]]
    }
    ret <- cbind(theta0,products)
    ret
}

# MH sampler for theta values
draw.thetas <- function(theta0, pars, fulldata, itemloc, cand.t.var, prior.t.var,
                        prior.mu, prodlist, OffTerm, CUSTOM.IND)
{
    N <- nrow(fulldata)
    J <- length(pars) - 1L
    unif <- runif(N)
    sigma <- if(ncol(theta0) == 1L) matrix(cand.t.var) else diag(rep(cand.t.var,ncol(theta0)))
    total_0 <- attr(theta0, 'log.lik_full')
    theta1 <- theta0 + mirt_rmvnorm(N, sigma = sigma)
    if(is.null(total_0)) theta1 <- theta0 #for intial draw
    log_den1 <- mirt_dmvnorm(theta1,prior.mu,prior.t.var,log=TRUE)
    if(length(prodlist) > 0L)
        theta1 <- prodterms(theta1,prodlist)
    itemtrace1 <- computeItemtrace(pars=pars, Theta=theta1, itemloc=itemloc,
                                   offterm=OffTerm, CUSTOM.IND=CUSTOM.IND)
    total_1 <- rowSums(fulldata * log(itemtrace1)) + log_den1
    if(!is.null(prodlist))
        theta1 <- theta1[ ,1L:(pars[[1L]]@nfact - pars[[1L]]@nfixedeffects -
                                   length(prodlist)), drop=FALSE]
    if(is.null(total_0)){ #for intial draw
        attr(theta1, 'log.lik_full') <- total_1
        return(theta1)
    }
    diff <- total_1 - total_0
    accept <- unif < exp(diff)
    theta1[!accept, ] <- theta0[!accept, ]
    total_1[!accept] <- total_0[!accept]
    log.lik <- sum(total_1)
    attr(theta1, "Proportion Accepted") <- sum(accept)/N
    attr(theta1, "log.lik") <- log.lik
    attr(theta1, 'log.lik_full') <- total_1
    return(theta1)
}

imputePars <- function(pars, covB, imputenums, constrain){
    shift <- mirt_rmvnorm(1L, sigma=covB)
    for(i in 1L:length(pars)){
        pn <- pars[[i]]@parnum
        pick2 <- imputenums %in% pn
        pick1 <- pn %in% imputenums
        pars[[i]]@par[pick1] <- pars[[i]]@par[pick1] + shift[pick2]
        if(is(pars[[i]], 'graded')){
            where <- (length(pars[[i]]@par) - pars[[i]]@ncat + 2L):length(pars[[i]]@par)
            if(!(all(sort(pars[[i]]@par[where], decreasing=TRUE) == pars[[i]]@par[where])))
                stop('Drawn values out of order')
        } else if(is(pars[[i]], 'grsm')){
            where <- (length(pars[[i]]@par) - pars[[i]]@ncat + 1L):(length(pars[[i]]@par)-1L)
            if(!(all(sort(pars[[i]]@par[where], decreasing=TRUE) == pars[[i]]@par[where])))
                stop('Drawn values out of order')
        }
    }
    if(length(constrain)){
        for(con in 1L:length(constrain)){
            tmp <- shift[imputenums %in% constrain[[con]][1L]]
            if(length(tmp)){
                for(i in 1L:length(pars)){
                    pick <- pars[[i]]@parnum %in% constrain[[con]][-1L]
                    pars[[i]]@par[pick] <- tmp + pars[[i]]@par[pick]
                }
            }
        }
    }
    return(pars)
}

# Rotation function
Rotate <- function(F, rotate, Target = NULL, ...)
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
	return(unclass(rotF))
}

# Gamma correlation, mainly for obtaining a sign
gamma.cor <- function(x)
{
	concordant <- function(x){
			mat.lr <- function(r, c){
				lr <- x[(r.x > r) & (c.x > c)]
				sum(lr)
			}
		r.x <- row(x)
		c.x <- col(x)
		sum(x * mapply(mat.lr, r = r.x, c = c.x))
	}
	discordant <- function(x){
		mat.ll <- function(r, c){
			ll <- x[(r.x > r) & (c.x < c)]
			sum(ll)
		}
		r.x <- row(x)
		c.x <- col(x)
		sum(x * mapply(mat.ll, r = r.x, c = c.x))
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
	nitems <- ncol(fulldata)
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

test_info <- function(pars, Theta, Alist, K){
    infolist <- list()
    for(cut in 1L:length(Alist)){
        A <- Alist[[cut]]
        info <- rep(0,nrow(Theta))
        for(j in 1L:length(K)){
            info <- info + ItemInfo(pars[[j]], A[j,], Theta)
        }
        infolist[[cut]] <- info
    }
    tmp <- 0
    for(i in 1L:length(infolist)){
        tmp <- tmp + infolist[[i]]
    }
    info <- tmp/length(infolist)
    info
}

Lambdas <- function(pars, Names, explor = FALSE, alpha = .05){
    J <- length(pars) - 1L
    lambdas <- lowerlambdas <- upperlambdas <-
        matrix(NA, J, length(ExtractLambdas(pars[[1L]])))
    gcov <- ExtractGroupPars(pars[[J+1L]])$gcov
    if(ncol(gcov) < ncol(lambdas)){
        tmpcov <- diag(ncol(lambdas))
        tmpcov[1L:ncol(gcov), 1L:ncol(gcov)] <- gcov
        gcov <- tmpcov
    }
    z <- qnorm(1 - alpha/2)
    rownames(lambdas) <- rownames(upperlambdas) <- rownames(lowerlambdas) <- Names
    for(i in 1L:J){
        tmp <- pars[[i]]
        lambdas[i,] <- ExtractLambdas(tmp) /1.702
        tmp@par <- pars[[i]]@par - z * pars[[i]]@SEpar
        lowerlambdas[i,] <- ExtractLambdas(tmp) /1.702
        tmp@par <- pars[[i]]@par + z * pars[[i]]@SEpar
        upperlambdas[i,] <- ExtractLambdas(tmp) /1.702
    }
    norm <- sqrt(1 + rowSums(lambdas^2))
    F <- as.matrix(lambdas/norm)
    if(!explor){
        norml <- sqrt(1 + rowSums(lowerlambdas^2, na.rm=TRUE))
        normh <- sqrt(1 + rowSums(upperlambdas^2, na.rm=TRUE))
        ret <- list(F=F, lower=as.matrix(lowerlambdas/norml),
                    upper=as.matrix(upperlambdas/normh))
    } else {
        ret <- list(F=F, lower=list(), upper=list())
    }
    ret
}

#change long pars for groups into mean in sigma
ExtractGroupPars <- function(x){
    nfact <- x@nfact
    gmeans <- x@par[1L:nfact]
    tmp <- x@par[-(1L:nfact)]
    gcov <- matrix(0, nfact, nfact)
    gcov[lower.tri(gcov, diag=TRUE)] <- tmp
    if(nfact != 1L)
        gcov <- gcov + t(gcov) - diag(diag(gcov))
    return(list(gmeans=gmeans, gcov=gcov))
}

reloadConstr <- function(par, constr, obj){
    par2 <- rep(NA, length(constr[[1L]]))
    notconstr <- rep(TRUE, length(par2))
    for(i in 1L:length(constr)){
        par2[constr[[i]]] <- par[i]
        notconstr[constr[[i]]] <- FALSE
    }
    par2[notconstr] <- par[(length(constr)+1L):length(par)]
    ind <- 1L
    for(i in 1L:length(obj)){
        obj[[i]]@par[obj[[i]]@est] <- par2[ind:(ind + sum(obj[[i]]@est) - 1L)]
        ind <- ind + sum(obj[[i]]@est)
    }
    return(obj)
}

bfactor2mod <- function(model, J){
    tmp <- tempfile('tempfile')
    unique <- sort(unique(model))
    index <- 1L:J
    tmp2 <- c()
    for(i in 1L:length(unique)){
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

updatePrior <- function(pars, Theta, Thetabetween, list, ngroups, nfact, J, N,
                        BFACTOR, sitems, cycles, rlist, prior, lrPars = list(), full=FALSE){
    Prior <- Priorbetween <- vector('list', ngroups)
    if(list$EH){
        Prior[[1L]] <- list$EHPrior[[1L]]
    } else {
        for(g in 1L:ngroups){
            gp <- ExtractGroupPars(pars[[g]][[J+1L]])
            if(BFACTOR){
                sel <- 1L:(nfact-ncol(sitems) + 1L)
                sel2 <- sel[-length(sel)]
                Priorbetween[[g]] <- mirt_dmvnorm(Thetabetween,
                                                      gp$gmeans[sel2], gp$gcov[sel2,sel2,drop=FALSE])
                Priorbetween[[g]] <- Priorbetween[[g]]/sum(Priorbetween[[g]])
                Prior[[g]] <- mirt_dmvnorm(Theta[ ,sel], gp$gmeans[sel], gp$gcov[sel,sel,drop=FALSE])
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
                next
            }
            if(full){
                Prior[[g]] <- mirt_dmvnorm(Theta[ ,1L:nfact,drop=FALSE], lrPars@mus, gp$gcov,
                                           quad=TRUE)
                Prior[[g]] <- Prior[[g]]/rowSums(Prior[[g]])
            } else {
                Prior[[g]] <- mirt_dmvnorm(Theta[ ,1L:nfact,drop=FALSE], gp$gmeans, gp$gcov)
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
            }
        }
    }
    if(list$EH){
        if(cycles > 1L){
            for(g in 1L:ngroups)
                Prior[[g]] <- rowSums(rlist[[g]][[1L]]) / sum(rlist[[g]][[1L]])
        } else {
            for(g in 1L:ngroups){
                Prior[[g]] <- mirt_dmvnorm(Theta, 0, matrix(1))
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
            }
        }
    } else if(!is.null(list$customPriorFun)){
        for(g in 1L:ngroups)
            Prior[[g]] <- list$customPriorFun(Theta, Etable=rlist[[g]][[1L]])
    }
    return(list(Prior=Prior, Priorbetween=Priorbetween))
}

UpdateConstrain <- function(pars, constrain, invariance, nfact, nLambdas, J, ngroups, PrepList,
                            method, itemnames, model, groupNames)
{
    if(!is.numeric(model[[1L]])){
        if(any(model[[1L]]$x[,1L] == 'CONSTRAIN')){
            groupNames <- as.character(groupNames)
            names(pars) <- groupNames
            input <- model[[1L]]$x[model[[1L]]$x[,1L] == 'CONSTRAIN', 2L]
            input <- gsub(' ', replacement='', x=input)
            elements <- strsplit(input, '\\),\\(')[[1L]]
            elements <- gsub('\\(', replacement='', x=elements)
            elements <- gsub('\\)', replacement='', x=elements)
            esplit <- strsplit(elements, ',')
            esplit <- lapply(esplit, function(x, groupNames)
                if(!(x[length(x)] %in% c(groupNames, 'all'))) c(x, 'all') else x,
                             groupNames=as.character(groupNames))
            esplit <- lapply(esplit, function(x){
                            newx <- c()
                            if(length(x) < 3L)
                                stop('PRIOR = ... has not been supplied enough arguments')
                            for(i in 1L:(length(x)-2L)){
                                if(grepl('-', x[i])){
                                    tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                                    newx <- c(newx, tmp[1L]:tmp[2L])
                                } else newx <- c(newx, x[i])
                            }
                            x <- c(newx, x[length(x)-1L], x[length(x)])
                            x
                        })
            for(i in 1L:length(esplit)){
                if(!(esplit[[i]][length(esplit[[i]])] %in% c(groupNames, 'all')))
                    stop('Invalid group name passed to CONSTRAIN = ... syntax.')
                if(esplit[[i]][length(esplit[[i]])] == 'all'){
                    for(g in 1L:ngroups){
                        constr <- c()
                        p <- pars[[g]]
                        sel <- suppressWarnings(
                            as.numeric(esplit[[i]][1L:(length(esplit[[i]])-1L)]))
                        picknames <- c(is.na(sel), FALSE)
                        sel <- na.omit(sel)
                        for(j in 1L:length(sel)){
                            pick <- p[[sel[j]]]@parnum[names(p[[sel[j]]]@est) %in%
                                                           esplit[[i]][picknames]]
                            if(!length(pick))
                                stop('CONSTRAIN = ... indexed a parameter that was not relavent for item ', sel[j])
                            constr <- c(constr, pick)
                        }
                        constrain[[length(constrain) + 1L]] <- constr
                    }
                } else {
                    constr <- c()
                    p <- pars[[esplit[[i]][length(esplit[[i]])]]]
                    sel <- as.numeric(esplit[[i]][1L:(length(esplit[[i]])-2L)])
                    for(j in 1L:length(sel)){
                        pick <- p[[sel[j]]]@parnum[names(p[[sel[j]]]@est) ==
                                                       esplit[[i]][length(esplit[[i]])-1L]]
                        if(!length(pick))
                            stop('CONSTRAIN = ... indexed a parameter that was not relavent for item ', sel[j])
                        constr <- c(constr, pick)
                    }
                    constrain[[length(constrain) + 1L]] <- constr
                }
            }
        }
        if(any(model[[1L]]$x[,1L] == 'CONSTRAINB')){
            if(length(unique(groupNames)) == 1L)
                stop('CONSTRAINB model argument not valid for single group models')
            groupNames <- as.character(groupNames)
            names(pars) <- groupNames
            input <- model[[1L]]$x[model[[1L]]$x[,1L] == 'CONSTRAINB', 2L]
            input <- gsub(' ', replacement='', x=input)
            elements <- strsplit(input, '\\),\\(')[[1L]]
            elements <- gsub('\\(', replacement='', x=elements)
            elements <- gsub('\\)', replacement='', x=elements)
            esplit <- strsplit(elements, ',')
            esplit <- lapply(esplit, function(x, groupNames)
                if(!(x[length(x)] %in% c(groupNames, 'all'))) c(x, 'all') else x,
                             groupNames=as.character(groupNames))
            esplit <- lapply(esplit, function(x){
                newx <- c()
                if(length(x) < 3L)
                    stop('PRIOR = ... has not been supplied enough arguments')
                for(i in 1L:(length(x)-2L)){
                    if(grepl('-', x[i])){
                        tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                        newx <- c(newx, tmp[1L]:tmp[2L])
                    } else newx <- c(newx, x[i])
                }
                x <- c(newx, x[length(x)-1L], x[length(x)])
                x
            })
            for(i in 1L:length(esplit)){
                sel <- as.numeric(esplit[[i]][1L:(length(esplit[[i]])-2L)])
                for(j in 1L:length(sel)){
                    constr <- c()
                    for(g in 1L:ngroups){
                        p <- pars[[g]]
                        pick <- p[[sel[j]]]@parnum[names(p[[sel[j]]]@est) ==
                                                       esplit[[i]][length(esplit[[i]])-1L]]
                        if(!length(pick))
                            stop('CONSTRAINB = ... indexed a parameter that was not relavent accross groups')
                        constr <- c(constr, pick)
                    }
                    constrain[[length(constrain) + 1L]] <- constr
                }
            }
        }
    }

    #within group item constraints only
    for(g in 1L:ngroups)
        if(length(PrepList[[g]]$constrain) > 0L)
            for(i in 1L:length(PrepList[[g]]$constrain))
                constrain[[length(constrain) + 1L]] <- PrepList[[g]]$constrain[[i]]
    if('covariances' %in% invariance){ #Fix covariance accross groups (only makes sense with vars = 1)
        tmpmat <- matrix(NA, nfact, nfact)
        low_tri <- lower.tri(tmpmat)
        tmp <- c()
        tmpmats <- tmpestmats <- matrix(NA, ngroups, nfact*(nfact+1L)/2)
        for(g in 1L:ngroups){
            tmpmats[g,] <- pars[[g]][[J + 1L]]@parnum[(nfact+1L):length(pars[[g]][[J + 1L]]@parnum)]
            tmpestmats[g,] <- pars[[g]][[J + 1L]]@est[(nfact+1L):length(pars[[g]][[J + 1L]]@est)]
        }
        select <- colSums(tmpestmats) == ngroups
        for(i in 1L:length(select))
            if(select[i])
                constrain[[length(constrain) + 1L]] <- tmpmats[1L:ngroups, i]

    }
    if(any(itemnames %in% invariance)){
        matched <- na.omit(match(invariance, itemnames))
        for(i in matched){
            jj <- sum(pars[[1L]][[i]]@est)
            stopifnot(jj > 0)
            for(j in 1L:jj){
                tmp <- c()
                for(g in 1L:ngroups)
                    tmp <- c(tmp, pars[[g]][[i]]@parnum[pars[[g]][[i]]@est][j])
                constrain[[length(constrain) + 1L]] <- tmp
            }
        }
    }
    if('slopes' %in% invariance){ #Equal factor loadings
        tmpmats <- tmpests <- list()
        for(g in 1L:ngroups)
            tmpmats[[g]] <- tmpests[[g]] <- matrix(NA, J, nLambdas)
        for(g in 1L:ngroups){
            for(i in 1L:J){
                tmpmats[[g]][i,] <- pars[[g]][[i]]@parnum[1L:nLambdas]
                tmpests[[g]][i,] <- pars[[g]][[i]]@est[1L:nLambdas]
            }
        }
        for(i in 1L:J){
            for(j in 1L:nLambdas){
                tmp <- c()
                for(g in 1L:ngroups){
                    if(tmpests[[1L]][[i, j]])
                        tmp <- c(tmp, tmpmats[[g]][i,j])
                }
                constrain[[length(constrain) + 1L]] <- tmp
            }
        }
    }
    if('intercepts' %in% invariance){ #Equal item intercepts (and all other item pars)
        tmpmats <- tmpests <- list()
        for(g in 1L:ngroups)
            tmpmats[[g]] <- tmpests[[g]] <- list()
        for(g in 1L:ngroups){
            for(i in 1L:J){
                ind <- (nLambdas+1L):length(pars[[g]][[i]]@parnum)
                if(is(pars[[g]][[i]], 'dich')) ind <- ind[1L:(length(ind)-2L)]
                if(is(pars[[g]][[i]], 'partcomp')) ind <- ind[1L:(length(ind)-1L)]
                tmpmats[[g]][[i]] <- pars[[g]][[i]]@parnum[ind]
                tmpests[[g]][[i]] <- pars[[g]][[i]]@est[ind]
            }
        }
        for(i in 1L:J){
            for(j in 1L:length(tmpmats[[1L]][[i]])){
                tmp <- c()
                for(g in 1L:ngroups){
                    if(tmpests[[1L]][[i]][j])
                        tmp <- c(tmp, tmpmats[[g]][[i]][j])
                }
                constrain[[length(constrain) + 1L]] <- tmp
            }
        }
    }
    #remove redundent constraints
    if(TRUE){
        redun <- rep(FALSE, length(constrain))
        if(length(constrain) > 0L){
            for(i in 1L:length(redun)){
                for(j in 1L:length(redun)){
                    if(j < i){
                        if(all(constrain[[i]] %in% constrain[[j]] ||
                                all(constrain[[j]] %in% constrain[[i]]))){
                            if(length(constrain[[i]]) < length(constrain[[j]])) redun[i] <- TRUE
                            else redun[j] <- TRUE
                        }
                    }
                }
            }
        }
        constrain[redun] <- NULL
    }
    return(constrain)
}

UpdatePrior <- function(PrepList, model, groupNames){
    if(!is.numeric(model[[1L]])){
        if(!length(model[[1L]]$x[model[[1L]]$x[,1L] == 'PRIOR', 2L])) return(PrepList)
        groupNames <- as.character(groupNames)
        ngroups <- length(groupNames)
        pars <- vector('list', length(PrepList))
        for(g in 1L:length(PrepList))
            pars[[g]] <- PrepList[[g]]$pars
        names(pars) <- groupNames
        input <- model[[1L]]$x[model[[1L]]$x[,1L] == 'PRIOR', 2L]
        input <- gsub(' ', replacement='', x=input)
        elements <- strsplit(input, '\\),\\(')[[1L]]
        elements <- gsub('\\(', replacement='', x=elements)
        elements <- gsub('\\)', replacement='', x=elements)
        esplit <- strsplit(elements, ',')
        esplit <- lapply(esplit, function(x, groupNames)
            if(!(x[length(x)] %in% c(groupNames, 'all'))) c(x, 'all') else x,
                         groupNames=as.character(groupNames))
        esplit <- lapply(esplit, function(x){
            newx <- c()
            if(length(x) < 5L)
                stop('PRIOR = ... has not been supplied enough arguments')
            for(i in 1L:(length(x)-5L)){
                if(grepl('-', x[i])){
                    tmp <- as.numeric(strsplit(x[i], '-')[[1L]])
                    newx <- c(newx, tmp[1L]:tmp[2L])
                } else newx <- c(newx, x[i])
            }
            x <- c(newx, x[(length(x)-4L):length(x)])
            x
        })
        for(i in 1L:length(esplit)){
            if(!(esplit[[i]][length(esplit[[i]])] %in% c(groupNames, 'all')))
                stop('Invalid group name passed to PRIOR = ... syntax.')
            if(esplit[[i]][length(esplit[[i]])] == 'all'){
                for(g in 1L:ngroups){
                    sel <- as.numeric(esplit[[i]][1L:(length(esplit[[i]])-5L)])
                    name <- esplit[[i]][length(esplit[[i]])-4L]
                    type <- esplit[[i]][length(esplit[[i]])-3L]
                    if(!(type %in% c('norm', 'beta', 'lnorm')))
                        stop('Prior type specified in PRIOR = ... not available')
                    type <- switch(type, norm=1L, lnorm=2L, beta=3L, 0L)
                    val1 <- as.numeric(esplit[[i]][length(esplit[[i]])-2L])
                    val2 <- as.numeric(esplit[[i]][length(esplit[[i]])-1L])
                    for(j in 1L:length(sel)){
                        which <- names(pars[[g]][[j]]@est) == name
                        if(!any(which)) stop('Parameter \'', name, '\' does not exist for item ', j)
                        pars[[g]][[sel[j]]]@any.prior <- TRUE
                        pars[[g]][[sel[j]]]@prior.type[which] <- type
                        pars[[g]][[sel[j]]]@prior_1[which] <- val1
                        pars[[g]][[sel[j]]]@prior_2[which] <- val2
                    }
                }
            } else {
                sel <- as.numeric(esplit[[i]][1L:(length(esplit[[i]])-5L)])
                gname <- esplit[[i]][length(esplit[[i]])]
                name <- esplit[[i]][length(esplit[[i]])-4L]
                type <- esplit[[i]][length(esplit[[i]])-3L]
                if(!(type %in% c('norm', 'beta', 'lnorm')))
                    stop('Prior type specified in PRIOR = ... not available')
                type <- switch(type, norm=1L, lnorm=2L, beta=3L, 0L)
                val1 <- as.numeric(esplit[[i]][length(esplit[[i]])-2L])
                val2 <- as.numeric(esplit[[i]][length(esplit[[i]])-1L])
                for(j in 1L:length(sel)){
                    which <- names(pars[[gname]][[j]]@est) == name
                    if(!any(which)) stop('Parameter \'', name, '\' does not exist for item ', j)
                    pars[[gname]][[sel[j]]]@any.prior <- TRUE
                    pars[[gname]][[sel[j]]]@prior.type[which] <- type
                    pars[[gname]][[sel[j]]]@prior_1[which] <- val1
                    pars[[gname]][[sel[j]]]@prior_2[which] <- val2
                }
            }
        }
        for(g in 1L:length(PrepList))
            PrepList[[g]]$pars <- pars[[g]]
    }
    return(PrepList)
}

ReturnPars <- function(PrepList, itemnames, random, lrPars, MG = FALSE){
    parnum <- par <- est <- item <- parname <- gnames <- class <-
        lbound <- ubound <- prior.type <- prior_1 <- prior_2 <- c()
    if(!MG) PrepList <- list(full=PrepList)
    for(g in 1L:length(PrepList)){
        tmpgroup <- PrepList[[g]]$pars
        for(i in 1L:length(tmpgroup)){
            if(i <= length(itemnames))
                item <- c(item, rep(itemnames[i], length(tmpgroup[[i]]@parnum)))
            class <- c(class, rep(class(tmpgroup[[i]]), length(tmpgroup[[i]]@parnum)))
            parname <- c(parname, names(tmpgroup[[i]]@est))
            parnum <- c(parnum, tmpgroup[[i]]@parnum)
            par <- c(par, tmpgroup[[i]]@par)
            est <- c(est, tmpgroup[[i]]@est)
            lbound <- c(lbound, tmpgroup[[i]]@lbound)
            ubound <- c(ubound, tmpgroup[[i]]@ubound)
            tmp <- sapply(as.character(tmpgroup[[i]]@prior.type),
                                 function(x) switch(x, '1'='norm', '2'='lnorm', '3'='beta', 'none'))
            prior.type <- c(prior.type, tmp)
            prior_1 <- c(prior_1, tmpgroup[[i]]@prior_1)
            prior_2 <- c(prior_2, tmpgroup[[i]]@prior_2)
        }
        item <- c(item, rep('GROUP', length(tmpgroup[[i]]@parnum)))
    }
    if(length(random) > 0L){
        for(i in 1L:length(random)){
            parname <- c(parname, names(random[[i]]@est))
            parnum <- c(parnum, random[[i]]@parnum)
            par <- c(par, random[[i]]@par)
            est <- c(est, random[[i]]@est)
            lbound <- c(lbound, random[[i]]@lbound)
            ubound <- c(ubound, random[[i]]@ubound)
            tmp <- sapply(as.character(random[[i]]@prior.type),
                          function(x) switch(x, '1'='norm', '2'='lnorm', '3'='beta', 'none'))
            prior.type <- c(prior.type, tmp)
            prior_1 <- c(prior_1, random[[i]]@prior_1)
            prior_2 <- c(prior_2, random[[i]]@prior_2)
            class <- c(class, rep('RandomPars', length(random[[i]]@parnum)))
            item <- c(item, rep('RANDOM', length(random[[i]]@parnum)))
        }
    }
    if(length(lrPars) > 0L){
        parname <- c(parname, names(lrPars@est))
        parnum <- c(parnum, lrPars@parnum)
        par <- c(par, lrPars@par)
        est <- c(est, lrPars@est)
        lbound <- c(lbound, lrPars@lbound)
        ubound <- c(ubound, lrPars@ubound)
        tmp <- sapply(as.character(lrPars@prior.type),
                      function(x) switch(x, '1'='norm', '2'='lnorm', '3'='beta', 'none'))
        prior.type <- c(prior.type, tmp)
        prior_1 <- c(prior_1, lrPars@prior_1)
        prior_2 <- c(prior_2, lrPars@prior_2)
        class <- c(class, rep('lrPars', length(lrPars@parnum)))
        item <- c(item, rep('BETA', length(lrPars@parnum)))
    }
    gnames <- rep(names(PrepList), each = length(est)/length(PrepList))
    par[parname %in% c('g', 'u')] <- antilogit(par[parname %in% c('g', 'u')])
    ret <- data.frame(group=gnames, item=item, class=class, name=parname, parnum=parnum, value=par,
                      lbound=lbound, ubound=ubound, est=est, prior.type=prior.type,
                      prior_1=prior_1, prior_2=prior_2)
    ret
}

UpdatePrepList <- function(PrepList, pars, random, lrPars = list(), MG = FALSE){
    currentDesign <- ReturnPars(PrepList, PrepList[[1L]]$itemnames, random=random,
                                lrPars=lrPars, MG = TRUE)
    if(nrow(currentDesign) != nrow(pars))
        stop('Rows in supplied and starting value data.frame objects do not match. Were the
             data or itemtype input arguments modified?')
    if(!all(as.matrix(currentDesign[,c('group', 'item', 'class', 'name', 'parnum')]) ==
                as.matrix(pars[,c('group', 'item', 'class', 'name', 'parnum')])))
        stop('Critical internal parameter labels do not match those returned from pars = \'values\'')
    if(!all(sapply(currentDesign, class) == sapply(pars, class)))
        stop('pars input does not contain the appropriate classes, which should match pars = \'values\'')
    if(!all(unique(pars$prior.type) %in% c('none', 'norm', 'beta', 'lnorm')))
        stop('prior.type input in pars contains invalid prior types')
    if(!MG) PrepList <- list(PrepList)
    len <- length(PrepList[[length(PrepList)]]$pars)
    maxparnum <- max(PrepList[[length(PrepList)]]$pars[[len]]@parnum)
    pars$value[pars$name %in% c('g', 'u')] <- logit(pars$value[pars$name %in% c('g', 'u')])
    if(PrepList[[1L]]$nfact > 1L)
        PrepList[[1L]]$exploratory <- all(pars$est[pars$name %in% paste0('a', 1L:PrepList[[1L]]$nfact)])
    ind <- 1L
    for(g in 1L:length(PrepList)){
        for(i in 1L:length(PrepList[[g]]$pars)){
            for(j in 1L:length(PrepList[[g]]$pars[[i]]@par)){
                PrepList[[g]]$pars[[i]]@par[j] <- pars[ind,'value']
                PrepList[[g]]$pars[[i]]@est[j] <- as.logical(pars[ind,'est'])
                PrepList[[g]]$pars[[i]]@lbound[j] <- pars[ind,'lbound']
                PrepList[[g]]$pars[[i]]@ubound[j] <- pars[ind,'ubound']
                tmp <- as.character(pars[ind,'prior.type'])
                PrepList[[g]]$pars[[i]]@prior.type[j] <-
                    switch(tmp, norm=1L, lnorm=2L, beta=3L, 0L)
                PrepList[[g]]$pars[[i]]@prior_1[j] <- pars[ind,'prior_1']
                PrepList[[g]]$pars[[i]]@prior_2[j] <- pars[ind,'prior_2']
                ind <- ind + 1L
            }
            if(is(PrepList[[g]]$pars[[i]], 'graded')){
                tmp <- ExtractZetas(PrepList[[g]]$pars[[i]])
                if(!all(tmp == sort(tmp, decreasing=TRUE)) || length(unique(tmp)) != length(tmp))
                    stop('Graded model intercepts for item ', i, ' in group ', g,
                         ' do not descend from highest to lowest. Please fix')
            }
            PrepList[[g]]$pars[[i]]@any.prior <- any(1L:3L %in%
                                                         PrepList[[g]]$pars[[i]]@prior.type)
        }
    }
    if(length(random) > 0L){
        for(i in 1L:length(random)){
            for(j in 1L:length(random[[i]]@par)){
                random[[i]]@par[j] <- pars[ind,'value']
                random[[i]]@est[j] <- as.logical(pars[ind,'est'])
                random[[i]]@lbound[j] <- pars[ind,'lbound']
                random[[i]]@ubound[j] <- pars[ind,'ubound']
                ind <- ind + 1L
            }
        }
        attr(PrepList, 'random') <- random
    }
    if(!MG) PrepList <- PrepList[[1L]]
    return(PrepList)
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
    if(any(x@prior.type %in% 3L)){ #beta
        ind <- x@prior.type %in% 3L
        val <- x@par[ind]
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
        for(i in 1L:length(val)){
            tmp <- dnorm(val[i], u[i], s[i], log=TRUE)
            LL <- LL + ifelse(tmp == -Inf, log(1e-100), tmp)
        }
    }
    if(any(x@prior.type %in% 2L)){
        ind <- x@prior.type %in% 2L
        val <- x@par[ind]
        val <- ifelse(val > 0, val, 1e-100)
        u <- x@prior_1[ind]
        s <- x@prior_2[ind]
        for(i in 1L:length(val)){
            tmp <- dlnorm(val[i], u[i], s[i], log=TRUE)
            LL <- LL + ifelse(tmp == -Inf, log(1e-100), tmp)
        }
    }
    if(any(x@prior.type %in% 3L)){
        ind <- x@prior.type %in% 3L
        val <- x@par[ind]
        a <- x@prior_1[ind]
        b <- x@prior_2[ind]
        for(i in 1L:length(val)){
            tmp <- dbeta(val[i], a[i], b[i], log=TRUE)
            LL <- LL + ifelse(tmp == -Inf, log(1e-100), tmp)
        }
    }
    return(LL)
}

ItemInfo <- function(x, Theta, cosangle, total.info = TRUE){
    P <- ProbTrace(x, Theta)
    dx <- DerivTheta(x, Theta)
    info <- matrix(0, nrow(Theta), ncol(P))
    cosanglefull <- matrix(cosangle, nrow(P), length(cosangle), byrow = TRUE)
    if(ncol(cosanglefull) < ncol(dx$grad[[1L]]))
        cosanglefull <- cbind(cosanglefull, matrix(1, nrow(cosanglefull),
                                                   ncol(dx$grad[[1L]]) - ncol(cosanglefull)))
    for(i in 1L:x@ncat)
        dx$grad[[i]] <- matrix(rowSums(dx$grad[[i]] * cosanglefull))
    for(i in 1L:x@ncat)
        info[,i] <- (dx$grad[[i]])^2 / P[ ,i]
    if(total.info) info <- matrix(rowSums(info))
    return(info)
}

ItemInfo2 <- function(x, Theta, total.info = TRUE){
    P <- ProbTrace(x, Theta)
    dx <- DerivTheta(x, Theta)
    info <- matrix(0, nrow(Theta), ncol(P))
    for(i in 1L:x@ncat)
        info[,i] <- (dx$grad[[i]])^2 / P[ ,i]
    if(total.info) info <- rowSums(info)
    return(info)
}

nameInfoMatrix <- function(info, correction, L, npars){
    #give info meaningful names for wald test
    parnames <- names(correction)
    tmp <- outer(1L:npars, rep(1L, npars))
    matind <- matrix(0, ncol(tmp), nrow(tmp))
    matind[lower.tri(matind, diag = TRUE)] <- tmp[lower.tri(tmp, diag = TRUE)]
    matind <- matind * L
    matind[matind == 0 ] <- NA
    matind[!is.na(matind)] <- tmp[!is.na(matind)]
    shortnames <- c()
    for(i in 1L:length(correction)){
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

maketabData <- function(stringfulldata, stringtabdata, group, groupNames, nitem, K, itemloc,
                        Names, itemnames, survey.weights){
    tabdata2 <- lapply(strsplit(stringtabdata, split='/'), as.integer)
    tabdata2 <- do.call(rbind, tabdata2)
    tabdata2[tabdata2 == 99999L] <- NA
    tabdata <- matrix(0L, nrow(tabdata2), sum(K))
    for(i in 1L:nitem){
        uniq <- sort(na.omit(unique(tabdata2[,i])))
        if(length(uniq) < K[i]) uniq <- 0L:(K[i]-1L)
        for(j in 1L:length(uniq))
            tabdata[,itemloc[i] + j - 1L] <- as.integer(tabdata2[,i] == uniq[j])
    }
    tabdata[is.na(tabdata)] <- 0L
    colnames(tabdata) <- Names
    colnames(tabdata2) <- itemnames
    groupFreq <- vector('list', length(groupNames))
    names(groupFreq) <- groupNames
    for(g in 1L:length(groupNames)){
        Freq <- integer(length(stringtabdata))
        tmpstringdata <- stringfulldata[group == groupNames[g]]
        if(!is.null(survey.weights)){
            mtc <- match(tmpstringdata, stringtabdata)
            Freq <- mySapply(1L:nrow(tabdata), function(x, std, tstd, w)
                sum(w[stringtabdata[x] == tstd]), std=stringtabdata, tstd=tmpstringdata,
                w=survey.weights)
        } else {
            Freq[stringtabdata %in% tmpstringdata] <- as.integer(table(
                match(tmpstringdata, stringtabdata)))
        }
        groupFreq[[g]] <- Freq
    }
    ret <- list(tabdata=tabdata, tabdata2=tabdata2, Freq=groupFreq)
    ret
}

makeLmats <- function(pars, constrain, random = list(), lrPars = list()){
    f <- function(k) (k+1) / (k*2)
    ngroups <- length(pars)
    J <- length(pars[[1L]]) - 1L
    L <- c()
    for(g in 1L:ngroups)
        for(i in 1L:(J+1L))
            L <- c(L, pars[[g]][[i]]@est)
    if(length(random))
        for(i in 1L:length(random))
            L <- c(L, random[[i]]@est)
    if(length(lrPars))
        L <- c(L, lrPars@est)
    L <- diag(as.numeric(L))
    redun_constr <- rep(FALSE, ncol(L))
    if(length(constrain) > 0L){
        for(i in 1L:length(constrain)){
            L[constrain[[i]], constrain[[i]]] <- 1L
            for(j in 2L:length(constrain[[i]]))
                redun_constr[constrain[[i]][j]] <- TRUE
        }
    }
    return(list(L=L, redun_constr=redun_constr))
}

updateGrad <- function(g, L) L %*% g
updateHess <- function(h, L) L %*% h %*% L

makeopts <- function(method = 'MHRM', draws = 2000L, calcLL = TRUE, quadpts = NULL,
                     rotate = 'varimax', Target = NaN, SE = FALSE, verbose = TRUE,
                     SEtol = .001, grsm.block = NULL, D = 1, TOL = NULL,
                     rsm.block = NULL, calcNull = TRUE, BFACTOR = FALSE,
                     technical = list(), use = 'pairwise.complete.obs',
                     SE.type = 'crossprod', large = NULL, accelerate = 'Ramsay', empiricalhist = FALSE,
                     optimizer = NULL, solnp_args = list(), alabama_args = list(), ...)
{
    opts <- list()
    if(method == 'MHRM' || method == 'MIXED') SE.type <- 'MHRM'
    if(!(method %in% c('MHRM', 'MIXED', 'BL', 'EM', 'QMCEM')))
        stop('method argument not supported')
    D <- 1
    opts$method = method
    if(draws < 1) stop('draws must be greater than 0')
    opts$draws = draws
    opts$calcLL = calcLL
    opts$quadpts = quadpts
    opts$rotate = rotate
    opts$Target = Target
    opts$SE = SE
    opts$SE.type = SE.type
    opts$verbose = verbose
    opts$SEtol = ifelse(is.null(technical$SEtol), .001, technical$SEtol)
    opts$grsm.block = grsm.block
    opts$D = D
    opts$rsm.block = rsm.block
    opts$calcNull = calcNull
    opts$customPriorFun = technical$customPriorFun
    opts$BFACTOR = BFACTOR
    opts$accelerate = accelerate
    opts$TOL <- ifelse(is.null(TOL), if(method == 'EM' || method == 'QMCEM') 1e-4 else
        if(method == 'BL') 1e-8 else 1e-3, TOL)
    if(SE.type == 'SEM' && SE){
        opts$accelerate <- 'none'
        if(is.null(TOL)) opts$TOL <- 1e-5
        if(is.null(technical$NCYCLES)) technical$NCYCLES <- 1000L
    }
    if(is.null(technical$symmetric_SEM)) technical$symmetric_SEM <- TRUE
    opts$warn <- if(is.null(technical$warn)) TRUE else technical$warn
    opts$message <- if(is.null(technical$message)) TRUE else technical$message
    opts$technical <- technical
    opts$use <- use
    opts$technical$parallel <- ifelse(is.null(technical$parallel), TRUE, technical$parallel)
    opts$MAXQUAD <- ifelse(is.null(technical$MAXQUAD), 10000L, technical$MAXQUAD)
    opts$NCYCLES <- ifelse(is.null(technical$NCYCLES), 2000L, technical$NCYCLES)
    if(opts$method %in% c('EM', 'QMCEM'))
        opts$NCYCLES <- ifelse(is.null(technical$NCYCLES), 500L, technical$NCYCLES)
    opts$BURNIN <- ifelse(is.null(technical$BURNIN), 150L, technical$BURNIN)
    opts$SEMCYCLES <- ifelse(is.null(technical$SEMCYCLES), 50, technical$SEMCYCLES)
    opts$KDRAWS  <- ifelse(is.null(technical$KDRAWS), 1L, technical$KDRAWS)
    opts$empiricalhist <- empiricalhist
    if(empiricalhist){
        if(opts$method != 'EM')
            stop('empirical histogram method only applicable when method = \'EM\' ')
        if(opts$TOL == 1e-4) opts$TOL <- 3e-5
        if(is.null(opts$quadpts)) opts$quadpts <- 199L
        if(opts$NCYCLES == 500L) opts$NCYCLES <- 2000L
    }
    if(method == 'QMCEM' && is.null(opts$quadpts)) opts$quadpts <- 2000L
    opts$MSTEPTOL <- ifelse(is.null(technical$MSTEPTOL), opts$TOL/1000, technical$MSTEPTOL)
    if(opts$method == 'MHRM' || opts$method =='MIXED' || SE.type == 'MHRM')
        set.seed(12345L)
    if(!is.null(technical$set.seed)) set.seed(technical$set.seed)
    opts$gain <- c(0.15,0.65)
    if(!is.null(technical$gain)){
        if(length(technical$gain) == 2L && is.numeric(technical$gain))
            opts$gain <- technical$gain
    }
    opts$NULL.MODEL <- ifelse(is.null(technical$NULL.MODEL), FALSE, TRUE)
    opts$USEEM <- ifelse(method == 'EM', TRUE, FALSE)
    opts$returnPrepList <- FALSE
    opts$PrepList <- NULL
    if(is.null(optimizer)){
        opts$Moptim <- if(method %in% c('EM','BL','QMCEM')) 'BFGS' else 'NR'
    } else {
        opts$Moptim <- optimizer
    }
    if(opts$Moptim == 'solnp'){
        if(is.null(solnp_args$control)) solnp_args$control <- list()
        if(is.null(solnp_args$control$trace)) solnp_args$control$trace <- 0
        if(method != 'EM') stop('solnp only supported for optimization with EM algorithm')
        opts$solnp_args <- solnp_args
    } else if(opts$Moptim == 'alabama'){
        if(method != 'EM') stop('alabama only supported for optimization with EM algorithm')
        if(is.null(alabama_args$control.outer)) alabama_args$control.outer <- list()
        if(is.null(alabama_args$control.optim)) alabama_args$control.optim <- list()
        if(is.null(alabama_args$control.outer$trace)) alabama_args$control.outer$trace <- FALSE
        opts$solnp_args <- alabama_args
    }
    if(!is.null(large)){
        if(is.logical(large))
            if(large) opts$returnPrepList <- TRUE
        if(is.list(large)) opts$PrepList <- large
    }
    if(!is.null(technical$customK)) opts$calcNull <- FALSE
    return(opts)
}

reloadPars <- function(longpars, pars, ngroups, J){
    if(FALSE){
        pars <- .Call('reloadPars', longpars, pars, ngroups, J)
    } else {
        #slower version for now till R 3.1.0 evaluation bug gets fixed. FIXME
        ind <- 1L
        for(g in 1L:ngroups){
            for(i in 1L:(J+1L)){
                tmp <- pars[[g]][[i]]@par
                pars[[g]][[i]]@par <- longpars[ind:(ind+length(tmp)-1L)]
                ind <- ind + length(tmp)
            }
        }
    }
    return(pars)
}

computeItemtrace <- function(pars, Theta, itemloc, offterm = matrix(0L, 1L, length(itemloc)-1L),
                             CUSTOM.IND){
    itemtrace <- .Call('computeItemTrace', pars, Theta, itemloc, offterm)
    if(length(CUSTOM.IND)){
        for(i in CUSTOM.IND)
            itemtrace[,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(pars[[i]], Theta=Theta)
    }
    return(itemtrace)
}

assignItemtrace <- function(pars, itemtrace, itemloc){
    for(i in 1L:(length(pars)-1L))
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
            warning('Could not invert information matrix; model likely is not identified.')
        ESTIMATE$fail_invert_info <- TRUE
        return(ESTIMATE)
    } else ESTIMATE$fail_invert_info <- FALSE
    SEtmp <- diag(solve(info))
    if(any(SEtmp < 0)){
        if(warn)
            warning("Negative SEs set to NaN.\n")
        SEtmp[SEtmp < 0 ] <- NaN
    }
    SEtmp <- sqrt(SEtmp)
    SE <- rep(NA, length(longpars))
    SE[ESTIMATE$estindex_unique[!isna]] <- SEtmp
    index <- 1L:length(longpars)
    if(length(constrain) > 0L)
        for(i in 1L:length(constrain))
            SE[index %in% constrain[[i]][-1L]] <- SE[constrain[[i]][1L]]
    ind1 <- 1L
    for(g in 1L:ngroups){
        for(i in 1L:(J+1L)){
            ind2 <- ind1 + length(pars[[g]][[i]]@par) - 1L
            pars[[g]][[i]]@SEpar <- SE[ind1:ind2]
            ind1 <- ind2 + 1L
        }
    }
    ESTIMATE$pars <- pars
    return(ESTIMATE)
}

make.randomdesign <- function(random, longdata, covnames, itemdesign, N){
    itemcovnames <- colnames(itemdesign)
    J <- nrow(itemdesign)
    ret <- vector('list', length(random))
    for(i in 1L:length(random)){
        f <- gsub(" ", "", as.character(random[[i]])[2L])
        splt <- strsplit(f, '\\|')[[1L]]
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

make.lrdesign <- function(df, formula, factorNames, EM=FALSE){
    nfact <- length(factorNames)
    if(is.list(formula)){
        if(!all(names(formula) %in% factorNames))
            stop('List of fixed effect names do not match factor names')
        estnames <- X <- vector('list', length(formula))
        for(i in 1L:length(formula)){
            X[[i]] <- model.matrix(formula[[i]], df)
            estnames[[i]] <- colnames(X[[i]])
        }
        X <- do.call(cbind, X)
        X <- X[,unique(colnames(X))]
    } else {
        X <- model.matrix(formula, df)
    }
    tXX <- t(X) %*% X
    if(ncol(X) > 1) inv_tXX <- solve(tXX)
    else inv_tXX <- matrix(0)
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
    est[1,] <- FALSE
    est <- as.logical(est)
    names(est) <- as.character(t(outer(factorNames, colnames(X),
                                     FUN = function(X, Y) paste(X,Y,sep="_"))))
    colnames(beta) <- factorNames
    rownames(beta) <- colnames(X)
    par <- as.numeric(beta)
    ret <- new('lrPars',
               par=par,
               SEpar=rep(NaN,length(par)),
               est=est,
               beta=beta,
               sigma=sigma,
               nfact=nfact,
               nfixed=ncol(X),
               df=df,
               X=X,
               tXX=tXX,
               inv_tXX=inv_tXX,
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


OffTerm <- function(random, J, N){
    ret <- numeric(N*J)
    for(i in 1L:length(random)){
        tmp <- rowSums(random[[i]]@gdesign*random[[i]]@drawvals[random[[i]]@mtch, ,drop=FALSE])
        ret <- ret + tmp
    }
    return(matrix(ret, N, J))
}

reloadRandom <- function(random, longpars, parstart){
    ind1 <- parstart
    for(i in 1L:length(random)){
        ind2 <- ind1 + length(random[[i]]@par) - 1L
        random[[i]]@par <- longpars[ind1:ind2]
        ind1 <- ind2 + 1L
    }
    random
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

assignInformationMG <- function(object){
    J <- ncol(object@Data$data)
    names <- colnames(object@information)
    spl_names <- strsplit(names, split="\\.")
    spl_names_par <- sapply(spl_names, function(x) x[1L])
    spl_names <- lapply(spl_names,
                        function(x) as.numeric(x[-1L]))
    spl_names <- do.call(rbind, spl_names)
    for(g in 1L:length(object@pars)){
        from <- object@pars[[g]]@pars[[1L]]@parnum[1L]
        to <- object@pars[[g]]@pars[[J+1L]]@parnum[length(
            object@pars[[g]]@pars[[J+1L]]@parnum)]
        pick <- spl_names[,g] >= from & spl_names[,g] <= to
        tmp <- object@information[pick,pick]
        colnames(tmp) <- rownames(tmp) <-
            paste(spl_names_par[pick], spl_names[pick,g], sep='.')
        object@pars[[g]]@information <- tmp
    }
    object
}

BL.LL <- function(p, est, longpars, pars, ngroups, J, Theta, PrepList, specific, sitems,
               CUSTOM.IND, EH, EHPrior, Data, BFACTOR, itemloc, theta){
    longpars[est] <- p
    pars2 <- reloadPars(longpars=longpars, pars=pars, ngroups=ngroups, J=J)
    gstructgrouppars <- prior <- Prior <- vector('list', ngroups)
    if(EH){
        Prior[[1L]] <- EHPrior[[1L]]
    } else {
        for(g in 1L:ngroups){
            gstructgrouppars[[g]] <- ExtractGroupPars(pars2[[g]][[J+1L]])
            if(BFACTOR){
                prior[[g]] <- dnorm(theta, 0, 1)
                prior[[g]] <- prior[[g]]/sum(prior[[g]])
                Prior[[g]] <- apply(expand.grid(prior[[g]], prior[[g]]), 1L, prod)
                next
            }
            Prior[[g]] <- mirt_dmvnorm(Theta,gstructgrouppars[[g]]$gmeans,
                                       gstructgrouppars[[g]]$gcov)
            Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
        }
    }
    LL <- 0
    for(g in 1L:ngroups){
        if(BFACTOR){
            expected <- Estep.bfactor(pars=pars2[[g]],
                                      tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                      Theta=Theta, prior=prior[[g]],
                                      specific=specific, sitems=sitems,
                                      itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)$expected
        } else {
            expected <- Estep.mirt(pars=pars2[[g]],
                                   tabdata=Data$tabdatalong, freq=Data$Freq[[g]],
                                   Theta=Theta, prior=Prior[[g]], itemloc=itemloc,
                                   CUSTOM.IND=CUSTOM.IND, full=FALSE)$expected
        }
        LL <- LL + sum(Data$Freq[[g]] * log(expected), na.rm = TRUE)
    }
    LL
}

select_quadpts <- function(nfact) switch(as.character(nfact),
                                         '1'=61, '2'=31, '3'=15, '4'=9, '5'=7, 3)

select_quadpts2 <- function(nfact) switch(as.character(nfact),
                                          '1'=41, '2'=21, '3'=11, '4'=7, '5'=5, 3)

mirt_rmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                         check = FALSE)
{
    # Version modified from mvtnorm::rmvnorm, version 0.9-9996, 19-April, 2014.
    if(check){
        if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE))
            stop("sigma must be a symmetric matrix")
        if (length(mean) != nrow(sigma))
            stop("mean and sigma have non-conforming size")
    }
    ev <- eigen(sigma, symmetric = TRUE)
    if(check)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1])))
            warning("sigma is numerically not positive definite")
    retval <- ev$vectors %*%  diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*%  retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}

mirt_dmvnorm <- function(x, mean, sigma, log = FALSE, quad = FALSE)
{
    if(quad && is.matrix(mean)){
        isigma <- solve(sigma)
        distval <- matrix(0, nrow(mean), nrow(x))
        for(i in 1L:nrow(mean)){
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
    logdet <- sum(log(eigen(sigma, symmetric=TRUE,
                            only.values=TRUE)$values))
    logretval <- -(ncol(x)*log(2*pi) + logdet + distval)/2
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

makeObstables <- function(dat, K){
    ret <- vector('list', ncol(dat))
    sumscore <- rowSums(dat)
    for(i in 1L:length(ret)){
        ret[[i]] <- matrix(0, sum(K-1L)+1L, K[i])
        colnames(ret[[i]]) <- paste0(1L:K[i]-1L)
        rownames(ret[[i]]) <- paste0(1L:nrow(ret[[i]])-1L)
        split <- by(sumscore, dat[,i], table)
        for(j in 1L:K[i]){
            m <- match(names(split[[j]]), rownames(ret[[i]]))
            ret[[i]][m,j] <- split[[j]]
        }
        ret[[i]] <- ret[[i]][-c(1L, nrow(ret[[i]])), ]
    }
    ret
}

collapseCells <- function(O, E, mincell = 1){
    for(i in 1L:length(O)){
        On <- O[[i]]
        En <- E[[i]]
        drop <- which(rowSums(is.na(En)) > 0)
        En[is.na(En)] <- 0

        #collapse known upper and lower sparce cells
        if(length(drop) > 0L){
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

        #drop 0's and 1's
        drop <- rowSums(On) == 0L
        On <- On[!drop,]
        En <- En[!drop,]
        L <- En < mincell
        drop <- c()
        for(j in 1L:(nrow(On)-1L)){
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
            for(j in 1L:nrow(En)){
                if(!any(L[j,])) next
                tmp <- En[j, ]
                tmp2 <- On[j, ]
                while(length(tmp) > 2L){
                    m <- min(tmp)
                    whc <- which(m == tmp)
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

        #drop columns if they are very rare

        En[is.na(En)] <- 0
        dropcol <- logical(ncol(En))
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

rmsea <- function(X2, df, N){
    ret <- ifelse((X2 - df) > 0,
                  sqrt(X2 - df) / sqrt(df * (N-1)), 0)
    ret <- ifelse(is.na(ret), NaN, ret)
    ret
}

controlCandVar <- function(PA, cand, min = .1, max = .6){
    if(PA > max) cand <- cand * 1.05
    else if(PA < min) cand <- cand * 0.9
    if(cand < .001) cand <- .001
    cand
}

mirtClusterEnv <- new.env()
mirtClusterEnv$ncores <- 1L

myApply <- function(X, MARGIN, FUN, ...){
    if(mirtClusterEnv$ncores > 1L){
        return(t(parallel::parApply(cl=mirtClusterEnv$MIRTCLUSTER, X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    } else {
        return(t(apply(X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    }
}

myLapply <- function(X, FUN, ...){
    if(mirtClusterEnv$ncores > 1L){
        return(parallel::parLapply(cl=mirtClusterEnv$MIRTCLUSTER, X=X, fun=FUN, ...))
    } else {
        return(lapply(X=X, FUN=FUN, ...))
    }
}

mySapply <- function(X, FUN, ...){
    if(mirtClusterEnv$ncores > 1L){
        return(t(parallel::parSapply(cl=mirtClusterEnv$MIRTCLUSTER, X=X, FUN=FUN, ...)))
    } else {
        return(t(sapply(X=X, FUN=FUN, ...)))
    }
}
