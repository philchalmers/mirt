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
    theta1 <- theta0 + mvtnorm::rmvnorm(N,prior.mu, sigma)
    log_den0 <- mvtnorm::dmvnorm(theta0,prior.mu,prior.t.var,log=TRUE)
    log_den1 <- mvtnorm::dmvnorm(theta1,prior.mu,prior.t.var,log=TRUE)
    if(length(prodlist) > 0L){
        theta0 <- prodterms(theta0,prodlist)
        theta1 <- prodterms(theta1,prodlist)
    }
    itemtrace0 <- computeItemtrace(pars=pars, Theta=theta0, itemloc=itemloc,
                                   offterm=OffTerm, CUSTOM.IND=CUSTOM.IND)
    itemtrace1 <- computeItemtrace(pars=pars, Theta=theta1, itemloc=itemloc,
                                   offterm=OffTerm, CUSTOM.IND=CUSTOM.IND)
    totals <- .Call('denRowSums', fulldata, itemtrace0, itemtrace1, log_den0, 
                    log_den1, mirtClusterEnv$ncores)
    total_0 <- totals[[1L]]
    total_1 <- totals[[2L]]
    diff <- total_1 - total_0
    accept <- diff > 0
    accept[unif < exp(diff)] <- TRUE
    theta1[!accept, ] <- theta0[!accept, ]
    total_1[!accept] <- total_0[!accept]
    log.lik <- sum(total_1)
    if(!is.null(prodlist))
        theta1 <- theta1[ ,1L:(pars[[1L]]@nfact - pars[[1L]]@nfixedeffects -
                                  length(prodlist)), drop=FALSE]
    attr(theta1, "Proportion Accepted") <- sum(accept)/N
    attr(theta1, "log.lik") <- log.lik
    return(theta1)
}

imputePars <- function(pars, covB, imputenums, constrain){
    shift <- mvtnorm::rmvnorm(1L, sigma=covB)
    for(i in 1L:length(pars)){
        pn <- pars[[i]]@parnum 
        pick2 <- imputenums %in% pn
        pick1 <- pn %in% imputenums
        pars[[i]]@par[pick1] <- pars[[i]]@par[pick1] + shift[pick2]
        if(is(pars[[i]], 'graded')){
            where <- (length(pars[[i]]@par) - pars[[i]]@ncat + 2L):length(pars[[i]]@par)
            pars[[i]]@par[where] <- sort(pars[[i]]@par[where], decreasing=TRUE)
        } else if(is(pars[[i]], 'grsm')){
            where <- (length(pars[[i]]@par) - pars[[i]]@ncat + 1L):(length(pars[[i]]@par)-1L)
            pars[[i]]@par[where] <- sort(pars[[i]]@par[where], decreasing=TRUE)
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
    ret <- log(x / (1 - x))
    ret <- ifelse(x == 0, -999, ret)
    ret <- ifelse(x == 1, 999, ret)
    ret
}

antilogit <- function(x){
    ret <- 1 / (1 + exp(-x))
    ret <- ifelse(x == -999, 0, ret)
    ret <- ifelse(x == 999, 1, ret)
    ret
}

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
    z <- qnorm(1 - alpha/2)
    rownames(lambdas) <- rownames(upperlambdas) <- rownames(lowerlambdas) <- Names
    for(i in 1L:J){
        tmp <- pars[[i]]
        lambdas[i,] <- ExtractLambdas(tmp)/1.702
        tmp@par <- pars[[i]]@par - z * pars[[i]]@SEpar
        lowerlambdas[i,] <- ExtractLambdas(tmp)/1.702
        tmp@par <- pars[[i]]@par + z * pars[[i]]@SEpar
        upperlambdas[i,] <- ExtractLambdas(tmp)/1.702
    }
    if(ncol(lambdas) > 1L){
        norm <- sqrt(1 + rowSums(lambdas^2))
    } else {
        norm <- as.matrix(sqrt(1 + lambdas[ ,1L]^2))
    }
    F <- as.matrix(lambdas/norm)
    if(!explor){
        if(ncol(lambdas) > 1L){
            norml <- sqrt(1 + rowSums(lowerlambdas^2))
            normh <- sqrt(1 + rowSums(upperlambdas^2))
        } else {
            norml <- as.matrix(sqrt(1 + lowerlambdas[ ,1L]^2))
            normh <- as.matrix(sqrt(1 + upperlambdas[ ,1L]^2))
        }
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

updatePrior <- function(pars, Theta, Thetabetween, list, ngroups, nfact, J, 
                        BFACTOR, sitems, cycles, rlist, prior){
    Prior <- Priorbetween <- vector('list', ngroups)
    if(list$EH){
        Prior[[1L]] <- list$EHPrior[[1L]]
    } else {
        for(g in 1L:ngroups){
            gp <- ExtractGroupPars(pars[[g]][[J+1L]])
            if(BFACTOR){
                sel <- 1L:(nfact-ncol(sitems) + 1L)
                sel2 <- sel[-length(sel)]
                Priorbetween[[g]] <- mvtnorm::dmvnorm(Thetabetween,
                                                      gp$gmeans[sel2], gp$gcov[sel2,sel2,drop=FALSE])
                Priorbetween[[g]] <- Priorbetween[[g]]/sum(Priorbetween[[g]])
                Prior[[g]] <- mvtnorm::dmvnorm(Theta[ ,sel], gp$gmeans[sel], gp$gcov[sel,sel,drop=FALSE])
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
                next
            }
            Prior[[g]] <- mvtnorm::dmvnorm(Theta[ ,1L:nfact,drop=FALSE], gp$gmeans, gp$gcov)
            Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
        }
    }
    if(list$EH){
        if(cycles > 1L){
            for(g in 1L:ngroups)
                Prior[[g]] <- rowSums(rlist[[g]][[1L]]) / sum(rlist[[g]][[1L]])
        } else {
            for(g in 1L:ngroups){
                Prior[[g]] <- mvtnorm::dmvnorm(Theta, 0, matrix(1))
                Prior[[g]] <- Prior[[g]]/sum(Prior[[g]])
            }
        }
    } else if(!is.null(list$customPriorFun)){
        for(g in 1L:ngroups)
            Prior[[g]] <- list$customPriorFun(Theta)
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

ReturnPars <- function(PrepList, itemnames, random, MG = FALSE){
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
    gnames <- rep(names(PrepList), each = length(est)/length(PrepList))
    par[parname %in% c('g', 'u')] <- antilogit(par[parname %in% c('g', 'u')])
    ret <- data.frame(group=gnames, item=item, class=class, name=parname, parnum=parnum, value=par,
                      lbound=lbound, ubound=ubound, est=est, prior.type=prior.type,
                      prior_1=prior_1, prior_2=prior_2)
    ret
}

UpdatePrepList <- function(PrepList, pars, random, MG = FALSE){
    currentDesign <- ReturnPars(PrepList, PrepList[[1L]]$itemnames, random=random, MG = TRUE)
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
    for(i in 1L:x@ncat){
        dx$grad[[i]] <- matrix(rowSums(dx$grad[[i]] * cosanglefull))
        dx$hess[[i]] <- matrix(rowSums(dx$hess[[i]] * cosanglefull))
    }
    for(i in 1L:x@ncat)
        info[,i] <- (dx$grad[[i]])^2 / P[ ,i] - dx$hess[[i]]
    if(total.info) info <- matrix(rowSums(info))
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
                        Names, itemnames){
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
    ret1 <- ret2 <- vector('list', length(groupNames))
    for(g in 1L:length(groupNames)){
        Freq <- integer(length(stringtabdata))
        tmpstringdata <- stringfulldata[group == groupNames[g]]
        Freq[stringtabdata %in% tmpstringdata] <- as.integer(table(match(tmpstringdata, stringtabdata)))
        ret1[[g]] <- cbind(tabdata, Freq)
        ret2[[g]] <- cbind(tabdata2, Freq)
    }
    ret <- list(tabdata=ret1, tabdata2=ret2)
    ret
}

makeLmats <- function(pars, constrain, random = NULL){
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
    L <- L2 <- L3 <- diag(as.numeric(L))
    redun_constr <- rep(FALSE, ncol(L))
    if(length(constrain) > 0L){
        for(i in 1L:length(constrain)){
            LC <- length(constrain[[i]])
            L[constrain[[i]], constrain[[i]]] <- 1/LC
            L2[constrain[[i]], constrain[[i]]] <- sqrt(1/LC)
            L3[constrain[[i]], constrain[[i]]] <- f(LC)
            for(j in 2L:length(constrain[[i]]))
                redun_constr[constrain[[i]][j]] <- TRUE
        }
    }
    return(list(L=L, L2=L2, L3=L3, redun_constr=redun_constr))
}

updateHess <- function(h, L2, L3){
    hess <- L3 %*% h %*% L3
    tmp <- L2 %*% h %*% L2
    diag(hess) <- diag(tmp)
    hess
}

makeopts <- function(method = 'MHRM', draws = 2000L, calcLL = TRUE, quadpts = NaN,
                     rotate = 'varimax', Target = NaN, SE = FALSE, verbose = TRUE,
                     SEtol = .0001, grsm.block = NULL, D = 1, TOL = NULL,
                     rsm.block = NULL, calcNull = TRUE, BFACTOR = FALSE,
                     technical = list(), use = 'pairwise.complete.obs',
                     SE.type = 'MHRM', large = NULL, accelerate = TRUE, empiricalhist = FALSE,
                     ...)
{
    opts <- list()
    if(method == 'MHRM' && SE.type == 'SEM') SE.type <- 'MHRM'
    D <- 1
    opts$method = method
    opts$draws = draws
    opts$calcLL = calcLL
    opts$quadpts = quadpts
    opts$rotate = rotate
    opts$Target = Target
    opts$SE = SE
    opts$SE.type = SE.type
    opts$verbose = verbose
    opts$SEtol = ifelse(is.null(technical$SEtol), .0001, technical$SEtol)
    opts$grsm.block = grsm.block
    opts$D = D
    opts$rsm.block = rsm.block
    opts$calcNull = calcNull
    opts$customPriorFun = technical$customPriorFun
    opts$BFACTOR = BFACTOR
    opts$accelerate = accelerate
    opts$TOL <- ifelse(is.null(TOL), if(method == 'EM') 1e-4 else 1e-3, TOL)
    if(SE.type == 'SEM' && SE){
        opts$accelerate <- FALSE
        if(is.null(TOL)) opts$TOL <- 1e-6
        if(is.null(technical$NCYCLES)) technical$NCYCLES <- 1000L
        opts$SEtol <- ifelse(is.null(technical$SEtol), .001, technical$SEtol)
    }
    if(BFACTOR && is.nan(quadpts)) opts$quadpts <- 21L
    if(is.null(technical$symmetric_SEM)) technical$symmetric_SEM <- TRUE
    opts$technical <- technical
    opts$use <- use
    opts$MAXQUAD <- ifelse(is.null(technical$MAXQUAD), 10000L, technical$MAXQUAD)
    opts$NCYCLES <- ifelse(is.null(technical$NCYCLES), 2000L, technical$NCYCLES)
    if(opts$method == 'EM')
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
    if(is.null(technical$Moptim)){
        opts$Moptim <- if(method == 'EM') 'BFGS' else 'NR'
    } else {
        opts$Moptim <- technical$Moptim
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

loadESTIMATEinfo <- function(info, ESTIMATE, constrain){
    longpars <- ESTIMATE$longpars
    pars <- ESTIMATE$pars
    ngroups <- length(pars)
    J <- length(pars[[1L]]) - 1L
    info <- nameInfoMatrix(info=info, correction=ESTIMATE$correction, L=ESTIMATE$L,
                           npars=length(longpars))
    isna <- rowSums(is.na(info)) > 0L
    info <- info[!isna, !isna, drop=FALSE]
    ESTIMATE$info <- info
    acov <- try(solve(info), TRUE)
    if(is(acov, 'try-error')){
        warning('Could not invert information matrix; model likely is not identified.')
        ESTIMATE$fail_invert_info <- TRUE
        return(ESTIMATE)
    } else ESTIMATE$fail_invert_info <- FALSE
    SEtmp <- diag(solve(info))
    if(any(SEtmp < 0)){
        warning("Negative SEs set to NaN.\n")
        SEtmp[SEtmp < 0 ] <- NaN
    }
    SEtmp <- sqrt(SEtmp)
    SE <- rep(NA, length(longpars))
    SE[ESTIMATE$estindex_unique[!isna]] <- SEtmp
    index <- 1L:ncol(info)
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
        if(colnames(gframe) %in% covnames){
            between <- TRUE
        } else if(colnames(gframe) %in% itemcovnames){
            between <- FALSE
        } else stop('grouping variable not in itemdesign or covdata')
        if(between){
            gframe <- gframe[1L:N, , drop=FALSE]
            sframe <- sframe[1L:N, , drop=FALSE]
        } else {
            gframe <- itemdesign[, which(colnames(gframe) == itemcovnames), drop=FALSE]
            sframe <- itemdesign[, which(colnames(sframe) == itemcovnames), drop=FALSE]
        }
        matpar <- diag(ncol(gframe) + ncol(sframe))
        estmat <- lower.tri(matpar, diag=TRUE)
        ndim <- ncol(matpar)
        if(strsplit(f, '+')[[1L]][[1L]] == '-')
            estmat[lower.tri(estmat)] <- FALSE
        fn <- paste0('COV_', c(colnames(gframe), colnames(sframe)))
        FNCOV <- outer(fn, c(colnames(gframe), colnames(sframe)), FUN=paste, sep='_')
        par <- matpar[lower.tri(matpar, diag=TRUE)]
        est <- estmat[lower.tri(estmat, diag=TRUE)]
        names(par) <- names(est) <- FNCOV[lower.tri(FNCOV, diag=TRUE)]
        drawvals <- matrix(0, length(unique(gframe)[[1L]]), ndim,
                           dimnames=list(unique(gframe)[[1L]], NULL))
        mtch <- match(gframe[[1L]], rownames(drawvals))
        gdesign <- matrix(1, nrow(gframe), 1L, dimnames = list(NULL, colnames(gframe)))
        if(ncol(sframe) != 0L)
            gdesign <- cbind(model.matrix(as.formula(paste0('~',splt[1L])), sframe), gdesign)
        tmp <- matrix(-Inf, ndim, ndim)
        diag(tmp) <- 1e-4
        lbound <- tmp[lower.tri(tmp, diag=TRUE)]
        ret[[i]] <- new('RandomPars',
                        par=par,
                        est=est,
                        ndim=ndim,
                        lbound=lbound,
                        ubound=rep(Inf, length(par)),
                        gframe=gframe,
                        gdesign=gdesign,
                        between=between,
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

OffTerm <- function(random, J, N){
    ret <- matrix(0, N, J)
    for(i in 1L:length(random)){
        if(random[[i]]@between){
            tmp <- rowSums(random[[i]]@gdesign*random[[i]]@drawvals[random[[i]]@mtch, ,drop=FALSE])
            for(j in 1L:J) ret[,j] <- ret[,j] + tmp
        } else {
            tmp <- matrix(rowSums(random[[i]]@gdesign *
                          random[[i]]@drawvals[random[[i]]@mtch, ,drop=FALSE]), nrow(ret), J,
                          byrow=TRUE)
            ret <- ret + tmp
        }
    }
    ret
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

assignInformationMG <- function(object){
    J <- ncol(object@data)
    names <- colnames(object@information)
    spl_names <- strsplit(names, split="\\.")
    spl_names_par <- sapply(spl_names, function(x) x[1L])
    spl_names <- lapply(spl_names, 
                        function(x) as.numeric(x[-1L]))
    spl_names <- do.call(rbind, spl_names)
    for(g in 1L:length(object@cmods)){
        from <- object@cmods[[g]]@pars[[1L]]@parnum[1L]
        to <- object@cmods[[g]]@pars[[J+1L]]@parnum[length(
            object@cmods[[g]]@pars[[J+1L]]@parnum)]
        pick <- spl_names[,g] >= from & spl_names[,g] <= to
        tmp <- object@information[pick,pick]
        colnames(tmp) <- rownames(tmp) <- 
            paste(spl_names_par[pick], spl_names[pick,g], sep='.')
        object@cmods[[g]]@information <- tmp
    }
    object
}

mirtClusterEnv <- new.env()
mirtClusterEnv$ncores <- 1L

myApply <- function(X, MARGIN, FUN, ...){
    if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
        return(t(parallel::parApply(cl=mirtClusterEnv$MIRTCLUSTER, X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    } else {
        return(t(apply(X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    }
}

myLapply <- function(X, FUN, ...){
    if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
        return(parallel::parLapply(cl=mirtClusterEnv$MIRTCLUSTER, X=X, fun=FUN, ...))
    } else {
        return(lapply(X=X, FUN=FUN, ...))
    }
}

mySapply <- function(X, FUN, ...){
    if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
        return(t(parallel::parSapply(cl=mirtClusterEnv$MIRTCLUSTER, X=X, FUN=FUN, ...)))
    } else {
        return(t(sapply(X=X, FUN=FUN, ...)))
    }
}
