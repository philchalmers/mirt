setMethod(
	f = "fscores.internal",
	signature = 'ExploratoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, 
	                      returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only)
	{
	    #local functions for apply
	    MAP <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist){
	        estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars, patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, hessian=TRUE))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
	        SEest <- try(sqrt(diag(solve(estimate$hessian))))
	        if(is(SEest, 'try-error')) SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    ML <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist){
	        if(any(scores[ID, ] %in% c(-Inf, Inf)))
                return(c(scores[ID, ], rep(NA, ncol(scores))))
	        estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars,patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, ML=TRUE, hessian=TRUE))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
	        SEest <- try(sqrt(diag(solve(estimate$hessian))))
	        if(is(SEest, 'try-error')) SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    WLE <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist){
	        estimate <- try(nlm(gradnorm.WLE,scores[ID, ],pars=pars,patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
            TI <- 0
            for(i in 1L:(length(itemloc)-1L))
                TI <- TI + iteminfo(pars[[i]], Theta=estimate$estimate)
	        SEest <- 1 / sqrt(TI)
	        return(c(estimate$estimate, SEest))
	    }
	    EAP <- function(ID, log_itemtrace, tabdata, ThetaShort, W){
	        L <- rowSums(log_itemtrace[ ,as.logical(tabdata[ID,]), drop = FALSE])
	        thetas <- colSums(ThetaShort * exp(L) * W / sum(exp(L) * W))
	        SE <- sqrt(colSums(t((t(ThetaShort) - thetas))^2 * exp(L) * W / sum(exp(L) * W)))
	        return(c(thetas, SE))
	    }

        if(!is.null(response.pattern)){
            drop <- FALSE
            if(!is.matrix(response.pattern)){
                response.pattern <- rbind(response.pattern, response.pattern)
                drop <- TRUE
            }
            nfact <- object@nfact
            sv <- mod2values(object)
            sv$est <- FALSE
            mins <- apply(object@data, 2, min, na.rm=TRUE)
            response.pattern <- response.pattern - matrix(mins, nrow(response.pattern),
                                                          ncol(response.pattern), byrow=TRUE)
            colnames(response.pattern) <- colnames(object@data)
            newmod <- mirt(response.pattern, nfact, itemtype = object@itemtype, pars=sv, calcNull=FALSE,
                           technical=list(customK=object@K))
            ret <- fscores(newmod, rotate=rotate, full.scores=full.scores, scores.only=scores.only,
                           method=method, quadpts=quadpts, verbose=FALSE,
                           response.pattern=NULL)
            if(!scores.only || !full.scores)
                ret[,1L:ncol(response.pattern)] <- ret[,1L:ncol(response.pattern)] +
                    matrix(mins, nrow(ret), ncol(response.pattern), byrow=TRUE)
            if(drop){
                if(full.scores){
                    ret <- ret[-1L, , drop=FALSE]
                } else ret[1L, ncol(response.pattern)+1L] <- 1L
            }
            return(ret)
        }
        pars <- object@pars
		K <- object@K
        J <- length(K)
        prodlist <- attr(pars, 'prodlist')
        nfact <- object@nfact
        nLambdas <- object@nfact
        itemloc <- object@itemloc
        gp <- ExtractGroupPars(object@pars[[length(itemloc)]])
        if(rotate != 'CONFIRMATORY'){
            so <- summary(object, rotate = rotate, verbose = FALSE)
            a <- rotateLambdas(so)
            for(i in 1L:J)
                pars[[i]]@par[1L:nfact] <- a[i, ]
            gp$gmeans <- rep(0, nfact)
            gp$gcov <- so$fcor
        }
        if(!is.null(gmean)) gp$gmeans <- gmean
        if(!is.null(gcov)) gp$gcov <- gcov
        if(method == 'EAPsum') return(EAPsum(object, full.scores=full.scores,
                                             quadpts=quadpts, gp=gp, verbose=verbose))
        if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))
		theta <- as.matrix(seq(-4,4,length.out = quadpts))
		ThetaShort <- Theta <- thetaComb(theta,nfact)
        if(length(prodlist) > 0L)
            Theta <- prodterms(Theta,prodlist)
		fulldata <- object@data
		tabdata <- object@tabdatalong
		tabdata <- tabdata[ ,-ncol(tabdata), drop = FALSE]
        USETABDATA <- TRUE
        if(full.scores && nrow(tabdata) > nrow(fulldata)/5){
            USETABDATA <- FALSE
            tabdata <- object@fulldata
        }
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)
		W <- mvtnorm::dmvnorm(ThetaShort,gp$gmeans,gp$gcov)
		W <- W/sum(W)
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
	    log_itemtrace <- log(itemtrace)
	    if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
	        tmp <- t(parallel::parApply(cl=mirtClusterEnv$MIRTCLUSTER, matrix(1:nrow(scores)), 1, EAP,
	                                    log_itemtrace=log_itemtrace, tabdata=tabdata,
                                        ThetaShort=ThetaShort, W=W))
	    } else {
	        tmp <- t(apply(matrix(1:nrow(scores)), 1, EAP, log_itemtrace=log_itemtrace, tabdata=tabdata,
	                       ThetaShort=ThetaShort, W=W))
	    }
	    scores <- tmp[ ,1:nfact, drop = FALSE]
	    SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
		if(method == "MAP"){
            if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
                tmp <- t(parallel::parApply(cl=mirtClusterEnv$MIRTCLUSTER, matrix(1:nrow(scores)), 1, MAP,
                                     scores=scores, pars=pars, tabdata=tabdata, itemloc=itemloc,
                                     gp=gp, prodlist=prodlist))
            } else {
                tmp <- t(apply(matrix(1:nrow(scores)), 1, MAP, scores=scores, pars=pars,
                             tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist))
            }
            scores <- tmp[ ,1:nfact, drop = FALSE]
            SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
		}
		if(method == "ML"){
            tabdata2 <- object@tabdata[,-ncol(object@tabdata)]
			tmp2 <- tabdata[,itemloc[-1L] - 1L, drop = FALSE]
            tmp2[is.na(tabdata2)] <- 1
			scores[rowSums(tmp2) == J,] <- Inf
            SEscores[rowSums(tmp2) == J,] <- NA
            tmp2 <- tabdata[,itemloc[-length(itemloc)], drop = FALSE]
            tmp2[is.na(tabdata2)] <- 1
            scores[rowSums(tmp2) == J,] <- -Inf
            SEscores[rowSums(tmp2) == J,] <- NA
			SEscores[is.na(scores[,1L]), ] <- rep(NA, nfact)
            if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
                tmp <- t(parallel::parApply(cl=mirtClusterEnv$MIRTCLUSTER, matrix(1:nrow(scores)), 1, ML,
                                            scores=scores, pars=pars, tabdata=tabdata, itemloc=itemloc,
                                            gp=gp, prodlist=prodlist))
            } else {
                tmp <- t(apply(matrix(1:nrow(scores)), 1, ML, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist))
            }
            scores <- tmp[ ,1:nfact, drop = FALSE]
            SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
		}
        if(method == 'WLE'){
            if(nfact > 1L)
                stop('WLE method only supported for unidimensional models')
            itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
            if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
                tmp <- t(parallel::parApply(cl=mirtClusterEnv$MIRTCLUSTER, matrix(1:nrow(scores)), 1, WLE,
                                            scores=scores, pars=pars, tabdata=tabdata, itemloc=itemloc,
                                            gp=gp, prodlist=prodlist))
            } else {
                tmp <- t(apply(matrix(1:nrow(scores)), 1, WLE, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist))
            }
            scores <- tmp[ ,1:nfact, drop = FALSE]
            SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
        }
		colnames(scores) <- paste('F', 1:ncol(scores), sep='')
		if (full.scores){
            if(USETABDATA){
                tabdata2 <- object@tabdatalong
                tabdata2 <- tabdata2[,-ncol(tabdata2)]
                stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
                sfulldata <- apply(object@fulldata, 1, paste, sep='', collapse = '/')
                scoremat <- scores[match(sfulldata, stabdata2), , drop = FALSE]
    			colnames(scoremat) <- colnames(scores)
            } else {
                scoremat <- scores
            }
            if(scores.only) return(scoremat)
			else return(cbind(fulldata,scoremat))
		} else {
            r <- object@tabdata[,ncol(object@tabdata)]
            T <- E <- matrix(NA, 1, ncol(scores))
            for(i in 1L:nrow(scores)){
                if(any(scores[i, ] %in% c(Inf, -Inf))) next
                T <- rbind(T, matrix(rep(scores[i, ], r[i]), ncol=ncol(scores), byrow = TRUE))
                E <- rbind(E, matrix(rep(SEscores[i, ], r[i]), ncol=ncol(scores), byrow = TRUE))
            }
            T <- na.omit(T)
            E <- na.omit(E)
            reliability <- diag(var(T)) / (diag(var(T)) + colMeans(E^2))
            names(reliability) <- colnames(scores)
            if(returnER) return(reliability)
			if(verbose){
                cat("\nMethod: ", method)
                cat("\n\nEmpirical Reliability:\n")
                print(round(reliability, 4L))
			}
			colnames(SEscores) <- paste('SE_', colnames(scores), sep='')
            ret <- cbind(object@tabdata,scores,SEscores)
            if(nrow(ret) > 1L) ret <- ret[do.call(order, as.data.frame(ret[,1L:J])), ]
			return(ret)
		}
	}
)

#------------------------------------------------------------------------------
setMethod(
	f = "fscores.internal",
	signature = 'ConfirmatoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP",
	                      quadpts = NULL, response.pattern = NULL, 
	                      returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only)
	{
        class(object) <- 'ExploratoryClass'
        ret <- fscores(object, rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, quadpts=quadpts,
                       response.pattern=response.pattern, returnER=returnER, verbose=verbose,
                       mean=gmean, cov=gcov, scores.only=scores.only)
        return(ret)
	}
)

#------------------------------------------------------------------------------
setMethod(
    f = "fscores.internal",
    signature = 'MultipleGroupClass',
    definition = function(object, rotate = '', full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, 
                          returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only)
    {
        cmods <- object@cmods
        ngroups <- length(cmods)
        for(g in 1L:ngroups)
            class(cmods[[g]]) <- 'ConfirmatoryClass'
        ret <- vector('list', length(cmods))
        for(g in 1L:ngroups)
            ret[[g]] <- fscores(cmods[[g]], rotate = 'CONFIRMATORY', full.scores=full.scores, method=method,
                           quadpts=quadpts, returnER=returnER, verbose=verbose,
                                mean=gmean[[g]], cov=gcov[[g]], scores.only=scores.only)
        names(ret) <- object@groupNames
        if(full.scores){
            id <- c()
            fulldata <- matrix(NA, 1, ncol(ret[[1]]))
            for(g in 1L:ngroups){
                id <- c(id, rownames(ret[[g]]))
                fulldata <- rbind(fulldata, ret[[g]])
            }
            if(!scores.only){
                fulldata <- fulldata[-1L, ]
                fulldata <- data.frame(id=as.numeric(id), fulldata)
                ret <- fulldata[order(fulldata$id), ]
                ret <- ret[ ,-1L]
            } else {
                ret <- fulldata[-1L, , drop = FALSE]
            }
        }
        return(ret)
    }
)

# MAP scoring for mirt
MAP.mirt <- function(Theta, pars, patdata, itemloc, gp, prodlist, ML=FALSE)
{
    ThetaShort <- Theta
    Theta <- matrix(Theta, nrow=1)
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta,prodlist)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
    L <- sum(log(itemtrace)[as.logical(patdata)])
    prior <- mvtnorm::dmvnorm(ThetaShort, gp$gmeans, gp$gcov)
    L <- ifelse(ML, -L, (-1)*(L + log(prior)))
    L
}

gradnorm.WLE <- function(Theta, pars, patdata, itemloc, gp, prodlist){
    ThetaShort <- Theta
    Theta <- matrix(Theta, nrow=1)
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta,prodlist)
    nfact <- ncol(Theta)
    itemtrace <- matrix(0, ncol=length(patdata), nrow=nrow(Theta))
    dP <- d2P <- vector('list', nfact)
    for(i in 1L:nfact)
        dP[[i]] <- d2P[[i]] <- itemtrace
    I <- numeric(1)
    dW <- dL <- numeric(nfact)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
    for (i in 1L:(length(itemloc)-1L)){
        deriv <- DerivTheta(x=pars[[i]], Theta=Theta)
        for(k in 1L:nfact){
            dPitem <- d2Pitem <- matrix(0, 1, length(deriv[[1L]]))
            for(j in 1L:length(deriv[[1L]])){
                dPitem[1, j] <- deriv$grad[[j]][ ,k]
                d2Pitem[1, j] <- deriv$hess[[j]][ ,k]
            }
            dP[[k]][ ,itemloc[i]:(itemloc[i+1L] - 1L)] <- dPitem
            d2P[[k]][ ,itemloc[i]:(itemloc[i+1L] - 1L)] <- d2Pitem
            dW[k] <- dW[k] + sum(dPitem * d2Pitem / itemtrace[ ,itemloc[i]:(itemloc[i+1L] - 1L)])
        }
        I <- I + iteminfo(x=pars[[i]], Theta=Theta)
    }
    dW <- 1/(2*I^2) * dW
    for(i in 1L:nfact)
        dL[i] <- sum(patdata * dP[[i]] / itemtrace)
    grad <- dL + dW*I
    MIN <- sum(grad^2)
    MIN
}

EAPsum <- function(x, full.scores = FALSE, quadpts = NULL, S_X2 = FALSE, gp, verbose){
    calcL1 <- function(itemtrace, K, itemloc){
        J <- length(K)
        L0 <- L1 <- matrix(1, sum(K-1L) + 1L, ncol(itemtrace))
        L0[1L:K[1L], ] <- itemtrace[1:K[1], ]
        nstar <- K[1L] + K[2L] - 3L
        Sum.Scores <- 1L:nrow(L0)-1L
        MAX.Scores <- cumsum(K-1L)
        for(i in 1L:(J-1L)){
            T <- itemtrace[itemloc[i+1L]:(itemloc[i+2L] - 1L), ]
            L1[1L, ] <- L0[1L, ] * T[1L, ]
            for(j in 1L:nstar+1L){
                sums <- 0
                for(k in 1L:K[i+1L]-1L)
                    if(Sum.Scores[j] >= k && (MAX.Scores[i] + k) >= Sum.Scores[j])
                        sums <- sums + L0[j - k, ] * T[1L + k, ]
                L1[j, ] <- sums
            }
            L1[j+1L, ] <- L0[j - k + 1L, ] * T[nrow(T), ]
            L0 <- L1
            if(i != J) nstar <- nstar + K[i+2L] - 1L
        }
        list(L1=L1, Sum.Scores=Sum.Scores)
    }
    if(x@nfact > 1L) stop('EAP sum score method only is applicable to unidimensional models')
    if(is.null(quadpts)) quadpts <- 40
    Theta <- as.matrix(seq(-4,4,length.out = quadpts))
    prior <- mvtnorm::dmvnorm(Theta,gp$gmeans,gp$gcov)
    prior <- prior/sum(prior)
    pars <- x@pars
    K <- x@K
    J <- length(K)
    itemloc <- x@itemloc
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
    itemtrace <- t(itemtrace)
    tmp <- calcL1(itemtrace=itemtrace, K=K, itemloc=itemloc)
    L1 <- tmp$L1
    Sum.Scores <- tmp$Sum.Scores
    if(S_X2){
        L1total <- L1 %*% prior
        Elist <- vector('list', J)
        for(i in 1L:J){
            KK <- K[-i]
            T <- itemtrace[c(itemloc[i]:(itemloc[i+1L]-1L)), ]
            itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i+1L]-1L)), ]
            if(i != J){
                itemloc2 <- itemloc[-i]
                itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
            } else itemloc2 <- itemloc[-(J+1)]
            tmp <- calcL1(itemtrace=itemtrace2, K=KK, itemloc=itemloc2)
            E <- matrix(NA, nrow(L1total), nrow(T))
            for(j in 1L:(nrow(T)))
                E[1L:nrow(tmp$L1)+j-1L,j] <- tmp$L1 %*% (T[j,] * prior) /
                    L1total[1L:nrow(tmp$L1)+j-1L, ]
            Elist[[i]] <- E[-c(1L, nrow(E)), ]
        }
        return(Elist)
    }
    thetas <- SEthetas <- numeric(nrow(L1))
    for(i in 1L:length(thetas)){
        thetas[i] <- sum(Theta * L1[i, ] * prior / sum(L1[i,] * prior))
        SEthetas[i] <- sqrt(sum((Theta - thetas[i])^2 * L1[i, ] * prior / sum(L1[i,] * prior)))
    }
    ret <- data.frame(Sum.Scores=Sum.Scores, Theta=thetas, SE.Theta=SEthetas)
    rownames(ret) <- ret$Sum.Scores
    if(full.scores){
        if(any(is.na(x@data))) stop('Full scores requires a complete dataset (no N\'s)')
        dat <- x@data
        adj <- apply(dat, 2, min)
        if(any(adj > 0L)) message('Data adjusted so that every item has a lowest score of 0')
        dat <- t(t(dat) - adj)
        scores <- rowSums(dat)
        EAPscores <- ret$Theta[match(scores, ret$Sum.Scores)]
        ret <- data.frame(Sum.Scores=scores, Theta=EAPscores)
    } else {
        dat <- x@data
        E <- L1 %*% prior * nrow(dat)
        adj <- apply(dat, 2, min)
        dat <- t(t(dat) - adj)
        Otmp <- matrix(table(sort(rowSums(dat))))
        got <- as.numeric(names(table(sort(rowSums(dat)))))
        if(min(got) == 0) got <- got + 1
        O <- matrix(0, nrow(E), 1)
        O[got, 1] <- Otmp
        keep <- O != 0
        ret$observed <- O
        ret$expected <- E
        O <- O[keep]
        E <- E[keep]
        df <- length(ret$observed) - 1
        X2 <- sum((ret$observed - ret$expected)^2 / ret$expected)
        G2 <- 2 * sum(O * log(O/E))
        attr(ret, 'fit') <- data.frame(df=df, X2=X2, p.X2 = pchisq(X2, df, lower.tail=FALSE),
                                       G2=G2, p.G2 = pchisq(G2, df, lower.tail=FALSE))
        if(verbose){
            print(attr(ret, 'fit'))
            cat('\n')
        }
    }
    ret
}
