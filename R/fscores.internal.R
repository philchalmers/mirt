setMethod(
	f = "fscores.internal",
	signature = 'ExploratoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.vector = NULL, degrees = NULL,
	                      returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only)
	{
        if(!is.null(response.vector)){
            if(!is.matrix(response.vector)) response.vector <- matrix(response.vector, nrow = 1)
            v <- response.vector
            newdata <- rbind(object@data, v)
            nfact <- object@nfact
            sv <- mod2values(object)
            sv$est <- FALSE
            newmod <- mirt(newdata, nfact, pars=sv, calcNull=FALSE)
            tmptabdata <- t(newmod@tabdata[, 1L:length(v)])
            tmptabdata[is.na(tmptabdata)] <- 9999
            v[is.na(v)] <- 9999
            ind <- which(colSums(tmptabdata == as.numeric(v)) == length(v))
            newmod@tabdata <- newmod@tabdata[ind, , drop = FALSE]
            newmod@tabdatalong <- newmod@tabdatalong[ind, , drop = FALSE]
            ret <- fscores(newmod, rotate=rotate, full.scores=FALSE, method=method,
                           quadpts=quadpts, verbose=FALSE, degrees=degrees, response.vector=NULL)
            ret <- ret[,colnames(ret) != 'Freq']
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
                                             quadpts=quadpts, gp=gp))
        if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))
		theta <- as.matrix(seq(-4,4,length.out = quadpts))
		ThetaShort <- Theta <- thetaComb(theta,nfact)
        if(length(prodlist) > 0L)
            Theta <- prodterms(Theta,prodlist)
		fulldata <- object@data
		tabdata <- object@tabdatalong
		tabdata <- tabdata[ ,-ncol(tabdata), drop = FALSE]
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)
		W <- mvtnorm::dmvnorm(ThetaShort,gp$gmeans,gp$gcov)
		W <- W/sum(W)
        itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)
		for (i in 1L:nrow(tabdata)){
			L <- rowSums(log(itemtrace)[ ,as.logical(tabdata[i,]), drop = FALSE])
			thetas <- colSums(ThetaShort * exp(L) * W / sum(exp(L) * W))
			SE <- sqrt(colSums(t((t(ThetaShort) - thetas))^2 * exp(L) * W / sum(exp(L) * W)))
			scores[i, ] <- thetas
			SEscores[i, ] <- SE
		}
		if(method == "MAP"){
			for (i in 1L:nrow(scores)){
				tmp <- scores[i, ]
                estimate <- try(nlm(MAP.mirt,tmp,pars=pars, patdata=tabdata[i, ],
                                itemloc=itemloc, gp=gp, prodlist=prodlist, hessian=TRUE))
                if(is(estimate, 'try-error')) {
                    scores[i, ] <- SEscores[i, ] <- NA
                    next
                }
				scores[i, ] <- estimate$estimate
				SEest <- try(sqrt(diag(solve(estimate$hessian))))
				if(is(SEest, 'try-error')) SEest <- rep(NA, nfact)
				SEscores[i, ] <- SEest
			}
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
			for (i in 1L:nrow(scores)){
				if(any(scores[i, ] %in% c(-Inf, Inf))) next
				Theta <- scores[i, ]
				estimate <- try(nlm(MAP.mirt,Theta,pars=pars,patdata=tabdata[i, ],
				    itemloc=itemloc, gp=gp, prodlist=prodlist, ML=TRUE, hessian = TRUE))
				if(is(estimate, 'try-error')) {
				    scores[i, ] <- SEscores[i, ] <- NA
				    next
				}
				scores[i, ] <- estimate$estimate
                SEest <- try(sqrt(diag(solve(estimate$hessian))))
                if(is(SEest, 'try-error')) SEest <- rep(NA, nfact)
				SEscores[i, ] <- SEest
			}
		}
        if(method == 'WLE'){
            tabdata2 <- object@tabdata[,-ncol(object@tabdata)]
            tmp2 <- tabdata[,itemloc[-1L] - 1L, drop = FALSE]
            tmp2[is.na(tabdata2)] <- 1
            scores[rowSums(tmp2) == J,] <- Inf
            tmp2 <- tabdata[,itemloc[-length(itemloc)], drop = FALSE]
            tmp2[is.na(tabdata2)] <- 1
            scores[rowSums(tmp2) == J,] <- -Inf
            SEscores[is.na(scores[,1L]), ] <- rep(NA, nfact)
            SEscores[!is.na(SEscores)] <- NA
            for (i in 1L:nrow(scores)){
                if(any(scores[i, ] %in% c(-Inf, Inf))) next
                Theta <- scores[i, ]
                estimate <- try(nlm(gradnorm.WLE,Theta,pars=pars,patdata=tabdata[i, ],
                                itemloc=itemloc, gp=gp, prodlist=prodlist, degrees=degrees,
                                hessian = TRUE))
                if(is(estimate, 'try-error')) {
                    scores[i, ] <- NA
                    next
                }
                scores[i, ] <- estimate$estimate
            }
        }
		colnames(scores) <- paste('F', 1:ncol(scores), sep='')
		if (full.scores){
            tabdata2 <- object@tabdatalong
            tabdata2 <- tabdata2[,-ncol(tabdata2)]
            stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
            sfulldata <- apply(object@fulldata, 1, paste, sep='', collapse = '/')
            scoremat <- scores[match(sfulldata, stabdata2), , drop = FALSE]
			colnames(scoremat) <- colnames(scores)
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
	                      quadpts = NULL, response.vector = NULL, degrees = NULL,
	                      returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only)
	{
        class(object) <- 'ExploratoryClass'
        ret <- fscores(object, rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, quadpts=quadpts,
                       response.vector=response.vector, degrees=degrees, returnER=returnER, verbose=verbose, 
                       mean=gmean, cov=gcov, scores.only=scores.only)
        return(ret)
	}
)

#------------------------------------------------------------------------------
setMethod(
    f = "fscores.internal",
    signature = 'MultipleGroupClass',
    definition = function(object, rotate = '', full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.vector = NULL, degrees = NULL,
                          returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only)
    {
        cmods <- object@cmods
        ngroups <- length(cmods)
        for(g in 1L:ngroups)
            class(cmods[[g]]) <- 'ConfirmatoryClass'
        ret <- vector('list', length(cmods))
        for(g in 1L:ngroups)
            ret[[g]] <- fscores(cmods[[g]], rotate = 'CONFIRMATORY', full.scores=full.scores, method=method,
                           quadpts=quadpts, degrees=degrees, returnER=returnER, verbose=verbose, 
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
    itemtrace <- matrix(0, ncol=length(patdata), nrow=nrow(Theta))
    for (i in 1L:(length(itemloc)-1L))
        itemtrace[ ,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(x=pars[[i]], Theta=Theta)
    itemtrace[itemtrace < 1e-8] <- 1e-8
    L <- sum(log(itemtrace)[as.logical(patdata)])
    prior <- mvtnorm::dmvnorm(ThetaShort, gp$gmeans, gp$gcov)
    L <- ifelse(ML, -L, (-1)*(L + log(prior)))
    L
}

gradnorm.WLE <- function(Theta, pars, patdata, itemloc, gp, prodlist, degrees){
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
    for (i in 1L:(length(itemloc)-1L)){
        itemtrace[ ,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(x=pars[[i]], Theta=Theta)
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
        I <- I + iteminfo(x=pars[[i]], Theta=Theta, degrees=degrees)
    }
    dW <- 1/(2*I^2) * dW
    for(i in 1L:nfact)
        dL[i] <- sum(patdata * dP[[i]] / itemtrace)
    grad <- dL - dW*I
    MIN <- sum(grad^2)
    MIN
}

EAPsum <- function(x, full.scores = FALSE, quadpts = NULL, S_X2 = FALSE, gp){
    calcL1 <- function(itemtrace, K){
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
    itemtrace <- matrix(0, ncol=ncol(x@tabdatalong)-1, nrow=nrow(Theta))
    for (i in 1L:J)
        itemtrace[ ,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(x=pars[[i]], Theta=Theta)
    itemtrace <- t(itemtrace)
    tmp <- calcL1(itemtrace=itemtrace, K=K)
    L1 <- tmp$L1
    Sum.Scores <- tmp$Sum.Scores
    if(S_X2){
        L1total <- L1 %*% prior
        Elist <- vector('list', J)
        for(i in 1L:J){
            KK <- K[-i]
            T <- itemtrace[c(itemloc[i]:(itemloc[i+1L]-1L)), ]
            itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i+1L]-1L)), ]
            tmp <- calcL1(itemtrace=itemtrace2, K=KK)
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
    }
    ret
}
