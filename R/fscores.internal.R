setMethod(
	f = "fscores.internal",
	signature = 'ExploratoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, theta_lim, MI, 
	                      returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only,
	                      full.scores.SE)
	{
	    #local functions for apply
	    MAP <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND){
	        estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars, patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, hessian=TRUE, 
                                CUSTOM.IND=CUSTOM.IND))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
	        SEest <- try(sqrt(diag(solve(estimate$hessian))))
	        if(is(SEest, 'try-error')) SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    ML <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND){
	        if(any(scores[ID, ] %in% c(-Inf, Inf)))
                return(c(scores[ID, ], rep(NA, ncol(scores))))
	        estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars,patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, ML=TRUE, hessian=TRUE, 
                                CUSTOM.IND=CUSTOM.IND))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
	        SEest <- try(sqrt(diag(solve(estimate$hessian))))
	        if(is(SEest, 'try-error')) SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    WLE <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND){
	        estimate <- try(nlm(gradnorm.WLE,scores[ID, ],pars=pars,patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, CUSTOM.IND=CUSTOM.IND))
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
        CUSTOM.IND <- object@CUSTOM.IND
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
                                             quadpts=quadpts, gp=gp, verbose=verbose, 
                                             CUSTOM.IND=CUSTOM.IND))
        if (is.null(quadpts)) quadpts <- ceiling(40/(nfact^1.5))
		theta <- as.matrix(seq(theta_lim[1L], theta_lim[2L], length.out=quadpts))
		fulldata <- object@data
		tabdata <- object@tabdatalong
		tabdata <- tabdata[ ,-ncol(tabdata), drop = FALSE]
		keep <- object@tabdata[,ncol(object@tabdata)] > 0L
		tabdata <- tabdata[keep, , drop=FALSE]
        USETABDATA <- TRUE
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)
        list_SEscores <- list_scores <- vector('list', MI)
        if(MI == 0) MI <- 1	    
        impute <- MI > 1
        opars <- pars
        if(impute){
            if(is(try(chol(object@information), silent=TRUE), 'try-error')){
                stop('Proper information matrix must be precomputed in model for MI estimation')
            } else {
                names <- colnames(object@information)
                imputenums <- as.numeric(sapply(names, function(x, split){
                    strsplit(x, split=split)[[1L]][2L]
                }, split='\\.'))
                covB <- solve(object@information)
            }
        }
        for(mi in 1L:MI){
            if(impute)
                pars <- imputePars(pars=opars, covB=covB, imputenums=imputenums, 
                                   constrain=object@constrain)
            if(nfact < 3 || method == 'EAP'){
                ThetaShort <- Theta <- thetaComb(theta,nfact)
                if(length(prodlist) > 0L)
                    Theta <- prodterms(Theta,prodlist)
                W <- mirt_dmvnorm(ThetaShort,gp$gmeans,gp$gcov)
                itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, 
                                              CUSTOM.IND=CUSTOM.IND)
                log_itemtrace <- log(itemtrace)
        	    tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=EAP, log_itemtrace=log_itemtrace,
                               tabdata=tabdata, ThetaShort=ThetaShort, W=W)
        	    scores <- tmp[ ,1:nfact, drop = FALSE]
        	    SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
            }
    		if(method == "EAP"){
                #do nothing
    		} else if(method == "MAP"){
                tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=MAP, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist, 
                               CUSTOM.IND=CUSTOM.IND)
                scores <- tmp[ ,1:nfact, drop = FALSE]
                SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
    		} else if(method == "ML"){
                isna <- apply(object@tabdata[,-ncol(object@tabdata)], 1L, 
                              function(x) sum(is.na(x)))[keep]
    			allzero <- (rowSums(tabdata[,itemloc[-length(itemloc)]]) + isna) == J
    			allmcat <- (rowSums(tabdata[,itemloc[-1L]-1L]) + isna) == J
    			scores[allmcat,] <- Inf
                SEscores[allmcat,] <- NA                
                scores[allzero,] <- -Inf
                SEscores[allzero,] <- NA
                tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=ML, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist,
                               CUSTOM.IND=CUSTOM.IND)
                scores <- tmp[ ,1:nfact, drop = FALSE]
                SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
    		} else if(method == 'WLE'){
                if(nfact > 1L)
                    stop('WLE method only supported for unidimensional models')
                itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                              CUSTOM.IND=CUSTOM.IND)
                tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=WLE, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist, 
                               CUSTOM.IND=CUSTOM.IND)
                scores <- tmp[ ,1:nfact, drop = FALSE]
                SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
            } else {
                stop('method not defined')
            }
    		colnames(scores) <- paste('F', 1:ncol(scores), sep='')
            if(impute){
                list_SEscores[[mi]] <- SEscores
                list_scores[[mi]] <- scores
            } 
        }
        if(impute){
            scores <- list_scores[[1L]]/MI
            Ubar <- list_SEscores[[1L]]^2 / MI
            for(i in 2L:MI){
                scores <- list_scores[[i]]/MI + scores
                Ubar <- list_SEscores[[i]]^2 / MI + Ubar
            }
            B <- matrix(0, nrow(scores), ncol(scores))
            for(i in 1L:MI)
                B <- B + (1 / (MI-1L)) * ((list_scores[[i]] - scores)^2)
            SEscores <- sqrt(Ubar + (1 + 1/MI) * B)
        }
		if (full.scores){
            if(USETABDATA){
                tabdata2 <- object@tabdatalong
                tabdata2 <- tabdata2[tabdata2[,ncol(tabdata2)] > 0L, -ncol(tabdata2)]
                stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
                sfulldata <- apply(object@fulldata, 1, paste, sep='', collapse = '/')
                scoremat <- scores[match(sfulldata, stabdata2), , drop = FALSE]
                SEscoremat <- SEscores[match(sfulldata, stabdata2), , drop = FALSE]
    			colnames(scoremat) <- colnames(scores)
    			colnames(SEscoremat) <- paste0('SE_',colnames(scores))
            } else {
                scoremat <- scores
                SEscoremat <- SEscores
                colnames(SEscoremat) <- paste0('SE_',colnames(scores))
            }
            if(full.scores.SE)
                scoremat <- cbind(scoremat, SEscoremat)
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
            ret <- cbind(object@tabdata[keep, ,drop=FALSE],scores,SEscores)
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
	                      quadpts = NULL, response.pattern = NULL, theta_lim, MI,
	                      returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only,
	                      full.scores.SE)
	{
        class(object) <- 'ExploratoryClass'
        ret <- fscores(object, rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, quadpts=quadpts,
                       response.pattern=response.pattern, returnER=returnER, verbose=verbose,
                       mean=gmean, cov=gcov, scores.only=scores.only, theta_lim=theta_lim, MI=MI,
                       full.scores.SE=full.scores.SE)
        return(ret)
	}
)

#------------------------------------------------------------------------------
setMethod(
    f = "fscores.internal",
    signature = 'MultipleGroupClass',
    definition = function(object, rotate = '', full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, theta_lim, MI,
                          returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only,
                          full.scores.SE)
    {
        cmods <- object@cmods
        ngroups <- length(cmods)
        for(g in 1L:ngroups)
            class(cmods[[g]]) <- 'ConfirmatoryClass'
        if(MI > 0){
            object <- assignInformationMG(object)
            cmods <- object@cmods
        }
        ret <- vector('list', length(cmods))
        for(g in 1L:ngroups)
            ret[[g]] <- fscores(cmods[[g]], rotate = 'CONFIRMATORY', full.scores=full.scores, method=method,
                           quadpts=quadpts, returnER=returnER, verbose=verbose, theta_lim=theta_lim,
                                mean=gmean[[g]], cov=gcov[[g]], scores.only=FALSE, MI=MI,
                           full.scores.SE=full.scores.SE)
        names(ret) <- object@groupNames
        if(full.scores){
            id <- c()
            fulldata <- matrix(NA, 1, ncol(ret[[1]]))
            for(g in 1L:ngroups){
                id <- c(id, rownames(ret[[g]]))
                fulldata <- rbind(fulldata, ret[[g]])
            }
            fulldata <- fulldata[-1L, ]
            fulldata <- data.frame(id=as.numeric(id), fulldata)
            ret <- fulldata[order(fulldata$id), ]
            ret <- ret[ ,-1L]
            if(scores.only)
                ret <- ret[ ,!(colnames(ret) %in% colnames(object@data)), drop=FALSE]
        }
        if(is.data.frame(ret))
            ret <- as.matrix(ret)
        return(ret)
    }
)

# MAP scoring for mirt
MAP.mirt <- function(Theta, pars, patdata, itemloc, gp, prodlist, CUSTOM.IND, ML=FALSE)
{
    Theta <- matrix(Theta, nrow=1L)
    ThetaShort <- Theta
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta,prodlist)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
    L <- sum(log(itemtrace)[as.logical(patdata)])
    prior <- mirt_dmvnorm(ThetaShort, gp$gmeans, gp$gcov)
    L <- ifelse(ML, -L, (-1)*(L + log(prior)))
    L
}

gradnorm.WLE <- function(Theta, pars, patdata, itemloc, gp, prodlist, CUSTOM.IND){
    Theta <- matrix(Theta, nrow=1L)
    ThetaShort <- Theta
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta,prodlist)
    nfact <- ncol(Theta)
    itemtrace <- matrix(0, ncol=length(patdata), nrow=nrow(Theta))
    dP <- d2P <- vector('list', nfact)
    for(i in 1L:nfact)
        dP[[i]] <- d2P[[i]] <- itemtrace
    I <- numeric(1)
    dW <- dL <- numeric(nfact)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
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

EAPsum <- function(x, full.scores = FALSE, quadpts = NULL, S_X2 = FALSE, gp, verbose, CUSTOM.IND){
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
    prior <- mirt_dmvnorm(Theta,gp$gmeans,gp$gcov)
    prior <- prior/sum(prior)
    pars <- x@pars
    K <- x@K
    J <- length(K)
    itemloc <- x@itemloc
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, 
                                  CUSTOM.IND=CUSTOM.IND)
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
