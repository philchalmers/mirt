setMethod(
	f = "fscores.internal",
	signature = 'ExploratoryClass',
	definition = function(object, rotate = '', full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, theta_lim, MI, 
	                      returnER = FALSE, verbose = TRUE, gmean, gcov, scores.only,
	                      full.scores.SE, return.acov = FALSE)
	{
	    #local functions for apply
	    MAP <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
	                    hessian, return.acov = FALSE){
	        estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars, patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, hessian=hessian, 
                                CUSTOM.IND=CUSTOM.IND))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
            if(hessian){
                vcov <- try(solve(estimate$hessian))
                if(return.acov) return(vcov)
    	        SEest <- try(sqrt(diag(vcov)))
    	        if(is(SEest, 'try-error')) SEest <- rep(NA, ncol(scores))
            } else SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    ML <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
	                   hessian, return.acov = FALSE){
	        if(any(scores[ID, ] %in% c(-Inf, Inf)))
                return(c(scores[ID, ], rep(NA, ncol(scores))))
	        estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars,patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, ML=TRUE, hessian=hessian, 
                                CUSTOM.IND=CUSTOM.IND))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
            if(hessian){
    	        vcov <- try(solve(estimate$hessian))
    	        if(return.acov) return(vcov)
    	        SEest <- try(sqrt(diag(vcov)))
    	        if(is(SEest, 'try-error')) SEest <- rep(NA, ncol(scores))
            } else SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    WLE <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
                        hessian, return.acov = FALSE){
	        estimate <- try(nlm(gradnorm.WLE,scores[ID, ],pars=pars,patdata=tabdata[ID, ], hessian=FALSE,
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, CUSTOM.IND=CUSTOM.IND))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
            if(hessian){
                TI <- 0
                for(i in 1L:(length(itemloc)-1L))
                    TI <- TI + iteminfo(pars[[i]], Theta=estimate$estimate)
                if(return.acov) return(1/TI)
    	        SEest <- 1 / sqrt(TI)
            } else SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    EAP <- function(ID, log_itemtrace, tabdata, ThetaShort, W, hessian, return.acov = FALSE){
            nfact <- ncol(ThetaShort)
	        L <- rowSums(log_itemtrace[ ,as.logical(tabdata[ID,]), drop = FALSE])
	        thetas <- colSums(ThetaShort * exp(L) * W / sum(exp(L) * W))
            if(hessian){
    	        thetadif <- t((t(ThetaShort) - thetas))
                Thetaprod <- matrix(0, nrow(ThetaShort), nfact * (nfact + 1L)/2L)
                ind <- 1L
                for(i in 1L:nfact){
                    for(j in 1L:nfact){
                        if(i <= j){
                            Thetaprod[,ind] <- thetadif[,i] * thetadif[,j]
                            ind <- ind + 1L
                        }
                    }
                }
                vcov <- matrix(0, nfact, nfact)
                vcov[lower.tri(vcov, TRUE)] <- colSums(Thetaprod * exp(L) * W / sum(exp(L) * W))
                if(nfact > 1L) vcov <- vcov + t(vcov) - diag(diag(vcov))
                if(return.acov) return(vcov)
    	        SE <- sqrt(diag(vcov))
            } else SE <- rep(NA, nfact)
	        return(c(thetas, SE))
	    }
        
        if(return.acov && MI != 0)
            stop('simultaneous impute and return.acov option not supported')
	    if(return.acov && returnER)
	        stop('simultaneous returnER and return.acov option not supported')
        if(!is.null(response.pattern)){
            if(is.data.frame(response.pattern))
                response.pattern <- as.matrix(response.pattern)
            if(!is.matrix(response.pattern))
                response.pattern <- matrix(response.pattern, 1L)
            nfact <- object@nfact
            sv <- mod2values(object)
            mins <- object@Data$mins
            if(!all(mins == 0L))
                response.pattern <- response.pattern - matrix(mins, nrow(response.pattern),
                                                          ncol(response.pattern), byrow=TRUE)
            colnames(response.pattern) <- colnames(object@Data$data)
            large <- suppressWarnings(mirt(response.pattern, nfact, technical=list(customK=object@K), 
                          large=TRUE))
            newmod <- object
            newmod@Data <- list(data=response.pattern, tabdata=large$tabdata2, 
                               tabdatalong=large$tabdata, Freq=large$Freq)
            ret <- fscores(newmod, rotate=rotate, full.scores=TRUE, scores.only=FALSE,
                           method=method, quadpts=quadpts, verbose=FALSE, full.scores.SE=TRUE,
                           response.pattern=NULL, return.acov=return.acov)
            if(return.acov) return(ret)
            if(!all(mins == 0L))
                ret[,1L:ncol(response.pattern)] <- ret[,1L:ncol(response.pattern)] +
                    matrix(mins, nrow(ret), ncol(response.pattern), byrow=TRUE)
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
                                             CUSTOM.IND=CUSTOM.IND, theta_lim=theta_lim))
		theta <- as.matrix(seq(theta_lim[1L], theta_lim[2L], length.out=quadpts))
		fulldata <- object@Data$data
		tabdata <- object@Data$tabdatalong
		keep <- object@Data$Freq[[1L]] > 0L
		tabdata <- tabdata[keep, , drop=FALSE]
        USETABDATA <- TRUE
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)
        list_SEscores <- list_scores <- vector('list', MI)
        if(MI == 0) MI <- 1	    
        impute <- MI > 1
        opars <- pars
        estHess <- !full.scores | return.acov | full.scores.SE
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
            if(impute){
                while(TRUE){
                    pars <- try(imputePars(pars=opars, covB=covB, imputenums=imputenums, 
                                       constrain=object@constrain), silent=TRUE)
                    if(!is(pars, 'try-error')) break
                }
            }
            if(nfact < 3 || method == 'EAP'){
                ThetaShort <- Theta <- thetaComb(theta,nfact)
                if(length(prodlist) > 0L)
                    Theta <- prodterms(Theta,prodlist)
                W <- mirt_dmvnorm(ThetaShort,gp$gmeans,gp$gcov)
                itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, 
                                              CUSTOM.IND=CUSTOM.IND)
                log_itemtrace <- log(itemtrace)
                if(method == 'EAP' && return.acov){
                    tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=EAP, log_itemtrace=log_itemtrace,
                                   tabdata=tabdata, ThetaShort=ThetaShort, W=W, return.acov=TRUE, hessian=TRUE)
                } else {
            	    tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=EAP, log_itemtrace=log_itemtrace,
                                   tabdata=tabdata, ThetaShort=ThetaShort, W=W, hessian=estHess && method == 'EAP')
            	    scores <- tmp[ ,1:nfact, drop = FALSE]
            	    SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]                
                }
            }
    		if(method == "EAP"){
                #do nothing
    		} else if(method == "MAP"){
                tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=MAP, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist, 
                               CUSTOM.IND=CUSTOM.IND, return.acov=return.acov, hessian=estHess)
    		} else if(method == "ML"){
                isna <- apply(object@Data$tabdata[,-ncol(object@Data$tabdata),drop=FALSE], 1L, 
                              function(x) sum(is.na(x)))[keep]
    			allzero <- (rowSums(tabdata[,itemloc[-length(itemloc)], drop=FALSE]) + isna) == J
    			allmcat <- (rowSums(tabdata[,itemloc[-1L]-1L, drop=FALSE]) + isna) == J
    			scores[allmcat,] <- Inf
                SEscores[allmcat,] <- NA                
                scores[allzero,] <- -Inf
                SEscores[allzero,] <- NA
                tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=ML, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist,
                               CUSTOM.IND=CUSTOM.IND, return.acov=return.acov, hessian=estHess)
    		} else if(method == 'WLE'){
                if(nfact > 1L)
                    stop('WLE method only supported for unidimensional models')
                itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                              CUSTOM.IND=CUSTOM.IND)
                tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=WLE, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist, 
                               CUSTOM.IND=CUSTOM.IND, hessian=estHess)
            } else {
                stop('method not defined')
            }
    		if(return.acov){
    		    scores <- tmp
                if(nrow(scores) < ncol(scores)) scores <- t(scores)
    		} else {
    		    scores <- tmp[ ,1:nfact, drop = FALSE]
    		    SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
    		    colnames(scores) <- paste('F', 1:ncol(scores), sep='')
    		    if(impute){
    		        list_SEscores[[mi]] <- SEscores
    		        list_scores[[mi]] <- scores
    		    }
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
                tabdata2 <- object@Data$tabdata[keep, , drop=FALSE]
                stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
                sfulldata <- apply(object@Data$data, 1, paste, sep='', collapse = '/')
                scoremat <- scores[match(sfulldata, stabdata2), , drop = FALSE]
                if(return.acov){
                    ret <- vector('list', nrow(scoremat))
                    for(i in 1L:nrow(scoremat))
                        ret[[i]] <- matrix(scoremat[i,], nfact, nfact)
                    names(ret) <- 1:nrow(scoremat)
                    return(ret)
                }
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
            if(return.acov){
                ret <- vector('list', nrow(scores))
                for(i in 1L:nrow(scores))
                    ret[[i]] <- matrix(scores[i,], nfact, nfact)
                names(ret) <- paste0('pattern_', 1:nrow(scores))
                return(ret)
            }
            r <- object@Data$Freq[[1L]]
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
            ret <- cbind(object@Data$tabdata[keep, ,drop=FALSE],scores,SEscores)
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
	                      full.scores.SE, return.acov = FALSE)
	{
        class(object) <- 'ExploratoryClass'
        ret <- fscores(object, rotate = 'CONFIRMATORY', full.scores=full.scores, method=method, quadpts=quadpts,
                       response.pattern=response.pattern, returnER=returnER, verbose=verbose,
                       mean=gmean, cov=gcov, scores.only=scores.only, theta_lim=theta_lim, MI=MI,
                       full.scores.SE=full.scores.SE, return.acov = return.acov)
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
                          full.scores.SE, return.acov = FALSE)
    {
        pars <- object@pars
        ngroups <- length(pars)
        for(g in 1L:ngroups)
            class(pars[[g]]) <- 'ConfirmatoryClass'
        if(MI > 0){
            object <- assignInformationMG(object)
            pars <- object@pars
        }
        ret <- vector('list', length(pars))
        for(g in 1L:ngroups){
            tmp <- pars[[g]]
            tmp@Data <- object@Data
            tmp@Data$data <- tmp@Data$data[tmp@Data$group == tmp@Data$groupName[g],,drop=FALSE]
            tmp@Data$Freq[[1L]] <- tmp@Data$Freq[[g]]
            ret[[g]] <- fscores(tmp, rotate = 'CONFIRMATORY', full.scores=full.scores, method=method,
                           quadpts=quadpts, returnER=returnER, verbose=verbose, theta_lim=theta_lim,
                                mean=gmean[[g]], cov=gcov[[g]], scores.only=FALSE, MI=MI,
                           full.scores.SE=full.scores.SE, return.acov=return.acov)
        }
        names(ret) <- object@Data$groupNames
        if(full.scores){
            if(return.acov){
                group <- object@Data$group
                groupNames <- object@Data$groupNames
                count <- numeric(length(groupNames))
                out <- vector('list', length(group))
                for(i in 1L:length(group)){
                    which <- which(groupNames %in% group[i])
                    count[which] <- count[which] + 1L
                    out[[i]] <- ret[[which]][[count[which]]]
                }
                names(out) <- 1L:length(out)
                return(out)
            }
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
                ret <- ret[ ,!(colnames(ret) %in% colnames(object@Data$data)), drop=FALSE]
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

EAPsum <- function(x, full.scores = FALSE, quadpts = NULL, S_X2 = FALSE, gp, verbose, CUSTOM.IND,
                   theta_lim){
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
    theta <- seq(theta_lim[1L],theta_lim[2L],length.out = quadpts)
    Theta <- thetaComb(theta, x@nfact)
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
    thetas <- SEthetas <- matrix(0, nrow(L1), x@nfact)
    for(i in 1L:nrow(thetas)){
        thetas[i,] <- colSums(Theta * L1[i, ] * prior / sum(L1[i,] * prior))
        SEthetas[i,] <- sqrt(colSums((t(t(Theta) - thetas[i,]))^2 * L1[i, ] * prior / 
                                         sum(L1[i,] * prior)))
    }
    ret <- data.frame(Sum.Scores=Sum.Scores, Theta=thetas, SE.Theta=SEthetas)
    rownames(ret) <- ret$Sum.Scores
    if(full.scores){
        if(any(is.na(x@Data$data))) stop('Full scores requires a complete dataset (no N\'s)')
        dat <- x@Data$data
        adj <- x@Data$min
        if(any(adj > 0L)) message('Data adjusted so that every item has a lowest score of 0')
        dat <- t(t(dat) - adj)
        scores <- rowSums(dat)
        EAPscores <- ret[match(scores, ret$Sum.Scores), -1L, drop=FALSE]
        ret <- EAPscores[,1L:x@nfact, drop=FALSE]
    } else {
        dat <- x@Data$data
        E <- L1 %*% prior * nrow(dat)
        adj <- x@Data$min
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
        tmp <- suppressWarnings(expand.table(cbind(ret[,2L:(ncol(ret)-1L)], ret$observed)))
        pick <- 1L:x@nfact
        rxx <- apply(tmp[,pick, drop=FALSE], 2L, var) /
            (apply(tmp[,pick, drop=FALSE], 2L, var) + apply(tmp[,pick, drop=FALSE], 2L, 
                                                            function(x) mean(x^2)))
        attr(ret, 'fit') <- data.frame(df=df, X2=X2, p.X2 = pchisq(X2, df, lower.tail=FALSE),
                                       G2=G2, p.G2 = pchisq(G2, df, lower.tail=FALSE),
                                       rxx=as.list(rxx))
        if(verbose){
            print(attr(ret, 'fit'))
            cat('\n')
        }
    }
    ret
}
