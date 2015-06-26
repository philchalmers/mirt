setMethod(
	f = "fscores.internal",
	signature = 'SingleGroupClass',
	definition = function(object, rotate, Target, full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, theta_lim, MI,
	                      returnER = FALSE, verbose = TRUE, gmean, gcov,
	                      plausible.draws, full.scores.SE, return.acov = FALSE,
                          QMC, custom_den = NULL, custom_theta = NULL, digits=4,
	                      min_expected, ...)
	{
        den_fun <- mirt_dmvnorm
        if(!is.null(custom_den)) den_fun <- custom_den

	    #local functions for apply
	    MAP <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
	                    hessian, mirtCAT = FALSE, return.acov = FALSE, den_fun, ...){
	        if(any(is.na(scores[ID, ])))
	            return(c(scores[ID, ], rep(NA, ncol(scores))))
            if(mirtCAT){
                estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars, patdata=tabdata[ID, ], den_fun=den_fun,
                                    itemloc=itemloc, gp=gp, prodlist=prodlist, hessian=hessian,
                                    CUSTOM.IND=CUSTOM.IND, ID=ID, iterlim=1, stepmax=1e-20, ...))
            } else {
    	        estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars, patdata=tabdata[ID, ], den_fun=den_fun,
    	                            itemloc=itemloc, gp=gp, prodlist=prodlist, hessian=hessian,
                                    CUSTOM.IND=CUSTOM.IND, ID=ID, ...))
            }
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
            if(hessian){
                vcov <- try(solve(estimate$hessian))
                if(return.acov) return(vcov)
    	        SEest <- try(sqrt(diag(vcov)))
            } else SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    ML <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
	                   hessian, return.acov = FALSE, den_fun, ...){
            if(any(scores[ID, ] %in% c(-Inf, Inf, NA)))
                return(c(scores[ID, ], rep(NA, ncol(scores))))
            estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars,patdata=tabdata[ID, ], den_fun=NULL,
    	                        itemloc=itemloc, gp=gp, prodlist=prodlist, ML=TRUE,
                                hessian=hessian, CUSTOM.IND=CUSTOM.IND, ID=ID, ...))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
	        est <- estimate$estimate
            if(hessian){
                pick <- diag(estimate$hessian) > 0
                if(!all(pick)){
    	            vcov_small <- try(solve(estimate$hessian[pick,pick,drop=FALSE]))
                    vcov <- matrix(0, length(pick), length(pick))
                    vcov[pick,pick] <- vcov_small
                } else vcov <- try(solve(estimate$hessian))
    	        if(return.acov) return(vcov)
    	        SEest <- try(sqrt(diag(vcov)))
                if(any(SEest > 30)){
                    est[SEest > 30] <- Inf * sign(est[SEest > 30])
                    SEest[SEest > 30] <- NA
                }
            } else SEest <- rep(NA, ncol(scores))
	        return(c(est, SEest))
	    }
	    WLE <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
	                    hessian, data, return.acov = FALSE, ...){
	        if(any(is.na(scores[ID, ])))
	            return(c(scores[ID, ], rep(NA, ncol(scores))))
	        estimate <- try(nlm(WLE.mirt, scores[ID, ], pars=pars, patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, data=data[ID, ],
	                            hessian=hessian, CUSTOM.IND=CUSTOM.IND, ID=ID, ...))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2))
	        if(hessian){
	            vcov <- try(solve(estimate$hessian))
	            if(return.acov) return(vcov)
	            SEest <- try(sqrt(diag(vcov)))
	        } else SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest))
	    }
	    EAP <- function(ID, log_itemtrace, tabdata, ThetaShort, W, hessian, scores, return.acov = FALSE){
	        if(any(is.na(scores[ID, ])))
	            return(c(scores[ID, ], rep(NA, ncol(scores))))
            nfact <- ncol(ThetaShort)
	        L <- rowSums(log_itemtrace[ ,as.logical(tabdata[ID,]), drop = FALSE])
            expLW <- if(is.matrix(W)) exp(L) * W[ID, ] else exp(L) * W
            maxL <- max(expLW)
	        thetas <- colSums(ThetaShort * expLW / (sum(expLW/maxL)*maxL))
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
                vcov[lower.tri(vcov, TRUE)] <- colSums(Thetaprod * expLW / (sum(expLW/maxL)*maxL))
                if(nfact > 1L) vcov <- vcov + t(vcov) - diag(diag(vcov))
                if(return.acov) return(vcov)
    	        SE <- sqrt(diag(vcov))
            } else SE <- rep(NA, nfact)
	        return(c(thetas, SE))
	    }

        if(plausible.draws > 0){
            fs <- fscores(object, rotate=rotate, Target=Target, full.scores = TRUE, method=method,
                          quadpts = quadpts, theta_lim=theta_lim, verbose=FALSE,
                          return.acov = FALSE, QMC=QMC, custom_den = NULL, ...)
            fs_acov <- fscores(object, rotate = rotate, Target=Target, full.scores = TRUE, method=method,
                          quadpts = quadpts, theta_lim=theta_lim, verbose=FALSE,
                          plausible.draws=0, full.scores.SE=FALSE,
                          return.acov = TRUE, QMC=QMC, custom_den = NULL, ...)
            ret <- vector('list', plausible.draws)
            for(i in 1L:plausible.draws){
                suppressWarnings(jit <- lapply(fs_acov, function(x) mirt_rmvnorm(1L, sigma = x)))
                jit <- do.call(rbind, jit)
                if(any(is.nan(jit)))
                    stop('Could not draw unique plausible values. Response pattern ACOVs may
                         not be positive definite')
                ret[[i]] <- fs + jit
            }
            if(plausible.draws == 1L) return(ret[[1L]])
            else return(ret)
        }
        if(return.acov && MI != 0)
            stop('simultaneous impute and return.acov option not supported', call.=FALSE)
	    if(return.acov && returnER)
	        stop('simultaneous returnER and return.acov option not supported', call.=FALSE)
        if(!is.null(response.pattern)){
            if(is.data.frame(response.pattern))
                response.pattern <- as.matrix(response.pattern)
            if(!is.matrix(response.pattern))
                response.pattern <- matrix(response.pattern, 1L)
            nfact <- object@nfact
            mins <- object@Data$mins
            if(!all(mins == 0L))
                response.pattern <- response.pattern - matrix(mins, nrow(response.pattern),
                                                          ncol(response.pattern), byrow=TRUE)
            colnames(response.pattern) <- colnames(object@Data$data)
            newmod <- object
            if(nrow(response.pattern) > 1L){
                large <- suppressWarnings(mirt(response.pattern, nfact, technical=list(customK=object@K),
                              large=TRUE))
                newmod@Data <- list(data=response.pattern, tabdata=large$tabdata2,
                                   tabdatalong=large$tabdata, Freq=large$Freq)
                ret <- fscores(newmod, rotate=rotate, Target=Target, full.scores=TRUE,
                               method=method, quadpts=quadpts, verbose=FALSE, full.scores.SE=TRUE,
                               response.pattern=NULL, return.acov=return.acov, theta_lim=theta_lim,
                               MI=MI, mean=gmean, cov=gcov, custom_den=custom_den,
                               custom_theta=custom_theta, ...)
                if(return.acov) return(ret)
                ret <- cbind(response.pattern, ret)
            } else {
                pick <- which(!is.na(response.pattern))
                rp <- response.pattern[,pick,drop=FALSE]
                large <- suppressWarnings(mirt(rp, nfact, large=TRUE,
                                            technical=list(customK=object@K[pick])))
                newmod@Data <- list(data=rp, tabdata=large$tabdata2,
                                    tabdatalong=large$tabdata, Freq=large$Freq)
                newmod@pars <- newmod@pars[c(pick, length(newmod@pars))]
                newmod@itemloc <- c(1L, 1L + cumsum(object@K[pick]))
                newmod@K <- object@K[pick]
                newmod@F <- newmod@F[pick, , drop=FALSE]
                newmod@h2 <- newmod@h2[pick]
                ret <- fscores(newmod, rotate=rotate, Target=Target, full.scores=TRUE,
                               method=method, quadpts=quadpts, verbose=FALSE, full.scores.SE=TRUE,
                               response.pattern=NULL, return.acov=return.acov, theta_lim=theta_lim,
                               MI=MI, mean=gmean, cov=gcov, custom_den=custom_den,
                               custom_theta=custom_theta, ...)
                if(return.acov) return(ret)
                ret <- cbind(response.pattern, ret)
            }
            if(return.acov) return(ret)
            if(!all(mins == 0L))
                ret[,1L:ncol(response.pattern)] <- ret[,1L:ncol(response.pattern)] +
                    matrix(mins, nrow(ret), ncol(response.pattern), byrow=TRUE)
            return(ret)
        }
        dots <- list(...)
        discrete <- FALSE
        if(method == 'Discrete' || method == 'DiscreteSum'){
            discrete <- TRUE
            method <- ifelse(method == 'Discrete', 'EAP', 'EAPsum')
        }
        mirtCAT <- FALSE
        if(!is.null(dots$mirtCAT)) mirtCAT <- TRUE
        pars <- object@pars
		K <- object@K
        J <- length(K)
        CUSTOM.IND <- object@CUSTOM.IND
        prodlist <- attr(pars, 'prodlist')
        nfact <- object@nfact
        nLambdas <- object@nfact
        itemloc <- object@itemloc
        gp <- ExtractGroupPars(object@pars[[length(itemloc)]])
        if(object@exploratory){
            so <- summary(object, rotate=rotate, Target=Target, verbose = FALSE, digits = Inf)
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
                                             CUSTOM.IND=CUSTOM.IND, theta_lim=theta_lim,
                                             discrete=discrete, QMC=QMC, den_fun=den_fun,
                                             min_expected=min_expected, ...))
		theta <- as.matrix(seq(theta_lim[1L], theta_lim[2L], length.out=quadpts))
		fulldata <- object@Data$data
		LR <- .hasSlot(object@lrPars, 'beta')
		USETABDATA <- TRUE
		if(LR){
            if(!(method %in% c('EAP', 'MAP')))
                warning('Latent regression information only used in MAP and EAP estimates',
                        call.=FALSE)
            full.scores <- TRUE
            USETABDATA <- FALSE
		    gp$gmeans <- fixef(object)
		    tabdata <- object@Data$fulldata[[1L]]
            keep <- rep(TRUE, nrow(tabdata))
		} else {
    		tabdata <- object@Data$tabdatalong
    		keep <- object@Data$Freq[[1L]] > 0L
    		tabdata <- tabdata[keep, , drop=FALSE]
		}
		SEscores <- scores <- matrix(0, nrow(tabdata), nfact)
		drop <- rowSums(tabdata) == 0L
		scores[drop, ] <- SEscores[drop, ] <- NA
        list_SEscores <- list_scores <- vector('list', MI)
        if(MI == 0) MI <- 1
        impute <- MI > 1
        opars <- pars
        estHess <- !full.scores | return.acov | full.scores.SE
        if(impute){
            if(length(object@information) == 1L)
                stop('Stop an information matrix must be computed for imputations', call.=FALSE)
            if(is(try(chol(object@information), silent=TRUE), 'try-error')){
                stop('Information matrix is not positive definite', call.=FALSE)
            } else {
                names <- colnames(object@information)
                imputenums <- as.numeric(sapply(names, function(x, split){
                    strsplit(x, split=split)[[1L]][2L]
                }, split='\\.'))
                covB <- solve(object@information)
            }
            pre.ev <- eigen(covB)
        }
        for(mi in 1L:MI){
            if(impute){
                while(TRUE){
                    pars <- try(imputePars(pars=opars, imputenums=imputenums,
                                       constrain=object@constrain, pre.ev=pre.ev), silent=TRUE)
                    if(!is(pars, 'try-error')) break
                }
            }
            if(nfact < 3 || method == 'EAP' && !mirtCAT){
                if(discrete){
                    ThetaShort <- Theta <- object@Theta
                    W <- object@Prior[[1L]]
                } else {
                    if(is.null(custom_theta)){
                        ThetaShort <- Theta <- if(QMC){
                            QMC_quad(npts=quadpts, nfact=nfact, lim=theta_lim)
                        } else thetaComb(theta,nfact)
                    } else {
                        if(ncol(custom_theta) != object@nfact)
                            stop('ncol(custom_theta) does not match model', call.=FALSE)
                        ThetaShort <- Theta <- custom_theta
                    }
                    if(length(prodlist) > 0L)
                        Theta <- prodterms(Theta,prodlist)
                    W <- den_fun(ThetaShort, mean=gp$gmeans, sigma=gp$gcov, quad=LR, ...)
                    W <- W/sum(W)
                }
                itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                              CUSTOM.IND=CUSTOM.IND)
                log_itemtrace <- log(itemtrace)
                if(method == 'EAP' && return.acov){
                    tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=EAP, log_itemtrace=log_itemtrace,
                                   tabdata=tabdata, ThetaShort=ThetaShort, W=W, return.acov=TRUE,
                                   scores=scores, hessian=TRUE)
                } else {
            	    tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=EAP, log_itemtrace=log_itemtrace,
                                   tabdata=tabdata, ThetaShort=ThetaShort, W=W, scores=scores,
                                   hessian=estHess && method == 'EAP')
            	    scores <- tmp[ ,1L:nfact, drop = FALSE]
            	    SEscores <- tmp[ ,-c(1L:nfact), drop = FALSE]
                }
            }
    		if(method == "EAP"){
                #do nothing
    		} else if(method == "MAP"){
                tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=MAP, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist, den_fun=den_fun,
                               CUSTOM.IND=CUSTOM.IND, return.acov=return.acov, hessian=estHess,
                               ...)
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
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist, den_fun=NULL,
                               CUSTOM.IND=CUSTOM.IND, return.acov=return.acov, hessian=estHess,
                               ...)
    		} else if(method == 'WLE'){
                tmp <- myApply(X=matrix(1L:nrow(scores)), MARGIN=1L, FUN=WLE, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist,
                               CUSTOM.IND=CUSTOM.IND, hessian=estHess, data=object@Data$tabdata, ...)
            } else {
                stop('method not defined', call.=FALSE)
            }
    		if(return.acov){
    		    scores <- tmp
                if(nrow(scores) < ncol(scores)) scores <- t(scores)
    		} else {
    		    scores <- tmp[ ,1:nfact, drop = FALSE]
    		    SEscores <- tmp[ ,-c(1:nfact), drop = FALSE]
    		    colnames(scores) <- paste('F', 1L:ncol(scores), sep='')
    		    if(impute){
    		        list_SEscores[[mi]] <- SEscores
    		        list_scores[[mi]] <- scores
    		    }
    		}
        }
        if(impute){
            tmp <- averageMI(list_scores, list_SEscores, as.data.frame=FALSE,
                             digits = 200)
            scores <- tmp[[1L]]
            SEscores <- tmp[[2L]]
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
                if(return.acov){
                    ret <- vector('list', nrow(scoremat))
                    for(i in 1L:nrow(scoremat))
                        ret[[i]] <- matrix(scoremat[i,], nfact, nfact)
                    names(ret) <- 1L:nrow(scoremat)
                    return(ret)
                } else colnames(SEscoremat) <- paste0('SE_',colnames(scores))
            }
            if(full.scores.SE)
                scoremat <- cbind(scoremat, SEscoremat)
            return(scoremat)
		} else {
            if(return.acov){
                ret <- vector('list', nrow(scores))
                for(i in 1L:nrow(scores))
                    ret[[i]] <- matrix(scores[i,], nfact, nfact)
                names(ret) <- paste0('pattern_', 1:nrow(scores))
                return(ret)
            }
            r <- object@Data$Freq[[1L]]
            T <- E <- matrix(NA, 1L, ncol(scores))
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
			if(verbose && !discrete){
                cat("\nMethod: ", method)
			    if(object@exploratory) cat("\nRotate: ", rotate)
                cat("\n\nEmpirical/Marginal Reliability:\n\n")
                print(round(reliability, digits))
			}
			colnames(SEscores) <- paste('SE_', colnames(scores), sep='')
            ret <- cbind(object@Data$tabdata[keep, ,drop=FALSE],scores,SEscores)
            if(nrow(ret) > 1L) ret <- ret[do.call(order, as.data.frame(ret[,1L:J])), ]
			return(round(ret, digits))
		}
	}
)

#------------------------------------------------------------------------------
setMethod(
    f = "fscores.internal",
    signature = 'DiscreteClass',
    definition = function(object, rotate, full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, theta_lim, MI,
                          returnER = FALSE, verbose = TRUE, gmean, gcov,
                          full.scores.SE, return.acov = FALSE, ...)
    {
        class(object) <- 'MultipleGroupClass'
        if(!any(method %in% c('EAP', 'EAPsum')))
            stop('Only EAP and EAPsum methods are supported for DiscreteClass objects', call.=FALSE)
        method <- ifelse(method == 'EAP', 'Discrete', 'DiscreteSum')
        ret <- fscores(object, full.scores=full.scores, method=method, quadpts=quadpts,
                       response.pattern=response.pattern, returnER=FALSE, verbose=verbose,
                       mean=gmean, cov=gcov, theta_lim=theta_lim, MI=MI,
                       full.scores.SE=FALSE, return.acov = FALSE, rotate='none', ...)
        if(!full.scores){
            if(method == 'Discrete'){
                nclass <- ncol(object@Theta)
                ret <- lapply(ret, function(x, nclass){
                  nx <- x[,1L:(ncol(x)-nclass)]
                  names <- colnames(x)
                  colnames(nx) <- c(names[1:(ncol(nx)-nclass)], paste0('Class_', 1L:nclass))
                  nx
                }, nclass=nclass)
            } else if(method == 'DiscreteSum'){
                names <- paste0('Class_', 1L:object@nfact)
                names2 <- paste0('SE.Theta.', 1L:object@nfact)
                ret <- lapply(ret, function(x, names, names2){
                    nx <- x[,!(colnames(x) %in% names2)]
                    colnames(nx) <- c('Sum.Scores', names, 'observed', 'expected')
                    nx
                }, names=names, names2=names2)
            }
        }
        if(length(ret) == 1L) ret <- ret[[1L]]
        return(ret)
    }
)

#------------------------------------------------------------------------------
setMethod(
    f = "fscores.internal",
    signature = 'MultipleGroupClass',
    definition = function(object, rotate, full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, theta_lim, MI,
                          returnER = FALSE, verbose = TRUE, gmean, gcov,
                          full.scores.SE, return.acov = FALSE, QMC, ...)
    {
        pars <- object@pars
        ngroups <- length(pars)
        for(g in 1L:ngroups)
            class(pars[[g]]) <- 'SingleGroupClass'
        if(MI > 0){
            object <- assignInformationMG(object)
            pars <- object@pars
        }
        ret <- vector('list', length(pars))
        for(g in 1L:ngroups){
            tmp_obj <- MGC2SC(object, g)
            ret[[g]] <- fscores(tmp_obj, rotate = rotate, full.scores=full.scores, method=method,
                           quadpts=quadpts, returnER=returnER, verbose=verbose, theta_lim=theta_lim,
                                mean=gmean[[g]], cov=gcov[[g]], MI=MI,
                           full.scores.SE=full.scores.SE, return.acov=return.acov, QMC=QMC, ...)
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
            out <- matrix(NA, nrow(object@Data$data), ncol(ret[[1L]]))
            colnames(out) <- colnames(ret[[1L]])
            for(g in 1L:object@Data$ngroups){
                wch <- which(object@Data$group == object@Data$groupNames[g])
                for(j in 1L:ncol(ret[[1L]]))
                    out[wch, j] <- ret[[g]][,j]
            }
            ret <- out
            rownames(ret) <- rownames(object@Data$data)
        }
        if(is.data.frame(ret))
            ret <- as.matrix(ret)
        return(ret)
    }
)

# MAP scoring for mirt
MAP.mirt <- function(Theta, pars, patdata, itemloc, gp, prodlist, CUSTOM.IND, ID,
                     ML=FALSE, den_fun, ...)
{
    Theta <- matrix(Theta, nrow=1L)
    ThetaShort <- Theta
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta,prodlist)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
    L <- sum(log(itemtrace)[as.logical(patdata)])
    mu <- if(is.matrix(gp$gmeans)) gp$gmeans[ID, ] else gp$gmeans
    if(!ML){
        prior <- den_fun(ThetaShort, mean=mu, sigma=gp$gcov, ...)
        if(prior < 1e-200) prior <- 1e-200
    }
    L <- ifelse(ML, -L, (-1)*(L + log(prior)))
    L
}

WLE.mirt <- function(Theta, pars, patdata, itemloc, gp, prodlist, CUSTOM.IND, ID, data)
{
    Theta <- matrix(Theta, nrow=1L)
    ThetaShort <- Theta
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta,prodlist)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
    L <- sum(log(itemtrace)[as.logical(patdata)])
    if(ncol(ThetaShort) == 1L){
        infos <- numeric(length(data))
        for(i in 1L:length(infos)){
            if(!is.na(data[i]))
                infos[i] <- ItemInfo2(x=pars[[i]], Theta=Theta, total.info=TRUE)
        }
        infos <- sum(infos)
    } else {
        infos <- matrix(0, ncol(Theta), ncol(Theta))
        for(i in 1L:length(data)){
            if(!is.na(data[i]))
                infos <- infos + ItemInfo2(x=pars[[i]], Theta=Theta, total.info=TRUE, MD=TRUE)
        }
        infos <- det(infos)
    }
    return(-(log(sqrt(infos)) + L))
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
                   theta_lim, discrete, QMC, den_fun, min_expected, ...){
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
    collapseTotals <- function(table, min_expected){

        E <- table$expected
        O <- table$observed
        while(TRUE){
            small <- E < min_expected
            if(!any(small)) break
            for(i in 1L:(length(E)-1L)){
                if(small[i]){
                    E[i+1L] <- E[i+1L] + E[i]
                    O[i+1L] <- O[i+1L] + O[i]
                    O[i] <- E[i] <- NA
                }
            }
            if(small[i+1L]){
                E[i] <- E[i+1L] + E[i]
                O[i] <- O[i+1L] + O[i]
                O[i+1L] <- E[i+1L] <- NA
            }
            O <- na.omit(O)
            E <- na.omit(E)
        }
        df <- length(O) - 1L
        X2 <- sum((O - E)^2 / E)
        list(X2=X2, df=df)
    }

    prodlist <- attr(x@pars, 'prodlist')
    if(discrete){
        Theta <- ThetaShort <- x@Theta
        prior <- x@Prior[[1L]]
    } else {
        nfact <- x@nfact
        ThetaShort <- Theta <- if(QMC){
            QMC_quad(npts=quadpts, nfact=nfact, lim=theta_lim)
        } else {
            theta <- seq(theta_lim[1L],theta_lim[2L],length.out = quadpts)
            thetaComb(theta,nfact)
        }
        prior <- den_fun(Theta, mean=gp$gmeans, sigma=gp$gcov, ...)
        prior <- prior/sum(prior)
        if(length(prodlist) > 0L)
            Theta <- prodterms(Theta, prodlist)
    }
    pars <- x@pars
    K <- x@K
    J <- length(K)
    itemloc <- x@itemloc
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
    itemtrace <- t(itemtrace)
    tmp <- calcL1(itemtrace=itemtrace, K=K, itemloc=itemloc)
    L1 <- tmp$L1
    maxLs <- apply(L1, 1L, max)
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
        thetas[i,] <- colSums(ThetaShort * L1[i, ] * prior / sum(L1[i,] * prior))
        SEthetas[i,] <- sqrt(colSums((t(t(ThetaShort) - thetas[i,]))^2 * L1[i, ] * prior /
                                         sum(L1[i,] * prior)))
    }
    ret <- data.frame(Sum.Scores=Sum.Scores + sum(x@Data$min), Theta=thetas, SE.Theta=SEthetas)
    rownames(ret) <- ret$Sum.Scores
    if(full.scores){
        if(any(is.na(x@Data$data)))
            stop('Full scores requires a complete dataset (no N\'s)', call.=FALSE)
        dat <- x@Data$data
        adj <- x@Data$min
        dat <- t(t(dat) - adj)
        scores <- rowSums(dat)
        EAPscores <- ret[match(scores, Sum.Scores), -1L, drop=FALSE]
        ret <- EAPscores[,1L:x@nfact, drop=FALSE]
        rownames(ret) <- NULL
    } else {
        dat <- x@Data$data
        if(any(is.na(dat)))
            stop('EAPsum scores are not meaningful when data contains missing values')
        E <- L1 %*% prior * nrow(dat)
        adj <- x@Data$min
        dat <- t(t(dat) - adj)
        Otmp <- matrix(table(sort(rowSums(dat))))
        got <- as.numeric(names(table(sort(rowSums(dat))))) + 1L
        O <- matrix(0, nrow(E), 1)
        O[got, 1] <- Otmp
        keep <- O != 0
        ret$observed <- O
        ret$expected <- E
        tmp <- collapseTotals(ret, min_expected)
        df <- tmp$df
        X2 <- tmp$X2
        tmp <- suppressWarnings(expand.table(cbind(ret[,2L:(ncol(ret)-1L)], ret$observed)))
        pick <- 1L:x@nfact
        rxx <- apply(tmp[,pick, drop=FALSE], 2L, var) /
            (apply(tmp[,pick, drop=FALSE], 2L, var) + apply(tmp[,pick+x@nfact, drop=FALSE], 2L,
                                                            function(x) mean(x^2)))
        names(rxx) <- paste0('rxx_Theta.', 1L:x@nfact)
        fit <- data.frame(df=df, X2=X2, p.X2 = pchisq(X2, df, lower.tail=FALSE))
        fit <- cbind(fit, t(as.data.frame(rxx)))
        rownames(fit) <- 'stats'
        attr(ret, 'fit') <- fit
        if(verbose && !discrete){
            print(attr(ret, 'fit'))
            cat('\n')
        }
    }
    ret
}
