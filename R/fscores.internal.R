setMethod(
	f = "fscores.internal",
	signature = 'SingleGroupClass',
	definition = function(object, rotate, Target, full.scores = FALSE, method = "EAP",
                          quadpts = NULL, response.pattern = NULL, append_response.pattern = TRUE,
	                      theta_lim, MI, pis=NULL, mixture=FALSE,
	                      returnER = FALSE, verbose = TRUE, gmean, gcov,
	                      plausible.draws, full.scores.SE, return.acov = FALSE,
                          QMC, custom_den = NULL, custom_theta = NULL,
	                      min_expected, converge_info, plausible.type, start,
	                      use_dentype_estimate, ...)
	{
        den_fun <- mirt_dmvnorm
        if(!is.null(custom_den)) den_fun <- custom_den
        if(use_dentype_estimate && !(method %in% c('EAP', 'EAPsum', 'plausible')))
            stop("use_dentype_estimate only supported for EAP, EAPsum, or plausible method", call.=FALSE)

	    #local functions for apply
	    MAP <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
	                    hessian, mirtCAT = FALSE, return.acov = FALSE, den_fun, max_theta, ...){
	        if(any(is.na(scores[ID, ])))
	            return(c(scores[ID, ], rep(NA, ncol(scores))))
            if(mirtCAT){
                estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars, patdata=tabdata[ID, ], den_fun=den_fun,
                                    itemloc=itemloc, gp=gp, prodlist=prodlist, max_theta=max_theta, hessian=hessian,
                                    CUSTOM.IND=CUSTOM.IND, ID=ID, iterlim=1, stepmax=1e-20, ...))
            } else {
    	        estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars, patdata=tabdata[ID, ], den_fun=den_fun,
    	                            itemloc=itemloc, gp=gp, prodlist=prodlist, max_theta=max_theta, hessian=hessian,
                                    CUSTOM.IND=CUSTOM.IND, ID=ID, ...))
            }
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2 + 1L))
            if(hessian){
                vcov <- try(solve(estimate$hessian))
                if(return.acov) return(vcov)
    	        SEest <- try(sqrt(diag(vcov)))
            } else SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest, estimate$code))
	    }
	    ML <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
	                   hessian, return.acov = FALSE, den_fun, max_theta, ...){
            if(any(scores[ID, ] %in% c(-Inf, Inf, NA)))
                return(c(scores[ID, ], rep(NA, ncol(scores) + 1L)))
            estimate <- try(nlm(MAP.mirt,scores[ID, ],pars=pars,patdata=tabdata[ID, ], den_fun=NULL,
    	                        itemloc=itemloc, gp=gp, prodlist=prodlist, ML=TRUE, max_theta=max_theta,
                                hessian=hessian, CUSTOM.IND=CUSTOM.IND, ID=ID, ...))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2 + 1L))
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
	        return(c(est, SEest, estimate$code))
	    }
	    WLE <- function(ID, scores, pars, tabdata, itemloc, gp, prodlist, CUSTOM.IND,
	                    hessian, data, DERIV, return.acov = FALSE, max_theta, ...){
	        if(any(is.na(scores[ID, ])))
	            return(c(scores[ID, ], rep(NA, ncol(scores))))
	        estimate <- try(nlm(WLE.mirt, scores[ID, ], pars=pars, patdata=tabdata[ID, ],
	                            itemloc=itemloc, gp=gp, prodlist=prodlist, data=data[ID, ], max_theta=max_theta,
	                            hessian=hessian, CUSTOM.IND=CUSTOM.IND, ID=ID, DERIV=DERIV, ...))
	        if(is(estimate, 'try-error'))
	            return(rep(NA, ncol(scores)*2 + 1L))
	        if(hessian){
	            vcov <- try(solve(estimate$hessian))
	            if(return.acov) return(vcov)
	            SEest <- try(sqrt(diag(vcov)))
	        } else SEest <- rep(NA, ncol(scores))
	        return(c(estimate$estimate, SEest, estimate$code))
	    }
	    EAP <- function(ID, log_itemtrace, tabdata, ThetaShort, W, hessian, scores, return.acov = FALSE,
	                    return_zeros = FALSE){
	        if(any(is.na(scores[ID, ])))
	            return(c(scores[ID, ], rep(NA, ncol(scores))))
            nfact <- ncol(ThetaShort)
	        L <- rowSums(log_itemtrace[ ,as.logical(tabdata[ID,]), drop = FALSE])
            expLW <- if(is.matrix(W)) exp(L) * W[ID, ] else exp(L) * W
            LW <- if(is.matrix(W)) L + log(W[ID, ]) else L + log(W)
            maxLW <- max(LW)
            nc <- sum(exp(LW - maxLW)) * exp(maxLW)
            if(nc == 0){
                if(return_zeros){
                    if(return.acov) return(matrix(NA, nfact, nfact))
                    return(numeric(nfact*2L + 1L))
                }
                warning(paste0('Unable to compute normalization constant for EAP estimates; ',
                               'consider using MAP estimates instead. Returning NaNs'),
                        call.=FALSE)
                return(c(rep(NaN, nfact*2), 0))
            }
	        thetas <- colSums(ThetaShort * expLW / nc)
            if(hessian){
    	        thetadif <- t((t(ThetaShort) - thetas))
                Thetaprod <- matrix(0, nrow(ThetaShort), nfact * (nfact + 1L)/2L)
                ind <- 1L
                for(i in seq_len(nfact)){
                    for(j in seq_len(nfact)){
                        if(i <= j){
                            Thetaprod[,ind] <- thetadif[,i] * thetadif[,j]
                            ind <- ind + 1L
                        }
                    }
                }
                vcov <- matrix(0, nfact, nfact)
                vcov[lower.tri(vcov, TRUE)] <- colSums(Thetaprod * expLW / nc)
                if(nfact > 1L) vcov <- vcov + t(vcov) - diag(diag(vcov))
                if(return.acov) return(vcov)
    	        SE <- sqrt(diag(vcov))
            } else SE <- rep(NA, nfact)
	        return(c(thetas, SE, 1))
	    }

        if(plausible.draws > 0){
            if(plausible.type == 'MH'){
                dots <- list(...)
                technical <- if(!is.null(dots$technical)) dots$technical else list()
                technical$plausible.draws <- plausible.draws
                formulas <- try(extract.mirt(object, 'lrformulas'), TRUE)
                if(!is(formulas, 'try-error'))
                    stop('MH plausible.type currently not supported for latent regression model', call.=FALSE)
                sv <- mod2values(object)
                sv$est <- FALSE
                ret <- mirt(extract.mirt(object, 'data'),
                            extract.mirt(object, 'model'),
                            extract.mirt(object, 'itemtype'),
                            pars=sv,
                            method='MHRM',
                            technical=technical)
                if(plausible.draws == 1L) return(ret[[1L]])
                else return(ret)
            } else if(plausible.type == 'normal'){
                fs <- fscores(object, rotate=rotate, Target=Target, full.scores = TRUE, method=method,
                              quadpts = quadpts, theta_lim=theta_lim, verbose=FALSE, cov=gcov,
                              return.acov = FALSE, QMC=QMC, custom_den=custom_den, converge_info=FALSE, ...)
                if(any(is.na(fs)))
                    stop('Plausible values cannot be drawn for completely empty response patterns.
                         Please remove these from your analysis.', call.=FALSE)
                fs_acov <- fscores(object, rotate = rotate, Target=Target, full.scores = TRUE, method=method,
                              quadpts = quadpts, theta_lim=theta_lim, verbose=FALSE,
                              plausible.draws=0, full.scores.SE=FALSE, cov=gcov,
                              return.acov = TRUE, QMC=QMC, custom_den=custom_den, converge_info=FALSE, ...)
                suppressWarnings(jit <- myLapply(seq_len(nrow(fs)), function(i, mu, sig)
                    mirt_rmvnorm(plausible.draws, mean = mu[i,], sigma = sig[[i]]),
                    mu=fs, sig=fs_acov))
                if(any(sapply(jit, is.nan)))
                    stop('Could not draw unique plausible values. Response pattern ACOVs may
                         not be positive definite')
                ret <- vector('list', plausible.draws)
                for(i in seq_len(plausible.draws)){
                    ret[[i]] <- matrix(NA, nrow(fs), ncol(fs))
                    for(j in seq_len(nrow(fs))) ret[[i]][j,] <- jit[[j]][i,]
                }
                if(plausible.draws == 1L) return(ret[[1L]])
                else return(ret)
            } else stop('plausible.type not supported', call.=FALSE)
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
            nfact <- object@Model$nfact
            mins <- extract.mirt(object, 'mins')
            if(!all(mins == 0L))
                response.pattern <- response.pattern - matrix(mins, nrow(response.pattern),
                                                          ncol(response.pattern), byrow=TRUE)
            colnames(response.pattern) <- colnames(object@Data$data)
            newmod <- object
            if(nrow(response.pattern) > 1L){
                large <- suppressWarnings(mirt(response.pattern, nfact, technical=list(customK=object@Data$K),
                              large=TRUE))
                newmod@Data <- list(data=response.pattern, tabdata=large$tabdata2,
                                   tabdatalong=large$tabdata, Freq=large$Freq,
                                   K=extract.mirt(object, 'K'), mins=rep(0L, ncol(response.pattern)))
                ret <- fscores(newmod, rotate=rotate, Target=Target, full.scores=TRUE,
                               method=method, quadpts=quadpts, verbose=FALSE, full.scores.SE=TRUE,
                               response.pattern=NULL, return.acov=return.acov, theta_lim=theta_lim,
                               MI=MI, mean=gmean, cov=gcov, custom_den=custom_den, QMC=QMC,
                               custom_theta=custom_theta, converge_info=converge_info,
                               start=start, use_dentype_estimate=use_dentype_estimate, ...)
                if(return.acov) return(ret)
                if(append_response.pattern) ret <- cbind(response.pattern, ret)
            } else {
                pick <- which(!is.na(response.pattern))
                rp <- response.pattern[,pick,drop=FALSE]
                large <- suppressWarnings(mirt(rp, nfact, large=TRUE,
                                            technical=list(customK=object@Data$K[pick])))
                newmod@Data <- list(data=rp, tabdata=large$tabdata2, K=object@Data$K[pick],
                                    tabdatalong=large$tabdata, Freq=large$Freq,
                                    mins=rep(0L, ncol(response.pattern))[pick])
                newmod@ParObjects$pars <- newmod@ParObjects$pars[c(pick, length(newmod@ParObjects$pars))]
                newmod@Model$itemloc <- c(1L, 1L + cumsum(object@Data$K[pick]))
                if(newmod@Options$exploratory)
                    stop('exploratory models not supported for single response pattern inputs', call.=FALSE)
                ret <- fscores(newmod, rotate=rotate, Target=Target, full.scores=TRUE,
                               method=method, quadpts=quadpts, verbose=FALSE, full.scores.SE=TRUE,
                               response.pattern=NULL, return.acov=return.acov, theta_lim=theta_lim,
                               MI=MI, mean=gmean, cov=gcov, custom_den=custom_den, QMC=QMC,
                               custom_theta=custom_theta, converge_info=converge_info,
                               start=start, use_dentype_estimate=use_dentype_estimate, ...)
                if(return.acov) return(ret)
                if(append_response.pattern) ret <- cbind(response.pattern, ret)
            }
            if(return.acov) return(ret)
            if(!all(mins == 0L) && append_response.pattern)
                ret[,seq_len(ncol(response.pattern))] <- ret[,seq_len(ncol(response.pattern))] +
                    matrix(mins, nrow(ret), ncol(response.pattern), byrow=TRUE)
            return(ret)
        }
        dots <- list(...)
        discrete <- FALSE
        if(object@Model$nfact > 3L && !QMC && method %in% c('EAP', 'EAPsum'))
            warning('High-dimensional models factor scores should use quasi-Monte Carlo integration. Pass QMC=TRUE',
                    call.=FALSE)
        if(method == 'Discrete' || method == 'DiscreteSum'){
            discrete <- TRUE
            method <- ifelse(method == 'Discrete', 'EAP', 'EAPsum')
        }
        if(mixture) discrete <- TRUE
        mirtCAT <- FALSE
        if(!is.null(dots$mirtCAT)) mirtCAT <- TRUE
        pars <- object@ParObjects$pars
		K <- object@Data$K
        J <- length(K)
        CUSTOM.IND <- object@Internals$CUSTOM.IND
        prodlist <- attr(pars, 'prodlist')
        nfact <- object@Model$nfact
        itemloc <- object@Model$itemloc
        gp <- if(!discrete) ExtractGroupPars(object@ParObjects$pars[[length(itemloc)]]) else list()
        if(object@Options$exploratory){
            so <- summary(object, rotate=rotate, Target=Target, verbose = FALSE)
            a <- rotateLambdas(so)
            for(i in seq_len(J))
                pars[[i]]@par[seq_len(nfact)] <- a[i, ]
            gp$gmeans <- rep(0, nfact)
            gp$gcov <- so$fcor
        }
        if(!is.null(gmean)) gp$gmeans <- gmean
        if(!is.null(gcov)) gp$gcov <- gcov
        if(method == 'EAPsum') return(EAPsum(object, full.scores=full.scores, full.scores.SE=full.scores.SE,
                                             quadpts=quadpts, gp=gp, verbose=verbose,
                                             CUSTOM.IND=CUSTOM.IND, theta_lim=theta_lim,
                                             discrete=discrete, QMC=QMC, den_fun=den_fun,
                                             min_expected=min_expected, pis=pis, mixture=mixture,
                                             use_dentype_estimate=use_dentype_estimate, ...))
		theta <- as.matrix(seq(theta_lim[1L], theta_lim[2L], length.out=quadpts))
		LR <- .hasSlot(object@Model$lrPars, 'beta')
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
            if(length(object@vcov) == 1L)
                stop('Stop an information matrix must be computed for imputations', call.=FALSE)
            if(!object@OptimInfo$secondordertest){
                stop('Information matrix is not positive definite', call.=FALSE)
            } else {
                names <- colnames(object@vcov)
                imputenums <- as.numeric(sapply(names, function(x, split){
                    strsplit(x, split=split)[[1L]][2L]
                }, split='\\.'))
                covB <- object@vcov
            }
            pre.ev <- eigen(covB)
        }
        for(mi in seq_len(MI)){
            if(impute){
                while(TRUE){
                    pars <- try(imputePars(pars=opars, imputenums=imputenums,
                                       constrain=object@Model$constrain, pre.ev=pre.ev), silent=TRUE)
                    if(!is(pars, 'try-error')) break
                }
            }
            if(nfact < 3 || method == 'EAP' && !mirtCAT){
                if(discrete){
                    ThetaShort <- Theta <- object@Model$Theta
                    W <- if(mixture) do.call(c, object@Internals$Prior) else object@Internals$Prior[[1L]]
                } else {
                    if(is.null(custom_theta)){
                        ThetaShort <- Theta <- if(QMC){
                            QMC_quad(npts=quadpts, nfact=nfact, lim=theta_lim)
                        } else thetaComb(theta,nfact)
                    } else {
                        if(ncol(custom_theta) != object@Model$nfact)
                            stop('ncol(custom_theta) does not match model', call.=FALSE)
                        ThetaShort <- Theta <- custom_theta
                    }
                    if(length(prodlist) > 0L)
                        Theta <- prodterms(Theta,prodlist)
                    W <- if(QMC) rep(1, nrow(Theta)) else
                        den_fun(ThetaShort, mean=gp$gmeans, sigma=gp$gcov, quad=LR, ...)
                    W <- W/sum(W)
                }
                if(use_dentype_estimate){
                    Theta <- ThetaShort <- object@Model$Theta
                    W <- object@Internals$Prior[[1L]]
                    den_fun <- pars[[J+1L]]@den
                }
                itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                              CUSTOM.IND=CUSTOM.IND, pis=pis)
                log_itemtrace <- log(itemtrace)
                if(mixture) ThetaShort <- thetaStack(ThetaShort, length(pis))
                if(method == 'EAP' && return.acov){
                    tmp <- myApply(X=matrix(seq_len(nrow(scores))), MARGIN=1L, FUN=EAP, log_itemtrace=log_itemtrace,
                                   tabdata=tabdata, ThetaShort=ThetaShort, W=W, return.acov=TRUE,
                                   scores=scores, hessian=TRUE)
                } else {
            	    tmp <- myApply(X=matrix(seq_len(nrow(scores))), MARGIN=1L, FUN=EAP, log_itemtrace=log_itemtrace,
                                   tabdata=tabdata, ThetaShort=ThetaShort, W=W, scores=scores,
                                   hessian=estHess && method == 'EAP', return_zeros=method != 'EAP')
                }
                scores <- tmp[ ,seq_len(nfact), drop = FALSE]
                SEscores <- tmp[ , seq_len(nfact) + nfact, drop = FALSE]
            }
            if(!is.null(start) && method != "EAP"){ #replace scores with start
                if(all(dim(scores) == dim(start)))
                    scores <- start
            }
    		if(method == "EAP"){
                #do nothing
    		} else if(method == "MAP"){
                tmp <- myApply(X=matrix(seq_len(nrow(scores))), MARGIN=1L, FUN=MAP, scores=scores, pars=pars,
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
                tmp <- myApply(X=matrix(seq_len(nrow(scores))), MARGIN=1L, FUN=ML, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist, den_fun=NULL,
                               CUSTOM.IND=CUSTOM.IND, return.acov=return.acov, hessian=estHess,
                               ...)
    		} else if(method == 'WLE'){
    		    DERIV <- vector('list', extract.mirt(object, 'nitems'))
    		    cls <- sapply(object@ParObjects$pars, class)
    		    for(i in seq_len(length(cls)-1L))
    		        DERIV[[i]] <- selectMethod(DerivTheta, c(cls[i], 'matrix'))
                tmp <- myApply(X=matrix(seq_len(nrow(scores))), MARGIN=1L, FUN=WLE, scores=scores, pars=pars,
                               tabdata=tabdata, itemloc=itemloc, gp=gp, prodlist=prodlist, DERIV=DERIV,
                               CUSTOM.IND=CUSTOM.IND, hessian=estHess, data=object@Data$tabdata, ...)
            } else {
                stop('method not defined', call.=FALSE)
            }
    		if(return.acov){
    		    scores <- tmp
    		    converge_info_vec <- rep(1, nrow(scores))
                if(nrow(scores) < ncol(scores)) scores <- t(scores)
    		} else {
    		    scores <- tmp[ ,seq_len(nfact), drop = FALSE]
    		    SEscores <- tmp[ , seq_len(nfact) + nfact, drop = FALSE]
    		    colnames(scores) <- paste('F', seq_len(ncol(scores)), sep='')
    		    converge_info_vec <- tmp[,ncol(tmp)]
    		    if(impute){
    		        list_SEscores[[mi]] <- SEscores
    		        list_scores[[mi]] <- scores
    		    }
    		}
        }
        if(impute){
            tmp <- averageMI(list_scores, list_SEscores, as.data.frame=FALSE)
            scores <- tmp[[1L]]
            SEscores <- tmp[[2L]]
        }
        if(any(is.na(scores) & !is.nan(scores)))
            warning('NAs returned for response patterns with no data. Consider removing',
                    call.=FALSE)
		if (full.scores){
            if(USETABDATA){
                tabdata2 <- object@Data$tabdata[keep, , drop=FALSE]
                stabdata2 <- apply(tabdata2, 1, paste, sep='', collapse = '/')
                sfulldata <- apply(object@Data$data, 1, paste, sep='', collapse = '/')
                scoremat <- scores[match(sfulldata, stabdata2), , drop = FALSE]
                if(return.acov){
                    ret <- vector('list', nrow(scoremat))
                    for(i in seq_len(nrow(scoremat)))
                        ret[[i]] <- matrix(scoremat[i,], nfact, nfact)
                    names(ret) <- seq_len(nrow(scoremat))
                    return(ret)
                }
                SEscoremat <- SEscores[match(sfulldata, stabdata2), , drop = FALSE]
                converge_info_mat <- converge_info_vec[match(sfulldata, stabdata2)]
    			colnames(scoremat) <- colnames(scores)
    			colnames(SEscoremat) <- paste0('SE_',colnames(scores))
            } else {
                scoremat <- scores
                SEscoremat <- SEscores
                converge_info_mat <- converge_info_vec
                if(return.acov){
                    ret <- vector('list', nrow(scoremat))
                    for(i in seq_len(nrow(scoremat)))
                        ret[[i]] <- matrix(scoremat[i,], nfact, nfact)
                    names(ret) <- seq_len(nrow(scoremat))
                    return(ret)
                } else colnames(SEscoremat) <- paste0('SE_',colnames(scores))
            }
            if(full.scores.SE)
                scoremat <- cbind(scoremat, SEscoremat)
            if(converge_info) scoremat <- cbind(scoremat, converged=converge_info_mat)
            return(scoremat)
		} else {
            if(return.acov){
                ret <- vector('list', nrow(scores))
                for(i in seq_len(nrow(scores)))
                    ret[[i]] <- matrix(scores[i,], nfact, nfact)
                names(ret) <- paste0('pattern_', 1:nrow(scores))
                return(ret)
            }
            r <- object@Data$Freq[[1L]]
            T <- E <- matrix(NA, 1L, ncol(scores))
            for(i in seq_len(nrow(scores))){
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
			    if(object@Options$exploratory) cat("\nRotate: ", rotate)
                cat("\n\nEmpirical Reliability:\n\n")
                print(round(reliability, 4))
			}
			colnames(SEscores) <- paste('SE_', colnames(scores), sep='')
            ret <- cbind(object@Data$tabdata[keep, ,drop=FALSE],scores,SEscores)
            if(converge_info) ret <- cbind(ret, converged=converge_info_vec)
            if(nrow(ret) > 1L) ret <- ret[do.call(order, as.data.frame(ret[,seq_len(J)])), ]
			return(ret)
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
                nclass <- ncol(object@Model$Theta)
                ret <- lapply(ret, function(x, nclass){
                  nx <- x[,seq_len(ncol(x)-nclass)]
                  names <- colnames(x)
                  colnames(nx) <- c(names[seq_len(ncol(nx)-nclass)], paste0('Class_', seq_len(nclass)))
                  nx
                }, nclass=nclass)
            } else if(method == 'DiscreteSum'){
                names <- paste0('Class_', seq_len(object@Model$nfact))
                names2 <- paste0('SE.Theta.', seq_len(object@Model$nfact))
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
                          full.scores.SE, return.acov = FALSE, QMC, plausible.draws, ...)
    {
        pars <- object@ParObjects$pars
        ngroups <- length(pars)
        for(g in seq_len(ngroups))
            class(pars[[g]]) <- 'SingleGroupClass'
        ret <- vector('list', length(pars))
        for(g in seq_len(ngroups)){
            tmp_obj <- MGC2SC(object, g)
            ret[[g]] <- fscores(tmp_obj, rotate = rotate, full.scores=full.scores, method=method,
                           quadpts=quadpts, returnER=returnER, verbose=verbose, theta_lim=theta_lim,
                                mean=gmean[[g]], cov=gcov[[g]], MI=MI, plausible.draws=plausible.draws,
                           full.scores.SE=full.scores.SE, return.acov=return.acov, QMC=QMC, ...)
        }
        names(ret) <- object@Data$groupNames
        if(plausible.draws > 0){
            pv <- plausible.draws
            out <- matrix(NA, nrow(object@Data$data), ncol(ret[[1L]][[1L]]))
            out2 <- vector('list', pv)
            colnames(out) <- colnames(ret[[1L]][[1L]])
            for(i in seq_len(pv)){
                for(g in seq_len(object@Data$ngroups)){
                    wch <- which(object@Data$group == object@Data$groupNames[g])
                    for(j in seq_len(ncol(ret[[1L]][[1L]])))
                        out[wch, j] <- ret[[g]][[i]][,j]
                }
                out2[[i]] <- out
            }
            return(out2)
        }
        if(full.scores){
            if(return.acov){
                group <- object@Data$group
                groupNames <- object@Data$groupNames
                count <- numeric(length(groupNames))
                out <- vector('list', length(group))
                for(i in seq_len(length(group))){
                    which <- which(groupNames %in% group[i])
                    count[which] <- count[which] + 1L
                    out[[i]] <- ret[[which]][[count[which]]]
                }
                names(out) <- seq_len(length(out))
                return(out)
            }
            out <- matrix(NA, nrow(object@Data$data), ncol(ret[[1L]]))
            colnames(out) <- colnames(ret[[1L]])
            for(g in seq_len(object@Data$ngroups)){
                wch <- which(object@Data$group == object@Data$groupNames[g])
                for(j in seq_len(ncol(ret[[1L]])))
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

#------------------------------------------------------------------------------
setMethod(
    f = "fscores.internal",
    signature = 'MixtureClass',
    definition = function(object, gmeans, gcov, ...)
    {
        class(object) <- 'SingleGroupClass'
        pis <- extract.mirt(object, 'pis')
        fscores.internal(object, mixture = TRUE,
                gmean=NULL, gcov=NULL, pis=pis, ...)
    }
)

#------------------------------------------------------------------------------

# MAP scoring for mirt
MAP.mirt <- function(Theta, pars, patdata, itemloc, gp, prodlist, CUSTOM.IND, ID,
                     ML=FALSE, den_fun, max_theta, ...)
{
    if(any(abs(Theta) > max_theta)) return(1e10)
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

WLE.mirt <- function(Theta, pars, patdata, itemloc, gp, prodlist, CUSTOM.IND, ID, data, DERIV,
                     max_theta)
{
    if(any(abs(Theta) > max_theta)) return(1e10)
    Theta <- matrix(Theta, nrow=1L)
    ThetaShort <- Theta
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta,prodlist)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
    L <- sum(log(itemtrace)[as.logical(patdata)])
    if(ncol(ThetaShort) == 1L){
        infos <- numeric(length(data))
        for(i in seq_len(length(infos))){
            if(!is.na(data[i]))
                infos[i] <- ItemInfo2(x=pars[[i]], Theta=Theta, total.info=TRUE, DERIV=DERIV[[i]],
                                      P=itemtrace[,itemloc[i]:(itemloc[i+1L]-1L),drop=FALSE])
        }
        infos <- sum(infos)
    } else {
        infos <- matrix(0, ncol(Theta), ncol(Theta))
        for(i in seq_len(length(data))){
            if(!is.na(data[i]))
                infos <- infos + ItemInfo2(x=pars[[i]], Theta=Theta, total.info=TRUE,
                                           MD=TRUE, DERIV=DERIV[[i]],
                                           P=itemtrace[,itemloc[i]:(itemloc[i+1L]-1L),drop=FALSE])
        }
        infos <- det(infos)
        if(closeEnough(infos, -1e-20, 1e-20))
            stop('Information matrix has a determinate of 0', call.=FALSE)
    }
    return(-(log(sqrt(infos)) + L))
}

gradnorm.WLE <- function(Theta, pars, patdata, itemloc, gp, prodlist, CUSTOM.IND){
    Theta <- matrix(Theta, nrow=1L)
    if(length(prodlist) > 0L)
        Theta <- prodterms(Theta,prodlist)
    nfact <- ncol(Theta)
    itemtrace <- matrix(0, ncol=length(patdata), nrow=nrow(Theta))
    dP <- d2P <- vector('list', nfact)
    for(i in seq_len(nfact))
        dP[[i]] <- d2P[[i]] <- itemtrace
    I <- numeric(1)
    dW <- dL <- numeric(nfact)
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND)
    for (i in seq_len(length(itemloc)-1L)){
        deriv <- DerivTheta(x=pars[[i]], Theta=Theta)
        for(k in seq_len(nfact)){
            dPitem <- d2Pitem <- matrix(0, 1, length(deriv[[1L]]))
            for(j in seq_len(length(deriv[[1L]]))){
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
    for(i in seq_len(nfact))
        dL[i] <- sum(patdata * dP[[i]] / itemtrace)
    grad <- dL + dW*I
    MIN <- sum(grad^2)
    MIN
}

EAPsum <- function(x, full.scores = FALSE, full.scores.SE = FALSE,
                   quadpts = NULL, S_X2 = FALSE, gp, verbose, CUSTOM.IND,
                   theta_lim, discrete, mixture, QMC, den_fun, min_expected,
                   which.items = 2:length(x@ParObjects$pars)-1,
                   use_dentype_estimate = FALSE, pis, ...){
    calcL1 <- function(itemtrace, K, itemloc){
        J <- length(K)
        L0 <- L1 <- matrix(1, sum(K-1L) + 1L, ncol(itemtrace))
        L0[1L:K[1L], ] <- itemtrace[1:K[1], ]
        nstar <- K[1L] + K[2L] - 3L
        Sum.Scores <- 1L:nrow(L0)-1L
        MAX.Scores <- cumsum(K-1L)
        for(i in seq_len(J-1L)){
            T <- itemtrace[itemloc[i+1L]:(itemloc[i+2L] - 1L), , drop=FALSE]
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
            for(i in seq_len(length(E)-1L)){
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

    prodlist <- attr(x@ParObjects$pars, 'prodlist')
    if(discrete){
        Theta <- ThetaShort <- x@Model$Theta
        prior <- if(mixture) do.call(c, x@Internals$Prior) else x@Internals$Prior[[1L]]
    } else {
        nfact <- x@Model$nfact
        ThetaShort <- Theta <- if(QMC){
            QMC_quad(npts=quadpts, nfact=nfact, lim=theta_lim)
        } else {
            theta <- seq(theta_lim[1L],theta_lim[2L],length.out = quadpts)
            thetaComb(theta,nfact)
        }
        prior <- if(QMC) rep(1, nrow(Theta)) else
            den_fun(Theta, mean=gp$gmeans, sigma=gp$gcov, ...)
        prior <- prior/sum(prior)
        if(length(prodlist) > 0L)
            Theta <- prodterms(Theta, prodlist)
    }
    if(use_dentype_estimate){
        Theta <- ThetaShort <- x@Model$Theta
        prior <- x@Internals$Prior[[1L]]
    }
    pars <- x@ParObjects$pars
    K <- x@Data$K
    J <- length(K)
    itemloc <- x@Model$itemloc
    itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc,
                                  CUSTOM.IND=CUSTOM.IND, pis=pis)
    itemtrace <- t(itemtrace)
    tmp <- calcL1(itemtrace=itemtrace, K=K, itemloc=itemloc)
    L1 <- tmp$L1
    Sum.Scores <- tmp$Sum.Scores
    if(S_X2){
        L1total <- L1 %*% prior
        Elist <- vector('list', J)
        for(i in which.items){
            KK <- K[-i]
            T <- itemtrace[c(itemloc[i]:(itemloc[i+1L]-1L)), , drop=FALSE]
            itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i+1L]-1L)), , drop=FALSE]
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
    if(mixture) ThetaShort <- thetaStack(ThetaShort, length(pis))
    thetas <- SEthetas <- matrix(0, nrow(L1), x@Model$nfact)
    for(i in seq_len(nrow(thetas))){
        expLW <- L1[i,] * prior
        LW <- log(L1[i,]) + log(prior)
        maxLW <- max(LW)
        nc <- sum(exp(LW - maxLW)) * exp(maxLW)
        if(nc == 0){
            warning('Unable to compute normalization constant for EAPsum estimates. Returning NaNs',
                     call.=FALSE)
            thetas[i, ] <- SEthetas[i, ] <- NaN
        } else {
            thetas[i, ] <- colSums(ThetaShort * expLW / nc)
            SEthetas[i, ] <- sqrt(colSums((t(t(ThetaShort) - thetas[i,]))^2 * expLW / nc))
        }
    }
    ret <- data.frame(Sum.Scores=Sum.Scores + sum(x@Data$min), Theta=thetas, SE.Theta=SEthetas)
    rownames(ret) <- ret$Sum.Scores
    if(full.scores){
        if(any(is.na(x@Data$data)))
            stop('Full scores requires a complete dataset (no NA\'s). If possible, pass na.rm=TRUE', call.=FALSE)
        dat <- x@Data$data
        adj <- extract.mirt(x, 'mins')
        dat <- t(t(dat) - adj)
        scores <- rowSums(dat)
        EAPscores <- ret[match(scores, Sum.Scores), -1L, drop=FALSE]
        pick <- if(full.scores.SE) seq_len(x@Model$nfact*2) else 1L:x@Model$nfact
        ret <- as.matrix(EAPscores[,pick, drop=FALSE])
        rownames(ret) <- NULL
    } else {
        dat <- x@Data$data
        if(any(is.na(dat)))
            stop('EAPsum scores are not meaningful when data contains missing values. If possible, pass na.rm=TRUE', call.=FALSE)
        E <- L1 %*% prior * nrow(dat)
        adj <- extract.mirt(x, 'mins')
        dat <- t(t(dat) - adj)
        Otmp <- matrix(table(sort(rowSums(dat))))
        got <- as.numeric(names(table(sort(rowSums(dat))))) + 1L
        O <- matrix(0, nrow(E), 1)
        O[got, 1] <- Otmp
        ret$observed <- O
        ret$expected <- E
        tmp <- collapseTotals(ret, min_expected)
        df <- tmp$df
        X2 <- tmp$X2
        tmp <- suppressWarnings(expand.table(cbind(ret[,2L:(ncol(ret)-1L)], ret$observed)))
        pick <- seq_len(x@Model$nfact)
        rxx <- apply(tmp[,pick, drop=FALSE], 2L, var) /
            (apply(tmp[,pick, drop=FALSE], 2L, var) + apply(tmp[,pick+x@Model$nfact, drop=FALSE], 2L,
                                                            function(x) mean(x^2)))
        names(rxx) <- paste0('rxx_Theta.', seq_len(x@Model$nfact))
        fit <- data.frame(df=df, X2=X2, p.X2 = suppressWarnings(pchisq(X2, df, lower.tail=FALSE)))
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
