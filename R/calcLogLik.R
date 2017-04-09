setMethod(
	f = "calcLogLik",
	signature = signature(object = 'SingleGroupClass'),
	definition = function(object, draws = 5000, G2 = TRUE, lrPars=NULL)
	{
        LLdraws <- function(LLDUMMY=NULL, nfact, N, grp, prodlist, fulldata, object, J, random, ot,
                            CUSTOM.IND, lrPars, pars, itemloc, lr.random){
            theta <- mirt_rmvnorm(N,grp$gmeans, grp$gcov)
            if(length(lr.random)){
                mus <- matrix(0, N, length(lr.random))
                for(i in seq_len(length(lr.random))){
                    mat <- matrix(0, lr.random[[i]]@ndim, lr.random[[i]]@ndim)
                    mat[lower.tri(mat, TRUE)] <- lr.random[[i]]@par
                    if(ncol(mat) > 1L)
                        mat <- mat + t(mat) - diag(diag(mat))
                    mus[,i] <- mirt_rmvnorm(nrow(lr.random[[i]]@drawvals), sigma = mat)[lr.random[[i]]@mtch]
                }
                mus <- rowSums(mus) + lrPars@mus
                theta <- theta + mus
            } else {
                if(length(lrPars)) theta <- theta + lrPars@mus
            }
            if(length(prodlist) > 0L)
                theta <- prodterms(theta,prodlist)
            if(length(random) > 0L){
                for(i in seq_len(length(random))){
                    mat <- matrix(0, random[[i]]@ndim, random[[i]]@ndim)
                    mat[lower.tri(mat, TRUE)] <- random[[i]]@par
                    if(ncol(mat) > 1L)
                        mat <- mat + t(mat) - diag(diag(mat))
                    random[[i]]@drawvals <- mirt_rmvnorm(nrow(random[[i]]@drawvals), sigma = mat)
                }
                ot <- OffTerm(random, J=J, N=N)
            }
            itemtrace <- computeItemtrace(pars=pars, Theta=theta, itemloc=itemloc, offterm=ot,
                                          CUSTOM.IND=CUSTOM.IND)
            return(rowSums(log(itemtrace)*fulldata))
        }
        pars <- object@ParObjects$pars
        fulldata <- object@Data$fulldata
        prodlist <- object@Model$prodlist
        itemloc <- object@Model$itemloc
        N <- nrow(fulldata)
	    J <- length(pars)-1L
	    nfact <- length(ExtractLambdas(pars[[1L]])) - length(object@Model$prodlist) - pars[[1L]]@nfixedeffects
        LL <- matrix(0, N, draws)
        grp <- ExtractGroupPars(pars[[length(pars)]])
        if(length(object@ParObjects$random) == 0L){
            ot <- matrix(0, 1L, J)
        } else ot <- OffTerm(object@ParObjects$random, J=J, N=N)
        LL <- t(myApply(X=LL, MARGIN=2L, FUN=LLdraws, nfact=nfact, lrPars=lrPars, pars=pars, itemloc=itemloc,
                        N=N, grp=grp, prodlist=prodlist, fulldata=fulldata, object=object, J=J,
                        random=object@ParObjects$random, ot=ot, CUSTOM.IND=object@Internals$CUSTOM.IND,
                        lr.random=object@ParObjects$lr.random))
        LL <- exp(LL)
        LL[is.nan(LL)] <- 0
        rwmeans <- rowMeans(LL)
        logLik <- sum(log(rwmeans))
        SElogLik <- sqrt(var(log(rowMeans(LL))) / draws)
        logLikpre <- 0
        if(G2 == 'return'){
            logLikpre <- logLik
            G2 <- TRUE
        }
		data <- object@Data$data
        tabdata <- cbind(object@Data$tabdata, object@Data$Freq[[1L]])
        r <- object@Data$Freq[[1L]]
        expected <- .Call('sumExpected', t(data), tabdata, rwmeans, J)
		tabdata <- cbind(tabdata,expected*N)
        object@Internals$Pl <- expected
        nestpars <- nconstr <- 0L
        for(i in seq_len(length(pars)))
            nestpars <- nestpars + sum(pars[[i]]@est)
        if(length(object@Model$constrain))
            for(i in seq_len(length(object@Model$constrain)))
                nconstr <- nconstr + length(object@Model$constrain[[i]]) - 1L
        nfact <- object@Model$nfact - length(prodlist)
        Fit <- object@Fit
		if(G2){
			if(any(is.na(data))){
			    Fit$G2 <- NaN
			} else {
                pick <- r != 0L
                r <- r[pick]
                expected <- expected[pick]
				G2 <- 2 * sum(r*log(r/(sum(r)*expected)))
				Fit$G2 <- G2
			}
		}
		Fit$logLik <- logLik
        if(logLikpre < 0)
            Fit$logLik <- logLikpre
		Fit$SElogLik <- SElogLik
		LP <- 0
		if(length(lrPars))
		    if(lrPars@any.prior)
		        LP <- LL.Priors(x=lrPars, LL=LP)
		if(length(object@ParObjects$random))
		    for(i in seq_len(length(object@ParObjects$random)))
		        if(object@ParObjects$random[[i]]@any.prior)
		            LP <- LL.Priors(x=object@ParObjects$random[[i]], LL=LP)
		for(i in seq_len(length(pars)))
		    if(pars[[i]]@any.prior)
		        LP <- LL.Priors(x=pars[[i]], LL=LP)
		Fit$logPrior <- unname(LP)
		object@Fit <- Fit
		return(object)
	}
)

setMethod(
    f = "calcLogLik",
    signature = signature(object = 'MixedClass'),
    definition = function(object, draws = 5000)
    {
        class(object) <- 'SingleGroupClass'
        ret <- calcLogLik(object, draws=draws, G2=FALSE, lrPars=object@Model$lrPars)
        class(ret) <- 'MixedClass'
        return(ret)

    }
)
