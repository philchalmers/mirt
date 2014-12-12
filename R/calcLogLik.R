setMethod(
	f = "calcLogLik",
	signature = signature(object = 'SingleGroupClass'),
	definition = function(object, draws = 5000, G2 = TRUE, lrPars=NULL)
	{
        LLdraws <- function(LLDUMMY=NULL, nfact, N, grp, prodlist, fulldata, object, J, random, ot,
                            CUSTOM.IND, lrPars){
            theta <- mirt_rmvnorm(N,grp$gmeans, grp$gcov)
            if(length(lrPars)) theta <- theta + lrPars@mus
            if(length(prodlist) > 0L)
                theta <- prodterms(theta,prodlist)
            if(length(random) > 0L){
                for(i in 1L:length(random)){
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
        pars <- object@pars
        fulldata <- object@Data$fulldata
        prodlist <- object@prodlist
        itemloc <- object@itemloc
        N <- nrow(fulldata)
	    J <- length(pars)-1L
	    nfact <- length(ExtractLambdas(pars[[1L]])) - length(object@prodlist) - pars[[1L]]@nfixedeffects
        LL <- matrix(0, N, draws)
        grp <- ExtractGroupPars(pars[[length(pars)]])
        if(length(object@random) == 0L){
            ot <- matrix(0, 1L, J)
        } else ot <- OffTerm(object@random, J=J, N=N)
        LL <- t(myApply(X=LL, MARGIN=2L, FUN=LLdraws, nfact=nfact, lrPars=lrPars,
                        N=N, grp=grp, prodlist=prodlist, fulldata=fulldata, object=object, J=J,
                        random=object@random, ot=ot, CUSTOM.IND=object@CUSTOM.IND))
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
        expected <- .Call('sumExpected', t(data), tabdata, rwmeans, J, mirtClusterEnv$ncores)
		tabdata <- cbind(tabdata,expected*N)
        object@Pl <- expected
        nestpars <- nconstr <- 0L
        for(i in 1L:length(pars))
            nestpars <- nestpars + sum(pars[[i]]@est)
        if(length(object@constrain) > 0L)
            for(i in 1L:length(object@constrain))
                nconstr <- nconstr + length(object@constrain[[i]]) - 1L
        nfact <- object@nfact - length(prodlist)
        nmissingtabdata <- sum(is.na(rowSums(object@Data$tabdata)))
		if(G2){
			if(any(is.na(data))){
			    object@G2 <- NaN
			} else {
                pick <- r != 0L
                r <- r[pick]
                expected <- expected[pick]
				G2 <- 2 * sum(r*log(r/(sum(r)*expected)))
                df <- object@df
				object@G2 <- G2
                if(logLikpre == 0){
    				null.mod <- object@null.mod
    				object@TLI <- (null.mod@G2 / null.mod@df - G2/df) / (null.mod@G2 / null.mod@df - 1)
                }
			}
		}
		object@logLik <- logLik
        if(logLikpre < 0)
            object@logLik <- logLikpre
		object@SElogLik <- SElogLik
		return(object)
	}
)

setMethod(
    f = "calcLogLik",
    signature = signature(object = 'MixedClass'),
    definition = function(object, draws = 5000)
    {
        class(object) <- 'SingleGroupClass'
        ret <- calcLogLik(object, draws=draws, G2=FALSE, lrPars=object@lrPars)
        class(ret) <- 'MixedClass'
        return(ret)

    }
)
