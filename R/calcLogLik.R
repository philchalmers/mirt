# Monte Carlo Log-Likelihood Calculation
#
# Calculates a new object that contain the Monte Carlo estimated observed
# log-likelihood values for mirt objects estimated with the MH-RM algorithm.
# Can be estimated in parallel automatically by defining a parallel object with
# \code{\link{mirtCluster}}.
#
# @name calcLogLik
# @usage
# calcLogLik(object, ...)
#
# \S4method{calcLogLik}{ExploratoryClass}(object,
#    draws = 5000, G2 = TRUE)
# \S4method{calcLogLik}{ConfirmatoryClass}(object,
#    draws = 5000, G2 = TRUE)
# \S4method{calcLogLik}{MixedClass}(object,
#    draws = 5000)
# @aliases calcLogLik-method calcLogLik,ExploratoryClass-method
#   calcLogLik,ConfirmatoryClass-method calcLogLik,MixedClass-method
# @param object a model of class \code{ConfirmatoryClass} or \code{ExploratoryClass}
# @param draws the number of Monte Carlo draws
# @param G2 logical; estimate the G2 model fit statistic?
# @param ... parameters that are passed
# @section Methods:
# \describe{ \item{calcLogLik}{\code{signature(object = "ConfirmatoryClass")},
#           \code{signature(object = "ExploratoryClass")}, \code{signature(object = "MixedClass")} }
# }
# @return Returns an object with the log-likelihood and Monte Carlo standard errors,
#   and (possibly) the G^2 and other model fit statistic if there is no missing data.
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
# @docType methods
# @rdname calcLogLik-method
# @export calcLogLik
# @seealso \code{\link{mirt}}, \code{\link{multipleGroup}}, \code{\link{mixedmirt}}
# @keywords calcLogLik
# @examples
#
# \dontrun{
#
# # no parallel
# mod1withLogLik <- calcLogLik(mod1, draws=5000)
#
# #with parallel using detected number of cores
# library(parallel)
# mirtCluster(detectCores())
# mod1withLogLik <- calcLogLik(mod1, draws=5000)
#
#   }
#
setMethod(
	f = "calcLogLik",
	signature = signature(object = 'ExploratoryClass'),
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
                    random[[i]]@drawvals <- DrawValues(x=random[[i]], Theta=theta, pars=pars,
                                                       fulldata=fulldata, itemloc=itemloc,
                                                       offterm0=ot, CUSTOM.IND=CUSTOM.IND)
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
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, draws = 5000, G2 = TRUE, lrPars=NULL)
    {
        class(object) <- 'ExploratoryClass'
        ret <- calcLogLik(object, draws=draws, G2=G2, lrPars=lrPars)
        class(ret) <- 'ConfirmatoryClass'
        return(ret)
    }
)

setMethod(
    f = "calcLogLik",
    signature = signature(object = 'MixedClass'),
    definition = function(object, draws = 5000)
    {
        class(object) <- 'ExploratoryClass'
        ret <- calcLogLik(object, draws=draws, G2=FALSE, lrPars=object@lrPars)
        class(ret) <- 'MixedClass'
        return(ret)

    }
)
