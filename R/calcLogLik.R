#' Monte Carlo Log-Likelihood Calculation
#'
#' Calculates a new object that contain the Monte Carlo estimated observed
#' log-likelihood values for mirt objects estimated with the MH-RM algorithm. 
#' Can be estimated in parallel automatically by defining a parallel object with
#' \code{\link{mirtCluster}}.
#'
#' @name calcLogLik
#' @usage
#' calcLogLik(object, ...)
#'
#' \S4method{calcLogLik}{ExploratoryClass}(object,
#'    draws = 5000, G2 = TRUE)
#' \S4method{calcLogLik}{ConfirmatoryClass}(object,
#'    draws = 5000, G2 = TRUE)
#' \S4method{calcLogLik}{MixedClass}(object,
#'    draws = 5000)
#' @aliases calcLogLik-method calcLogLik,ExploratoryClass-method
#' calcLogLik,ConfirmatoryClass-method calcLogLik,MixedClass-method
#' @param object a model of class \code{ConfirmatoryClass} or \code{ExploratoryClass}
#' @param draws the number of Monte Carlo draws
#' @param G2 logical; estimate the G2 model fit statistic?
#' @param ... parameters that are passed
#' @section Methods:
#' \describe{ \item{calcLogLik}{\code{signature(object = "ConfirmatoryClass")},
#' \code{signature(object = "ExploratoryClass")}, \code{signature(object = "MixedClass")} }
#' }
#' @return Returns an object with the log-likelihood and Monte Carlo standard errors, 
#' and (possibly) the G^2 and other model fit statistic if there is no missing data.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @docType methods
#' @rdname calcLogLik-methods
#' @export calcLogLik
#' @seealso
#' \code{\link{mirt}}, \code{\link{multipleGroup}}
#' @keywords calcLogLik
#' @examples
#'
#' \dontrun{
#'
#' # no parallel
#' mod1withLogLik <- calcLogLik(mod1, draws=5000)
#'
#' #with parallel using detected number of cores
#' library(parallel)
#' mirtCluster(detectCores())
#' mod1withLogLik <- calcLogLik(mod1, draws=5000)
#'
#'   }
#'
setMethod(
	f = "calcLogLik",
	signature = signature(object = 'ExploratoryClass'),
	definition = function(object, draws = 5000, G2 = TRUE)
	{
        LLdraws <- function(LLDUMMY=NULL, nfact, N, grp, prodlist, fulldata, object, J, random, ot,
                            PROBTRACE){
            theta <- mvtnorm::rmvnorm(N,grp$gmeans, grp$gcov)            
            if(length(prodlist) > 0L)
                theta <- prodterms(theta,prodlist)            
            if(length(random) > 0L){                                
                for(i in 1L:length(random)){
                    random[[i]]@drawvals <- DrawValues(x=random[[i]], Theta=theta, pars=pars, 
                                                fulldata=fulldata, itemloc=itemloc, offterm0=ot)
                }
                ot <- OffTerm(random, J=J, N=N)                
            }
            itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)            
            for (i in 1L:J) itemtrace[ ,itemloc[i]:(itemloc[i+1L] - 1L)] <-
                PROBTRACE[[i]](x=pars[[i]], Theta=theta, ot=ot[,i])
            return(exp(rowSums(log(itemtrace)*fulldata)))
        }           
        pars <- object@pars
	    tol <- .Machine$double.eps
        fulldata <- object@fulldata
        prodlist <- object@prodlist
        itemloc <- object@itemloc
        N <- nrow(fulldata)
	    J <- length(pars)-1L
	    nfact <- length(ExtractLambdas(pars[[1L]])) - length(object@prodlist) - pars[[1L]]@nfixedeffects
        LL <- matrix(0, N, draws)
        grp <- ExtractGroupPars(pars[[length(pars)]])
        if(length(object@random) == 0L){
            ot <- matrix(0, 1, J)
        } else ot <- OffTerm(object@random, J=J, N=N)
        PROBTRACE <- vector('list', J)
        for(i in 1L:J)
            PROBTRACE[[i]] <- selectMethod(ProbTrace, c(class(pars[[i]]), 'matrix'))        
        if(!is.null(globalenv()$MIRTCLUSTER)){
            LL <- parallel::parApply(cl=globalenv()$MIRTCLUSTER, LL, MARGIN=1, FUN=LLdraws, nfact=nfact, 
                                     N=N, grp=grp, prodlist=prodlist, fulldata=fulldata, object=object, J=J, 
                                     random=object@random, ot=ot, PROBTRACE=PROBTRACE)
        } else for(draw in 1L:draws)
            LL[ ,draw] <- LLdraws(nfact=nfact, N=N, grp=grp, prodlist=prodlist, PROBTRACE=PROBTRACE,
                                  fulldata=fulldata, object=object, J=J, random=object@random, ot=ot)
        LL[is.nan(LL)] <- 0
        rwmeans <- rowMeans(LL)
        logLik <- sum(log(rwmeans))
        SElogLik <- sqrt(var(log(rowMeans(LL))) / draws)
        logLikpre <- 0
        if(G2 == 'return'){
            logLikpre <- logLik
            G2 <- TRUE
        }
		data <- object@data
        tabdata <- object@tabdata
        r <- tabdata[,ncol(tabdata)]
		expected <- rep(0,nrow(tabdata))
		for (j in 1L:nrow(tabdata)){
			TFvec <- colSums(ifelse(t(data) == tabdata[j,1L:J],1,0)) == J
			TFvec[is.na(TFvec)] <- FALSE
			expected[j] <- mean(rwmeans[TFvec])
		}
		expected[is.nan(expected)] <- NA
		tabdata <- cbind(tabdata,expected*N)
        object@Pl <- expected
		logN <- 0
		logr <- rep(0,length(r))
		for (i in 1L:N) logN <- logN + log(i)
		for (i in 1L:length(r))
			for (j in 1L:r[i])
				logr[i] <- logr[i] + log(j)
		if(sum(logr) != 0)
			logLik <- logLik + logN/sum(logr)
        nestpars <- nconstr <- 0L
        for(i in 1:length(pars))
            nestpars <- nestpars + sum(pars[[i]]@est)
        if(length(object@constrain) > 0L)
            for(i in 1:length(object@constrain))
                nconstr <- nconstr + length(object@constrain[[i]]) - 1L
        nfact <- object@nfact - length(prodlist)
        nmissingtabdata <- sum(is.na(rowSums(object@tabdata)))
        df <- length(r) - nestpars + nconstr + nfact*(nfact - 1)/2 - 1 - nmissingtabdata
		AIC <- (-2) * logLik + 2 * (length(r) - df - 1)
		BIC <- (-2) * logLik + (length(r) - df - 1)*log(N)
		if(G2){
			if(any(is.na(data))){
			    object@G2 <- object@X2 <- NaN
			} else {
                r <- r[!is.na(expected)]
                expected <- expected[!is.na(expected)]
				G2 <- 2 * sum(r*log(r/(sum(r)*expected)))
                X2 <- sum((r - sum(r)*expected)^2 / (sum(r)*expected))
				object@G2 <- G2
                object@X2 <- X2
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
		object@AIC <- AIC
		object@BIC <- BIC
		object@df <- as.integer(df)
		return(object)
	}
)

setMethod(
    f = "calcLogLik",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, draws = 2000, G2 = TRUE)
    {
        class(object) <- 'ExploratoryClass'
        ret <- calcLogLik(object, draws=draws, G2=G2)
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
        ret <- calcLogLik(object, draws=draws, G2=FALSE)
        class(ret) <- 'MixedClass'
        return(ret)
        
    }
)