#' Extract various elements from estimated model objects
#'
#' A generic function to extract the internal objects from any estimated model. If possible, the
#' function will returns a vector/scalar containing the desired elements, otherwise it will return a list.
#'
#' Objects which can be extracted from mirt objects include:
#'
#' \describe{
#'   \item{logLik}{observed log-likelihood}
#'   \item{logPrior}{log term contributed by prior parameter distributions}
#'   \item{G2}{goodness of fit statistic}
#'   \item{df}{degrees of freedom}
#'   \item{p}{p-value for G2 statistic}
#'   \item{RMSEA}{root mean-square error of approximation based on G2}
#'   \item{CFI}{CFI fit statistic}
#'   \item{TLI}{TLI fit statistic}
#'   \item{AIC}{AIC}
#'   \item{AICc}{corrected AIC}
#'   \item{BIC}{BIC}
#'   \item{SABIC}{sample size adjusted BIC}
#'   \item{DIC}{DIC}
#'   \item{F}{unrotated standardized loadings matrix}
#'   \item{h2}{factor communality estimates}
#'   \item{parvec}{vector containing uniquely estimated paramters}
#'   \item{vcov}{parameter covariance matrix (associated with parvec)}
#'   \item{nest}{number of freely estimated parameters}
#'   \item{LLhistory}{EM log-likelihood history}
#'   \item{exp_resp}{expected probability of the unique response patterns}
#' }
#'
#' @aliases extract.mirt
#' @export extract.mirt
#' @param x mirt model of class 'SingleGroupClass', 'MultipleGroupClass', 'MixedClass' or
#'   'DiscreteGroupClass'
#' @param what a character vector indicating what to extract. Can contain more than one element
#'
#' @keywords extract
#' @seealso \code{\link{extract.group}}
#' @export extract.item
#' @examples
#'
#' \dontrun{
#' mod <- mirt(Science, 1)
#'
#' extract.mirt(mod, 'logLik')
#' extract.mirt(mod, c('G2', 'df', 'p'))
#' extract.mirt(mod, c('F', 'h2'))
#'
#' }
extract.mirt <- function(x, what){
    ret <- vector('list', length(what))
    for(i in 1:length(ret)){
        ret[[i]] <- switch(what[i],
                      G2 = x@Fit$G2,
                      logLik = x@Fit$logLik,
                      p = x@Fit$p,
                      TLI = x@Fit$TLI,
                      CFI = x@Fit$CFI,
                      RMSEA = x@Fit$RMSEA,
                      df = x@Fit$df,
                      AIC = x@Fit$AIC,
                      AICc = x@Fit$AICc,
                      BIC = x@Fit$BIC,
                      SABIC = x@Fit$SABIC,
                      DIC = x@Fit$DIC,
                      logPrior = x@Fit$logPrior,
                      F = x@Fit$F,
                      h2 = x@Fit$h2,
                      parvec = x@Internals$shortpars,
                      vcov = x@vcov,
                      nest = x@Model$nest,
                      LLhistory = x@Internals$collectLL,
                      exp_resp = x@Internals$Pl,
                      NULL
        )
    }
    names(ret) <- what
    if(length(ret) == 1L) return(ret[[1L]])
    if(!any(what %in% c('F', 'h2', 'vcov', 'parvec', 'exp_resp')))
        ret <- do.call(c, ret)
    ret
}
