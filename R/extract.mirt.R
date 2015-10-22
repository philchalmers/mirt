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
#'   \item{LLhistory}{EM log-likelihood history}
#'   \item{tabdata}{a tabular version of the raw response data input. Frequencies are stored
#'     in \code{freq}}
#'   \item{freq}{frequencies associated with \code{tabdata}}
#'   \item{K}{an integer vector indicating the number of unique elements for each item}
#'   \item{tabdatalong}{similar to \code{tabdata}, however the responses have been transformed into
#'     dummy coded variables}
#'   \item{fulldatalong}{analogous to \code{tabdatafull}, but for the raw input data instead of the
#'     tabulated frequencies}
#'   \item{exp_resp}{expected probability of the unique response patterns}
#'   \item{converged}{a logical value indicating whether the model terminated within
#'     the convergence criteria}
#'   \item{iterations}{number of iterations it took to reach the convergence criteria}
#'   \item{nest}{number of freely estimated parameters}
#'   \item{parvec}{vector containing uniquely estimated parameters}
#'   \item{vcov}{parameter covariance matrix (associated with parvec)}
#'   \item{condnum}{the condition number of the Hessian (if computed). Otherwise NA}
#'   \item{Prior}{prior density distribution for the latent traits}
#'   \item{secondordertest}{a logical indicating whether the model passed the second-order test
#'     based on the Hessian matrix. Indicates whether model is a potential local maximum solution}
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
                      K = x@Data$K,
                      parvec = x@Internals$shortpars,
                      vcov = x@vcov,
                      nest = x@Model$nest,
                      iterations = x@OptimInfo$iter,
                      LLhistory = x@Internals$collectLL,
                      exp_resp = x@Internals$Pl,
                      converged = x@OptimInfo$converged,
                      condnum = x@OptimInfo$condnum,
                      Prior = x@Internals$Prior,
                      secondordertest = x@OptimInfo$secondordertest,
                      tabdata = x@Data$tabdata,
                      freq = x@Data$Freq,
                      tabdatalong = x@Data$tabdatalong,
                      fulldatalong = x@Data$fulldata,
                      NULL
        )
    }
    names(ret) <- what
    if(length(ret) == 1L) return(ret[[1L]])
    if(!any(what %in% c('F', 'h2', 'vcov', 'parvec', 'exp_resp', 'Prior', 'freq',
                        'tabdata', 'tabdatalong', 'fulldatalong')))
        ret <- do.call(c, ret)
    ret
}
