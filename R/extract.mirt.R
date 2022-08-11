#' Extract various elements from estimated model objects
#'
#' A generic function to extract the internal objects from estimated models.
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
#'   \item{BIC}{BIC}
#'   \item{SABIC}{sample size adjusted BIC}
#'   \item{HQ}{HQ}
#'   \item{F}{unrotated standardized loadings matrix}
#'   \item{h2}{factor communality estimates}
#'   \item{LLhistory}{EM log-likelihood history}
#'   \item{tabdata}{a tabular version of the raw response data input. Frequencies are stored
#'     in \code{freq}}
#'   \item{freq}{frequencies associated with \code{tabdata}}
#'   \item{K}{an integer vector indicating the number of unique elements for each item}
#'   \item{mins}{an integer vector indicating the lowest category found in the input \code{data}}
#'   \item{model}{input model syntax}
#'   \item{method}{estimation method used}
#'   \item{itemtype}{a vector of item types for each respective item (e.g., 'graded', '2PL', etc)}
#'   \item{itemnames}{a vector of item names from the input data}
#'   \item{factorNames}{a vector of factor names from the model definition}
#'   \item{rowID}{an integer vector indicating all valid row numbers used in the model estimation
#'    (when all cases are used this will be \code{1:nrow(data)})}
#'   \item{data}{raw input data of item responses}
#'   \item{covdata}{raw input data of data used as covariates}
#'   \item{tabdatalong}{similar to \code{tabdata}, however the responses have been transformed into
#'     dummy coded variables}
#'   \item{fulldatalong}{analogous to \code{tabdatafull}, but for the raw input data instead of the
#'     tabulated frequencies}
#'   \item{EMhistory}{if saved, extract the EM iteration history}
#'   \item{exp_resp}{expected probability of the unique response patterns}
#'   \item{survey.weights}{if supplied, the vector of survey weights used during estimation (NULL if missing)}
#'   \item{converged}{a logical value indicating whether the model terminated within
#'     the convergence criteria}
#'   \item{iterations}{number of iterations it took to reach the convergence criteria}
#'   \item{nest}{number of freely estimated parameters}
#'   \item{parvec}{vector containing uniquely estimated parameters}
#'   \item{vcov}{parameter covariance matrix (associated with parvec)}
#'   \item{condnum}{the condition number of the Hessian (if computed). Otherwise NA}
#'   \item{constrain}{a list of item parameter constraints to indicate which item parameters were equal
#'     during estimation}
#'   \item{Prior}{prior density distribution for the latent traits}
#'   \item{key}{if supplied, the data scoring key}
#'   \item{nfact}{number of latent traits/factors}
#'   \item{nitems}{number of items}
#'   \item{ngroups}{number of groups}
#'   \item{groupNames}{character vector of unique group names}
#'   \item{group}{a character vector indicating the group membership}
#'   \item{invariance}{a character vector indicating \code{invariance} input from \code{\link{multipleGroup}}}
#'   \item{secondordertest}{a logical indicating whether the model passed the second-order test
#'     based on the Hessian matrix. Indicates whether model is a potential local maximum solution}
#'   \item{SEMconv}{logical; check whether the supplemented EM information matrix converged. Will be \code{NA}
#'     if not applicable}
#'   \item{time}{estimation time, broken into different sections}
#' }
#'
#' @aliases extract.mirt
#' @export extract.mirt
#' @param x mirt model of class 'SingleGroupClass', 'MultipleGroupClass', 'MixedClass' or
#'   'DiscreteGroupClass'
#' @param what a string indicating what to extract
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords extract
#' @seealso \code{\link{extract.group}}, \code{\link{extract.item}}, \code{\link{mod2values}}
#' @export extract.item
#' @examples
#'
#' \dontrun{
#' mod <- mirt(Science, 1)
#'
#' extract.mirt(mod, 'logLik')
#' extract.mirt(mod, 'F')
#'
#' #multiple group model
#' grp <- rep(c('G1', 'G2'), each = nrow(Science)/2)
#' mod2 <- multipleGroup(Science, 1, grp)
#'
#' grp1 <- extract.group(mod2, 1) #extract single group model
#' extract.mirt(mod2, 'parvec')
#' extract.mirt(grp1, 'parvec')
#'
#' }
extract.mirt <- function(x, what){
    if(what == 'DIF_coefficients')
        return(attr(x, 'DIF_coefficients'))
    ret <- switch(what,
                  G2 = x@Fit$G2,
                  logLik = x@Fit$logLik,
                  p = x@Fit$p,
                  TLI = x@Fit$TLI,
                  CFI = x@Fit$CFI,
                  RMSEA = x@Fit$RMSEA,
                  df = x@Fit$df,
                  AIC = x@Fit$AIC,
                  HQ = x@Fit$HQ,
                  BIC = x@Fit$BIC,
                  SABIC = x@Fit$SABIC,
                  method = x@Options$method,
                  logPrior = x@Fit$logPrior,
                  F = x@Fit$F,
                  h2 = x@Fit$h2,
                  K = x@Data$K,
                  mins = x@Data$mins,
                  itemtype =  x@Model$itemtype,
                  itemnames = colnames(x@Data$data),
                  parvec = x@Internals$shortpars,
                  vcov = x@vcov,
                  nest = x@Model$nest,
                  constrain = x@Model$constrain,
                  iterations = x@OptimInfo$iter,
                  LLhistory = x@Internals$collectLL,
                  exp_resp = x@Internals$Pl,
                  converged = x@OptimInfo$converged,
                  condnum = x@OptimInfo$condnum,
                  key = x@Internals$key,
                  nfact = x@Model$nfact,
                  nitems = ncol(x@Data$data),
                  ngroups = x@Data$ngroups,
                  model = x@Model$model,
                  group = x@Data$group,
                  Prior = x@Internals$Prior,
                  secondordertest = x@OptimInfo$secondordertest,
                  SEMconv = x@OptimInfo$SEMconv,
                  data = x@Data$data,
                  covdata = x@Data$covdata,
                  tabdata = x@Data$tabdata,
                  freq = x@Data$Freq,
                  tabdatalong = x@Data$tabdatalong,
                  fulldatalong = x@Data$fulldata,
                  groupNames = x@Data$groupNames,
                  time = x@time,
                  survey.weights = x@Internals$survey.weights,
                  rowID = x@Data$rowID,
                  factorNames=x@Model$factorNames,
                  EMhistory=x@Internals$EMhistory,
                  invariance=x@Model$invariance,
                  # undocumented
                  factorNames = x@Model$factorNames,
                  parprior = x@Model$parprior,
                  pars = x@ParObjects$pars,
                  lrPars = x@ParObjects$lrPars,
                  random = x@ParObjects$random,
                  formulas = x@Model$formulas,
                  lrformulas=x@ParObjects$lrPars@formula,
                  itemdesign = x@Data$itemdesign,
                  itemloc = x@Model$itemloc,
                  CUSTOM.IND = x@Internals$CUSTOM.IND,
                  dentype = x@Options$dentype,
                  pis = x@Model$pis,
                  nestpars=x@Model$nestpars,
                  prodlist=x@Model$prodlist,
                  completely_missing=x@Data$completely_missing,
                  stop(sprintf("Could not extract element \'%s\'", what), call.=FALSE))
        ret
}
