#' Extract Empirical Estimating Functions
#'
#' A function for extracting the empirical estimating functions of a fitted
#' \code{\link{mirt}}, \code{\link{multipleGroup}} or \code{\link{bfactor}}
#' model. This is the derivative of the log-likelihood with respect to the
#' parameter vector, evaluated at the observed (case-wise) data. In other
#' words, this function returns the case-wise scores, evaluated at the fitted
#' model parameters. Currently, models fitted via the \code{EM} or \code{BL}
#' method are supported. For the computations, the internal \code{Theta} grid of
#' the model is being used which was already used during the estimation of
#' the model itself along with its matching normalized density.
#'
#' @return An n x k matrix corresponding to n observations and k parameters
#'
#' @aliases estfun.AllModelClass
#' @param x a fitted model object of class \code{SingleGroupClass} or
#'   \code{MultipleGroupClass}
#' @param weights by default, the \code{survey.weights} which were (optionally)
#'   specified when fitting the model are included to calculate the scores.
#'   If specified by the user, this should be a numeric vector of length equal
#'   to the total sample size. Note that if not all cases were weighted equally
#'   when fitting the model, the weights must be corrected by taking their
#'   square root if the scores are being used to compute the outer product of
#'   gradients (OPG) estimate of the variance-covariance matrix (see examples
#'   below).
#' @author Lennart Schneider \email{lennart.sch@@web.de}
#' @keywords scores
#' @seealso \code{\link{mirt}}, \code{\link{multipleGroup}},
#'   \code{\link{bfactor}}
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # fit a 2PL on the LSAT7 data and get the scores
#' mod1 <- mirt(expand.table(LSAT7), 1, SE = TRUE, SE.type = "crossprod")
#' sc1 <- estfun.AllModelClass(mod1)
#' # get the gradient
#' colSums(sc1)
#' # calculate the OPG estimate of the variance-covariance matrix "by hand"
#' vc1 <- vcov(mod1)
#' all.equal(crossprod(sc1), chol2inv(chol(vc1)), check.attributes = FALSE)
#' 
#' # fit a multiple group 2PL and do the same as above
#' group <- rep(c("G1", "G2"), 500)
#' mod2 <- multipleGroup(expand.table(LSAT7), 1, group, SE = TRUE,
#'   SE.type = "crossprod")
#' sc2 <- estfun.AllModelClass(mod2)
#' colSums(sc2)
#' vc2 <- vcov(mod2)
#' all.equal(crossprod(sc2), chol2inv(chol(vc2)), check.attributes = FALSE)
#'
#' # fit a bifactor model with 2 specific factors and do the same as above
#' mod3 <- bfactor(expand.table(LSAT7), c(2, 2, 1, 1, 2), SE = TRUE,
#'   SE.type = "crossprod")
#' sc3 <- estfun.AllModelClass(mod3)
#' colSums(sc3)
#' vc3 <- vcov(mod3)
#' all.equal(crossprod(sc3), chol2inv(chol(vc3)), check.attributes = FALSE)
#'
#' # fit a 2PL not weighting all cases equally
#' survey.weights <- c(rep(2, sum(LSAT7$freq) / 2), rep(1, sum(LSAT7$freq) / 2))
#' survey.weights <- survey.weights / sum(survey.weights) * sum(LSAT7$freq)
#' mod4 <- mirt(expand.table(LSAT7), 1, SE = TRUE, SE.type = "crossprod",
#'   survey.weights = survey.weights)
#' sc4 <- estfun.AllModelClass(mod4,
#'   weights = extract.mirt(mod4, "survey.weights"))
#' # get the gradient
#' colSums(sc4)
#' # to calculate the OPG estimate of the variance-covariance matrix "by hand",
#' # the weights must be adjusted by taking their square root
#' sc4_crp <- estfun.AllModelClass(mod4,
#'   weights = sqrt(extract.mirt(mod4, "survey.weights")))
#' vc4 <- vcov(mod4)
#' all.equal(crossprod(sc4_crp), chol2inv(chol(vc4)), check.attributes = FALSE)
#'
#'}

estfun.AllModelClass <- function(x, weights = extract.mirt(x, "survey.weights"))
{
  ## check class
  stopifnot(class(x) %in% c("SingleGroupClass", "MultipleGroupClass"))
  ## check estimation method
  stopifnot(x@Options$method %in% c("EM", "BL"))
  ## check latent regression
  if(length(x@Model$lrPars)) {
    stop("Scores computations currently not supported for latent regression estimates.")
  }
  ## check items
  CUSTOM.IND <- x@Internals$CUSTOM.IND
  SLOW.IND <- x@Internals$SLOW.IND
  whichitems <- unique(c(CUSTOM.IND, SLOW.IND))
  if(length(whichitems)) {
    stop("Scores computations currently not supported for at least one of the supplied items.")
  }
  ## get relevant model info
  constrain <- x@Model$constrain
  Data <- x@Data
  group <- x@Data$group
  groupNames <- x@Data$groupNames
  itemloc <- x@Model$itemloc
  nitems <- x@Data$nitems
  ngroups <- x@Data$ngroups
  Prior <- x@Internals$Prior
  prior <- x@Internals$bfactor$prior
  Priorbetween <- x@Internals$bfactor$Priorbetween
  sitems <- x@Internals$bfactor$sitems
  sw <- x@Internals$survey.weights
  Theta <- x@Model$Theta
  isbifactor <- length(Priorbetween[[1L]]) > 0L
  ## check if not bifactor
  if(!isbifactor) {
    sitems <- as.matrix(0)
    prior <- Priorbetween <- list(matrix(0))
  }
  pars <-
  if(ngroups == 1L) {
      list(x@ParObjects$pars)
  } else {
      lapply(x@ParObjects$pars, function(pr) pr@ParObjects$pars)
  }
  Ls <- makeLmats(pars, constrain)
  L <- Ls$L
  redun_constr <- Ls$redun_constr
  epars <- mod2values(x)
  epars$group <- factor(epars$group, levels = groupNames)
  eparsgroup <- split(epars, epars$group)
  sel <- lapply(eparsgroup, function(grp) grp$parnum[grp$est])
  ## constrains: cb = between groups, cw = within groups
  if(length(constrain)) {
    if(ngroups > 1L) {
      constraingroup <- lapply(constrain, function(cst) epars$group[cst])
      cb <- which(sapply(constraingroup, function(cst) !any(duplicated(cst))))
      if(length(cb) == 0L) {
        cw <- seq_along(constrain)
      } else {
        cw <- seq_along(constrain)[-cb]
      }
    } else {
      cb <- NULL
      cw <- seq_along(constrain)
    }
  } else {
    cb <- NULL
    cw <- NULL
  }
  ## get relevant response patterns
  datcollpsd <-
  sapply(seq_len(dim(Data$data)[1L]), function(rw) {
    paste(Data$data[rw, ], collapse = "")
  })
  tabdatacollpsd <-
  sapply(seq_len(dim(Data$tabdata)[1L]), function(tb) {
    paste(Data$tabdata[tb, ], collapse = "")
  })
  datpats <- match(datcollpsd, tabdatacollpsd)
  ## compute gradient for each pattern
  gitemtrace <- vector("list", ngroups)
  for(g in seq_len(ngroups)) {
    gitemtrace[[g]] <- computeItemtrace(pars = pars[[g]],
      Theta = Theta, itemloc = itemloc, CUSTOM.IND = CUSTOM.IND)
    gp <- ExtractGroupPars(pars[[g]][[nitems + 1L]])
    pars[[g]][[nitems + 1L]]@mu <- gp$gmeans
    pars[[g]][[nitems + 1L]]@sig <- gp$gcov
    pars[[g]][[nitems + 1L]]@invsig <- solve(gp$gcov)
    pars[[g]][[nitems + 1L]]@meanTheta <- colMeans(Theta)
  }
  npars <- ncol(L)
  gPrior <- t(do.call(rbind, Prior))
  gradient <- .Call("computeGradient", pars, Theta, gPrior, prior,
    do.call(rbind, Priorbetween), Data$tabdatalong, sitems, itemloc,
    gitemtrace, npars, isbifactor)
  ## select matching gradient
  scores <- gradient[datpats, ]
  ## apply constrain handling
  if(ngroups > 1L) {
    for(g in seq_len(ngroups)) {
      scores[group == groupNames[g], -sel[[g]]] <- 0
    }
  }
  if(length(cw)) {
    for(cst in cw) {
      scores[, constrain[[cst]][1L]] <-
      scores[, constrain[[cst]][1L]] +
      rowSums(as.matrix(scores[, constrain[[cst]][-1L]]))
    }
  }
  if(length(cb)) {
    for(cst in cb) {
      for(g in seq_len(ngroups)) {
        tmp <- which(epars$group[constrain[[cst]]] == groupNames[g])
        if(length(tmp)) {
          scores[group == groupNames[g], constrain[[cst]][1L]] <-
          scores[group == groupNames[g], constrain[[cst]][tmp]]
        }
      }
    }
  }
  ## select uniquely estimated model parameters
  sel <- unlist(sel)
  scores <- scores[, sel[!(sel %in% which(redun_constr))]]
  ## handle survey weights
  if(is.null(weights)) {
    weights <- rep.int(1L, Data$N)
  }
  scores <- scores * weights
  ## apply proper naming
  colnames(scores) <-
  if(x@Options$SE) {
    colnames(x@vcov)
  } else {
    NULL
  }
  return(scores)
}
