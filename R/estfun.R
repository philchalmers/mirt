#' Extract Empirical Estimating Functions
#'
#' A function for extracting the empirical estimating functions of a fitted
#' \code{\link{mirt}} or \code{\link{multipleGroup}} model. This is the
#' derivative of the log-likelihood with respect to the parameter vector,
#' evaluated at the observed (case-wise) data. In other words, this function
#' returns the case-wise scores, evaluated at the fitted model parameters.
#'
#' @return An n x k matrix corresponding to n observations and k parameters.
#'
#' @aliases estfun.AllModelClass
#' @param object a computed model object of class \code{SingleGroupClass} or
#'   \code{MultipleGroupClass}
#' @author Lennart Schneider \email{lennart.sch@@web.de}
#' @keywords empirical estimating functions
#' @seealso \code{\link{mirt}}, \code{\link{multipleGroup}}
#' @export estfun.AllModelClass
#'
#' @examples
#'
#' \dontrun{
#' mod <- mirt(expand.table(LSAT7), 1, SE = TRUE, SE.type = "crossprod")
#' (sc <- estfun.AllModelClass(mod))
#' colSums(sc)
#' vc <- mod@vcov
#' all.equal(crossprod(sc), chol2inv(chol(vc)), check.attributes = FALSE)
#'
#'}

estfun.AllModelClass <- function(object) {
  ### only allow for EM and BL?
  ### only allow for models of class SingleGroupClass and MultipleGroupClass?
  ### get relevant model info
  constrain <- object@Model$constrain
  CUSTOM.IND <- object@Internals$CUSTOM.IND
  Data <- object@Data
  full <- object@Options$full
  group <- object@Data$group
  groupNames <- object@Data$groupNames
  itemloc <- object@Model$itemloc
  nitems <- object@Data$nitems
  ngroups <- object@Data$ngroups
  Prior <- object@Internals$Prior
  prior <- object@Internals$bfactor$prior
  Priorbetween <- object@Internals$bfactor$Priorbetween
  sitems <- object@Internals$bfactor$sitems
  SLOW.IND <- object@Internals$SLOW.IND
  tabdata <- Data$tabdatalong
  Theta <- object@Model$Theta
  isbifactor <- length(Priorbetween[[1L]]) > 0L
  ### check items
  whichitems <- unique(c(CUSTOM.IND, SLOW.IND))
  if(length(whichitems)) {
    stop("Scores computations currently not supported for at least one of the supplied items.")
  }
  ### check if not bifactor
  if(!isbifactor) {
    sitems <- as.matrix(0)
    prior <- Priorbetween <- list(matrix(0))
  }
  ### get rest of relevant model info
  pars <-
  if(ngroups == 1L) {
      list(object@ParObjects$pars)
  } else {
      lapply(object@ParObjects$pars, function(x) x@ParObjects$pars)
  }
  Ls <- makeLmats(pars, constrain)
  L <- Ls$L
  redun_constr <- Ls$redun_constr
  epars <- mod2values(object)
  eparsgroup <- split(epars, epars$group)
  sel <- lapply(eparsgroup, function(x) x$parnum[x$est])
  ### handle constrains; cb = between groups, cw = within groups
  if(length(constrain)) {
    if(ngroups) {
      constraingroup <- lapply(constrain, function(x) epars$group[x])
      cb <- which(sapply(constraingroup, function(x) !any(duplicated(x))))
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
  ### get relevant response patterns
  datcollpsd <-
  sapply(seq_len(dim(Data$data)[1L]), function(x) {
    paste(Data$data[x, ], collapse = "")
  })
  tabdatacollpsd <-
  sapply(seq_len(dim(Data$tabdata)[1L]), function(x) {
    paste(Data$tabdata[x, ], collapse = "")
  })
  datpats <- match(datcollpsd, tabdatacollpsd)
  ### compute gradient for each pattern
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
  rs <- do.call(rbind, Data$Freq)
  gradient <- .Call("computeGradient", pars, Theta, gPrior, prior,
    do.call(rbind, Priorbetween), Data$tabdatalong, rs, sitems, itemloc,
    gitemtrace, npars, isbifactor)
  ### select matching gradient
  scores <- gradient[datpats, ]
  ### apply constrain handling
  if(ngroups > 1L) {
    for(g in seq_len(ngroups)) {
      scores[group == groupNames[g], -sel[[g]]] <- 0
    }
  }
  if(length(cw)) {
    for(x in cw) {
      scores[, constrain[[x]][1L]] <-
      scores[, constrain[[x]][1L]] +
      rowSums(as.matrix(scores[, constrain[[x]][-1L]]))
    }
  }
  if(length(cb)) {
    for(x in cb) {
      for(g in seq_len(ngroups)) {
        tmp <- which(epars$group[constrain[[x]]] == groupNames[g])
          if(length(tmp)) {
            scores[group == groupNames[g], constrain[[x]][1L]] <-
            scores[group == groupNames[g], constrain[[x]][tmp]]
        }
      }
    }
  }
  ### select uniquely estimated model parameters
  sel <- unlist(sel)
  scores <- scores[, sel[!(sel %in% which(redun_constr))]]
  colnames(scores) <-
  if(all(!is.na(object@vcov))) {
    colnames(object@vcov)
  } else {
    NULL
  }
  return(scores)
}

