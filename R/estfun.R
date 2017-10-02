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
#' @param object a fitted model object of class \code{SingleGroupClass} or
#'   \code{MultipleGroupClass}
#' @author Lennart Schneider \email{lennart.sch@@web.de}
#' @keywords scores
#' @seealso \code{\link{mirt}}, \code{\link{multipleGroup}},
#'   \code{\link{bfactor}}
#' @export
#'
#' @examples
#'
#' \dontrun{
#' mod1 <- mirt(expand.table(LSAT7), 1, SE = TRUE, SE.type = "crossprod")
#' (sc1 <- estfun.AllModelClass(mod1))
#' colSums(sc1)
#' vc1 <- vcov(mod1)
#' all.equal(crossprod(sc1), chol2inv(chol(vc1)), check.attributes = FALSE)
#' 
#' group <- rep(c("G1", "G2"), 500)
#' mod2 <- multipleGroup(expand.table(LSAT7), 1, group, SE = TRUE,
#'   SE.type = "crossprod")
#' (sc2 <- estfun.AllModelClass(mod2))
#' colSums(sc2)
#' vc2 <- vcov(mod2)
#' all.equal(crossprod(sc2), chol2inv(chol(vc2)), check.attributes = FALSE)
#'
#'}

estfun.AllModelClass <- function(object) {
  ### check class
  stopifnot(class(object) %in% c("SingleGroupClass", "MultipleGroupClass"))
  ### check estimation method
  stopifnot(object@Options$method  %in% c("EM", "BL"))
  ### check latent regression
  if(length(object@Model$lrPars)) {
    stop("Scores computations currently not supported for latent regression estimates.")
  }
  ### check items
  CUSTOM.IND <- object@Internals$CUSTOM.IND
  SLOW.IND <- object@Internals$SLOW.IND
  whichitems <- unique(c(CUSTOM.IND, SLOW.IND))
  if(length(whichitems)) {
    stop("Scores computations currently not supported for at least one of the supplied items.")
  }
  ### get relevant model info
  constrain <- object@Model$constrain
  Data <- object@Data
  group <- object@Data$group
  groupNames <- object@Data$groupNames
  itemloc <- object@Model$itemloc
  nitems <- object@Data$nitems
  ngroups <- object@Data$ngroups
  Prior <- object@Internals$Prior
  prior <- object@Internals$bfactor$prior
  Priorbetween <- object@Internals$bfactor$Priorbetween
  sitems <- object@Internals$bfactor$sitems
  Theta <- object@Model$Theta
  isbifactor <- length(Priorbetween[[1L]]) > 0L
  ### check if not bifactor
  if(!isbifactor) {
    sitems <- as.matrix(0)
    prior <- Priorbetween <- list(matrix(0))
  }
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
  ### constrains; cb = between groups, cw = within groups
  if(length(constrain)) {
    if(ngroups > 1L) {
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
  if(object@Options$SE) {
    colnames(object@vcov)
  } else {
    NULL
  }
  return(scores)
}
