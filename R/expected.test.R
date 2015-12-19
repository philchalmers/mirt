#' Function to calculate expected test score
#'
#' Given an estimated model compute the expected test score. Returns the expected values in the
#' same form as the data used to estimate the model.
#'
#' @aliases expected.test
#' @param x an estimated mirt object
#' @param Theta a matrix of latent trait values
#' @param group a number signifying which group the item should be extracted from (applies to
#'   'MultipleGroupClass' objects only)
#' @param mins logical; include the minimum value constants in the dataset. If FALSE, the
#'   expected values for each item are determined from the scoring 0:(ncat-1)
#' @keywords expected score
#' @seealso \code{\link{expected.item}}
#' @export expected.test
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(deAyala)
#' model <- 'F = 1-5
#'           CONSTRAIN = (1-5, a1)'
#' mod <- mirt(dat, model)
#'
#' Theta <- matrix(seq(-6,6,.01))
#' tscore <- expected.test(mod, Theta)
#' tail(cbind(Theta, tscore))
#'
#' }
expected.test <- function(x, Theta, group = NULL, mins = TRUE){
    if(missing(x)) missingMsg('x')
    if(missing(Theta)) missingMsg('Theta')
    pars <- if(is(x, 'MultipleGroupClass')) x@ParObjects$pars[[group]]@ParObjects$pars else x@ParObjects$pars
    trace <- computeItemtrace(pars, Theta, itemloc=x@Model$itemloc, CUSTOM.IND=x@Internals$CUSTOM.IND)
    score <- do.call(c, lapply(x@Data$K, function(x) 0L:(x-1L)))
    ret <- as.numeric(score %*% t(trace))
    if(mins) ret <- ret + sum(x@Data$mins)
    return(ret)
}
