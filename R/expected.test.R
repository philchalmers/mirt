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
#' @param individual logical; return tracelines for individual items?
#' @param which.items an integer vector indicating which items to include in the expected test score. Default
#'   uses all possible items
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
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
#' # use only first two items (i.e., a bundle)
#' bscore <- expected.test(mod, Theta, which.items = 1:2)
#' tail(cbind(Theta, bscore))
#'
#' }
expected.test <- function(x, Theta, group = NULL, mins = TRUE, individual = FALSE, which.items = NULL){
    if(missing(x)) missingMsg('x')
    if(missing(Theta)) missingMsg('Theta')
    pars <- if(is(x, 'MultipleGroupClass')) x@ParObjects$pars[[group]]@ParObjects$pars else x@ParObjects$pars
    K <- extract.mirt(x, 'K')
    if(is.null(which.items) || length(x@Internals$CUSTOM.IND)){
        pick <- 1L:(length(K)+1L)
        which.items <- 1L:length(K)
    } else pick <- c(which.items, length(pars))
    stopifnot(all(pick %in% 1L:(length(K)+1L)))
    pars <- pars[pick]
    itemloc <- c(1L, 1L + cumsum(K[which.items]))
    MINS <- x@Data$mins[which.items]
    trace <- computeItemtrace(pars, Theta, itemloc=itemloc, CUSTOM.IND=x@Internals$CUSTOM.IND)
    if(individual){
        ret <- sapply(1L:length(MINS), function(item, trace, itemloc){
            index <- score <- itemloc[item]:(itemloc[item+1L]-1L)
            score <- score - min(score)
            score %*% t(trace[ ,index])
        }, trace=trace, itemloc=itemloc)
        if(mins) ret <- t(t(ret) + MINS)
    } else {
        score <- do.call(c, lapply(K[which.items], function(x) 0L:(x-1L)))
        ret <- as.numeric(score %*% t(trace))
        if(mins) ret <- ret + sum(MINS)
    }
    return(ret)
}
