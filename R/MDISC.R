#' Compute multidimensional discrimination index
#'
#' Returns a vector containing the MDISC values for each item in the model input object (Reckase, 2009).
#'
#' @aliases MDISC
#' @param x an object of class 'SingleGroupClass', or an object of class 'MultipleGroupClass' if a suitable
#'   \code{group} input were supplied
#' @param group group argument to pass to \code{\link{extract.group}} function. Required when the input object is
#'   a multiple-group model
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Reckase, M. D. (2009). Multidimensional Item Response Theory. Springer.
#'
#' @seealso \code{\link{extract.group}}
#'
#' @keywords discrimination
#' @export MDISC
#' @examples
#' \donttest{
#'
#' mod <- mirt(Science, 2)
#' MDISC(mod)
#'
#' }
MDISC <- function(x, group = NULL){
    if(missing(x)) missingMsg('x')
    if(is(x, 'MultipleGroupClass') && is.null(group))
        stop('Input must be a SingleGroupClass object or a MultipleGroupClass object with a suitable group input',
             call.=FALSE)
    if(!is.null(group) && is(x, 'MultipleGroupClass'))
        x <- extract.group(x=x, group=group)
    stopifnot(class(x) == 'SingleGroupClass')
    ret <- numeric(extract.mirt(x, 'nitems'))
    for(i in seq_len(length(ret))){
        item <- extract.item(x, i)
        as <- ExtractLambdas(item)
        ret[i] <- sqrt(as %*% as)
    }
    names(ret) <- extract.mirt(x, 'itemnames')
    ret
}
