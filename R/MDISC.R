#' Compute multidimensional discrimination index
#'
#' Returns a vector containing the MDSIC values (Reckase, 2009).
#'
#' @aliases MDISC
#' @param x an object of class 'SingleGroupClass'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Reckase, M. D. (2009). Multidimensional Item Response Theory. Springer.
#'
#' @seealso \code{\link{extract.group}}
#'
#' @keywords discrimination
#' @export MDISC
#' @examples
#' \dontrun{
#'
#' mod <- mirt(Science, 2)
#' MDISC(mod)
#'
#' }
MDISC <- function(x){
    if(missing(x)) missingMsg('x')
    stopifnot(class(x) == 'SingleGroupClass')
    ret <- numeric(extract.mirt(x, 'nitems'))
    for(i in 1L:length(ret)){
        item <- extract.item(x, i)
        as <- ExtractLambdas(item)
        ret[i] <- as %*% as
    }
    names(ret) <- extract.mirt(x, 'itemnames')
    ret
}
