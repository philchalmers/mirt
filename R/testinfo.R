#' Function to calculate test information
#'
#' Given an estimated model compute the test information.
#'
#' @aliases testinfo
#' @param x an object of class 'SingleGroupClass', or an object of class 'MultipleGroupClass' if a suitable
#'   \code{group} input were supplied
#' @param Theta a matrix of latent trait values
#' @param degrees a vector of angles in degrees that are between 0 and 90.
#'   Only applicable when the input object is multidimensional
#' @param individual logical; return a data.frame of information traceline for each item?
#' @param which.items an integer vector indicating which items to include in the expected information function.
#'   Default uses all possible items
#' @param group group argument to pass to \code{\link{extract.group}} function. Required when the input object is
#'   a multiple-group model
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords information
#' @export testinfo
#' @examples
#'
#' dat <- expand.table(deAyala)
#' (mirt(dat, 1, '2PL', pars = 'values'))
#' mod <- mirt(dat, 1, '2PL', constrain = list(c(1,5,9,13,17)))
#'
#' Theta <- matrix(seq(-4,4,.01))
#' tinfo <- testinfo(mod, Theta)
#' plot(Theta, tinfo, type = 'l')
#'
#' \dontrun{
#'
#' #compare information loss between two tests
#' tinfo_smaller <- testinfo(mod, Theta, which.items = 3:5)
#'
#' #removed item informations
#' plot(Theta, iteminfo(extract.item(mod, 1), Theta), type = 'l')
#' plot(Theta, iteminfo(extract.item(mod, 2), Theta), type = 'l')
#'
#' #most loss of info around -1 when removing items 1 and 2; expected given item info functions
#' plot(Theta, tinfo_smaller - tinfo, type = 'l')
#'
#'
#' }
testinfo <- function(x, Theta, degrees = NULL, group = NULL, individual = FALSE,
                     which.items = 1:extract.mirt(x, 'nitems')){
    if(missing(x)) missingMsg('x')
    if(is(x, 'MultipleGroupClass') && is.null(group))
        stop('Input must be a SingleGroupClass object or a MultipleGroupClass object with a suitable group input',
             call.=FALSE)
    if(!is.null(group) && is(x, 'MultipleGroupClass'))
        x <- extract.group(x=x, group=group)
    if(missing(Theta)) missingMsg('Theta')
    if(!is.matrix(Theta)) Theta <- as.matrix(Theta)
    J <- length(extract.mirt(x, 'K'))
    if(is.null(which.items)) which.items <- 1L:J
    stopifnot(all(which.items <= J & which.items > 0L))
    info <- matrix(0, nrow(Theta), J)
    for(i in which.items){
        item <- extract.item(x, i, group=group)
        info[,i] <- iteminfo(item, Theta=Theta, degrees=degrees)
    }
    if(!individual){
        info <- rowSums(info)
    } else info <- info[,which.items, drop=FALSE]
    info
}
