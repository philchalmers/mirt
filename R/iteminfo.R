#' Function to calculate item information
#'
#' Given an internal mirt item object extracted by using \code{\link{extract.item}},
#' compute the item information.
#'
#' @aliases iteminfo
#' @param x an extracted internal mirt object containing item information
#' @param Theta a vector (unidimensional) or matrix (multidimensional) of latent trait values
#' @param degrees a vector of angles in degrees that are between 0 and 90 that jointly sum to 90.
#'   Only applicable when the input object is multidimensional
#' @param total.info logical; return the total information curve for the item? If \code{FALSE},
#'   information curves for each category are returned as a matrix
#' @param use_degrees logical; use \code{degrees} argument for multidimensional information?
#'   If \code{FALSE}, information will be calculated without reference to any angle
#'
#' @keywords information
#' @seealso
#' \code{\link{extract.item}}
#' @export iteminfo
#' @examples
#'
#' \dontrun{
#' mod <- mirt(Science, 1)
#' extr.2 <- extract.item(mod, 2)
#' Theta <- matrix(seq(-4,4, by = .1))
#' info.2 <- iteminfo(extr.2, Theta)
#'
#' #do something with the info?
#' plot(Theta, info.2, type = 'l', main = 'Item information')
#'
#' #category information curves
#' cat.info <- iteminfo(extr.2, Theta, total.info = FALSE)
#' plot(Theta, cat.info[,1], type = 'l', ylim = c(0, max(cat.info)),
#'      ylab = 'info', main = 'Category information')
#' for(i in 2:ncol(cat.info))
#'    lines(Theta, cat.info[,i], col = i)
#'
#' ## Customized test information plot
#' T1 <- T2 <- 0
#' dat <- expand.table(LSAT7)
#' mod1 <- mirt(dat, 1)
#' mod2 <- mirt(dat, 1, 'Rasch')
#' for(i in 1:5){
#'   T1 <- T1 + iteminfo(extract.item(mod1, i), Theta)
#'   T2 <- T2 + iteminfo(extract.item(mod2, i), Theta)
#' }
#' plot(Theta, T2/T1, type = 'l', ylab = 'Relative Test Information', las = 1)
#' lines(Theta, T1/T1, col = 'red')
#'
#' }
iteminfo <- function(x, Theta, degrees = NULL, total.info = TRUE, use_degrees = TRUE){
    if(is(Theta, 'vector')) Theta <- as.matrix(Theta)
    if(!is.matrix(Theta)) stop('Theta input must be a matrix', call.=FALSE)
    if(is.null(degrees) && ncol(Theta) == 1L) degrees <- 0
    if(is.null(degrees) && ncol(Theta) != 1L)
        stop('Multidimensional information requires prespecified angles in degrees that sum to 90',
             call.=FALSE)
    if(ncol(Theta) != x@nfact)
        stop('Theta does not have the correct number of dimensions', call.=FALSE)
    cosangle <- cos(d2r(degrees))
    info <- if(use_degrees){
        stopifnot(closeEnough(sum(cosangle^2), 1-1e-10, 1+1e-10))
        ItemInfo(x=x, Theta=Theta, cosangle=cosangle, total.info=total.info)
    } else ItemInfo2(x=x, Theta=Theta, total.info=total.info)
    info
}
