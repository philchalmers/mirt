#' Function to calculate item information
#'
#' Given an internal mirt item object extracted by using \code{\link{extract.item}},
#' compute the item information.
#'
#' @aliases iteminfo
#' @param x an extracted internal mirt object containing item information (see \code{\link{extract.item}})
#' @param Theta a vector (unidimensional) or matrix (multidimensional) of latent trait values
#' @param degrees a vector of angles in degrees that are between 0 and 90.
#'   Only applicable when the input object is multidimensional
#' @param total.info logical; return the total information curve for the item? If \code{FALSE},
#'   information curves for each category are returned as a matrix
#' @param multidim_matrix logical; compute the information matrix for each row in \code{Theta}? If \code{Theta}
#'   contains more than 1 row then a list of matrices will be returned, otherwise if \code{Theta} has exactly
#'   one row then a matrix will be returned
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords information
#' @seealso
#' \code{\link{extract.item}}
#' @export iteminfo
#' @examples
#'
#' mod <- mirt(Science, 1)
#' extr.2 <- extract.item(mod, 2)
#' Theta <- matrix(seq(-4,4, by = .1))
#' info.2 <- iteminfo(extr.2, Theta)
#'
#' #do something with the info?
#' plot(Theta, info.2, type = 'l', main = 'Item information')
#'
#' \dontrun{
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
#' # multidimensional
#' mod <- mirt(dat, 2, TOL=1e-2)
#' ii <- extract.item(mod, 1)
#' Theta <- as.matrix(expand.grid(-4:4, -4:4))
#'
#' iteminfo(ii, Theta, degrees=c(45,45)) # equal angle
#' iteminfo(ii, Theta, degrees=c(90,0)) # first dimension only
#'
#' # information matrices
#' iteminfo(ii, Theta, multidim_matrix = TRUE)
#' iteminfo(ii, Theta[1, , drop=FALSE], multidim_matrix = TRUE)
#'
#' }
iteminfo <- function(x, Theta, degrees = NULL, total.info = TRUE, multidim_matrix = FALSE){
    if(is(Theta, 'vector')) Theta <- as.matrix(Theta)
    if(!is.matrix(Theta)) Theta <- as.matrix(Theta)
    if(ncol(Theta) == 1L) degrees <- 0
    use_degrees <- ncol(Theta) > 1L
    if(ncol(Theta) != x@nfact)
        stop('Theta does not have the correct number of dimensions', call.=FALSE)
    if(multidim_matrix){
        ret <- lapply(1L:nrow(Theta), function(ind, x, Theta, total.info){
            ItemInfo2(x=x, Theta=Theta[ind, , drop=FALSE], total.info=total.info, MD=TRUE)
        }, x=x, Theta=Theta, total.info=total.info)
        if(length(ret) == 1L) ret <- ret[[1L]]
        return(ret)
    }
    if(is.null(degrees) && ncol(Theta) != 1L)
        stop('Multidimensional information requires pre-specified angles in degrees',
             call.=FALSE)
    cosangle <- cos(d2r(degrees))
    info <- if(use_degrees){
        if(ncol(Theta) != length(cosangle))
            stop('length of the degrees vector not equal to the number of factors', call.=FALSE)
        ItemInfo(x=x, Theta=Theta, cosangle=cosangle, total.info=total.info)
    } else ItemInfo2(x=x, Theta=Theta, total.info=total.info)
    info
}
