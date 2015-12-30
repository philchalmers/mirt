#' Function to calculate test information
#'
#' Given an estimated model compute the test information.
#'
#' @aliases testinfo
#' @param x an estimated mirt object
#' @param Theta a matrix of latent trait values
#' @param degrees a vector of angles in degrees that are between 0 and 90.
#'   Only applicable when the input object is multidimensional
#' @param group a number signifying which group the item should be extracted from (applies to
#'   'MultipleGroupClass' objects only)
#' @keywords information
#' @export testinfo
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(deAyala)
#' (mirt(dat, 1, '2PL', pars = 'values'))
#' mod <- mirt(dat, 1, '2PL', constrain = list(c(1,5,9,13,17)))
#'
#' Theta <- matrix(seq(-4,4,.01))
#' tinfo <- testinfo(mod, Theta)
#' plot(Theta, tinfo, type = 'l')
#'
#' #compare information loss between two tests
#' dat.smaller <- dat[,-c(1,2)]
#' mod2 <- mirt(dat.smaller, 1, '2PL', constrain = list(c(1,5,9)))
#' tinfo2 <- testinfo(mod2, Theta)
#'
#' #removed item informations
#' plot(Theta, iteminfo(extract.item(mod, 1), Theta), type = 'l')
#' plot(Theta, iteminfo(extract.item(mod, 2), Theta), type = 'l')
#'
#' #most loss of info around -1 when removing items 1 and 2; expected given item info functions
#' plot(Theta, tinfo2 - tinfo, type = 'l')
#'
#'
#' }
testinfo <- function(x, Theta, degrees = NULL, group = NULL){
    if(missing(x)) missingMsg('x')
    if(missing(Theta)) missingMsg('Theta')
    if(!is.matrix(Theta)) Theta <- as.matrix(Theta)
    J <- extract.mirt(x, 'nitems')
    info <- 0
    for(i in 1L:J){
        item <- extract.item(x, i, group=group)
        info <- info + iteminfo(item, Theta=Theta, degrees=degrees)
    }
    return(info)
}
