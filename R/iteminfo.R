#' Function to calculate item information
#'
#' Given an internal mirt object extracted from an estimated model
#' compute the item information.
#' 
#' @aliases iteminfo
#' @param x an extracted internal mirt object containing item information
#' @param Theta a matrix of latent trait values
#' @param degrees a vector of angles in degrees that are between 0 and 90 that jointly sum to 90. 
#' Only applicable when the input object is multidimensional
#' @keywords information
#' @export iteminfo
#' @examples 
#' 
#' \dontrun{
#' mod <- mirt(Science, 1, SE = FALSE)
#' extr.2 <- extract.item(mod, 2)
#' Theta <- matrix(seq(-4,4, by = .1))
#' info.2 <- iteminfo(extr.2, Theta)
#' 
#' #do something with the info?
#' plot(Theta, info.2, type = 'l', Main = 'Item information')
#' } 
iteminfo <- function(x, Theta, degrees = NULL){    
    if(!is.matrix(Theta)) stop('Theta input must be a matrix')
    if(is.null(degrees) && ncol(Theta) == 1) degrees <- 0
    if(is.null(degrees) && ncol(Theta) != 1)
        stop('Multidimensional information requires prespecified angles in degrees that sum to 90')        
    cosangle <- cos(d2r(degrees))
    info <- ItemInfo(x=x, Theta=Theta, cosangle=cosangle)   
    info    
}
