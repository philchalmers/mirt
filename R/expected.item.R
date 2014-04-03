#' Function to calculate expected value of item
#'
#' Given an internal mirt object extracted from an estimated model
#' compute the expected value for an item given the ability parameter(s).
#'
#' @aliases expected.item
#' @param x an extracted internal mirt object containing item information
#' @param Theta a vector (unidimensional) or matrix (multidimensional) of latent trait values
#' @param min a constant value added to the expected values indicating the lowest theoretical
#'   category. Default is 0
#' @keywords expected value
#' @export expected.item
#' @seealso \code{\link{extract.item}}, , \code{\link{expected.test}}
#' @examples
#'
#' \dontrun{
#' mod <- mirt(Science, 1)
#' extr.2 <- extract.item(mod, 2)
#' Theta <- matrix(seq(-6,6, length.out=200))
#' expected <- expected.item(extr.2, Theta, min(Science[,1])) #min() of first item
#' head(data.frame(expected, Theta=Theta))
#'
#' }
expected.item <- function(x, Theta, min = 0){
    if(is(Theta, 'vector')) Theta <- as.matrix(Theta)
    if(!is.matrix(Theta)) stop('Theta input must be a matrix')
    if(ncol(Theta) != x@nfact)
        stop('Theta does not have the correct number of dimensions')
    P <- ProbTrace(x=x, Theta=Theta)
    Emat <- matrix(0:(x@ncat-1), nrow(P), ncol(P), byrow = TRUE)
    E <- rowSums(P * Emat) + min
    E
}
