#' Function to calculate the marginal reliability
#'
#' Given an estimated model and a prior density function, compute the marginal reliability. This is only
#' available for unidimensional tests.
#'
#' @aliases marginal_rxx
#' @param mod an object of class \code{'SingleGroupClass'}
#' @param density a density function to use for integration. Default assumes the latent traits are from a
#'   normal (Gaussian) distribution
#' @param theta_lim a vector containing the range of integration
#' @param ... additional arguments passed to the density function
#' @keywords reliability
#' @export marginal_rxx
#' @seealso \code{\link{extract.group}}, \code{\link{testinfo}}
#' @examples
#'
#' \dontrun{
#'
#' dat <- expand.table(deAyala)
#' mod <- mirt(dat, 1)
#'
#' # marginal estimate
#' marginal_rxx(mod)
#'
#' # empirical estimate (assuming the same prior)
#' fscores(mod, returnER = TRUE)
#'
#' }
marginal_rxx <- function(mod, density = dnorm, theta_lim = c(-6,6), ...){
    stopifnot(mod@nfact == 1L)
    stopifnot(is(mod, 'SingleGroupClass'))
    Theta <- seq(theta_lim[1L], theta_lim[2L], length.out = 1000L)
    TI <- testinfo(mod, Theta)
    den <- density(Theta, ...)
    den <- den / sum(den)
    ret <- sum(TI / (TI + 1) * den)
    ret
}
