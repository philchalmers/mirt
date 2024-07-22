#' Function to calculate the marginal reliability
#'
#' Given an estimated model and a prior density function, compute the marginal reliability
#' (Thissen and Wainer, 2001). This is only available for unidimensional tests.
#'
#' @aliases marginal_rxx
#' @param mod an object of class \code{'SingleGroupClass'}
#' @param density a density function to use for integration. Default assumes the latent traits are from a
#'   normal (Gaussian) distribution
#' @param ... additional arguments passed to the density function
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Thissen, D. and Wainer, H. (2001). Test Scoring. Lawrence Erlbaum Associates.
#'
#' @keywords reliability
#' @export marginal_rxx
#' @seealso \code{\link{empirical_rxx}}, \code{\link{extract.group}}, \code{\link{testinfo}}
#' @examples
#'
#'
#' dat <- expand.table(deAyala)
#' mod <- mirt(dat)
#'
#' # marginal estimate treating item parameters as known
#' marginal_rxx(mod)
#'
#' # compare to alpha
#' itemstats(dat)$overall$alpha
#'
#' \dontrun{
#'
#' # empirical estimate (assuming the same prior)
#' fscores(mod, returnER = TRUE)
#'
#' # empirical rxx the alternative way, given theta scores and SEs
#' fs <- fscores(mod, full.scores.SE=TRUE)
#' head(fs)
#' empirical_rxx(fs)
#'
#' }
marginal_rxx <- function(mod, density = dnorm, ...){
    stopifnot(extract.mirt(mod, 'nfact') == 1L)
    stopifnot(is(mod, 'SingleGroupClass'))
    cfs <- coef(mod, simplify=TRUE)
    var_theta <- as.vector(cfs$cov)
    fn <- function(theta, mod, den, ...){
        TI <- testinfo(mod, matrix(theta))
        TI / (TI + var_theta) * den(theta, ...)
    }
    integrate(fn, lower = -Inf, upper=Inf, mod=mod, den=density, ...)$value
}
