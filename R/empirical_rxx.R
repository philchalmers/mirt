#' Function to calculate the empirical (marginal) reliability
#'
#' Given secondary latent trait estimates and their associated standard errors
#' returned from \code{\link{fscores}}, compute the empirical reliability.
#'
#' @aliases empirical_rxx
#' @param Theta_SE a matrix of latent trait estimates returned from \code{\link{fscores}} with the options
#'   \code{full.scores = TRUE} and \code{full.scores.SE = TRUE}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords reliability
#' @export empirical_rxx
#' @seealso \code{\link{fscores}}, \code{\link{marginal_rxx}}
#' @examples
#'
#' \dontrun{
#'
#' dat <- expand.table(deAyala)
#' mod <- mirt(dat, 1)
#'
#' theta_se <- fscores(mod, full.scores.SE = TRUE)
#' empirical_rxx(theta_se)
#'
#' theta_se <- fscores(mod, full.scores.SE = TRUE, method = 'ML')
#' empirical_rxx(theta_se)
#'
#' }
empirical_rxx <- function(Theta_SE){
    stopifnot(is.matrix(Theta_SE))
    stopifnot(ncol(Theta_SE) %% 2 == 0)
    nms <- colnames(Theta_SE)
    is_SE <- grepl('SE_', nms)
    if(!any(is_SE))
        stop('Must use full.scores.SE = TRUE in fscores() when computing estimates', call.=FALSE)
    Theta <- Theta_SE[ ,!is_SE, drop=FALSE]
    SE <- Theta_SE[ ,is_SE, drop=FALSE]
    pick <- !(rowSums(Theta) %in% c(NA, Inf, -Inf))
    T <- Theta[pick, , drop=FALSE]
    E <- SE[pick, , drop=FALSE]
    reliability <- diag(var(T)) / (diag(var(T)) + colMeans(E^2))
    reliability
}
