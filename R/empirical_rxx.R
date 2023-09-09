#' Function to calculate the empirical (marginal) reliability
#'
#' Given secondary latent trait estimates and their associated standard errors
#' returned from \code{\link{fscores}}, compute the empirical reliability.
#'
#' @aliases empirical_rxx
#' @param Theta_SE a matrix of latent trait estimates returned from \code{\link{fscores}} with the options
#'   \code{full.scores = TRUE} and \code{full.scores.SE = TRUE}
#'
#' @param T_as_X logical; should the observed variance be equal to
#'   \code{var(X) = var(T) + E(E^2)} or \code{var(X) = var(T)} when computing
#'   empirical reliability estimates? Default (\code{FALSE}) uses the former
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
#' itemstats(dat)
#' mod <- mirt(dat)
#'
#' theta_se <- fscores(mod, full.scores.SE = TRUE)
#' empirical_rxx(theta_se)
#'
#' theta_se <- fscores(mod, full.scores.SE = TRUE, method = 'ML')
#' empirical_rxx(theta_se)
#' empirical_rxx(theta_se, T_as_X = TRUE)
#'
#' }
empirical_rxx <- function(Theta_SE, T_as_X = FALSE){
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
    if(T_as_X)
        reliability <- 1 - colMeans(E^2) / (diag(var(T)))
    else
        reliability <- diag(var(T)) / (diag(var(T)) + colMeans(E^2))
    names(reliability) <- colnames(T)
    reliability
}
