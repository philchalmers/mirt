#' Re-scale continuous variables to [0,1] bound
#'
#' Transforms select continuous variables to fall within [0,1] boundary. Requires
#' upper and lower inputs for each continuous variable input.
#'
#' @param dat dataset containing continuous variables to be rescaled
#' @param L a vector of theoretical lower values (e.g., 0) for the respective input variables. Must
#'   be either a single value to applied to each variable, or a vector equal to the number
#'   of elements being transformed
#' @param U a vector of theoretical upper values (e.g., 100) for the respective
#'   input variables. Has the same input structure as \code{L}
#' @param select similar to \code{\link{subset}} to extract specific columns in \code{dat},
#'   however the complete \code{dat} object will be returned with transformations applied
#'   only to these variables. If omitted all variables will be selected for transformation
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' @seealso \code{\link{mirt}}
#' @keywords rescaling
#' @export
#' @examples
#'
#' # dataset with scores all from 0-30
#' dat <- as.data.frame(matrix(runif(5000, 0, 30), ncol=5))
#' head(dat)
#'
#' # rescaled to 0-1 for all columns (default)
#' rescaled <- rescaleContinuous(dat, L = 0, U = 30)
#' head(rescaled)
#'
#' # if U/L differ across variables provide suitable vectors
#' rescaled <- rescaleContinuous(dat,
#'                               L = rep(0, 5), U = rep(30, 5))
#' head(rescaled)
#'
#' # mix of continuous and discrete, rescaling only the former
#' LSAT6 <- expand.table(LSAT6)
#' dat.mix <- data.frame(LSAT6, dat)
#' head(dat.mix)
#'
#' rescaled <- rescaleContinuous(dat.mix, L = 0, U = 30,
#'                               select = V1:V5)
#' head(rescaled)
#'
rescaleContinuous <- function(dat, L, U, select){
    stopifnot(!missing(dat) || !missing(L) || !missing(U))
    vars <- if (missing(select))
        rep_len(TRUE, ncol(dat))
    else {
        nl <- as.list(seq_along(dat))
        names(nl) <- names(dat)
        eval(substitute(select), nl, parent.frame())
    }
    if(length(L) == 1L) L <- rep(L, length(vars))
    if(length(U) == 1L) U <- rep(U, length(vars))
    stopifnot(length(L) == length(vars))
    stopifnot(length(U) == length(vars))
    dat[,vars] <- t((t(dat[,vars,drop=FALSE]) - L) / (U - L))
    if(!all(dat[,vars] > 0 & dat[,vars] < 1))
        warning("rescaled data are not within [0,1]. Consider fixing.", call.=FALSE)
    dat
}
