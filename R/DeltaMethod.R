#' Numerical derivative version of delta method
#'
#' Delta method using numerical derivatives
#' (via \code{\link{numerical_deriv}} for the provided function.
#' Convenient when function is easier to automatic grammatically
#' rather than using an explicit formula or math expression, though
#' is also useful for checking results.
#'
#' @param fn a function or list of functions. Must be of the form
#'   \code{fn(par, ...)}
#' @param par numerical vector (typically MLEs)
#' @param acov numeric matrix of the MLE ACO
#' @param ... additional arguments passed to fn
#'
#' @export
#' @return returns a list of the transformed parameters, ACOV,
#'   and SEs. If fn was a list return will be a list of lists
#'
#' @examples
#'
#' # Slightly modified example from ?msm::deltamethod
#' # Multiple linear regression, E(y) = alpha + beta1 x + beta2 g
#' x <- 1:100
#' g <- rep(0:1, each=50)
#' y <- rnorm(100, 4*x, 5)
#' toy.lm <- lm(y ~ x + g)
#' estmean <- coef(toy.lm)
#' estvar <- summary(toy.lm)$cov.unscaled * summary(toy.lm)$sigma^2
#'
#' ## Estimate of (1 / (alphahat + betahat))
#' 1 / (estmean[1] + estmean[2])
#' 1 / (estmean[1] + estmean[2] + estmean[3])
#'
#' ## Approximate standard error (uncomment to check)
#' # msm::deltamethod (~ 1 / (x1 + x2), estmean, estvar)
#' # msm::deltamethod (~ 1 / (x1 + x2 + x3), estmean, estvar)
#'
#' # with DeltaMethod
#' fn <- function(par) 1 / sum(par[1:2])
#' DeltaMethod(fn, estmean, estvar)$se
#'
#' # index argument for more flexibility
#' fn <- function(par, index) 1 / sum(par[index])
#' DeltaMethod(fn, estmean, estvar, index=1:2)$se
#' DeltaMethod(fn, estmean, estvar, index=1:3)$se
#'
#' # as list of functions
#' fn1 <- function(par) 1 / sum(par[1:2])
#' fn2 <- function(par) 1 / sum(par[1:3])
#' out <- DeltaMethod(list(fn1, fn2), estmean, estvar)
#' c(out[[1]]$se, out[[2]]$se)
#'
DeltaMethod <- function(fn, par, acov, ...){
    if(!is.list(fn)) fn <- list(fn)
    ret <- lapply(1:length(fn), function(i, par){
        vals <- unname(fn[[i]](par, ...))
        dfn <- matrix(numerical_deriv(par, fn[[i]], ...), 1)
        acov_vals <- dfn %*% acov %*% t(dfn)
        SEs <- sqrt(acov_vals)
        list(fn_par=vals, acov=acov_vals, se=sqrt(diag(acov_vals)))
    }, par=par)
    if(length(ret) == 1) ret <- ret[[1]]
    ret
}
