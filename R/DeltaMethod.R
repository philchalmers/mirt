#' Numerical derivative version of delta method
#'
#' Delta method using numerical derivatives
#' (via \code{\link{numerical_deriv}}) for the provided function.
#' Convenient when target transformation function is
#' easier to automate programmatically
#' instead of using explicit formula or math expressions. Can
#' also be useful for checking analytic results.
#'
#' @param fn a function specifying the type of
#'   transformation to make for each new parameter of interest.
#'   Must be of the form \code{fn(par, ...)}, or more simply
#'   \code{fn(par)}, and return a numeric vector with one element
#' @param par numerical vector passed to \code{fn(par)}
#'   (typically vector of MLEs)
#' @param acov numeric matrix for the ACOV of the MLEs
#' @param ... additional arguments passed to fn
#'
#' @export
#' @return returns a list of the transformed parameters, ACOV,
#'   and SEs
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
#' estvar <- vcov(toy.lm)
#'
#' # Estimate of (1 / (b0 + b1)) and (1 / (b0 + b1 + b2))
#' 1 / (estmean[1] + estmean[2])
#' 1 / (estmean[1] + estmean[2] + estmean[3])
#'
#' ## Approximate standard error (uncomment to check)
#' # msm::deltamethod (~ 1 / (x1 + x2), estmean, estvar)
#' # msm::deltamethod (~ 1 / (x1 + x2 + x3), estmean, estvar)
#'
#' # with DeltaMethod
#' fn <- function(par) 1 / sum(par[1:2])
#' fn2 <- function(par) 1 / sum(par[1:3])
#' DeltaMethod(fn, estmean, estvar)$se
#' DeltaMethod(fn2, estmean, estvar)$se
#'
#' # index argument for easier flexibility
#' fn <- function(par, index) 1 / sum(par[index])
#' DeltaMethod(fn, estmean, estvar, index=1:2)$se
#' DeltaMethod(fn, estmean, estvar, index=1:3)$se
#'
DeltaMethod <- function(fn, par, acov, ...){
    stopifnot(length(par) == ncol(acov))
    vals <- unname(fn(par, ...))
    dfn <- matrix(numerical_deriv(par, fn, ...), 1)
    acov_vals <- dfn %*% acov %*% t(dfn)
    ret <- list(fn_par=vals, acov=acov_vals,
                se=sqrt(diag(acov_vals)))
    ret
}
