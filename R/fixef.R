#' Compute latent regression fixed effect expected values
#'
#' Create expected values for fixed effects parameters in latent regression models.
#'
#' @aliases fixef
#' @param x an estimated model object from the \code{\link{mixedmirt}} or \code{\link{mirt}}
#'   function
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{mirt}}, \code{\link{mixedmirt}}
#' @keywords fixed effects
#' @export fixef
#' @examples
#' \dontrun{
#'
#' #simulate data
#' set.seed(1234)
#' N <- 1000
#'
#' # covariates
#' X1 <- rnorm(N); X2 <- rnorm(N)
#' covdata <- data.frame(X1, X2)
#' Theta <- matrix(0.5 * X1 + -1 * X2 + rnorm(N, sd = 0.5))
#'
#' #items and response data
#' a <- matrix(1, 20); d <- matrix(rnorm(20))
#' dat <- simdata(a, d, 1000, itemtype = 'dich', Theta=Theta)
#'
#' #conditional model using X1 and X2 as predictors of Theta
#' mod1 <- mirt(dat, 1, 'Rasch', covdata=covdata, formula = ~ X1 + X2)
#'
#' #latent regression fixed effects (i.e., expected values)
#' fixef(mod1)
#'
#' # with mixedmirt()
#' mod1b <- mixedmirt(dat, covdata, 1, lr.fixed = ~ X1 + X2, fixed = ~ 0 + items)
#' fixef(mod1b)
#'
#' }
fixef <- function(x){
    if(missing(x)) missingMsg('x')
    if(!(is(x, 'MixedClass') || is(x, 'SingleGroupClass')))
        stop('Only applicable to MixedClass and SingleGroupClass objects', call.=FALSE)
    if(!length(x@lrPars))
        stop('No latent regression parameters were defined in the supplied model', call.=FALSE)
    ret <- x@lrPars@X %*% x@lrPars@beta
    ret
}
