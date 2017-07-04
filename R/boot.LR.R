#' Parametric bootstrap likleihood-ratio test
#'
#' Given two fitted models, compute a parametric bootstrap test to determine whether
#' the less restrictive models fits significantly better than the more restricted model.
#' Note that this hypothesis test also works when prior parameter distributions are included for
#' either model. Function can be run in parallel after using a stuitable \code{\link{mirtCluster}}
#' definition.
#'
#' @aliases boot.LR
#' @param mod an estimated model object
#' @param mod2 an estimated model object
#' @param R number of parametric bootstraps to use.
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @return a p-value evaluating whether the more restrictive model fits significantly worse
#'   than the less restrictive model
#' @keywords parametric bootstrap
#' @export boot.LR
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @examples
#'
#' \dontrun{
#'
#' #standard
#' dat <- expand.table(LSAT7)
#' mod1 <- mirt(dat, 1)
#' mod2 <- mirt(dat, 1, '3PL')
#'
#' # standard LR test
#' anova(mod1, mod2)
#'
#' # bootstrap LR test (run in parallel to save time)
#' mirtCluster()
#' boot.LR(mod1, mod2, R=200)
#'
#' }
boot.LR <- function(mod, mod2, R = 1000){
    stopifnot(is(mod, 'SingleGroupClass'))
    df1 <- extract.mirt(mod, 'df')
    df2 <- extract.mirt(mod2, 'df')
    if(df1 < df2){
        tmp <- mod
        mod <- mod2
        mod2 <- tmp
    }
    LR <- anova(mod, mod2, verbose=FALSE)$X2[2L]
    results <- mySapply(1L:R, function(r, mod, mod2){
        while(TRUE){
            dat <- simdata(model=mod)
            m0 <- mirt(dat, extract.mirt(mod, 'model'),
                       itemtype=extract.mirt(mod, 'itemtype'),
                       verbose=FALSE, technical=list(warn=FALSE, message=FALSE))
            if(!extract.mirt(m0, 'converged')) next
            m1 <- mirt(dat, extract.mirt(mod2, 'model'),
                       itemtype=extract.mirt(mod2, 'itemtype'),
                       verbose=FALSE, technical=list(warn=FALSE, message=FALSE))
            if(!extract.mirt(m1, 'converged')) next
            lr0 <- anova(m0, m1, verbose=FALSE)$X2[2L]
            if(lr0 <= 0) next
            break
        }
        lr0
    }, mod=mod, mod2=mod2)
    p <- (1 + sum(LR < results, na.rm = TRUE)) / (1 + R)
    p
}
