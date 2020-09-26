#' Generalized item difficulty summaries
#'
#' Function provides the four generalized item difficulty representations
#' for polytomous response models described by Ali, Chang, and Anderson (2015).
#' These estimates are used to gauge how difficult a polytomous item may be.
#'
#' @param mod a single factor model estimated by \code{\link{mirt}}
#' @param type type of generalized difficulty parameter to report.
#'   Can be \code{'IRF'} to use the item response function (default),
#'   \code{'mean'} to find the average of the difficulty estimates,
#'   \code{'median'} the median of the difficulty estimates, and
#'   \code{'trimmed'} to find the trimmed mean after removing the first
#'   and last difficulty estimates
#' @param interval interval range to search for \code{'IRF'} type
#' @param ... additional arguments to pass to \code{\link{uniroot}}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Ali, U. S., Chang, H.-H., & Anderson, C. J. (2015). \emph{Location indices for ordinal
#' polytomous items based on item response theory} (Research Report No. RR-15-20).
#' Princeton, NJ: Educational Testing Service. http://dx.doi.org/10.1002/ets2.12065
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' @export
#' @examples
#'
#' \dontrun{
#'
#' mod <- mirt(Science, 1)
#' coef(mod, simplify=TRUE, IRTpars = TRUE)$items
#'
#' gen.difficulty(mod)
#' gen.difficulty(mod, type = 'mean')
#'
#' # also works for dichotomous items (though this is unnecessary)
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1)
#' coef(mod, simplify=TRUE, IRTpars = TRUE)$items
#'
#' gen.difficulty(mod)
#' gen.difficulty(mod, type = 'mean')
#'
#' }
gen.difficulty <- function(mod, type = "IRF", interval = c(-30, 30), ...){

    LImean <- function(mod, type){
        cfs <- coef(mod, IRTpars = TRUE)
        cfs <- cfs[2:length(cfs) - 1L]
        ret <- sapply(cfs, function(x){
            bs <- x[colnames(x) %in% c('b', paste0('b', 2L:length(x) - 1))]
            mb <- ifelse(length(bs) > 0,
                         if(type == 'mean') mean(bs)
                         else if(type == 'median') median(bs)
                         else if(type == 'trimmed')
                             mean(sort(bs)[c(-1, -length(bs))]),
                         NA)
            mb
        })
        ret
    }

    LIIRF_1 <- function(item, ...){
        m <- attr(item, 'm')
        fn <- function(theta) expected.item(item, theta) - m
        uniroot(fn, interval=interval, ...)$root
    }

    if(type %in% c("mean", "median", "trimmed")){
        ret <- LImean(mod, type=type)
    } else if(type == "IRF"){
        m <- (extract.mirt(mod, 'K') - 1)/2
        items <- sapply(1:extract.mirt(mod, 'nitems'),
                        function(x){
                           out <- extract.item(mod, x)
                           attr(out, 'm') <- m[x]
                           out
                        })
        ret <- sapply(items, LIIRF_1, ...)
    } else stop('type not supported')
    names(ret) <- extract.mirt(mod, 'itemnames')
    ret
}
