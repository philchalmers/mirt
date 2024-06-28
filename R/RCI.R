#' Model-based Reliable Change Index
#'
#' Computes an IRT version of the "reliable change index" (RCI) proposed by
#' Jacobson and Traux (1991) but modified to use IRT information about scores
#' and measurement error. Main benefit of the IRT approach is the inclusion
#' of response pattern information in the pre/post data score estimates, as well
#' as conditional standard error of measurement information.
#'
#' @param mod single-group model fitted by \code{\link{mirt}}
#' @param predat a vector (if one individual) or matrix/data.frame
#'   of response data to be scored, where each individuals' responses are
#'   included in exactly one row
#' @param postdat same as \code{predat}, but with respect to the post/follow-up
#'   measurement
#' @param cutoffs optional vector of length 2 indicating the type of cut-offs to
#'   report (e.g., \code{c(-1.96, 1.96)} reflects the 95 percent z-score type cut-off)
#' @param ... additional arguments passed to \code{\link{fscores}}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Jacobson, N. S., & Truax, P. (1991). Clinical significance: A statistical approach
#' to defining meaningful change in psychotherapy research. Journal
#' of Consulting and Clinical Psychology, 59, 12-19.
#' @keywords reliable change index
#' @export
#' @examples
#'
#' \dontrun{
#'
#' mod <- mirt(Science, 1)
#'
#' # single response pattern change using EAP information
#' RCI(mod, predat = c(1,2,3,2), postdat = c(1,2,2,1))
#'
#' # WLE estimator
#' RCI(mod, predat = c(1,2,3,2), postdat = c(1,2,2,1), method = 'WLE')
#'
#' # multiple respondents
#' RCI(mod, predat = Science[1:5,], postdat = Science[2:6,])
#'
#' # include large-sample z-type cutoffs
#' RCI(mod, predat = Science[1:5,], postdat = Science[2:6,],
#'     cutoffs = c(-1.96, 1.96))
#'
#' ############################
#' # Example where individuals take completely different item set pre-post
#' #   but prior calibration has been performed to equate the items
#'
#' dat <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#'
#' mod <- mirt(dat, 1)
#'
#' # with N=5 individuals under investigation
#' predat <- postdat <- dat[1:5,]
#' predat[, 17:32] <- NA
#' postdat[, 1:16] <- NA
#'
#' head(predat)
#' head(postdat)
#'
#' RCI(mod, predat, postdat)
#'
#' }
RCI <- function(mod, predat, postdat, cutoffs = NULL, ...){
    stopifnot(extract.mirt(mod, 'nfact') == 1L)
    if(!is.null(cutoffs))
        stopifnot(length(cutoffs) == 2)
    fs_pre <- fscores(mod, response.pattern = predat, ...)
    fs_post <- fscores(mod, response.pattern = postdat, ...)

    diff <- fs_post[,1] - fs_pre[,1]
    pse <- sqrt(fs_pre[,2]^2 + fs_post[,2]^2)
    ret <- data.frame(pre.score=fs_pre[,1], post.score=fs_post[,1], diff,
                      pooled_SEM=pse, z=diff/pse)
    rownames(ret) <- NULL
    if(!is.null(cutoffs)){
        ret$cut_decision <- 'unchanged'
        ret$cut_decision[ret$z > max(cutoffs)] <- 'increased'
        ret$cut_decision[ret$z < min(cutoffs)] <- 'decreased'
    }
    ret <- as.mirt_df(ret)
    ret
}
