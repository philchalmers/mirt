#' Model-based Reliable Change Index
#'
#' Computes an IRT version of the "reliable change index" (RCI) proposed by
#' Jacobson and Traux (1991) but modified to use IRT information about scores
#' and measurement error. Main benefit of the IRT approach is the inclusion
#' of response pattern information in the pre/post data score estimates, as well
#' as conditional standard error of measurement information.
#'
#' @param mod single-group model fitted by \code{\link{mirt}}. If not supplied the
#'  information will be extracted from the data input objects to compute the classical
#'  test theory version of the RCI statistics
#' @param predat a vector (if one individual) or matrix/data.frame
#'   of response data to be scored, where each individuals' responses are
#'   included in exactly one row
#' @param postdat same as \code{predat}, but with respect to the post/follow-up
#'   measurement
#' @param cutoffs optional vector of length 2 indicating the type of cut-offs to
#'   report (e.g., \code{c(-1.96, 1.96)} reflects the 95 percent z-score type cut-off)
#'
#' @param rxx.pre CTT reliability of pretest. If not supplied will be computed using coefficient
#'  alpha from \code{predat}
#' @param rxx.post same as \code{rxx.pre}, but for post-test data
#' @param SD.pre standard deviation of pretest. If not supplied will be computed from \code{predat}
#' @param SD.post same as \code{SD.pre}, but for the post-test data
#' @param rxx.method for CTT version of RCI, which method to use for pooling the reliability
#'   information. Currently supports \code{'pooled'} to pool the pre-post
#'   reliability estimates (default) or \code{'pre'} for using just the pre-test
#'
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
#' # simulate some data
#' N <- 1000
#' J <- 20     # number of items
#' a <- matrix(rlnorm(J,.2,.3))
#' d <- rnorm(J)
#'
#' theta <- matrix(rnorm(N))
#' dat_pre <- simdata(a, d, itemtype = '2PL', Theta = theta)
#'
#' # first 3 cases decrease by 1/2
#' theta2 <- theta - c(1/2, 1/2, 1/2, numeric(N-3))
#' dat_post <- simdata(a, d, itemtype = '2PL', Theta = theta2)
#'
#' mod <- mirt(dat_pre)
#'
#' # all changes using fitted model from pre data
#' RCI(mod, predat=dat_pre, postdat=dat_post)
#'
#' # single response pattern change using EAP information
#' RCI(mod, predat=dat_pre[1,], postdat=dat_post[1,])
#'
#' # WLE estimator
#' RCI(mod, predat = dat_pre[1,], postdat = dat_post[1,], method = 'WLE')
#'
#' # multiple respondents
#' RCI(mod, predat = dat_pre[1:6,], postdat = dat_post[1:6,])
#'
#' # include large-sample z-type cutoffs
#' RCI(mod, predat = dat_pre[1:6,], postdat = dat_post[1:6,],
#'     cutoffs = c(-1.96, 1.96))
#'
#' ######
#' # CTT version by omitting IRT model (easiest to use complete dataset)
#' RCI(predat = dat_pre, postdat = dat_post)
#'
#' # CTT version with pre-computed information
#' RCI(predat = dat_pre[1:6,], postdat = dat_post[1:6,],
#'     rxx.pre=.6, rxx.post=.6, SD.pre=2, SD.post=3,
#'     cutoffs = c(-1.96, 1.96))
#'
#' # just pre-test rxx
#' RCI(predat = dat_pre[1:6,], postdat = dat_post[1:6,],
#'     rxx.pre=.6, SD.pre=2, rxx.method = 'pre')
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
RCI <- function(mod, predat, postdat, cutoffs = NULL,
                rxx.method = 'pooled',
                rxx.pre = NULL, rxx.post = NULL,
                SD.pre = NULL, SD.post = NULL, ...){
    if(!is.null(cutoffs))
        stopifnot(length(cutoffs) == 2)
    if(missing(mod)){
        stopifnot(rxx.method %in% c('pooled', 'pre'))
        if(is.vector(predat))
            predat <- matrix(predat, 1L)
        if(is.vector(postdat))
            postdat <- matrix(postdat, 1L)
        TS_pre <- rowSums(predat)
        TS_post <- rowSums(postdat)
        if(is.null(SD.pre)) SD.pre <- sd(TS_pre)
        if(is.null(SD.post)) SD.post <- sd(TS_post)
        if(is.null(rxx.pre))
            rxx.pre <- itemstats(predat)$overall$alpha
        if(rxx.method == 'pre')
            rxx.post <- rxx.pre
        if(is.null(rxx.post) && rxx.method == 'pooled')
            rxx.post <- itemstats(postdat)$overall$alpha
        SEM_pre <- as.numeric(SD.pre * sqrt(1 - rxx.pre))
        SEM_post <- as.numeric(SD.post * sqrt(1 - rxx.post))
        if(rxx.method == 'pooled'){
            SEM <- sqrt(SEM_pre^2 + SEM_post^2)
        } else if(rxx.method == 'pre'){
            SEM <- sqrt(2*SEM_pre^2)
        }
        diff <- TS_post - TS_pre
        z_JCI <- diff / SEM
        ret <- data.frame(pre.score=TS_pre, post.score=TS_post, diff,
                          SEM=SEM, z=z_JCI)
    } else {
        if(extract.mirt(mod, 'nfact') == 1L){
            fs_pre <- fscores(mod, response.pattern = predat, ...)
            fs_post <- fscores(mod, response.pattern = postdat, ...)

            diff <- fs_post[,1] - fs_pre[,1]
            pse <- sqrt(fs_pre[,2]^2 + fs_post[,2]^2)
            ret <- data.frame(pre.score=fs_pre[,1], post.score=fs_post[,1], diff,
                              pooled_SEM=pse, z=diff/pse)
        } else {
            stop('not yet supported')
            browser()
        }
    }
    rownames(ret) <- NULL
    if(!is.null(cutoffs) && ncol(ret) == 5L){ # only for unidim
        ret$cut_decision <- 'unchanged'
        ret$cut_decision[ret$z > max(cutoffs)] <- 'increased'
        ret$cut_decision[ret$z < min(cutoffs)] <- 'decreased'
    }
    ret <- as.mirt_df(ret)
    ret
}
