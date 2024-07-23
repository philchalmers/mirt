#' Model-based Reliable Change Index
#'
#' Computes an IRT version of the "reliable change index" (RCI) proposed by
#' Jacobson and Traux (1991) but modified to use IRT information about scores
#' and measurement error (see Jabrayilov, Emons, and Sijtsma (2016).
#' Main benefit of the IRT approach is the inclusion
#' of response pattern information in the pre/post data score estimates, as well
#' as conditional standard error of measurement information.
#'
#' @param mod_pre single-group model fitted by \code{\link{mirt}}. If not supplied the
#'  information will be extracted from the data input objects to compute the classical
#'  test theory version of the RCI statistics
#' @param mod_post (optional) IRT model for post-test if different from pre-test;
#'  otherwise, the pre-test model will be used
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
#' @param rxx.post (optional) same as \code{rxx.pre}, but for post-test data. Using this
#'   will create a pooled version of the SEM
#' @param SD.pre standard deviation of pretest. If not supplied will be computed from \code{predat}.
#'   Required when \code{rxx.pre} is specified
#' @param SD.post (optional) same as \code{SD.pre}, but for the post-test data
#' @param SEM.pre standard error of measurement for the pretest. This can be used instead of
#'   \code{rxx.pre} and \code{SD.pre}
#' @param SEM.post (optional) standard error of measurement for the post-test. This can be used instead of
#'   \code{rxx.post} and \code{SD.post}. Using this will create a pooled version of the SEM
#' @param Fisher logical; use the Fisher/expected information function to compute the
#'   SE terms? If \code{FALSE} the SE information will be extracted from the select
#'   \code{\link{fscores}} method (default). Only applicable for unidimensional models
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
#'
#' Jabrayilov, R. , Emons, W. H. M., & Sijtsma, K. (2016). Comparison of
#' Classical Test Theory and Item Response Theory in Individual Change Assessment.
#' \emph{Applied Psychological Measurement, 40} (8), 559-572.
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
#' # WLE estimator with Fisher information for SE (see Jabrayilov et al. 2016)
#' RCI(mod, predat = dat_pre[1,], postdat = dat_post[1,],
#'     method = 'WLE', Fisher = TRUE)
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
#' # equivalent, but using SEM instead
#' SEM.pre <- 2 * sqrt(1 - .6)
#' SEM.post <- 3 * sqrt(1 - .6)
#' RCI(predat = dat_pre[1:6,], postdat = dat_post[1:6,],
#'     SEM.pre=SEM.pre, SEM.post=SEM.post,
#'     cutoffs = c(-1.96, 1.96))
#'
#' # use just pre-test rxx for computing SEM
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
#' mod <- mirt(dat)
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
RCI <- function(mod_pre, predat, postdat,
                mod_post = mod_pre, cutoffs = NULL,
                SEM.pre = NULL, SEM.post = NULL,
                rxx.pre = NULL, rxx.post = NULL,
                SD.pre = NULL, SD.post = NULL, Fisher = FALSE, ...){
    if(!is.null(cutoffs))
        stopifnot(length(cutoffs) == 2)
    nfact <- 1L
    if(missing(mod_pre)){
        if(is.vector(predat))
            predat <- matrix(predat, 1L)
        if(is.vector(postdat))
            postdat <- matrix(postdat, 1L)
        TS_pre <- rowSums(predat)
        TS_post <- rowSums(postdat)
        if(is.null(SEM.post)) SEM.post <- SEM.pre
        if(is.null(SEM.pre)){
            if(is.null(rxx.pre)) rxx.pre <- CA(predat)
            else if(is.null(SD.pre)) stop('rxx.pre and SD.pre must both be included')
            if(is.null(SD.pre)) SD.pre <- sd(TS_pre)
            if(is.null(SD.post)) SD.post <- sd(TS_post)
            if(is.null(rxx.post)) rxx.post <- CA(postdat)
            SEM.pre <- as.numeric(SD.pre * sqrt(1 - rxx.pre))
            SEM.post <- as.numeric(SD.post * sqrt(1 - rxx.post))
        } else {
            stopifnot(is.numeric(SEM.pre) && length(SEM.pre) == 1L)
            if(!is.null(rxx.pre))
                stop('Please use either SEM or rxx/SD inputs, not both')
        }
        SEM <- sqrt(SEM.pre^2 + SEM.post^2)
        diff <- TS_post - TS_pre
        z_JCI <- diff / SEM
        ret <- data.frame(pre.score=TS_pre, post.score=TS_post, diff,
                          SEM=SEM, z=z_JCI,
                          p=pnorm(abs(z_JCI), lower.tail = FALSE)*2)
    } else {
        if(is.null(mod_post)) mod_post <- mod_pre
        nfact <- extract.mirt(mod_pre, 'nfact')
        if(nfact == 1L){
            fs_pre <- fscores(mod_pre, response.pattern = predat, ...)
            fs_post <- fscores(mod_post, response.pattern = postdat, ...)
            diff <- fs_post[,1] - fs_pre[,1]
            if(Fisher){
                fs_pre[,2] <- 1/sqrt(testinfo(mod_pre, Theta = fs_pre[,1]))
                fs_post[,2] <- 1/sqrt(testinfo(mod_post, Theta = fs_post[,1]))
            }
            pse <- sqrt(fs_pre[,2]^2 + fs_post[,2]^2)
            z <- diff/pse
            ret <- data.frame(pre.score=fs_pre[,1], post.score=fs_post[,1], diff,
                              SEM=pse, z=z,
                              p=pnorm(abs(z), lower.tail = FALSE)*2)
        } else {
            fs_pre <- fscores(mod_pre, response.pattern=predat, ...)
            fs_post <- fscores(mod_post, response.pattern=postdat, ...)
            fs_acov <- fs_pre_acov <- fscores(mod_pre, response.pattern = predat,
                                              return.acov=TRUE, ...)
            fs_post_acov <- fscores(mod_post, response.pattern=postdat,
                                    return.acov=TRUE, ...)
            for(i in 1L:length(fs_acov))
                fs_acov[[i]] <- fs_pre_acov[[i]] + fs_post_acov[[i]]
            diff <- fs_post[,1L:nfact] - fs_pre[,1L:nfact]
            joint <- vector('list', nrow(diff))
            for(i in 1:nrow(diff))
                joint[[i]] <- wald.test(diff[i,], covB = fs_acov[[i]],
                                        L = diag(ncol(diff)))
            joint <- do.call(rbind, joint)
            SEs <- do.call(rbind, lapply(fs_acov, \(x) sqrt(diag(x))))
            z <- diff/SEs
            ret <- data.frame(diff=diff, joint, z=z,
                              p=pnorm(abs(z), lower.tail = FALSE)*2)
        }
    }

    rownames(ret) <- NULL
    if(!is.null(cutoffs) && nfact == 1L){
        ret$cut_decision <- 'unchanged'
        ret$cut_decision[ret$z > max(cutoffs)] <- 'increased'
        ret$cut_decision[ret$z < min(cutoffs)] <- 'decreased'
    }
    ret <- as.mirt_df(ret)
    ret
}
