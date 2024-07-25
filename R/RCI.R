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
#' @param SEM.pre standard error of measurement for the pretest. This can be used instead of
#'   \code{rxx.pre} and \code{SD.pre}
#' @param SEM.post (optional) standard error of measurement for the post-test.
#'   Using this will create a pooled version of the SEM; otherwise, \code{SEM.post = SEM.pre}
#' @param Fisher logical; use the Fisher/expected information function to compute the
#'   SE terms? If \code{FALSE} the SE information will be extracted from the select
#'   \code{\link{fscores}} method (default). Only applicable for unidimensional models
#' @param verbose logical; include extra information in the printout?
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
#' # CTT version by omitting IRT model
#'     # Requires either sample or population SEM's as input
#' (istats <- itemstats(dat_pre)$overall)
#' SEM.alpha <- istats$SEM.alpha    # SEM estimate of dat_pre
#'
#' # assumes SEM.post = SEM.pre
#' RCI(predat = dat_pre, postdat = dat_post, SEM.pre=SEM.alpha)
#'
#' # include cutoffs
#' RCI(predat = dat_pre, postdat = dat_post, SEM.pre=SEM.alpha,
#'     cutoffs=c(-1.96, 1.96))
#'
#' # allows SEM.post != SEM.pre
#' (istats.post <- itemstats(dat_post)$overall)
#' SEM.alpha.post <- istats.post$SEM.alpha
#'
#' RCI(predat = dat_pre, postdat = dat_post,
#'    SEM.pre=SEM.alpha, SEM.post=SEM.alpha.post)
#'
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
                Fisher = FALSE, ...){
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
        if(is.null(SEM.pre))
            stop('Must include SEM.pre')
        SEM.pre <- unname(SEM.pre)
        stopifnot(is.numeric(SEM.pre) && length(SEM.pre) == 1L)
        if(is.null(SEM.post)) SEM.post <- SEM.pre
        SEM.post <- unname(SEM.post)
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
                fs_pre[,2L] <- 1/sqrt(testinfo(mod_pre, Theta = fs_pre[,1]))
                fs_post[,2L] <- 1/sqrt(testinfo(mod_post, Theta = fs_post[,1]))
            }
            pse <- sqrt(fs_pre[,2]^2 + fs_post[,2]^2)
            z <- diff/pse
            converge_pre <- converge_post <- rep(TRUE, length(z))
            converge_pre[attr(fs_pre, 'failed2converge')] <- FALSE
            converge_post[attr(fs_post, 'failed2converge')] <- FALSE
            ret <- data.frame(pre.score=fs_pre[,1], post.score=fs_post[,1],
                              converged=converge_pre & converge_post, diff,
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
            converge_pre <- converge_post <- rep(TRUE, length(z))
            converge_pre[attr(fs_pre, 'failed2converge')] <- FALSE
            converge_post[attr(fs_post, 'failed2converge')] <- FALSE
            ret <- data.frame(diff=diff, converged=converge_pre & converge_post,
                              joint, z=z,
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
