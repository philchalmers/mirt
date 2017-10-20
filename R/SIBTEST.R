#' Simultaneous Item Bias Test (SIBTEST)
#'
#' Classical test theory approach to detecting unidirectional and bidirectional (with one
#' crossing location) DIF. This family of statistics is intended for unidimensional tests,
#' and applies a regression-corrected matched-total score approach to quantify the response
#' bias between two groups. Can be used for DIF, DBF, and DTF testing.
#'
#' SIBTEST is similar to the Mantel-Haenszel approach for detecting DIF but uses a regression
#' correction based on the KR-20/coefficient alpha reliability index to correct the observed
#' differences when the latent trait distributions are not equal.
#' Function supports the standard SIBTEST for dichotomous and polytomous data (compensatory) and
#' supports crossing DIF testing (i.e., non-compensatory/non-uniform) using the asymptotic sampling
#' distribution version of the Crossing-SIBTEST (CSIBTEST) statistic described by
#' Chalmers (accepted).
#'
#' @param dat integer dataset to be tested containing dichotomous or polytomous responses
#' @param group a vector indicating group membership
#' @param match_set an integer vector indicating which items to use as the items which are matched
#'   (i.e., contain no DIF). These are analogous to 'achor' items in the likelihood method to locate
#'   DIF. If missing, all items other than the items found in the suspect_set will be used
#' @param focal_name name of the focal group; e.g., 'focal'. If not specified then one will be
#'   selected automatically
#' @param suspect_set an integer vector indicating which items to inspect with SIBTEST. Including only
#'   one value will perform a DIF test, while including more than one will perform a simultaneous
#'   bundle test (DBF); including all non-matched items will perform DTF.
#'   If missing, a simultaneous test using all the items not listed in match_set
#'   will be used (i.e., DTF)
#' @param guess_correction a vector of numbers from 0 to 1 indicating how much to correct the items
#'   for guessing. It's length should be the same as ncol(dat)
#' @param Jmin the minimum number of observations required when splitting the data into focal and
#'   reference groups conditioned on the matched set
#' @param pk_focal logical; using the group weights from the focal group instead of the total
#'   sample? Default is FALSE as per Shealy and Stout's recommendation
#' @param correction logical; apply the composite correction for the difference between focal
#'   composite scores using the true-score regression technique? Default is \code{TRUE},
#'   reflecting Shealy and Stout's method
#' @param details logical; return a data.frame containing the details required to compute SIBTEST?
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords SIBTEST, crossing SIBTEST
#' @aliases SIBTEST
#' @export SIBTEST
#'
#' @references
#'
#' Chalmers, R. P. (accepted). Improving the Crossing-SIBTEST statistic for
#' detecting non-uniform DIF. \emph{Psychometrika}.
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Chang, H. H., Mazzeo, J. & Roussos, L. (1996). DIF for Polytomously Scored Items: An Adaptation of the
#'   SIBTEST Procedure. \emph{Journal of Educational Measurement, 33}, 333-353.
#'
#' Li, H.-H. & Stout, W. (1996). A new procedure for detection of crossing DIF.
#' \emph{Psychometrika, 61}, 647-677.
#'
#' Shealy, R. & Stout, W. (1993). A model-based standardization approach that separates true
#'   bias/DIF from group ability differences and detect test bias/DTF as well as item bias/DIF.
#'   \emph{Psychometrika, 58}, 159-194.
#'
#' @examples
#' \dontrun{
#'
#' set.seed(1234)
#' n <- 30
#' N <- 500
#' a <- matrix(1, n)
#' d <- matrix(rnorm(n), n)
#' group <- c(rep('reference', N), rep('focal', N*2))
#'
#' ## -------------
#' # groups completely equal
#' dat1 <- simdata(a, d, N, itemtype = 'dich')
#' dat2 <- simdata(a, d, N*2, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#'
#' #DIF (all other items as anchors)
#' SIBTEST(dat, group, suspect_set = 6)
#'
#' #DIF (specific anchors)
#' SIBTEST(dat, group, match_set = 1:5, suspect_set = 6)
#'
#' # DBF (all and specific anchors, respectively)
#' SIBTEST(dat, group, suspect_set = 11:30)
#' SIBTEST(dat, group, match_set = 1:5, suspect_set = 11:30)
#'
#' #DTF
#' SIBTEST(dat, group, suspect_set = 11:30)
#' SIBTEST(dat, group, match_set = 1:10) #equivalent
#'
#' # different hyper pars
#' dat1 <- simdata(a, d, N, itemtype = 'dich')
#' dat2 <- simdata(a, d, N*2, itemtype = 'dich', mu = .5, sigma = matrix(1.5))
#' dat <- rbind(dat1, dat2)
#' SIBTEST(dat, group, 6:30)
#' SIBTEST(dat, group, 11:30)
#'
#' #DIF testing with anchors 1 through 5
#' SIBTEST(dat, group, 6, match_set = 1:5)
#' SIBTEST(dat, group, 7, match_set = 1:5)
#' SIBTEST(dat, group, 8, match_set = 1:5)
#'
#' #DIF testing with all other items as anchors
#' SIBTEST(dat, group, 6)
#' SIBTEST(dat, group, 7)
#' SIBTEST(dat, group, 8)
#'
#' ## -------------
#' ## systematic differing slopes and intercepts (clear DTF)
#' dat1 <- simdata(a, d, N, itemtype = 'dich')
#' dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)), d + c(numeric(15), rnorm(n-15, 1, 1)),
#'   N*2, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' SIBTEST(dat, group, 6:30)
#' SIBTEST(dat, group, 11:30)
#'
#' #DIF testing using valid anchors
#' SIBTEST(dat, group, suspect_set = 6, match_set = 1:5)
#' SIBTEST(dat, group, suspect_set = 7, match_set = 1:5)
#' SIBTEST(dat, group, suspect_set = 30, match_set = 1:5)
#'
#' }
SIBTEST <- function(dat, group, suspect_set, match_set, focal_name,
                    guess_correction = 0, Jmin = 2,
                    pk_focal = FALSE, correction = TRUE, details = FALSE){

    CA <- function(dat, guess_correction = rep(0, ncol(dat))){
        n <- ncol(dat)
        C <- cov(dat)
        d <- diag(C)
        dich <- apply(dat, 2, function(x) length(unique(x)) == 2L)
        for(i in seq_len(ncol(dat))){
            if(dich[i] & guess_correction[i] > 0){
                U <- mean(dat[,i]) - min(dat[,i])
                v <- max((U - guess_correction[i]) / (1 - guess_correction[i]), 0)
                d[i] <- v * (1 - v)
            }
        }
        (n/(n - 1)) * (1 - sum(d)/sum(C))
    }
    find_intersection <- function(diff, weight, use, scores){
        k <- scores[use]
        diff <- diff[use]
        weight <- weight[use]
        mod <- lm(diff ~ k, weights = weight)
        cfs <- coef(mod)
        ks <- -cfs[1L]/cfs[2L]
        scores > signif(ks, 1L)
    }

    if(any(is.na(dat)))
        stop('SIBTEST does not support datasets with missing values.')
    if(missing(focal_name))
        focal_name <- unique(group)[2L]
    stopifnot(focal_name %in% group)
    stopifnot(nrow(dat) == length(group))
    group <- ifelse(group == focal_name, 'reference', 'focal')
    stopifnot(!(missing(suspect_set) && missing(match_set)))
    index <- 1L:ncol(dat)
    if(missing(match_set)) match_set <- index[-suspect_set]
    else if(missing(suspect_set)) suspect_set <- index[-match_set]
    if(length(guess_correction) > 1L){
        stopifnot(length(guess_correction) == ncol(dat))
    } else guess_correction <- rep(guess_correction, ncol(dat))
    if(missing(suspect_set)){
        suspect_set <- index[!(index %in% match_set)]
    } else suspect_set <- suspect_set
    stopifnot(!any(match_set %in% suspect_set))
    N <- nrow(dat)
    focal_dat <- dat[group == 'focal',]
    focal_match_scores <- rowSums(focal_dat[,match_set, drop=FALSE])
    focal_suspect_scores <- rowSums(focal_dat[,suspect_set, drop=FALSE])
    ref_dat <- dat[group == 'reference',]
    ref_match_scores <- rowSums(ref_dat[,match_set, drop=FALSE])
    ref_suspect_scores <- rowSums(ref_dat[,suspect_set, drop=FALSE])
    tab_scores <- table(rowSums(dat[ ,match_set, drop=FALSE]))
    pkstar <- tab_scores / sum(tab_scores)
    scores <- as.integer(names(pkstar))

    # selection
    tab_focal <- tab_ref <- numeric(length(tab_scores))
    II <- tab_scores > Jmin
    tab1 <- table(focal_match_scores)
    match <- match(names(II), names(tab1), nomatch=0)
    II[match] <- tab1 > Jmin
    tab_focal[match] <- tab1
    tab2 <- table(ref_match_scores)
    match <- match(names(II), names(tab2), nomatch=0)
    II[match] <- II[match] & tab2 > Jmin
    tab_ref[match] <- tab2
    II <- II & scores != min(scores) & scores != max(scores)

    n <- length(match_set)
    Xbar_ref <- mean(ref_match_scores)
    Xbar_focal <- mean(focal_match_scores)
    b_focal <- CA(focal_dat[,match_set, drop=FALSE], guess_correction[match_set])
    b_ref <- CA(ref_dat[,match_set, drop=FALSE], guess_correction[match_set])
    V_ref <- V_focal <- Ybar_ref <- Ybar_focal <- sigma_ref <- sigma_focal <-
        numeric(length(tab_scores))
    for(kk in seq_len(length(tab_scores))){
        k <- scores[kk]
        pick_ref <- k == ref_match_scores
        pick_focal <- k == focal_match_scores
        V_ref[kk] <- 1/n * (Xbar_ref + b_ref * (k - Xbar_ref))
        V_focal[kk] <- 1/n * (Xbar_focal + b_focal * (k - Xbar_focal))
        sigma_focal[kk] <- var(focal_suspect_scores[pick_focal])
        sigma_ref[kk] <- var(ref_suspect_scores[pick_ref])
        Ybar_ref[kk] <- mean(ref_suspect_scores[pick_ref])
        Ybar_focal[kk] <- mean(focal_suspect_scores[pick_focal])
    }
    V <- (V_ref + V_focal) / 2
    sigma_focal <- ifelse(is.na(sigma_focal), 0, sigma_focal)
    sigma_ref <- ifelse(is.na(sigma_ref), 0, sigma_ref)
    Ybar_ref <- ifelse(is.nan(Ybar_ref), 0, Ybar_ref)
    Ybar_focal <- ifelse(is.nan(Ybar_focal), 0, Ybar_focal)
    II <- II & sigma_ref != 0 & sigma_focal != 0
    II[scores < mean(guess_correction*ncol(dat))] <- FALSE
    if(pk_focal){
        tmp <- table(rowSums(focal_dat[ ,match_set, drop=FALSE]))
        tmp <- tmp / sum(tmp)
        match <- match(names(II), names(tmp), nomatch=0)
        pkstar[] <- 0
        pkstar[match] <- tmp
    }
    pkstar[!II] <- 0
    pkstar <- pkstar / sum(pkstar)
    tab_focal[tab_focal == 0] <- NA
    tab_ref[tab_ref == 0] <- NA
    sigma_uni <- sqrt(sum(pkstar^2 * (sigma_focal/tab_focal + sigma_ref/tab_ref), na.rm = TRUE))
    ystar_ref_vec <- ystar_focal_vec <- numeric(length(II))
    for(kk in seq_len(length(tab_scores))){
        if(!II[kk]) next
        if(correction){
            M_ref <- (Ybar_ref[kk+1] - Ybar_ref[kk-1]) / (V_ref[kk+1] - V_ref[kk-1])
            M_focal <- (Ybar_focal[kk+1] - Ybar_focal[kk-1]) / (V_focal[kk+1] - V_focal[kk-1])
            ystar_ref_vec[kk] <- Ybar_ref[kk] + M_ref * (V[kk] - V_ref[kk])
            ystar_focal_vec[kk] <- Ybar_focal[kk] + M_focal * (V[kk] - V_focal[kk])
        } else {
            ystar_ref_vec[kk] <- Ybar_ref[kk]
            ystar_focal_vec[kk] <- Ybar_focal[kk]
        }
    }
    crossvec <- find_intersection(ystar_ref_vec - ystar_focal_vec, pmax(tab_ref, tab_focal),
                                  use = pmax(tab_ref, tab_focal)/N > .01, scores=scores)

    # compute stats
    beta_uni <- beta_cross <- 0
    for(kk in seq_len(length(tab_scores))){
        if(!II[kk]) next
        beta_uni <- beta_uni + pkstar[kk] * (ystar_focal_vec[kk] - ystar_ref_vec[kk])
        if(!crossvec[kk]) beta_cross <- beta_cross + pkstar[kk] * (ystar_ref_vec[kk] - ystar_focal_vec[kk])
        else beta_cross <- beta_cross + pkstar[kk] * (ystar_focal_vec[kk] - ystar_ref_vec[kk])
    }
    beta_cross <- abs(beta_cross)
    X2_uni <- (beta_uni / sigma_uni)^2
    p_uni <- (1 - pchisq(X2_uni, 1L))
    beta1 <- sum(pkstar[crossvec] * (ystar_focal_vec[crossvec] - ystar_ref_vec[crossvec]))
    beta2 <- sum(pkstar[!crossvec] * (ystar_focal_vec[!crossvec] - ystar_ref_vec[!crossvec]))
    sigma1 <- sqrt(sum(pkstar[crossvec]^2 * (sigma_focal[crossvec]/tab_focal[crossvec] +
                                                 sigma_ref[crossvec]/tab_ref[crossvec]), na.rm = TRUE))
    sigma2 <- sqrt(sum(pkstar[!crossvec]^2 * (sigma_focal[!crossvec]/tab_focal[!crossvec] +
                                                 sigma_ref[!crossvec]/tab_ref[!crossvec]), na.rm = TRUE))
    df <- 0L
    if(sigma1 > 0) df <- df + 1L else sigma1 <- NA
    if(sigma2 > 0) df <- df + 1L else sigma2 <- NA
    X2_cross <- sum((beta1/sigma1)^2, (beta2/sigma2)^2, na.rm = TRUE)
    p_cross <- pchisq(X2_cross, df, lower.tail = FALSE)
    ret <- data.frame(focal_group=focal_name, n_matched_set=length(match_set),
                      n_suspect_set = length(suspect_set),
                      beta = c(beta_uni, beta_cross), SE=c(sigma_uni, NA),
                      X2=c(X2_uni, X2_cross),
                      df=c(1, df), p = c(p_uni, p_cross))
    rownames(ret) <- c('SIBTEST', 'CSIBTEST')
    class(ret) <- c('mirt_df', 'data.frame')
    if(details){
        ret <- data.frame(pkstar=unname(as.numeric(pkstar)),
                          sigma_focal=sigma_focal, sigma_ref=sigma_ref,
                          Y_focal=Ybar_focal, Y_ref=Ybar_ref,
                          Ystar_focal=ystar_focal_vec, Ystar_ref=ystar_ref_vec,
                          row.names = names(pkstar))
    }
    ret
}
