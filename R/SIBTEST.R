#' Simultaneous Item Bias Test (SIBTEST)
#'
#' Classical test theory approach to detecting DIF for unidimensional tests
#' by applying a regression-corrected matched-total score approach.
#' SIBTEST is similar to the Mantel-Haenszel approach for detecting DIF but uses a regression
#' correction based on the KR-20/coefficient alpha reliability index to correct the observed
#' differences when the latent trait distributions are not equal.
#' Function supports the standard SIBTEST for dichotomous and poltomous data (compensatory) and
#' also supports crossed DIF testing (i.e., non-compensatory).
#'
#' @param dat integer dataset to be tested containing dichotomous or polytomous responses
#' @param group a vector indicating group membership
#' @param match_set an integer vector indicating which items to use as the items which are matched
#'   (i.e., contain no DIF). These are analogous to 'achor' items in the likelihood method to locate
#'   DIF. If missing, all items other than the items found in the focal_set will be used
#' @param focal_name name of the focal group; e.g., 'focal'. If not specified then one will be
#'   selected automatically
#' @param focal_set an integer vector indicating which items to inspect with SIBTEST. Including only
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
#' @param cross logical; perform the crossing test for non-compensatory bias? Default is \code{FALSE}
#' @param permute number of permutations to perform when \code{cross = TRUE}. Default is 1000
#' @param correction logical; apply the composite correction for the difference between focal
#'   composite scores using the true-score regression technique? Default is \code{TRUE},
#'   reflecting Shealy and Stout's method
#' @param details logical; return a data.frame containing the details required to compute SIBTEST?
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords SIBTEST, crossed-SIBTEST
#' @aliases SIBTEST
#' @export SIBTEST
#'
#' @references
#'
#' Chang, H. H., Mazzeo, J. & Roussos, L. (1996). DIF for Polytomously Scored Items: An Adaptation of the
#'   SIBTEST Procedure. Journal of Educational Measurement, 33, 333-353.
#'
#' Li, H.-H. & Stout, W. (1996). A new procedure for detetion of crossing DIF. Psychometrika, 61, 647-677.
#'
#' Shealy, R. & Stout, W. (1993). A model-based standardization approach that separates true
#'   bias/DIF from group ability differences and ddetect test bias/DTF as well as item bias/DIF.
#'   Psychometrika, 58, 159-194.
#'
#' @examples
#' \dontrun{
#'
#' library(mirt)
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
#' dat1 <- simdata(a, d, N, itemtype = '2PL')
#' dat2 <- simdata(a, d, N*2, itemtype = '2PL')
#' dat <- rbind(dat1, dat2)
#'
#' #DIF (all other items as anchors)
#' SIBTEST(dat, group, focal_set = 6)
#'
#' #DIF (specific anchors)
#' SIBTEST(dat, group, match_set = 1:5, focal_set = 6)
#'
#' # DBF (all and specific anchors, respectively)
#' SIBTEST(dat, group, focal_set = 11:30)
#' SIBTEST(dat, group, match_set = 1:5, focal_set = 11:30)
#'
#' #DTF
#' SIBTEST(dat, group, focal_set = 11:30)
#' SIBTEST(dat, group, match_set = 1:10) #equivalent
#'
#' # different hyper pars
#' dat1 <- simdata(a, d, N, itemtype = '2PL')
#' dat2 <- simdata(a, d, N*2, itemtype = '2PL', mu = .5, sigma = matrix(1.5))
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
#' #crossed SIBTEST
#' SIBTEST(dat, group, 6, match_set = 1:5, cross=TRUE)
#' SIBTEST(dat, group, 7, match_set = 1:5, cross=TRUE)
#' SIBTEST(dat, group, 8, match_set = 1:5, cross=TRUE)
#'
#' ## -------------
#' ## systematic differing slopes and intercepts (clear DTF)
#' dat1 <- simdata(a, d, N, itemtype = '2PL')
#' dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)), d + c(numeric(15), rnorm(n-15, 1, 1)),
#'   N*2, itemtype = '2PL')
#' dat <- rbind(dat1, dat2)
#' SIBTEST(dat, group, 6:30)
#' SIBTEST(dat, group, 11:30)
#'
#' #DIF testing using valid anchors
#' SIBTEST(dat, group, focal_set = 6, match_set = 1:5)
#' SIBTEST(dat, group, focal_set = 7, match_set = 1:5)
#' SIBTEST(dat, group, focal_set = 30, match_set = 1:5)
#'
#' SIBTEST(dat, group, focal_set = 11, match_set = 1:10, cross=TRUE)
#' SIBTEST(dat, group, focal_set = 30, match_set = 1:15, cross=TRUE)
#'
#' }
SIBTEST <- function(dat, group, focal_set, match_set, focal_name,
                    guess_correction = 0, Jmin = 2, cross = FALSE, permute = 1000,
                    pk_focal = FALSE, correction = TRUE, details = FALSE){

    CA <- function(dat, guess_correction = rep(0, ncol(dat))){
        n <- ncol(dat)
        C <- cov(dat)
        d <- diag(C)
        dich <- apply(dat, 2, function(x) length(unique(x)) == 2L)
        for(i in 1L:ncol(dat)){
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
        ret <- scores > signif(ks, 1L)
        pick <- which(ret)
        if(length(pick)){
            ret[min(pick)] <- NA
        } else {
            ret[length(ret)] <- NA
        }
        ret
    }

    if(any(is.na(dat)))
        stop('SIBTEST does not support datasets with missing values.')
    if(missing(focal_name))
        focal_name <- unique(group)[2L]
    stopifnot(focal_name %in% group)
    group <- ifelse(group == focal_name, 'focal', 'reference')
    stopifnot(!(missing(focal_set) && missing(match_set)))
    index <- 1L:ncol(dat)
    if(missing(match_set)) match_set <- index[-focal_set]
    else if(missing(focal_set)) focal_set <- index[-match_set]
    if(length(guess_correction) > 1L){
        stopifnot(length(guess_correction) == ncol(dat))
    } else guess_correction <- rep(guess_correction, ncol(dat))
    if(missing(focal_set)){
        suspect_set <- index[!(index %in% match_set)]
    } else suspect_set <- focal_set
    stopifnot(!any(match_set %in% focal_set))
    N <- nrow(dat)
    focal_dat <- dat[group == 'focal',]
    focal_match_scores <- rowSums(focal_dat[,match_set, drop=FALSE])
    focal_suspect_scores <- rowSums(focal_dat[,suspect_set, drop=FALSE])
    ref_dat <- dat[group == 'reference',]
    ref_match_scores <- rowSums(ref_dat[,match_set, drop=FALSE])
    ref_suspect_scores <- rowSums(ref_dat[,suspect_set, drop=FALSE])
    if(pk_focal){
        tab_scores <- table(rowSums(focal_dat[ ,match_set, drop=FALSE]))
    } else {
        tab_scores <- table(rowSums(dat[ ,match_set, drop=FALSE]))
    }
    pk <- pkstar <- tab_scores / sum(tab_scores)
    scores <- as.integer(names(pk))

    # selection
    tab_focal <- tab_ref <- numeric(length(tab_scores))
    II <- tab_scores > Jmin
    tab1 <- table(focal_match_scores)
    tab2 <- table(ref_match_scores)
    match <- match(names(II), names(tab1), nomatch=0)
    II[match] <- tab1 > Jmin
    tab_focal[match] <- tab1
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
    for(kk in 1L:length(tab_scores)){
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
    pkstar[!II] <- 0
    pkstar <- pkstar / sum(pkstar)
    sigma_uni <- sqrt(sum(pkstar^2 * (sigma_focal/tab_focal + sigma_ref/tab_ref), na.rm = TRUE))
    crossvec <- logical(length(II))
    ystar_ref_vec <- ystar_focal_vec <- numeric(length(II))
    for(kk in 1L:length(tab_scores)){
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
    if(cross)
        crossvec <- find_intersection(ystar_ref_vec - ystar_focal_vec, pmax(tab_ref, tab_focal),
                                      use = pmax(tab_ref, tab_focal)/N > .01, scores=scores)

    # compute stats
    beta_uni <- 0
    for(kk in 1L:length(tab_scores)){
        if(!II[kk] || is.na(crossvec[kk])) next
        if(!crossvec[kk]) beta_uni <- beta_uni + pkstar[kk] * (ystar_ref_vec[kk] - ystar_focal_vec[kk])
        else beta_uni <- beta_uni + pkstar[kk] * (ystar_focal_vec[kk] - ystar_ref_vec[kk])
    }
    if(cross){
        sigma_uni <- sqrt(sum((pkstar^2 * (sigma_focal/tab_focal + sigma_ref/tab_ref))[!is.na(crossvec)],
                              na.rm = TRUE))
        B <- abs(beta_uni/sigma_uni)
        B_vec <- numeric(permute)
        for(p in 1L:permute){
            diff <- sample(c(-1,1), length(ystar_ref_vec), replace = TRUE) *
                (ystar_ref_vec - ystar_focal_vec)
            crossvec <- find_intersection(diff, pmax(tab_ref, tab_focal),
                                          use = pmax(tab_ref, tab_focal)/N > .01, scores=scores)
            sigma_uni <- sqrt(sum((pkstar^2 * (sigma_focal/tab_focal + sigma_ref/tab_ref))[!is.na(crossvec)],
                                  na.rm = TRUE))
            beta <- 0
            for(kk in 1L:length(tab_scores)){
                if(!II[kk] || is.na(crossvec[kk])) next
                if(!crossvec[kk]) beta <- beta + pkstar[kk] * (diff[kk])
                else beta <- beta + pkstar[kk] * (diff[kk])
            }
            B_vec[p] <- beta/sigma_uni
        }
        z <- B * sign(beta_uni)
        p <- mean(abs(B_vec) >= B)
    } else {
        z <- beta_uni / sigma_uni
        p <- (1 - pnorm(abs(z))) * 2
    }
    ret <- data.frame(focal_group=focal_name, n_matched_set=length(match_set),
                      n_focal_set = length(focal_set),
                      beta = beta_uni, z, p = p)
    name <- ifelse(cross, 'Crossed_SIBTEST', 'SIBTEST')
    rownames(ret) <- name
    if(details){
        ret <- data.frame(pkstar=unname(as.numeric(pkstar)),
                          sigma_focal=sigma_focal, sigma_ref=sigma_ref,
                          Y_focal=Ybar_focal, Y_ref=Ybar_ref,
                          Ystar_focal=ystar_focal_vec, Ystar_ref=ystar_ref_vec,
                          row.names = names(pkstar))
    }
    ret
}
