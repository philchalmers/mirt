#' (Generalized) Simultaneous Item Bias Test (SIBTEST)
#'
#' Classical test theory approach to detecting unidirectional and bidirectional (with one
#' crossing location) DIF. This family of statistics is intended for unidimensional tests,
#' and applies a regression-corrected matched-total score approach to quantify the response
#' bias between two or more groups. Can be used for DIF, DBF, and DTF testing with two or more
#' discrete groups.
#'
#' SIBTEST is similar to the Mantel-Haenszel approach for detecting DIF but uses a regression
#' correction based on the KR-20/coefficient alpha reliability index to correct the observed
#' differences when the latent trait distributions are not equal.
#' Function supports the standard SIBTEST for dichotomous and polytomous data (compensatory) and
#' supports crossing DIF testing (i.e., non-compensatory/non-uniform) using the asymptotic sampling
#' distribution version of the Crossing-SIBTEST (CSIBTEST) statistic described by
#' Chalmers (2018) and the permutation method described by Li and Stout (1996). This
#' function also supports the multi-group generalizations (GSIBTEST and GCSIBTEST)
#' proposed by Chalmers and Zheng (in press), where users may specify alternative
#' contrast matrices to evaluate specific comparisons between groups as well as
#' perform joint hypothesis tests.
#'
#' @param dat integer-based dataset to be tested, containing dichotomous or polytomous responses
#' @param group a vector indicating group membership with the same length as the number of rows in
#'   \code{dat}
#' @param match_set an integer vector indicating which items to use as the items which are matched
#'   (i.e., contain no DIF). These are analogous to 'anchor' items in the likelihood method to locate
#'   DIF. If missing, all items other than the items found in the \code{suspect_set} will be used
#' @param focal_name name of the focal group; e.g., \code{'focal'}. If not specified then one will be
#'   selected automatically using \code{unique(group)[2]}
#' @param suspect_set an integer vector indicating which items to inspect with SIBTEST. Including only
#'   one value will perform a DIF test, while including more than one will perform a simultaneous
#'   bundle test (DBF); including all non-matched items will perform DTF.
#'   If missing, a simultaneous test using all the items not listed in match_set
#'   will be used (i.e., DTF)
#' @param guess_correction a vector of numbers from 0 to 1 indicating how much to correct the items
#'   for guessing. It's length should be the same as ncol(dat)
#' @param na.rm logical; remove rows in \code{dat} with any missing values? If \code{TRUE},
#'   rows with missing data will be removed, as well as the corresponding elements in the \code{group}
#'   input
#' @param randomize logical; perform the crossing test for non-compensatory bias
#'   using Li and Stout's (1996) permutation approach? Default is \code{FALSE}, which uses the
#'   ad-hoc mixed degrees of freedom method suggested by Chalmers (2018)
#' @param permute number of permutations to perform when \code{randomize = TRUE}. Default is 1000
#' @param p.adjust.method a character input dictating which \code{method} to use in \code{\link{p.adjust}}.
#'   when studying more than two groups. Default does not present any p-value adjustments
#' @param C a contrast matrix to use for pooled testing with more than two groups. Default uses an
#'   effects coding approach, where the last group (last column of the matrix) is treated as the reference
#'   group, and each column is associated with the respective name via \code{unique(group)} (i.e., the first
#'   column is the coefficient for \code{unique(group)[1]}, second column for \code{unique(group)[2]}, and
#'   so on)
#' @param Jmin the minimum number of observations required when splitting the data into focal and
#'   reference groups conditioned on the matched set
#' @param pk_focal logical; using the group weights from the focal group instead of the total
#'   sample? Default is FALSE as per Shealy and Stout's recommendation
#' @param correction logical; apply the composite correction for the difference between focal
#'   composite scores using the true-score regression technique? Default is \code{TRUE},
#'   reflecting Shealy and Stout's linear extrapolation method
#' @param remove_cross logical; remove the subtest information associated with the approximate
#'   crossing location? If TRUE this reflects the CSIBTEST definition of Li and Stout (1996);
#'   if FALSE, this reflects the version of CSIBTEST utilized by Chalmers (2018). Only applicable
#'   in two-group settings (in multi-group this is fixed to FALSE)
#' @param details logical; return a data.frame containing the details required to compute SIBTEST?
#' @param plot a character input indicating the type of plot to construct. Options are \code{'none'}
#'   (default), \code{'observed'} for the scaled focal subtest scores against the matched subtest
#'   scores, \code{'weights'} for the proportion weights used (i.e., the proportion of observations at
#'   each matched score), \code{'difference'} for the difference between the scaled focal subtest scores
#'   against the matched subtest scores, and \code{'wdifference'} for the conditional differences multiplied
#'   by each respective weight. Note that the last plot reflects the components used in SIBTEST, and therefore
#'   the sum of these plotted observations will equal the beta coefficient for SIBTEST
#' @param ... additional plotting arguments to be passed
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords SIBTEST, crossing SIBTEST, generalized SIBTEST
#' @aliases SIBTEST
#' @export SIBTEST
#'
#' @references
#'
#' Chalmers, R. P. (2018). Improving the Crossing-SIBTEST statistic for
#' detecting non-uniform DIF. \emph{Psychometrika, 83}, 2, 376-386.
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Chalmers, R. P. and Zheng, G. (accepted). Multi-group Generalizations of
#' SIBTEST and Crossing-SIBTEST. \emph{Applied Measurement in Education}.
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
#'
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
#' # DIF (all other items as anchors)
#' SIBTEST(dat, group, suspect_set = 6)
#'
#' # Some plots depicting the above tests
#' SIBTEST(dat, group, suspect_set = 6, plot = 'observed')
#' SIBTEST(dat, group, suspect_set = 6, plot = 'weights')
#' SIBTEST(dat, group, suspect_set = 6, plot = 'wdifference')
#'
#' # Include CSIBTEST with randomization method
#' SIBTEST(dat, group, suspect_set = 6, randomize = TRUE)
#'
#' # remove crossing-location (identical to Li and Stout 1996 definition of CSIBTEST)
#' SIBTEST(dat, group, suspect_set = 6, randomize = TRUE, remove_cross=TRUE)
#'
#' # DIF (specific anchors)
#' SIBTEST(dat, group, match_set = 1:5, suspect_set = 6)
#' SIBTEST(dat, group, match_set = 1:5, suspect_set = 6, randomize=TRUE)
#'
#' # DBF (all and specific anchors, respectively)
#' SIBTEST(dat, group, suspect_set = 11:30)
#' SIBTEST(dat, group, match_set = 1:5, suspect_set = 11:30)
#'
#' # DTF
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
#' # DIF testing with anchors 1 through 5
#' SIBTEST(dat, group, 6, match_set = 1:5)
#' SIBTEST(dat, group, 7, match_set = 1:5)
#' SIBTEST(dat, group, 8, match_set = 1:5)
#'
#' # DIF testing with all other items as anchors
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
#' # Some plots depicting the above tests
#' SIBTEST(dat, group, suspect_set = 11:30, plot = 'observed')
#' SIBTEST(dat, group, suspect_set = 11:30, plot = 'weights')
#' SIBTEST(dat, group, suspect_set = 11:30, plot = 'wdifference')
#'
#' # DIF testing using valid anchors
#' SIBTEST(dat, group, suspect_set = 6, match_set = 1:5)
#' SIBTEST(dat, group, suspect_set = 7, match_set = 1:5)
#' SIBTEST(dat, group, suspect_set = 30, match_set = 1:5)
#'
#' # randomization method is fairly poor when smaller matched-set used
#' SIBTEST(dat, group, suspect_set = 30, match_set = 1:5, randomize=TRUE)
#' SIBTEST(dat, group, suspect_set = 30, randomize=TRUE)
#'
#' ## ----------------------------------
#' # three group SIBTEST test
#' set.seed(1234)
#' n <- 30
#' N <- 1000
#' a <- matrix(1, n)
#' d <- matrix(rnorm(n), n)
#' group <- c(rep('group1', N), rep('group2', N), rep('group3', N))
#'
#' # groups completely equal
#' dat1 <- simdata(a, d, N, itemtype = 'dich')
#' dat2 <- simdata(a, d, N, itemtype = 'dich')
#' dat3 <- simdata(a, d, N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2, dat3)
#'
#' # omnibus test using effects-coding contrast matrix (default)
#' SIBTEST(dat, group, suspect_set = 6)
#' SIBTEST(dat, group, suspect_set = 6, randomize=TRUE)
#'
#' # explicit contrasts
#' SIBTEST(dat, group, suspect_set = 6, randomize=TRUE,
#'         C = matrix(c(1,-1,0), 1))
#'
#' ## systematic differing slopes and intercepts
#' dat2 <- simdata(a + c(numeric(15), .5,.5,.5,.5,.5, numeric(10)),
#'         d + c(numeric(15), 0,.6,.7,.8,.9, numeric(10)),
#'         N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2, dat3)
#'
#' SIBTEST(dat, group, suspect_set = 16)
#' SIBTEST(dat, group, suspect_set = 16, randomize=TRUE)
#'
#' SIBTEST(dat, group, suspect_set = 19)
#' SIBTEST(dat, group, suspect_set = 19, randomize=TRUE)
#'
#'
#' }
SIBTEST <- function(dat, group, suspect_set, match_set, focal_name = unique(group)[2],
                    guess_correction = 0, Jmin = 5, na.rm = FALSE, randomize = FALSE,
                    C = cbind(1, -diag(length(unique(group)) - 1L)),
                    p.adjust.method = 'none', permute = 1000, pk_focal = FALSE,
                    correction = TRUE, remove_cross = FALSE, details = FALSE, plot = 'none', ...){

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
    find_intersection <- function(diff, weight, use, scores, remove_cross,
                                  tab_match = NULL, C = NULL){
        if(is.matrix(diff)){
            Ystar <- diff
            ks <- Ns <- numeric(ncol(diff)*(ncol(diff)-1)/2)
            ind <- 1L
            for(i in 1L:ncol(Ystar)){
                for(j in 1L:ncol(Ystar)){
                    if(j > i){
                        weight <- pmax(tab_match[[i]], tab_match[[j]])
                        N <- sum(tab_match[[i]], na.rm=TRUE) +
                            sum(tab_match[[j]], na.rm=TRUE)
                        use <- weight/N > .01
                        use[is.na(use)] <- FALSE
                        k <- scores[use]
                        diff <- Ystar[use, i] - Ystar[use, j]
                        weight <- weight[use]
                        mod <- lm(diff ~ k, weights = weight)
                        cfs <- coef(mod)
                        ks[ind] <- -cfs[1L]/cfs[2L]
                        Ns[ind] <- N
                        ind <- ind + 1L
                    }
                }
            }
            # ks <- sum(ks * (Ns / sum(Ns)))
            ksmat <- Nsmat <- matrix(NA, ncol(C), ncol(C))
            ksmat[lower.tri(ksmat)] <- ks
            Nsmat[lower.tri(ksmat)] <- Ns
            ret <- matrix(FALSE, length(scores), nrow(C))

            # computed weighted ks per contrast
            for(j in 1:nrow(C)){
                whcC <- which(C[j,] != 0)
                wks <- sum(ksmat[whcC, whcC] * Nsmat[whcC, whcC] /
                               sum(Nsmat[whcC, whcC], na.rm = TRUE), na.rm=TRUE)
                ret[,j] <- scores > signif(wks, 1L)
                ## leave crossing location in for easier math
            }
        } else {
            k <- scores[use]
            diff <- diff[use]
            weight <- weight[use]
            mod <- lm(diff ~ k, weights = weight)
            cfs <- coef(mod)
            ks <- -cfs[1L]/cfs[2L]
            ret <- scores > signif(ks, 1L)
            if(remove_cross){
                pick <- which(ret)
                if(length(pick)){
                    ret[min(pick)] <- NA
                } else {
                    ret[length(ret)] <- NA
                }
            }
        }
        ret
    }

    if(na.rm){
        is_na <- is.na(rowSums(dat))
        dat <- dat[!is_na, , drop=FALSE]
        group <- group[!is_na]
    }
    if(any(is.na(dat)))
        stop(c('SIBTEST does not support datasets with missing values. \n  ',
               'Pass na.rm=TRUE to remove missing rows'), call.=FALSE)
    ngroups <- length(unique(group))
    N <- nrow(dat)
    stopifnot(ngroups >= 2L)
    stopifnot(focal_name %in% group)
    stopifnot(N == length(group))
    stopifnot(Jmin >= 2)
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

    if(ngroups == 2L){
        focal_dat <- dat[group == focal_name,]
        focal_match_scores <- rowSums(focal_dat[,match_set, drop=FALSE])
        focal_suspect_scores <- rowSums(focal_dat[,suspect_set, drop=FALSE])
        ref_dat <- dat[group != focal_name,]
        ref_match_scores <- rowSums(ref_dat[,match_set, drop=FALSE])
        ref_suspect_scores <- rowSums(ref_dat[,suspect_set, drop=FALSE])
        tab_scores <- table(rowSums(dat[ ,match_set, drop=FALSE]))
        pkstar <- tab_scores / sum(tab_scores)
        scores <- as.integer(names(pkstar))

        # selection
        tab_focal <- tab_ref <- numeric(length(tab_scores))
        II <- tab_scores > Jmin
        tab1 <- table(focal_match_scores)
        match <- names(II) %in% names(tab1)
        II[match] <- tab1 > Jmin
        tab_focal[match] <- tab1
        tab2 <- table(ref_match_scores)
        match <- names(II) %in% names(tab2)
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
            match <- names(II) %in% names(tmp)
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
                                      use = pmax(tab_ref, tab_focal)/N > .01, scores=scores,
                                      remove_cross=remove_cross)

        # compute stats
        beta_uni <- beta_cross <- 0
        for(kk in seq_len(length(tab_scores))){
            if(!II[kk]) next
            beta_uni <- beta_uni + pkstar[kk] * (ystar_ref_vec[kk]- ystar_focal_vec[kk])
            if(is.na(crossvec[kk])) next
            if(!crossvec[kk])
                beta_cross <- beta_cross + pkstar[kk] * (ystar_ref_vec[kk] - ystar_focal_vec[kk])
            else
                beta_cross <- beta_cross + pkstar[kk] * (ystar_focal_vec[kk] - ystar_ref_vec[kk])
        }
        X2_uni <- (beta_uni / sigma_uni)^2
        p_uni <- (1 - pchisq(X2_uni, 1L))
        beta1 <- sum(pkstar[crossvec] * (ystar_focal_vec[crossvec] - ystar_ref_vec[crossvec]))
        beta2 <- sum(pkstar[!crossvec] * (ystar_focal_vec[!crossvec] - ystar_ref_vec[!crossvec]))
        sigma1 <- sqrt(sum(pkstar[crossvec]^2 * (sigma_focal[crossvec]/tab_focal[crossvec] +
                                                     sigma_ref[crossvec]/tab_ref[crossvec]),
                           na.rm = TRUE))
        sigma2 <- sqrt(sum(pkstar[!crossvec]^2 * (sigma_focal[!crossvec]/tab_focal[!crossvec] +
                                                      sigma_ref[!crossvec]/tab_ref[!crossvec]),
                           na.rm = TRUE))
        df <- 0L
        if(sigma1 > 0) df <- df + 1L else sigma1 <- NA
        if(sigma2 > 0) df <- df + 1L else sigma2 <- NA
        X2_cross <- sum((beta1/sigma1)^2, (beta2/sigma2)^2, na.rm = TRUE)
        p_cross <- pchisq(X2_cross, df, lower.tail = FALSE)
        B_vec <- numeric(permute)
        sigma_cross <- NA
        ret <- data.frame(focal_group=focal_name, n_matched_set=length(match_set),
                          n_suspect_set = length(suspect_set),
                          beta = c(beta_uni, beta_cross), SE=c(sigma_uni, NA),
                          X2=c(X2_uni, X2_cross),
                          df=c(1, df), p = c(p_uni, p_cross))
        rownames(ret) <- c('SIBTEST', 'CSIBTEST')
        if(randomize){
            crossvec <- find_intersection(ystar_ref_vec - ystar_focal_vec,
                                          pmax(tab_ref, tab_focal),
                                          use = pmax(tab_ref, tab_focal)/N > .01,
                                          scores=scores,
                                          remove_cross=remove_cross)
            beta_cross2 <- 0
            for(kk in seq_len(length(tab_scores))){
                if(!II[kk] || is.na(crossvec[kk])) next
                if(!crossvec[kk])
                    beta_cross2 <- beta_cross2 +
                        pkstar[kk] * (ystar_ref_vec[kk] - ystar_focal_vec[kk])
                else
                    beta_cross2 <- beta_cross2 +
                        pkstar[kk] * (ystar_focal_vec[kk] - ystar_ref_vec[kk])
            }
            sigma_cross <- sqrt(sum((pkstar^2 * (sigma_focal/tab_focal +
                                                     sigma_ref/tab_ref))[!is.na(crossvec)],
                                    na.rm = TRUE))
            B <- abs(beta_cross2/sigma_cross)
            for(p in 1L:permute){
                diff <- sample(c(-1,1), length(ystar_ref_vec), replace = TRUE) *
                    (ystar_ref_vec - ystar_focal_vec)
                crossvec <- find_intersection(diff, pmax(tab_ref, tab_focal),
                                              use = pmax(tab_ref, tab_focal)/N > .01,
                                              scores=scores,
                                              remove_cross=remove_cross)
                beta <- 0
                for(kk in 1L:length(tab_scores)){
                    if(!II[kk] || is.na(crossvec[kk])) next
                    if(!crossvec[kk]) beta <- beta + pkstar[kk] * (diff[kk])
                    else beta <- beta + pkstar[kk] * (-diff[kk])
                }
                B_vec[p] <- beta/sigma_cross
            }
            p_cross2 <- mean(abs(B_vec) >= B)
            ret$X2[2] <- ret$df[2] <- NA
            ret$SE[2] <- sigma_cross
            ret$p[2] <- p_cross2
        }
        class(ret) <- c('mirt_df', 'data.frame')
        if(details){
            ret <- data.frame(pkstar=unname(as.numeric(pkstar)),
                              sigma_focal=sigma_focal, sigma_ref=sigma_ref,
                              Y_focal=Ybar_focal, Y_ref=Ybar_ref,
                              Ystar_focal=ystar_focal_vec,
                              Ystar_ref=ystar_ref_vec,
                              row.names = names(pkstar))
            if(randomize) attr(ret, "B_vec") <- B_vec
        }
        if(plot != 'none'){
            ret <- data.frame(total_score = 1:length(pkstar) - 1,
                              pkstar=unname(as.numeric(pkstar)),
                              Ystar=c(ystar_focal_vec, ystar_ref_vec),
                              Ystar_diff =ystar_focal_vec - ystar_ref_vec,
                              Ystar_diff_pkstar =
                                  unname((ystar_focal_vec-ystar_ref_vec)*as.numeric(pkstar)),
                              group = rep(c('focal', 'reference'), each=length(pkstar)))
            if(plot != "freq")
                for(i in 1L:nrow(ret))
                    if(ret$pkstar[i] == 0) ret[i, ] <- NA
            if(plot == "observed")
                return(lattice::xyplot(Ystar ~ total_score, data=ret,
                                       groups = group, xlab='Matched subtest',
                                       ylab = 'Scaled focal subtest', type='b',
                                       auto.key=list(space='right'), ...))
            if(plot == "weights")
                return(lattice::xyplot(pkstar~total_score,
                                       data=subset(ret, group=='focal'),
                                       xlab='Matched subtest',
                                       ylab = 'Proportion', type = 'b', ...))
            if(plot == "difference")
                return(lattice::xyplot(Ystar_diff~total_score,
                                       data=subset(ret, group=='focal'),
                                       xlab='Matched subtest',
                                       ylab = 'Focal subtest difference',
                                       panel = function(x, y, ...){
                                           panel.xyplot(x, y, type = 'b', ...)
                                           panel.abline(h=0, lty=2, col='red')
                                       }, ...))
            if(plot == "wdifference")
                return(lattice::xyplot(Ystar_diff_pkstar~total_score,
                                       data=subset(ret, group=='focal'),
                                       xlab='Matched subtest',
                                       ylab = 'Weighted focal subtest difference',
                                       panel = function(x, y, ...){
                                           panel.xyplot(x, y, type = 'b', ...)
                                           panel.abline(h=0, lty=2, col='red')
                                       }, ...))
            stop('plot argument not supported', call.=FALSE)
        }
    } else {
        # Multi-group

        if(plot != 'none')
            stop('Multi-group plots not currently supported')

        groupnms <- unique(group)
        if(ncol(C) != ngroups)
            stop(sprintf('The C matrix must have %i columns', ngroups), call.=FALSE)
        splt_dat_match <-
            lapply(groupnms, function(x) dat[group == x, match_set, drop=FALSE])
        splt_dat_suspect <-
            lapply(groupnms, function(x) dat[group == x, suspect_set, drop=FALSE])
        splt_match_scores <- lapply(splt_dat_match, rowSums)
        splt_suspect_scores <- lapply(splt_dat_suspect, rowSums)
        tab_scores <- table(rowSums(dat[ ,match_set, drop=FALSE]))
        pkstar <- tab_scores / sum(tab_scores)
        scores <- as.integer(names(pkstar))

        # selection
        splt_tab_match <- lapply(splt_match_scores, function(x){
            tab_full <- numeric(length(tab_scores))
            tab <- table(x)
            match <- names(tab_scores) %in% names(tab)
            tab_full[match] <- tab
            tab_full
        })
        II <- vector('list', nrow(C))
        for(cont in 1L:nrow(C)){
            pick <- ifelse(C[cont,] != 0, TRUE, FALSE)
            II[[cont]] <- colSums(do.call(rbind,
                                          lapply(splt_tab_match,
                                                 function(x) x > Jmin)) * pick) == sum(pick)
            II[[cont]] <- II[[cont]] & scores != min(scores) & scores != max(scores)
        }

        n <- length(match_set)
        Xbars <- sapply(splt_match_scores, mean)
        bs <- sapply(splt_dat_match, CA,
                     guess_correction=guess_correction[match_set])
        cond_stats_list <- lapply(1L:ngroups, function(ind){
            Vg <- Ybar <- sigma <- numeric(length(tab_scores))
            for(kk in seq_len(length(tab_scores))){
                k <- scores[kk]
                pick <- k == splt_match_scores[[ind]]
                Vg[kk] <- 1/n * (Xbars[ind] + bs[ind] * (k - Xbars[ind]))
                sigma[kk] <- var(splt_suspect_scores[[ind]][pick])
                Ybar[kk] <- mean(splt_suspect_scores[[ind]][pick])
            }
            sigma <- ifelse(is.na(sigma), 0, sigma)
            Ybar <- ifelse(is.nan(Ybar), 0, Ybar)
            list(Vg=Vg, sigma=sigma, Ybar=Ybar)
        })
        V <- rowMeans(do.call(cbind, lapply(cond_stats_list, function(x) x$Vg)))
        sigmas <- do.call(cbind, lapply(cond_stats_list, function(x) x$sigma))
        for(cont in 1L:nrow(C)){
            pick <- ifelse(C[cont,] != 0, TRUE, FALSE)
            pick_sigmas <- apply(sigmas, 1, function(x) all(x[pick] > 0))
            II[[cont]] <- II[[cont]] & pick_sigmas
            II[[cont]][scores < mean(guess_correction*ncol(dat))] <- FALSE
        }
        if(any(sapply(II, sum) < 3))
            stop('Too few sum scores used. Consider reducing Jmin argument
                 or adjusting guess_correction (if necessary)',
                 call.=FALSE)
        pkstar[!apply(do.call(cbind, II), 1L, any)] <- 0
        pkstar <- pkstar / sum(pkstar)
        splt_tab_match <- lapply(splt_tab_match, function(x){
            x[x == 0] <- NA
            x
        })
        Sigma_list <- lapply(1L:length(tab_scores), function(i){
            vals <- numeric(ngroups)
            for(g in 1L:ngroups)
                vals[g] <- cond_stats_list[[g]]$sigma[i] / splt_tab_match[[g]][i]
            C %*% diag(vals) %*% t(C) * pkstar[i]^2
        })
        Ystar <- sapply(1L:ngroups, function(ind){
            ystar_vec <- numeric(length(II[[1L]]))
            ystar_vec[] <- NA
            for(kk in seq_len(length(tab_scores))){
                if(!any(sapply(II, function(x) x[kk]))) next
                if(correction){
                    M <- with(cond_stats_list[[ind]], (Ybar[kk+1] - Ybar[kk-1]) /
                                  (Vg[kk+1] - Vg[kk-1]))
                    ystar_vec[kk] <- with(cond_stats_list[[ind]], Ybar[kk] +
                                              M * (V[kk] - Vg[kk]))
                } else {
                    ystar_vec[kk] <- with(cond_stats_list[[ind]], Ybar[kk])
                }
            }
            ystar_vec
        })

        crossmat <- find_intersection(Ystar, scores=scores, remove_cross=FALSE,
                                      tab_match=splt_tab_match, C=C)
        # compute stats
        Ystar[is.na(Ystar)] <- 0
        beta_uni <- beta_cross <- C %*% t(Ystar) %*% pkstar
        for(j in 1L:length(beta_cross))
            beta_cross[j] <- C[j, ,drop=FALSE] %*%(t(Ystar[!crossmat[,j], , drop=FALSE]) %*%
                                                       pkstar[!crossmat[,j]] -
                                                       t(Ystar[crossmat[,j], , drop=FALSE]) %*% pkstar[crossmat[,j]])

        sigma_uni <- matrix(0, nrow(C), nrow(C))
        for(i in seq_len(length(tab_scores)))
            if(all(!is.na(Sigma_list[[i]]))) sigma_uni <- sigma_uni + Sigma_list[[i]]
        slv_sigma_uni <- solve(sigma_uni)
        X2_uni <- t(beta_uni) %*% slv_sigma_uni %*% beta_uni
        p_uni <- (1 - pchisq(X2_uni, nrow(C)))

        crossmat_org <- find_intersection(Ystar, scores=scores,
                                          remove_cross=FALSE,
                                          tab_match=splt_tab_match, C=C)
        not_crossmat <- !crossmat_org & !is.na(crossmat_org)
        crossmat <- crossmat_org & !is.na(crossmat_org)
        beta_cross <- numeric(nrow(C))
        for(j in 1L:length(beta_cross))
            beta_cross[j] <- C[j, ,drop=FALSE] %*%
            (t(Ystar[!crossmat[,j], , drop=FALSE]) %*%
                 pkstar[!crossmat[,j]] - t(Ystar[crossmat[,j], , drop=FALSE]) %*%
                 pkstar[crossmat[,j]])
        X2_cross <- t(beta_cross) %*% slv_sigma_uni %*% beta_cross

        if(nrow(C) == 1L){
            ret <- data.frame(n_matched_set=length(match_set),
                              n_suspect_set = length(suspect_set),
                              beta = c(beta_uni, beta_cross),
                              X2=c(X2_uni, X2_cross),
                              df=c(1, NA), p = c(p_uni, NA))
            rownames(ret) <- c('GSIBTEST', 'GCSIBTEST')
            class(ret) <- c('mirt_df', 'data.frame')
        } else {
            mat1 <- rbind(as.vector(beta_uni), beta_cross)
            colnames(mat1) <- paste0('beta_', 1:nrow(C))
            ret <- data.frame(n_matched_set=length(match_set),
                              n_suspect_set = length(suspect_set),
                              mat1, X2=c(X2_uni, X2_cross),
                              df=c(nrow(C), NA), p = c(p_uni, NA))
            rownames(ret) <- c('GSIBTEST', 'GCSIBTEST')
            class(ret) <- c('mirt_df', 'data.frame')
        }
        if(randomize){
            X2_cross_vec <- numeric(permute)
            for(p in 1L:permute){
                Ystar_p <-  sample(c(-1,1), nrow(Ystar), replace = TRUE) * Ystar
                crossmat_org <- find_intersection(Ystar_p, scores=scores,
                                                  remove_cross=FALSE,
                                                  tab_match=splt_tab_match, C=C)
                not_crossmat <- !crossmat_org & !is.na(crossmat_org)
                crossmat <- crossmat_org & !is.na(crossmat_org)
                beta_cross_p <- numeric(nrow(C))
                for(j in 1L:length(beta_cross_p))
                    beta_cross_p[j] <- C[j, ,drop=FALSE] %*%
                    (t(Ystar_p[!crossmat[,j], , drop=FALSE]) %*%
                         pkstar[!crossmat[,j]] - t(Ystar_p[crossmat[,j], , drop=FALSE]) %*%
                         pkstar[crossmat[,j]])
                X2_cross_vec[p] <- t(beta_cross_p) %*% slv_sigma_uni %*% beta_cross_p
            }
            p_cross2 <- mean(abs(X2_cross_vec) >= X2_cross[1L])
            ret$p[2] <- p_cross2
        }
        if(details){
            ret <- data.frame(pkstar=unname(as.numeric(pkstar)),
                              sigma_focal=sigma_focal, sigma_ref=sigma_ref,
                              Y_focal=Ybar_focal, Y_ref=Ybar_ref,
                              Ystar_focal=ystar_focal_vec,
                              Ystar_ref=ystar_ref_vec,
                              row.names = names(pkstar))
        }
        return(ret)
    }
    ret
}
