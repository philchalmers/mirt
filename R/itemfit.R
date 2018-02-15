#' Item fit statistics
#'
#' Computes item-fit statistics for a variety of unidimensional and multidimensional models.
#' Poorly fitting items should be inspected with the empirical plots/tables
#' for unidimensional models, otherwise \code{\link{itemGAM}} can be used to diagnose
#' where the functional form of the IRT model was misspecified, or models can be refit using
#' more flexible semi-parametric response models (e.g., \code{itemtype = 'spline'}).
#'
#' @aliases itemfit
#' @param x a computed model object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{DiscreteClass}
#' @param fit_stats a character vector indicating which fit statistics should be computed.
#'   Supported inputs are:
#'
#' \itemize{
#'   \item \code{'S_X2'} : Orlando and Thissen (2000, 2003) and
#'     Kang and Chen's (2007) signed chi-squared test (default)
#'   \item \code{'Zh'} : Drasgow, Levine, & Williams (1985) Zh
#'   \item \code{'X2'} : Bock's (1972) chi-squared method.
#'     The default inputs compute Yen's (1981) Q1 variant of the X2 statistic
#'     (i.e., uses a fixed \code{group.bins = 10}). However, Bock's group-size variable
#'     median-based method can be computed by passing \code{group.fun = median} and
#'     modifying the \code{group.size} input to the desired number of bins
#'   \item \code{'G2'} : McKinley & Mills (1985) G2 statistic (similar method to Q1,
#'     but with the likelihood-ratio test).
#'   \item \code{'PV_Q1'} : Chalmers and Ng's (2017) plausible-value variant
#'     of the Q1 statistic.
#'   \item \code{'PV_Q1*'} : Chalmers and Ng's (2017) plausible-value variant
#'     of the Q1 statistic that uses parametric bootstrapping to obtain a suitable empirical
#'     distribution.
#'   \item \code{'X2*'} : Stone's (2000) fit statistics that require parametric
#'     bootstrapping
#'   \item \code{'X2*_df'} : Stone's (2000) fit statistics that require parametric
#'     bootstrapping to obtain scaled versions of the X2* and degrees of freedom
#'   \item \code{'infit'} : (Unidimensional Rasch model only) compute the
#'     infit and outfit statistics. Ignored if models are not from the Rasch family
#' }
#'
#' Note that 'infit', 'S_X2', and 'Zh' cannot be computed when there are missing response data
#' (i.e., will require multiple-imputation techniques).
#'
#' @param which.items an integer vector indicating which items to test for fit.
#'   Default tests all possible items
#' @param mincell the minimum expected cell size to be used in the S-X2 computations. Tables will be
#'   collapsed across items first if polytomous, and then across scores if necessary
#' @param mincell.X2 the minimum expected cell size to be used in the X2 computations. Tables will be
#'   collapsed if polytomous, however if this condition can not be met then the group block will
#'   be ommited in the computations
#' @param S_X2.tables logical; return the tables in a list format used to compute the S-X2 stats?
#' @param group.size approximate size of each group to be used in calculating the \eqn{\chi^2}
#'   statistic. The default \code{NA}
#'   disables this command and instead uses the \code{group.bins} input to try and construct
#'   equally sized bins
#' @param group.bins the number of bins to use for X2 and G2. For example,
#'   setting \code{group.bins = 10} will will compute Yen's (1981) Q1 statistic when \code{'X2'} is
#'   requested
#' @param group.fun function used when \code{'X2'} or \code{'G2'} are computed. Determines the central
#'   tendancy measure within each partitioned group. E.g., setting \code{group.fun = median} will
#'   obtain the median of each respective ability estimate in each subgroup (this is what was used
#'   by Bock, 1972)
#' @param empirical.plot a single numeric value or character of the item name indicating which
#'   item to plot (via \code{itemplot}) and overlay with the empirical \eqn{\theta} groupings (see
#'   \code{empirical.CI}). Useful for plotting the expected bins based on the \code{'X2'} or
#'   \code{'G2'} method
#' @param empirical.table a single numeric value or character of the item name indicating which
#'   item table of expected values should be returned. Useful for visualizing the
#'   expected bins based on the \code{'X2'} or \code{'G2'} method
#' @param empirical.CI a numeric value indicating the width of the empirical confidence interval
#'   ranging between 0 and 1 (default of 0 plots not interval). For example, a 95% confidence
#'   interval would be plotted when \code{empirical.CI = .95}. Only applicable to dichotomous items
#' @param method type of factor score estimation method. See \code{\link{fscores}} for more detail
#' @param Theta a matrix of factor scores for each person used for statistics that require
#'   empirical estimates. If supplied, arguments typically passed to \code{fscores()} will be
#'   ignored and these values will be used instead. Also required when estimating statistics
#'   with missing data via imputation
#' @param pv_draws number of plausible-value draws to obtain for PV_Q1 and PV_Q1*
#' @param boot number of parametric bootstrap samples to create for PV_Q1* and X2*
#' @param boot_dfapprox number of parametric bootstrap samples to create for the X2*_df statistic
#'   to approximate the scaling factor for X2* as well as the scaled degrees of freedom estimates
#' @param ETrange rangone of integration nodes for Stone's X2* statistic
#' @param ETpoints number of integration nodes to use for Stone's X2* statistic
#' @param impute a number indicating how many imputations to perform (passed to
#'   \code{\link{imputeMissing}}) when there are missing data present.
#'   Will return a data.frame object with the mean estimates
#'   of the stats and their imputed standard deviations
#' @param par.strip.text plotting argument passed to \code{\link{lattice}}
#' @param par.settings plotting argument passed to \code{\link{lattice}}
#' @param ... additional arguments to be passed to \code{fscores()} and \code{\link{lattice}}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords item fit
#' @export itemfit
#'
#' @seealso
#' \code{\link{personfit}}, \code{\link{itemGAM}}
#'
#' @references
#'
#' Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored
#' in two or more nominal categories. \emph{Psychometrika, 37}, 29-51.
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Chalmers, R. P. & Ng, V. (2017). Plausible-Value Imputation Statistics for Detecting
#' Item Misfit. \emph{Applied Psychological Measurement, 41}, 372-387.
#' \doi{10.1177/0146621617692079}
#'
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with
#' polychotomous item response models and standardized indices.
#' \emph{British Journal of Mathematical and Statistical Psychology, 38}, 67-86.
#'
#' Kang, T. & Chen, Troy, T. (2007). An investigation of the performance of the generalized
#' S-X2 item-fit index for polytomous IRT models. ACT
#'
#' McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit statistics.
#' Applied Psychological Measurement, 9, 49-57.
#'
#' Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for dichotomous item
#' response theory models. \emph{Applied Psychological Measurement, 24}, 50-64.
#'
#' Reise, S. P. (1990). A comparison of item- and person-fit methods of assessing model-data fit
#' in IRT. \emph{Applied Psychological Measurement, 14}, 127-137.
#'
#' Stone, C. A. (2000). Monte Carlo Based Null Distribution for an Alternative Goodness-of-Fit
#' Test Statistics in IRT Models. \emph{Journal of Educational Measurement, 37}, 58-75.
#'
#' Wright B. D. & Masters, G. N. (1982). \emph{Rating scale analysis}. MESA Press.
#'
#' Yen, W. M. (1981). Using simulation results to choose a latent trait model.
#' \emph{Applied Psychological Measurement, 5}, 245-262.
#'
#' @examples
#'
#' \dontrun{
#'
#' P <- function(Theta){exp(Theta^2 * 1.2 - 1) / (1 + exp(Theta^2 * 1.2 - 1))}
#'
#' #make some data
#' set.seed(1234)
#' a <- matrix(rlnorm(20, meanlog=0, sdlog = .1),ncol=1)
#' d <- matrix(rnorm(20),ncol=1)
#' Theta <- matrix(rnorm(2000))
#' items <- rep('2PL', 20)
#' ps <- P(Theta)
#' baditem <- numeric(2000)
#' for(i in 1:2000)
#'    baditem[i] <- sample(c(0,1), 1, prob = c(1-ps[i], ps[i]))
#' data <- cbind(simdata(a,d, 2000, items, Theta=Theta), baditem=baditem)
#'
#' x <- mirt(data, 1)
#' raschfit <- mirt(data, 1, itemtype='Rasch')
#' fit <- itemfit(x)
#' fit
#'
#' itemfit(x)
#' itemfit(x, 'X2') # just X2
#' itemfit(x, c('S_X2', 'X2')) #both S_X2 and X2
#' itemfit(x, group.bins=15, empirical.plot = 1) #empirical item plot with 15 points
#' itemfit(x, group.bins=15, empirical.plot = 21)
#'
#' # PV and X2* statistics (parametric bootstrap stats not run to save time)
#' itemfit(x, 'PV_Q1')
#'
#' # mirtCluster() # improve speed of bootstrap samples by running in parallel
#' # itemfit(x, 'PV_Q1*')
#' # itemfit(x, 'X2*') # Stone's 1993 statistic
#' # itemfit(x, 'X2*_df') # Stone's 2000 scaled statistic with df estimate
#'
#' #empirical tables
#' itemfit(x, empirical.table=1)
#' itemfit(x, empirical.table=21)
#'
#' #infit/outfit statistics. method='ML' agrees better with eRm package
#' itemfit(raschfit, 'infit', method = 'ML') #infit and outfit stats
#'
#' #same as above, but inputting ML estimates instead
#' Theta <- fscores(raschfit, method = 'ML')
#' itemfit(raschfit, 'infit', Theta=Theta)
#'
#' # fit a new more flexible model for the mis-fitting item
#' itemtype <- c(rep('2PL', 20), 'spline')
#' x2 <- mirt(data, 1, itemtype=itemtype)
#' itemfit(x2)
#' itemplot(x2, 21)
#' anova(x2, x)
#'
#' #------------------------------------------------------------
#'
#' #similar example to Kang and Chen 2007
#' a <- matrix(c(.8,.4,.7, .8, .4, .7, 1, 1, 1, 1))
#' d <- matrix(rep(c(2.0,0.0,-1,-1.5),10), ncol=4, byrow=TRUE)
#' dat <- simdata(a,d,2000, itemtype = rep('graded', 10))
#' head(dat)
#'
#' mod <- mirt(dat, 1)
#' itemfit(mod)
#' itemfit(mod, 'X2') #pretty much useless given inflated Type I error rates
#' itemfit(mod, empirical.plot = 1)
#'
#' # collapsed tables (see mincell.X2) for X2 and G2
#' itemfit(mod, empirical.table = 1)
#'
#' mod2 <- mirt(dat, 1, 'Rasch')
#' itemfit(mod2, 'infit')
#'
#' #massive list of tables
#' tables <- itemfit(mod, S_X2.tables = TRUE)
#'
#' #observed and expected total score patterns for item 1 (post collapsing)
#' tables$O[[1]]
#' tables$E[[1]]
#'
#' # fit stats with missing data (run in parallel using all cores)
#' data[sample(1:prod(dim(data)), 500)] <- NA
#' raschfit <- mirt(data, 1, itemtype='Rasch')
#'
#' mirtCluster() # run in parallel
#' itemfit(raschfit, c('S_X2', 'infit'), impute = 10)
#'
#' #alternative route: use only valid data, and create a model with the previous parameter estimates
#' data2 <- na.omit(data)
#' raschfit2 <- mirt(data2, 1, itemtype = 'Rasch', pars=mod2values(raschfit), TOL=NaN)
#' itemfit(raschfit2, 'infit')
#'
#' # note that X2, G2, PV-Q1, and X2* do not require complete datasets
#' itemfit(raschfit, c('X2', 'G2'))
#' itemfit(raschfit, empirical.plot=1)
#' itemfit(raschfit, empirical.table=1)
#'
#'}
#'
itemfit <- function(x, fit_stats = 'S_X2', which.items = 1:extract.mirt(x, 'nitems'),
                    group.bins = 10, group.size = NA, group.fun = mean,
                    mincell = 1, mincell.X2 = 2, S_X2.tables = FALSE,
                    pv_draws = 30, boot = 1000, boot_dfapprox = 200,
                    ETrange = c(-2,2), ETpoints = 11,
                    empirical.plot = NULL, empirical.CI = .95, empirical.table = NULL,
                    method = 'EAP', Theta = NULL, impute = 0,
                    par.strip.text = list(cex = 0.7),
                    par.settings = list(strip.background = list(col = '#9ECAE1'),
                                        strip.border = list(col = "black")), ...){

    fn <- function(ind, Theta, obj, vals, ...){
        tmpobj <- obj
        tmpdat <- imputeMissing(obj, Theta[[ind]], warn=FALSE)
        tmpmod <- mirt(tmpdat, model=1, TOL=NA,
                       technical=list(customK=obj@Data$K, message=FALSE, warn=FALSE))
        tmpobj@Data <- tmpmod@Data
        whc <- 1L:length(Theta)
        return(itemfit(tmpobj, Theta=Theta[[sample(whc[-ind], 1L)]], ...))
    }
    PV_itemfit <- function(mod, which.items = 1:extract.mirt(mod, 'nitems'),
                           draws = 100, ...){
        pv <- fscores(mod, plausible.draws = draws, ...)
        draws <- length(pv)
        df.X2 <- Q1 <- matrix(NA, length(which.items), draws)
        for (i in seq_len(draws)){
            tmp <- itemfit(mod, fit_stats='X2', which.items=which.items,
                           Theta = pv[[i]], ...)
            Q1[,i] <- tmp$X2
            df.X2[,i] <- tmp$df.X2
        }
        Q1_m <- rowMeans(Q1)
        df.X2_m <- rowMeans(df.X2)
        p.Q1 <- pchisq(Q1_m, df.X2_m, lower.tail = FALSE)
        p.Q1 <- ifelse(df.X2_m == 0, NaN, p.Q1)
        ret <- data.frame(PV_Q1=Q1_m, df.PV_Q1=df.X2_m, p.PV_Q1=p.Q1)
        ret
    }
    boot_PV <- function(mod, org, is_NA, which.items = 1:extract.mirt(mod, 'nitems'),
                        itemtype, boot = 1000, draws = 30, verbose = FALSE, ...){
        pb_fun <- function(ind, mod, N, sv, which.items, draws, itemtype, ...){
            count <- 0L
            while(TRUE){
                count <- count + 1L
                if(count == 20)
                    stop('20 consecutive parametric bootstraps failed for PV_Q1*', call.=FALSE)
                dat <- simdata(model=mod, N=N)
                dat[is_NA] <- NA
                mod2 <- mirt(dat, model, itemtype=itemtype,
                             verbose=FALSE, pars=sv, technical=list(warn=FALSE))
                if(!extract.mirt(mod2, 'converged')) next
                tmp <- PV_itemfit(mod2, which.items=which.items, draws=draws, ...)
                ret <- tmp$p.PV_Q1
                if(any(is.nan(ret) | is.na(ret))) next
                break
            }
            ret
        }
        N <- nrow(extract.mirt(mod, 'data'))
        retQ1 <- matrix(NA, boot, length(which.items))
        stopifnot(nrow(org) == length(which.items))
        model <- extract.mirt(mod, 'model')
        sv <- mod2values(mod)
        retQ1 <- mySapply(1L:boot, pb_fun, mod=mod, N=N, sv=sv, itemtype=itemtype,
                          which.items=which.items, draws=draws, ...)
        if(nrow(retQ1) == 1L) retQ1 <- t(retQ1)
        Q1 <- (1 + rowSums(org$p.PV_Q1 > t(retQ1), na.rm = TRUE)) / (1 + boot)
        ret <- data.frame("p.PV_Q1_star"=Q1)
        ret
    }
    StoneFit <- function(mod, is_NA, which.items = 1:extract.mirt(mod, 'nitems'), itemtype,
                         dfapprox = FALSE, boot = 1000, ETrange = c(-2,2), ETpoints = 11,
                         verbose = FALSE, ...){
        X2star <- function(mod, which.items, ETrange, ETpoints, itemtype, ...){
            sv <- mod2values(mod)
            sv$est <- FALSE
            Theta <- matrix(seq(ETrange[1L], ETrange[2L], length.out=ETpoints))
            dat <- extract.mirt(mod, 'data')
            Emod <- mirt(dat, 1, itemtype=itemtype,
                         pars=sv, verbose=FALSE,
                         technical=list(storeEtable=TRUE, customTheta=Theta))
            Etable <- Emod@Internals$Etable[[1]]$r1
            itemloc <- extract.mirt(mod, 'itemloc')
            X2 <- rep(NA, ncol(dat))
            for(i in seq_len(length(which.items))){
                pick <- itemloc[which.items[i]]:(itemloc[which.items[i]+1L] - 1L)
                O <- Etable[ ,pick]
                item <- extract.item(mod, which.items[i])
                E <- probtrace(item, Theta) * rowSums(O)
                X2[which.items[i]] <- sum((O - E)^2 / E, na.rm = TRUE)
            }
            X2[which.items]
        }
        pb_fun <- function(ind, is_NA, mod, N, model, itemtype, sv, which.items, ETrange,
                           ETpoints, ...){
            count <- 0L
            while(TRUE){
                count <- count + 1L
                if(count == 20)
                    stop('20 consecutive parametric bootstraps failed for X2*', call.=FALSE)
                dat <- simdata(model=mod, N=N)
                dat[is_NA] <- NA
                mod2 <- mirt(dat, model, itemtype=itemtype, verbose=FALSE, pars=sv,
                             technical=list(warn=FALSE))
                if(!extract.mirt(mod2, 'converged')) next
                ret <- X2star(mod2, which.items=which.items, ETrange=ETrange,
                              ETpoints=ETpoints, itemtype=itemtype, ...)
                if(any(is.nan(ret) | is.na(ret))) next
                break
            }
            ret
        }

        N <- nrow(extract.mirt(mod, 'data'))
        X2bs <- matrix(NA, boot, length(which.items))
        org <- X2star(mod, which.items=which.items, itemtype=itemtype,
                      ETrange=ETrange, ETpoints=ETpoints, ...)
        stopifnot(length(org) == length(which.items))
        sv <- mod2values(mod)
        model <- extract.mirt(mod, 'model')
        X2bs <- mySapply(1L:boot, pb_fun, mod=mod, N=N, model=model, is_NA=is_NA,
                         itemtype=itemtype, sv=sv, which.items=which.items,
                         ETrange=ETrange, ETpoints=ETpoints, ...)
        if(nrow(X2bs) == 1L) X2bs <- t(X2bs)
        if(dfapprox){
            M <- colMeans(X2bs)
            V <- apply(X2bs, 2, var)
            upsilon <- 2 * M^2 / V
            gamma <- M / upsilon
            df <- upsilon
            for(i in seq_len(length(which.items))){
                item <- extract.item(mod, which.items[i])
                df[i] <- upsilon[i] - sum(item@est)
            }
            ret <- data.frame(X2_star_scaled=org/gamma, df.X2_star_scaled=df,
                              p.X2_star_scaled=pchisq(org/gamma, df, lower.tail=FALSE))
        } else {
            p <- apply(t(X2bs) > org, 1, mean)
            ret <- data.frame(X2_star=org, p.X2_star=p)
        }
        ret
    }

    if(missing(x)) missingMsg('x')
    if(is(x, 'MixedClass'))
        stop('MixedClass objects are not supported', call.=FALSE)
    if(!is.null(empirical.plot) && !is.null(empirical.table))
        stop('Please select empirical.plot or empirical.table, not both', call.=FALSE)
    if(!all(fit_stats %in% c('S_X2', 'Zh', 'X2', 'G2', 'infit', 'PV_Q1', 'PV_Q1*', 'X2*', 'X2*_df')))
        stop('Unsupported fit_stats element requested', call.=FALSE)
    if(any(c('X2', 'G2', 'PV_Q1', 'PV_Q1*') %in% fit_stats) && extract.mirt(x, 'nfact') > 1L)
        stop('X2, G2, PV_Q1, or PV_Q1* are for unidimensional models only', call.=FALSE)
    S_X2 <- 'S_X2' %in% fit_stats
    Zh <- 'Zh' %in% fit_stats
    X2 <- 'X2' %in% fit_stats
    G2 <- 'G2' %in% fit_stats
    infit <- 'infit' %in% fit_stats
    if(!is.null(empirical.plot))
        which.items <- empirical.plot
    if(!is.null(empirical.table))
        which.items <- empirical.table
    which.items <- sort(which.items)
    if(!is.null(empirical.plot) || !is.null(empirical.table)){
        Zh <- FALSE
        if(length(which.items) > 1L)
            stop('Plots and tables only supported for 1 item at a time', call.=FALSE)
    }
    stopifnot(is.numeric(empirical.CI))
    if(!is.null(empirical.table) || !is.null(empirical.plot)){
        X2 <- TRUE
        S_X2 <- Zh <- infit <- G2 <- FALSE
    }
    J <- ncol(x@Data$data)
    if(any(is.na(x@Data$data)) && (Zh || S_X2 || infit) && impute == 0)
        stop('Only X2, G2, PV_Q1, PV_Q1*, X2*, and X2*_df can be computed with missing data.
             Consider using imputed datasets', call.=FALSE)

    if(is(x, 'MultipleGroupClass') || is(x, 'DiscreteClass')){
        discrete <- is(x, 'DiscreteClass')
        if(discrete)
            for(g in seq_len(x@Data$ngroups))
                x@ParObjects$pars[[g]]@ParObjects$pars[[J+1L]]@est[] <- FALSE
        ret <- vector('list', x@Data$ngroups)
        if(is.null(Theta))
            Theta <- fscores(x, method=method, full.scores=TRUE, plausible.draws=impute, ...)
        for(g in seq_len(x@Data$ngroups)){
            if(impute > 0L){
                tmpTheta <- vector('list', impute)
                for(i in seq_len(length(tmpTheta)))
                    tmpTheta[[i]] <- Theta[[i]][x@Data$groupNames[g] == x@Data$group, , drop=FALSE]
            } else tmpTheta <- Theta[x@Data$groupNames[g] == x@Data$group, , drop=FALSE]
            tmp_obj <- MGC2SC(x, g)
            ret[[g]] <- itemfit(tmp_obj, fit_stats=fit_stats, group.size=group.size, group.bins=group.bins,
                                group.fun=group.fun, mincell=mincell, mincell.X2=mincell.X2,
                                S_X2.tables=S_X2.tables, empirical.plot=empirical.plot,
                                empirical.table=empirical.table,
                                Theta=tmpTheta, empirical.CI=empirical.CI, method=method,
                                impute=impute, discrete=discrete, ...)
        }
        names(ret) <- x@Data$groupNames
        if(extract.mirt(x, 'ngroups') == 1L) return(ret[[1L]])
        return(ret)
    }
    dots <- list(...)
    discrete <- dots$discrete
    discrete <- ifelse(is.null(discrete), FALSE, discrete)
    if(impute != 0 && !is(x, 'MultipleGroupClass')){
        if(impute == 0)
            stop('Fit statistics cannot be computed when there are missing data. Pass a suitable
                 impute argument to compute statistics following multiple data
                 imputations', call.=FALSE)
        if(sum(is.na(x@Data$data)) / length(x@Data$data) > .10)
            warning('Imputations for large amounts of missing data may be overly conservative', call.=FALSE)
        stopifnot(impute > 1L)
        if(is.null(Theta))
            Theta <- fscores(x, plausible.draws = impute, method = ifelse(method == 'MAP', 'MAP', 'EAP'), ...)
        collect <- vector('list', impute)
        vals <- mod2values(x)
        vals$est <- FALSE
        collect <- myLapply(1L:impute, fn, Theta=Theta, obj=x, vals=vals, fit_stats=fit_stats,
                            group.size=group.size, group.bins=group.bins,
                            mincell=mincell, mincell.X2=mincell.X2,
                            S_X2.tables=S_X2.tables, empirical.plot=empirical.plot,
                            empirical.CI=empirical.CI, empirical.table=empirical.table,
                            method=method, impute=0, discrete=discrete, ...)
        ave <- SD <- collect[[1L]]
        pick1 <- 1:nrow(ave)
        pick2 <- sapply(ave, is.numeric)
        ave[pick1, pick2] <- SD[pick1, pick2] <- 0
        for(i in seq_len(impute))
            ave[pick1, pick2] <- ave[pick1, pick2] + collect[[i]][pick1, pick2]
        ave[pick1, pick2] <- ave[pick1, pick2]/impute
        for(i in seq_len(impute))
            SD[pick1, pick2] <- SD[pick1, pick2] + (ave[pick1, pick2] - collect[[i]][pick1, pick2])^2
        SD[pick1, pick2] <- sqrt(SD[pick1, pick2]/impute)
        SD$item <- paste0('SD_', SD$item)
        SD <- rbind(NA, SD)
        ret <- rbind(ave, SD)
        class(ret) <- c('mirt_df', 'data.frame')
        return(ret)
    }
    if(S_X2.tables || discrete) Zh <- X2 <- FALSE
    ret <- data.frame(item=colnames(x@Data$data)[which.items])
    itemloc <- x@Model$itemloc
    pars <- x@ParObjects$pars
    if(all(x@Model$itemtype %in% c('Rasch', '2PL', 'rsm', 'gpcm')) && infit){
        infit <- FALSE
        oneslopes <- rep(FALSE, length(x@Model$itemtype))
        slope <- x@ParObjects$pars[[1L]]@par[1L]
        for(i in seq_len(length(x@Model$itemtype)))
            oneslopes[i] <- closeEnough(x@ParObjects$pars[[i]]@par[1L], slope-1e-10, slope+1e-10)
        if(all(oneslopes)) infit <- TRUE
    } else infit <- FALSE
    if(Zh || infit){
        if(is.null(Theta))
            Theta <- fscores(x, verbose=FALSE, full.scores=TRUE, method=method, ...)
        prodlist <- attr(pars, 'prodlist')
        nfact <- x@Model$nfact + length(prodlist)
        fulldata <- x@Data$fulldata[[1L]]
        if(any(Theta %in% c(Inf, -Inf))){
            for(i in 1L:ncol(Theta)){
                tmp <- Theta[,i]
                tmp[tmp %in% c(-Inf, Inf)] <- NA
                Theta[Theta[,i] == Inf, i] <- max(tmp, na.rm=TRUE) + .1
                Theta[Theta[,i] == -Inf, i] <- min(tmp, na.rm=TRUE) - .1
            }
        }
        N <- nrow(Theta)
        itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)
        for (i in which.items)
            itemtrace[ ,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(x=pars[[i]], Theta=Theta)
        log_itemtrace <- log(itemtrace)
        LL <- log_itemtrace * fulldata
        Lmatrix <- matrix(LL[as.logical(fulldata)], N, J)
        mu <- sigma2 <- rep(0, J)
        for(item in which.items){
            P <- itemtrace[ ,itemloc[item]:(itemloc[item+1L]-1L)]
            log_P <- log_itemtrace[ ,itemloc[item]:(itemloc[item+1L]-1L)]
            mu[item] <- sum(P * log_P)
            for(i in seq_len(ncol(P)))
                for(j in seq_len(ncol(P)))
                    if(i != j)
                        sigma2[item] <- sigma2[item] + sum(P[,i] * P[,j] *
                                                               log_P[,i] * log(P[,i]/P[,j]))
        }
        tmp <- (colSums(Lmatrix) - mu) / sqrt(sigma2)
        if(Zh) ret$Zh <- tmp[which.items]
        #if all Rasch models, infit and outfit
        if(infit){
            attr(x, 'inoutfitreturn') <- TRUE
            pf <- personfit(x, method=method, Theta=Theta)
            z2 <- pf$resid^2 / pf$W
            outfit <- colSums(z2) / N
            q.outfit <- sqrt(colSums((pf$C / pf$W^2) / N^2) - 1 / N)
            q.outfit[q.outfit > 1.4142] <- 1.4142
            z.outfit <- (outfit^(1/3) - 1) * (3/q.outfit) + (q.outfit/3)
            infit <- colSums(pf$W * z2) / colSums(pf$W)
            q.infit <- sqrt(colSums(pf$C - pf$W^2) / colSums(pf$W)^2)
            q.infit[q.infit > 1.4142] <- 1.4142
            z.infit <- (infit^(1/3) - 1) * (3/q.infit) + (q.infit/3)
            ret$outfit <- outfit[which.items]
            ret$z.outfit <- z.outfit[which.items]
            ret$infit <- infit[which.items]
            ret$z.infit <- z.infit[which.items]
        }
    }
    if(( (X2 || G2) || !is.null(empirical.plot) || !is.null(empirical.table)) && x@Model$nfact == 1L){
        if(is.null(Theta))
            Theta <- fscores(x, verbose=FALSE, full.scores=TRUE, method=method, ...)
        nfact <- ncol(Theta)
        prodlist <- attr(pars, 'prodlist')
        fulldata <- x@Data$fulldata[[1L]]
        if(any(Theta %in% c(Inf, -Inf))){
            for(i in seq_len(ncol(Theta))){
                tmp <- Theta[,i]
                tmp[tmp %in% c(-Inf, Inf)] <- NA
                Theta[Theta[,i] == Inf, i] <- max(tmp, na.rm=TRUE) + .1
                Theta[Theta[,i] == -Inf, i] <- min(tmp, na.rm=TRUE) - .1
            }
        }
        ord <- order(Theta[,1L])
        fulldata <- fulldata[ord,]
        pick <- !is.na(x@Data$data)
        pick <- pick[ord, ]
        Theta <- Theta[ord, , drop = FALSE]
        den <- dnorm(Theta, 0, .5)
        den <- den / sum(den)
        cumTheta <- cumsum(den)
        if(!is.na(group.size)){
            Groups <- rep(20, length(ord))
            ngroups <- ceiling(nrow(fulldata) / group.size)
            weight <- 1/ngroups
            for(i in seq_len(length(Groups)))
                Groups[round(cumTheta,2) >= weight*(i-1) & round(cumTheta,2) < weight*i] <- i
        } else {
            ngroups <- group.bins
            Groups <- rep(1:group.bins, each = floor(length(ord) / ngroups))
            if(length(ord) %% ngroups > 0L){
                c1 <- length(ord) %% ngroups
                Groups <- c(rep(1, floor(c1/2)), Groups)
                Groups <- c(Groups, rep(ngroups, c1 - floor(c1/2)))
            }
        }
        X2.value <- G2.value <- df.G2 <- df.X2 <- numeric(J)
        if(!is.null(empirical.plot)){
            if(nfact > 1L) stop('Cannot make empirical plot for multidimensional models', call.=FALSE)
            if(!is.numeric(empirical.plot)){
                inames <- colnames(x@Data$data)
                ind <- 1L:length(inames)
                empirical.plot <- ind[inames == empirical.plot]
            }
            empirical.plot_points <- matrix(NA, length(unique(Groups)), x@Data$K[empirical.plot] + 2L)
        }
        if(!is.null(empirical.table)){
            Etable <- vector('list', ngroups)
            mtheta_nms <- numeric(ngroups)
        }
        for (i in which.items){
            for(j in unique(Groups)){
                dat <- fulldata[Groups == j & pick[,i], itemloc[i]:(itemloc[i+1] - 1), drop = FALSE]
                colnames(dat) <- paste0("cat_", sort(unique(extract.mirt(x, "data")[,i])))
                if(nrow(dat) <= 1L) next
                r <- colSums(dat)
                N <- nrow(dat)
                mtheta <- matrix(group.fun(Theta[Groups == j & pick[,i],]), nrow=1)
                if(!is.null(empirical.plot)){
                    tmp <- r/N
                    empirical.plot_points[j, ] <- c(mtheta, N, tmp)
                }
                P <- ProbTrace(x=pars[[i]], Theta=mtheta)
                if(is.null(empirical.table) && any(N * P < mincell.X2)){
                    while(TRUE){
                        wch <- which(N * P < mincell.X2)
                        if(!length(wch) || length(r) == 1L) break
                        for(p in wch){
                            if(p == 1L){
                                r[2L] <- r[1L] + r[2L]
                                P[2L] <- P[1L] + P[2L]
                            } else {
                                r[p-1L] <- r[p] + r[p-1L]
                                P[p-1L] <- P[p] + P[p-1L]
                            }
                        }
                        r <- r[-wch]
                        P <- P[-wch]
                    }
                    if(length(r) == 1L) next
                }
                E <- N*P
                if(!is.null(empirical.table)){
                    Etable[[j]] <- data.frame(Observed=r, Expected=as.vector(E),
                                              z.Residual=as.vector(sqrt((r - E)^2 / E) * sign(r-E)))
                    mtheta_nms[j] <- mtheta
                } else {
                    X2.value[i] <- X2.value[i] + sum((r - E)^2 / E)
                    df.X2[i] <- df.X2[i] + length(r) - 1L
                    tmp <- r * log(r/E)
                    tmp <- tmp[is.finite(tmp)]
                    if(length(tmp) > 1L){
                        G2.value[i] <- G2.value[i] + 2*sum(tmp)
                        df.G2[i] <- df.G2[i] + length(tmp) - 1L
                    }
                }
            }
            if(!is.null(empirical.table)){
                names(Etable) <- paste0('theta = ', round(mtheta_nms, 4))
                return(Etable)
            }
            df.X2[i] <- df.X2[i] - sum(pars[[i]]@est)
            df.G2[i] <- df.G2[i] - sum(pars[[i]]@est)
        }
        X2.value[X2.value == 0] <- NA
        G2.value[G2.value == 0] <- NA
        if(!is.null(empirical.plot)){
            K <- x@Data$K[empirical.plot]
            EPCI.lower <- EPCI.upper <- NULL
            if(K == 2 && empirical.CI != 0){
                p.L <- function(x, alpha) if (x[1] == 0) 0 else qbeta(alpha, x[1], x[2] - x[1] + 1)
                p.U <- function(x, alpha) if (x[1] == x[2]) 1 else
                    qbeta(1 - alpha, x[1] + 1, x[2] - x[1])
                N <- empirical.plot_points[,2]
                O <- empirical.plot_points[,ncol(empirical.plot_points)] * N
                EPCI.lower <- apply(cbind(O, N), 1, p.L, (1-empirical.CI)/2)
                EPCI.upper <- apply(cbind(O, N), 1, p.U, (1-empirical.CI)/2)
            }
            theta <- seq(-4,4, length.out=max(c(50, ngroups)))
            ThetaFull <- thetaComb(theta, nfact)
            empirical.plot_P <- ProbTrace(pars[[empirical.plot]], ThetaFull)
            empirical.plot_points <- empirical.plot_points[,-2]
            colnames(empirical.plot_points) <- c('theta', paste0('p.', 1:K))
            while(nrow(empirical.plot_points) < nrow(empirical.plot_P))
                empirical.plot_points <- rbind(empirical.plot_points,
                                               rep(NA, length(empirical.plot_points[1,])))
            plt.1 <- data.frame(id = 1:nrow(ThetaFull), Theta=ThetaFull, P=empirical.plot_P)
            plt.1 <- reshape(plt.1, varying = 3:ncol(plt.1), direction = 'long', timevar = 'cat')
            plt.2 <- data.frame(id = 1:nrow(empirical.plot_points), empirical.plot_points)
            plt.2 <- reshape(plt.2, varying = 3:ncol(plt.2), direction = 'long', timevar = 'cat')
            plt <- cbind(plt.1, plt.2)
            if(K == 2) plt <- plt[plt$cat != 1, ]
            plt$cat <- factor(paste0('Category ', plt$cat - 1 + extract.mirt(x, 'mins')[which.items]))
            return(xyplot(if(K == 2) P~Theta else P ~ Theta|cat, plt, groups = cat,
                          main = paste('Empirical plot for item', empirical.plot),
                            ylim = c(-0.1,1.1), xlab = expression(theta), ylab=expression(P(theta)),
                          EPCI.lower=EPCI.lower, EPCI.upper=EPCI.upper,
                          panel = function(x, y, groups, subscripts, EPCI.lower, EPCI.upper, ...){
                              panel.xyplot(x=x, y=y, groups=groups, type='l',
                                           subscripts=subscripts, ...)
                              panel.points(cbind(plt$theta[subscripts], plt$p[subscripts]),
                                           col='black', ...)
                              if(!is.null(EPCI.lower)){
                                  theta <- na.omit(plt$theta)
                                  for(i in 1:length(theta))
                                      panel.lines(c(theta[i], theta[i]), c(EPCI.lower[i],
                                                                           EPCI.upper[i]),
                                                  lty = 2, col = 'red')
                              }
                          },
                          par.strip.text=par.strip.text, par.settings=par.settings, ...))
        }
        if(X2){
            ret$X2 <- X2.value[which.items]
            ret$df.X2 <- df.X2[which.items]
            ret$p.X2 <- suppressWarnings(pchisq(ret$X2, ret$df.X2, lower.tail=FALSE))
            ret$df.X2[ret$df.X2 <= 0] <- 0
            ret$p.X2[ret$df.X2 <= 0] <- NaN
        }
        if(G2){
            ret$G2 <- G2.value[which.items]
            ret$df.G2 <- df.G2[which.items]
            ret$p.G2 <- suppressWarnings(pchisq(ret$G2, ret$df.G2, lower.tail=FALSE))
            ret$df.G2[ret$df.G2 <= 0] <- 0
            ret$p.G2[ret$df.G2 <= 0] <- NaN
        }
    }
    if(S_X2){
        dat <- x@Data$data
        adj <- x@Data$mins
        dat <- t(t(dat) - adj)
        S_X2 <- df.S_X2 <- rep(NA, J)
        O <- makeObstables(dat, x@Data$K, which.items=which.items)
        Nk <- rowSums(O[[which(sapply(O, is.matrix))[1L]]])
        dots <- list(...)
        QMC <- ifelse(is.null(dots$QMC), FALSE, dots$QMC)
        quadpts <- dots$quadpts
        if(is.null(quadpts) && QMC) quadpts <- 5000L
        if(is.null(quadpts)) quadpts <- select_quadpts(x@Model$nfact)
        if(x@Model$nfact > 3L && !QMC && method %in% c('EAP', 'EAPsum') && !discrete)
            warning('High-dimensional models should use quasi-Monte Carlo integration. Pass QMC=TRUE',
                    call.=FALSE)
        theta_lim <- dots$theta_lim
        if(is.null(theta_lim)) theta_lim <- c(-6,6)
        gp <- ExtractGroupPars(pars[[length(pars)]])
        E <- EAPsum(x, S_X2 = TRUE, gp = gp, CUSTOM.IND=x@Internals$CUSTOM.IND, den_fun=mirt_dmvnorm,
                    quadpts=quadpts, theta_lim=theta_lim, discrete=discrete, QMC=QMC,
                    which.items=which.items)
        for(i in which.items)
            E[[i]] <- E[[i]] * Nk
        coll <- collapseCells(O, E, mincell=mincell)
        if(S_X2.tables) return(list(O.org=O, E.org=E, O=coll$O, E=coll$E))
        O <- coll$O
        E <- coll$E
        for(i in seq_len(J)){
            if (is.null(dim(O[[i]])) || is.null(E[[i]])) next
            S_X2[i] <- sum((O[[i]] - E[[i]])^2 / E[[i]], na.rm = TRUE)
            df.S_X2[i] <- sum(!is.na(E[[i]])) - nrow(E[[i]]) - sum(pars[[i]]@est)
        }
        df.S_X2 <- df.S_X2 - sum(pars[[J+1L]]@est)
        df.S_X2[df.S_X2 < 0] <- 0
        S_X2[df.S_X2 == 0] <- NaN
        ret$S_X2 <- S_X2[which.items]
        ret$df.S_X2 <- df.S_X2[which.items]
        ret$p.S_X2 <- suppressWarnings(pchisq(ret$S_X2, ret$df.S_X2, lower.tail=FALSE))
    }
    itemtype <- extract.mirt(x, 'itemtype')
    if(any(c('PV_Q1', 'PV_Q1*') %in% fit_stats)){
        tmp <- PV_itemfit(x, which.items=which.items, draws=pv_draws, itemtype=itemtype, ...)
        ret <- cbind(ret, tmp)
    }
    is_NA <- is.na(x@Data$data)
    if('PV_Q1*' %in% fit_stats){
        tmp <- boot_PV(x, is_NA=is_NA, org=tmp, which.items=which.items,
                       itemtype=itemtype, boot=boot, draws=pv_draws, ...)
        ret <- cbind(ret, tmp)
    }
    if('X2*' %in% fit_stats){
        tmp <- StoneFit(x, is_NA=is_NA, which.items=which.items, boot=boot, dfapprox=FALSE,
                        itemtype=itemtype, ETrange=ETrange, ETpoints=ETpoints, ...)
        ret <- cbind(ret, tmp)
    }
    if('X2*_df' %in% fit_stats){
        tmp <- StoneFit(x, is_NA=is_NA, which.items=which.items, boot=boot_dfapprox, dfapprox=TRUE,
                        itemtype=itemtype, ETrange=ETrange, ETpoints=ETpoints, ...)
        ret <- cbind(ret, tmp)
    }
    class(ret) <- c('mirt_df', 'data.frame')
    return(ret)
}
