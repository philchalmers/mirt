#' Item fit statistics
#'
#' \code{itemfit} calculates the Zh values from Drasgow, Levine and Williams (1985),
#' \eqn{\chi^2} and \eqn{G^2} values for unidimensional models, and S-X2 statistics for unidimensional and
#' multidimensional models (Kang & Chen, 2007; Orlando & Thissen, 2000).
#' For Rasch, partial credit, and rating scale models infit and outfit statistics are
#' also produced. Poorly fitting items should be inspected with the empirical plots/tables
#' for unidimensional models, otherwise \code{\link{itemGAM}} can be used to diagnose
#' where the functional form of the IRT model was misspecified, or models can be refit using
#' more flexible semi-parametric response models (e.g., \code{itemtype = 'spline'}).
#' For discrete models, only the S-X2 statistic will be calculated.
#'
#' @aliases itemfit
#' @param x a computed model object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{DiscreteClass}
#' @param which.items an integer vector indicating which items to test for fit.
#'   Default tests all possible items
#' @param S_X2 logical; calculate the S_X2 statistic?
#' @param Zh logical; calculate Zh?
#' @param X2 logical; calculate the X2 statistic for unidimensional models?
#' @param G2 logical; calculate the G2 statistic for unidimensional models?
#' @param infit logical; calculate the infit/outfit values for unidimensional Rasch models? If the models
#'   are not from the Rasch family then this will be ignored
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
#' @param group.bins the number of bins to use when \code{X2 = TRUE}. For example,
#'   setting \code{group.bins = 10} will will compute Yen's (1981) Q1 statistic
#' @param group.fun function used when \code{X2} or \code{G2} are computed. Determines the central
#'   tendancy measure within each partitioned group. E.g., setting \code{group.fun = median} will
#'   obtain the median of each respective ability estimate in each subgroup (this is what was used
#'   by Bock, 1972)
#' @param empirical.plot a single numeric value or character of the item name indicating which
#'   item to plot (via \code{itemplot}) and overlay with the empirical \eqn{\theta} groupings (see
#'   \code{empirical.CI}). Useful for plotting the expected bins after passing \code{X2 = TRUE} or
#'   \code{G2 = TRUE}
#' @param empirical.table a single numeric value or character of the item name indicating which
#'   item table of expected values should be returned. Useful for plotting the expected bins after
#'   passing \code{X2 = TRUE} or \code{G2 = TRUE}
#' @param empirical.CI a numeric value indicating the width of the empirical confidence interval
#'   ranging between 0 and 1 (default of 0 plots not interval). For example, a 95% confidence
#'   interval would be plotted when \code{empirical.CI = .95}. Only applicable to dichotomous items
#' @param method type of factor score estimation method. See \code{\link{fscores}} for more detail
#' @param Theta a matrix of factor scores for each person used for statistics that require
#'   empirical estimates. If supplied, arguments typically passed to \code{fscores()} will be
#'   ignored and these values will be used instead. Also required when estimating statistics
#'   with missing data via imputation
#' @param impute a number indicating how many imputations to perform (passed to
#'   \code{\link{imputeMissing}}) when there are missing data present.
#'   Will return a data.frame object with the mean estimates
#'   of the stats and their imputed standard deviations
#' @param digits number of digits to round result to. Default is 4
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
#' in two or more nominal categories. Psychometrika, 37, 29-51.
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
#' Wright B. D. & Masters, G. N. (1982). \emph{Rating scale analysis}. MESA Press.
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
#' items <- rep('dich', 20)
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
#' itemfit(x, X2=TRUE)
#' itemfit(x, group.bins=15, empirical.plot = 1) #empirical item plot with 15 points
#' itemfit(x, group.bins=15, empirical.plot = 21)
#'
#' #empirical tables
#' itemfit(x, empirical.table=1)
#' itemfit(x, empirical.table=21)
#'
#' #infit/outfit statistics. method='ML' agrees better with eRm package
#' itemfit(raschfit, method = 'ML', infit = TRUE) #infit and outfit stats
#'
#' #same as above, but inputting ML estimates instead
#' Theta <- fscores(raschfit, method = 'ML')
#' itemfit(raschfit, Theta=Theta, infit = TRUE)
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
#' itemfit(mod, X2 = TRUE)
#' itemfit(mod, empirical.plot = 1)
#'
#' # collapsed tables (see mincell.X2) for X2 and G2
#' itemfit(mod, empirical.table = 1)
#'
#' mod2 <- mirt(dat, 1, 'Rasch')
#' itemfit(mod2, infit = TRUE)
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
#' itemfit(raschfit, impute = 10, infit = TRUE)
#'
#' #alternative route: use only valid data, and create a model with the previous parameter estimates
#' data2 <- na.omit(data)
#' raschfit2 <- mirt(data2, 1, itemtype = 'Rasch', pars=mod2values(raschfit), TOL=NaN)
#' itemfit(raschfit2, infit = TRUE)
#'
#' # note that X2 and G2 do not require complete datasets
#' itemfit(raschfit, X2=TRUE, G2 = TRUE, S_X2=FALSE)
#' itemfit(raschfit, empirical.plot=1)
#' itemfit(raschfit, empirical.table=1)
#'
#'}
#'
itemfit <- function(x, which.items = 1:extract.mirt(x, 'nitems'),
                    S_X2 = TRUE, Zh = FALSE, X2 = FALSE, G2 = FALSE, infit = FALSE,
                    group.bins = 10, group.size = NA, group.fun = mean,
                    mincell = 1, mincell.X2 = 2, S_X2.tables = FALSE,
                    empirical.plot = NULL, empirical.CI = .95, empirical.table = NULL,
                    method = 'EAP', Theta = NULL, impute = 0, digits = 4,
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
        return(itemfit(tmpobj, Theta=Theta[[sample(whc[-ind], 1L)]], digits = Inf, ...))
    }

    if(missing(x)) missingMsg('x')
    if(is(x, 'MixedClass'))
        stop('MixedClass objects are not supported', call.=FALSE)
    discrete <- FALSE
    if(is(x, 'DiscreteClass')){
        if(is.null(Theta))
            Theta <- fscores(x, method=method, full.scores=TRUE, plausible.draws=impute, ...)
        class(x) <- 'MultipleGroupClass'
        discrete <- TRUE
    }
    if(!is.null(empirical.plot) && !is.null(empirical.table))
        stop('Please select empirical.plot or empirical.table, not both', call.=FALSE)
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

    stopifnot(Zh || X2 || S_X2 || infit || G2)
    stopifnot(is.numeric(empirical.CI))
    if(!is.null(empirical.table) || !is.null(empirical.plot)){
        X2 <- TRUE
        S_X2 <- Zh <- infit <- G2 <- FALSE
    }
    if(any(is.na(x@Data$data)) && (Zh || S_X2 || infit) && impute == 0)
        stop('Only X2 and G2 can be computed with missing data. Consider using imputed datasets', call.=FALSE)

    if(impute != 0 && !is(x, 'MultipleGroupClass')){
        if(impute == 0)
            stop('Fit statistics cannot be computed when there are missing data. Pass a suitable
                 impute argument to compute statistics following multiple data
                 imputations', call.=FALSE)
        if(sum(is.na(x@Data$data)) / length(x@Data$data) > .10)
            warning('Imputations for large amounts of missing data may be overly conservative', call.=FALSE)
        stopifnot(impute > 1L)
        Theta <- fscores(x, plausible.draws = impute, method = ifelse(method == 'MAP', 'MAP', 'EAP'), ...)
        collect <- vector('list', impute)
        vals <- mod2values(x)
        vals$est <- FALSE
        collect <- myLapply(1L:impute, fn, Theta=Theta, obj=x, vals=vals, S_X2=S_X2,
                            Zh=Zh, X2=X2, G2=G2, group.size=group.size, group.bins=group.bins,
                            mincell=mincell, mincell.X2=mincell.X2, infit=infit,
                            S_X2.tables=S_X2.tables, empirical.plot=empirical.plot,
                            empirical.CI=empirical.CI, empirical.table=empirical.table,
                            method=method, impute=0, discrete=discrete, ...)
        ave <- SD <- collect[[1L]]
        pick1 <- 1:nrow(ave)
        pick2 <- sapply(ave, is.numeric)
        ave[pick1, pick2] <- SD[pick1, pick2] <- 0
        for(i in 1L:impute)
            ave[pick1, pick2] <- ave[pick1, pick2] + collect[[i]][pick1, pick2]
        ave[pick1, pick2] <- ave[pick1, pick2]/impute
        for(i in 1L:impute)
            SD[pick1, pick2] <- SD[pick1, pick2] + (ave[pick1, pick2] - collect[[i]][pick1, pick2])^2
        SD[pick1, pick2] <- sqrt(SD[pick1, pick2]/impute)
        SD$item <- paste0('SD_', SD$item)
        SD <- rbind(NA, SD)
        ret <- rbind(ave, SD)
        ret[,sapply(ret, class) == 'numeric'] <- round(ret[,sapply(ret, class) == 'numeric'], digits)
        return(ret)
    }
    if(is(x, 'MultipleGroupClass')){
        ret <- vector('list', x@Data$ngroups)
        if(is.null(Theta))
            Theta <- fscores(x, method=method, full.scores=TRUE, plausible.draws=impute, ...)
        for(g in 1L:x@Data$ngroups){
            if(impute > 0L){
                tmpTheta <- vector('list', impute)
                for(i in 1L:length(tmpTheta))
                    tmpTheta[[i]] <- Theta[[i]][x@Data$groupNames[g] == x@Data$group, , drop=FALSE]
            } else tmpTheta <- Theta[x@Data$groupNames[g] == x@Data$group, , drop=FALSE]
            tmp_obj <- MGC2SC(x, g)
            ret[[g]] <- itemfit(tmp_obj, Zh=Zh, X2=X2, group.size=group.size, group.bins=group.bins,
                                group.fun=group.fun, mincell=mincell, mincell.X2=mincell.X2, infit=infit,
                                S_X2.tables=S_X2.tables, empirical.plot=empirical.plot,
                                empirical.table=empirical.table, G2=G2,
                                Theta=tmpTheta, empirical.CI=empirical.CI, method=method,
                                impute=impute, discrete=discrete, digits=digits, S_X2=S_X2, ...)
        }
        names(ret) <- x@Data$groupNames
        if(extract.mirt(x, 'ngroups') == 1L) return(ret[[1L]])
        return(ret)
    }
    dots <- list(...)
    discrete <- dots$discrete
    discrete <- ifelse(is.null(discrete), FALSE, discrete)
    if(S_X2.tables || discrete) Zh <- X2 <- FALSE
    ret <- data.frame(item=colnames(x@Data$data)[which.items])
    J <- ncol(x@Data$data)
    itemloc <- x@Model$itemloc
    pars <- x@ParObjects$pars
    if(all(x@Model$itemtype %in% c('Rasch', 'rsm', 'gpcm')) && infit){
        infit <- FALSE
        oneslopes <- rep(FALSE, length(x@Model$itemtype))
        slope <- x@ParObjects$pars[[1L]]@par[1L]
        for(i in 1L:length(x@Model$itemtype))
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
            for(i in 1L:ncol(P))
                for(j in 1L:ncol(P))
                    if(i != j)
                        sigma2[item] <- sigma2[item] + sum(P[,i] * P[,j] *
                                                               log_P[,i] * log(P[,i]/P[,j]))
        }
        tmp <- (colSums(Lmatrix) - mu) / sqrt(sigma2)
        if(Zh) ret$Zh <- tmp[which.items]
        #if all Rasch models, infit and outfit
        if(all(x@Model$itemtype %in% c('Rasch', 'rsm', 'gpcm'))){
            oneslopes <- rep(FALSE, length(x@Model$itemtype))
            for(i in 1L:length(x@Model$itemtype))
                oneslopes[i] <- closeEnough(x@ParObjects$pars[[i]]@par[1L], 1-1e-10, 1+1e-10)
            if(all(oneslopes)){
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
    }
    if(( (X2 || G2) || !is.null(empirical.plot) || !is.null(empirical.table)) && x@Model$nfact == 1L){
        if(is.null(Theta))
            Theta <- fscores(x, verbose=FALSE, full.scores=TRUE, method=method, ...)
        nfact <- ncol(Theta)
        prodlist <- attr(pars, 'prodlist')
        fulldata <- x@Data$fulldata[[1L]]
        if(any(Theta %in% c(Inf, -Inf))){
            for(i in 1L:ncol(Theta)){
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
            for(i in 1L:length(Groups))
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
                if(nrow(dat) <= 1L) next
                r <- colSums(dat)
                N <- nrow(dat)
                mtheta <- matrix(group.fun(Theta[Groups == j & pick[,i],]), nrow=1)
                if(!is.null(empirical.plot)){
                    tmp <- r/N
                    empirical.plot_points[j, ] <- c(mtheta, N, tmp)
                }
                P <- ProbTrace(x=pars[[i]], Theta=mtheta)
                if(any(N * P < mincell.X2)){
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
                X2.value[i] <- X2.value[i] + sum((r - E)^2 / E)
                df.X2[i] <- df.X2[i] + length(r) - 1L
                tmp <- r * log(r/E)
                tmp <- tmp[is.finite(tmp)]
                if(length(tmp) > 1L){
                    G2.value[i] <- G2.value[i] + 2*sum(tmp)
                    df.G2[i] <- df.G2[i] + length(tmp) - 1L
                }
                if(!is.null(empirical.table)){
                    Etable[[j]] <- data.frame(Observed=r, Expected=as.vector(E),
                                              z.Residual=as.vector(sqrt((r - E)^2 / E) * sign(r-E)))
                    rownames(Etable[[j]]) <- 1:nrow(Etable[[j]]) + extract.mirt(x, 'mins')[which.items] - 1
                    mtheta_nms[j] <- mtheta
                }
            }
            if(!is.null(empirical.table)){
                Etable <- lapply(Etable, round, digits=digits)
                names(Etable) <- paste0('theta = ', round(mtheta_nms, digits))
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
            ret$p.X2 <- 1 - suppressWarnings(pchisq(ret$X2, ret$df.X2))
            ret$p.X2[ret$df.X2 <= 0] <- NaN
        }
        if(G2){
            ret$G2 <- G2.value[which.items]
            ret$df.G2 <- df.G2[which.items]
            ret$p.G2 <- 1 - suppressWarnings(pchisq(ret$G2, ret$df.G2))
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
        if(is.null(quadpts) && QMC) quadpts <- 15000L
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
        for(i in 1L:J){
            if (is.null(dim(O[[i]])) || is.null(E[[i]])) next
            S_X2[i] <- sum((O[[i]] - E[[i]])^2 / E[[i]], na.rm = TRUE)
            df.S_X2[i] <- sum(!is.na(E[[i]])) - nrow(E[[i]]) - sum(pars[[i]]@est)
        }
        S_X2[df.S_X2 <= 0] <- NaN
        ret$S_X2 <- na.omit(S_X2)
        ret$df.S_X2 <- na.omit(df.S_X2)
        ret$p.S_X2 <- 1 - suppressWarnings(pchisq(ret$S_X2, ret$df.S_X2))
    }
    ret[,sapply(ret, class) == 'numeric'] <- round(ret[,sapply(ret, class) == 'numeric'], digits)
    return(ret)
}
