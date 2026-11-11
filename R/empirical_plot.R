#' Function to generate empirical unidimensional item and test plots
#'
#' Given a dataset containing item responses this function will construct empirical graphics
#' using the observed responses, potentially conditioned on the (reduced)
#' total/composite score. When individual
#' item plots are requested then the total score will be formed without the item of interest
#' (i.e., the total score without that item).
#'
#' Note that some of these plot types should only be used for unidimensional
#' tests with monotonically increasing item
#' response functions. If monotonicity is not true for all items, however, then these plots may
#' serve as a visual diagnostic tool so long as the majority of items are indeed monotonic.
#'
#' @aliases empirical_plot
#' @param data a \code{data.frame} or \code{matrix} of item responses (see \code{\link{mirt}}
#'   for typical input)
#' @param group (optional) grouping variable to be paired with \code{data} input to create comparisons
#'   across discrete population subsets
#' @param org.data identical to \code{data}, but contains the "unscored" response options (e.g.,
#'   the original coding in a multiple-choice test). Used in various item-level plots for
#'   diagnostic purposes, such as in distractor analyses. This will also automatically add
#'   size and linetype changes to highlight the detected scored categories
#' @param which.items a numeric vector indicating which items to plot in a faceted image plot.
#'   If NULL then empirical test plots will be constructed instead
#' @param smooth logical; include a GAM smoother instead of the raw proportions? Default is FALSE
#' @param type character vector specifying type of plot to draw, some of which change as a function of
#'   the \code{which.item} input.
#'
#'   \describe{
#'    \item{'prop' (default)}{cumulative total-score proportions. If  \code{which.items}
#'      specified then item-level conditional proportions are instead plotted
#'      against the reduced total scores}
#'    \item{'hist'}{histogram of total scores}
#'    \item{'cor'}{heatmap of correlations. \code{...} arguments are passed to \code{\link{plot_cors}}}
#'    \item{'cordiff'}{heatmap of correlations differences (limited to cases when \code{group} has exactly
#'      two levels). \code{...} arguments are passed to \code{\link{plot_cors}}}
#'    \item{'bar'}{barplot of total scores}
#'    \item{'freq'}{item response frequencies (supports \code{which.item})}
#'    \item{'discrim'}{reduced item-total correlations to visualize discrimination effects}
#'    \item{'difficulty'}{mean/proportion of each item}
#'    \item{'discrim_diff'}{reduced item-total correlation against item difficulty}
#'    \item{'alpha_rm'}{effect on coefficient alpha if item were removed}
#'    \item{'boxplot'}{conditional boxplots of reduced total scores (supports \code{which.items})}
#'    \item{'bubble'}{bivariate frequency  bubble plots (requires that \code{which.items}
#'      is exactly of length two)}
#'   }
#'
#' @param sort logical; when applicable, sort the items first (e.g., in discrimination plot)?
#' @param formula formula used for the GAM smoother
#' @param main the main title for the plot. If NULL an internal default will be used
#' @param auto.key plotting argument passed to \code{\link[lattice]{lattice}}
#' @param par.strip.text plotting argument passed to \code{\link[lattice]{lattice}}
#' @param par.settings plotting argument passed to \code{\link[lattice]{lattice}}
#' @param discrim.cut horizontal cut-off line to use when \code{type = 'discrim'}. Default
#'   is .2 (to omit, use \code{NA})
#' @param ts_fun similar to the same argument in \code{\link{itemstats}}
#' @param ... additional arguments to be passed to \code{\link[lattice]{lattice}} and \code{coef()}
#' @keywords empirical plots
#' @export empirical_plot
#' @references
#' Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @seealso \code{\link{itemstats}}, \code{\link{itemplot}}, \code{\link{itemGAM}}, \code{\link{plot_cors}}
#' @examples
#'
#' \donttest{
#'
#' SAT12[SAT12 == 8] <- NA
#' data <- key2binary(SAT12,
#'    key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#'
#' # test plot
#' empirical_plot(data)
#' empirical_plot(data, type = 'cor')
#' empirical_plot(data, type = 'hist')
#' empirical_plot(data, type = 'hist', ts_fun=rowMeans)
#' empirical_plot(data, type = 'hist', breaks=20)
#' empirical_plot(data, type = 'bar')
#' empirical_plot(data, type = 'discrim')
#' empirical_plot(data, type = 'discrim', sort=TRUE)
#' empirical_plot(data, type = 'difficulty')
#' empirical_plot(data, type = 'difficulty', sort=TRUE)
#' empirical_plot(data, type = 'discrim_diff')
#' empirical_plot(data, type = 'alpha_rm')
#' empirical_plot(data, type = 'freq')
#'
#'
#' # items 1, 2 and 5
#' empirical_plot(data, c(1, 2, 5), type = 'freq')
#' empirical_plot(data, c(1, 2, 5))
#' empirical_plot(data, c(1, 2, 5), smooth = TRUE)
#' empirical_plot(data, c(1, 2, 5), type = 'boxplot')
#' empirical_plot(data, c(1, 2), type = 'bubble')
#' empirical_plot(data, c(1, 5), type = 'bubble')
#'
#' # replace weird looking items with unscored versions for diagnostics
#' empirical_plot(data, 32)
#' data2 <- data
#' data2[,32] <- SAT12[,32]
#' empirical_plot(data2, 32)
#' empirical_plot(data2, 32, smooth = TRUE)
#'
#' # alternatively, distractor analyses using original dataset
#' empirical_plot(data, which.items=32, org.data=SAT12)
#' empirical_plot(data, which.items=32, org.data=SAT12, smooth=TRUE)
#' empirical_plot(data, which.items=1:12, org.data=SAT12, smooth=TRUE)
#' empirical_plot(data, which.items=13:32, org.data=SAT12, smooth=TRUE)
#'
#' # example where item 32 has two correct scoring keys, and item 31
#' # is incorrectly scored (perhaps miskeyed)
#' key <- cbind(c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,1,5),
#'              c(rep(NA, 31), 3))
#' data2 <- key2binary(SAT12, key=key)
#'
#' # item 30 is correct, but 31 proportions do not climb. Item 32 has two
#' # highlighted correct answers
#' empirical_plot(data2, which.items=30:32, org.data=SAT12)
#' empirical_plot(data2, which.items=30:32, org.data=SAT12, smooth=TRUE)
#'
#'
#' #################
#' # polytomous response data
#' empirical_plot(Science)
#' empirical_plot(Science, type = 'cor')
#' empirical_plot(cbind(Science,
#'    Science + sample(-1:1, prod(dim(Science)), replace=TRUE)), type = 'cor')
#' empirical_plot(Science, type = 'hist', breaks=20)
#' empirical_plot(Science, type = 'freq')
#' empirical_plot(Science, type = 'difficulty')
#' empirical_plot(Science, type = 'discrim_diff')
#' empirical_plot(Science, type = 'boxplot')
#'
#' # item-level
#' empirical_plot(Science, c(1, 2), type = 'bubble')
#' empirical_plot(Science, c(1, 3), type = 'bubble')
#' empirical_plot(Science, which.items = 1:4, type = 'prop')
#' empirical_plot(Science, which.items = 1:4, type = 'prop', smooth=TRUE)
#'
#' # last plot very similar to model-based approach (though conditioned
#' #   on reduced total scores rather than scaled latent trait)
#' mod <- mirt(Science)
#' plot(mod, type='trace')
#'
#' # when missing values present
#' Science[1:3, 1] <- NA
#' Science[6:8, 2] <- NA
#' empirical_plot(Science, type = 'freq')
#'
#'
#' #####################
#' # Using the group input
#'
#' # single factor
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('2PL', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' data <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#'
#' # test plot
#' empirical_plot(data, group=group)
#' empirical_plot(data, type = 'cor', group=group)
#'
#' # this is likely more useful
#' empirical_plot(data, type = 'cordiff', group=group)
#'
#' empirical_plot(data, type = 'hist', group=group)
#' empirical_plot(data, type = 'hist', ts_fun=rowMeans, group=group)
#' empirical_plot(data, type = 'hist', breaks=20, group=group)
#' empirical_plot(data, type = 'bar', group=group)
#' empirical_plot(data, type = 'discrim', group=group)
#' empirical_plot(data, type = 'discrim', sort=TRUE, group=group)
#' empirical_plot(data, type = 'difficulty', group=group)
#' empirical_plot(data, type = 'difficulty', sort=TRUE, group=group)
#' empirical_plot(data, type = 'discrim_diff', group=group)
#' empirical_plot(data, type = 'alpha_rm', group=group)
#' empirical_plot(data, type = 'freq', group=group)
#'
#' }
empirical_plot <- function(data, which.items = NULL, type = 'prop', group=NULL,
                           smooth = FALSE, sort=FALSE, formula = resp ~ s(TS, k = 5),
                           org.data = NULL, ts_fun = rowSums,
                           discrim.cut = .2, main = NULL, par.strip.text = list(cex = 0.7),
                           par.settings = list(strip.background = list(col = '#9ECAE1'),
                                               strip.border = list(col = "black")),
                           auto.key = list(space = 'right', points=FALSE, lines=TRUE), ...){
    stopifnot(is.matrix(data) || is.data.frame(data))
    if(is.null(which.items) && type %in% c("freq", 'boxplot', 'cor', 'cordiff'))
        which.items <- 1:ncol(data)
    stopifnot(type %in% c('prop', 'hist', 'boxplot', 'bubble', 'alpha_rm', 'bar', 'cor', 'cordiff',
                          'discrim', 'difficulty', 'discrim_diff', 'freq'))
    if(is.null(which.items))
        stopifnot("Must specify which.items"=type %in%
                      c('prop', 'hist', 'alpha_rm', 'discrim', 'difficulty', 'discrim_diff', 'bar'))
    if(type == 'cor' || type == 'cordiff'){
        if(!is.null(group)){
            ugroups <- unique(group)
            ngroups <- length(unique(group))
            ret <- vector('list', ngroups)
            names(ret) <- ugroups
            if(type == 'cordiff'){
                stopifnot("cordiff only supported when there are two groups in group"=ngroups == 2)
                R1 <- cor(subset(data, group == ugroups[1], select=which.items), use = 'pairwise.complete')
                R2 <- cor(subset(data, group == ugroups[2], select=which.items), use = 'pairwise.complete')
                ret <- plot_cors(R1 - R2,
                                 main = sprintf('Difference between %s and %s',
                                                ugroups[1], ugroups[2]), ...)
            } else {
                for(i in 1:ngroups)
                    ret[[i]] <- plot_cors(subset(data, group == ugroups[i],
                                                select=which.items),
                                          main = paste0('group: ', ugroups[i]), ...)
            }
            return(ret)
        }
        return(plot_cors(subset(data, select=which.items), ...))
    }
    if(type %in% c('boxplot', 'bubble')) smooth <- FALSE
    if(!(type %in% c('freq')))
        data <- na.omit(as.matrix(data))
    if(is.null(org.data))
        org.data <- data
    if(!(type %in% c('freq')))
        org.data <- na.omit(as.matrix(org.data))
    stopifnot("dimensions of org.data do not match data" =
                  all(dim(data) == dim(org.data)))
    K <- apply(org.data, 2, function(x) length(unique(x)))
    if(all(K == 2L)) auto.key <- FALSE
    key <- NULL
    if(!identical(data, org.data)){
        key <- lapply(1:ncol(data), \(i){
            na.omit(unique(org.data[data[,i] == 1, i]))
        })
        names(key) <- colnames(data)
        key <- key[which.items]
    }
    TS <- ts_fun(data)
    if(is.null(group)){
        ord <- order(TS)
        data <- data[ord,]
        org.data <- org.data[ord,]
        TS <- TS[ord]
    }
    tab <- table(TS)
    if(type %in% c('discrim', 'difficulty', 'discrim_diff', 'freq', 'alpha_rm')){
        isummary <- itemstats(data, ts_fun=ts_fun, itemfreq = 'proportions', group=group)
        if(is.null(group)){
            is <- isummary$itemstats
            is$item <- factor(rownames(is), levels=colnames(data))
            if(type == 'discrim'){
                if(sort){
                    is <- is[order(is$cor_if_rm),]
                    is$item <- factor(as.character(is$item), levels=as.character(is$item))
                }
                plt <- lattice::xyplot(cor_if_rm ~ item, is,
                                       pch = 16,
                                       panel = function(x, y, ...) {
                                           panel.xyplot(x, y, type = 'b', lty = 'dashed', ...)
                                           panel.abline(h = discrim.cut, col='red', lty=2)
                                           panel.abline(h = 0, col='black', lty=1)
                                       },
                                       main = if(is.null(main)) "Reduced Item-Total Correlation" else main,
                                       xlab = 'Item', ylab='Correlation',
                                       scales = list(x = list(rot = 90)), ...)
                return(plt)
            }
            if(type == 'difficulty'){
                if(sort){
                    is <- is[order(is$mean),]
                    is$item <- factor(as.character(is$item), levels=as.character(is$item))
                }
                plt <- lattice::xyplot(mean ~ item, is, type = 'b', lty = 'dashed',
                                       pch = 16,
                                       main = if(is.null(main)) "Item Difficulty" else main,
                                       xlab = 'Item', ylab='Mean',
                                       scales = list(x = list(rot = 90)), ...)
                return(plt)
            }
            if(type == 'discrim_diff'){
                plt <- lattice::xyplot(cor_if_rm ~ mean, is,
                                       pch = 16,
                                       panel = function(x, y, subscripts, ...) {
                                           panel.xyplot(x, y, ...)
                                           panel.text(x, y, labels = is$item[subscripts],
                                                      pos = 3, offset = 0.8, cex = 0.75, col = "black")
                                           panel.abline(h = discrim.cut, col='red', lty=2)
                                           panel.abline(h = 0, col='black', lty=1)
                                       },
                                       main = if(is.null(main)) "Difficulty by Discrimination" else main,
                                       xlab = 'Mean', ylab='Reduced item-total correlation', ...)
                return(plt)
            }
            if(type == 'freq'){
                Freq <- cbind(item=is$item, isummary$proportions * isummary$overall$N)[which.items,]
                mlt <- reshape(Freq, varying=list(which(colnames(Freq) != 'item')),
                               direction='long', v.names='freq', timevar = 'cat',
                               times=colnames(isummary$proportions))
                mlt$cat[is.na(mlt$cat)] <- '<NA>'
                mlt$cat <- factor(mlt$cat)
                col <- rep('#9ECAE1', length(levels(mlt$cat)))
                if(any(mlt$cat == '<NA>')) col[length(col)] <- 'red'
                plt <- lattice::barchart(freq ~ cat|item, mlt, horizontal = FALSE,
                                         ylab = 'Frequency', xlab = 'Category', col = col,
                                  main = if(is.null(main)) "Item Response Frequency" else main, ...)
                return(plt)
            }
            if(type == 'alpha_rm'){
                if(sort){
                    is <- is[order(is$alpha_if_rm),]
                    is$item <- factor(as.character(is$item), levels=as.character(is$item))
                }
                plt <- lattice::xyplot(alpha_if_rm ~ item, is,
                                       pch = 16,
                                       panel = function(x, y, ...) {
                                           panel.xyplot(x, y, type = 'b', lty = 'dashed', ...)
                                           panel.abline(h = isummary$overall$alpha, lty=2)
                                       },
                                       main = if(is.null(main)) "Alpha if Item Removed" else main,
                                       xlab = 'Item', ylab=expression(alpha),
                                       ylim = c(min(is$alpha_if_rm) - .05, max(is$alpha_if_rm) + .05),
                                       scales = list(x = list(rot = 90)), ...)
                return(plt)
            }
        } else { # multiple-group
            is <- lapply(isummary, \(x) x$itemstats)
            is <- lapply(is, \(x){
                x$item <- factor(rownames(x), levels=colnames(data))
                x
            })
            if(type == 'discrim'){
                if(sort){
                    ord <- order(is[[1]]$cor_if_rm)
                    lvls <- as.character(is[[1]]$item[ord])
                    is <- lapply(is, \(x){
                        x <- x[ord, ]
                        x$item <- factor(as.character(x$item), levels=lvls)
                        x
                    })
                }
                grouplvl <- rep(names(is), each=nrow(is[[1]]))
                isbind <- do.call(rbind, is)
                isbind$group <- grouplvl
                plt <- lattice::xyplot(cor_if_rm ~ item, isbind, group=group,
                                       pch = 16,
                                       panel = function(x, y, ...) {
                                           panel.xyplot(x, y, type = 'b', lty = 'dashed', ...)
                                           panel.abline(h = discrim.cut, col='red', lty=2)
                                           panel.abline(h = 0, col='black', lty=1)
                                       },
                                       main = if(is.null(main)) "Reduced Item-Total Correlation" else main,
                                       xlab = 'Item', ylab='Correlation',
                                       scales = list(x = list(rot = 90)),
                                       auto.key=TRUE, ...)
                return(plt)
            }
            if(type == 'difficulty'){
                if(sort){
                    ord <- order(is[[1]]$mean)
                    lvls <- as.character(is[[1]]$item[ord])
                    is <- lapply(is, \(x){
                        x <- x[ord, ]
                        x$item <- factor(as.character(x$item), levels=lvls)
                        x
                    })
                }
                grouplvl <- rep(names(is), each=nrow(is[[1]]))
                isbind <- do.call(rbind, is)
                isbind$group <- grouplvl
                plt <- lattice::xyplot(mean ~ item, isbind, type = 'b', lty = 'dashed',
                                       pch = 16, group=group,
                                       main = if(is.null(main)) "Item Difficulty" else main,
                                       xlab = 'Item', ylab='Mean',
                                       scales = list(x = list(rot = 90)),
                                       auto.key=TRUE, ...)
                return(plt)
            }
            if(type == 'discrim_diff'){
                grouplvl <- rep(names(is), each=nrow(is[[1]]))
                isbind <- do.call(rbind, is)
                isbind$group <- grouplvl
                plt <- lattice::xyplot(cor_if_rm ~ mean|group, isbind,
                                       pch = 16, group=isbind$group,
                                       panel = function(x, y, subscripts, groups, ...) {
                                           panel.xyplot(x, y, subscripts=subscripts, group=groups, ...)
                                           panel.text(x, y, labels = isbind$item[subscripts],
                                                      pos = 3, offset = 0.8, cex = 0.75, col = "black")
                                           panel.abline(h = discrim.cut, col='red', lty=2)
                                           panel.abline(h = 0, col='black', lty=1)
                                       },
                                       main = if(is.null(main)) "Difficulty by Discrimination" else main,
                                       xlab = 'Mean', ylab='Reduced item-total correlation',
                                       auto.key=TRUE, ...)
                return(plt)
            }
            if(type == 'alpha_rm'){
                if(sort){
                    ord <- order(is[[1]]$alpha_if_rm)
                    lvls <- as.character(is[[1]]$item[ord])
                    is <- lapply(is, \(x){
                        x <- x[ord, ]
                        x$item <- factor(as.character(x$item), levels=lvls)
                        x
                    })
                }
                grouplvl <- rep(names(is), each=nrow(is[[1]]))
                isbind <- do.call(rbind, is)
                isbind$group <- grouplvl
                tpar <- trellis.par.get()
                plt <- lattice::xyplot(alpha_if_rm ~ item, isbind, group=group,
                                       pch = 16,
                                       panel = function(x, y, groups, ...) {
                                           panel.xyplot(x, y, type = 'b', lty = 'dashed', group=groups, ...)
                                           alphas <- sapply(isummary, \(x) x$overall$alpha)
                                           panel.abline(h = alphas,
                                                        col=tpar$superpose.symbol$col[1:length(unique(group))], lty=2)
                                       },
                                       main = if(is.null(main)) "Alpha if Item Removed" else main,
                                       xlab = 'Item', ylab=expression(alpha),
                                       ylim = c(min(isbind$alpha_if_rm) - .05, max(isbind$alpha_if_rm) + .05),
                                       scales = list(x = list(rot = 90)), ...)
                return(plt)
            }
            if(type == 'freq'){
                ngroups <- length(unique(group))
                Freq <- do.call(rbind, lapply(1:ngroups, \(i){
                    cbind(item=is[[i]]$item,
                                  isummary[[i]]$proportions * isummary[[i]]$overall$N,
                          group=names(is)[i])[which.items,]
                }))
                mlt <- reshape(Freq, varying=list(which(colnames(Freq) != c('item', 'group'))),
                               direction='long', v.names='freq', timevar = 'cat',
                               times=colnames(isummary[[1]]$proportions))
                mlt$cat[is.na(mlt$cat)] <- '<NA>'
                mlt$cat <- factor(mlt$cat)
                plt <- lattice::barchart(freq ~ cat|item, mlt, horizontal = FALSE, group=group,
                                         ylab = 'Frequency', xlab = 'Category',
                                         main = if(is.null(main)) "Item Response Frequency" else main,
                                         auto.key=TRUE, ...)
                return(plt)
            }

        }
    }
    if(type == 'bubble'){
        stopifnot("which.items must be of length two"=length(which.items) == 2)
        dat.sub <- as.data.frame(data[,which.items, drop=FALSE])
        colnames(dat.sub) <- c('x', 'y')
        tab <- as.data.frame(table(dat.sub))
        cfs <- coef(lm(y ~ x, dat.sub))
        if(min(dat.sub$y) == 0) cfs[1] <- cfs[1] + 1
        plt <- lattice::xyplot(y ~ x, tab,
                               panel = function(x, y, subscripts, cex, cfs, ...) {
                                   panel.xyplot(x, y, cex = cex[subscripts], ...)
                                   panel.abline(coef=cfs, col='red', lty=2, lwd=2)
                               },
                               pch = 16, cfs=cfs,
                               cex=sqrt(tab$Freq) / max(sqrt(tab$Freq)) * 3,
                               main = if(is.null(main))
                                   paste0('Correlation = ', round(cor(dat.sub)[1,2], 2)) else main,
                               xlab = colnames(data)[which.items[1]],
                               ylab = colnames(data)[which.items[2]])
        return(plt)
    }
    if(is.null(which.items)){
        if(type == 'prop'){
            prop <- cumsum(tab) / nrow(data)
            df <- data.frame(TS=as.integer(names(tab)), P=prop)
            plt <- lattice::xyplot(P ~ TS, df, type = 'b',
                                   main = if(is.null(main)) 'Cumulative Composite Score' else main,
                                   xlab = 'Composite Score', ylab = 'Cumulative Proportion',
                                   ylim = c(-.1, 1.1), ...)
        } else if(type == 'hist'){
            plt <- lattice::histogram(TS,
                                      main = if(is.null(main)) 'Composite Scores' else main,
                                      xlab = 'Composite Score', ylab = 'Frequency', type = 'count', ...)
        } else if(type == 'bar'){
            plt <- lattice::barchart(table(TS), horizontal = FALSE,
                                     main = if(is.null(main)) 'Composite Scores' else main,
                                     xlab = 'Composite Score', ylab = 'Frequency', ...)
        }
    } else {
        stopifnot(all(which.items >= 1L & which.items <= ncol(data)))
        pltdat <- vector('list', length(which.items))
        nms <- colnames(data)
        for(i in 1:length(which.items)){
            item <- org.data[, which.items[i]]
            TS <- ts_fun(data[ ,-which.items[i]])
            ord <- order(TS)
            item <- item[ord]
            TS <- TS[ord]
            if(smooth){
                uniq <- sort(unique(item))
                if(length(uniq) == 2L) uniq <- uniq[2L]
                splt <- as.list(uniq)
                for(j in 1:length(uniq)){
                    df <- data.frame(resp = item==uniq[j], TS=TS)
                    props <- fitted(gam(formula, df, family = binomial()))
                    splt[[j]] <- data.frame(item=nms[which.items[i]], TS=TS,
                                            props=as.numeric(props), cat=uniq[j])
                }
                pltdat[[i]] <- do.call(rbind, splt)
            } else {
                tab <- table(TS)
                tmptab <- tab
                splt <- split(TS, item)
                if(type != "boxplot")
                    if(length(splt) == 2L) splt[[1L]] <- NULL
                for(j in 1:length(splt)){
                    tmptab[] <- NA
                    tab2 <- table(splt[[j]])
                    tmptab[names(tab) %in% names(tab2)] <- unname(tab2)
                    props <- tmptab / tab
                    names(props) <- names(tab)
                    splt[[j]] <- data.frame(item=nms[which.items[i]], TS=as.integer(names(tab)),
                                            props=as.numeric(props), freq=as.numeric(tab),
                                            cat=names(splt)[j])
                }
                pltdat[[i]] <- do.call(rbind, splt)
            }
        }
        df <- na.omit(do.call(rbind, pltdat))
        df$cat <- factor(df$cat)
        df$item <- factor(df$item, levels=colnames(data)[which.items])
        if(type == "boxplot"){
            plt <- lattice::bwplot(TS ~ cat | item, df,
                                   main = if(is.null(main)) "Item Category by Composite" else main,
                                   xlab = 'Item Category', ylab = 'Reduced Composite Score',
                                   par.strip.text=par.strip.text, par.settings=par.settings,
                                   auto.key=auto.key, ...)
        } else if(type == 'prop'){
            plt <- lattice::xyplot(props ~ TS|item, df, groups = cat,
                                   type = ifelse(smooth, 'l', 'b'),
                                   panel = function(x, y, groups, ...) {
                                       current_panel <- panel.number()
                                       lwd <- lty <- rep(1, length(unique(groups)))
                                       if(is.list(key)){
                                           lty[] <- 2
                                           item <- levels(df$item)[current_panel]
                                           pick <- levels(groups) %in% key[[item]]
                                           lwd[pick] <- 2
                                           lty[pick] <- 1
                                       }
                                       panel.xyplot(x, y, groups=groups, lwd=lwd, lty=lty, ...)
                                   },
                                   main = if(is.null(main)) "Item-Total Plot" else main,
                                   xlab = 'Reduced Composite Score', ylab = 'Proportion',
                                   par.strip.text=par.strip.text, par.settings=par.settings,
                                   auto.key=auto.key, ylim = c(-.1, 1.1), ...)
        }
    }
    plt
}


#' Correlation matrix plot
#'
#' Heatmap version of correlation matrix to visualize bivariate correlations.
#'
#' @param r square matrix (correlations, or any matrix with values on a common
#'              scale), or a data.frame/non-square matrix of raw data, in which
#'              case \code{cor(r, use = "pairwise.complete.obs")} is computed
#' @param numbers    print the (rounded) values in the cells
#' @param colors     if FALSE, a grey scale is used
#' @param n          number of colour bins
#' @param main       title
#' @param zlim       range of the colour scale (values are clamped to it)
#' @param show.legend draw the colour key
#' @param labels     variable names (default: dimnames of r)
#' @param select     indices/names of the variables to show
#' @param pval       optional matrix of p-values (same dim as r). If supplied, the
#'              upper triangle shows the p-values as plain text with no colour
#'              fill, and the lower triangle/diagonal show the coloured
#'              correlations
#' @param cuts       p-value cutpoints for stars, e.g. c(.001, .01) -> *** / ** / *
#'              (a final single * is used for p < .05)
#' @param stars      optionally draw significance stars after the lower-triangle
#'             correlations (requires pval; FALSE by default, since the
#'              p-values themselves are shown in the upper triangle)
#' @param pdigits    decimals for printed p-values (smaller values print as <.001)
#' @param scale      if TRUE, text is drawn at size cex for every cell (diagonal
#'              included) whose |r| >= thresh, and shrinks linearly with |r|
#'              below thresh
#' @param thresh     \code{|r|} below which the text starts to shrink
#' @param min.cex    smallest relative text size (reached as |r| -> 0)
#' @param fit        if \code{TRUE}, text is shrunk (never enlarged) to fit the drawn cells;
#'              the size is measured at draw time, so it adapts to the device.
#' @param cex   A common shrink factor is set from the bulk of the labels (not the
#'              single widest one), then any label that still overflows its own
#'              cell is capped individually
#' @param fill       proportion of a cell width/height the typical label may fill
#' @param drop.zero  drop the leading zero of the printed values (.45, -.45)
#' @param upper      draw the upper triangle?
#' @param diag       draw the diagonal?
#' @param symmetric  if \code{TRUE} and \code{r} is square, row i is plotted in the same position
#'              as column i
#' @param  xaxis       1 = x labels at the bottom, 3 = at the top
#' @param  xlas       rotation (degrees) of x labels
#' @param  gr         colour-ramp function, e.g., \code{colorRampPalette(c("red","white","blue"))}
#' @param  digits     digits shown in cells
#' @param  alpha      transparency of the cell colours
#' @param  border     colour of cell borders (NA/"transparent" for none)
#' @param  ...        passed on to \code{lattice::levelplot()}
#' @export
#'
plot_cors <- function(r, numbers = TRUE, colors = TRUE, n = 51,
                            main = NULL, zlim = c(-1, 1), show.legend = TRUE,
                            labels = NULL, select = NULL, pval = NULL,
                            cuts = c(.001, .01), stars = FALSE, pdigits = 3, scale = TRUE,
                            cex = 1, thresh = .2, min.cex = .6, drop.zero = TRUE,
                            fit = TRUE, fill = .9,
                            upper = TRUE, diag = TRUE,
                            symmetric = TRUE, xaxis = 1, xlas = 90,
                            gr = NULL, digits = 2, alpha = .75,
                            border = "white", ...) {
    stopifnot(requireNamespace("lattice", quietly = TRUE))
    if (is.data.frame(r)) r <- as.matrix(r)
    if (!is.matrix(r)) stop("'r' must be a matrix or data.frame")
    if (nrow(r) != ncol(r)) r <- stats::cor(r, use = "pairwise.complete.obs")

    ## names
    if (is.null(labels)) labels <- colnames(r)
    if (is.null(labels)) labels <- paste0("V", seq_len(ncol(r)))
    if (is.null(rownames(r))) rownames(r) <- labels
    colnames(r) <- labels
    if (!is.null(pval)) { dimnames(pval) <- dimnames(r) }

    ## subset
    if (!is.null(select)) {
        r <- r[select, select, drop = FALSE]
        if (!is.null(pval)) pval <- pval[select, select, drop = FALSE]
    }
    p <- ncol(r)
    nms <- colnames(r)

    ## mask triangles
    keep <- matrix(TRUE, p, p)
    if (!upper) keep[upper.tri(keep)] <- FALSE
    if (!diag)  diag(keep) <- FALSE

    ## long format; row 1 plotted at the top
    df <- expand.grid(row = seq_len(p), col = seq_len(p))
    df$z <- as.vector(r)               # row varies fastest = column-major, matching as.vector()
    df$z[!as.vector(keep)] <- NA
    ## when p-values are supplied the upper triangle holds p-values: those
    ## cells get no fill, so they are excluded from the colour variable
    has.p <- !is.null(pval)
    is.up <- df$row < df$col
    is.lo <- df$row > df$col
    if (has.p) df$z[is.up] <- NA
    df$zc <- pmin(pmax(df$z, zlim[1]), zlim[2])   # clamp to colour range

    ## p-value matrix, symmetrised (use the other triangle if one is missing)
    if (has.p) {
        pm <- pval
        pm[is.na(pm)] <- t(pm)[is.na(pm)]
        pv <- as.vector(pm)
    }
    if (stars && !has.p) stop("'stars = TRUE' requires 'pval'")

    ## cell text
    txt <- rep("", nrow(df))
    if (numbers) {
        txt <- formatC(df$z, digits = digits, format = "f")
        ## unit diagonal printed as "1" rather than "1.00" to save space
        is.unit <- df$row == df$col & !is.na(df$z) & abs(df$z - 1) < 1e-8
        txt[is.unit] <- "1"
        ## drop the leading zero (0.45 -> .45, -0.45 -> -.45) to save space
        if (drop.zero) txt <- sub("^(-?)0\\.", "\\1.", txt)
    }
    txt[is.na(df$z)] <- ""
    if (stars) {                        # stars go on the lower-triangle correlations
        cuts <- sort(cuts)
        st <- ifelse(pv < cuts[1], "***",
                     ifelse(pv < cuts[2], "**",
                            ifelse(pv < .05,     "*", "")))
        st[!is.lo | is.na(st)] <- ""
        txt <- paste0(txt, st)
    }
    if (has.p) {                        # p-values in the upper triangle
        up <- is.up & keep[cbind(df$row, df$col)] & !is.na(pv)
        ptxt <- ifelse(pv < 10^-pdigits, paste0("<.", strrep("0", pdigits - 1), "1"),
                       formatC(pv, digits = pdigits, format = "f"))
        if (drop.zero) ptxt <- sub("^0\\.", ".", ptxt)
        txt[up] <- ptxt[up]
    }
    df$txt <- txt
    df$cex <- cex
    if (scale) {
        ok <- !is.na(df$zc)
        df$cex[ok] <- cex * pmax(min.cex, pmin(1, abs(df$zc[ok]) / thresh))
    }

    ## colours
    if (colors) {
        if (is.null(gr)) gr <- grDevices::colorRampPalette(c("red", "white", "blue"))
        cols <- gr(n)
    } else {
        cols <- grDevices::gray(seq(0, 1, length.out = n))
    }
    cols <- grDevices::adjustcolor(cols, alpha.f = alpha)
    brks <- seq(zlim[1], zlim[2], length.out = n + 1)

    ## legend
    ck <- if (show.legend) list(space = "right", width = 1.2, tck = .6) else FALSE

    ## axes: y runs top-to-bottom, so reverse ylim
    sc <- list(x = list(at = seq_len(p), labels = nms, rot = xlas, tck = 0,
                        alternating = if (xaxis == 3) 2 else 1),
               y = list(at = seq_len(p), labels = rownames(r), rot = 0, tck = 0,
                        alternating = 1))

    lattice::levelplot(
        zc ~ col * row, data = df,
        main = main, xlab = NULL, ylab = NULL,
        at = brks, col.regions = cols, colorkey = ck,
        scales = sc, xlim = c(.5, p + .5), ylim = c(p + .5, .5),
        aspect = "iso", border = border,
        panel = function(x, y, z, subscripts, ...) {
            lattice::panel.levelplot(x, y, z, subscripts, ...)
            if (numbers || stars || has.p) {
                s  <- subscripts[nzchar(df$txt[subscripts])]
                cx <- df$cex[s]
                if (fit && length(s)) {
                    fs <- lattice::trellis.par.get("fontsize")$text
                    wi <- grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE) / p
                    hi <- grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE) / p
                    str <- df$txt[s]
                    u   <- unique(str)
                    gp  <- grid::gpar(fontsize = fs)
                    tw1 <- vapply(u, function(z)
                        grid::convertWidth(grid::grobWidth(grid::textGrob(z, gp = gp)),
                                           "inches", valueOnly = TRUE), numeric(1))
                    th1 <- grid::convertHeight(grid::grobHeight(grid::textGrob("0", gp = gp)),
                                               "inches", valueOnly = TRUE)
                    tw  <- tw1[match(str, u)]
                    ratio <- pmin(fill * wi / tw, fill * hi / th1) / cx
                    f   <- min(1, stats::quantile(ratio, .1, names = FALSE))
                    cx  <- pmin(cx * f, 0.98 * wi / tw)
                }
                lattice::panel.text(df$col[s], df$row[s], df$txt[s],
                                    cex = cx, col = "black")
            }
        },
        ...
    )
}
