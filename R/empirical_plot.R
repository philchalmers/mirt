#' Function to generate empirical unidimensional item and test plots
#'
#' Given a dataset containing item responses this function will construct empirical graphics
#' using the observed responses, potentially conditioned on the (reduced) total score. When individual
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
#' @param ... additional arguments to be passed to \code{\link[lattice]{lattice}} and \code{coef()}
#' @keywords empirical plots
#' @export empirical_plot
#' @references
#' Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @seealso \code{\link{itemstats}}, \code{\link{itemplot}}, \code{\link{itemGAM}}
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
#' empirical_plot(data, type = 'hist')
#' empirical_plot(data, type = 'hist', breaks=20)
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
#'
#' #################
#' # polytomous response data
#' empirical_plot(Science)
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
#' }
empirical_plot <- function(data, which.items = NULL, type = 'prop',
                           smooth = FALSE, sort=FALSE, formula = resp ~ s(TS, k = 5),
                           org.data = NULL,
                           discrim.cut = .2, main = NULL, par.strip.text = list(cex = 0.7),
                           par.settings = list(strip.background = list(col = '#9ECAE1'),
                                               strip.border = list(col = "black")),
                           auto.key = list(space = 'right', points=FALSE, lines=TRUE), ...){
    stopifnot(is.matrix(data) || is.data.frame(data))
    if(is.null(which.items) && type %in% c("freq", 'boxplot'))
        which.items <- 1:ncol(data)
    stopifnot(type %in% c('prop', 'hist', 'boxplot', 'bubble', 'alpha_rm',
                          'discrim', 'difficulty', 'discrim_diff', 'freq'))
    if(is.null(which.items))
        stopifnot("Must specify which.items"=type %in%
                      c('prop', 'hist', 'alpha_rm', 'discrim', 'difficulty', 'discrim_diff'))
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
            pick <- na.omit(unique(org.data[data[,i] == 1, i]))
            which(sort(unique(org.data[,i])) == pick) })
        names(key) <- colnames(data)
        key <- key[which.items]
    }
    TS <- rowSums(data)
    ord <- order(TS)
    data <- data[ord,]
    org.data <- org.data[ord,]
    TS <- TS[ord]
    tab <- table(TS)
    if(type %in% c('discrim', 'difficulty', 'discrim_diff', 'freq', 'alpha_rm')){
        isummary <- itemstats(data)
        is <- isummary$itemstats
        is$item <- factor(rownames(is), levels=colnames(data))
        if(type == 'discrim'){
            if(sort){
                is <- is[order(is$total.r_if_rm),]
                is$item <- factor(as.character(is$item), levels=as.character(is$item))
            }
            plt <- lattice::xyplot(total.r_if_rm ~ item, is,
                                   pch = 16,
                                   panel = function(x, y, ...) {
                                       panel.xyplot(x, y, ...)
                                       panel.abline(h = discrim.cut, col='red', lty=2)
                                   },
                                   main = if(is.null(main)) "Reduced Item-total Correlation" else main,
                                   xlab = 'Item', ylab='Correlation',
                                   scales = list(x = list(rot = 90)), ...)
            return(plt)
        }
        if(type == 'difficulty'){
            if(sort){
                is <- is[order(is$mean),]
                is$item <- factor(as.character(is$item), levels=as.character(is$item))
            }
            plt <- lattice::xyplot(mean ~ item, is,
                                   pch = 16,
                                   main = if(is.null(main)) "Item Difficulty" else main,
                                   xlab = 'Item', ylab='Mean',
                                   scales = list(x = list(rot = 90)), ...)
            return(plt)
        }
        if(type == 'discrim_diff'){
            plt <- lattice::xyplot(total.r_if_rm ~ mean, is,
                                   pch = 16,
                                   panel = function(x, y, subscripts, ...) {
                                       panel.xyplot(x, y, ...)
                                       panel.text(x, y, labels = is$item[subscripts],
                                                  pos = 3, offset = 0.8, cex = 0.75, col = "black")
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
                                       panel.xyplot(x, y, ...)
                                       panel.abline(h = isummary$overall$alpha, col='red', lty=2)
                                   },
                                   main = if(is.null(main)) "Alpha if Item Removed" else main,
                                   xlab = 'Item', ylab=expression(alpha),
                                   ylim = c(min(is$alpha_if_rm) - .05, max(is$alpha_if_rm) + .05),
                                   scales = list(x = list(rot = 90)), ...)
            return(plt)
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
                               xlab = paste0('Item ', which.items[1]),
                               ylab=paste0('Item ', which.items[2]))
        return(plt)
    }
    if(is.null(which.items)){
        if(type == 'prop'){
            prop <- cumsum(tab) / nrow(data)
            df <- data.frame(TS=as.integer(names(tab)), P=prop)
            plt <- lattice::xyplot(P ~ TS, df, type = 'b',
                                   main = if(is.null(main)) 'Cumulative Total Score' else main,
                                   xlab = 'Total Score', ylab = 'Cumulative Proportion',
                                   ylim = c(-.1, 1.1), ...)
        } else if(type == 'hist'){
            df <- data.frame(TS=as.integer(names(tab)), freq=as.integer(tab))
            plt <- lattice::histogram(TS,
                                      main = if(is.null(main)) 'Total Scores' else main,
                                      xlab = 'Total Score', ylab = 'Frequency', ...)
        }
    } else {
        stopifnot(all(which.items >= 1L & which.items <= ncol(data)))
        pltdat <- vector('list', length(which.items))
        nms <- colnames(data)
        for(i in 1:length(which.items)){
            item <- org.data[, which.items[i]]
            TS <- rowSums(data[ ,-which.items[i]])
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
                                   xlab = 'Item Category', ylab = 'Reduced Total Score',
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
                                           lwd[key[[current_panel]]] <- 2
                                           lty[key[[current_panel]]] <- 1
                                       }
                                       panel.xyplot(x, y, groups=groups, lwd=lwd, lty=lty, ...)
                                   },
                                   main = if(is.null(main)) "Item-total Plot" else main,
                                   xlab = 'Reduced Total Score', ylab = 'Proportion',
                                   par.strip.text=par.strip.text, par.settings=par.settings,
                                   auto.key=auto.key, ylim = c(-.1, 1.1), ...)
        }
    }
    plt
}
