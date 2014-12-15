#' Differential item functioning tests
#'
#' This function runs the Wald and likelihood-ratio approaches for testing differential
#' item functioning (DIF). This is primarily a convenience wrapper to the
#' \code{\link{multipleGroup}} function for performing standard DIF procedures. Models can be
#' estimated in parallel automatically by defining a parallel object with \code{\link{mirtCluster}},
#' which will help to decrease the runtime. For best results, the baseline model should contain
#' a set of 'anchor' items and have freely estimated hyper-parameters in the focal groups.
#'
#' Generally, the precomputed baseline model should have been
#' configured with two estimation properties: 1) a set of 'anchor' items exist,
#' where the anchor items have various parameters that have been constrainted to be equal
#' across the groups, and 2) contain freely estimated latent mean and variance terms in
#' all but one group. These two properties help to fix the metric of the groups so that
#' item parameter estimates do not contain latent distribution characteristics.
#'
#' @aliases DIF
#' @param MGmodel an object returned from \code{\link{multipleGroup}} to be used as the reference
#'   model
#' @param which.par a character vector containing the parameter names which will be inspected for
#'   DIF
#' @param Wald logical; perform Wald tests for DIF instead of likelihood ratio test?
#' @param items2test a numeric vector, or character vector containing the item names, indicating
#'   which items will be tested for DIF. In models where anchor items are known, omit them from
#'   this vector. For example, if items 1 and 2 are anchors in a 10 item test, then
#'   \code{items2test = 3:10} would work for testing the remaining items (important to remember
#'   when using sequential schemes)
#' @param return_models logical; return estimated model objects for further analysis?
#'   Default is FALSE
#' @param scheme type of DIF analysis to perform, either by adding or dropping constraints across
#'   groups. These can be:
#' \describe{
#'   \item{'add'}{parameters in \code{which.par} will be constrained each item one at a time for
#'     items that are specified in \code{items2test}. This is beneficial when examining DIF from a
#'     model with parameters freely estimated across groups, and when inspecting differences via
#'     the Wald test}
#'   \item{'drop'}{parameters in \code{which.par} will be freely estimated for items that are
#'     specified in \code{items2test}. This is useful when supplying an overly restrictive model
#'     and attempting to detect DIF with a slightly less restrictive model}
#'   \item{'add_sequential'}{sequentially loop over the items being tested, and at the end of the
#'     loop treat DIF tests that satisfy the \code{seq_stat} criteria as invariant. The loop is
#'     then re-run on the remaining invariant items to determine if they are now displaying DIF in
#'     the less constrained model, and when no new invariant item is found the algorithm stops and
#'     returns the items that displayed DIF}
#'   \item{'drop_sequential'}{sequentially loop over the items being tested, and at the end of the
#'     loop treat items that violate the \code{seq_stat} criteria as demonstrating DIF. The loop is
#'     then re-run, leaving the items that previously demonstrated DIF as variable across groups,
#'     and the remaining test items that previously showed invariance are re-tested. The algorithm
#'     stops when no more items showing DIF are found and returns the items that displayed DIF}
#' }
#' @param seq_stat select a statistic to test for in the sequential schemes. Potential values are
#'   (in descending order of power) \code{'AIC'}, \code{'AICc'}, \code{'SABIC'}, and \code{'BIC'}.
#'   If a numeric value is input that ranges between 0 and 1, the 'p' value will be tested
#'   (e.g., \code{seq_stat = .05} will test for the difference of p < .05 in the add scheme,
#'   or p > .05 in the drop scheme), along with the specified \code{p.adjust} input
#' @param max_run a number indicating the maximum number of cycles to perform in sequential
#'   searches. The default is to perform search until no further DIF is found
#' @param plotdif logical; create itemplots for items that are displaying DIF according to the
#'   \code{seq_stat} criteria? Only available for 'add' type schemes
#' @param type the \code{type} of plot argument passed to \code{plot()}. Default is 'trace', though
#'   another good option is 'infotrace'. For ease of viewing, the \code{facet_item} argument to
#'   mirt's \code{plot()} function is set to \code{TRUE}
#' @param p.adjust string to be passed to the \code{\link{p.adjust}} function to adjust p-values.
#'   Adjustments are located in the \code{adj_pvals} element in the returned list
#' @param verbose logical print extra information to the console?
#' @param ... additional arguments to be passed to \code{\link{multipleGroup}} and \code{plot}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords DIF
#' @seealso \code{\link{multipleGroup}}, \code{\link{DTF}}
#' @export DIF
#' @examples
#' \dontrun{
#'
#' #simulate data where group 2 has a smaller slopes and more extreme intercepts
#' set.seed(12345)
#' a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
#' a2[1:2, ] <- a1[1:2, ]/3
#' d1[c(1,3), ] <- d2[c(1,3), ]/4
#' head(data.frame(a.group1 = a1, a.group2 = a2, d.group1 = d1, d.group2 = d2))
#' itemtype <- rep('dich', nrow(a1))
#' N <- 1000
#' dataset1 <- simdata(a1, d1, N, itemtype)
#' dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#'
#' #### no anchors, all items tested for DIF by adding item constrains one item at a time.
#' # define a parallel cluster (optional) to help speed up internal functions
#' mirtCluster()
#'
#' #  Information matrix with crossprod (not controlling for latent group differences)
#' model <- multipleGroup(dat, 1, group, SE = TRUE)
#'
#' #test whether adding slopes and intercepts constraints results in DIF. Plot items showing DIF
#' resulta1d <- DIF(model, c('a1', 'd'), plotdif = TRUE)
#' resulta1d
#'
#' #same as above, but using Wald tests with Benjamini & Hochberg adjustment
#' resulta1dWald <- DIF(model, c('a1', 'd'), Wald = TRUE, p.adjust = 'fdr')
#' resulta1dWald
#' round(resulta1dWald$adj_pvals, 4)
#'
#' #test whether adding only slope constraints results in DIF for all items
#' resulta1 <- DIF(model, 'a1')
#' resulta1
#'
#' #following up on resulta1d, to determine whether it's a1 or d parameter causing DIF
#' (a1s <- DIF(model, 'a1', items2test = 1:3))
#' (ds <- DIF(model, 'd', items2test = 1:3))
#'
#' #### using items 4 to 15 as anchors
#' itemnames <- colnames(dat)
#' model_anchor <- multipleGroup(dat, model = 1, group = group,
#'   invariance = c(itemnames[4:15], 'free_means', 'free_var'))
#' anchor <- DIF(model_anchor, c('a1', 'd'), items2test = 1:3)
#' anchor
#'
#' ### drop down approach (freely estimating parameters accross groups) when
#' ### specifying a highly constrained model with estimated latent parameters
#' model_constrained <- multipleGroup(dat, 1, group,
#'   invariance = c(colnames(dat), 'free_means', 'free_var'))
#' dropdown <- DIF(model_constrained, 'd', scheme = 'drop')
#' dropdown
#'
#' ### sequential searches using SABIC as the selection criteria
#' # starting from completely different models
#' model <- multipleGroup(dat, 1, group)
#' stepup <- DIF(model, c('a1', 'd'), scheme = 'add_sequential')
#' stepup
#'
#' #step down procedure for highly constrained model
#' model <- multipleGroup(dat, 1, group, invariance = itemnames)
#' stepdown <- DIF(model, c('a1', 'd'), scheme = 'drop_sequential')
#' stepdown
#' }
DIF <- function(MGmodel, which.par, scheme = 'add', items2test = 1:ncol(MGmodel@Data$data),
                seq_stat = 'SABIC', Wald = FALSE, p.adjust = 'none', return_models = FALSE,
                max_run = Inf, plotdif = FALSE, type = 'trace', verbose = TRUE, ...){

    loop_test <- function(item, model, which.par, values, Wald, itemnames, invariance, drop,
                          return_models, ...)
    {
        constrain <- model@constrain
        parnum <- list()
        for(i in 1L:length(which.par))
            parnum[[i]] <- values$parnum[values$name == which.par[i] &
                                             values$item == itemnames[item]]
        for(i in length(parnum):1L)
            if(!length(parnum[[i]])) parnum[[i]] <- NULL
        if(!length(parnum))
            stop('Item ', item, ' does not contain any of the parameters defined in which.par.
                 Consider removing it from the item2test input or adding relevant parameters
                 to which.par')
        if(Wald){
            wv <- wald(model)
            infoname <- wv[1L, ]
            L <- matrix(0, length(parnum), length(infoname))
            for(i in 1L:length(parnum)){
                L[i, paste0(which.par[i], '.', parnum[[i]][1L]) == infoname] <- 1
                L[i, paste0(which.par[i], '.', parnum[[i]][2L]) == infoname] <- -1
            }
            res <- wald(model, L)
            return(res)
        }
        if(drop){
            for(j in 1L:length(parnum)){
                for(i in length(constrain):1L){
                    if(all(parnum[[j]] == sort(constrain[[i]])))
                        constrain[[i]] <- NULL
                }
            }
        } else {
            for(i in 1L:length(parnum))
                constrain[[length(constrain) + 1L]] <- parnum[[i]]
        }
        newmodel <- multipleGroup(model@Data$data, model@model[[1L]], group=model@Data$group,
                                  invariance = invariance, constrain=constrain,
                                  verbose = FALSE, ...)
        aov <- anova(newmodel, model, verbose = FALSE)
        attr(aov, 'parnum') <- parnum
        if(return_models) aov <- newmodel
        return(aov)
    }

    bfactorlist <- MGmodel@bfactor
    if(!is.null(bfactorlist$Priorbetween[[1L]]))
        stop('bifactor models are currently not supported in this function')
    itemnames <- colnames(MGmodel@Data$data)
    if(!any(scheme %in% c('add', 'drop', 'add_sequential', 'drop_sequential')))
        stop('scheme input is not valid')
    if(return_models){
        if(Wald)
            stop('return_models argument only valid for likelihood ratio tests')
    }
    if(Wald){
        if(scheme != 'add')
            stop('Wald tests are only appropriate when add scheme is used')
        if(length(MGmodel@information) == 1)
            stop('Information matrix was not calculated')
    }
    pval <- 0
    if(is.numeric(seq_stat)){
        pval <- seq_stat
        seq_stat <- 'p'
    } else if(!any(seq_stat %in% c('p', 'AIC', 'AICc', 'SABIC', 'BIC'))){
        stop('Invalid seq_stat input')
    }
    if(is.character(items2test)) items2test <- which(items2test %in% itemnames)
    data <- MGmodel@Data$data
    invariance <- MGmodel@invariance
    values <- mod2values(MGmodel)
    drop <- scheme == 'drop' || scheme == 'drop_sequential'
    invariance <- MGmodel@invariance[MGmodel@invariance %in%
                                         c('free_means', 'free_var', 'free_varcov', 'free_cov')]
    if(!length(invariance)) invariance <- ''
    res <- myLapply(X=items2test, FUN=loop_test, model=MGmodel, which.par=which.par, values=values,
                    Wald=Wald, drop=drop, itemnames=itemnames, invariance=invariance,
                    return_models=return_models, ...)
    names(res) <- itemnames[items2test]
    if(scheme %in% c('add_sequential', 'drop_sequential')){
        lastkeep <- rep(TRUE, length(res))
        updatedModel <- MGmodel
        run_number <- 2L
        if(run_number > max_run)
            stop('max_run number must be greater than 1 for sequential searches')
        while(TRUE){
            statdiff <- do.call(c, lapply(res, function(x, stat){
                if(stat == 'p') return(x[2L, 'p'])
                return(x[1L, stat] - x[2L, stat])
                }, stat = seq_stat))
            if(seq_stat == 'p'){
                statdiff <- p.adjust(statdiff, p.adjust)
                keep <- statdiff < pval
            } else {
                keep <- statdiff < 0
            }
            if(all(keep == lastkeep)) break
            lastkeep <- keep
            if(verbose)
                cat(sprintf('\rChecking for DIF in %d more items',
                            ifelse(drop, sum(keep), sum(!keep))))
            if(ifelse(drop, sum(keep), sum(!keep)) == 0) break
            constrain <- updatedModel@constrain
            for(j in 1L:length(keep)){
                parnum <- list()
                for(i in 1L:length(which.par))
                    parnum[[i]] <- values$parnum[values$name == which.par[i] &
                                                     values$item == itemnames[j]]
                for(i in length(parnum):1L)
                    if(!length(parnum[[i]])) parnum[[i]] <- NULL
                if(!length(parnum)) break
                if(keep[j] && !drop){
                    for(i in 1L:length(parnum))
                        constrain[[length(constrain) + 1L]] <- parnum[[i]]
                } else if(!keep[j] && drop){
                    for(j in 1L:length(parnum)){
                        for(i in length(constrain):1L){
                            if(all(parnum[[j]] == sort(constrain[[i]])))
                                constrain[[i]] <- NULL
                        }
                    }
                }
            }
            updatedModel <- multipleGroup(MGmodel@Data$data, MGmodel@model[[1L]],
                                          group=MGmodel@Data$group,
                                          invariance = invariance, constrain=constrain,
                                          verbose = FALSE, ...)
            pick <- !keep
            if(drop) pick <- !pick
            tmp <- myLapply(X=items2test[pick], FUN=loop_test, model=updatedModel,
                            which.par=which.par, values=values, Wald=Wald, drop=drop,
                            itemnames=itemnames, invariance=invariance, return_models=FALSE, ...)
            names(tmp) <- itemnames[pick]
            for(i in names(tmp))
                res[[i]] <- tmp[[i]]
            if(run_number == max_run) break
            run_number <- run_number + 1L
        }
        if(verbose)
            cat('\nComputing final DIF estimates...')
        res <- myLapply(X=items2test[!keep], FUN=loop_test, model=updatedModel,
                        which.par=which.par, values=values, Wald=Wald, drop=FALSE,
                        itemnames=itemnames, invariance=invariance, return_models=return_models,
                        ...)
        names(res) <- itemnames[!keep]
    }

    for(i in 1L:length(res))
        attr(res[[i]], 'parnum') <- NULL
    if(return_models) return(res)
    if(p.adjust != 'none'){
        if(Wald){
            ps <- do.call(c, lapply(res, function(x) x$p))
        } else {
            ps <- do.call(c, lapply(res, function(x, stat){
                if(stat == 'p') return(x[2L, 'p'])
                return(x[1L, stat] - x[2L, stat])
            }, stat = 'p'))
        }
        ps <- p.adjust(ps, p.adjust)
        res$adj_pvals <- ps
    }
    if(plotdif && any(scheme %in% c('add', 'add_sequential'))){
        if(seq_stat != 'p'){
            statdiff <- do.call(c, lapply(res, function(x, stat){
                if(stat == 'p') return(x[2L, 'p'])
                return(x[1L, stat] - x[2L, stat])
            }, stat = seq_stat))
            keep <- statdiff < 0
        } else {
            statdiff <- res$adj_pvals
            if(is.null(statdiff)){
                if(Wald){
                    statdiff <- do.call(c, lapply(res, function(x) x$p))
                } else {
                    statdiff <- do.call(c, lapply(res, function(x, stat){
                        if(stat == 'p') return(x[2L, 'p'])
                        return(x[1L, stat] - x[2L, stat])
                    }, stat = 'p'))
                }
            }
            statdiff <- p.adjust(statdiff, p.adjust)
            keep <- statdiff > pval
        }
        which.item <- which(!keep)
        print(plot(MGmodel, type = type, which.items=which.item, facet_items=TRUE, ...))
    }
    return(res)
}
