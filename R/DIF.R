#' Differential item functioning tests
#'
#' This function runs the Wald and likelihood-ratio approaches for testing differential
#' item functioning (DIF). This is primarily a convenience wrapper to the \code{\link{multipleGroup}}
#' function for performing standard DIF procedures. Models can be estimated in
#' parallel automatically by defining a parallel object with \code{\link{mirtCluster}}, which will help
#' to decrease the runtime.
#'
#' @aliases DIF
#' @param MGmodel an object returned from \code{\link{multipleGroup}} to be used as the reference model
#' @param which.par a character vector containing the parameter names which will be inspected for DIF
#' @param Wald logical; perform Wald tests for DIF instead of likelihood ratio test?
#' @param items2test a numeric vector, or character vector containing the item names, indicating which items
#'   will be tested for DIF. In models where anchor items are known, omit them from this vector. For example,
#'   if items 1 and 2 are anchors in a 10 item test, then \code{items2test = 3:10} would work for testing the
#'   remaining items (important to remember when using sequential schemes)
#' @param scheme type of DIF analysis to perform, either by adding or dropping constraints across groups.
#'   These can be:
#' \describe{
#'   \item{'add'}{parameters in \code{which.par} will be constrained each item one at a time for items that are
#'     specified in \code{items2test}. This is beneficial when examining DIF from a model with parameters
#'     freely estimated across groups, and when inspecting differences via the Wald test}
#'   \item{'drop'}{parameters in \code{which.par} will be freely estimated for items that are
#'     specified in \code{items2test}. This is useful when supplying an overly restrictive model and attempting to
#'     detect DIF with a slightly less restrictive model}
#'   \item{'add_sequential'}{sequentially loop over the items being tested, and at the end of the loop treat
#'     DIF tests that satisfy the \code{seq_stat} criteria as invariant. The loop is then re-run on the remaining
#'     invariant items to determine if they are now displaying DIF in the less constrained model,
#'     and when no new invariant item is found the algorithm stops and returns the items that displayed DIF}
#'   \item{'drop_sequential'}{sequentially loop over the items being tested, and at the end of the loop treat
#'     items that violate the \code{seq_stat} criteria as demonstrating DIF. The loop is then re-run, leaving the items
#'     that previously demonstrated DIF as variable across groups, and the remaining test items that previously showed
#'     invariance are re-tested. The algorithm stops when no more items showing DIF are found and returns the items that
#'     displayed DIF}
#' }
#' @param seq_stat select a statistic to test for in the sequential schemes. Potential values are
#'   (in descending order of power) \code{'AIC'}, \code{'AICc'}, \code{'SABIC'}, and \code{'BIC'}.
#'   If a numeric value is input that ranges between 0 and 1, the 'p' value will be tested
#'   (e.g., \code{seq_stat = .05} will test for the difference of p < .05 in the add scheme,
#'   or p > .05 in the drop scheme), along with the specified \code{p.adjust} input
#' @param p.adjust string to be passed to the \code{\link{p.adjust}} function to adjust p-values.
#'   Adjustments are located in the \code{adj_pvals} element in the returned list
#' @param verbose logical print extra information to the console?
#' @param ... additional arguments to be passed to \code{\link{multipleGroup}}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords DIF
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
#' #  Information matrix with S-EM, therefore drop TOL for better accuracy
#' model <- multipleGroup(dat, 1, group, SE = TRUE, technical = list(TOL = 1e-6))
#'
#' #test whether freeing slopes and intercepts results in DIF
#' resulta1d <- DIF(model, c('a1', 'd'))
#' resulta1d
#'
#' #same as above, but using Wald tests with Benjamini & Hochberg adjustment
#' resulta1dWald <- DIF(model, c('a1', 'd'), Wald = TRUE, p.adjust = 'fdr')
#' resulta1dWald
#'
#' #test whether freeing only slopes results in DIF for all items
#' resulta1 <- DIF(model, 'a1')
#' resulta1
#'
#' #following up on resulta1d, to determine whether it's a1 or d parameter causing DIF
#' itemplot(model, 1)
#' itemplot(model, 2)
#' itemplot(model, 3)
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
#' ### drop down approach when specifying a highly constrained model
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
#' #step down procedure
#' model <- multipleGroup(dat, 1, group, invariance = itemnames)
#' stepdown <- DIF(model, c('a1', 'd'), scheme = 'drop_sequential')
#' stepdown
#'
#' }
DIF <- function(MGmodel, which.par, scheme = 'add', items2test = 1:ncol(MGmodel@data),
                seq_stat = 'SABIC', Wald = FALSE, p.adjust = 'none', verbose = TRUE, ...){

    loop_test <- function(item, model, which.par, values, Wald, itemnames, invariance, drop, ...)
    {
        constrain <- model@constrain
        parnum <- list()
        for(i in 1L:length(which.par))
            parnum[[i]] <- values$parnum[values$name == which.par[i] &
                                             values$item == itemnames[item]]
        if(!length(parnum))
            stop('Item ', item, ' does not contain any of the parameters defined in which.par.
                 Consider removing it from the item2test input or adding relevant parameters to which.par')
        if(Wald){
            wv <- wald(model)
            L <- matrix(0, length(parnum), length(wv))
            for(i in 1L:length(parnum)){
                L[i, paste0(which.par[i], '.', parnum[[i]][1L]) == wv] <- 1
                L[i, paste0(which.par[i], '.', parnum[[i]][2L]) == wv] <- -1
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
        newmodel <- multipleGroup(model@data, model@model[[1L]], group=model@group,
                                  invariance = invariance, constrain=constrain,
                                  verbose = FALSE, ...)
        aov <- anova(newmodel, model, verbose = FALSE)
        attr(aov, 'parnum') <- parnum
        return(aov)
    }

    itemnames <- colnames(MGmodel@data)
    if(!any(scheme %in% c('add', 'drop', 'add_sequential', 'drop_sequential')))
        stop('scheme input is not valid')
    if(Wald){
        if(scheme != 'add')
            stop('Wald test are only appropriate when add scheme is used')
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
    data <- MGmodel@data
    invariance <- MGmodel@invariance
    values <- mod2values(MGmodel)
    drop <- scheme == 'drop' || scheme == 'drop_sequential'
    invariance <- MGmodel@invariance[MGmodel@invariance %in% c('free_means', 'free_var', 'free_varcov', 'free_cov')]
    if(!length(invariance)) invariance <- ''

    if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
        res <- parLapply(cl=mirtClusterEnv$MIRTCLUSTER, X=items2test, fun=loop_test,
                         model=MGmodel, drop=drop, which.par=which.par, values=values,
                         Wald=Wald, itemnames=itemnames, invariance=invariance, ...)
    } else {
        res <- lapply(X=items2test, FUN=loop_test, model=MGmodel,
                      which.par=which.par, values=values, Wald=Wald, drop=drop,
                      itemnames=itemnames, invariance=invariance, ...)
    }
    names(res) <- itemnames[items2test]
    if(scheme %in% c('add_sequential', 'drop_sequential')){
        lastkeep <- rep(TRUE, length(res))
        updatedModel <- MGmodel
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
                cat(sprintf('\rChecking for DIF in %d more items', ifelse(drop, sum(keep), sum(!keep))))
            constrain <- updatedModel@constrain
            for(j in 1L:length(keep)){
                parnum <- list()
                for(i in 1L:length(which.par))
                    parnum[[i]] <- values$parnum[values$name == which.par[i] &
                                                     values$item == itemnames[j]]
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
            updatedModel <- multipleGroup(MGmodel@data, MGmodel@model[[1L]], group=MGmodel@group,
                                          invariance = invariance, constrain=constrain,
                                          verbose = FALSE, ...)
            pick <- !keep
            if(drop) pick <- !pick
            if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
                tmp <- parLapply(cl=mirtClusterEnv$MIRTCLUSTER, X=items2test[pick], fun=loop_test,
                                 model=updatedModel, drop=drop, which.par=which.par, values=values,
                                 Wald=Wald, itemnames=itemnames, invariance=invariance, ...)
            } else {
                tmp <- lapply(X=items2test[pick], FUN=loop_test, model=updatedModel,
                              which.par=which.par, values=values, Wald=Wald, drop=drop,
                              itemnames=itemnames, invariance=invariance, ...)
            }
            names(tmp) <- itemnames[pick]
            for(i in names(tmp))
                res[[i]] <- tmp[[i]]
        }
        if(verbose)
            cat('\nComputing final DIF estimates...')
        if(!is.null(mirtClusterEnv$MIRTCLUSTER)){
            res <- parLapply(cl=mirtClusterEnv$MIRTCLUSTER, X=items2test[!keep], fun=loop_test,
                             model=updatedModel, drop=FALSE, which.par=which.par, values=values,
                             Wald=Wald, itemnames=itemnames, invariance=invariance, ...)
        } else {
            res <- lapply(X=items2test[!keep], FUN=loop_test, model=updatedModel,
                          which.par=which.par, values=values, Wald=Wald, drop=FALSE,
                          itemnames=itemnames, invariance=invariance, ...)
        }
        names(res) <- itemnames[!keep]
    }

    for(i in 1L:length(res))
        attr(res[[i]], 'parnum') <- NULL
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
    return(res)
}