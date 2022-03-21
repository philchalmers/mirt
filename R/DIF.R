#' Differential item functioning statistics
#'
#' This function runs the Wald and likelihood-ratio approaches for testing differential
#' item functioning (DIF). This is primarily a convenience wrapper to the
#' \code{\link{multipleGroup}} function for performing standard DIF procedures. Independent
#' models can be estimated in parallel by defining a parallel object with \code{\link{mirtCluster}},
#' which will help to decrease the runtime. For best results, the baseline model should contain
#' a set of 'anchor' items and have freely estimated hyper-parameters in the focal groups.
#'
#' Generally, the precomputed baseline model should have been
#' configured with two estimation properties: 1) a set of 'anchor' items,
#' where the anchor items have various parameters that have been constrained to be equal
#' across the groups, and 2) contain freely estimated latent mean and variance terms in
#' all but one group (the so-called 'reference' group).
#' These two properties help to fix the metric of the groups so that
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
#' @param return_seq_model logical; on the last iteration of the sequential schemes, return
#'   the fitted multiple-group model containing the freely estimated parameters indicative of
#'   DIF? This is generally only useful when \code{scheme = 'add_sequential'}. Default is FALSE
#' @param simplify logical; simplify the output by returning a data.frame object with
#'   the differences between AIC, BIC, etc, as well as the chi-squared test (X2) and associated
#'   df and p-values
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
#'     returns the items that displayed DIF. Note that the DIF statistics are relative to this final,
#'     less constrained model which includes the DIF effects}
#'   \item{'drop_sequential'}{sequentially loop over the items being tested, and at the end of the
#'     loop treat items that violate the \code{seq_stat} criteria as demonstrating DIF. The loop is
#'     then re-run, leaving the items that previously demonstrated DIF as variable across groups,
#'     and the remaining test items that previously showed invariance are re-tested. The algorithm
#'     stops when no more items showing DIF are found and returns the items that displayed DIF.
#'     Note that the DIF statistics are relative to this final,
#'     less constrained model which includes the DIF effects}
#' }
#' @param seq_stat select a statistic to test for in the sequential schemes. Potential values are
#'   (in descending order of power) \code{'AIC'}, \code{'SABIC'}, \code{'HQ'}, and \code{'BIC'}.
#'   If a numeric value is input that ranges between 0 and 1, the 'p' value will be tested
#'   (e.g., \code{seq_stat = .05} will test for the difference of p < .05 in the add scheme,
#'   or p > .05 in the drop scheme), along with the specified \code{p.adjust} input
#' @param max_run a number indicating the maximum number of cycles to perform in sequential
#'   searches. The default is to perform search until no further DIF is found
#' @param plotdif logical; create item plots for items that are displaying DIF according to the
#'   \code{seq_stat} criteria? Only available for 'add' type schemes
#' @param type the \code{type} of plot argument passed to \code{plot()}. Default is 'trace', though
#'   another good option is 'infotrace'. For ease of viewing, the \code{facet_item} argument to
#'   mirt's \code{plot()} function is set to \code{TRUE}
#' @param p.adjust string to be passed to the \code{\link{p.adjust}} function to adjust p-values.
#'   Adjustments are located in the \code{adj_pvals} element in the returned list
#' @param verbose logical print extra information to the console?
#' @param ... additional arguments to be passed to \code{\link{multipleGroup}} and \code{plot}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Chalmers, R. P., Counsell, A., and Flora, D. B. (2016). It might not
#'   make a big DIF: Improved Differential Test Functioning statistics that account for
#'   sampling variability. \emph{Educational and Psychological Measurement, 76}, 114-140.
#'   \doi{10.1177/0013164415584576}
#' @keywords differential item functioning
#' @seealso \code{\link{multipleGroup}}, \code{\link{DRF}}
#, \code{\link{DTF}}
#' @export DIF
#' @examples
#' \dontrun{
#'
#' # simulate data where group 2 has a smaller slopes and more extreme intercepts
#' set.seed(12345)
#' a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
#' a2[1:2, ] <- a1[1:2, ]/3
#' d1[c(1,3), ] <- d2[c(1,3), ]/4
#' head(data.frame(a.group1 = a1, a.group2 = a2, d.group1 = d1, d.group2 = d2))
#' itemtype <- rep('2PL', nrow(a1))
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
#' # Information matrix with Oakes' identity (not controlling for latent group differences)
#' # NOTE: Without properly equating the groups the following example code is not testing for DIF,
#'      # but instead reflects a combination of DIF + latent-trait distribution effects
#' model <- multipleGroup(dat, 1, group, SE = TRUE)
#'
#' # Likelihood-ratio test for DIF (as well as model information)
#' DIF(model, c('a1', 'd'))
#' DIF(model, c('a1', 'd'), simplify=FALSE) # return list output
#'
#' # same as above, but using Wald tests with Benjamini & Hochberg adjustment
#' DIF(model, c('a1', 'd'), Wald = TRUE, p.adjust = 'fdr')
#'
#' # equate the groups by assuming the last 5 items have no DIF
#' itemnames <- colnames(dat)
#' model <- multipleGroup(dat, 1, group, SE = TRUE,
#'    invariance = c(itemnames[11:ncol(dat)], 'free_means', 'free_var'))
#'
#' # test whether adding slopes and intercepts constraints results in DIF. Plot items showing DIF
#' resulta1d <- DIF(model, c('a1', 'd'), plotdif = TRUE, items2test=1:10)
#' resulta1d
#'
#' # test whether adding only slope constraints results in DIF for all items
#' DIF(model, 'a1', items2test=1:10)
#'
#' # Determine whether it's a1 or d parameter causing DIF (could be joint, however)
#' (a1s <- DIF(model, 'a1', items2test = 1:3))
#' (ds <- DIF(model, 'd', items2test = 1:3))
#'
#' ### drop down approach (freely estimating parameters across groups) when
#' ### specifying a highly constrained model with estimated latent parameters
#' model_constrained <- multipleGroup(dat, 1, group,
#'   invariance = c(colnames(dat), 'free_means', 'free_var'))
#' dropdown <- DIF(model_constrained, c('a1', 'd'), scheme = 'drop')
#' dropdown
#'
#' ### sequential schemes (add constraints)
#'
#' ### sequential searches using SABIC as the selection criteria
#' # starting from completely different models
#' stepup <- DIF(model, c('a1', 'd'), scheme = 'add_sequential',
#'               items2test=1:10)
#' stepup
#'
#' # step down procedure for highly constrained model
#' stepdown <- DIF(model_constrained, c('a1', 'd'), scheme = 'drop_sequential')
#' stepdown
#'
#' # view final MG model (only useful when scheme is 'add_sequential')
#' updated_mod <- DIF(model, c('a1', 'd'), scheme = 'add_sequential',
#'                return_seq_model=TRUE)
#' plot(updated_mod, type='trace')
#'
#' }
DIF <- function(MGmodel, which.par, scheme = 'add', items2test = 1:extract.mirt(MGmodel, 'nitems'),
                seq_stat = 'SABIC', Wald = FALSE, p.adjust = 'none', return_models = FALSE,
                return_seq_model = FALSE, max_run = Inf, plotdif = FALSE, type = 'trace',
                simplify = TRUE, verbose = TRUE, ...){

    loop_test <- function(item, model, which.par, values, Wald, itemnames, invariance, drop,
                          return_models, technical = list(), ...)
    {
        constrain <- model@Model$constrain
        mirt_model <- model@Model$model
        if(is(mirt_model, 'mirt.model'))
            mirt_model$x <- mirt_model$x[mirt_model$x[,"Type"] != 'CONSTRAINB',
                                         , drop=FALSE]
        technical$omp <- FALSE
        parnum <- list()
        for(i in seq_len(length(which.par)))
            parnum[[i]] <- values$parnum[values$name == which.par[i] &
                                             values$item == itemnames[item]]
        for(i in length(parnum):1L)
            if(!length(parnum[[i]])) parnum[[i]] <- NULL
        if(!length(parnum))
            stop('Item ', item, ' does not contain any of the parameters defined in which.par.
                 Consider removing it from the item2test input or adding relevant parameters
                 to which.par', call.=FALSE)
        if(Wald){
            wv <- wald(model)
            infoname <- names(wv)
            L <- matrix(0, length(parnum), length(infoname))
            for(i in seq_len(length(parnum))){
                L[i, paste0(which.par[i], '.', parnum[[i]][1L]) == infoname] <- 1
                L[i, paste0(which.par[i], '.', parnum[[i]][2L]) == infoname] <- -1
            }
            res <- wald(model, L)
            return(res)
        }
        if(drop){
            sv <- values
            for(j in seq_len(length(parnum))){
                for(i in length(constrain):1L){
                    if(all(parnum[[j]] %in% sort(constrain[[i]])))
                        constrain[[i]] <- NULL
                }
            }
        } else {
            sv <- NULL
            for(i in seq_len(length(parnum)))
                constrain[[length(constrain) + 1L]] <- parnum[[i]]
        }
        newmodel <- multipleGroup(model@Data$data, mirt_model, group=model@Data$group,
                                  invariance = invariance, constrain=constrain, pars=sv,
                                  itemtype = model@Model$itemtype, verbose=FALSE, technical=technical,
                                  ...)
        aov <- anova(newmodel, model, verbose = FALSE)
        attr(aov, 'parnum') <- parnum
        attr(aov, 'converged') <- extract.mirt(newmodel, 'converged')
        if(return_models) aov <- newmodel
        return(aov)
    }

    if(missing(MGmodel)) missingMsg('MGmodel')
    if(missing(which.par)) missingMsg('which.par')
    if(!is(MGmodel, 'MultipleGroupClass'))
        stop('Input model must be fitted by multipleGroup()', call.=FALSE)
    aov <- anova(MGmodel)
    has_priors <- !is.null(aov$logPost)
    if(has_priors && is.numeric(seq_stat))
        stop('p-value seq_stat for models fitted with Bayesian priors in not meaningful. Please select alternative',
             call.=FALSE)

    if(!any(sapply(MGmodel@ParObjects$pars, function(x, pick) x@ParObjects$pars[[pick]]@est,
                   pick = MGmodel@Data$nitems + 1L)))
        message(paste('No hyper-parameters were estimated in the DIF model. For effective',
                'DIF testing, freeing the focal group hyper-parameters is recommended.'))
    bfactorlist <- MGmodel@Internals$bfactor
    if(!is.null(bfactorlist$Priorbetween[[1L]]))
        stop('bifactor models are currently not supported in this function', call.=FALSE)
    itemnames <- colnames(MGmodel@Data$data)
    if(!any(scheme %in% c('add', 'drop', 'add_sequential', 'drop_sequential')))
        stop('scheme input is not valid', call.=FALSE)
    if(return_models){
        if(Wald)
            stop('return_models argument only valid for likelihood ratio tests', call.=FALSE)
    }
    if(Wald){
        if(scheme != 'add')
            stop('Wald tests are only appropriate when add scheme is used', call.=FALSE)
        if(!MGmodel@Options$SE)
            stop('Information matrix was not calculated', call.=FALSE)
    }
    if(plotdif && any(scheme %in% c('drop', 'drop_sequential')))
        stop('plotdif not supported for dropping schemes', call.=FALSE)
    pval <- 0
    if(is.numeric(seq_stat)){
        pval <- seq_stat
        seq_stat <- 'p'
    } else if(!any(seq_stat %in% c('p', 'AIC', 'SABIC', 'BIC', 'DIC', 'HQ'))){
        stop('Invalid seq_stat input', call.=FALSE)
    }
    if(is.character(items2test)) items2test <- which(items2test %in% itemnames)
    invariance <- MGmodel@Model$invariance
    values <- mod2values(MGmodel)
    drop <- scheme == 'drop' || scheme == 'drop_sequential'
    invariance <- MGmodel@Model$invariance[MGmodel@Model$invariance %in%
                                         c('free_means', 'free_var')]
    if(!length(invariance)) invariance <- ''
    res <- myLapply(X=items2test, FUN=loop_test, progress=verbose,
                    model=MGmodel, which.par=which.par, values=values,
                    Wald=Wald, drop=drop, itemnames=itemnames, invariance=invariance,
                    return_models=return_models, ...)
    names(res) <- itemnames[items2test]
    if(scheme %in% c('add_sequential', 'drop_sequential')){
        lastkeep <- rep(TRUE, length(res))
        updatedModel <- MGmodel
        run_number <- 2L
        if(run_number > max_run)
            stop('max_run number must be greater than 1 for sequential searches', call.=FALSE)
        while(TRUE){
            statdiff <- do.call(c, lapply(res, function(x, stat){
                if(stat == 'p') return(x[2L, 'p'])
                return(x[1L, stat] - x[2L, stat])
                }, stat = seq_stat))
            if(seq_stat == 'p'){
                statdiff <- p.adjust(statdiff, p.adjust)
                keep <- statdiff >= pval
            } else {
                keep <- statdiff <= 0
            }
            if(run_number == 2L && (all(!keep) || all(keep))){
                if(verbose)
                    message('sequential scheme not required; all/no items contain DIF on first iteration')
                if(return_seq_model) return(MGmodel)
                if(scheme == 'add_sequential'){
                    scheme <- 'add'
                    if(all(!keep)) break
                }
                if(scheme == 'drop_sequential'){
                    scheme <- 'drop'
                    if(all(!keep)) break
                }
                ret <- data.frame()
                ret <- as.mirt_df(ret)
                return(ret)
            }
            if(all(keep == lastkeep)) break
            if(drop && run_number > 2L){
                lastkeep <- keep | lastkeep
            } else lastkeep <- keep
            if(verbose)
                cat(sprintf('\rChecking for DIF in %d more items', if(drop) sum(keep) else sum(!keep)))
            if(ifelse(drop, sum(keep), sum(!keep)) == 0) break
            constrain <- updatedModel@Model$constrain
            for(j in seq_len(length(keep))){
                parnum <- list()
                for(i in seq_len(length(which.par)))
                    parnum[[i]] <- values$parnum[values$name == which.par[i] &
                                                     values$item == itemnames[j]]
                for(i in length(parnum):1L)
                    if(!length(parnum[[i]])) parnum[[i]] <- NULL
                if(!length(parnum)) break
                if(keep[j] && !drop){
                    for(i in 1L:length(parnum))
                        constrain[[length(constrain) + 1L]] <- parnum[[i]]
                } else if(!keep[j] && drop){
                    for(j in seq_len(length(parnum))){
                        for(i in length(constrain):1L){
                            if(all(parnum[[j]] == sort(constrain[[i]])))
                                constrain[[i]] <- NULL
                        }
                    }
                }
            }
            updatedModel <- multipleGroup(MGmodel@Data$data, MGmodel@Model$model,
                                          group=MGmodel@Data$group, itemtype=MGmodel@Model$itemtype,
                                          invariance = invariance, constrain=constrain,
                                          verbose = FALSE, ...)
            pick <- !keep
            if(drop) pick <- !pick
            tmp <- myLapply(X=items2test[pick], FUN=loop_test, progress=verbose, model=updatedModel,
                            which.par=which.par, values=values, Wald=Wald, drop=drop,
                            itemnames=itemnames, invariance=invariance, return_models=FALSE, ...)
            names(tmp) <- itemnames[items2test][pick]
            for(i in names(tmp))
                res[[i]] <- tmp[[i]]
            if(run_number == max_run) break
            run_number <- run_number + 1L
        }

        if(verbose && !(scheme %in% c('add', 'drop')))
            cat('\nComputing final DIF estimates...\n')
        pick <- !lastkeep
        if(return_seq_model) return(updatedModel)
        if(!(scheme %in% c('add', 'drop'))){ # will equal 'add/drop' if all items on first loop have DIF
            res <- myLapply(X=items2test[pick], FUN=loop_test, progress=verbose, model=updatedModel,
                            which.par=which.par, values=values, Wald=Wald, drop=FALSE,
                            itemnames=itemnames, invariance=invariance, return_models=return_models,
                            ...)
            names(res) <- itemnames[items2test][pick]
        }
    }
    converged <- logical(length(res))
    for(i in seq_len(length(res))){
        if(!Wald) converged[i] <- attr(res[[i]], 'converged')
        attr(res[[i]], 'converged') <- attr(res[[i]], 'parnum') <- NULL
    }
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
            keep <- statdiff <= 0
        } else {
            statdiff <- res$adj_pvals
            if(is.null(statdiff)){
                if(Wald){
                    statdiff <- do.call(c, lapply(res, function(x) x$p))
                } else {
                    statdiff <- do.call(c, lapply(res, function(x, stat){
                        if(stat == 'p') return(x[2L, 'p'])
                        return(x[1L, stat] - x[2L, stat])
                    }, stat = seq_stat))
                }
            }
            if(seq_stat == 'p' || Wald){
                statdiff <- p.adjust(statdiff, p.adjust)
                keep <- statdiff >= pval
            } else {
                keep <- !(statdiff <= 0 & !sapply(statdiff, closeEnough, low=-1e-4, up=1e-4))
            }
        }
        which.item <- which(!keep)
        if(length(which.item)){
            print(plot(MGmodel, type = type, which.items=which.item, facet_items=TRUE, ...))
        } else {
            message('No DIF items were detected for plotting.')
        }
    }
    pick <- names(res)
    pick <- pick[pick != 'adj_pvals']
    if(Wald && !return_models){
        adj_pvals <- res$adj_pvals
        res <- do.call(rbind, res[pick])
        res$adj_pvals <- adj_pvals
        res <- as.mirt_df(res)
        return(res)
    }
    if(simplify && !return_models){
        adj_pvals <- res$adj_pvals
        out <- lapply(res[pick], function(x){
             r <- x[2L, ] - x[1L, ]
             if(!has_priors)
                r[,c("X2", 'df', 'p')] <- x[2L, c("X2", 'df', 'p')]
             else r[,c('df')] <- x[2L, c('df')]
             r$logLik <- NULL
             r
         })
         res <- cbind(converged, do.call(rbind, out))
         if(!has_priors) res$adj_pvals <- adj_pvals
         res <- as.mirt_df(res)
    }
    return(res)
}
