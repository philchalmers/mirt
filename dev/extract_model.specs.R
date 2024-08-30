#' Extract model elements for fitting same model to different data
#'
#'
#' @param mod model fitted by the \code{\link{mirt}} function
#'
#' @return list of model information to fit the same single-group model
#'   to different data
#'
#' @examples
#' library(mirt)
#'
#' # fitted model
#'  mod <- mirt(expand.table(LSAT7),
#'              model="F = 1-5
#'              CONSTRAIN = (3, 5, g)",
#'              itemtype = c('2PL', '2PL', '3PL', '2PL', '3PL'))
#'
#' coef(mod, simplify=TRUE)
#' model_specs <- extract_model.specs(mod)
#'
#' # check
#' same.mod <- mirt(expand.table(LSAT7), model_specs=model_specs)
#' coef(same.mod, simplify=TRUE)
#' anova(mod, same.mod)
#'
#' # use with new data
#' new.mod <- mirt(expand.table(LSAT6), model_specs=model_specs)
#' coef(new.mod, simplify=TRUE)
#'
extract_model.specs <- function(mod){
    itemtype <- extract.mirt(mod, 'itemtype')
    model <- extract.mirt(mod, 'model')
    constrain <- extract.mirt(mod, 'constrain')
    monopoly.k <- extract.mirt(mod, 'monopoly.k')
    gpcm_mats <- extract.mirt(mod, 'gpcm_mats')
    grsm.block <- extract.mirt(mod, "grsm.block")
    rsm.block <- extract.mirt(mod, "rsm.block")
    key <- extract.mirt(mod, "key")
    dentype <- extract.mirt(mod, 'dentype')

    # maybe formula?
    ret <- list(itemtype=itemtype, model=model, constrain=constrain,
                monopoly.k=monopoly.k, gpcm_mats=gpcm_mats, dentype=dentype,
                grsm.block=grsm.block, rsm.block=rsm.block, key=key)
    ret
}

replace_model.specs <- function(lst, env){
    if(!is.null(lst)){
        nms <- names(lst)
        for(i in 1L:length(lst))
            assign(nms[i], lst[[i]], envir = env)
    }
    invisible(NULL)
}

# in mirt() use
# dots <- list(...)
# replace_model.specs(dots$model_specs, env=environment())