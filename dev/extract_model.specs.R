#' Extract model elements
#'
#'
#' @param mod model fitted by mirt package
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
#' model_specs <- extract_model.specs(prev.mirt.model)
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
    browser()

    # maybe formula?
    ret <- list(itemtype=itemtype, model=model, constrain=constrain,
                monopoly.k=monopoly.k, gpcm_mats=gpcm_mats,
                grsm.block=grsm.block, rsm.block=rsm.block, key=key)
    ret
}
