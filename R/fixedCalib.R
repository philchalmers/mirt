#' Fixed-item calibration method
#'
#' Implements the set of fixed-item calibration methods described by Kim (2006). The initial
#' calibrated model must be fitted via \code{\link{mirt}}, is currently limited to
#' unidimensional models only, and should only be utilized when the new set of responses
#' are obtained from a population with similar distributional characteristics in the latent traits.
#' For more flexible calibration of items, including a fixed-item calibration variant involving
#' anchor items for equating, see \code{\link{multipleGroup}}.
#'
#' @param old_mod a model of class SingleGroupClass fitted using \code{\link{mirt}}
#'
#' @param data new data to be used for calibration. Note that to be consistent
#'   with the \code{mod} object, observed responses/NA placeholders must be included
#'   to link the item names used in the original \code{mod} definition
#'   (i.e., \code{extract.mirt(mod, what = 'itemnames')})
#'
#' @param PAU prior ability update (PAU) approach. Supports none (\code{"NWU"}),
#'   one (\code{"OWU"}), and many (\code{"MWU"})
#'
#' @param NEM number of EM cycles (NEMC) to use for the to-be-estimated parameters.
#'   Supports one (\code{"OEM"}) and many (\code{"MEM"})
#'
#' @param technical list of technical estimation arguments
#'   (see \code{\link{mirt}} for details)
#'
#' @param ... additional arguments to pass to \code{\link{mirt}}
#'
#' @seealso \code{\link{mirt}}, \code{\link{multipleGroup}}
#'
#' @export
#'
#' @references
#'
#' Kim, S. (2006). A comparative study of IRT fixed parameter calibration methods.
#' \emph{Journal of Educational Measurement, 4}(43), 355–381.
#'
#' @examples
#' \dontrun{
#'
#' # single factor
#' set.seed(12345)
#' J <- 50
#' a <- matrix(abs(rnorm(J,1,.3)), ncol=1)
#' d <- matrix(rnorm(J,0,.7),ncol=1)
#' itemtype <- rep('2PL', nrow(a))
#'
#' # calibration data theta ~ N(0,1)
#' N <- 3000
#' dataset1 <- simdata(a, d, N = N, itemtype=itemtype)
#'
#' # new data (again, theta ~ N(0,1))
#' dataset2 <- simdata(a, d, N = 1000, itemtype=itemtype)
#'
#' # last 40% of experimental items not given to calibration group
#' #     (unobserved; hence removed)
#' dataset1 <- dataset1[,-c(J:(J*.6))]
#' head(dataset1)
#'
#' # assume first 60% of items not given to new group
#' dataset2[,colnames(dataset1)] <- NA
#' head(dataset2)
#'
#' #--------------------------------------
#'
#' # calibrated model from dataset1 only
#' mod <- mirt(dataset1, model = 1)
#' coef(mod, simplify=TRUE)
#'
#' # No Prior Weights Updating and One EM Cycle (NWU-OEM)
#' NWU_OEM <- fixedCalib(dataset2, model = 1, old_mod = mod, PAU = 'NWU', NEMC = 'OEM')
#' coef(NWU_OEM, simplify=TRUE)
#' data.frame(coef(NWU_OEM, simplify=TRUE)$items[,c('a1','d')], pop_a1=a, pop_d=d)
#' plot(NWU_OEM, type = 'empiricalhist')
#'
#' # No Prior Weights Updating and Multiple EM Cycles (NWU-MEM)
#' NWU_MEM <- fixedCalib(dataset2, model = 1, old_mod = mod, PAU = 'NWU')
#' coef(NWU_MEM, simplify=TRUE)
#' data.frame(coef(NWU_MEM, simplify=TRUE)$items[,c('a1','d')], pop_a1=a, pop_d=d)
#' plot(NWU_MEM, type = 'empiricalhist')
#'
#' # One Prior Weights Updating and One EM Cycle (OWU-OEM)
#' OWU_OEM <- fixedCalib(dataset2, model = 1, old_mod = mod, PAU = 'OWU', NEMC = "OEM")
#' coef(OWU_OEM, simplify=TRUE)
#' data.frame(coef(OWU_OEM, simplify=TRUE)$items[,c('a1','d')], pop_a1=a, pop_d=d)
#' plot(OWU_OEM, type = 'empiricalhist')
#'
#' # One Prior Weights Updating and Multiple EM Cycles (OWU-MEM)
#' OWU_MEM <- fixedCalib(dataset2, model = 1, old_mod = mod, PAU = 'OWU')
#' coef(OWU_MEM, simplify=TRUE)
#' data.frame(coef(OWU_MEM, simplify=TRUE)$items[,c('a1','d')], pop_a1=a, pop_d=d)
#' plot(OWU_MEM, type = 'empiricalhist')
#'
#' # Multiple Prior Weights Updating and Multiple EM Cycles (MWU-MEM)
#' MWU_MEM <- fixedCalib(dataset2, model = 1, old_mod = mod)
#' coef(MWU_MEM, simplify=TRUE)
#' data.frame(coef(MWU_MEM, simplify=TRUE)$items[,c('a1','d')], pop_a1=a, pop_d=d)
#' plot(MWU_MEM, type = 'empiricalhist')
#'
#' }
fixedCalib <- function(data, model = 1, old_mod,  PAU = 'MWU', NEMC = "MEM",
                       technical = list(), ...){

    stopifnot(PAU %in% c('NWU', 'OWU', 'MWU'))
    stopifnot(NEMC %in% c("OEM", "MEM"))

    global_prior <- NULL
    custom_den_const <- function(obj, Theta) global_prior
    custom_den_EH <- function(obj, Theta) {
        ret <- if(length(obj@rr) == 0L) # very first EM iteration use global info
            global_prior
        else obj@rr / sum(obj@rr) # prior from normalized E-table of counts
        ret
    }

    nfact <- extract.mirt(old_mod, 'nfact')
    stopifnot(nfact == 1L)
    old_itemnames <- extract.mirt(old_mod, 'itemnames')
    stopifnot(all(old_itemnames %in% colnames(data)))

    olddata <- extract.mirt(old_mod, 'data')
    tmp <- t(data.frame(rep(NA, ncol(data) - ncol(olddata))))
    rownames(tmp) <- NULL
    olddata <- data.frame(olddata, tmp)
    colnames(olddata) <- c(old_itemnames, colnames(data)[!(colnames(data) %in% old_itemnames)])
    fulldata <- rbind(olddata, data)

    sv <- mod2values(old_mod)
    sv$est <- FALSE
    sv2 <- mirt(fulldata, model, pars='values')
    for(item in c(old_itemnames, 'GROUP')){
        pick1 <- sv$item == item
        pick2 <- sv2$item == item
        sv2$value[pick2] <- sv$value[pick1]
        sv2$est[pick2] <- FALSE
    }

    if(is.null(technical$NCYCLES))
        if(NEMC == 'OEM') technical$NCYCLES <- 1L
    technical$message <- FALSE

    if(PAU %in% c('OWU', 'MWU')){
        FC_mod_den <- mirt(data=extract.mirt(old_mod, 'data'),
                           model=extract.mirt(old_mod, 'model'),
                           pars=sv, dentype = 'EH', verbose=FALSE,
                           technical = list(NCYCLES = 1L, message = FALSE))
        global_prior <- FC_mod_den@Internals$Prior[[1L]]
        den <- if(PAU == 'OWU') custom_den_const else custom_den_EH
        par <- c(MEAN_1 = 0, COV_11 = 1) # names preserved to avoid new starting value issues
        est <- c(FALSE, FALSE)
        grp <- createGroup(par, est, den, nfact = 1L)
    }

    mod <- if(PAU == "NWU"){
        mirt(fulldata, model, pars=sv2, technical=technical, ...)
    } else if(PAU %in% c("OWU", "MWU")){
        mirt(fulldata, model, pars=sv2, customGroup=grp, dentype='EH',
             technical=technical, ...)
    }
    sv_final <- mod2values(mod)
    FC_mod_den <- mirt(fulldata, model, pars=sv_final, dentype='EH', TOL=NaN, ...)
    FC_mod_den
}