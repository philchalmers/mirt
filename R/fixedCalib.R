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
#' @param model type of model to fit for the complete dataset (not that for the fixed items
#'   in \code{old_mod} the factor loadings/constraints specified by the potential \code{\link{mirt.model}}
#'   specification is not relevant)
#'
#' @param data new data to be used for calibration. Note that to be consistent
#'   with the \code{mod} object, observed responses/NA placeholders must be included
#'   to link the item names used in the original \code{mod} definition
#'   (i.e., \code{extract.mirt(mod, what = 'itemnames')})
#'
#' @param PAU prior ability update (PAU) approach. Supports none (\code{"NWU"}),
#'   one (\code{"OWU"}), and many (\code{"MWU"})
#'
#' @param NEMC number of EM cycles (NEMC) to use for the to-be-estimated parameters.
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
#' \emph{Journal of Educational Measurement, 4}(43), 355-381.
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
#' #--------------------------------------
#'
#' # calibrated model from dataset1 only
#' mod <- mirt(dataset1, model = 1)
#' coef(mod, simplify=TRUE)
#'
#' # No Prior Weights Updating and One EM Cycle (NWU-OEM)
#' NWU_OEM <- fixedCalib(dataset2, model=1, old_mod=mod, PAU='NWU', NEMC='OEM')
#' coef(NWU_OEM, simplify=TRUE)
#' data.frame(coef(NWU_OEM, simplify=TRUE)$items[,c('a1','d')],
#'            pop_a1=a, pop_d=d)
#' plot(NWU_OEM, type = 'empiricalhist')
#'
#' # No Prior Weights Updating and Multiple EM Cycles (NWU-MEM)
#' NWU_MEM <- fixedCalib(dataset2, model = 1, old_mod = mod, PAU = 'NWU')
#' coef(NWU_MEM, simplify=TRUE)
#' data.frame(coef(NWU_MEM, simplify=TRUE)$items[,c('a1','d')],
#'            pop_a1=a, pop_d=d)
#' plot(NWU_MEM, type = 'empiricalhist')
#'
#' # One Prior Weights Updating and One EM Cycle (OWU-OEM)
#' OWU_OEM <- fixedCalib(dataset2, model=1, old_mod=mod, PAU='OWU', NEMC="OEM")
#' coef(OWU_OEM, simplify=TRUE)
#' data.frame(coef(OWU_OEM, simplify=TRUE)$items[,c('a1','d')], pop_a1=a, pop_d=d)
#' plot(OWU_OEM, type = 'empiricalhist')
#'
#' # One Prior Weights Updating and Multiple EM Cycles (OWU-MEM)
#' OWU_MEM <- fixedCalib(dataset2, model = 1, old_mod = mod, PAU = 'OWU')
#' coef(OWU_MEM, simplify=TRUE)
#' data.frame(coef(OWU_MEM, simplify=TRUE)$items[,c('a1','d')],
#'            pop_a1=a, pop_d=d)
#' plot(OWU_MEM, type = 'empiricalhist')
#'
#' # Multiple Prior Weights Updating and Multiple EM Cycles (MWU-MEM)
#' MWU_MEM <- fixedCalib(dataset2, model = 1, old_mod = mod)
#' coef(MWU_MEM, simplify=TRUE)
#' data.frame(coef(MWU_MEM, simplify=TRUE)$items[,c('a1','d')],
#'            pop_a1=a, pop_d=d)
#' plot(MWU_MEM, type = 'empiricalhist')
#'
#' # factor scores distribution check
#' fs <- fscores(MWU_MEM)
#' hist(fs)
#' c(mean_calib=mean(fs[1:N, ]), sd_calib=sd(fs[1:N, ]))
#' c(mean_exper=mean(fs[-c(1:N), ]), sd_exper=sd(fs[-c(1:N), ]))
#'
#'
#' ############################
#' ## Item length constraint example for each participant in the experimental
#' ## items group. In this example, all participants were forced to have a test
#' ## length of J=30, though the item pool had J=50 total items.
#'
#' # new experimental data (relatively extreme, theta ~ N(.5,1.5))
#' dataset2 <- simdata(a, d, N = 1000, itemtype=itemtype,
#'     mu=.5, sigma=matrix(1.5))
#'
#' # Add missing values to each participant in new dataset where individuals
#' # were randomly administered 10 experimental items, subject to the constraint
#' # that each participant received a test with J=30 items.
#' dataset2 <- t(apply(dataset2, 1, function(x){
#'    NA_precalib <- sample(1:30, 10)
#'    NA_experimental <- sample(31:50, 10)
#'    x[c(NA_precalib, NA_experimental)] <- NA
#'    x
#' }))
#' head(dataset2)
#'
#' # check that all individuals had 30 items
#' all(rowSums(!is.na(dataset2)) == 30)
#'
#' #' Multiple Prior Weights Updating and Multiple EM Cycles (MWU-MEM)
#' MWU_MEM <- fixedCalib(dataset2, model = 1, old_mod = mod)
#' coef(MWU_MEM, simplify=TRUE)
#' data.frame(coef(MWU_MEM, simplify=TRUE)$items[,c('a1','d')],
#'            pop_a1=a, pop_d=d)
#' plot(MWU_MEM, type = 'empiricalhist')
#'
#' ## factor scores check
#' fs <- fscores(MWU_MEM)
#' hist(fs)
#' c(mean_calib=mean(fs[1:N, ]), sd_calib=sd(fs[1:N, ]))
#'
#' ## shrinkage, but generally different from calibrated sample
#' c(mean_exper=mean(fs[-c(1:N), ]), sd_exper=sd(fs[-c(1:N), ]))
#'
#'
#' }
fixedCalib <- function(data, model = 1, old_mod, PAU = 'MWU', NEMC = "MEM",
                       technical = list(), ...){

    stopifnot(PAU %in% c('NWU', 'OWU', 'MWU'))
    stopifnot(NEMC %in% c("OEM", "MEM"))

    global_prior <- NULL
    custom_den_const <- function(obj, Theta) global_prior
    custom_den_EH <- function(obj, Theta) {
        ret <- if(length(obj@rr) == 0L) # first EM iteration use global info
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
    colnames(olddata) <- c(old_itemnames,
                           colnames(data)[!(colnames(data) %in% old_itemnames)])
    fulldata <- rbind(olddata, data)

    sv <- mod2values(old_mod)
    sv$est <- FALSE
    sv2 <- mirt(fulldata, model, pars='values', ...)
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
                           technical = list(NCYCLES = 1L, message = FALSE),
                           ...)
        global_prior <- FC_mod_den@Internals$Prior[[1L]]
        den <- if(PAU == 'OWU') custom_den_const else custom_den_EH
        # names preserved to avoid new starting value issues
        par <- c(MEAN_1 = 0, COV_11 = 1)
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
    FC_mod_den <- mirt(fulldata, model, pars=sv_final,
                       dentype='EH', TOL=NaN, ...)
    FC_mod_den
}