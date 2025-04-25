#' RMSD effect size statistic to quantify category-level DIF
#'
#' This function computes a set of RMSD "badness-of-fit" statistics when investing
#' DIF across a set of grouping variables. In a first step, a (potentially highly constrained)
#' multiple group model is fitted, while in a second step the item (and person) parameters
#' are estimated based on all examines across all groups. Category level DIF is assessed
#' based on how well the pseudo-table of counts match the (constrained) probability functions
#' implied by the original multiple group model (while also weighing across the implied density
#' function of the latent traits). If the RSMD fit is poor, indicating non-ignorable DIF,
#' then the multiple-group model should be adjusted to better account for the large response bias
#' due to using a pooled model. See Lee and von Davier (2020) and Buchholz and Hartig (2019) for details.
#'
#' This function is generally designed to work with multiple-group models, however if a single-group
#' model is passed this will return the RMSD information regarding the observed vs expected response functions.
#'
#' @param pooled_mod a multiple-group model (used to compute the model-implied
#'   probability in the goodness-of-fit test) or a single-group model for computing RMSD goodness-of-fit
#'   values relative to the model-implied information
#' @param flag a numeric value used as a cut-off to help flag larger RMSD values
#'   (e.g., \code{flag = .03} will highlight only categories with RMSD values greater than
#'   .03)
#' @param probfun logical; use probability functions to compute RMSD? If FALSE, the expected score
#'   functions will be integrated instead, which may be useful for collapsing across the
#'   categories in polytomous items
#' @param dentype density to use for the latent trait.
#'   Can be \code{'norm'} to use a normal Gaussian density where the mean/variance are extracted
#'   from the model object(default), \code{'snorm'} for a standard normal distribution,
#'   or \code{'empirical'} to use the density estimate obtained via the E-table
#'
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#'
#' @seealso \code{\link{DIF}}, \code{\link{DRF}}, \code{\link{multipleGroup}}, \code{\link{empirical_ES}}
#'
#' @references
#'
#' Buchholz, J., and Hartig, J. (2019). Comparing Attitudes Across Groups: An IRT-Based Item-Fit Statistic
#'   for the Analysis of Measurement Invariance. \emph{Applied Psychological Measurement, 43}(3), 241-250.
#'   \doi{https://doi.org/10.1177/0146621617748323}
#'
#' Lee, S. S., and von Davier, M. (2020). Improving measurement properties of the PISA home
#'   possessions scale through partial invariance modeling.
#'   \emph{Psychological test and assessment modeling}, 62(1):55-83.
#'
#' @export
#' @examples
#'
#' \donttest{
#'
#' #----- generate some data
#' set.seed(12345)
#' a <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
#'
#' # item 1 has DIF
#' d2[1] <- d[1] - .5
#' a2[1] <- a[1] + 1
#'
#' itemtype <- rep('2PL', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a2, d2, N, itemtype)
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#'
#' #-----
#'
#' # fully pooled model
#' pooled_mod <- multipleGroup(dat, 1, group=group,
#'    invariance = c(colnames(dat), 'free_mean', 'free_var'))
#' coef(pooled_mod, simplify=TRUE)
#'
#' RMSD_DIF(pooled_mod)
#' RMSD_DIF(pooled_mod, dentype = 'empirical')
#' RMSD_DIF(pooled_mod, flag = .03)
#'
#' # more freely estimated model (item 1 has 2 parameters estimated)
#' MGmod <- multipleGroup(dat, 1, group=group,
#'                        invariance = c(colnames(dat)[-1], 'free_mean', 'free_var'))
#' coef(MGmod, simplify=TRUE)
#'
#' # RMSD in item.1 now reduced (MG model accounts for DIF)
#' RMSD_DIF(MGmod)
#' RMSD_DIF(MGmod, flag = .03)
#'
#' #################
#' # NA placeholders included when groups do not respond to specific items
#'
#' a1 <- a2 <- rlnorm(20)
#' d <- d2 <- rnorm(20)
#' # item 5 contains DIF
#' a2[5] <- a1[5] + 1
#' d2[5] <- d[5] - 1/2
#' g <- rbeta(20, 5, 17)
#'
#' dat1 <- simdata(a1, d, guess = g, N=1000, itemtype = '3PL')
#' dat1[, 11:13] <- NA  # items 11:13 items NA for g1
#' dat2 <- simdata(a2, d2, guess = g, N=1000, itemtype = '3PL',
#'    mu=1/4, sigma=matrix(.75))
#' dat2[,1:3] <- NA # items 1:3 items NA for g2
#' dat <- rbind(dat1, dat2)
#' group <- c(rep('g1', 1000), rep('g2', 1000))
#'
#' mod <- multipleGroup(dat, "Theta = 1-20
#'                             PRIOR = (1-20, g, norm, -1, 0.5)",
#'                      group=group, itemtype='3PL',
#'                      invariance = c(colnames(dat), 'free_mean', 'free_var'))
#' coef(mod, simplify = TRUE)
#'
#' RMSD_DIF(mod)
#' RMSD_DIF(mod, flag = .03)
#'
#' #################
#' # Single group model (GOF only)
#' mod <- mirt(dat)
#' RMSD_DIF(mod)
#'
#' #################
#' # polytomous example
#' set.seed(12345)
#' a <- a2 <- matrix(rlnorm(20,.2,.3))
#'
#' # for the graded model, ensure that there is enough space between the intercepts,
#' # otherwise closer categories will not be selected often (minimum distance of 0.3 here)
#' diffs <- t(apply(matrix(runif(20*4, .3, 1), 20), 1, cumsum))
#' diffs <- -(diffs - rowMeans(diffs))
#' d <- d2 <- diffs + rnorm(20)
#'
#' # item 1 has slope + dif for first intercept parameter
#' d2[1] <- d[1] - .5
#' a2[1] <- a[1] + 1
#'
#' itemtype <- rep('graded', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a2, d2, N, itemtype)
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#'
#' #-----
#'
#' # fully pooled model
#' pooled_mod <- multipleGroup(dat, 1, group=group,
#'          invariance = c(colnames(dat), 'free_mean', 'free_var'))
#' coef(pooled_mod, simplify=TRUE)
#'
#' # Item_1 fits poorly in several categories (RMSD > .05)
#' RMSD_DIF(pooled_mod)
#' RMSD_DIF(pooled_mod, flag = .05)
#' RMSD_DIF(pooled_mod, flag = .1, probfun = FALSE) # use expected score function
#'
#' # more freely estimated model (item 1 has more parameters estimated)
#' MGmod <- multipleGroup(dat, 1, group=group,
#'                        invariance = c(colnames(dat)[-1], 'free_mean', 'free_var'))
#' coef(MGmod, simplify=TRUE)
#'
#' # RMSDs in Item_1 now reduced (MG model better accounts for DIF)
#' RMSD_DIF(MGmod)
#' RMSD_DIF(MGmod, flag = .05)
#' RMSD_DIF(MGmod, probfun = FALSE, flag = .1) # use expected score function
#'
#' }
#'
RMSD_DIF <- function(pooled_mod, flag = 0, probfun = TRUE, dentype = 'norm'){
    dat <- extract.mirt(pooled_mod, 'data')
    group <- extract.mirt(pooled_mod, "group")
    which.groups <- extract.mirt(pooled_mod, 'groupNames')
    ret <- vector('list', length(which.groups))
    names(ret) <- which.groups
    for(which.group in which.groups){
        if(is(pooled_mod, 'SingleGroupClass'))
            smod <- pooled_mod else smod <- extract.group(pooled_mod, which.group)
        nfact <- extract.mirt(smod, 'nfact')
        stopifnot(nfact == 1L)
        sv <- mod2values(smod)
        sv$est <- FALSE

        # make wider grid for better numerical integration
        Theta <- matrix(seq(-6,6,length.out = 201))

        tmpdat <- subset(dat, group == which.group)
        tmpdat <- tmpdat[,apply(tmpdat, 2, function(x) all(is.na(x)))]
        mod_g <- mirt(subset(dat, group == which.group), nfact,
                      itemtype = extract.mirt(smod, 'itemtype'),
                      pars = sv, technical = list(storeEtable=TRUE, customTheta=Theta,
                                                  customK=extract.mirt(pooled_mod, 'K')))
        Etable <- mod_g@Internals$Etable[[1]]$r1

        # standard normal dist for theta
        if(dentype %in% c('norm', 'snorm')){
            mu <- sv$value[sv$name == "MEAN_1"]
            sigma2 <- sv$value[sv$name == "COV_11"]
            if(dentype == 'snorm'){
                mu <- 0
                sigma2 <- 1
            }
            f_theta <- dnorm(Theta, mean = mu, sd = sqrt(sigma2))
            f_theta <- as.vector(f_theta / sum(f_theta))
        } else if(dentype == 'empirical'){
            f_theta <- rowSums(Etable) / sum(Etable)
        } else stop('dentype not supported', call.=FALSE)

        itemloc <- extract.mirt(mod_g, 'itemloc')
        which.items <- 1L:ncol(dat)
        ret2 <- vector('list', ncol(dat))
        names(ret2) <- extract.mirt(mod_g, 'itemnames')
        for(i in seq_len(length(which.items))){
            pick <- itemloc[which.items[i]]:(itemloc[which.items[i]+1L] - 1L)
            O <- Etable[ ,pick]
            P_o <- O / rowSums(O)
            item <- extract.item(smod, which.items[i])
            P_e <- probtrace(item, Theta)
            ret2[[i]] <- if(probfun){
                sqrt(colSums((P_o - P_e)^2 * f_theta))
            } else {
                S <- 1L:ncol(P_o)-  1L
                c("S(theta)" = sqrt(sum(( colSums(S*t(P_o - P_e)) )^2 * f_theta)))
            }
        }
        nms <- lapply(ret2, names)
        unms <- unique(do.call(c, nms))
        items <- matrix(NA, length(ret2), length(unms))
        rownames(items) <- names(ret2)
        colnames(items) <- unms
        for(i in seq_len(nrow(items)))
            items[i, nms[[i]]] <- ret2[[i]]
        if(flag > 0)
            items[items < flag] <- NA
        items[is.nan(items)] <- NA
        items <- as.data.frame(items)
        items <- as.mirt_df(items)

        ret[[which.group]] <- items
    }
    if(length(ret) == 1) ret <- ret[[1]]
    ret
}
