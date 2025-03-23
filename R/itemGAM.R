#' Parametric smoothed regression lines for item response probability functions
#'
#' This function uses a generalized additive model (GAM) to estimate response curves for items that
#' do not seem to fit well in a given model. Using a stable axillary model, traceline functions for
#' poorly fitting dichotomous or polytomous items can be inspected using point estimates
#' (or plausible values) of the latent trait. Plots of the tracelines and their associated standard
#' errors are available to help interpret the misfit. This function may also be useful when adding
#' new items to an existing, well established set of items, especially when the parametric form of
#' the items under investigation are unknown.
#'
#' @aliases itemGAM
#' @param item a single poorly fitting item to be investigated. Can be a vector or matrix
#' @param Theta a list or matrix of latent trait estimates typically returned from \code{\link{fscores}}
#' @param formula an R formula to be passed to the \code{gam} function. Default fits a spline model
#'   with 10 nodes. For multidimensional models, the traits are assigned the names 'Theta1', 'Theta2',
#'   ..., 'ThetaN'
#' @param CI a number ranging from 0 to 1 indicating the confidence interval range. Default provides the
#'   95 percent interval
#' @param theta_lim range of latent trait scores to be evaluated
#' @param return.models logical; return a list of GAM models for each category? Useful when the GAMs
#'   should be inspected directly, but also when fitting multidimensional models (this is set to
#'   TRUE automatically for multidimensional models)
#' @param ... additional arguments to be passed to \code{gam} or \code{lattice}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords item fit
#' @seealso \code{\link{itemfit}}
#' @export itemGAM
#' @examples
#'
#' \donttest{
#' set.seed(10)
#' N <- 1000
#' J <- 30
#'
#' a <- matrix(1, J)
#' d <- matrix(rnorm(J))
#' Theta <- matrix(rnorm(N, 0, 1.5))
#' dat <- simdata(a, d, N, itemtype = '2PL', Theta=Theta)
#'
#' # make a bad item
#' ps <- exp(Theta^2 + Theta) / (1 + exp(Theta^2 + Theta))
#' item1 <- sapply(ps, function(x) sample(c(0,1), size = 1, prob = c(1-x, x)))
#'
#' ps2 <- exp(2 * Theta^2 + Theta + .5 * Theta^3) / (1 + exp(2 * Theta^2 + Theta + .5 * Theta^3))
#' item2 <- sapply(ps2, function(x) sample(c(0,1), size = 1, prob = c(1-x, x)))
#'
#' # how the actual item looks in the population
#' plot(Theta, ps, ylim = c(0,1))
#' plot(Theta, ps2, ylim = c(0,1))
#'
#' baditems <- cbind(item1, item2)
#' newdat <- cbind(dat, baditems)
#'
#' badmod <- mirt(newdat, 1)
#' itemfit(badmod) #clearly a bad fit for the last two items
#' mod <- mirt(dat, 1) #fit a model that does not contain the bad items
#' itemfit(mod)
#'
#' #### Pure non-parametric way of investigating the items
#' library(KernSmoothIRT)
#' ks <- ksIRT(newdat, rep(1, ncol(newdat)), 1)
#' plot(ks, item=c(1,31,32))
#' par(ask=FALSE)
#'
#' # Using point estimates from the model
#' Theta <- fscores(mod)
#' IG0 <- itemGAM(dat[,1], Theta) #good item
#' IG1 <- itemGAM(baditems[,1], Theta)
#' IG2 <- itemGAM(baditems[,2], Theta)
#' plot(IG0)
#' plot(IG1)
#' plot(IG2)
#'
#' # same as above, but with plausible values to obtain the standard errors
#' set.seed(4321)
#' ThetaPV <- fscores(mod, plausible.draws=10)
#' IG0 <- itemGAM(dat[,1], ThetaPV) #good item
#' IG1 <- itemGAM(baditems[,1], ThetaPV)
#' IG2 <- itemGAM(baditems[,2], ThetaPV)
#' plot(IG0)
#' plot(IG1)
#' plot(IG2)
#'
#' ## for polytomous test items
#' SAT12[SAT12 == 8] <- NA
#' dat <- key2binary(SAT12,
#'                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' dat <- dat[,-32]
#' mod <- mirt(dat, 1)
#'
#' # Kernal smoothing is very sensitive to which category is selected as 'correct'
#' # 5th category as correct
#' ks <- ksIRT(cbind(dat, SAT12[,32]), c(rep(1, 31), 5), 1)
#' plot(ks, items = c(1,2,32))
#'
#' # 3rd category as correct
#' ks <- ksIRT(cbind(dat, SAT12[,32]), c(rep(1, 31), 3), 1)
#' plot(ks, items = c(1,2,32))
#'
#' # splines approach
#' Theta <- fscores(mod)
#' IG <- itemGAM(SAT12[,32], Theta)
#' plot(IG)
#'
#' set.seed(1423)
#' ThetaPV <- fscores(mod, plausible.draws=10)
#' IG2 <- itemGAM(SAT12[,32], ThetaPV)
#' plot(IG2)
#'
#' # assuming a simple increasing parametric form (like in a standard IRT model)
#' IG3 <- itemGAM(SAT12[,32], Theta, formula = resp ~ Theta)
#' plot(IG3)
#' IG3 <- itemGAM(SAT12[,32], ThetaPV, formula = resp ~ Theta)
#' plot(IG3)
#'
#' ### multidimensional example by returning the GAM objects
#' mod2 <- mirt(dat, 2)
#' Theta <- fscores(mod2)
#' IG4 <- itemGAM(SAT12[,32], Theta, formula = resp ~ s(Theta1, k=10) + s(Theta2, k=10),
#'    return.models=TRUE)
#' names(IG4)
#' plot(IG4[[1L]], main = 'Category 1')
#' plot(IG4[[2L]], main = 'Category 2')
#' plot(IG4[[3L]], main = 'Category 3')
#'
#' }
itemGAM <- function(item, Theta, formula = resp ~ s(Theta, k = 10), CI = .95,
                    theta_lim = c(-3,3), return.models = FALSE, ...){
    if(return.models && is.list(Theta))
        stop('Models will only be returned when Theta is a matrix')
    Theta2 <- seq(theta_lim[1L], theta_lim[2L], length.out=1000)
    z <- qnorm((1 - CI) / 2 + CI)
    keep <- !is.na(item)
    mm <- model.matrix(~ 0 + factor(item))
    if(is.list(Theta)){
        stopifnot(is.matrix(Theta[[1L]]))
        pvfit <- pvfit_se <- vector('list', length(Theta))
        fit <- vector('list', ncol(mm))
        for(j in seq_len(ncol(mm))){
            for(pv in seq_len(length(Theta))){
                tmpdat <- data.frame(mm[,j], Theta[[pv]][keep,])
                colnames(tmpdat) <- c("resp", "Theta")
                out <- gam(formula, tmpdat, family=binomial, ...)
                fitlist <- predict(out, data.frame(Theta=Theta2), se.fit = TRUE) #from mgcv
                pvfit[[pv]] <- fitlist$fit
                pvfit_se[[pv]] <- fitlist$se.fit
            }
            MI <- averageMI(pvfit, pvfit_se)
            fit_high <- MI$par + z * MI$SEpar
            fit_low <- MI$par - z * MI$SEpar
            fit[[j]] <- data.frame(Theta=Theta2, Prob=plogis(MI$par), Prob_low=plogis(fit_low),
                                   Prob_high=plogis(fit_high), stringsAsFactors=FALSE)
        }
        cat <- rep(paste0('cat_', 1:ncol(mm)), each=nrow(fit[[1L]]))
        ret <- cbind(do.call(rbind, fit), cat)
    } else {
        stopifnot(is.matrix(Theta))
        nfact <- ncol(Theta)
        if(nfact > 1L && !return.models){
            return.models <- TRUE
            message('return.models is always set to TRUE for multidimensional models')
        }
        fit <- vector('list', ncol(mm))
        names(fit) <- paste0('cat_', seq_len(ncol(mm)))
        for(j in seq_len(ncol(mm))){
            tmpdat <- data.frame(mm[,j], Theta[keep,])
            if(nfact == 1L){
                colnames(tmpdat) <- c("resp", "Theta")
            } else {
                colnames(tmpdat) <- c("resp", paste0("Theta", seq_len(nfact)))
            }
            out <- gam(formula, tmpdat, family=binomial, ...)
            if(return.models){
                fit[[j]] <- out
            } else {
                fit[[j]] <- predict(out, data.frame(Theta=Theta2), se.fit = TRUE)
            }
        }
        if(return.models) return(fit)
        cat <- rep(names(fit), each=length(fit[[1L]]$fit))
        pred <- do.call(c, lapply(fit, function(x) x$fit))
        se.fit <- do.call(c, lapply(fit, function(x) x$se.fit))
        alphahalf <- (1 - CI)/2
        fit_high <- pred + qnorm(CI + alphahalf) * se.fit
        fit_low <- pred + qnorm(alphahalf) * se.fit
        ret <- data.frame(Theta=Theta2, cat=cat, Prob=plogis(pred),
                          Prob_high=plogis(fit_high), Prob_low=plogis(fit_low),
                          stringsAsFactors = FALSE)
    }
    class(ret) <- 'itemGAM'
    ret
}

#' @rdname itemGAM
#' @method plot itemGAM
#' @param x an object of class 'itemGAM'
#' @param y a \code{NULL} value ignored by the plotting function
#' @param auto.key plotting argument passed to \code{\link[lattice]{lattice}}
#' @param par.strip.text plotting argument passed to \code{\link[lattice]{lattice}}
#' @param par.settings plotting argument passed to \code{\link[lattice]{lattice}}
#' @export
plot.itemGAM <- function(x, y = NULL,
                         par.strip.text = list(cex = 0.7),
                         par.settings = list(strip.background = list(col = '#9ECAE1'),
                                             strip.border = list(col = "black")),
                         auto.key = list(space = 'right', points=FALSE, lines=TRUE), ...){
    class(x) <- 'data.frame'
    if(length(unique(x$cat)) == 2L){
        x <- subset(x, cat == 'cat_2')
        auto.key <- FALSE
    }
    if(!is.list(auto.key)){
        return(xyplot(Prob ~ Theta, data=x, groups = x$cat,
                  upper=x$Prob_high, lower=x$Prob_low, alpha = .2, fill = 'darkgrey',
                  panel = function(x, y, lower, upper, fill, alpha, ...){
                      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                    col = fill, border = FALSE, alpha=alpha, ...)
                      panel.xyplot(x, y, type='l', lty=1, col = 'black', ...)
                  },
                  ylim = c(-0.1,1.1), par.strip.text=par.strip.text,
                  par.settings=par.settings, auto.key=auto.key,
                  ylab = expression(P(theta)), xlab = expression(theta),
                  main = 'GAM item probability curves', ...))
    } else {
        return(xyplot(Prob ~ Theta | cat, data=x, ylim = c(-0.1,1.1), alpha = .2, fill = 'darkgrey',
                      upper=x$Prob_high, lower=x$Prob_low, par.strip.text=par.strip.text,
                      par.settings=par.settings, auto.key=auto.key,
                      panel = function(x, y, lower, upper, fill, alpha, subscripts, ...){
                          panel.polygon(c(x, rev(x)), c(upper[subscripts], rev(lower[subscripts])),
                                        col = fill, border = FALSE, alpha=alpha, ...)
                          panel.xyplot(x, y, type='l', lty=1, col = 'black', ...)
                      },
                      ylab = expression(P(theta)), xlab = expression(theta),
                      main = 'GAM item probability curves', ...))
    }
}
