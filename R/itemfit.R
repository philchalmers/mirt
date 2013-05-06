#' Item fit statistics
#'
#' \code{itemfit} calculates the Zh values from Drasgow, Levine and Williams (1985),
#' \eqn{\chi^2} values for unidimensional models, and S-X2 statistics for unidimensional models
#' (Kang & Chen, 2007; Orlando & Thissen, 2000). For Rasch, partial credit, and rating scale models
#' infit and outfit statistics are also produced.
#'
#' @aliases itemfit
#' @param x a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
#' \code{MultipleGroupClass}
#' @param Zh logical; calculate Zh and associated statistics (infit/outfit)? Disable this is you are only
#' interested in computing the S-X2 quickly
#' @param X2 logical; calculate the X2 statistic for unidimensional models?
#' @param mincell the minimum expected cell size to be used in the S-X2 computations. Tables will be
#' collapsed across items first if polytomous, and then across scores if necessary
#' @param S_X2.tables logical; return the tables in a list format used to compute the S-X2 stats?
#' @param group.size approximate size of each group to be used in calculating the \eqn{\chi^2} statistic
#' @param empirical.plot a single numeric value or character of the item name  indicating which item to plot
#'  (via \code{itemplot}) and
#' overlay with the empirical \eqn{\theta} groupings. Only applicable when \code{type = 'X2'}.
#' The default is \code{NULL}, therefore no plots are drawn
#' @param method type of factor score estimation method. Can be expected
#' a-posteriori (\code{"EAP"}), Bayes modal (\code{"MAP"}), weighted likelihood estimation
#' (\code{"WLE"}), or maximum likelihood (\code{"ML"})
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords item fit
#' @export itemfit
#'
#' @seealso
#' \code{\link{personfit}}
#'
#' @references
#'
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with
#' polychotomous item response models and standardized indices.
#' \emph{Journal of Mathematical and Statistical Psychology, 38}, 67-86.
#'
#' Kang, T. & Chen, Troy, T. (2007). An investigation of the performance of the generalized
#' S-X2 item-fit index for polytomous IRT models. ACT
#'
#' Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for dichotomous item response theory
#' models. \emph{Applied Psychological Measurement, 24}, 50-64.
#'
#' Reise, S. P. (1990). A comparison of item- and person-fit methods of assessing model-data fit
#' in IRT. \emph{Applied Psychological Measurement, 14}, 127-137.
#'
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' #make some data
#' set.seed(1234)
#' a <- matrix(rlnorm(20, meanlog=0, sdlog = .1),ncol=1)
#' d <- matrix(rnorm(20),ncol=1)
#' items <- rep('dich', 20)
#' data <- simdata(a,d, 2000, items)
#'
#' x <- mirt(data, 1)
#' raschfit <- mirt(data, 1, itemtype='Rasch')
#' fit <- itemfit(x)
#' fit
#'
#' itemfit(x, empirical.plot = 1) #empirical item plot
#' itemfit(raschfit, method = 'ML') #infit and outfit stats (method='ML' agrees better with eRm package)
#'
#' #similar example to Kang and Chen 2007
#' a <- matrix(c(.8,.4,.7, .8, .4, .7, 1, 1, 1, 1))
#' d <- matrix(rep(c(2.0,0.0,-1,-1.5),10), ncol=4, byrow=TRUE)
#' dat <- simdata(a,d,2000, itemtype = rep('graded', 10)) - 1
#' head(dat)
#'
#' mod <- mirt(dat, 1)
#' itemfit(mod)
#'
#' mod2 <- mirt(dat, 1, 'Rasch')
#' itemfit(mod2)
#'
#' #massive list of tables
#' tables <- itemfit(mod, S_X2.tables = TRUE)
#'
#' #observed and expected total score patterns for item 1 (post collapsing)
#' tables$O[[1]]
#' tables$E[[1]]
#'
#'
#'   }
#'
itemfit <- function(x, Zh = TRUE, X2 = FALSE, group.size = 150, mincell = 1, S_X2.tables = FALSE,
                    empirical.plot = NULL, method = 'EAP'){
    if(any(is.na(x@data)))
        stop('Fit statistics cannot be computed when there are missing data.')
    if(is(x, 'MultipleGroupClass')){
        ret <- list()
        for(g in 1:length(x@cmods)){
            x@cmods[[g]]@itemtype <- x@itemtype
            ret[[g]] <- itemfit(x@cmods[[g]], group.size=group.size, mincell = 1,
                                S_X2.tables = FALSE)
        }
        names(ret) <- x@groupNames
        return(ret)
    }
    if(S_X2.tables) Zh <- X2 <- FALSE
    ret <- data.frame(item=colnames(x@data))
    J <- ncol(x@data)
    itemloc <- x@itemloc
    pars <- x@pars
    if(Zh || X2){
        sc <- fscores(x, verbose = FALSE, full.scores = TRUE)
        prodlist <- attr(pars, 'prodlist')
        nfact <- x@nfact + length(prodlist)
        fulldata <- x@fulldata
        Theta <- sc[ ,ncol(sc):(ncol(sc) - nfact + 1), drop = FALSE]
        N <- nrow(Theta)
        itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)
        for (i in 1:J)
            itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)
        LL <- itemtrace * fulldata
        LL[LL < .Machine$double.eps] <- 1
        Lmatrix <- matrix(log(LL[as.logical(fulldata)]), N, J)
        mu <- sigma2 <- rep(0, J)
        log_itemtrace <- log(itemtrace)
        for(item in 1:J){
            P <- itemtrace[ ,itemloc[item]:(itemloc[item+1]-1)]
            log_P <- log_itemtrace[ ,itemloc[item]:(itemloc[item+1]-1)]
            mu[item] <- sum(P * log_P)
            for(i in 1:ncol(P))
                for(j in 1:ncol(P))
                    if(i != j)
                        sigma2[item] <- sigma2[item] + sum(P[,i] * P[,j] * log_P[,i] * log(P[,i]/P[,j]))
        }
        Zh <- (colSums(Lmatrix) - mu) / sqrt(sigma2)
        #if all Rasch models, infit and outfit
        if(all(x@itemtype %in% c('Rasch', 'rsm', 'gpcm'))){
            oneslopes <- rep(FALSE, length(x@itemtype))
            for(i in 1:length(x@itemtype))
                oneslopes[i] <- closeEnough((x@pars[[i]]@par[1] * x@pars[[1]]@D), 1-1e-10, 1+1e-10)
            if(all(oneslopes)){
                attr(x, 'inoutfitreturn') <- TRUE
                pf <- personfit(x, method=method)
                z2 <- pf$resid^2 / pf$W
                outfit <- colSums(z2) / N
                q.outfit <- sqrt(colSums((pf$C / pf$W^2) / N^2) - 1 / N)
                z.outfit <- (outfit^(1/3) - 1) * (3/q.outfit) + (q.outfit/3)
                infit <- colSums(pf$W * z2) / colSums(pf$W)
                q.infit <- sqrt(colSums(pf$C - pf$W^2) / colSums(pf$W)^2)
                z.infit <- (infit^(1/3) - 1) * (3/q.infit) + (q.infit/3)
                ret$outfit <- outfit
                ret$z.outfit <- z.outfit
                ret$infit <- infit
                ret$z.infit <- z.infit
                ret$Zh <- Zh
            } else {
                ret$Zh <- Zh
            }
        }
    }
    if((X2 || !is.null(empirical.plot)) && nfact == 1){
        ord <- order(Theta[,1])
        fulldata <- fulldata[ord,]
        Theta <- Theta[ord, , drop = FALSE]
        den <- dnorm(Theta, 0, .5)
        den <- den / sum(den)
        cumTheta <- cumsum(den)
        Groups <- rep(20, length(ord))
        ngroups <- ceiling(nrow(fulldata) / group.size)
        weight <- 1/ngroups
        for(i in 1:20)
            Groups[round(cumTheta,2) >= weight*(i-1) & round(cumTheta,2) < weight*i] <- i
        n.uniqueGroups <- length(unique(Groups))
        X2 <- df <- RMSEA <- rep(0, J)
        if(!is.null(empirical.plot)){
            if(nfact > 1) stop('Cannot make empirical plot for multidimensional models')
            theta <- seq(-4,4, length.out=40)
            ThetaFull <- thetaComb(theta, nfact)
            if(!is.numeric(empirical.plot)){
                inames <- colnames(x@data)
                ind <- 1:length(inames)
                empirical.plot <- ind[inames == empirical.plot]
            }
            P <- ProbTrace(pars[[empirical.plot]], ThetaFull)
            plot(ThetaFull, P[,1], type = 'l', ylim = c(0,1), las = 1,
                 main =paste('Item', empirical.plot), ylab = expression(P(theta)), xlab = expression(theta))
            for(i in 2:ncol(P))
                lines(ThetaFull, P[,i], col = i)
        }
        for (i in 1:J){
            if(!is.null(empirical.plot) && i != empirical.plot) next
            for(j in 1:n.uniqueGroups){
                dat <- fulldata[Groups == j, itemloc[i]:(itemloc[i+1] - 1), drop = FALSE]
                r <- colSums(dat)
                N <- nrow(dat)
                mtheta <- matrix(mean(Theta[Groups == j,]), nrow=1)
                if(!is.null(empirical.plot)){
                    tmp <- r/N
                    col <- 2:length(tmp)
                    points(rep(mtheta, length(tmp) - 1), tmp[-1], col = col, pch = col+2)
                }
                P <- ProbTrace(x=pars[[i]], Theta=mtheta)
                if(any(N * P < 2)){
                    df[i] <- df[i] - 1
                    next
                }
                X2[i] <- X2[i] + sum((r - N*P)^2 / N*P)
            }
            df[i] <- df[i] + n.uniqueGroups*(length(r) - 1) - sum(pars[[i]]@est)
        }
        X2[X2 == 0] <- NA
        if(!is.null(empirical.plot)) return(invisible(NULL))
        ret$X2 <- X2
        ret$df <- df
        ret$p.X2 <- round(1 - pchisq(X2, df), 4)
    }
    makeObstables <- function(dat, K){
        ret <- vector('list', ncol(dat))
        sumscore <- rowSums(dat)
        for(i in 1:length(ret)){
            ret[[i]] <- matrix(0, sum(K-1)+1, K[i])
            colnames(ret[[i]]) <- paste0(1:K[i]-1)
            rownames(ret[[i]]) <- paste0(1:nrow(ret[[i]])-1)
            split <- by(sumscore, dat[,i], table)
            for(j in 1:K[i]){
                m <- match(names(split[[j]]), rownames(ret[[i]]))
                ret[[i]][m,j] <- split[[j]]
            }
            ret[[i]] <- ret[[i]][-c(1, nrow(ret[[i]])), ]
        }
        ret
    }
    collapseCells <- function(O, E, mincell = 1){
        for(i in 1:length(O)){
            On <- O[[i]]
            En <- E[[i]]
            drop <- which(rowSums(is.na(En)) > 0)
            En[is.na(En)] <- 0
            #collapse known upper and lower sparce cells
            if(length(drop) > 0){
                up <- drop[1]:drop[length(drop)/2]
                low <- drop[length(drop)/2 + 1]:drop[length(drop)]
                En[max(up)+1, ] <- colSums(En[c(up, max(up)+1), , drop = FALSE])
                On[max(up)+1, ] <- colSums(On[c(up, max(up)+1), , drop = FALSE])
                En[min(low)-1, ] <- colSums(En[c(low, min(low)-1), , drop = FALSE])
                On[min(low)-1, ] <- colSums(On[c(low, min(low)-1), , drop = FALSE])
                En[c(up, low), ] <- On[c(up, low), ] <- NA
                En <- na.omit(En)
                On <- na.omit(On)
            }
            #collapse accross
            if(ncol(En) > 2){
                for(j in 1:(ncol(On)-1)){
                    L <- En < mincell
                    sel <- L[,j]
                    if(!any(sel)) next
                    On[sel, j+1]  <- On[sel, j] + On[sel, j+1]
                    En[sel, j+1]  <- En[sel, j] + En[sel, j+1]
                    On[sel, j] <- En[sel, j] <- NA
                }
                sel <- L[,j+1]
                sel[rowSums(is.na(En[, 1:j])) == (ncol(En)-1)] <- FALSE
                put <- apply(En[sel, 1:j], 1, function(x) max(which(!is.na(x))))
                put2 <- which(sel)
                for(k in 1:length(put)){
                    En[put2[k], put[k]] <- En[put2[k], put[k]] + En[put2[k], j+1]
                    En[put2[k], j+1] <- On[put2[k], j+1] <- NA
                }
            }
            L <- En < mincell
            L[is.na(L)] <- FALSE
            while(any(L)){
                drop <- c()
                for(j in 1:(nrow(On)-1)){
                    if(any(L[j,])) {
                        On[j+1, L[j,]] <- On[j+1, L[j,]] + On[j, L[j,]]
                        En[j+1, L[j,]] <- En[j+1, L[j,]] + En[j, L[j,]]
                        drop <- c(drop, j)
                        break
                    }
                }
                for(j in nrow(On):2){
                    if(any(L[j,])) {
                        On[j-1, L[j,]] <- On[j-1, L[j,]] + On[j, L[j,]]
                        En[j-1, L[j,]] <- En[j-1, L[j,]] + En[j, L[j,]]
                        drop <- c(drop, j)
                        break
                    }
                }
                if(nrow(On) > 4){
                    for(j in 2:(nrow(On)-1)){
                        if(any(L[j,])){
                            On[j+1, L[j,]] <- On[j+1, L[j,]] + On[j, L[j,]]
                            En[j+1, L[j,]] <- En[j+1, L[j,]] + En[j, L[j,]]
                            drop <- c(drop, j)
                            break
                        }
                    }
                }
                #drop
                if(!is.null(drop)){
                    En <- En[-drop, ]
                    On <- On[-drop, ]
                }
                L <- En < mincell
                L[is.na(L)] <- FALSE
            }
            E[[i]] <- En
            O[[i]] <- On
        }
        return(list(O=O, E=E))
    }
    if(x@nfact == 1){
        dat <- x@data
        adj <- apply(dat, 2, min)
        if(any(adj > 0))
            message('Data adjusted so that the lowest category score for every item is 0')
        dat <- t(t(dat) - adj)
        S_X2 <- df.S_X2 <- numeric(J)
        O <- makeObstables(dat, x@K)
        Nk <- rowSums(O[[1]])
        E <- EAPsum(x, S_X2 = TRUE, gp = list(gmeans=0, gcov=matrix(1)))
        for(i in 1:J)
            E[[i]] <- E[[i]] * Nk
        coll <- collapseCells(O, E, mincell=mincell)
        if(S_X2.tables) return(list(O.org=O, E.org=E, O=coll$O, E=coll$E))
        O <- coll$O
        E <- coll$E
        for(i in 1:J){
            S_X2[i] <- sum((O[[i]] - E[[i]])^2 / E[[i]], na.rm = TRUE)
            df.S_X2[i] <- (ncol(O[[i]])-1) * nrow(O[[i]]) - sum(pars[[i]]@est) - sum(is.na(E[[i]]))
        }
        ret$S_X2 <- S_X2
        ret$df.S_X2 <- df.S_X2
        ret$p.S_X2 <- round(1 - pchisq(S_X2, df.S_X2), 4)
    }
    return(ret)
}
