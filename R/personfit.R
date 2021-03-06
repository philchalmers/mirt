#' Person fit statistics
#'
#' \code{personfit} calculates the Zh values from Drasgow, Levine and Williams (1985) for
#' unidimensional and multidimensional models, as well as the infit and outfit statistics.
#' The returned object is a \code{data.frame}
#' consisting either of the tabulated data or full data with the statistics appended to the
#' rightmost columns.
#'
#'
#' @aliases personfit
#' @param x a computed model object of class \code{SingleGroupClass} or \code{MultipleGroupClass}
#' @param method type of factor score estimation method. See \code{\link{fscores}} for more detail
#' @param Theta a matrix of factor scores used for statistics that require empirical estimates. If
#'   supplied, arguments typically passed to \code{fscores()} will be ignored and these values will
#'   be used instead
#' @param stats.only logical; return only the person fit statistics without their associated
#'   response pattern?
#' @param ... additional arguments to be passed to \code{fscores()}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords person fit
#' @export personfit
#'
#' @seealso
#' \code{\link{itemfit}}
#'
#' @references
#'
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness measurement with
#' polychotomous item response models and standardized indices.
#' \emph{British Journal of Mathematical and Statistical Psychology, 38}, 67-86.
#'
#' Reise, S. P. (1990). A comparison of item- and person-fit methods of assessing model-data fit
#' in IRT. \emph{Applied Psychological Measurement, 14}, 127-137.
#'
#' Wright B. D. & Masters, G. N. (1982). \emph{Rating scale analysis}. MESA Press.
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' #make some data
#' set.seed(1)
#' a <- matrix(rlnorm(20),ncol=1)
#' d <- matrix(rnorm(20),ncol=1)
#' items <- rep('2PL', 20)
#' data <- simdata(a,d, 2000, items)
#'
#' x <- mirt(data, 1)
#' fit <- personfit(x)
#' head(fit)
#'
#' #using precomputed Theta
#' Theta <- fscores(x, method = 'MAP', full.scores = TRUE)
#' head(personfit(x, Theta=Theta))
#'
#' #multiple group Rasch model example
#' set.seed(12345)
#' a <- matrix(rep(1, 15), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('dich', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#' models <- 'F1 = 1-15'
#' mod_Rasch <- multipleGroup(dat, models, itemtype = 'Rasch', group = group)
#' coef(mod_Rasch, simplify=TRUE)
#' pf <- personfit(mod_Rasch, method='MAP')
#' head(pf)
#'
#'   }
#'
personfit <- function(x, method = 'EAP', Theta = NULL, stats.only = TRUE, ...){
    if(missing(x)) missingMsg('x')
    if(is(x, 'DiscreteClass'))
        stop('Discrete latent structures not yet supported', call.=FALSE)
    if(is(x, 'MixtureClass'))
        stop('Mixture latent structures not yet supported', call.=FALSE)
    if(any(is.na(x@Data$data)))
        stop('Fit statistics cannot be computed when there are missing data.', call.=FALSE)
    if(is(x, 'MultipleGroupClass')){
        ret <- vector('list', x@Data$ngroups)
        if(is.null(Theta))
            Theta <- fscores(x, method=method, full.scores=TRUE, ...)
        for(g in seq_len(x@Data$ngroups)){
            pick <- x@Data$groupNames[g] == x@Data$group
            tmp_obj <- MGC2SC(x, g)
            ret[[g]] <- personfit(tmp_obj, method=method, stats.only=stats.only,
                                  Theta=Theta[pick, , drop=FALSE], ...)
        }
        ret2 <- matrix(0, nrow(x@Data$data), ncol(ret[[1L]]))
        for(g in seq_len(x@Data$ngroups)){
            pick <- x@Data$groupNames[g] == x@Data$group
            ret2[pick, ] <- as.matrix(ret[[g]])
        }
        colnames(ret2) <- colnames(ret[[1L]])
        rownames(ret2) <- rownames(x@Data$data)
        return(as.data.frame(ret2))
    }
    if(is.null(Theta))
        Theta <- fscores(x, verbose=FALSE, full.scores=TRUE, method=method, ...)
    J <- ncol(x@Data$data)
    itemloc <- x@Model$itemloc
    pars <- x@ParObjects$pars
    fulldata <- x@Data$fulldata[[1L]]
    if(nrow(fulldata) < nrow(Theta))
        Theta <- Theta[extract.mirt(x, 'rowID'), , drop=FALSE]
    for(i in seq_len(ncol(Theta))){
        tmp <- Theta[,i]
        tmp[tmp %in% c(-Inf, Inf)] <- NA
        Theta[Theta[,i] == Inf, i] <- max(tmp, na.rm=TRUE) + .1
        Theta[Theta[,i] == -Inf, i] <- min(tmp, na.rm=TRUE) - .1
    }
    N <- nrow(Theta)
    itemtrace <- matrix(0, ncol=ncol(fulldata), nrow=N)
    for (i in seq_len(J))
        itemtrace[ ,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(x=pars[[i]], Theta=Theta)
    LL <- itemtrace * fulldata
    LL[LL < .Machine$double.eps] <- 1
    LL <- rowSums(log(LL))
    Zh <- rep(0, length(LL))
    mu <- sigma2 <- numeric(N)
    log_itemtrace <- log(itemtrace)
    for(item in seq_len(J)){
        P <- itemtrace[ ,itemloc[item]:(itemloc[item+1L]-1L)]
        log_P <- log_itemtrace[ ,itemloc[item]:(itemloc[item+1L]-1L)]
        mu <- mu + rowSums(P * log_P)
        for(i in seq_len(ncol(P)))
            for(j in seq_len(ncol(P)))
                if(i != j)
                    sigma2 <- sigma2 + P[,i] * P[,j] * log_P[,i] * log(P[,i]/P[,j])
    }
    Zh <- (LL - mu) / sqrt(sigma2)
    W <- resid <- C <- matrix(0, ncol=J, nrow=N)
    K <- x@Data$K
    for (i in seq_len(J)){
        P <- ProbTrace(x=pars[[i]], Theta=Theta)
        Emat <- matrix(0:(K[i]-1), nrow(P), ncol(P), byrow = TRUE)
        dat <- fulldata[ ,itemloc[i]:(itemloc[i+1] - 1)]
        item <- extract.item(x, i)
        resid[, i] <- rowSums(dat*Emat) - rowSums(Emat * P)
        W[ ,i] <- rowSums((Emat - rowSums(Emat * P))^2 * P)
        C[ ,i] <- rowSums((Emat - rowSums(Emat * P))^4 * P)
    }
    W[W^2 < 1e-5] <- sqrt(1e-5)
    if(!is.null(attr(x, 'inoutfitreturn'))) return(list(resid=resid, W=W, C=C))
    outfit <- rowSums(resid^2/W) / J
    q.outfit <- sqrt(rowSums((C / W^2) / J^2) - 1 / J)
    q.outfit[q.outfit > 1.4142] <- 1.4142
    z.outfit <- (outfit^(1/3) - 1) * (3/q.outfit) + (q.outfit/3)
    infit <- rowSums(resid^2) / rowSums(W)
    q.infit <- sqrt(rowSums(C - W^2) / rowSums(W)^2)
    q.infit[q.infit > 1.4142] <- 1.4142
    z.infit <- (infit^(1/3) - 1) * (3/q.infit) + (q.infit/3)
    ret <- data.frame(x@Data$data, outfit=outfit, z.outfit=z.outfit,
                      infit=infit, z.infit=z.infit, Zh=Zh)
    rownames(ret) <- NULL
    if(stats.only)
        ret <- ret[, !(colnames(ret) %in% colnames(x@Data$data)), drop=FALSE]
    completely_missing <- extract.mirt(x, 'completely_missing')
    ret <- addMissing(ret, whc=completely_missing)
    return(ret)
}

