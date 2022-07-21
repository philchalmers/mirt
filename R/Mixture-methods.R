# Methods
setMethod(
    f = "print",
    signature = signature(x = 'MixtureClass'),
    definition = function(x)
    {
        class(x) <- 'SingleGroupClass'
        print(x)
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'MixtureClass'),
    definition = function(object) {
        print(object)
    }
)

setMethod(
    f = "coef",
    signature = 'MixtureClass',
    definition = function(object, ...)
    {
        dots <- list(...)
        simplify <- if(is.null(dots$simplify)) FALSE else dots$simplify
        ngroups <- object@Data$ngroups
        allPars <- vector('list', ngroups)
        names(allPars) <- object@Data$groupNames
        for(g in 1:ngroups){
            tmp <- object@ParObjects$pars[[g]]
            tmp@Model$lrPars <- object@ParObjects$lrPars
            tmp@Data$data <- object@Data$data[1L, , drop=FALSE]
            allPars[[g]] <- coef(tmp, ...)
            if(simplify)
                allPars[[g]]$class_proportion <- data.frame(pi=extract.mirt(object, 'pis')[g],
                                                            row.names = '')
        }
        return(allPars)
    }
)

setMethod(
    f = "summary",
    signature = signature(object = 'MixtureClass'),
    definition = function(object, rotate = 'oblimin', verbose = TRUE, ...) {
        ngroups <- object@Data$ngroups
        ret <- list()
        for(g in 1:ngroups){
            if(verbose) cat('\n----------\nGROUP:', as.character(object@Data$groupNames[g]), '\n')
            ret[[g]] <- summary(object@ParObjects$pars[[g]], verbose=verbose,
                                rotate = rotate, ...)
            ret[[g]]$class_proportion <- object@Model$pis[g]
            if(verbose)
                cat("\nClass proportion: ", round(object@Model$pis[g], 3), "\n")
        }
        invisible(ret)
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'MixtureClass'),
    definition = function(object, object2, ...)
    {
        class(object) <- 'SingleGroupClass'
        anova(object, object2, ..., frame = 2)
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'MixtureClass', y = 'missing'),
    definition = function(x, collapse = FALSE, ...)
    {
        if(collapse){
            ret <- plot_mixture(x, ...)
        } else {
            class(x) <- 'MultipleGroupClass'
            ret <- plot(x, ...)
        }
        ret
    }
)

plot_mixture <- function(x, y, type = 'score', npts = 200, degrees = 45,
                         theta_lim = c(-6,6), which.items = 1:extract.mirt(x, 'nitems'),
                         MI = 0, CI = .95, rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                         facet_items = TRUE, main = NULL,
                         drape = TRUE, colorkey = TRUE, ehist.cut = 1e-10, add.ylab2 = TRUE,
                         par.strip.text = list(cex = 0.7),
                         par.settings = list(strip.background = list(col = '#9ECAE1'),
                                             strip.border = list(col = "black")),
                         auto.key = list(space = 'right', points=FALSE, lines=TRUE),
                         profile = FALSE, ...)
{
    dots <- list(...)
    pis <- extract.mirt(x, 'pis')
    if(!(type %in% c('info', 'SE', 'infoSE', 'trace', 'score', 'itemscore',
                     'infocontour', 'infotrace', 'scorecontour')))
        stop('type supplied is not supported')
    if (any(degrees > 90 | degrees < 0))
        stop('Improper angle specified. Must be between 0 and 90.', call.=FALSE)
    rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
    nfact <- x@Model$nfact
    if(nfact > 3) stop("Can't plot high dimensional solutions.", call.=FALSE)
    J <- x@Data$nitems
    theta <- seq(theta_lim[1L],theta_lim[2L],length.out=npts/(nfact^2))
    ThetaFull <- Theta <- thetaComb(theta, nfact)
    prodlist <- attr(x@ParObjects$pars, 'prodlist')
    if(all(x@Data$K[which.items] == 2L) && facet_items) auto.key <- FALSE
    if(length(prodlist))
        ThetaFull <- prodterms(Theta,prodlist)
    if(length(degrees) == 1L) degrees <- rep(degrees, ncol(ThetaFull))
    info <- numeric(nrow(ThetaFull))
    if(type %in% c('info', 'infocontour', 'rxx', 'SE', 'infoSE', 'infotrace')){
        info <- pis[1L] * testinfo(x@ParObjects$pars[[1L]], ThetaFull, degrees = degrees, which.items=which.items)
        for(g in 2L:length(pis))
            info <- pis[g] * testinfo(x@ParObjects$pars[[g]], ThetaFull, degrees = degrees, which.items=which.items)
    }
    mins <- x@Data$mins
    maxs <- extract.mirt(x, 'K') + mins - 1
    rotate <- if(is.null(dots$rotate)) 'none' else dots$rotate
    if (x@Options$exploratory){
        if(!is.null(dots$rotate)){
            so <- summary(x, verbose=FALSE, digits=5, ...)
            a <- rotateLambdas(so) * 1.702
            for(i in 1:J)
                x@ParObjects$pars[[i]]@par[1:nfact] <- a[i, ]
        }
    }
    itemtrace <- computeItemtrace(x@ParObjects$pars, ThetaFull, x@Model$itemloc,
                                  CUSTOM.IND=x@Internals$CUSTOM.IND, pis=pis)
    score <- c()
    for(i in 1:J)
        score <- c(score, (0:(x@Data$K[i]-1) + mins[i]) * (i %in% which.items))
    score <- matrix(score, nrow(itemtrace), ncol(itemtrace), byrow = TRUE)
    plt <- data.frame(cbind(info,score=rowSums(score*itemtrace),Theta=Theta))
    bundle <- length(which.items) != J
    mins <- mins[which.items]
    maxs <- maxs[which.items]
    ybump <- (max(maxs) - min(mins))/15
    ybump_full <- (sum(maxs) - sum(mins))/15
    if(nfact == 3){
        colnames(plt) <- c("info", "score", "Theta1", "Theta2", "Theta3")
        plt$SE <- 1 / sqrt(plt$info)
        if(type == 'infocontour'){
            if(is.null(main)){
                main <- paste("Test Information Contour")
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(contourplot(info ~ Theta1 * Theta2 | Theta3, data = plt,
                               main = main, xlab = expression(theta[1]),
                               ylab = expression(theta[2]),
                               par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'scorecontour'){
            if(is.null(main)){
                main <- paste("Expected Score Contour")
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(contourplot(score ~ Theta1 * Theta2 | Theta3, data = plt,
                               ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                               main = main, xlab = expression(theta[1]),
                               ylab = expression(theta[2]), ylim=c(sum(mins)-.1, sum(maxs)+.1),
                               par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'info'){
            if(is.null(main)){
                main <- "Test Information"
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(wireframe(info ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                             zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                             scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                             par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'SEcontour'){
            if(is.null(main)){
                main <- "Test Standard Errors"
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(contourplot(score ~ Theta1 * Theta2 | Theta3, data = plt,
                               main = main, xlab = expression(theta[1]),
                               ylab = expression(theta[2]),
                               par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'score'){
            if(is.null(main)){
                main <- if(bundle) "Expected Bundle Score" else "Expected Total Score"
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(wireframe(score ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                             ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                             zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                             scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                             par.strip.text=par.strip.text, par.settings=par.settings,
                             ylim=c(sum(mins)-.1, sum(maxs)+.1), ...))
        } else if(type == 'SE'){
            if(is.null(main)){
                main <- "Test Standard Errors"
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(wireframe(SE ~ Theta1 + Theta2 | Theta3, data = plt, main = main,
                             zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                             scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                             par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else {
            stop('plot type not supported for three dimensional model', call.=FALSE)
        }
    } else if(nfact == 2){
        colnames(plt) <- c("info", "score", "Theta1", "Theta2")
        plt$SE <- 1 / sqrt(plt$info)
        if(type == 'infocontour'){
            if(is.null(main)){
                main <- paste("Test Information Contour")
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(contourplot(info ~ Theta1 * Theta2, data = plt,
                               main = main, xlab = expression(theta[1]),
                               ylab = expression(theta[2]),
                               par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'scorecontour'){
            if(is.null(main)){
                main <- paste("Expected Score Contour")
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(contourplot(score ~ Theta1 * Theta2, data = plt,
                               ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                               main = main, xlab = expression(theta[1]),
                               ylab = expression(theta[2]),
                               par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'info'){
            if(is.null(main)){
                main <- "Test Information"
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(wireframe(info ~ Theta1 + Theta2, data = plt, main = main,
                             zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                             scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                             par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'SEcontour'){
            if(is.null(main)){
                main <- "Test Standard Errors"
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(contourplot(score ~ Theta1 * Theta2, data = plt,
                               main = main, xlab = expression(theta[1]),
                               ylab = expression(theta[2]),
                               par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'score'){
            if(is.null(main)){
                main <- if(bundle) "Expected Bundle Score" else "Expected Total Score"
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(wireframe(score ~ Theta1 + Theta2, data = plt, main = main,
                             zlim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                             zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                             scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                             par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else if(type == 'infoangle'){
            if(is.null(main)){
                main <- 'Information across different angles'
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            graphics::symbols(plt$Theta1, plt$Theta2, circles = sqrt(plt$info/pi), inches = .35, fg='white', bg='blue',
                              xlab = expression(theta[1]), ylab = expression(theta[2]),
                              main = main)
        } else if(type == 'SE'){
            if(is.null(main)){
                main <- "Test Standard Errors"
                if(x@Options$exploratory) main <- paste0(main, ' (rotate = \'', rotate, '\')')
            }
            return(wireframe(SE ~ Theta1 + Theta2, data = plt, main = main,
                             zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                             scales = list(arrows = FALSE), screen = rot, colorkey = colorkey, drape = drape,
                             par.strip.text=par.strip.text, par.settings=par.settings, ...))
        } else {
            stop('plot type not supported for two dimensional model', call.=FALSE)
        }
    } else {
        colnames(plt) <- c("info", "score", "Theta")
        plt$SE <- 1 / sqrt(plt$info)
        if(type == 'info'){
            if(is.null(main))
                main <- 'Test Information'
            if(MI > 0){
                return(xyplot(info ~ Theta, data=plt,
                              upper=plt$CIinfoupper, lower=plt$CIinfolower,
                              panel = function(x, y, lower, upper, ...){
                                  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                col="#E6E6E6", border = FALSE, ...)
                                  panel.xyplot(x, y, type='l', lty=1,...)
                              },
                              main = main, ylim=c(min(plt$CIinfolower), max(plt$CIinfoupper)),
                              ylab = expression(I(theta)), xlab = expression(theta),
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                return(xyplot(info~Theta, plt, type='l', main = main,
                              xlab = expression(theta), ylab=expression(I(theta)),
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            }
        } else if(type == 'rxx'){
            if(is.null(main))
                main <- 'Reliability'
            if(MI > 0){
                return(xyplot(rxx ~ Theta, data=plt,
                              upper=plt$CIrxxupper, lower=plt$CIrxxlower,
                              panel = function(x, y, lower, upper, ...){
                                  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                col="#E6E6E6", border = FALSE, ...)
                                  panel.xyplot(x, y, type='l', lty=1,...)
                              },
                              main = main, ylim=c(-0.1, 1.1),
                              ylab = expression(r[xx](theta)), xlab = expression(theta),
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                return(xyplot(rxx~Theta, plt, type='l', main = main, ylim=c(-0.1, 1.1),
                              xlab = expression(theta), ylab=expression(r[xx](theta)),
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            }
        } else if(type == 'SE'){
            if(is.null(main))
                main <- 'Test Standard Errors'
            if(MI > 0){
                return(xyplot(SE ~ Theta, data=plt,
                              upper=plt$CISEupper, lower=plt$CISElower,
                              panel = function(x, y, lower, upper, ...){
                                  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                col="#E6E6E6", border = FALSE, ...)
                                  panel.xyplot(x, y, type='l', lty=1,...)
                              },
                              main = main, ylim=c(min(plt$CISElower), max(plt$CISEupper)),
                              ylab = expression(I(theta)), xlab = expression(theta),
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                return(xyplot(SE~Theta, plt, type='l', main = main,
                              xlab = expression(theta), ylab=expression(SE(theta)),
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            }
        } else if(type == 'infoSE'){
            if(is.null(main))
                main <- 'Test Information and Standard Errors'
            obj1 <- xyplot(info~Theta, plt, type='l', main = main,
                           xlab = expression(theta), ylab=expression(I(theta)),
                           par.strip.text=par.strip.text, par.settings=par.settings)
            obj2 <- xyplot(SE~Theta, plt, type='l', ylab=expression(SE(theta)),
                           par.strip.text=par.strip.text, par.settings=par.settings)
            if(requireNamespace("latticeExtra", quietly = TRUE)){
                return(latticeExtra::doubleYScale(obj1, obj2, add.ylab2 = add.ylab2, ...))
            } else {
                stop('latticeExtra package is not available. Please install.', call.=FALSE)
            }
        } else if(type == 'trace'){
            if(is.null(main))
                main <- 'Item trace lines'
            P <- vector('list', length(which.items))
            names(P) <- colnames(x@Data$data)[which.items]
            ind <- 1L
            for(i in which.items){
                tmp <- probtrace(extract.item(x, i), ThetaFull)
                if(ncol(tmp) == 2L) tmp <- tmp[,2, drop=FALSE]
                tmp2 <- data.frame(P=as.numeric(tmp), cat=gl(ncol(tmp), k=nrow(Theta),
                                                             labels=paste0('P', seq_len(ncol(tmp)))))
                P[[ind]] <- tmp2
                ind <- ind + 1L
            }
            nrs <- sapply(P, nrow)
            Pstack <- do.call(rbind, P)
            names <- c()
            for(i in seq_len(length(nrs)))
                names <- c(names, rep(names(P)[i], nrs[i]))
            plotobj <- data.frame(Pstack, item=names, Theta=Theta)
            plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
            if(facet_items){
                return(xyplot(P ~ Theta|item, plotobj, ylim = c(-0.1,1.1), groups = cat,
                              xlab = expression(theta), ylab = expression(P(theta)),
                              auto.key = auto.key, type = 'l', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                return(xyplot(P ~ Theta, plotobj, groups=plotobj$item, ylim = c(-0.1,1.1),
                              xlab = expression(theta), ylab = expression(P(theta)),
                              auto.key = auto.key, type = 'l', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            }
        } else if(type == 'itemscore'){
            if(is.null(main))
                main <- 'Expected item scoring function'
            S <- vector('list', length(which.items))
            names(S) <- colnames(x@Data$data)[which.items]
            ind <- 1L
            for(i in which.items){
                S[[ind]] <- pis[1L] * expected.item(extract.item(x, i, group = 1), ThetaFull, mins[i])
                for(g in 2L:length(pis))
                    S[[ind]] <- S[[ind]] + pis[g] * expected.item(extract.item(x, i, group = 1), ThetaFull, mins[i])
                ind <- ind + 1L
            }
            Sstack <- do.call(c, S)
            names <- rep(names(S), each = nrow(ThetaFull))
            plotobj <- data.frame(S=Sstack, item=names, Theta=Theta)
            plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
            if(facet_items){
                return(xyplot(S ~ Theta|item, plotobj, ylim=c(min(mins)-ybump, max(maxs)+ybump),
                              xlab = expression(theta), ylab = expression(S(theta)),
                              auto.key = auto.key, type = 'l', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                return(xyplot(S ~ Theta, plotobj, groups=plotobj$item, ylim=c(min(mins)-.1, max(maxs)+.1),
                              xlab = expression(theta), ylab = expression(S(theta)),
                              auto.key = auto.key, type = 'l', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            }
        } else if(type == 'infotrace'){
            if(is.null(main))
                main <- 'Item information trace lines'
            I <- matrix(NA, nrow(Theta), J)
            for(i in which.items)
                I[,i] <- iteminfo(extract.item(x, i), ThetaFull)
            I <- t(na.omit(t(I)))
            items <- rep(colnames(x@Data$data)[which.items], each=nrow(Theta))
            plotobj <- data.frame(I = as.numeric(I), Theta=Theta, item=items)
            plotobj$item <- factor(plotobj$item, levels = colnames(x@Data$data)[which.items])
            if(facet_items){
                return(xyplot(I ~ Theta|item, plotobj,
                              xlab = expression(theta), ylab = expression(I(theta)),
                              auto.key = auto.key, type = 'l', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                return(xyplot(I ~ Theta, plotobj, groups = plotobj$item,
                              xlab = expression(theta), ylab = expression(I(theta)),
                              auto.key = auto.key, type = 'l', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            }
        } else if(type == 'score'){
            if(is.null(main))
                main <- if(bundle) "Expected Bundle Score" else "Expected Total Score"
            if(MI > 0){
                return(xyplot(score ~ Theta, data=plt,
                              ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                              upper=plt$CIscoreupper, lower=plt$CIscorelower,
                              panel = function(x, y, lower, upper, ...){
                                  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                col="#E6E6E6", border = FALSE, ...)
                                  panel.xyplot(x, y, type='l', lty=1,...)
                              },
                              main = main,
                              ylab = expression(T(theta)), xlab = expression(theta),
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            } else {
                return(xyplot(score ~ Theta, plt, ylim=c(sum(mins)-ybump_full, sum(maxs)+ybump_full),
                              xlab = expression(theta), ylab = expression(T(theta)),
                              type = 'l', main = main,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            }
        } else {
            stop('plot not supported for unidimensional models', call.=FALSE)
        }
    }
}


setMethod(
    f = "residuals",
    signature = signature(object = 'MixtureClass'),
    definition = function(object, ...)
    {
        stop('residuals() not supported for mixture IRT models yet', call.=FALSE)
        pis <- ExtractMixtures(object)
        class(object) <- 'SingleGroupClass'
        residuals(object, pis=pis, mixture=TRUE, ...)
    }
)

# Methods
setMethod(
    f = "vcov",
    signature = signature(object = 'MixtureClass'),
    definition = function(object)
    {
        class(object) <- 'SingleGroupClass'
        vcov(object)
    }
)

setMethod(
    f = "logLik",
    signature = signature(object = 'MixtureClass'),
    definition = function(object){
        extract.mirt(object, 'logLik')
    }
)
