# Methods
setMethod(
    f = "print",
    signature = signature(x = 'MultipleGroupClass'),
    definition = function(x)
    {
        class(x) <- 'SingleGroupClass'
        print(x)
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object) {
        print(object)
    }
)

setMethod(
    f = "coef",
    signature = 'MultipleGroupClass',
    definition = function(object, ...)
    {
        ngroups <- object@Data$ngroups
        allPars <- vector('list', ngroups)
        names(allPars) <- object@Data$groupNames
        for(g in 1:ngroups){
            tmp <- object@ParObjects$pars[[g]]
            tmp@Data$data <- object@Data$data[1L, , drop=FALSE]
            allPars[[g]] <- coef(tmp, ...)
        }
        return(allPars)
    }
)

setMethod(
    f = "summary",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, digits = 3, rotate = 'oblimin', verbose = TRUE, ...) {
        ngroups <- object@Data$ngroups
        ret <- list()
        for(g in 1:ngroups){
            if(verbose) cat('\n----------\nGROUP:', as.character(object@Data$groupNames[g]), '\n')
            ret[[g]] <- summary(object@ParObjects$pars[[g]], digits=digits, verbose=verbose,
                                rotate = rotate, ...)
        }
        invisible(ret)
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, object2, ...)
    {
        class(object) <- 'SingleGroupClass'
        anova(object, object2, ...)
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'MultipleGroupClass', y = 'missing'),
    definition = function(x, y, type = 'score', npts = 50, degrees = 45,
                          which.items = 1:extract.mirt(x, 'nitems'),
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                          facet_items = TRUE,
                          theta_lim = c(-6,6),
                          par.strip.text = list(cex = 0.7),
                          par.settings = list(strip.background = list(col = '#9ECAE1'),
                                              strip.border = list(col = "black")),
                          auto.key = list(space = 'right'), ...)
    {
        if (!type %in% c('info','infocontour', 'SE', 'RE', 'score', 'empiricalhist', 'trace',
                         'itemscore', 'infotrace'))
            stop(type, " is not a valid plot type.", call.=FALSE)
        if (any(degrees > 90 | degrees < 0))
            stop('Improper angle specifed. Must be between 0 and 90.', call.=FALSE)
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
        ngroups <- x@Data$ngroups
        J <- x@Data$nitems
        nfact <- x@Model$nfact
        if(nfact > 2) stop("Can't plot high dimensional solutions.", call.=FALSE)
        if(nfact == 1) degrees <- 0
        theta <- seq(theta_lim[1],theta_lim[2],length.out=npts)
        ThetaFull <- Theta <- thetaComb(theta, nfact)
        prodlist <- attr(x@ParObjects$pars, 'prodlist')
        if(length(prodlist) > 0)
            ThetaFull <- prodterms(Theta,prodlist)
        infolist <- vector('list', ngroups)
        for(g in 1:ngroups)
            infolist[[g]] <- testinfo(extract.group(x, g), ThetaFull, degrees = degrees)
        if(type == 'RE') infolist <- lapply(infolist, function(x) x / infolist[[1]])
        info <- do.call(c, infolist)
        Theta <- ThetaFull
        groups <- gl(ngroups, nrow(ThetaFull), labels=x@Data$groupNames)
        adj <- x@Data$mins
        gscore <- c()
        for(g in 1:ngroups){
            itemtrace <- computeItemtrace(x@ParObjects$pars[[g]]@ParObjects$pars, ThetaFull, x@Model$itemloc,
                                          CUSTOM.IND=x@Internals$CUSTOM.IND)
            score <- c()
            for(i in 1:J)
                score <- c(score, 0:(x@Data$K[i]-1) + adj[i])
            score <- matrix(score, nrow(itemtrace), ncol(itemtrace), byrow = TRUE)
            gscore <- c(gscore, rowSums(score * itemtrace))
        }
        plt <- data.frame(info=info, score=gscore, Theta, group=groups)
        if(nfact == 2){
            colnames(plt) <- c("info", "score", "Theta1", "Theta2", "group")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'infocontour')
                return(contourplot(info ~ Theta1 * Theta2|group, data = plt,
                                   main = paste("Test Information Contour"), xlab = expression(theta[1]),
                                   ylab = expression(theta[2]),
                                   par.strip.text=par.strip.text, par.settings=par.settings, ...))
            if(type == 'info')
                return(wireframe(info ~ Theta1 + Theta2|group, data = plt, main = "Test Information",
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = auto.key, par.strip.text=par.strip.text, par.settings=par.settings,
                                 ...))
            if(type == 'RE')
                return(wireframe(info ~ Theta1 + Theta2|group, data = plt, main = "Relative Efficiency",
                                 zlab=expression(RE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = auto.key, par.strip.text=par.strip.text, par.settings=par.settings,
                                 ...))
            if(type == 'SE')
                return(wireframe(SE ~ Theta1 + Theta2|group, data = plt, main = "Test Standard Errors",
                                 zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = auto.key, par.strip.text=par.strip.text, par.settings=par.settings,
                                 ...))
            if(type == 'score')
                return(wireframe(score ~ Theta1 + Theta2|group, data = plt, main = "Expected Total Score",
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = auto.key, par.strip.text=par.strip.text, par.settings=par.settings,
                                 ...))
        } else {
            colnames(plt) <- c("info", "score", "Theta", "group")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'info')
                return(xyplot(info~Theta, plt, type='l', groups=plt$group, main = 'Test Information',
                              xlab = expression(theta), ylab=expression(I(theta)), auto.key = auto.key,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            if(type == 'RE')
                return(xyplot(info~Theta, plt, type='l', groups=plt$group, main = 'Relative Efficiency',
                              xlab = expression(theta), ylab=expression(RE(theta)), auto.key = auto.key,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            if(type == 'infocontour')
                cat('No \'contour\' plots for 1-dimensional models\n')
            if(type == 'SE')
                return(xyplot(SE~Theta, plt, type='l', groups=plt$group, main = 'Test Standard Errors',
                              xlab = expression(theta), ylab=expression(SE(theta)), auto.key = auto.key,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            if(type == 'score')
                return(xyplot(score~Theta, plt, type='l', groups=plt$group, main = 'Expected Total Score',
                              xlab = expression(theta), ylab=expression(Total(theta)), auto.key = auto.key,
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            if(type == 'empiricalhist'){
                Prior <- Theta <- pltfull <- vector('list', ngroups)
                for(g in 1L:ngroups){
                    Theta[[g]] <- as.matrix(seq(-(.8 * sqrt(x@Options$quadpts)), .8 * sqrt(x@Options$quadpts),
                                           length.out = x@Options$quadpts))
                    Prior[[g]] <- x@Internals$Prior[[g]] * nrow(x@Data$data)
                    cuts <- cut(Theta[[g]], floor(npts/2))
                    Prior[[g]] <- do.call(c, lapply(split(Prior[[g]], cuts), mean))
                    Theta[[g]] <- do.call(c, lapply(split(Theta[[g]], cuts), mean))
                    keep1 <- min(which(Prior[[g]] > 1e-10))
                    keep2 <- max(which(Prior[[g]] > 1e-10))
                    plt <- data.frame(Theta=Theta[[g]], Prior=Prior[[g]], group=x@Data$groupNames[g])
                    plt <- plt[keep1:keep2, , drop=FALSE]
                    pltfull[[g]] <- plt
                }
                plt <- do.call(rbind, pltfull)
                return(xyplot(Prior ~ Theta, plt, groups=plt$group, auto.key = auto.key,
                              xlab = expression(theta), ylab = 'Expected Frequency',
                              type = 'b', main = 'Empirical Histogram',
                              par.strip.text=par.strip.text, par.settings=par.settings, ...))
            }
            if(type == 'trace'){
                plt <- vector('list', ngroups)
                P <- vector('list', length(which.items))
                for(g in 1L:ngroups){
                    names(P) <- colnames(x@Data$data)[which.items]
                    count <- 1
                    for(i in which.items){
                        tmp <- probtrace(extract.item(x, i, group=g), ThetaFull)
                        if(ncol(tmp) == 2L) tmp <- tmp[,2, drop=FALSE]
                        tmp2 <- data.frame(P=as.numeric(tmp), cat=gl(ncol(tmp), k=nrow(ThetaFull),
                                                                     labels=paste0('cat', 1L:ncol(tmp))))
                        P[[count]] <- tmp2
                        count <- count + 1
                    }
                    nrs <- sapply(P, nrow)
                    Pstack <- do.call(rbind, P)
                    names <- c()
                    for(i in 1L:length(nrs))
                        if(!is.null(nrs[i]))
                            names <- c(names, rep(names(P)[i], nrs[i]))
                    plotobj <- data.frame(Pstack, item=names, Theta=ThetaFull, group=x@Data$groupNames[g])
                    plt[[g]] <- plotobj
                }
                plt <- do.call(rbind, plt)
                plt$item <- factor(plt$item, levels = colnames(x@Data$data)[which.items])
                if(facet_items){
                    return(xyplot(P ~ Theta|item, plt, groups = plt$cat:plt$group, ylim = c(-0.1,1.1),
                           xlab = expression(theta), ylab = expression(P(theta)),
                           auto.key = auto.key, type = 'l', main = 'Item trace lines',
                           par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(P ~ Theta|group, plt, groups = plt$cat:plt$item, ylim = c(-0.1,1.1),
                                  xlab = expression(theta), ylab = expression(P(theta)),
                                  auto.key = auto.key, type = 'l', main = 'Item trace lines',
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            }
            if(type == 'itemscore'){
                plt <- vector('list', ngroups)
                S <- vector('list', length(which.items))
                mins <- extract.mirt(x, 'mins')
                for(g in 1L:ngroups){
                    names(S) <- colnames(x@Data$data)[which.items]
                    count <- 1
                    for(i in which.items){
                        S[[count]] <- expected.item(extract.item(x, i, group=g), ThetaFull, min = mins[i])
                        count <- count + 1
                    }
                    Sstack <- do.call(c, S)
                    names <- rep(names(S), each = nrow(ThetaFull))
                    plotobj <- data.frame(S=Sstack, item=names, Theta=ThetaFull, group=x@Data$groupNames[g])
                    plt[[g]] <- plotobj
                }
                plt <- do.call(rbind, plt)
                plt$item <- factor(plt$item, levels = colnames(x@Data$data)[which.items])
                if(facet_items){
                    return(xyplot(S ~ Theta|item, plt, groups = plt$group,
                                  xlab = expression(theta), ylab = expression(S(theta)),
                                  auto.key = auto.key, type = 'l', main = 'Item score function',
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(S ~ Theta|group, plt, groups = plt$item,
                                  xlab = expression(theta), ylab = expression(S(theta)),
                                  auto.key = auto.key, type = 'l', main = 'Item score function',
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            }
            if(type == 'infotrace'){
                plt <- vector('list', ngroups)
                for(g in 1L:ngroups){
                    I <- matrix(NA, nrow(ThetaFull), J)
                    for(i in which.items)
                        I[,i] <- iteminfo(extract.item(x, i, group=g), ThetaFull)
                    I <- t(na.omit(t(I)))
                    items <- rep(colnames(x@Data$data)[which.items], each=nrow(Theta))
                    plotobj <- data.frame(I = as.numeric(I), Theta=ThetaFull, item=items)
                    plt[[g]] <- plotobj
                }
                plt <- do.call(rbind, plt)
                plt$item <- factor(plt$item, levels = colnames(x@Data$data)[which.items])
                plt$group <- rep(x@Data$groupNames, each = nrow(ThetaFull)*length(which.items))
                if(facet_items){
                    return(xyplot(I ~ Theta | item, plt, groups = plt$group,
                                  xlab = expression(theta), ylab = expression(I(theta)),
                                  auto.key = auto.key, type = 'l', main = 'Item information trace lines',
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                } else {
                    return(xyplot(I ~ Theta | group, plt, groups = plt$item,
                                  xlab = expression(theta), ylab = expression(I(theta)),
                                  auto.key = auto.key, type = 'l', main = 'Item information trace lines',
                                  par.strip.text=par.strip.text, par.settings=par.settings, ...))
                }
            }
        }
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, ...)
    {
        ret <- vector('list', length(object@Data$groupNames))
        names(ret) <- object@Data$groupNames
        for(g in 1L:length(ret)){
            cmod <- object@ParObjects$pars[[g]]
            cmod@Data <- object@Data
            cmod@Data$data <- object@Data$data[object@Data$group == object@Data$groupName[g], ]
            cmod@Data$Freq[[1L]] <- cmod@Data$Freq[[g]]
            cmod@Options$quadpts <- object@Options$quadpts
            cmod@Internals$bfactor <- object@Internals$bfactor
            ret[[g]] <- residuals(cmod, verbose = FALSE, ...)
        }
        ret
    }
)

# Methods
setMethod(
    f = "vcov",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object)
    {
        class(object) <- 'SingleGroupClass'
        vcov(object)
    }
)

setMethod(
    f = "logLik",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object){
        extract.mirt(object, 'logLik')
    }
)
