# Methods
setMethod(
    f = "print",
    signature = signature(x = 'MultipleGroupClass'),
    definition = function(x)
    {
        class(x) <- 'ExploratoryClass'
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
        ngroups <- length(object@pars)
        allPars <- vector('list', ngroups)
        names(allPars) <- object@Data$groupNames
        for(g in 1:ngroups){
            tmp <- object@pars[[g]]
            tmp@Data$data <- object@Data$data[1L, , drop=FALSE]
            allPars[[g]] <- coef(tmp, ...)
        }
        return(allPars)
    }
)

setMethod(
    f = "summary",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, digits = 3, verbose = TRUE, ...) {
        ngroups <- length(object@pars)
        groupind <- length(object@pars[[1]]@pars)
        nfact <- object@nfact
        ret <- list()
        coeflist <- coef(object)
        for(g in 1:ngroups){
            if(verbose) cat('\n----------\nGROUP:', as.character(object@Data$groupNames[g]), '\n')
            ret[[g]] <- summary(object@pars[[g]], digits=digits, verbose=verbose, ...)
            if(is(coeflist[[g]][[groupind]], 'matrix'))
                ret[[g]]$mean <- coeflist[[g]][[groupind]][1, 1:nfact]
            else ret[[g]]$mean <- coeflist[[g]][[groupind]][1:nfact]
            names(ret[[g]]$mean) <- colnames(ret[[g]]$fcor)
            if(verbose){
                cat('\nFactor means:\n')
                print(round(ret[[g]]$mean, digits))
            }
        }
        invisible(ret)
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, object2, ...)
    {
        class(object) <- 'ExploratoryClass'
        anova(object, object2, ...)
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'MultipleGroupClass', y = 'missing'),
    definition = function(x, y, type = 'info', npts = 50, theta_angle = 45,
                          which.items = 1:ncol(x@Data$data),
                          rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                          facet_items = TRUE, auto.key = TRUE,
                          theta_lim = c(-6,6), ...)
    {
        if (!type %in% c('info','infocontour', 'SE', 'RE', 'score', 'empiricalhist', 'trace', 'infotrace'))
            stop(type, " is not a valid plot type.")
        if (any(theta_angle > 90 | theta_angle < 0))
            stop('Improper angle specifed. Must be between 0 and 90.')
        if(length(theta_angle) > 1) stop('No info-angle plot is available')
        rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
        ngroups <- length(x@pars)
        J <- length(x@pars[[1]]@pars) - 1L
        nfact <- x@nfact
        if(nfact > 2) stop("Can't plot high dimensional solutions.")
        if(nfact == 1) theta_angle <- 0
        pars <- x@pars
        theta <- seq(theta_lim[1],theta_lim[2],length.out=npts)
        ThetaFull <- Theta <- thetaComb(theta, nfact)
        prodlist <- attr(x@pars, 'prodlist')
        if(length(prodlist) > 0)
            ThetaFull <- prodterms(Theta,prodlist)
        infolist <- vector('list', ngroups)
        for(g in 1:ngroups){
            info <- 0
            for(i in 1:J){
                tmp <- extract.item(x, i, g)
                info <- info + iteminfo(tmp, Theta=ThetaFull, degrees=theta_angle)
            }
            infolist[[g]] <- info
        }
        if(type == 'RE') infolist <- lapply(infolist, function(x) x / infolist[[1]])
        info <- do.call(rbind, infolist)
        Theta <- ThetaFull
        for(g in 2:ngroups) Theta <- rbind(Theta, ThetaFull)
        groups <- gl(ngroups, nrow(ThetaFull), labels=x@Data$groupNames)
        adj <- apply(x@Data$data, 2, min)
        gscore <- c()
        for(g in 1:ngroups){
            itemtrace <- computeItemtrace(x@pars[[g]]@pars, ThetaFull, x@itemloc, 
                                          CUSTOM.IND=x@CUSTOM.IND)
            score <- c()
            for(i in 1:J)
                score <- c(score, 0:(x@K[i]-1) + adj[i])
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
                                   ylab = expression(theta[2]), ...))
            if(type == 'info')
                return(wireframe(info ~ Theta1 + Theta2|group, data = plt, main = "Test Information",
                                 zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = TRUE, ...))
            if(type == 'RE')
                return(wireframe(info ~ Theta1 + Theta2|group, data = plt, main = "Relative Efficiency",
                                 zlab=expression(RE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = TRUE, ...))
            if(type == 'SE')
                return(wireframe(SE ~ Theta1 + Theta2|group, data = plt, main = "Test Standard Errors",
                                 zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = TRUE, ...))
            if(type == 'score')
                return(wireframe(score ~ Theta1 + Theta2|group, data = plt, main = "Expected Total Score",
                                 zlab=expression(Total(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen = rot, colorkey = TRUE, drape = TRUE,
                                 auto.key = TRUE, ...))
        } else {
            colnames(plt) <- c("info", "score", "Theta", "group")
            plt$SE <- 1 / sqrt(plt$info)
            if(type == 'info')
                return(xyplot(info~Theta, plt, type='l', group=group, main = 'Test Information',
                              xlab = expression(theta), ylab=expression(I(theta)), auto.key = TRUE, ...))
            if(type == 'RE')
                return(xyplot(info~Theta, plt, type='l', group=group, main = 'Relative Efficiency',
                              xlab = expression(theta), ylab=expression(RE(theta)), auto.key = TRUE, ...))
            if(type == 'infocontour')
                cat('No \'contour\' plots for 1-dimensional models\n')
            if(type == 'SE')
                return(xyplot(SE~Theta, plt, type='l', group=group, main = 'Test Standard Errors',
                              xlab = expression(theta), ylab=expression(SE(theta)), auto.key = TRUE, ...))
            if(type == 'score')
                return(xyplot(score~Theta, plt, type='l', group=group, main = 'Expected Total Score',
                              xlab = expression(theta), ylab=expression(Total(theta)), auto.key = TRUE, ...))
            if(type == 'empiricalhist'){
                if(!length(x@Prior)) stop('Empirical histogram was not estimated for this object')
                Prior <- Theta <- pltfull <- vector('list', ngroups)
                for(g in 1L:ngroups){
                    Theta[[g]] <- as.matrix(seq(-(.8 * sqrt(x@quadpts)), .8 * sqrt(x@quadpts),
                                           length.out = x@quadpts))
                    Prior[[g]] <- x@Prior[[g]] * nrow(x@Data$data)
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
                return(xyplot(Prior ~ Theta, plt, group=group, auto.key = TRUE,
                              xlab = expression(theta), ylab = 'Expected Frequency',
                              type = 'b', main = 'Empirical Histogram', ...))
            }
            if(type == 'trace'){
                plt <- vector('list', ngroups)
                P <- vector('list', length(which.items))
                for(g in 1L:ngroups){
                    names(P) <- colnames(x@Data$data)[which.items]
                    count <- 1
                    for(i in which.items){
                        tmp <- probtrace(extract.item(x, i, group=x@Data$groupNames[g]), ThetaFull)
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
                if(facet_items){
                    return(xyplot(P ~ Theta|item, plt, group = cat:group, ylim = c(-0.1,1.1),
                           xlab = expression(theta), ylab = expression(P(theta)),
                           auto.key = auto.key, type = 'l', main = 'Item trace lines', ...))
                } else {
                    return(xyplot(P ~ Theta|group, plt, group = cat:item, ylim = c(-0.1,1.1),
                                  xlab = expression(theta), ylab = expression(P(theta)),
                                  auto.key = auto.key, type = 'l', main = 'Item trace lines', ...))
                }
            }
            if(type == 'infotrace'){
                plt <- vector('list', ngroups)
                for(g in 1L:ngroups){
                    I <- matrix(NA, nrow(ThetaFull), J)
                    for(i in which.items)
                        I[,i] <- iteminfo(extract.item(x, i, group=x@Data$groupNames[g]), ThetaFull)
                    I <- t(na.omit(t(I)))
                    items <- gl(n=length(unique(which.items)), k=nrow(ThetaFull),
                                labels = paste('Item', which.items))
                    plotobj <- data.frame(I = as.numeric(I), Theta=ThetaFull, item=items, 
                                          group=x@Data$groupNames[g])
                    plt[[g]] <- plotobj
                }
                plt <- do.call(rbind, plt)
                if(facet_items){
                    return(xyplot(I ~ Theta | item, plt, group = group,
                                  xlab = expression(theta), ylab = expression(I(theta)),
                                  auto.key = auto.key, type = 'l', main = 'Item information trace lines', ...))
                } else {
                    return(xyplot(I ~ Theta | group, plt, group = item,
                                  xlab = expression(theta), ylab = expression(I(theta)),
                                  auto.key = auto.key, type = 'l', main = 'Item information trace lines', ...))
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
            cmod <- object@pars[[g]]
            cmod@Data <- object@Data
            cmod@Data$data <- object@Data$data[object@Data$group == object@Data$groupName[g], ]
            cmod@Data$Freq[[1L]] <- cmod@Data$Freq[[g]]
            cmod@quadpts <- object@quadpts
            cmod@bfactor <- object@bfactor
            ret[[g]] <- residuals(cmod, verbose = FALSE, ...)            
        }
        ret
    }
)
