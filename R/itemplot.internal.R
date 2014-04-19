setMethod(
	f = "itemplot.internal",
	signature = signature(object = 'ExploratoryClass'),
	definition = function(object, ...)
	{
		x <- itemplot.main(object, ...)
		return(invisible(x))
	}
)

#------------------------------------------------------------------------------
setMethod(
	f = "itemplot.internal",
	signature = signature(object = 'ConfirmatoryClass'),
	definition = function(object, ...)
	{
	    x <- itemplot.main(object, ...)
	    return(invisible(x))
	}
)

#------------------------------------------------------------------------------
setMethod(
    f = "itemplot.internal",
    signature = signature(object = 'list'),
    definition = function(object, ...)
    {
        newobject <- new('MultipleGroupClass', cmods=object, nfact=object[[1]]@nfact,
                         groupNames=factor(names(object)))
        x <- itemplot.internal(newobject, ...)
        return(invisible(x))
    }
)

#------------------------------------------------------------------------------
setMethod(
    f = "itemplot.internal",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, item, type, degrees, CE, CEalpha, CEdraws, rot,
        auto.key = TRUE, main = NULL, ...)
    {
        Pinfo <- list()
        gnames <- object@groupNames
        nfact <- object@nfact
        K <- object@cmods[[1]]@pars[[item]]@ncat
        for(g in 1:length(gnames)){
            object@cmods[[g]]@information <- object@information
            Pinfo[[g]] <- itemplot.main(object@cmods[[g]], item=item, type='RETURN',
                                        degrees=degrees, CE=FALSE, CEalpha=CEalpha,
                                        CEdraws=CEdraws, rot=rot, ...)
            Pinfo[[g]]$group <- rep(gnames[g], nrow(Pinfo[[g]]))
        }
        if(type == 'RE'){
            for(g in length(gnames):1)
                Pinfo[[g]]$info <- Pinfo[[g]]$info / Pinfo[[1]]$info
        }
        dat <- Pinfo[[1]]
        for(g in 2:length(gnames))
            dat <- rbind(dat, Pinfo[[g]])
        if(K == 2){
            K <- 1
            dat <- dat[ , -1]
        }
        Plist <- unclass(dat[, 1:K, drop = FALSE])
        P <- c()
        dat2 <- dat[, (K+1):ncol(dat)]
        for(i in 1:length(Plist))
            P <- c(P, Plist[[i]])
        if(length(Plist) > 1)
            for(i in 2:length(Plist))
                dat2 <- rbind(dat2, dat[, (K+1):ncol(dat)])
        dat2$P <- P
        dat2$cat <- rep(as.character(1:(length(Plist))), each = nrow(dat))
        if(all(dat2$cat == '0')) dat2$cat <- rep('1', length(dat2$cat))
        if(nfact == 1){
            if(type == 'info'){
                if(is.null(main))
                    main <- paste('Information for item', item)
                return(xyplot(info ~ Theta, dat, group=group, type = 'l',
                                       auto.key = auto.key, main = main,
                                       ylab = expression(I(theta)), xlab = expression(theta), ...))
            } else if(type == 'trace'){
                if(is.null(main))
                    main <- paste("Item", item, "Trace")
                return(xyplot(P  ~ Theta | cat, dat2, group=group, type = 'l',
                            auto.key = auto.key, main = main, ylim = c(-0.1,1.1),
                            ylab = expression(P(theta)), xlab = expression(theta), ...))
            } else if(type == 'RE'){
                if(is.null(main))
                    main <- paste('Relative efficiency for item', item)
                return(xyplot(info ~ Theta, dat, group=group, type = 'l',
                                       auto.key = auto.key, main = main,
                                       ylab = expression(RE(theta)), xlab = expression(theta), ...))
            } else {
                stop('Plot type not supported for unidimensional model')
            }
        }
        if(nfact == 2){
            Names <- colnames(dat)
            Names[c(length(Names) - 2,length(Names) - 1)] <- c('Theta1', 'Theta2')
            Names2 <- colnames(dat2)
            Names2[2:3] <- c('Theta2', 'Theta1')
            colnames(dat) <- Names
            colnames(dat2) <- Names2
            if(type == 'info'){
                if(is.null(main))
                    main <- paste("Item", item, "Information")
                return(wireframe(info ~ Theta1 + Theta2, data = dat, group=group, main=main,
                                          zlab=expression(I(theta)), xlab=expression(theta[1]),
                                          ylab=expression(theta[2]), screen=rot,
                                          scales = list(arrows = FALSE),
                                          auto.key = auto.key, ...))
            } else if(type == 'trace'){
                if(is.null(main))
                    main <- paste("Item", item, "Trace")
                return(wireframe(P ~ Theta1 + Theta2|cat, data = dat2, group = group, main = main,
                                          zlab=expression(P(theta)),
                                          xlab=expression(theta[1]),
                                          ylab=expression(theta[2]), zlim = c(-0.1,1.1),
                                          scales = list(arrows = FALSE), screen=rot,
                                          auto.key = auto.key, ...))
            } else if(type == 'RE'){
                if(is.null(main))
                    main <- paste("Relative efficiency for item", item)
                return(wireframe(info ~ Theta1 + Theta2, data = dat, group=group, main=main,
                                          zlab=expression(RE(theta)), xlab=expression(theta[1]),
                                          ylab=expression(theta[2]),
                                          scales = list(arrows = FALSE), screen=rot,
                                          auto.key = auto.key, ...))
            } else {
                stop('Plot type not supported for 2 dimensional model')
            }
        }
    }
)


itemplot.main <- function(x, item, type, degrees, CE, CEalpha, CEdraws, drop.zeros, rot,
                          theta_lim, cuts = 30, colorkey = TRUE, auto.key = TRUE, main = NULL,
                          add.ylab2 = TRUE, drape = TRUE, ...){
    if(drop.zeros) x@pars[[item]] <- extract.item(x, item, drop.zeros=TRUE)
    nfact <- min(x@pars[[item]]@nfact, x@nfact)
    if(nfact > 3) stop('Can not plot high dimensional models')
    if(nfact == 2 && is.null(degrees)) stop('Please specify a vector of angles that sum to 90')
    theta <- seq(theta_lim[1L],theta_lim[2L], length.out=40)
    if(nfact == 3) theta <- seq(theta_lim[1L],theta_lim[2L], length.out=20)
    prodlist <- attr(x@pars, 'prodlist')
    if(length(prodlist) > 0){
        Theta <- thetaComb(theta, x@nfact)
        ThetaFull <- prodterms(Theta,prodlist)
    } else Theta <- ThetaFull <- thetaComb(theta, nfact)
    if(is(x, 'ExploratoryClass')){
        cfs <- coef(x, ..., verbose=FALSE, rawug=TRUE)
        x@pars[[item]]@par <- as.numeric(cfs[[item]])
    }
    P <- ProbTrace(x=x@pars[[item]], Theta=ThetaFull)
    K <- x@pars[[item]]@ncat
    info <- 0
    if(is(x@pars[[item]], 'custom') && any(type %in% c('info', 'infocontour')))
        stop('Unable to compute information for custom items')
    if(!class(x@pars[[item]]) %in% c('custom')){
        if(nfact == 3){
            if(length(degrees) != 3 && any(type %in% 'info', 'SE')){
                warning('Information plots require the degrees input to be of length 3')
            } else {
                info <- iteminfo(x=x@pars[[item]], Theta=ThetaFull, degrees=degrees)
            }
        }
        if(nfact == 2){
            for(i in 1:length(degrees))
                info <- info + iteminfo(x=x@pars[[item]], Theta=ThetaFull, degrees=c(degrees[i],
                                                                                 90 - degrees[i]))
        } else {
            info <- iteminfo(x=x@pars[[item]], Theta=ThetaFull, degrees=0)
        }
    } else message('Information functions could not be computed')
    CEinfoupper <- CEinfolower <- info
    CEprobupper <- CEproblower <- P
    if(CE && nfact != 3){
        tmpitem <- x@pars[[item]]
        if(length(tmpitem@SEpar) == 0) stop('Must calculate the information matrix first.')
        splt <- strsplit(colnames(x@information), '\\.')
        parnums <- as.numeric(do.call(rbind, splt)[,2])
        tmp <- x@pars[[item]]@parnum[x@pars[[item]]@est]
        if(length(x@constrain) > 0)
            for(i in 1:length(x@constrain))
                if(any(tmp %in% x@constrain[[i]]))
                    tmp[tmp %in% x@constrain[[i]]] <- x@constrain[[i]][1]
        tmp <- parnums %in% tmp
        mu <- tmpitem@par[x@pars[[item]]@est]
        smallinfo <- solve(x@information[tmp, tmp])
        #make symetric
        smallinfo <-(smallinfo + t(smallinfo))/2
        delta <- mirt_rmvnorm(CEdraws, mean=mu, sigma=smallinfo)
        tmp <- mirt_dmvnorm(delta, mu, smallinfo)
        sorttmp <- sort(tmp)
        lower <- sorttmp[floor(length(tmp) * CEalpha/2)]
        upper <- sorttmp[ceiling(length(tmp) * (1-CEalpha/2))]
        delta <- delta[tmp < upper & tmp > lower, , drop=FALSE]
        tmpitem@par[tmpitem@est] <- delta[1, ]
        degrees <- if(nfact == 2) c(degrees[i], 90 - degrees[i]) else 0
        CEinfoupper <- CEinfolower <- iteminfo(tmpitem, ThetaFull, degrees=degrees)
        CEprobupper <- CEproblower <- ProbTrace(tmpitem, ThetaFull)
        for(i in 2:nrow(delta)){
            tmpitem@par[tmpitem@est] <- delta[i, ]
            CEinfo <- iteminfo(tmpitem, ThetaFull, degrees=degrees)
            CEprob <- ProbTrace(tmpitem, ThetaFull)
            CEinfoupper <- apply(cbind(CEinfoupper, CEinfo), 1, max)
            CEinfolower <- apply(cbind(CEinfolower, CEinfo), 1, min)
            for(j in 1:ncol(CEprobupper)){
                CEprobupper[,j] <- apply(cbind(CEprobupper[,j], CEprob[,j]), 1, max)
                CEproblower[,j] <- apply(cbind(CEproblower[,j], CEprob[,j]), 1, min)
            }
        }
    }
    if(type == 'RETURN') return(data.frame(P=P, info=info, Theta=Theta))
    score <- matrix(0:(ncol(P) - 1), nrow(Theta), ncol(P), byrow = TRUE)
    score <- rowSums(score * P)
    if(class(x@pars[[item]]) %in% c('nominal', 'graded', 'rating'))
        score <- score + 1
    if(ncol(P) == 2){
        P <- P[ ,-1, drop = FALSE]
        CEprobupper <- CEprobupper[ ,-1, drop = FALSE]
        CEproblower <- CEproblower[ ,-1, drop = FALSE]
    }
    if(nfact == 1){
        plt <- data.frame(info = info, Theta = Theta)
        plt2 <- data.frame(P = P, Theta = Theta)
        colnames(plt2) <- c(paste("P", 1:ncol(P), sep=''), "Theta")
        plt2 <- reshape(plt2, direction='long', varying = paste("P", 1:ncol(P), sep=''), v.names = 'P',
                        times = paste("P", 1:ncol(P), sep=''))
        colnames(plt) <- c("info", "Theta")
        plt$score <- score
        plt$SE <- 1/sqrt(plt$info)
        plt$CEinfoupper <- CEinfoupper
        plt$CEinfolower <- CEinfolower
        plt2$upper <- as.numeric(CEprobupper)
        plt2$lower <- as.numeric(CEproblower)
        if(type == 'trace'){
            if(is.null(main))
                main <- paste('Trace lines for item', item)
            if(CE){
                return(xyplot(P + upper + lower ~ Theta|time, plt2, type = 'l',
                                col = c('black', 'red', 'red'), lty = c(1,2,2),
                                main = main, ylim = c(-0.1,1.1),
                                ylab = expression(P(theta)), xlab = expression(theta), ... ))
            }
            else
                return(xyplot(P ~ Theta, plt2, group = time, type = 'l', auto.key = auto.key,
                                main = main, ylim = c(-0.1,1.1),
                                ylab = expression(P(theta)), xlab = expression(theta), ... ))
        } else if(type == 'info'){
            if(is.null(main))
                main <- paste('Information for item', item)
            if(CE){
                return(xyplot(info + CEinfoupper + CEinfolower ~ Theta, plt, type = 'l',
                                col = c('black', 'red', 'red'), lty = c(1,2,2),
                                main = main, ylab = expression(I(theta)), xlab = expression(theta), ... ))
            } else
                return(xyplot(info ~ Theta, plt, type = 'l',
                                auto.key = auto.key, main = main,
                                ylab = expression(I(theta)), xlab = expression(theta), ...))
        } else if(type == 'score'){
            if(is.null(main))
                main <- paste('Expected score for item', item)
            return(xyplot(score ~ Theta, plt, type = 'l',
                            auto.key = auto.key, main = main,
                            ylab = expression(E(theta)), xlab = expression(theta), ...))
        } else if(type == 'SE'){
            if(is.null(main))
                main <- paste('Standard error plot for item', item)
            return(xyplot(SE ~ Theta, plt, type = 'l',
                                   auto.key = auto.key, main = main,
                                   ylab = expression(SE(theta)), xlab = expression(theta), ...))
        } else if(type == 'infoSE'){
            if(is.null(main))
                main <- paste('Item information and standard errors for item', item)
            obj1 <- xyplot(info~Theta, plt, type='l',
                           main = main, xlab = expression(theta), ylab=expression(I(theta)))
            obj2 <- xyplot(SE~Theta, plt, type='l', ylab=expression(SE(theta)))
            if(!require(latticeExtra)) require(latticeExtra)
            return(doubleYScale(obj1, obj2, add.ylab2 = add.ylab2))
        } else if(type == 'infotrace'){
            if(is.null(main))
                main <- paste('Trace lines and information for item', item)
            obj1 <- xyplot(P ~ Theta, plt2, type = 'l', lty = c(1:K), group=time, main = main,
                           ylim = c(-0.1,1.1), ylab = expression(P(theta)), xlab = expression(theta), ... )
            obj2 <- xyplot(info~Theta, plt, type='l', xlab = expression(theta), ylab=expression(I(theta)),
                           ylim = c(-0.1,max(plt$info) + .5))
            if(!require(latticeExtra)) require(latticeExtra)
            return(doubleYScale(obj1, obj2, add.ylab2 = add.ylab2))
        } else {
            stop('Plot type not supported for unidimensional model')
        }
    } else if(nfact == 2){
        plt <- data.frame(info = info, SE = 1/sqrt(info), Theta1 = Theta[,1], Theta2 = Theta[,2])
        plt2 <- data.frame(P = P, Theta1 = Theta[,1], Theta2 = Theta[,2])
        colnames(plt2) <- c(paste("P", 1:ncol(P), sep=''), "Theta1", "Theta2")
        plt2 <- reshape(plt2, direction='long', varying = paste("P", 1:ncol(P), sep=''), v.names = 'P',
                times = paste("P", 1:ncol(P), sep=''))
        plt$score <- score
        plt$CEinfoupper <- CEinfoupper
        plt$CEinfolower <- CEinfolower
        plt2$upper <- as.numeric(CEprobupper)
        plt2$lower <- as.numeric(CEproblower)
        if(type == 'infocontour'){
            if(is.null(main))
                main <- paste("Item", item, "Information Contour")
            return(contourplot(info ~ Theta1 * Theta2, data = plt,
                               main = main, xlab = expression(theta[1]),
                               ylab = expression(theta[2]), ...))
        } else if(type == 'tracecontour'){
            if(is.null(main))
                main <- paste("Item", item, "Probabiliy Contour")
            return(contourplot(P ~ Theta1 * Theta2, data = plt2,
                               main = main, xlab = expression(theta[1]),
                               ylab = expression(theta[2]), cuts=cuts, ...))
        } else if(type == 'info'){
            if(is.null(main))
                main <- paste("Item", item, "Information")
            if(CE)
                return(wireframe(info + CEinfolower + CEinfoupper ~ Theta1 + Theta2, data = plt,
                                   main = main, col = c('black', 'red', 'red'),
                                   zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                   scales=list(arrows = FALSE), colorkey=colorkey, drape=drape, screen=rot, ...))
            else
                return(wireframe(info ~ Theta1 + Theta2, data = plt, main = main,
                             zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                             scales = list(arrows = FALSE), colorkey=colorkey, drape = drape, screen=rot, ...))
        } else if(type == 'trace'){
            if(is.null(main))
                main <- paste("Item", item, "Trace")
            if(CE)
                return(wireframe(P + upper + lower ~ Theta1 + Theta2 | time, data = plt2,
                                          main = main, zlim = c(-0.1,1.1),
                                          zlab=expression(P(theta)), xlab=expression(theta[1]),
                                          ylab=expression(theta[2]), col = c('black', 'red', 'red'),
                                          scales=list(arrows = FALSE), colorkey=colorkey, drape=drape, screen=rot, ...))

            else
                return(wireframe(P ~ Theta1 + Theta2, data = plt2, group = time,
                             main = main, zlim = c(-0.1,1.1),
                             zlab=expression(P(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                             scales = list(arrows = FALSE), colorkey = colorkey, drape = drape, screen=rot, ...))
        } else if(type == 'score'){
            if(is.null(main))
                main <- paste("Item", item, "Expected Scores")
            return(wireframe(score ~ Theta1 + Theta2, data = plt, main = main,
                                      zlab=expression(E(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                      zlim = c(min(floor(plt$score)), max(ceiling(plt$score))),scales = list(arrows = FALSE),
                                      colorkey = colorkey, drape = drape, screen=rot, ...))
        } else if(type == 'SE'){
            if(is.null(main))
                main <- paste("Item", item, "Standard Errors")
            return(wireframe(SE ~ Theta1 + Theta2, data = plt, main = main,
                                      zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                      scales = list(arrows = FALSE),
                                      colorkey = colorkey, drape = drape, screen=rot, ...))
        } else {
            stop('Plot type not supported for 2 dimensional model')
        }
    } else {
        plt <- data.frame(info = info, SE = 1/sqrt(info), Theta1 = Theta[,1], Theta2 = Theta[,2],
                          Theta3 = Theta[,3])
        plt2 <- data.frame(P = P, Theta1 = Theta[,1], Theta2 = Theta[,2], Theta3 = Theta[,3])
        colnames(plt2) <- c(paste("P", 1:ncol(P), sep=''), "Theta1", "Theta2", "Theta3")
        plt2 <- reshape(plt2, direction='long', varying = paste("P", 1:ncol(P), sep=''), v.names = 'P',
                        times = paste("P", 1:ncol(P), sep=''))
        plt$score <- score
        if(type == 'trace'){
            if(is.null(main))
                main <- paste("Item", item, "Trace")
            return(wireframe(P ~ Theta1 + Theta2|Theta3, data = plt2, group = time,
                                      main = main, zlim = c(-0.1,1.1),
                                      zlab=expression(P(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                      scales = list(arrows = FALSE), colorkey = colorkey, drape = drape, screen=rot, ...))
        } else if(type == 'score'){
            if(is.null(main))
                main <- paste("Item", item, "Expected Scores")
            return(wireframe(score ~ Theta1 + Theta2|Theta3, data = plt, main = main,
                                      zlab=expression(E(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                      zlim = c(min(floor(plt$score)), max(ceiling(plt$score))),scales = list(arrows = FALSE),
                                      colorkey = colorkey, drape = drape, screen=rot, ...))
        } else if(type == 'info'){
            if(is.null(main))
                main <- paste("Item", item, "Information")
            return(wireframe(info ~ Theta1 + Theta2|Theta3, data = plt, main = main,
                                      zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                      scales = list(arrows = FALSE), colorkey = colorkey, drape = drape, screen=rot, ...))
        } else if(type == 'SE'){
            if(is.null(main))
                main <- paste("Item", item, "Standard Errors")
            return(wireframe(SE ~ Theta1 + Theta2|Theta3, data = plt, main = main,
                                      zlab=expression(SE(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]),
                                      scales = list(arrows = FALSE),
                                      colorkey = colorkey, drape = drape, screen=rot, ...))
        } else {
            stop('Plot type not supported for 3 dimensional model')
        }

    }
}
