setMethod(
	f = "itemplot.internal",
	signature = signature(object = 'SingleGroupClass'),
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
        Data <- object[[1L]]@Data
        Data$groupNames <- factor(names(object))
        Model <- object[[1L]]@Model
        ParObjects <- object[[1L]]@ParObjects
        ParObjects$pars <- object
        newobject <- new('MultipleGroupClass', Data=Data, Model=Model, ParObjects=ParObjects)
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
        gnames <- object@Data$groupNames
        nfact <- object@Model$nfact
        K <- object@ParObjects$pars[[1L]]@ParObjects$pars[[item]]@ncat
        for(g in 1L:length(gnames)){
            object@ParObjects$pars[[g]]@vcov <- object@vcov
            Pinfo[[g]] <- itemplot.main(object@ParObjects$pars[[g]], item=item, type='RETURN',
                                        degrees=degrees, CE=FALSE, CEalpha=CEalpha,
                                        CEdraws=CEdraws, rot=rot, ...)
            Pinfo[[g]]$group <- rep(gnames[g], nrow(Pinfo[[g]]))
        }
        if(type == 'RE'){
            for(g in length(gnames):1L)
                Pinfo[[g]]$info <- Pinfo[[g]]$info / Pinfo[[1L]]$info
        }
        dat <- Pinfo[[1]]
        mins <- extract.mirt(object, 'mins')[item]
        score <- vector('list', g)
        for(g in 1L:length(gnames))
            score[[g]] <- colSums(t(Pinfo[[g]][,1L:K]) * (0L:(K-1L) + mins))
        for(g in 2L:length(gnames))
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
        ymin_score <- extract.mirt(object, 'mins')[item]
        ymax_score <- extract.mirt(object, 'K')[item] + ymin_score - 1
        ybump <- (ymax_score - ymin_score)/15
        if(nfact == 1){
            if(type == 'info'){
                if(is.null(main))
                    main <- paste('Information for item', item)
                return(xyplot(info ~ Theta, dat, groups=dat$group, type = 'l',
                                       auto.key = auto.key, main = main,
                                       ylab = expression(I(theta)), xlab = expression(theta), ...))
            } else if(type == 'trace'){
                if(is.null(main))
                    main <- paste("Item", item, "Trace")
                return(xyplot(P  ~ Theta | cat, dat2, groups=dat2$group, type = 'l',
                            auto.key = auto.key, main = main, ylim = c(-0.1,1.1),
                            ylab = expression(P(theta)), xlab = expression(theta), ...))
            } else if(type == 'score'){
                if(is.null(main))
                    main <- paste("Expected score for item", item)
                dat$score <- do.call(c, score)
                return(xyplot(score ~ Theta, dat, groups=dat$group, type = 'l',
                              ylim=c(ymin_score-ybump, ymax_score+ybump),
                              auto.key = auto.key, main = main,
                              ylab = expression(S(theta)), xlab = expression(theta), ...))
            } else if(type == 'RE'){
                if(is.null(main))
                    main <- paste('Relative efficiency for item', item)
                return(xyplot(info ~ Theta, dat, groups=dat$group, type = 'l',
                                       auto.key = auto.key, main = main,
                                       ylab = expression(RE(theta)), xlab = expression(theta), ...))
            } else {
                stop('Plot type not supported for unidimensional model', call.=FALSE)
            }
        }
        if(nfact == 2){
            Names <- colnames(dat)
            Names[c(length(Names) - 2,length(Names) - 1)] <- c('Theta1', 'Theta2')
            Names2 <- colnames(dat2)
            Names2[2:3] <- c('Theta1', 'Theta2')
            colnames(dat) <- Names
            colnames(dat2) <- Names2
            if(type == 'info'){
                if(is.null(main))
                    main <- paste("Item", item, "Information")
                return(wireframe(info ~ Theta1 + Theta2, data = dat, group=dat$group, main=main,
                                          zlab=expression(I(theta)), xlab=expression(theta[1]),
                                          ylab=expression(theta[2]), screen=rot,
                                          scales = list(arrows = FALSE),
                                          auto.key = auto.key, ...))
            } else if(type == 'trace'){
                if(is.null(main))
                    main <- paste("Item", item, "Trace")
                return(wireframe(P ~ Theta1 + Theta2|cat, data = dat2, group = dat2$group, main = main,
                                          zlab=expression(P(theta)),
                                          xlab=expression(theta[1]),
                                          ylab=expression(theta[2]), zlim = c(-0.1,1.1),
                                          scales = list(arrows = FALSE), screen=rot,
                                          auto.key = auto.key, ...))
            } else if(type == 'score'){
                if(is.null(main))
                    main <- paste("Expected score for item", item)
                dat$score <- do.call(c, score)
                return(wireframe(score ~ Theta1 + Theta2, data = dat, group=dat$group, main=main,
                                 ylim=c(ymin_score-ybump, ymax_score+ybump),
                                 zlab=expression(S(theta)), xlab=expression(theta[1]),
                                 ylab=expression(theta[2]),
                                 scales = list(arrows = FALSE), screen=rot,
                                 auto.key = auto.key, ...))
            } else if(type == 'RE'){
                if(is.null(main))
                    main <- paste("Relative efficiency for item", item)
                return(wireframe(info ~ Theta1 + Theta2, data = dat, group=dat$group, main=main,
                                          zlab=expression(RE(theta)), xlab=expression(theta[1]),
                                          ylab=expression(theta[2]),
                                          scales = list(arrows = FALSE), screen=rot,
                                          auto.key = auto.key, ...))
            } else {
                stop('Plot type not supported for 2 dimensional model', call.=FALSE)
            }
        }
    }
)


itemplot.main <- function(x, item, type, degrees, CE, CEalpha, CEdraws, drop.zeros, rot,
                          theta_lim, cuts = 30, colorkey = TRUE, auto.key = TRUE, main = NULL,
                          add.ylab2 = TRUE, drape = TRUE, ...){
    if(drop.zeros){
        if(x@Options$exploratory) stop('Cannot drop zeros in exploratory models', call.=FALSE)
        x@ParObjects$pars[[item]] <- extract.item(x, item, drop.zeros=TRUE)
    }
    nfact <- min(x@ParObjects$pars[[item]]@nfact, x@Model$nfact)
    if(nfact > 3) stop('Can not plot high dimensional models', call.=FALSE)
    theta <- seq(theta_lim[1L],theta_lim[2L], length.out=40)
    if(nfact == 3) theta <- seq(theta_lim[1L],theta_lim[2L], length.out=20)
    prodlist <- attr(x@ParObjects$pars, 'prodlist')
    if(length(prodlist) > 0){
        Theta <- thetaComb(theta, x@Model$nfact)
        ThetaFull <- prodterms(Theta,prodlist)
    } else Theta <- ThetaFull <- thetaComb(theta, nfact)
    if(length(degrees) == 1) degrees <- rep(degrees, ncol(ThetaFull))
    if(is(x, 'SingleGroupClass') && x@Options$exploratory){
        cfs <- coef(x, ..., verbose=FALSE, rawug=TRUE)
        x@ParObjects$pars[[item]]@par <- as.numeric(cfs[[item]][1L,])
    }
    P <- ProbTrace(x=x@ParObjects$pars[[item]], Theta=ThetaFull)
    K <- x@ParObjects$pars[[item]]@ncat
    info <- numeric(nrow(ThetaFull))
    if(K == 2L) auto.key <- FALSE
    if(type %in% c('info', 'SE', 'infoSE', 'infotrace', 'RE', 'infocontour', 'RETURN')){
        if(nfact == 1){
            info <- iteminfo(x=x@ParObjects$pars[[item]], Theta=ThetaFull, degrees=0)
        } else {
            info <- iteminfo(x=x@ParObjects$pars[[item]], Theta=ThetaFull, degrees=degrees)
        }
    }
    CEinfoupper <- CEinfolower <- info
    CEprobupper <- CEproblower <- P
    if(CE && nfact != 3){
        tmpitem <- x@ParObjects$pars[[item]]
        if(length(tmpitem@SEpar) == 0) stop('Must calculate the information matrix first.', call.=FALSE)
        splt <- strsplit(colnames(x@vcov), '\\.')
        parnums <- as.numeric(do.call(rbind, splt)[,2])
        tmp <- tmpitem@parnum[tmpitem@est]
        constrain <- x@Model$constrain
        if(length(constrain))
            for(i in 1:length(constrain))
                if(any(tmp %in% constrain[[i]]))
                    tmp[tmp %in% constrain[[i]]] <- constrain[[i]][1L]
        tmp <- parnums %in% tmp
        mu <- tmpitem@par[tmpitem@est]
        vcov <- extract.mirt(x, 'vcov')
        smallinfo <- vcov[tmp, tmp]
        smallinfo <-(smallinfo + t(smallinfo))/2 #make symetric
        delta <- mirt_rmvnorm(CEdraws, mean=mu, sigma=smallinfo)
        tmp <- mirt_dmvnorm(delta, mu, smallinfo)
        sorttmp <- sort(tmp)
        lower <- sorttmp[floor(length(tmp) * CEalpha/2)]
        upper <- sorttmp[ceiling(length(tmp) * (1-CEalpha/2))]
        delta <- delta[tmp < upper & tmp > lower, , drop=FALSE]
        tmpitem@par[tmpitem@est] <- delta[1, ]
        CEinfoupper <- CEinfolower <- iteminfo(tmpitem, ThetaFull, degrees=degrees)
        CEprobupper <- CEproblower <- ProbTrace(tmpitem, ThetaFull)
        CEscoreupper <- CEscorelower <- expected.item(tmpitem, ThetaFull, min = x@Data$mins[item])
        for(i in 2:nrow(delta)){
            tmpitem@par[tmpitem@est] <- delta[i, ]
            CEinfo <- iteminfo(tmpitem, ThetaFull, degrees=degrees)
            CEprob <- ProbTrace(tmpitem, ThetaFull)
            CEscore <- expected.item(tmpitem, ThetaFull, min = x@Data$mins[item])
            CEinfoupper <- apply(cbind(CEinfoupper, CEinfo), 1, max)
            CEinfolower <- apply(cbind(CEinfolower, CEinfo), 1, min)
            CEscoreupper <- apply(cbind(CEscoreupper, CEscore), 1, max)
            CEscorelower <- apply(cbind(CEscorelower, CEscore), 1, min)
            for(j in 1:ncol(CEprobupper)){
                CEprobupper[,j] <- apply(cbind(CEprobupper[,j], CEprob[,j]), 1, max)
                CEproblower[,j] <- apply(cbind(CEproblower[,j], CEprob[,j]), 1, min)
            }
        }
    }
    if(type == 'RETURN') return(data.frame(P=P, info=info, Theta=Theta))
    score <- expected.item(x@ParObjects$pars[[item]], Theta=ThetaFull, min=x@Data$mins[item])
    if(ncol(P) == 2){
        P <- P[ ,-1, drop = FALSE]
        CEprobupper <- CEprobupper[ ,-1, drop = FALSE]
        CEproblower <- CEproblower[ ,-1, drop = FALSE]
    }
    ymin_score <- extract.mirt(x, 'mins')[item]
    ymax_score <- extract.mirt(x, 'K')[item] + ymin_score - 1
    ybump <- (ymax_score - ymin_score)/15
    if(nfact == 1){
        plt <- data.frame(info = info, Theta = Theta)
        plt2 <- data.frame(P = P, Theta = Theta)
        colnames(plt2) <- c(paste("P", 1:ncol(P), sep=''), "Theta")
        plt2 <- reshape(plt2, direction='long', varying = paste("P", 1:ncol(P), sep=''), v.names = 'P',
                        times = paste("P", 1:ncol(P), sep=''))
        colnames(plt) <- c("info", "Theta")
        plt$score <- score
        plt$SE <- 1/sqrt(plt$info)
        if(CE){
            plt$CEinfoupper <- CEinfoupper
            plt$CEinfolower <- CEinfolower
            plt$CEscoreupper <- CEscoreupper
            plt$CEscorelower <- CEscorelower
            plt2$upper <- as.numeric(CEprobupper)
            plt2$lower <- as.numeric(CEproblower)
        }
        if(type == 'trace'){
            if(is.null(main))
                main <- paste('Trace lines for item', item)
            if(CE){
                return(xyplot(P ~ Theta|time, data=plt2,
                       upper=plt2$upper, lower=plt2$lower,
                       panel = function(x, y, lower, upper, subscripts, ...){
                           upper <- upper[subscripts]
                           lower <- lower[subscripts]
                           panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                         col="#E6E6E6", border = FALSE, ...)
                           panel.xyplot(x, y, type='l', lty=1,...)
                       },
                       main = main, ylim = c(-0.1,1.1), auto.key = auto.key,
                       ylab = expression(P(theta)), xlab = expression(theta), ...))
            } else {
                return(xyplot(P ~ Theta, plt2, groups = time, type = 'l', auto.key = auto.key,
                                main = main, ylim = c(-0.1,1.1),
                                ylab = expression(P(theta)), xlab = expression(theta), ... ))
            }
        } else if(type == 'info'){
            if(is.null(main))
                main <- paste('Information for item', item)
            if(CE){
                return(xyplot(info ~ Theta, data=plt, auto.key = auto.key,
                              upper=plt$CEinfoupper, lower=plt$CEinfolower,
                              panel = function(x, y, lower, upper, ...){
                                  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                col="#E6E6E6", border = FALSE, ...)
                                  panel.xyplot(x, y, type='l', lty=1,...)
                              },
                              main = main, ylim=c(min(plt$CEinfolower), max(plt$CEinfoupper)),
                              ylab = expression(I(theta)), xlab = expression(theta), ...))
            } else {
                return(xyplot(info ~ Theta, plt, type = 'l',
                                auto.key = auto.key, main = main,
                                ylab = expression(I(theta)), xlab = expression(theta), ...))
            }
        } else if(type == 'score'){
            if(is.null(main))
                main <- paste('Expected score for item', item)
            if(CE){
                return(xyplot(score ~ Theta, data=plt, auto.key = auto.key,
                              ylim=c(ymin_score-ybump, ymax_score+ybump),
                              upper=plt$CEscoreupper, lower=plt$CEscorelower,
                              panel = function(x, y, lower, upper, ...){
                                  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                col="#E6E6E6", border = FALSE, ...)
                                  panel.xyplot(x, y, type='l', lty=1,...)
                              },
                              main = main, ylim=c(min(plt$CEscorelower), max(plt$CEscoreupper)),
                              ylab = expression(S(theta)), xlab = expression(theta), ...))
            } else {
                return(xyplot(score ~ Theta, plt, type = 'l',
                              ylim=c(ymin_score-ybump, ymax_score+ybump),
                                auto.key = auto.key, main = main,
                                ylab = expression(E(theta)), xlab = expression(theta), ...))
            }
        } else if(type == 'SE'){
            if(is.null(main))
                main <- paste('Standard error plot for item', item)
            if(CE){
                plt$CESEupper <- 1/sqrt(CEinfolower)
                plt$CESElower <- 1/sqrt(CEinfoupper)
                return(xyplot(SE ~ Theta, data=plt, auto.key = auto.key,
                              upper=plt$CESEupper, lower=plt$CESElower,
                              panel = function(x, y, lower, upper, ...){
                                  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                                col="#E6E6E6", border = FALSE, ...)
                                  panel.xyplot(x, y, type='l', lty=1,...)
                              },
                              main = main, ylim=c(min(plt$CESElower), max(plt$CESEupper)),
                              ylab = expression(SE(theta)), xlab = expression(theta), ...))
            } else {
                return(xyplot(SE ~ Theta, plt, type = 'l',
                                       auto.key = auto.key, main = main,
                                       ylab = expression(SE(theta)), xlab = expression(theta), ...))
            }
        } else if(type == 'infoSE'){
            if(is.null(main))
                main <- paste('Item information and standard errors for item', item)
            obj1 <- xyplot(info~Theta, plt, type='l',
                           main = main, xlab = expression(theta), ylab=expression(I(theta)))
            obj2 <- xyplot(SE~Theta, plt, type='l', ylab=expression(SE(theta)))
            if(requireNamespace("latticeExtra", quietly = TRUE)){
                return(latticeExtra::doubleYScale(obj1, obj2, add.ylab2 = add.ylab2))
            } else {
                stop('latticeExtra package is not available. Please install.', call.=FALSE)
            }
        } else if(type == 'infotrace'){
            if(is.null(main))
                main <- paste('Trace lines and information for item', item)
            obj1 <- xyplot(P ~ Theta, plt2, type = 'l', lty = c(1:K), groups=time, main = main,
                           ylim = c(-0.1,1.1), ylab = expression(P(theta)), xlab = expression(theta), ... )
            obj2 <- xyplot(info~Theta, plt, type='l', xlab = expression(theta), ylab=expression(I(theta)),
                           ylim = c(-0.1,max(plt$info) + .5))
            if(requireNamespace("latticeExtra", quietly = TRUE)){
                return(latticeExtra::doubleYScale(obj1, obj2, add.ylab2 = add.ylab2))
            } else {
                stop('latticeExtra package is not available. Please install.', call.=FALSE)
            }
        } else {
            stop('Plot type not supported for unidimensional model', call.=FALSE)
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
            stop('Plot type not supported for 2 dimensional model', call.=FALSE)
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
            stop('Plot type not supported for 3 dimensional model', call.=FALSE)
        }

    }
}
