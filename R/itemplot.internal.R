setMethod(
	f = "itemplot.internal",
	signature = signature(object = 'ExploratoryClass'),
	definition = function(object, item, type, degrees, CE, CEalpha, CEdraws, ...)
	{  			
		x <- itemplot.main(object, item, type, degrees, CE, CEalpha, CEdraws, ...)		        
		return(invisible(x))
	}
)

#------------------------------------------------------------------------------
setMethod(
	f = "itemplot.internal",
	signature = signature(object = 'ConfirmatoryClass'),
	definition = function(object, item, type, degrees, CE, CEalpha, CEdraws, ...)
	{
	    x <- itemplot.main(object, item, type, degrees, CE, CEalpha, CEdraws, ...)    	
	    return(invisible(x))
	}
)

#------------------------------------------------------------------------------
setMethod(
    f = "itemplot.internal",
    signature = signature(object = 'list'),
    definition = function(object, item, type, degrees, CE, CEalpha, CEdraws, ...)
    {        
        newobject <- new('MultipleGroupClass', cmods=object, nfact=object[[1]]@nfact, 
                         groupNames=factor(names(object)))        
        x <- itemplot.internal(newobject, item, type, degrees, CE, CEalpha, CEdraws, ...)    	
        return(invisible(x))
    }
)

#------------------------------------------------------------------------------
setMethod(
    f = "itemplot.internal",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, item, type, degrees, CE, CEalpha, CEdraws, ...)
    {           
        Pinfo <- list()        
        gnames <- object@groupNames
        nfact <- object@nfact        
        K <- object@cmods[[1]]@pars[[item]]@ncat
        for(g in 1:length(gnames)){
            object@cmods[[g]]@information <- object@information
            Pinfo[[g]] <- itemplot.main(object@cmods[[g]], item=item, type='RETURN', 
                                        degrees=degrees, CE=FALSE, CEalpha=CEalpha, 
                                        CEdraws=CEdraws, ...)
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
            if(type == 'info')            
                return(lattice::xyplot(info ~ Theta, dat, group=group, type = 'l', 
                                       auto.key = TRUE, main = paste('Information for item', item), 
                                       ylab = expression(I(theta)), xlab = expression(theta), ...))            
            if(type == 'trace')
                return(lattice::xyplot(P  ~ Theta | cat, dat2, group=group, type = 'l', 
                            auto.key = TRUE, main = paste("Item", item, "Trace"), 
                            ylab = expression(P(theta)), xlab = expression(theta), ...))            
            if(type == 'RE')
                return(lattice::xyplot(info ~ Theta, dat, group=group, type = 'l', 
                                       auto.key = TRUE, 
                                       main = paste('Relative efficiency for item', item), 
                                       ylab = expression(RE(theta)), xlab = expression(theta), ...))
        }
        if(nfact == 2){
            Names <- colnames(dat)
            Names[c(length(Names) - 2,length(Names) - 1)] <- c('Theta1', 'Theta2')
            Names2 <- colnames(dat2)
            Names2[2:3] <- c('Theta2', 'Theta1')
            colnames(dat) <- Names
            colnames(dat2) <- Names2            
            if(type == 'info')            
                return(lattice::wireframe(info ~ Theta1 + Theta2, data = dat, group=group, 
                                          main=paste("Item", item, "Information"), 
                                          zlab=expression(I(theta)), xlab=expression(theta[1]), 
                                          ylab=expression(theta[2]), 
                                          scales = list(arrows = FALSE), 
                                          auto.key = TRUE, ...))            
            if(type == 'trace')
                return(lattice::wireframe(P ~ Theta1 + Theta2|cat, data = dat2, group = group, 
                                          main = paste("Item", item, "Trace"), 
                                          zlab=expression(P(theta)), 
                                          xlab=expression(theta[1]), 
                                          ylab=expression(theta[2]), 
                                          scales = list(arrows = FALSE), 
                                          auto.key = TRUE, ...))   
            if(type == 'RE')            
                return(lattice::wireframe(info ~ Theta1 + Theta2, data = dat, group=group, 
                                          main=paste("Relative efficiency for item", item), 
                                          zlab=expression(RE(theta)), xlab=expression(theta[1]), 
                                          ylab=expression(theta[2]), 
                                          scales = list(arrows = FALSE), 
                                          auto.key = TRUE, ...))
        }
    }
)


itemplot.main <- function(x, item, type, degrees, CE, CEalpha, CEdraws, ...){      
    nfact <- x@nfact
    if(nfact > 2) stop('Can not plot high dimensional models')
    if(nfact == 2 && is.null(degrees)) stop('Please specify a vector of angles that sum to 90')        
    theta <- seq(-4,4, length.out=40)
    Theta <- ThetaFull <- thetaComb(theta, nfact)   
    prodlist <- attr(x@pars, 'prodlist')
    if(length(prodlist) > 0)        
        ThetaFull <- prodterms(Theta,prodlist)
    P <- ProbTrace(x=x@pars[[item]], Theta=ThetaFull)         
    K <- x@pars[[item]]@ncat      
    info <- 0 
    if(nfact == 2){
        for(i in 1:length(degrees))
            info <- info + iteminfo(x=x@pars[[item]], Theta=ThetaFull, degrees=c(degrees[i], 
                                                                             90 - degrees[i]))
    } else {
        info <- iteminfo(x=x@pars[[item]], Theta=ThetaFull, degrees=0)
    }
    CEinfoupper <- CEinfolower <- info
    CEprobupper <- CEproblower <- P
    if(CE){                    
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
        delta <- mvtnorm::rmvnorm(CEdraws, mean=mu, sigma=smallinfo)
        tmp <- mvtnorm::dmvnorm(delta, mu, smallinfo)
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
    if(nfact > 1){        
        P2 <- P
        CEu2 <- CEl2 <- CEprobupper 
        for(i in 1:K){
            P2[,i] <- P[ ,ncol(P) + 1 - i]
            CEu2[,i] <- CEprobupper[ ,ncol(P) + 1 - i]
            CEl2[,i] <- CEproblower[ ,ncol(P) + 1 - i]
        }
        P <- P2
        CEprobupper <- CEu2
        CEproblower <- CEl2        
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
        plt$CEinfoupper <- CEinfoupper
        plt$CEinfolower <- CEinfolower        
        plt2$upper <- as.numeric(CEprobupper)
        plt2$lower <- as.numeric(CEproblower)
        if(type == 'trace'){
            if(CE){
                return(lattice::xyplot(P + upper + lower ~ Theta|time, plt2, type = 'l', 
                                col = c('black', 'red', 'red'), lty = c(1,2,2),
                                main = paste('Trace lines for item', item), ylim = c(-0.1,1.1),
                                ylab = expression(P(theta)), xlab = expression(theta), ... ))                
            }
            else
                return(lattice::xyplot(P ~ Theta, plt2, group = time, type = 'l', auto.key = TRUE,
                                main = paste('Trace lines for item', item), ylim = c(-0.1,1.1),
                                ylab = expression(P(theta)), xlab = expression(theta), ... ))               
        }
        if(type == 'info'){                        
            if(CE){                       
                return(lattice::xyplot(info + CEinfoupper + CEinfolower ~ Theta, plt, type = 'l', 
                                col = c('black', 'red', 'red'), lty = c(1,2,2),
                                main = paste('Information for item', item), 
                                ylab = expression(I(theta)), xlab = expression(theta), ... ))                                  
            } else 
                return(lattice::xyplot(info ~ Theta, plt, type = 'l', 
                                auto.key = TRUE, main = paste('Information for item', item), 
                                ylab = expression(I(theta)), xlab = expression(theta), ...))
        }
        if(type == 'score'){       
            return(lattice::xyplot(score ~ Theta, plt, type = 'l', 
                            auto.key = TRUE, main = paste('Expected score for item', item), 
                            ylab = expression(E(theta)), xlab = expression(theta), ...))                     
        }
        if(type == 'infocontour') stop('Cannot draw contours for 1 factor models')        
    } else {        
        plt <- data.frame(info = info, Theta1 = Theta[,1], Theta2 = Theta[,2])
        plt2 <- data.frame(P = P, Theta1 = Theta[,1], Theta2 = Theta[,2])
        colnames(plt2) <- c(paste("P", 1:ncol(P), sep=''), "Theta1", "Theta2")
        plt2 <- reshape(plt2, direction='long', varying = paste("P", 1:ncol(P), sep=''), v.names = 'P', 
                times = paste("P", 1:ncol(P), sep=''))
        colnames(plt) <- c("info", "Theta1", "Theta2")  
        plt$score <- score
        plt$CEinfoupper <- CEinfoupper
        plt$CEinfolower <- CEinfolower        
        plt2$upper <- as.numeric(CEprobupper)
        plt2$lower <- as.numeric(CEproblower)
        if(type == 'infocontour')												
            return(contourplot(info ~ Theta1 * Theta2, data = plt, 
                               main = paste("Item", item, "Information Contour"), xlab = expression(theta[1]), 
                               ylab = expression(theta[2]), ...))
        if(type == 'info')
            if(CE) 
                return(lattice::wireframe(info + CEinfolower + CEinfoupper ~ Theta1 + Theta2, data = plt, 
                                   main = paste("Item", item, "Information"), col = c('black', 'red', 'red'),
                                   zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                                   scales = list(arrows = FALSE), colorkey = TRUE, drape = TRUE, ...))
            else 
                return(lattice::wireframe(info ~ Theta1 + Theta2, data = plt, 
                             main = paste("Item", item, "Information"), 
                             zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                             scales = list(arrows = FALSE), colorkey = TRUE, drape = TRUE, ...))
        if(type == 'trace'){
            if(CE) 
                return(lattice::wireframe(P + upper + lower ~ Theta1 + Theta2 | time, data = plt2, 
                                          main = paste("Item", item, "Trace"), 
                                          zlab=expression(P(theta)), xlab=expression(theta[1]), 
                                          ylab=expression(theta[2]), col = c('black', 'red', 'red'),
                                          scales = list(arrows = FALSE), colorkey = TRUE, drape = TRUE, ...)) 
                
            else
                return(lattice::wireframe(P ~ Theta1 + Theta2, data = plt2, group = time, 
                             main = paste("Item", item, "Trace"), 
                             zlab=expression(P(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                             scales = list(arrows = FALSE), colorkey = TRUE, drape = TRUE, ...))            
        } 
        if(type == 'score'){            
            return(lattice::wireframe(score ~ Theta1 + Theta2, data = plt, main = paste("Item", item, "Expected Score"), 
                                      zlab=expression(E(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                                      zlim = c(min(floor(plt$score)), max(ceiling(plt$score))),scales = list(arrows = FALSE), 
                                      colorkey = TRUE, drape = TRUE, ...))
        }
    }    
}
