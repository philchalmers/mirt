# theta combinations
thetaComb <- function(theta, nfact)
{
	if (nfact == 1) Theta <- matrix(theta)
	else if (nfact == 2) Theta <- expand.grid(theta,theta)   
	else if (nfact == 3) Theta <- expand.grid(theta,theta,theta)  
	else if (nfact == 4) Theta <- expand.grid(theta,theta,theta,theta)
	else if (nfact == 5) Theta <- expand.grid(theta,theta,theta,theta,theta)        	
	else if (nfact == 6) Theta <- expand.grid(theta,theta,theta,theta,theta,theta)
	if(nfact > 6) stop('Are you crazy?!?!? That is way too many factors for this quadrature method.
                       Try using confmirt() instead for better numerical accuracy')
	Theta <- as.matrix(Theta)	
	return(Theta)     
}

# Product terms in confmirt
prodterms <- function(theta0, prodlist)
{
    products <- matrix(1, ncol = length(prodlist), nrow = nrow(theta0))
    for(i in 1:length(prodlist)){
        tmp <- prodlist[[i]]
        for(j in 1:length(tmp)) 
            products[ ,i] <- products[ ,i] * theta0[ ,tmp[j]]	
    }	
    ret <- cbind(theta0,products)
    ret
}

# MH sampler for theta values
draw.thetas <- function(theta0, pars, fulldata, itemloc, cand.t.var, prior.t.var, 
                        prior.mu, prodlist, mixedlist = NULL, debug) 
{         
    if(debug == 'draw.thetas') browser()
    tol <- 1e-8
    N <- nrow(fulldata)
    J <- length(pars) - 1
    nfact <- ncol(theta0)					
    unif <- runif(N)
    if(nfact > 1){		
        theta1 <- theta0 + mvtnorm::rmvnorm(N,prior.mu, 
                                            diag(rep(cand.t.var,ncol(theta0)))) 
    } else
        theta1 <- theta0 + rnorm(N,prior.mu,sqrt(cand.t.var))							
    log_den0 <- mvtnorm::dmvnorm(theta0,prior.mu,prior.t.var,log=TRUE)
    log_den1 <- mvtnorm::dmvnorm(theta1,prior.mu,prior.t.var,log=TRUE)		
    if(length(prodlist) > 0){
        theta0 <- prodterms(theta0,prodlist)
        theta1 <- prodterms(theta1,prodlist)	
    }	    
    itemtrace0 <- itemtrace1 <- matrix(0, ncol=ncol(fulldata), nrow=nrow(theta0))    
    fixed.design0 <- fixed.design1 <- NULL    
    if(!is.null(mixedlist)){
        colnames(theta0) <- colnames(theta1) <- mixedlist$factorNames
        fixed.design0 <- designMats(covdata=mixedlist$covdata, fixed=mixedlist$fixed, 
                                    Thetas=theta0, nitems=J, 
                                    itemdesign=mixedlist$itemdesign, 
                                    fixed.identical=mixedlist$fixed.identical)
        fixed.design1 <- designMats(covdata=mixedlist$covdata, fixed=mixedlist$fixed, 
                                    Thetas=theta1, nitems=J, 
                                    itemdesign=mixedlist$itemdesign, 
                                    fixed.identical=mixedlist$fixed.identical)
    }    
    for (i in 1:J){
        itemtrace0[ ,itemloc[i]:(itemloc[i+1] - 1)] <- 
            ProbTrace(x=pars[[i]], Theta=theta0, fixed.design=fixed.design0[[i]])
        itemtrace1[ ,itemloc[i]:(itemloc[i+1] - 1)] <- 
            ProbTrace(x=pars[[i]], Theta=theta1, fixed.design=fixed.design1[[i]])        
    }    
    tmp0 <- itemtrace0*fulldata
    tmp1 <- itemtrace1*fulldata
    tmp0[tmp0 < tol] <- tmp1[tmp1 < tol] <- 1    
    total_0 <- rowSums(log(tmp0)) + log_den0
    total_1 <- rowSums(log(tmp1)) + log_den1    
    diff <- total_1 - total_0
    diff[is.nan(diff)] <- -50
    accept <- diff > 0
    accept[unif < exp(diff)] <- TRUE    
    theta1[!accept, ] <- theta0[!accept, ]
    total_1[!accept] <- total_0[!accept]
    log.lik <- sum(total_1)	
    if(!is.null(prodlist)) 
        theta1 <- theta1[ ,1:(pars[[1]]@nfact - pars[[1]]@nfixedeffects - 
                                  length(prodlist)), drop=FALSE]
    attr(theta1, "Proportion Accepted") <- sum(accept)/N 				
    attr(theta1, "log.lik") <- log.lik	
    return(theta1) 
}

# Rotation function
Rotate <- function(F, rotate, Target = NULL, ...)
{	
    if(ncol(F) == 1) rotF <- list()    
    if(rotate == 'none') rotF <- list(loadings=F, Phi=diag(ncol(F)), orthogonal=TRUE)            
	if(rotate == 'promax'){
        mypromax <- function (x, m = 4) {
                #borrowed and modified from stats::promax on Febuary 13, 2013
                if (ncol(x) < 2) 
                    return(x)
                dn <- dimnames(x)
                xx <- varimax(x)
                x <- xx$loadings
                Q <- x * abs(x)^(m - 1)
                U <- lm.fit(x, Q)$coefficients
                d <- diag(solve(t(U) %*% U))
                U <- U %*% diag(sqrt(d))
                dimnames(U) <- NULL
                z <- x %*% U
                U <- xx$rotmat %*% U
                ui <- solve(U)
                Phi <- ui %*% t(ui)
                dimnames(z) <- dn
                class(z) <- "loadings"
                list(loadings = z, rotmat = U, Phi = Phi, orthogonal = FALSE)
            }               
        rotF <- mypromax(F, ...)                
	}    
    if(rotate == 'oblimin') rotF <- GPArotation::oblimin(F, ...)     
	if(rotate == 'quartimin') rotF <- GPArotation::quartimin(F, ...)
	if(rotate == 'targetT') rotF <- GPArotation::targetT(F, Target = Target, ...)
	if(rotate == 'targetQ') rotF <- GPArotation::targetQ(F, Target = Target, ...)
	if(rotate == 'pstT') rotF <- GPArotation::pstT(F, Target = Target, ...)
	if(rotate == 'pstQ') rotF <- GPArotation::pstQ(F, Target = Target, ...)
	if(rotate == 'oblimax') rotF <- GPArotation::oblimax(F, ...)
	if(rotate == 'entropy') rotF <- GPArotation::entropy(F, ...)
	if(rotate == 'quartimax') rotF <- GPArotation::quartimax(F, ...)
	if(rotate == 'varimax') rotF <- GPArotation::Varimax(F, ...)
	if(rotate == 'simplimax') rotF <- GPArotation::simplimax(F, ...)
	if(rotate == 'bentlerT') rotF <- GPArotation::bentlerT(F, ...)
	if(rotate == 'bentlerQ') rotF <- GPArotation::bentlerQ(F, ...)
	if(rotate == 'tandemI') rotF <- GPArotation::tandemI(F, ...)
	if(rotate == 'tandemII') rotF <- GPArotation::tandemII(F, ...)
	if(rotate == 'geominT') rotF <- GPArotation::geominT(F, ...)
	if(rotate == 'geominQ') rotF <- GPArotation::geominQ(F, ...)
	if(rotate == 'cfT') rotF <- GPArotation::cfT(F, ...)
	if(rotate == 'cfQ') rotF <- GPArotation::cfQ(F, ...)
	if(rotate == 'infomaxT') rotF <- GPArotation::infomaxT(F, ...)
	if(rotate == 'infomaxQ') rotF <- GPArotation::infomaxQ(F, ...)
	if(rotate == 'mccammon') rotF <- GPArotation::mccammon(F, ...)
	if(rotate == 'bifactorT') rotF <- GPArotation::bifactorT(F, ...)
	if(rotate == 'bifactorQ') rotF <- GPArotation::bifactorQ(F, ...)
	return(unclass(rotF))
}  

# Gamma correlation, mainly for obtaining a sign
gamma.cor <- function(x)
{ 
	concordant <- function(x){ 	  
			mat.lr <- function(r, c){ 
				lr <- x[(r.x > r) & (c.x > c)] 
				sum(lr) 
			} 	  
		r.x <- row(x) 
		c.x <- col(x) 	  
		sum(x * mapply(mat.lr, r = r.x, c = c.x)) 
	} 	
	discordant <- function(x){ 	  
		mat.ll <- function(r, c){ 
			ll <- x[(r.x > r) & (c.x < c)] 
			sum(ll) 
		} 	  
		r.x <- row(x) 
		c.x <- col(x) 	  
		sum(x * mapply(mat.ll, r = r.x, c = c.x)) 
	} 
	c <- concordant(x) 
	d <- discordant(x) 
	gamma <- (c - d) / (c + d) 
	gamma 
} 

# Beta prior for grad and hess
betaprior <- function(g,a,b)
{	
    if(g < 1e-8) return(list(grad=0, hess=0))
	grad <- ((a-1) * g^(a-1) * (1-g)^(b-1) - (b-1)*g^(a-1)*(1-g)^(b-1))/ 
		(g^(a-1) * (1-g)^(b-1))
	hess <- -((g^(a-1)*(a-1)^2*(1-g)^(b-1)/g^2 - g^(a-1)*(a-1)*(1-g)^(b-1)/g^2 
		- 2*g^(a-1)*(a-1)*(1-g)^(b-1)*(b-1)/(g*(1-g)) + g^(a-1)*(1-g)^(b-1)*(b-1)^2/(1-g)^2 
		- g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g)^2)/(g^(a-1)*(1-g)^(b-1))	
		- ((g^(a-1)*(a-1)*(1-g)^(b-1)/g-g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g))*(a-1)/(g^(a-1)*(1-g)^(b-1)*g))
		+ ((g^(a-1)*(a-1)*(1-g)^(b-1)/g-g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g))*(b-1)/(g^(a-1)*(1-g)^(b-1)*(1-g))))
	return(list(grad=grad, hess=hess))	
}

# Approximation to polychoric matrix for initial values
cormod <- function(fulldata, K, guess, smooth = TRUE, use = 'pairwise.complete.obs') 
{  
	fulldata <- as.matrix(fulldata) 
	nitems <- ncol(fulldata)                 
	cormat <- cor(fulldata, use=use)      	
	cormat <- abs(cormat)^(1/1.15) * sign(cormat)  
	if(smooth){  
		eig <- eigen(cormat)
		negvalues <- eig$values < 0 
		if (any(negvalues)) {
			negeig <- sum(abs(eig$value[eig$value < 0])) 
			eig$value[eig$value < 0] <- 0
			L <- nitems/sum(eig$value)*eig$value[!negvalues]
			V <- eig$vector[ ,!negvalues] 
			cormat <- V %*% diag(L) %*% t(V)    
		}      
	}	
	cormat
}  

# Ramsey rate acceleration adjustment for EM
rateChange <- function(pars, listpars, lastpars1, lastpars2)
{   
	p <- unlist(listpars)	
	lp1 <- unlist(lastpars1)
	lp2 <- unlist(lastpars2)
	rate <- rep(0, length(p))
	d1 <- lp1 - p
	d2 <- lp2 - p
	rate <- ifelse(abs(d1) > 0.001 & (d1*d2 > 0.0) & (d1/d2 < 1.0),
		(1 - (1 - rate) * (d1/d2)),
		0)	    		
	rate[p > 4] <- 0
	rate[p < -4] <- 0    
	p <- lp1*rate*(-2) + (1 - rate*(-2))*p
    ind <- 1
    for(i in 1:length(pars)){
        pars[[i]]@par <- p[ind:(ind + length(pars[[i]]@par) - 1)]
        ind <- ind + length(pars[[i]]@par)
    }	
	pars
}

# Rotate lambda coefficients
rotateLambdas <- function(so){    
    F <- so$rotF
    h2 <- so$h2
    h <- matrix(rep(sqrt(1 - h2), ncol(F)), ncol = ncol(F))
    a <- F / h
    a    
}

d2r <-function(d) pi*d/180

test_info <- function(pars, Theta, Alist, K){
    infolist <- list()    
    for(cut in 1:length(Alist)){
        A <- Alist[[cut]]
        info <- rep(0,nrow(Theta))
        for(j in 1:length(K)){
            info <- info + ItemInfo(pars[[j]], A[j,], Theta)
        }
        infolist[[cut]] <- info
    }    
    tmp <- 0
    for(i in 1:length(infolist)){
        tmp <- tmp + infolist[[i]]
    }
    info <- tmp/length(infolist)
    info
}

Lambdas <- function(pars){    
    lambdas <- list()
    J <- ifelse(is(pars[[length(pars)]], 'GroupPars'), length(pars)-1, length(pars))
    for(i in 1:J)    
        lambdas[[i]] <- ExtractLambdas(pars[[i]])    
    lambdas <- do.call(rbind,lambdas)
    lambdas
}

#change long pars for groups into mean in sigma
ExtractGroupPars <- function(x){
    nfact <- x@nfact
    gmeans <- x@par[1:nfact]    
    tmp <- x@par[-(1:nfact)]    
    gcov <- matrix(0, nfact, nfact)
    gcov[lower.tri(gcov, diag=TRUE)] <- tmp
    if(nfact != 1)
        gcov <- gcov + t(gcov) - diag(diag(gcov))
    return(list(gmeans=gmeans, gcov=gcov))    
}

reloadConstr <- function(par, constr, obj){
    par2 <- rep(NA, length(constr[[1]]))         
    notconstr <- rep(TRUE, length(par2))
    for(i in 1:length(constr)){
        par2[constr[[i]]] <- par[i]           
        notconstr[constr[[i]]] <- FALSE
    }
    par2[notconstr] <- par[(length(constr)+1):length(par)]
    ind <- 1    
    for(i in 1:length(obj)){
        obj[[i]]@par[obj[[i]]@est] <- par2[ind:(ind + sum(obj[[i]]@est) - 1)]
        ind <- ind + sum(obj[[i]]@est)                   
    }
    return(obj)
}

bfactor2mod <- function(model, J){        
    tmp <- tempfile('tempfile')
    unique <- sort(unique(model))
    index <- 1:J
    tmp2 <- sprintf(c('G =', paste('1-', J, sep='')))
    for(i in 1:length(unique)){
        ind <- index[model == unique[i]]
        comma <- rep(',', 2*length(ind))
        TF <- rep(c(TRUE,FALSE), length(ind))
        comma[TF] <- ind
        comma[length(comma)] <- ""
        tmp2 <- c(tmp2, c(paste('\nF', i, ' =', sep=''), comma))
    }
    cat(tmp2, file=tmp)
    model <- confmirt.model(file=tmp, quiet = TRUE)        
    unlink(tmp)
    return(model)
}

calcEMSE <- function(object, data, model, itemtype, fitvalues, constrain, parprior, verbose){          
    if(is(model, 'numeric') && length(model) > 1)
        model <- bfactor2mod(model, data)
    pars <- confmirt(data, model, itemtype=itemtype, pars=fitvalues, constrain=constrain, 
                     parprior=parprior, 
                     technical = list(BURNIN = 1, SEMCYCLES = 5, TOL = .01, 
                                    EMSE = TRUE), verbose = verbose)    
    for(i in 1:length(pars$pars))
        object@pars[[i]]@SEpar <- pars$pars[[i]]@SEpar                            
    object@information <- pars$info
    object@longpars <- pars$longpars
    return(object)
}

UpdateConstrain <- function(pars, constrain, invariance, nfact, nLambdas, J, ngroups, PrepList,
                            mixedlist, method)
{        
    #within group item constraints only
    for(g in 1:ngroups)  
        if(length(PrepList[[g]]$constrain) > 0)
            for(i in 1:length(PrepList[[g]]$constrain))
                constrain[[length(constrain) + 1]] <- PrepList[[g]]$constrain[[i]]        
    if('covariances' %in% invariance){ #Fix covariance accross groups (only makes sense with vars = 1)
        tmpmat <- matrix(NA, nfact, nfact)
        low_tri <- lower.tri(tmpmat)
        tmp <- c()                
        tmpmats <- tmpestmats <- matrix(NA, ngroups, nfact*(nfact+1)/2)
        for(g in 1:ngroups){
            tmpmats[g,] <- pars[[g]][[J + 1]]@parnum[(nfact+1):length(pars[[g]][[J + 1]]@parnum)] 
            tmpestmats[g,] <- pars[[g]][[J + 1]]@est[(nfact+1):length(pars[[g]][[J + 1]]@est)]
        }
        select <- colSums(tmpestmats) == ngroups
        for(i in 1:length(select))
            if(select[i])
                constrain[[length(constrain) + 1]] <- tmpmats[1:ngroups, i]
        
    }    
    if('slopes' %in% invariance){ #Equal factor loadings
        tmpmats <- tmpests <- list()
        for(g in 1:ngroups)
            tmpmats[[g]] <- tmpests[[g]] <- matrix(NA, J, nLambdas)                
        for(g in 1:ngroups){            
            for(i in 1:J){
                tmpmats[[g]][i,] <- pars[[g]][[i]]@parnum[1:nLambdas]
                tmpests[[g]][i,] <- pars[[g]][[i]]@est[1:nLambdas]
            }
        }
        for(i in 1:J){
            for(j in 1:nLambdas){
                tmp <- c()
                for(g in 1:ngroups){                    
                    if(tmpests[[1]][[i, j]])
                        tmp <- c(tmp, tmpmats[[g]][i,j])
                }
                constrain[[length(constrain) + 1]] <- tmp
            }
        }        
    }
    if('intercepts' %in% invariance){ #Equal item intercepts (and all other item pars)
        tmpmats <- tmpests <- list()
        for(g in 1:ngroups)
            tmpmats[[g]] <- tmpests[[g]] <- list() 
        for(g in 1:ngroups){            
            for(i in 1:J){
                ind <- (nLambdas+1):length(pars[[g]][[i]]@parnum)
                if(is(pars[[g]][[i]], 'dich')) ind <- ind[1:(length(ind)-2)]
                if(is(pars[[g]][[i]], 'partcomp')) ind <- ind[1:(length(ind)-1)]
                tmpmats[[g]][[i]] <- pars[[g]][[i]]@parnum[ind]
                tmpests[[g]][[i]] <- pars[[g]][[i]]@est[ind]
            }
        }
        for(i in 1:J){
            for(j in 1:length(tmpmats[[1]][[i]])){
                tmp <- c()
                for(g in 1:ngroups){                    
                    if(tmpests[[1]][[i]][j])
                        tmp <- c(tmp, tmpmats[[g]][[i]][j])
                }
                constrain[[length(constrain) + 1]] <- tmp
            }
        }
    }
    #accross item constraints for mixedmirt
    if(method == 'MIXED' && mixedlist$fixed.constrain){                    
        for(i in 1:pars[[1]][[1]]@nfixedeffects){
            tmp <- c()
            for(j in 1:J)
                tmp <- c(tmp, pars[[1]][[j]]@parnum[i])
            constrain[[length(constrain) + 1]] <- tmp           
        }        
    }
    return(constrain)
}

ReturnPars <- function(PrepList, itemnames, MG = FALSE){    
    parnum <- par <- est <- item <- parname <- gnames <- itemtype <- 
        lbound <- ubound <- c()                                    
    if(!MG) PrepList <- list(full=PrepList)                        
    for(g in 1:length(PrepList)){
        tmpgroup <- PrepList[[g]]$pars                                
        for(i in 1:length(tmpgroup)){
            if(i <= length(itemnames))
                item <- c(item, rep(itemnames[i], length(tmpgroup[[i]]@parnum)))
            parname <- c(parname, names(tmpgroup[[i]]@parnum))
            parnum <- c(parnum, tmpgroup[[i]]@parnum) 
            par <- c(par, tmpgroup[[i]]@par)
            est <- c(est, tmpgroup[[i]]@est)                    
            lbound <- c(lbound, tmpgroup[[i]]@lbound)
            ubound <- c(ubound, tmpgroup[[i]]@ubound)
        }
        item <- c(item, rep('GROUP', length(tmpgroup[[i]]@parnum)))                                
    }
    gnames <- rep(names(PrepList), each = length(est)/length(PrepList))
    ret <- data.frame(group=gnames, item = item, name=parname, parnum=parnum, value=par, 
                      lbound=lbound, ubound=ubound, est=est)
    ret
}

UpdatePrepList <- function(PrepList, pars, MG = FALSE){
    if(!MG) PrepList <- list(PrepList)    
    ind <- 1    
    for(g in 1:length(PrepList)){
        for(i in 1:length(PrepList[[g]]$pars)){ 
            for(j in 1:length(PrepList[[g]]$pars[[i]]@par)){
                PrepList[[g]]$pars[[i]]@par[j] <- pars[ind,'value']
                PrepList[[g]]$pars[[i]]@est[j] <- as.logical(pars[ind,'est'])                
                PrepList[[g]]$pars[[i]]@lbound[j] <- pars[ind,'lbound']
                PrepList[[g]]$pars[[i]]@ubound[j] <- pars[ind,'ubound']                
                ind <- ind + 1
            }
        }
    }    
    if(!MG) PrepList <- PrepList[[1]]
    return(PrepList)
}

#new gradient and hessian with priors
DerivativePriors <- function(x, grad, hess){
    if(any(!is.nan(x@n.prior.mu))){
        ind <- !is.na(x@n.prior.mu)            
        val <- x@par[ind]
        mu <- x@n.prior.mu[ind]
        s <- x@n.prior.sd[ind]
        h <- g <- rep(0, length(val))
        for(i in 1:length(val)){
            g[i] <- -(val[i] - mu[i])/(s[i]^2)
            h[i] <- -1/(s[i]^2)
        }
        grad[ind] <- grad[ind] + g
        if(length(val) == 1) hess[ind, ind] <- hess[ind, ind] + h
        else diag(hess[ind, ind]) <- diag(hess[ind, ind]) + h           
    }
    if(any(!is.nan(x@b.prior.alpha))){
        ind <- !is.na(x@b.prior.alpha)
        val <- x@par[ind]
        a <- x@b.prior.alpha[ind]
        b <- x@b.prior.beta[ind]
        bphess <- bpgrad <- rep(0, length(val))
        for(i in 1:length(val)){
            tmp <- betaprior(val[i], a[i], b[i])
            bpgrad[i] <- tmp$grad
            bphess[i] <- tmp$hess
        }
        if(length(val) == 1) hess[ind, ind] <- hess[ind, ind] + bpgrad
        else diag(hess[ind, ind]) <- diag(hess[ind, ind]) + bphess                
    }
    return(list(grad=grad, hess=hess))    
}

#new likelihood with priors
LL.Priors <- function(x, LL){
    if(any(!is.nan(x@n.prior.mu))){
        ind <- !is.nan(x@n.prior.mu)
        val <- x@par[ind]
        u <- x@n.prior.mu[ind]
        s <- x@n.prior.sd[ind]
        for(i in 1:length(val))            
            LL <- LL - log(dnorm(val[i], u[i], s[i]))
    }
    if(any(!is.nan(x@b.prior.alpha))){
        ind <- !is.nan(x@b.prior.alpha)
        val <- x@par[ind]
        a <- x@b.prior.alpha[ind]
        b <- x@b.prior.beta[ind]
        for(i in 1:length(val)){            
            tmp <- dbeta(val[i], a[i], b[i])
            LL <- LL - log(ifelse(tmp == 0, 1, tmp))
        }        
    }
    return(LL)
}

ItemInfo <- function(x, Theta, cosangle){    
    P <- ProbTrace(x, Theta)    
    dx <- DerivTheta(x, Theta)    
    info <- 0     
    cosanglefull <- matrix(cosangle, nrow(P), length(cosangle), byrow = TRUE)
    if(ncol(cosanglefull) < ncol(dx$grad[[1]]))
        cosanglefull <- cbind(cosanglefull, matrix(1, nrow(cosanglefull),
                                                   ncol(dx$grad[[1]]) - ncol(cosanglefull)))    
    for(i in 1:x@ncat){        
        dx$grad[[i]] <- matrix(rowSums(dx$grad[[i]] * cosanglefull))
        dx$hess[[i]] <- matrix(rowSums(dx$hess[[i]] * cosanglefull))
    }
    for(i in 1:x@ncat)        
        info <- info + ( (dx$grad[[i]])^2 / P[ ,i] - dx$hess[[i]])    
    return(info)
}

designMats <- function(covdata, fixed, Thetas, nitems, itemdesign = NULL, random = NULL, 
                       fixed.identical = FALSE){     
    fixed.design.list <- vector('list', nitems)     
    for(item in 1:nitems){
        if(item > 1 && fixed.identical){
            fixed.design.list[[item]] <- fixed.design.list[[1]]
            next            
        }                
        if(colnames(itemdesign)[1] != 'InTeRnAlUsElESsNaMe2'){
            dat <- data.frame(matrix(itemdesign[item, ], nrow(covdata), ncol(itemdesign), byrow=TRUE), 
                                   covdata, Thetas)
            colnames(dat) <- c(colnames(itemdesign), colnames(covdata), colnames(Thetas))
        } else dat <- data.frame(covdata, Thetas)        
        if(fixed == ~ 1) {
            fixed.design <- NULL
        } else fixed.design <- model.matrix(fixed, dat)[ ,-1, drop = FALSE]
        cn <- colnames(Thetas)
        CN <- colnames(fixed.design)
        drop <- rep(FALSE, length(CN))
        for(i in 1:ncol(Thetas))
            drop <- drop | CN == cn[i]        
        fixed.design.list[[item]] <- fixed.design[ , !drop, drop = FALSE]          
    }
    return(fixed.design.list)    
}

nameInfoMatrix <- function(info, correction, L, npars){
    #give info meaningful names for wald test    
    parnames <- names(correction)
    tmp <- outer(1:npars, rep(1, npars))
    matind <- matrix(0, ncol(tmp), nrow(tmp))
    matind[lower.tri(matind, diag = TRUE)] <- tmp[lower.tri(tmp, diag = TRUE)]    
    matind <- matind * L
    matind[matind == 0 ] <- NA  
    shortnames <- c()
    for(i in 1:length(correction)){
        keep <- is.na(matind[ , 1])
        while(all(keep)){            
            matind <- matind[-1, -1, drop = FALSE]        
            keep <- is.na(matind[ , 1])
        }            
        tmp <- paste0(parnames[i], paste0('.', matind[!keep, 1], collapse=''))
        shortnames <- c(shortnames, tmp)
        matind <- matind[keep, keep, drop = FALSE]        
    }
    colnames(info) <- rownames(info) <- shortnames
    return(info)
}

maketabData <- function(stringfulldata, stringtabdata, group, groupNames, nitem, K, itemloc,
                        Names, itemnames){    
    tabdata2 <- lapply(strsplit(stringtabdata, split='/'), as.numeric)
    tabdata2 <- do.call(rbind, tabdata2)
    tabdata2[tabdata2 == 99999] <- NA
    tabdata <- matrix(0, nrow(tabdata2), sum(K))
    for(i in 1:nitem){
        uniq <- sort(na.omit(unique(tabdata2[,i])))        
        for(j in 1:length(uniq))
            tabdata[,itemloc[i] + j - 1] <- as.numeric(tabdata2[,i] == uniq[j])        
    }        
    colnames(tabdata) <- Names
    colnames(tabdata2) <- itemnames
    ret1 <- ret2 <- vector('list', length(groupNames))    
    for(g in 1:length(groupNames)){
        tmpstringdata <- stringfulldata[group == groupNames[g]]
        Freq <- numeric(nrow(tabdata))
        for(i in 1:length(stringtabdata))
            Freq[i] <- sum(tmpstringdata == stringtabdata[i])
        ret1[[g]] <- cbind(tabdata, Freq)
        ret2[[g]] <- cbind(tabdata2, Freq)
    }
    ret <- list(tabdata=ret1, tabdata2=ret2)
    ret
}

makeopts <- function(method = 'MHRM', draws = 2000, calcLL = TRUE, quadpts = NaN, 
                     rotate = 'varimax', Target = NaN, SE = TRUE, verbose = TRUE, 
                     SEtol = .01, nested.mod = NULL, grsm.block = NULL, D = 1.702, 
                     rsm.block = NULL, calcNull = TRUE, cl = NULL, BFACTOR = FALSE, 
                     technical = list(), use = 'pairwise.complete.obs', debug = FALSE)
{    
    opts <- list()
    opts$method = method
    opts$draws = draws
    opts$calcLL = calcLL
    opts$quadpts = quadpts 
    opts$rotate = rotate 
    opts$Target = Target
    opts$SE = SE 
    opts$verbose = verbose 
    opts$SEtol = SEtol
    opts$nested.mod = nested.mod
    opts$grsm.block = grsm.block
    opts$D = D 
    opts$rsm.block = rsm.block
    opts$calcNull = calcNull
    opts$cl = cl 
    opts$BFACTOR = BFACTOR 
    opts$debug = debug    
    opts$technical <- technical
    opts$use <- use
    opts$MAXQUAD <- ifelse(is.null(technical$MAXQUAD), 10000, technical$MAXQUAD)
    opts$MSTEPMAXIT <- ifelse(is.null(technical$MSTEPMAXIT), 15, technical$MSTEPMAXIT)        
    opts$NCYCLES <- ifelse(is.null(technical$NCYCLES), 2000, technical$NCYCLES)
    if(opts$method == 'EM')
        opts$NCYCLES <- ifelse(is.null(technical$NCYCLES), 300, technical$NCYCLES)
    opts$BURNIN <- ifelse(is.null(technical$BURNIN), 150, technical$BURNIN)
    opts$SEMCYCLES <- ifelse(is.null(technical$SEMCYCLES), 50, technical$SEMCYCLES)
    opts$KDRAWS  <- ifelse(is.null(technical$KDRAWS), 1, technical$KDRAWS)
    opts$TOL <- ifelse(is.null(technical$TOL), .001, technical$TOL)
    set.seed(12345)
    if(!is.null(technical$set.seed)) set.seed(technical$set.seed)    
    opts$gain <- c(0.05,0.5,0.004)
    if(!is.null(technical$gain)){
        if(length(technical$gain) == 3 && is.numeric(technical$gain))
            opts$gain <- technical$gain
    }	 
    opts$NULL.MODEL <- ifelse(is.null(technical$NULL.MODEL), FALSE, TRUE)
    opts$USEEM <- ifelse(method == 'EM', TRUE, FALSE)    
    if(!is.null(cl)) require(parallel)    
    return(opts)
}
