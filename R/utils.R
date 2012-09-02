# theta combinations
thetaComb <- function(theta, nfact)
{
	if (nfact == 1) Theta <- matrix(theta)
	else if (nfact == 2) Theta <- expand.grid(theta,theta)   
	else if (nfact == 3) Theta <- expand.grid(theta,theta,theta)  
	else if (nfact == 4) Theta <- expand.grid(theta,theta,theta,theta)
	else if (nfact == 5) Theta <- expand.grid(theta,theta,theta,theta,theta)        	
	else if (nfact == 6) Theta <- expand.grid(theta,theta,theta,theta,theta,theta)
	if(nfact > 6) stop('Are you crazy?!?!? That\'s way too many factors for this quandrature method.
                       Try using confmirt() instead for better accuracy')
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
                        prior.mu, prodlist, debug) 
{         
    if(debug == 'draw.thetas') browser()
    tol <- 1e-8
    N <- nrow(fulldata)
    J <- length(pars) - 1
    nfact <- ncol(theta0)					
    unif <- runif(N)
    if(nfact > 1)		
        theta1 <- theta0 + mvtnorm::rmvnorm(N,prior.mu, 
                                            diag(rep(cand.t.var,ncol(theta0)))) 
    else
        theta1 <- theta0 + rnorm(N,prior.mu,sqrt(cand.t.var))							
    log_den0 <- mvtnorm::dmvnorm(theta0,prior.mu,prior.t.var,log=TRUE)
    log_den1 <- mvtnorm::dmvnorm(theta1,prior.mu,prior.t.var,log=TRUE)		
    if(length(prodlist) > 0){
        theta0 <- prodterms(theta0,prodlist)
        theta1 <- prodterms(theta1,prodlist)	
    }	
    itemtrace0 <- itemtrace1 <- matrix(0, ncol=ncol(fulldata), nrow=nrow(theta0))    
    for (i in 1:J){
        itemtrace0[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=theta0)
        itemtrace1[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=theta1)        
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
        theta1 <- theta1[ ,1:(pars[[1]]@nfact - length(prodlist)), drop=FALSE]
    attr(theta1, "Proportion Accepted") <- sum(accept)/N 				
    attr(theta1, "log.lik") <- log.lik	
    return(theta1) 
}

# Rotation function
Rotate <- function(F, rotate, Target = NULL, ...)
{	
    if(ncol(F) == 1) rotF <- list()
	if(rotate == 'promax'){
        rotF <- psych::Promax(F)
        rotF$orthogonal <- FALSE
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
cormod <- function(fulldata, K, guess, smooth = TRUE) 
{  
	fulldata <- as.matrix(fulldata) 
	nitems <- ncol(fulldata)         
	cormat <- cor(fulldata)      
	if (any(guess > 0)){
		for (i in 1:nitems){
			for (j in 1:nitems){
				if (i < j & K[i] == 2 & K[j] == 2 & guess[i]!= 0 ){         
					g1 <- guess[i]
					g2 <- guess[j]
					tabp <- tab <- table(fulldata[ ,i],fulldata[ ,j])/length(fulldata[ ,i])
					w1 <- (1 - g1)
					w2 <- (1 - g2)
					tabp[1,1] <- tab[1,1]/(w1*w2)
					tabp[1,2] <- (w2*tab[1,2] - g2*tab[1,1])/(w1*w2)
					tabp[2,1] <- (w1*tab[2,1] - g1*tab[1,1])/(w1*w2)
					tabp[2,2] <- 1 - tabp[1,1] - tabp[1,2] - tabp[2,1]
					tabp <- round(tabp*length(fulldata[ ,i]))
					if(any(tabp < 0)) next	
					cormat[i,j] <- cormat[j,i] <- 
						abs(psych::phi(tabp,6))^(1/1.15)*sign(psych::phi(tabp,6))          		  
				} 
				if(i < j & K[i] == 2 & K[j] > 2 & guess[i]!= 0) 
					cormat[i,j] <- cormat[j,i] <- abs(cormat[i,j])^(1/1.15) * sign(cormat[i,j])
			}	
		}      
	} 
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
    F <- so$rotF %*% t(chol(so$fcor))
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

calcEMSE <- function(object, data, model, constrain, parprior, verbose){  
    startvalues <- coef(object, allpars = T)  
    if(is(model, 'numeric') && length(model) > 1){
        J <- ncol(data)
        tmp <- tempfile('tempfile')
        unique <- unique(model)
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
        model <- confmirt.model(tmp, quiet = TRUE)        
        unlink(tmp)
    }
    pars <- confmirt(data, model, startvalues=startvalues, constrain=constrain, parprior=parprior, 
                    technical = list(BURNIN = 1, SEMCYCLES = 5, TOL = .01, 
                                    EMSE = TRUE), verbose = verbose)
    for(i in 1:length(pars))
        object@pars[[i]]@SEpar <- pars[[i]]@SEpar    
    return(object)
}

UpdateConstrain <- function(pars, constrain, invariance, nfact, nLambdas, J, ngroups){
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
    return(constrain)
}
