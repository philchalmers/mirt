# theta combinations
thetaComb <- function(theta, nfact)
{
	if (nfact == 1) Theta <- matrix(theta)
	else if (nfact == 2) Theta <- expand.grid(theta,theta)   
	else if (nfact == 3) Theta <- expand.grid(theta,theta,theta)  
	else if (nfact == 4) Theta <- expand.grid(theta,theta,theta,theta)
	else if (nfact == 5) Theta <- expand.grid(theta,theta,theta,theta,theta)        
	else if (nfact == 6) Theta <- expand.grid(theta,theta,theta,theta,theta,theta)        	
	if(nfact > 6) stop('Are you crazy?!?!? That\'s way too many factors for this quandrature method.')
	Theta <- as.matrix(Theta)	
	return(Theta)     
}

# start values
start.values <- function(fulldata, guess, Rpoly, nfact=2, bfactor = FALSE, nowarn = TRUE)
{	  	
	if (bfactor){ 
		suppressWarnings(FA <- psych::fa(Rpoly, 1, warnings = !nowarn))
		loads <- unclass(FA$load)
		cs <- sqrt(abs(FA$u))      
		dstart <- qnorm(colMeans(fulldata))/cs
		astart <- loads/cs
		startvalues <- cbind(astart,astart/2,dstart)
	} else {    
		suppressWarnings(FA <- psych::fa(Rpoly,nfact,rotate = 'none', warnings= !nowarn))	
		loads <- unclass(loadings(FA))
		u <- FA$unique
		u[u < .001 ] <- .2
		cs <- sqrt(u)
		dstart <- qnorm(colMeans(fulldata))/cs
		astart <- loads/cs
		startvalues <- cbind(astart,dstart)
	}  	
	startvalues
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

# MAP scoring for mirt
MAP.mirt <- function(Theta, a, d, guess, upper, patdata, itemloc, ML=FALSE)
{	
	itemtrace <- rep(0, ncol=length(patdata))
	Theta <- matrix(Theta, 1)
	for (i in 1:length(guess)){
		if(length(d[[i]]) == 1){
			itemtrace[itemloc[i]] <- P.mirt(a[i, ], d[[i]], Theta, guess[i], upper[i]) 
			itemtrace[itemloc[i] + 1] <- 1.0 - itemtrace[itemloc[i]]
		} else {
			itemtrace[itemloc[i]:(itemloc[i+1] - 1)] <- 
				P.poly(a[i, ], d[[i]], Theta, TRUE)	
		}
	}		
	L <- sum(log(itemtrace)[as.logical(patdata)])
	mu <- 0
	sigma <- 1
    L <- ifelse(ML, -L, (-1)*(L + sum(log(exp(-0.5*((Theta - mu)/sigma)^2)))))
	L  
}  

# MAP scoring for bfactor
MAP.bfactor <- function(Theta, a, d, guess, upper, patdata, logicalfact, itemloc, ML=FALSE)
{	
	itemtrace <- rep(0, ncol=length(patdata))
	Theta <- matrix(Theta, 1)
	for (i in 1:length(guess)){
		if(length(d[[i]]) == 1){
			itemtrace[itemloc[i]] <- P.mirt(a[i, logicalfact[i, ]], d[[i]], Theta, guess[i], upper[i]) 
			itemtrace[itemloc[i] + 1] <- 1.0 - itemtrace[itemloc[i]]
		} else {
			itemtrace[itemloc[i]:(itemloc[i+1] - 1)] <- 
				P.poly(a[i, logicalfact[i, ]], d[[i]], Theta, TRUE)	
		}
	}		
	L <- sum(log(itemtrace)[as.logical(patdata)])
	mu <- 0
	sigma <- 1
    L <- ifelse(ML, -L, (-1)*(L + sum(log(exp(-0.5*((Theta - mu)/sigma)^2)))))
	L  
}  

# Estep for mirt
Estep.mirt <- function(pars, tabdata, Theta, prior, itemloc) 
{       
	nfact <- ncol(Theta)
	nquad <- nrow(Theta)	
	r <- tabdata[ ,ncol(tabdata)]
	X <- tabdata[ ,1:(ncol(tabdata) - 1)]	
	itemtrace <- matrix(0, ncol=ncol(X), nrow=nrow(Theta))	
	for (i in 1:length(pars))
	    itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)
    retlist <- .Call("Estep", itemtrace, prior, X, nfact, r)    
	return(retlist)
} 

# Estep for bfactor
Estep.bfactor <- function(pars, tabdata, Theta, prior, specific, sitems, itemloc) 
{	    
	nfact <- pars[[1]]@nfact
	J <- length(pars)
	nquad <- nrow(Theta)		
	r <- tabdata[ ,ncol(tabdata)]
	X <- tabdata[ ,1:(ncol(tabdata) - 1)]	
	itemtrace <- matrix(0, ncol=ncol(X), nrow=nrow(Theta))	
	for (i in 1:J)
	    itemtrace[ ,itemloc[i]:(itemloc[i+1] - 1)] <- ProbTrace(x=pars[[i]], Theta=Theta)			
	retlist <- .Call("Estepbfactor", itemtrace, prior, X, r, sitems)	
	r1 <- matrix(0, nrow(Theta), ncol(X))	
	for (i in 1:J)
	    r1[ ,itemloc[i]:(itemloc[i+1]-1)] <- 		
			retlist$r1[ ,itemloc[i]:(itemloc[i+1]-1) + (specific[i] - 1)*ncol(X) ]		
	return(list(r1=r1, expected=retlist$expected))	
}      

Mstep.mirt <- function(par, obj, Theta, prior, constr = list()){ 
    if(length(constr) < 1){
        obj@par[obj@est] <- par    
        ret <- LogLik(x=obj, Theta=Theta)                
    } else {        
        obj <- reloadConstr(par=par, constr=constr, obj=obj)        
        ret <- 0
        for(i in 1:length(obj))            
            ret <- ret + LogLik(x=obj[[i]], Theta=Theta)               
    }
    return(ret)
}

# MH sampler for theta values
draw.thetas <- function(theta0, pars, fulldata, itemloc, cand.t.var, prior.t.var, 
                        prior.mu, prodlist) 
{ 	    
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
	if(!is.null(prodlist)){
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
betaprior <- function(a,b,g,W=20)
{
	a <- a + (1-g)*W
	b <- b + g*W
	grad <- ((a-1) * g^(a-1) * (1-g)^(b-1) - (b-1)*g^(a-1)*(1-g)^(b-1))/ 
		(g^(a-1) * (1-g)^(b-1))
	hess <- -((g^(a-1)*(a-1)^2*(1-g)^(b-1)/g^2 - g^(a-1)*(a-1)*(1-g)^(b-1)/g^2 
		- 2*g^(a-1)*(a-1)*(1-g)^(b-1)*(b-1)/(g*(1-g)) + g^(a-1)*(1-g)^(b-1)*(b-1)^2/(1-g)^2 
		- g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g)^2)/(g^(a-1)*(1-g)^(b-1))	
		- ((g^(a-1)*(a-1)*(1-g)^(b-1)/g-g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g))*(a-1)/(g^(a-1)*(1-g)^(b-1)*g))
		+ ((g^(a-1)*(a-1)*(1-g)^(b-1)/g-g^(a-1)*(1-g)^(b-1)*(b-1)/(1-g))*(b-1)/(g^(a-1)*(1-g)^(b-1)*(1-g))))
	return(list(g=grad, h=hess))	
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

# Extract model matricies and values for user specified confmirt.model()
model.elements <- function(model, factorNames, itemtype, nfactNames, nfact, J, K, fulldata, 
                           itemloc, data, N, guess, upper, itemnames, exploratory, constrain, 
                           startvalues, freepars, parprior, parnumber)
{       
    hasProdTerms <- ifelse(nfact == nfactNames, FALSE, TRUE)
    prodlist <- NULL
    if(hasProdTerms){
        tmp <- factorNames[grepl('\\(',factorNames)]
        tmp2 <- factorNames[!grepl('\\(',factorNames)] 
        tmp <- gsub("\\(","",tmp)	
        tmp <- gsub("\\)","",tmp)
        tmp <- gsub(" ","",tmp)
        prodlist <- strsplit(tmp,"\\*")
        for(j in 1:length(prodlist)){
            for(i in 1:nfact)
                prodlist[[j]][prodlist[[j]] == tmp2[[i]]] <- i		
            prodlist[[j]] <- as.numeric(prodlist[[j]])	
        }		
    }  
    #slopes specification
    estlam <- matrix(FALSE, ncol = nfactNames, nrow = J)	
    for(i in 1:nfactNames){
        tmp <- model[model[ ,1] == factorNames[i],2]
        if(any(regexpr(",",tmp)))
            tmp <- strsplit(tmp,",")[[1]]
        popout <- c()	
        for(j in 1:length(tmp)){
            if(regexpr("-",tmp[j]) > 1){
                popout <- c(popout,j)
                tmp2 <- as.numeric(strsplit(tmp[j],"-")[[1]])
                tmp2 <- as.character(tmp2[1]:tmp2[2])
                tmp <- c(tmp,tmp2)
            }
        }
        if(length(popout != 0))	
            estlam[as.numeric(tmp[-popout]),i] <- TRUE
        else 
            estlam[as.numeric(tmp),i] <- TRUE
    }
    lambdas <- ifelse(estlam, .5, 0)	  
    #INT
    cs <- sqrt(abs(1-rowSums(lambdas^2)))	
    zetas <- list()
    loc <- 1	
    for(i in 1:J){        
        if(K[i] == 2){
            div <- ifelse(cs[i] > .25, cs[i], .25)		
            zetas[[i]] <- qnorm(mean(fulldata[,itemloc[i]]))/div            
        } else {			
            temp <- table(data[,i])[1:(K[i]-1)]/N
            temp <- cumsum(temp)
            div <- ifelse(cs[i] > .25, cs[i], .25)		
            zetas[[i]] <- qnorm(1 - temp)/div	            
        }		
    }    
    estzetas <- list()        
    for(i in 1:J)
        estzetas[[i]] <- length(zetas[[i]])                
    #COV
    find <- 1:nfact
    estgcov <- matrix(FALSE,nfact,nfact)    
    if(any(model[,1] == 'COV')){
        tmp <- model[model[,1] == 'COV',2]		
        tmp <- strsplit(tmp,",")[[1]]
        tmp <- gsub(" ","",tmp)
        for(i in 1:length(tmp)){            
            tmp2 <- strsplit(tmp[i],"*",fixed=TRUE)[[1]]				
            ind1 <- find[tmp2[1] == factorNames]
            ind2 <- find[tmp2[2] == factorNames]
            estgcov[ind2,ind1] <- TRUE            	
        }
    }
    gcov <- ifelse(estgcov,.1,0) 
    diag(gcov) <- 1	  
    #MEAN
    gmeans <- rep(0, nfact)
    estgmeans <- rep(FALSE, nfact)
#   #PRIOR, 1 == norm, 2== beta
#   parpriors <- list()
#   parpriorscount <- 1
#   if(sum(estGuess) > 0){
#     for(i in 1:J){
#       if(estGuess[i]){
#         a <- guess[i] * guess.prior.n[i]
#         b <- (1 - guess[i]) * guess.prior.n
#         parpriors[[parpriorscount]] <- c(2,guessind[i],a,b)						
#         parpriorscount <- parpriorscount + 1			
#       }
#     }
#   }		
    ret <- LoadPars(itemtype=itemtype, itemloc=itemloc, lambdas=lambdas, zetas=zetas, guess=guess, upper=upper,
                 fulldata=fulldata, J=J, K=K, nfact=nfact, constrain=constrain, nfactNames=nfactNames,
                    startvalues=startvalues, freepars=freepars, parprior=parprior, parnumber=parnumber,
                    estLambdas=estlam)  
    parnumber <- ret[[length(ret)]]@parnum[length(ret[[length(ret)]]@parnum)]
    ret[[length(ret) + 1]] <- LoadGroupPars(gmeans=gmeans, gcov=gcov, estgmeans=estgmeans, 
                                            estgcov=estgcov, parnumber=parnumber+1)
    attr(ret, 'prodlist') <- prodlist    
    return(ret)    
}

# Take long parameter form and return list of pars for polymirt (obsolete)
sortPars <- function(pars, indlist, nfact, estGuess)
{
	lambdas <- matrix(pars[indlist$lamind],ncol=nfact,byrow=TRUE)	
	J <- nrow(lambdas)		
	zetas <- list()
	for(i in 1:J)
		zetas[[i]] <- pars[indlist$zetaind[[i]]]
	guess <- upper <- rep(0,J)
	guess[estGuess] <- pars[indlist$gind]		
	
	return(list(lambdas=lambdas, zetas=zetas, guess=guess))
}

# Take long parameter form and return list of pars for confmirt
sortParsConfmirt <- function(pars, indlist, nfact, nfactNames)
{
	J <- length(indlist$guessind)
	lambdas <- matrix(pars[indlist$lamind],J,nfactNames,byrow=TRUE)
	zetas <- list()
	for(i in 1:J)
		zetas[[i]] <- pars[indlist$zetaindlist[[i]]]
	guess <- pars[indlist$guessind]
    upper <- pars[indlist$upperind]
	mu <- pars[indlist$meanind]
	sig <- matrix(0, nfact, nfact)
	sig[lower.tri(sig, diag=TRUE)] <- pars[indlist$covind]
	if(nfact > 1)
		sig <- sig + t(sig) - diag(diag(sig))							
	
	return(list(lambdas=lambdas, zetas=zetas, guess=guess, upper=upper, mu=mu, sig=sig))
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

# Rebuild parameters given a list
rebuildPars <- function(p, pars)
{
	names(p) <- NULL
	pars2 <- pars
	pars2$lambdas <- matrix(p[1:length(pars$lambdas)], ncol=ncol(pars$lambdas), 
		nrow=nrow(pars$lambdas)) 
	ind1 <- length(pars$lambdas) + 1
	for(i in 1:length(pars$zetas) ){
		ind2 <- ind1 + length(pars$zetas[[i]]) - 1
		pars2$zetas[[i]] <- p[ind1:ind2]
		ind1 <- ind1 + length(pars$zetas[[i]])
	}			
	return(pars2)
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
    
test_info <- function(a, d, Theta, Alist, guess, upper, K){
    infolist <- list()    
    for(cut in 1:length(Alist)){
        A <- Alist[[cut]]
        info <- rep(0,nrow(Theta))
        for(j in 1:length(K)){
            if(K[j] > 2){
                P <- P.poly(a[j,], d[[j]], Theta, itemexp = FALSE)    	
                for(i in 1:K[j]){
                    w1 <- P[,i]*(1-P[,i])*A[j]
                    w2 <- P[,i+1]*(1-P[,i+1])*A[j]
                    I <- ((w1 - w2)^2) / (P[,i] - P[,i+1]) * P[,i]
                    info <- info + I
                }
            } else {
                P <- P.mirt(a[j,], d[[j]], Theta, guess[j], upper[j])
                Pstar <- P.mirt(a[j,], d[[j]], Theta, 0)
                info <- info + A[j]^2 * P * (1-P) * Pstar/P ###FIXME: might need new 4PL info
            }			
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

LoadPars <- function(itemtype, itemloc, lambdas, zetas, guess, upper, fulldata, J, K, nfact, 
                     constrain, startvalues, freepars, parprior, parnumber, bfactor = NULL, 
                     nfactNames = NULL, estLambdas=NULL){    
    pars <- list()   
    if(is.null(nfactNames)) nfactNames <- nfact
    if(is.null(estLambdas)) estLambdas <- matrix(TRUE, J, nfactNames)
    BFACTOR <- FALSE
    if(!is.null(bfactor)){
        estLambdas <- bfactor
        BFACTOR <- TRUE
    }                
    constr <- c()
    if(!is.null(constrain) && is.list(constrain)) 
        for(i in 1:length(constrain))
            constr <- c(constr, constrain[[i]])
    constr <- unique(constr)
    for(i in 1:J){
        tmp <- c(itemloc[i]:(itemloc[i+1] - 1)) #item location 
        if(itemtype[i] == 'NullModel' && K[i] == 2) 
            pars[[i]] <- new('dich', par=c(0,zetas[[i]],0,1), nfact=1, bfactor=BFACTOR,
                             dat=fulldata[ ,tmp], est=c(FALSE,TRUE,FALSE,FALSE), constr=FALSE)
        
        if(itemtype[i] == 'NullModel' && K[i] > 2) 
            pars[[i]] <- new('graded', par=c(0,zetas[[i]]), nfact=1, ncat=K[i], bfactor=BFACTOR,
                             dat=fulldata[ ,tmp], est=c(FALSE,rep(TRUE,K[i]-1)), constr=FALSE)
        
        if(any(itemtype[i] == c('2PL', '3PL', '3PLu', '4PL'))){ 
            pars[[i]] <- new('dich', par=c(lambdas[i,], zetas[[i]], guess[i], upper[i]),
                             nfact=nfactNames, dat=fulldata[ ,tmp], constr=FALSE, bfactor=BFACTOR)                    
            estpars <- c(estLambdas[i, ], TRUE, FALSE, FALSE) 
            if(any(itemtype[i] == c('3PL', '4PL'))) estpars[length(estpars)-1] <- TRUE
            if(any(itemtype[i] == c('3PLu', '4PL'))) estpars[length(estpars)] <- TRUE            
            pars[[i]]@est <- estpars
            tmp2 <- parnumber:(parnumber + length(estpars) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE
            names(tmp2) <- c(paste('a', 1:nfact, sep=''), 'd', 'g','u')
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(estpars)            
        }
        
        if(itemtype[i] == 'graded'){
            pars[[i]] <- new('graded', par=c(lambdas[i,], zetas[[i]]), nfact=nfactNames, ncat=K[i],
                             dat=fulldata[ ,tmp], constr=FALSE, bfactor=BFACTOR)            
            estpars <- c(estLambdas[i, ], rep(TRUE, K[i]-1))
            pars[[i]]@est <- estpars
            tmp2 <- parnumber:(parnumber + length(estpars) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE
            names(tmp2) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:(K[i]-1), sep=''))
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(estpars)
        }
        
        if(itemtype[i] == 'gpcm'){            
            pars[[i]] <- new('gpcm', par=c(lambdas[i,], zetas[[i]]), nfact=nfactNames, ncat=K[i],
                             dat=fulldata[ ,tmp], constr=FALSE, bfactor=BFACTOR)
            estpars <- c(estLambdas[i, ], rep(TRUE, K[i]))
            #identifiction constraints
            estpars[nfact+1] <- FALSE
            pars[[i]]@par[nfact+1] <- 0
            pars[[i]]@est <- estpars
            tmp2 <- parnumber:(parnumber + length(estpars) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE
            names(tmp2) <- c(paste('a', 1:nfact, sep=''), paste('d', 0:(K[i]-1), sep=''))
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(estpars)
        }        
        
        if(itemtype[i] == 'nominal'){
            pars[[i]] <- new('nominal', par=c(rep(.5, nfact), 0, rep(.5, K[i] - 2), K[i]-1, rep(0, K[i])), 
                             nfact=nfactNames, ncat=K[i], dat=fulldata[ ,tmp], constr=FALSE, bfactor=BFACTOR)
            estpars <- c(estLambdas[i, ], rep(TRUE, length(pars[[i]]@par) - nfact))
            #identifiction constraints
            estpars[c(nfact+1, nfact+ K[i], nfact + K[i] + 1)] <- FALSE
            pars[[i]]@par[c(nfact + 1, nfact + K[i] + 1)] <- 0
            pars[[i]]@par[nfact + K[i]] <- K[i] - 1
            pars[[i]]@est <- estpars
            tmp2 <- parnumber:(parnumber + length(estpars) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE
            names(tmp2) <- c(paste('a', 1:nfact, sep=''), paste('ak', 0:(K[i]-1), sep=''), 
                             paste('d', 0:(K[i]-1), sep=''))
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(estpars)
        } 
        
        if(any(itemtype[i] == c('PC2PL','PC3PL'))){
            pars[[i]] <- new('partcomp', par=c(rep(.5, nfact), rep(-1, nfact), 0, 1), 
                             nfact=nfactNames, dat=fulldata[ ,tmp], constr=FALSE, bfactor=BFACTOR)
            estpars <- c(estLambdas[i, ], estLambdas[i, ], FALSE, FALSE)
            if(itemtype[i] == 'PC3PL') estpars[length(estpars) - 1] <- TRUE
            pars[[i]]@est <- estpars
            tmp2 <- parnumber:(parnumber + length(estpars) - 1)
            if(length(intersect(tmp2, constr)) > 0 ) pars[[i]]@constr <- TRUE
            names(tmp2) <- c(paste('a', 1:nfact, sep=''), paste('d', 1:nfact, sep=''), 'g','u')
            pars[[i]]@parnum <- tmp2
            parnumber <- parnumber + length(estpars)
        }
    }
    if(!is.null(startvalues)) {
        if(startvalues != 'index')
            for(i in 1:J)
                pars[[i]]@par <- startvalues[[i]]        
    }
    if(!is.null(freepars)) {
        if(freepars != 'index')
            for(i in 1:J)
                pars[[i]]@est <- freepars[[i]]        
    }
    attr(pars, 'uniqueconstr') <- constr 
    return(pars)
}

LoadGroupPars <- function(gmeans, gcov, estgmeans, estgcov, parnumber){
    nfact <- length(gmeans)
    fn <- paste('COV_', 1:nfact, sep='')
    FNCOV <- outer(fn, 1:nfact, FUN=paste, sep='')
    FNMEANS <- paste('MEAN_', 1:nfact, sep='')  
    tri <- lower.tri(gcov, diag=TRUE)
    par <- c(gmeans, gcov[tri])
    parnum <- parnumber:(parnumber + length(par) - 1)
    names(parnum) <- c(FNMEANS,FNCOV[tri])
    ret <- new('GroupPars', par=par, est=c(estgmeans,estgcov[tri]), nfact=nfact, 
               parnum=parnum)
    ret    
}

#change long pars for groups into mean in sigma
ExtractGroupPars <- function(x){
    nfact <- x@nfact
    gmeans <- x@par[1:nfact]
    gmeans <- x@par[1:nfact]
    tmp <- x@par[-(1:nfact)]
    gcov <- matrix(0, nfact, nfact)
    gcov[lower.tri(gcov, diag=TRUE)] <- tmp
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
