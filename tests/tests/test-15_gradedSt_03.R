# FOLLOWING is copied from modified package mirt

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'grmSt', Theta = 'matrix'),
    definition = function(x, Theta){
        #if(nrow(x@fixed.design) > 1L && useDesign)
        #    Theta <- cbind(x@fixed.design, Theta)
        #SOURCE from ProbTrace(graded)

        #P.grm <- P.poly2(x@par[-x@nfact], Theta=Theta[,-x@nfact, drop=F], itemexp=itemexp, ot=ot)

        th1 = Theta[,1]; th2 = Theta[,2]

        # Differ from P.egrm
        ncat = x@ncat/2
        #ncat = x@ncat
        a = x@par[1]
        d = x@par[3:(ncat+1)]
        #
        if (all(d == sort(d))) {

            D.star = matrix(d, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
            # a.xi included
            TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
            A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
            P = 1/(1+exp(-1*(A*(TH1-D.star))))
            P.star=cbind(1, P)-cbind(P, 0)

            a.st <- x@par[x@nfact]
            p.st <- matrix(th2, ncol=ncat, nrow=nrow(Theta))
            ##p.repeat <- 2*abs(1/(1+exp(-a.st*p.st))-1/2)
            #p.st2 <- exp(-a.st*p.st)
            #p.repeat <- (p.st2-1)/(p.st2+1)
            p.repeat <- tanh(-a.st*p.st/2)
            q.repeat <- ifelse(p.repeat<0, -p.repeat, 0)
            p.repeat = ifelse(p.repeat > 0, p.repeat, 0)


            P <- cbind(P.star * (1 - p.repeat) * (1 - q.repeat) + q.repeat/(ncat-1),
                       P.star * (1 - p.repeat) * (1 - q.repeat) + p.repeat)
            #P <- cbind(P.star * q.repeat, P.star * q.repeat + p.repeat)

            P <- ifelse(P < 1e-20, 1e-20, P)
            P <- ifelse(P > (1 - 1e-20), (1 - 1e-20), P)


        } else {

            P <- matrix(0, nrow=nrow(Theta), ncol=ncat*2)}




        return(P)
    }
)


# FOLLOWING is copied from test-13_egrm10a.R

createR=function(n.subj, n.item, ncat, Theta) {

    #

    a = runif(n.item, 1, 2)
    d = matrix(runif(n.item*(ncat-1), -3, 3), ncol=ncat-1)
    d = t(apply(d, 1, sort))

        #P.grm from P2.egrm Theta -> th1
    P.grm <- function(par, th1, ncat) {
        #th1 = Theta[,1]; xi1 = Theta[,2];
        a = par[1]
        d = par[2:ncat]
        d.mean=mean(d);
        D.star = matrix(exp(0), nrow=length(th1), ncol=ncat-1) *
            matrix((d - d.mean) + d.mean, nrow=length(th1), ncol=ncat-1, byrow=T)
        TH1 = matrix(th1, nrow=length(th1), ncol=ncat-1)
        A = matrix(a, nrow=length(th1), ncol=ncat-1)
        P = 1/(1+exp(-1*(A*(TH1-D.star))))

        return(P)
    }

    # P.grmSt for creating Items adopted from ABOVE!
    P.grmSt = function(par, Theta, prev, ncat){
        #if(nrow(x@fixed.design) > 1L && useDesign)
        #    Theta <- cbind(x@fixed.design, Theta)
        #SOURCE from ProbTrace(graded)

        #P.grm <- P.poly2(x@par[-x@nfact], Theta=Theta[,-x@nfact, drop=F], itemexp=itemexp, ot=ot)

        stopifnot(length(prev)==nrow(Theta))

        th1 = Theta[,1]; th2 = Theta[,2]



        # Differ from P.egrm
        #ncat = ncat/2
        ncat=ncat

        a = par[1]; a.st <- par[2]
        d = par[3:(ncat+1)]
        #
        #if (all(d == sort(d))) {
        D.star = matrix(d, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
        # a.xi included
        TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
        A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
        P = 1/(1+exp(-1*(A*(TH1-D.star))))
        P.star=cbind(1, P)-cbind(P, 0)


        p.st <- matrix(th2, ncol=ncat, nrow=nrow(Theta))
        ##p.repeat <- 2*abs(1/(1+exp(-a.st*p.st))-1/2)
        #p.st2 <- exp(-a.st*p.st)
        #p.repeat <- (p.st2-1)/(p.st2+1)
        p.repeat <- tanh(-a.st*p.st/2)
        q.repeat <- ifelse(p.repeat<0, -p.repeat, 0)
        p.repeat = ifelse(p.repeat > 0, p.repeat, 0)

        P <- cbind(P.star * (1 - p.repeat) * (1 - q.repeat) + q.repeat/(ncat-1), P.star * (1 - p.repeat) * (1 - q.repeat) + p.repeat)
        #P <- cbind(P.star * q.repeat, P.star * q.repeat + p.repeat)


        pos.item = matrix(NA, nrow(Theta), ncat)
        for (i in 1:length(prev)) {
            pos.item[i,] = 1:ncat + ncat*(1:ncat == prev[i])
            # for prev[i] == 1, (1,0,0,0), for prev[i] == 2, (0,1,0,0)
        }

        P2 <- matrix(NA, nrow(Theta), ncat)
        for (i in 1:length(prev)) {
            P2[i,]=P[i,pos.item[i,]]
        }

        P2.prob <- 1-t(apply(P2, 1, cumsum))
        P2.prob <- P2.prob[,-ncol(P2.prob)]

        return(P2.prob)
    }


    #P.grmSt(c(1,1/3,-1,0,1), matrix(c(1,1,0,0,-1,-1), byrow=T, ncol=2), prev=c(1,2,3), ncat=4)

    #par=c(1,1/3,-1,0,1)
    #Theta=matrix(c(1,1,0,0,-1,-1), byrow=T, ncol=2)
    #prev=c(1,2,3)
    #ncat=4

    R = matrix(0, nrow=n.subj, ncol=n.item)
    R2 = matrix(0, nrow=n.subj, ncol=n.item)

    i.item = 1
    u <- runif(n.subj, 0, 1)
    U <- matrix(u, ncol=ncat, nrow=n.subj)
    P <- cbind(1, P.grm(c(a[i.item], d[i.item,]), Theta[,1],  ncat))
    R[, i.item]=apply(U-P <= 0 , 1, function(x) max(which(x)))

    for (i.item in 2:n.item) {
        u <- runif(n.subj, 0, 1)
        U <- matrix(u, ncol=ncat, nrow=n.subj)
        P <- cbind(1, P.grmSt(c(a[i.item], 1, d[i.item,]), Theta, R[, i.item-1], ncat=4))
        # par=c(a[i.item], 1/3, d[i.item,]); Theta=Theta; prev=R[, i.item-1]; ncat=4
        R[, i.item]=apply(U-P <= 0 , 1, function(x) max(which(x)))
        R2[, i.item]=apply(U-P <=0, 1, sum)  # in this case, all categories should be number 1,2,3,4,...
    }

    return(list(R,a,d))
}


# Data(R) Generation
n.subj = 300; n.item = 10; ncat=4; # we need this for following scripts
#Theta = matrix(rnorm(n.subj*2), ncol=2)
#Theta = matrix(c(rnorm(n.subj), runif(n.subj, -4, 4)), ncol=2)
# I thought runif( , -4, 4) would discriminate response all the more.
# But latent variable has normal prior.
Theta = matrix(c(rnorm(n.subj), rnorm(n.subj)), ncol=2)

list.RespData <- createR(n.subj, n.item, ncat=ncat, Theta)

R <- list.RespData[[1]]
a <- list.RespData[[2]]
d <- list.RespData[[3]]

# Estimation using a.st=1
apply(R, 2, function(x) length(unique(x)))

R.repeat <- t(apply(R, 1, diff) == 0)*1
#R.repeat1 <- (runif(nrow(R), 0, 1) > 1/2) *1

R2 <- cbind(0, R.repeat*ncat)  + R
apply(R2, 2, function(x) length(unique(x)))

#save(file="testing_grmSt.RData", R, R2, a, d, Theta)
#load(file="testing_grmSt.RData")

#save(file="testing_grmSt2.RData", R, R2, a, d, Theta)  # subj=1000, item=50
#load(file="testing_grmSt2.RData")

#save(file="testing_grmSt3.RData", R, R2, a, d, Theta)  # subj=300, item=10, Theta = matrix(c(rnorm(n.subj), runif(n.subj, -4, 4)), ncol=2)
#load(file="testing_grmSt3.RData")

#FOR CURRENT R
#save(file="testing_grmSt4.RData", R, R2, a, d, Theta)  # subj=300, item=30, Theta = matrix(c(rnorm(n.subj), runif(n.subj, -4, 4)), ncol=2)
#load(file="testing_grmSt4.RData")

model <- paste("F1 = 1-10
               F2 = 2-10
               ", sep="")

model <- paste("F1 = 1-30
               F2 = 2-30
               ", sep="")

n.item=ncol(R2)

#mirt(R, model=mirt.model(s.mirt.model.a), itemtype=rep('egrm10a', n.item), pars='values')
mod2 <- mirt(R2, model=mirt.model(model), itemtype=c('graded', rep('grmSt', n.item-1)), pars='values')
system.time(mod2e5 <- mirt(R2,
                         TOL=1e-5,
                         model=mirt.model(model), itemtype=c('graded', rep('grmSt', n.item-1)), technical=list(NCYCLES=1000)))
system.time(mod2e6 <- mirt(R2,
                         TOL=1e-6,
                         model=mirt.model(model), itemtype=c('graded', rep('grmSt', n.item-1)), technical=list(NCYCLES=2000)))
#mod2e5 <- mod2


# transforming normal to jeffrey's prior (beta(1/2,1/2)))
hist((tanh(2*rnorm(10000))+1)/2, freq=F)
hist(norm2p(rnorm(10000)), freq=F)
curve(dbeta(x,1/2,1/2), xlim=c(0,1), add=T)

norm2p = function(x) {
    (tanh(2*x)+1)/2
}


#================= TESTING ESTIMATES
mod = mod2e5
mod = mod2
mod = mod2MHRM

par(mfcol=c(2,1))
hist(tanh(Theta[,2]/2))
hist(tanh(fscores(mod, full.scores=T)[,2]/2), xlim=c(-1,1))

hist(Theta[,2])
hist(fscores(mod, full.scores=T)[,2])


a.est=c()
for (i in 1:n.item)
    a.est[i] = coef(mod)[[i]]["par","a1"]
plot(a.est,a); abline(a=0, b=1)

plot(a.est,a, pch="")
text(a.est,a)
cor(a.est,a)

d.est = list()
for (i in 1:n.item)
    d.est[[i]] = coef(mod)[[i]]["par", c(-1, -2), drop=T]

d.est[[1]] = sort(d.est[[1]])
names(d.est[[1]])=c("d1","d2","d3")

source("R/ParamDifficulty.r")
d.new = paramDiff(d, R)

d.est.unlist <- unlist(d.est); d.new.unlist <- unlist(d.new)
cor(unlist(d.est), unlist(d.new))

d1.est <- d.est.unlist[names(d.est.unlist)=="d1"]
d2.est <- d.est.unlist[names(d.est.unlist)=="d2"]
d3.est <- d.est.unlist[names(d.est.unlist)=="d3"]

d1.new <- d.new.unlist[names(d.new.unlist)=="d1"]
d2.new <- d.new.unlist[names(d.new.unlist)=="d2"]
d3.new <- d.new.unlist[names(d.new.unlist)=="d3"]

cor(d1.est, d1.new)
cor(d2.est, d2.new)
cor(d3.est, d3.new)

plot(unlist(d.est), unlist(d.new))
plot(d1.est, d1.new, pch="")
text(d1.est, d1.new)
plot(d2.est, d2.new, pch="")
text(d2.est, d2.new)
plot(d3.est, d3.new, pch="")
text(d3.est, d3.new)

cor(fscores(mod, full.scores=T)[,1], Theta[,1])
cor(fscores(mod2e5, full.scores=T)[,1], Theta[,1]);
cor(fscores(mod2e6, full.scores=T)[,1], Theta[,1]);
modgrm <- mirt(R, 1)
cor(fscores(modgrm, full.scores=T)[,1], Theta[,1]);

plot(fscores(mod, full.scores=T)[,1], Theta[,1])
plot(fscores(modgrm, full.scores=T)[,1], Theta[,1])
cor(fscores(mod, full.scores=T)[,2], Theta[,2])
plot(fscores(mod, full.scores=T)[,2], Theta[,2])
#plot(abs(fscores(mod2, full.scores=T)[,2]), abs(Theta[,2]))

#for n.subj=300, n.item=10, th1, cor=0.72 -> 0.76
#for n.subj=300, n.item=30, th1, cor=0.63 -> 0.82

modgrm <- mirt(R,1)
cor(fscores(modgrm, full.scores=T)[,1], Theta[,1])
plot(fscores(mod, full.scores=T)[,1], Theta[,1])
plot(fscores(modgrm, full.scores=T)[,1], Theta[,1]); points(fscores(mod, full.scores=T)[,1], Theta[,1], col="red")

str(mod@Data$fulldata)
mod@Data$Freq
mod@Data$tabdata

pp1 <- which(apply(R2,1,function(x) all(x==c(2,6,6,6,6,6,6,6,6,6))))
fscores(mod, full.scores=T)[pp1,]
Theta[pp1,]
apply(Theta[pp1,], 2, mean)

pp2 <- which(fscores(mod, full.scores=T)[,2] < -2.106 & fscores(mod, full.scores=T)[,2] > -2.107)
Theta[pp2,]
apply(Theta[pp2,], 2, mean)

str(R[2,])
str(c(2,6,6,6,6,6,6,6,6,6))
