# GRM with Straight Responses,
# Simulation Study

logis=function(x) {
    1/(1+exp(-x))
}


prob.st <- function(p.st, a.st=1/3) { 2*abs(1/(1+exp(-a.st*p.st))-1/2) }

createR=function(n.subj, n.item, ncat, Theta) {

    #n.subj=10; n.item=5; ncat=5; Theta=matrix(rnorm(n.subj*2), n.subj, 2 )
    #covariate=matrix(sample(ncat, n.subj, replace=T), n.subj)

    a = runif(n.item, 1, 2)
    d = matrix(runif(n.item*(ncat-1), -3, 3), ncol=ncat-1)
    d = t(apply(d, 1, sort))

    #P2.egrm for creating response data R
    P2.grmSt <- function(par, Theta, covariate, ncat) {
        th1 = Theta[,1]; pr = prob.st(Theta[,2]);
        a = par[1]
        d = par[2:ncat]
        #d.mean=mean(d);
        #D.star = matrix(exp(Theta[,2]), nrow=nrow(Theta), ncol=ncat-1) *
        #    matrix((d - d.mean) + d.mean, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
        D = matrix(d, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
        TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
        A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
        P_ = 1/(1+exp(-1*(A*(TH1-D))))


        if (!is.null(covariate)) {
            # covariate be n.subj*1 matrix
            PR = matrix(pr, nrow=n.subj, ncol=ncat)
            PRcat = t(apply(covariate, 1, function(x) {1:ncat == x})) * PR

            PR2 = PRcat

            for (i in ncat:1) {
                PR2[,i] = apply(PRcat[,i:ncat, drop=F], 1, sum)
            }
            PR2 = PR2[,-1]
        } else {
            PR = matrix(0, nrow=n.subj, ncol=ncat)
            PR2 = matrix(0, nrow=n.subj, ncol=ncat-1)
        }


        return((1-PR)[,-1]*P_+PR2)
    }

    # P.egrm for creating Items
    P.egrm <- function(par, Theta, ncat) {
        th1 = Theta[,1]; pr = logis(Theta[,2]);
        a = par[1]
        d = par[2:ncat]
        #d.mean=mean(d);
        #D.star = matrix(exp(Theta[,2]), nrow=nrow(Theta), ncol=ncat-1) *
        #    matrix((d - d.mean) + d.mean, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
        D = matrix(d, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
        TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
        A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
        P_ = 1/(1+exp(-1*(A*(TH1-D))))


        TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
        A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
        P = 1/(1+exp(-1*(A*(TH1-D.star))))

        if (!is.null(covariate)) {
            # covariate be n.subj*1 matrix
            PR = matrix(pr, nrow=n.subj, ncol=ncat)
            PRcat = t(apply(covariate, 1, function(x) {1:ncat == x})) * PR

            PR2 = PRcat

            for (i in ncat:1) {
                PR2[,i] = apply(PRcat[,i:ncat, drop=F], 1, sum)
            }
            PR2 = PR2[,-1]
        } else {
            PR = matrix(0, nrow=n.subj, ncol=ncat)
            PR2 = matrix(0, nrow=n.subj, ncol=ncat-1)
        }










        P.star=cbind(1, P)-cbind(P, 0)

        # Is this correct or justifiable?
        P.star <- ifelse(P.star < 1e-20, 1e-20, P.star)
        P.star <- ifelse(P.star > (1 - 1e-20), (1 - 1e-20), P.star)
        #

        return(P.star)
    }

    R = matrix(0, nrow=n.subj, ncol=n.item)
    R2 = matrix(0, nrow=n.subj, ncol=n.item)
    for (i.item in 1:n.item) {
        u <- runif(n.subj, 0, 1)
        U <- matrix(u, ncol=ncat, nrow=n.subj)
        if (i.item>1) {
            P <- cbind(1, P2.grmSt(c(a[i.item], d[i.item,]), Theta, covariate=R[, i.item-1, drop=F], ncat))
        } else {
            P <- cbind(1, P2.grmSt(c(a[i.item], d[i.item,]), Theta, covariate=NULL, ncat)) }

        R[, i.item]=apply(U-P <= 0 , 1, function(x) max(which(x)))
        R2[, i.item]=apply(U-P <=0, 1, sum)  # in this case, all categories should be number 1,2,3,4,...
    }

    return(list(R,a,d))
}

# Data(R) Generation
n.subj = 300; n.item = 10; ncat=4; # we need this for following scripts
Theta = matrix(rnorm(n.subj*2), ncol=2)
#Theta[,2] = Theta[,2]*.1
list.RespData <- createR(n.subj = n.subj, n.item = n.item, ncat=ncat, Theta=Theta)

R <- list.RespData[[1]]
a <- list.RespData[[2]]
d <- list.RespData[[3]]

#save(list.RespData, n.subj, n.item, ncat, file= "test_egrm.RData")
#load("test_egrm.RData")

apply(R, 2, function(x) length(unique(x)))

R.repeat <- t(apply(R, 1, diff) == 0)*1
R.repeat1 <- (runif(nrow(R), 0, 1) > 1/2) *1

R2 <- cbind(R.repeat1, R.repeat)  * ncat + R
apply(R2, 2, function(x) length(unique(x)))

#save(file="testing_gradedSt.RData", R, R2, a, d, Theta)
#load(file="testing_gradedSt.RData")

#Checking built-in itemtype "egrm10" : which 1-dimensional factor model, with slope for factor 1 and xi
system.time(mod2 <- mirt(R, 2, itemtype=rep("gradedSt",n.item)))
system.time(mod2 <- mirt(R2, 2, itemtype=rep("gradedSt",n.item)))
system.time(mod1 <- mirt(R, 1, itemtype=rep("graded",n.item)))



#============================
# Using Model; first item doesn't need second factor

model <- paste("F1 = 1-10
                        F2 = 2-10
                        ", sep="")

#mirt(R, model=mirt.model(s.mirt.model.a), itemtype=rep('egrm10a', n.item), pars='values')
mirt(R, model=mirt.model(model), itemtype=c('graded', rep('gradedSt', n.item-1)), pars='values')
system.time(mod2a <- mirt(R, model=mirt.model(model), itemtype=c('graded', rep('gradedSt', n.item-1))))

system.time(mod2a <- mirt(R, model=mirt.model(model), itemtype=c('graded', rep('gradedSt', n.item-1))), )

# sd of a.xi estimates
coef(mod3a)[[1]]["par","a.xi"]
sd(coef(mod3a)[[1]]["par","a.xi"]*fscores(mod3a, full.scores=T)[,2])
sd(fscores(mod3a, full.scores=T)[,2])
sd(Theta[,2])

#=============================

coef(mod1)


mirt(R, 2, itemtype=rep("gradedSt",n.item), pars='values')
mirt(R2, 2, itemtype=rep("gradedSt",n.item), pars='values')




# How accurate are the estimates?

cor(fscores(mod2, full.scores=T)[,1], Theta[,1])
cor(fscores(mod2, full.scores=T)[,2], Theta[,2])

mod = mod2
a.est=c()
for (i in 1:n.item)
    a.est[i] = coef(mod)[[i]]["par","a1"]
plot(a.est,a); abline(a=0, b=1)
cor(a.est,a)

d.est = list()
for (i in 1:n.item)
    d.est[[i]] = coef(mod)[[i]]["par", c(-1, -2), drop=T]

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
plot(d1.est, d1.new, pch=1:length(d1.est))
plot(d2.est, d2.new, pch=1:length(d1.est))
plot(d3.est, d3.new, pch=1:length(d1.est))

cor(fscores(mod2, full.scores=T)[,1], Theta[,1])
plot(fscores(mod2, full.scores=T)[,1], Theta[,1])
cor(fscores(mod2, full.scores=T)[,2], Theta[,2])
plot(abs(fscores(mod2, full.scores=T)[,2]), abs(Theta[,2]))

plot(prob.st(fscores(mod2, full.scores=T)[,2]), prob.st(Theta[,2]))
