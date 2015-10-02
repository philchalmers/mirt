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
        th1 = Theta[,1];
        pr = logis(Theta[,2])
        #pr = prob.st(Theta[,2]);
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

# Data Generation: Toy Example
RR <- createR(7, 5, 4, matrix(c(-3,-2,-1,0,1,2,3,0,0,0,10,0,0,0), ncol=2))
R <- RR[[1]]


# Data(R) Generation
n.subj = 300; n.item = 10; ncat=4; # we need this for following scripts
#Theta = matrix(rnorm(n.subj*2), ncol=2)
Theta = matrix(c(rnorm(n.subj), runif(n.subj, -4, 4)), ncol=2)
#Theta[,2] = Theta[,2]
hist(logis(Theta[,2]))
list.RespData <- createR(n.subj = n.subj, n.item = n.item, ncat=ncat, Theta=Theta)

R <- list.RespData[[1]]
a <- list.RespData[[2]]
d <- list.RespData[[3]]

#save(list.RespData, n.subj, n.item, ncat, file= "test_egrm.RData")
#load("test_egrm.RData")

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

model <- paste("F1 = 1-10
                        F2 = 2-10
               ", sep="")

model <- paste("F1 = 1-50
               F2 = 2-50
               ", sep="")

n.item=ncol(R2)

#mirt(R, model=mirt.model(s.mirt.model.a), itemtype=rep('egrm10a', n.item), pars='values')
mod2 <- mirt(R2, model=mirt.model(model), itemtype=c('graded', rep('grmSt', n.item-1)), pars='values')
system.time(mod2 <- mirt(R2, model=mirt.model(model), itemtype=c('graded', rep('grmSt', n.item-1))))

system.time(mod2MHRM <- mirt(R2, model=mirt.model(model), itemtype=c('graded', rep('grmSt', n.item-1)), method="MHRM"))
system.time(mod2MHRM <- mirt(R2, model=mirt.model(model), itemtype=c('graded', rep('grmSt', n.item-1)), method="QMCEM"))

mirt(R, model=mirt.model(model), itemtype=rep('graded', n.item), pars='values')
mirt(R, model=mirt.model(model), itemtype=rep('grmSt', n.item), pars='values')

system.time(mod2a <- mirt(R, model=mirt.model(model), itemtype=c('graded', rep('grmSt', n.item-1))))


th2 <- rnorm(1000)/3
hist(prob.st(th2))

#================= TESTING ESTIMATES
mod = mod2
mod = mod2MHRM

hist(prob.st(fscores(mod, full.scores=T)[,2]))
hist(logis(1/3*fscores(mod, full.scores=T)[,2]))
hist(fscores(mod, full.scores=T)[,2])

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
plot(d1.est, d1.new, pch="")
text(d1.est, d1.new)
plot(d2.est, d2.new, pch="")
text(d2.est, d2.new)
plot(d1.est, d1.new, pch="")
text(d1.est, d1.new)

cor(fscores(mod, full.scores=T)[,1], Theta[,1])
plot(fscores(mod, full.scores=T)[,1], Theta[,1])
cor(fscores(mod, full.scores=T)[,2], Theta[,2])
plot(fscores(mod, full.scores=T)[,2], Theta[,2])
#plot(abs(fscores(mod2, full.scores=T)[,2]), abs(Theta[,2]))

probSt = function(p.st, a.st=1/3) {
  p.st2 <- exp(-a.st*p.st)
  p.repeat <- (p.st2-1)/(p.st2+1)
  p.repeat }

plot(prob.st(fscores(mod, full.scores=T)[,2]), prob.st(Theta[,2]))
plot(logis(1/3*fscores(mod, full.scores=T)[,2]), logis(1/3*Theta[,2]))

plot(tanh(-1*1/3*fscores(mod, full.scores=T)[,2]/2),  logis(Theta[,2]))

cor(probSt(fscores(mod, full.scores=T)[,2]), prob.st(Theta[,2]))
plot(probSt(fscores(mod, full.scores=T)[,2]), prob.st(Theta[,2]))

plot(probSt(fscores(mod, full.scores=T)[,2]), logis(1/3*Theta[,2]))

which(probSt(fscores(mod, full.scores=T)[,2]) < -0.1 & logis(1/3*Theta[,2]) > 0.6)
# Look into the case of 299, 4888888824
fscores(mod, full.scores=T)[299,1]
prob.st(fscores(mod, full.scores=T)[299,2])

R2[299,]; Theta[299,1]

th <- fscores(mod, full.scores=T)[299,]

item3 <- mod@pars[[3]]

pb(mod@pars[[3]], Theta[299,,drop=F])
pb(mod@pars[[3]], matrix(th, ncol=2))
pb(mod@pars[[3]], matrix(c(2.228, 0.1), ncol=2))
pb(mod@pars[[3]], matrix(c(2.228, 0.2), ncol=2))
pb(mod@pars[[3]], matrix(c(2.228, 0.3), ncol=2))
pb(mod@pars[[3]], matrix(c(2.228, 0.4), ncol=2))
pb(mod@pars[[3]], matrix(c(2.228, 0.5), ncol=2))
pb(mod@pars[[3]], matrix(c(2.228, 0.6), ncol=2))
pb(mod@pars[[3]], matrix(c(2.228, 3), ncol=2))
pb(mod@pars[[3]], matrix(c(2.228, 10), ncol=2))

pb=function(x, Theta){
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
        p.repeat <- 2*abs(1/(1+exp(-a.st*p.st))-1/2)
        q.repeat <- 1- p.repeat

        P <- cbind(P.star * q.repeat, P.star * q.repeat + p.repeat)

        P <- ifelse(P < 1e-20, 1e-20, P)
        P <- ifelse(P > (1 - 1e-20), (1 - 1e-20), P)


    } else {

        P <- matrix(0, nrow=nrow(Theta), ncol=ncat*2)}

    return(P)
}
