# From test-15_gradedSt_03.R
# from test-15_gradedSt_04.R
# changing grmSt to egrmSt

norm2p = function(x) {
    (tanh(2*x)+1)/2
}

norm2p.gradedSt_03 = function(x) {
    tanh(-x/2)
}


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

        stopifnot(length(prev) == nrow(Theta));
        stopifnot(length(par) == ncat + 2);

        th1 = Theta[,1];
        xi1 = Theta[,2]; st1 = Theta[,3]

        # Differ from P.egrm
        #ncat = ncat/2
        ncat=ncat

        a = par[1]; a.xi <- par[2]; a.st <- par[3]
        d = par[4:(ncat+2)]
        d.mean = mean(d)
        #
        #if (all(d == sort(d))) {
        D.star = matrix(exp(Theta[,2]), nrow=nrow(Theta), ncol=ncat-1) *
            matrix((d-d.mean)+d.mean, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
        # a.xi included
        TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
        A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
        P = 1/(1+exp(-1*(A*(TH1-D.star))))
        P.star=cbind(1, P)-cbind(P, 0)

        p.st <- matrix(Theta[,3], ncol=ncat, nrow=nrow(Theta))
        ##p.repeat <- 2*abs(1/(1+exp(-a.st*p.st))-1/2)
        #p.st2 <- exp(-a.st*p.st)
        #p.repeat <- (p.st2-1)/(p.st2+1)

        #p.repeat <- tanh(-a.st*p.st/2)
        p.repeat <- norm2p(a.st*p.st)

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
        P <- cbind(1, P.grmSt(c(a[i.item], 1, 1, d[i.item,]), Theta, R[, i.item-1], ncat))
        # par=c(a[i.item], 1/3, d[i.item,]); Theta=Theta; prev=R[, i.item-1]; ncat=4
        R[, i.item]=apply(U-P <= 0 , 1, function(x) max(which(x)))
        R2[, i.item]=apply(U-P <=0, 1, sum)  # in this case, all categories should be number 1,2,3,4,...
    }

    return(list(R,a,d))
}


createRgrm=function(n.subj, n.item, ncat, Theta) {

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
        P <- cbind(1, P.grm(c(a[i.item], d[i.item,]), Theta[,1], ncat=4))
        # par=c(a[i.item], 1/3, d[i.item,]); Theta=Theta; prev=R[, i.item-1]; ncat=4
        R[, i.item]=apply(U-P <= 0 , 1, function(x) max(which(x)))
        R2[, i.item]=apply(U-P <=0, 1, sum)  # in this case, all categories should be number 1,2,3,4,...
    }

    return(list(R,a,d))
}


# Data(R) Generation
n.subj = 200; n.item = 10; ncat=4; # we need this for following scripts
#Theta = matrix(rnorm(n.subj*2), ncol=2)
#Theta = matrix(c(rnorm(n.subj), runif(n.subj, -4, 4)), ncol=2)
# I thought runif( , -4, 4) would discriminate response all the more.
# But latent variable has normal prior.
Theta = matrix(rnorm(n.subj*3), ncol=3)

list.RespData <- createR(n.subj, n.item, ncat=ncat, Theta)

#list.RespData <- createRgrm(n.subj, n.item, ncat=ncat, Theta)

R <- list.RespData[[1]]
a <- list.RespData[[2]]
d <- list.RespData[[3]]

# Estimation using a.st=1
apply(R, 2, function(x) length(unique(x)))

R.repeat <- t(apply(R, 1, diff) == 0)*1
#R.repeat1 <- (runif(nrow(R), 0, 1) > 1/2) *1

R2 <- cbind(0, R.repeat*ncat)  + R
apply(R2, 2, function(x) length(unique(x)))



#FOR CURRENT R, grm data
#save(file="testing_egrmSt.RData", R, R2, a, d, Theta)  # subj=500, item=15, ncat=5
#load(file="testing_egrmSt.RData")

#save(file="testing_egrmSt2.RData", R, R2, a, d, Theta)  # subj=500, item=20, ncat=4
#load(file="testing_egrmSt2.RData")

#save(file="testing_egrmSt3.RData", R, R2, a, d, Theta)  # subj=200, item=10, ncat=4
#load(file="testing_egrmSt3.RData")
#save(file="testing_egrmSt3_RESULT,RData",fit.grm, fit.Egrm, fit.grmSt, fit.EgrmSt )


model.grmSt <- paste("F1 = 1-",n.item,"
               F2 = 2-",n.item,"
               ", sep="")

model.EgrmSt <- paste("F1 = 1-",n.item,"
               F2 = 1-",n.item,"
               F3 = 2-",n.item,"
               ", sep="")

n.item=ncol(R2)

system.time(fit.grm <- mirt(R, 1, technical=list(NCYCLES=2000)))

system.time(fit.grm2 <- mirt(R, 2, technical=list(NCYCLES=5000)))
# default NCYLCES=500,
#    fit.grm2@converge = 1 (CONVERGED)
#    fit.grm2@converge = 0 (NOT CONVERGED)

system.time(fit.Egrm <- mirt(R, 2, itemtype=rep('egrm10', n.item),
                             technical=list(NCYCLES=2000)))

system.time(fit.grmSt <- mirt(R2, model=mirt.model(model.grmSt), itemtype=c('graded', rep('grmSt2', n.item-1)),
                               technical=list(NCYCLES=2000)))

system.time(fit.EgrmSt <- mirt(R2, model=mirt.model(model.EgrmSt), itemtype=c('egrm10', rep('egrmSt2', n.item-1)),
                               technical=list(NCYCLES=4000)))





# checking parameters
#mirt(R2, model=mirt.model(model), itemtype=c('egrm10', rep('egrmSt2', n.item-1)), pars='values')

system.time(fitMHRM.EgrmSt <- mirt(R2, model=mirt.model(model.EgrmSt), itemtype=c('egrm10', rep('egrmSt2', n.item-1)),
                              method="MHRM"))
# serious degeneration from EM!! WHY???


mirt.fit <- mirt(R2, model=mirt.model(model), itemtype=c('egrm10', rep('egrmSt2', n.item-1)),
                 TOL=1e-7,
                 optimizer = "nlminb")
# IT must be about optimizer?!!
# or should have set lbound, ubound?!
# TOL and parameter scale?



# transforming normal to jeffrey's prior (beta(1/2,1/2)))
hist((tanh(2*rnorm(10000))+1)/2, freq=F)
hist(norm2p(rnorm(10000)), freq=F)
curve(dbeta(x,1/2,1/2), xlim=c(0,1), add=T)



#================= TESTING ESTIMATES
# From looking at the specific data and model,
# I can find problems and the way to fix it.
# For instance, for EgrmSt model, d3 parameters tend to get very high... why?
#   Possible solution : NORMAL PRIOR?
#   Do xi's compensate pr's?
fit = fit.Egrm

par(mfcol=c(2,1))
hist(Theta[,2])
hist(fscores(fit, full.scores=T)[,2])

hist(norm2p(Theta[,3]), xlim=c(0,1))
hist(norm2p(fscores(fit, full.scores=T)[,3]), xlim=c(0,1))

a.est=c()
for (i in 1:n.item)
    a.est[i] = coef(fit)[[i]]["par","a1"]
plot(a.est,a); abline(a=0, b=1)

plot(a.est,a, pch="")
text(a.est,a)
cor(a.est,a)

d.est = list()
d.est[[1]] = coef(fit)[[1]]["par", c(-1, -2), drop=T]
for (i in 2:n.item)
    d.est[[i]] = coef(fit)[[i]]["par", c(-1, -2, -3), drop=T]

#d.est[[1]]=-d.est[[1]]
#d.est[[1]] = sort(d.est[[1]])
#names(d.est[[1]])=c("d1","d2","d3")

source("R/ParamDifficulty.r")
d.new = paramDiff(d, R)

d.est.unlist <- unlist(d.est); d.new.unlist <- unlist(d.new)
cor(unlist(d.est), unlist(d.new))

d1.est <- d.est.unlist[names(d.est.unlist)=="d1"]
d2.est <- d.est.unlist[names(d.est.unlist)=="d2"]
d3.est <- d.est.unlist[names(d.est.unlist)=="d3"]
d4.est <- d.est.unlist[names(d.est.unlist)=="d4"]

d1.new <- d.new.unlist[names(d.new.unlist)=="d1"]
d2.new <- d.new.unlist[names(d.new.unlist)=="d2"]
d3.new <- d.new.unlist[names(d.new.unlist)=="d3"]
d4.new <- d.new.unlist[names(d.new.unlist)=="d4"]

cor(d1.est, d1.new)
cor(d2.est, d2.new)
cor(d3.est, d3.new)
cor(d4.est, d4.new)

plot(unlist(d.est), unlist(d.new))
plot(d1.est, d1.new, pch="")
text(d1.est, d1.new)
plot(d2.est, d2.new, pch="")
text(d2.est, d2.new)
plot(d3.est, d3.new, pch="")
text(d3.est, d3.new)
plot(d4.est, d4.new, pch="")
text(d4.est, d4.new)

cor(fscores(fit, full.scores=T)[,1], Theta[,1])
#cor(fscores(fit2e5, full.scores=T)[,1], Theta[,1]);
#cor(fscores(fit2e6, full.scores=T)[,1], Theta[,1]);
fitgrm <- mirt(R, 1);
fitgrm2 <- mirt(R,2, technical=list(NCYCLES=5000));
fitgrm3 <- mirt(R,3, technical=list(NCYCLES=5000));
cor(fscores(fitgrm, full.scores=T)[,1], Theta[,1]);

plot(fscores(fit, full.scores=T)[,1], Theta[,1])
fs1 <- fscores(fit, full.scores=T, full.scores.SE=T)
fs2 <- fscores(fit, full.scores=T, full.scores.SE=T, rotate="none")
all.equal(fs1, fs2) # TRUE

# COMPARING SEVERAL ESTIMATES... EAP SEEMS BEST?!

fsML <- fscores(fit, full.scores=T, full.scores.SE=T, method="ML")
fsWLE <- fscores(fit, full.scores=T, full.scores.SE=T, method="WLE")
fsMAP <- fscores(fit, full.scores=T, full.scores.SE=T, method="MAP")
fsEAP <- fscores(fit, full.scores=T, full.scores.SE=T, method="EAP")
plot(fsEAP, fsWLE); plot(fsEAP[,1], fsWLE[,1])
plot(fsEAP, fsML); plot(fsEAP[,1], fsML[,1])
plot(fsEAP, fsMAP); plot(fsEAP[,1], fsMAP[,1])

cor(Theta[,1], fsEAP[,1])
cor(Theta[,1], fsWLE[,1])
cor(Theta[,1], fsMAP[,1])
cor(Theta[,1], fsML[,1], use="complete.obs")

plot(fs[,1], Theta[,1])
plot(fs1[,1], Theta[,1])
plot(fs1[,1], fs[,1])

fsgrm <- fscores(fitgrm, full.scores=T, full.scores.SE=T)
head(fs)
plot(fs[,1],fs[,3], ylim=c(0,1))
plot(fsgrm[,1],fsgrm[,2], ylim=c(0,1))
plot(fs[,2],fs[,4], ylim=c(0,1))
plot(Theta[,2],fs[,4], ylim=c(0,1))


# ESTIMATES OF F1, F2, F3, ...

fit <- fit.EgrmSt

par(mfcol=c(1,2))

fs <- fscores(fit, full.scores=T)[,"F1"]
sdev <- fscores(fit, full.scores=T, full.scores.SE=T)[,"SE_F1"]
th <- Theta[,1]

cor(fs, th);
plot(fs ~ th, pch=19, xlim=c(-3,3), ylim=c(-3,3)); abline(a=0, b=1, col="red", lty="dotted")
arrows(th, fs-sdev, th, fs+sdev, length=0.05, angle=90, code=3)


fs <- fscores(fit, full.scores=T)[,"F2"]
sdev <- fscores(fit, full.scores=T, full.scores.SE=T)[,"SE_F2"]
th <- Theta[,2]

cor(fs, th);
plot(fs ~ th, pch=19, xlim=c(-3,3), ylim=c(-3,3)); abline(a=0, b=1, col="red", lty="dotted")
arrows(th, fs-sdev, th, fs+sdev, length=0.05, angle=90, code=3)


fs <- fscores(fit, full.scores=T)[,"F3"]
sdev <- fscores(fit, full.scores=T, full.scores.SE=T)[,"SE_F3"]
th <- Theta[,3]

cor(fs, th);
plot(fs ~ th, pch=19, xlim=c(-3,3), ylim=c(-3,3)); abline(a=0, b=1, col="red", lty="dotted")
arrows(th, fs-sdev, th, fs+sdev, length=0.05, angle=90, code=3)



# COMPARISON: ESTIMATES OF F1, grm egrm, grmSt, egrmSt

par(mfcol=c(2,2))

fs <- fscores(fit.grm, full.scores=T)[,"F1"]
sdev <- fscores(fit.grm, full.scores=T, full.scores.SE=T)[,"SE_F1"]
th <- Theta[,1]

c = cor(fs, th);
plot(fs ~ th, pch=19, xlim=c(-3,3), ylim=c(-3,3), main="grm", xlab=paste("th :", round(c,2)));
abline(a=0, b=1, col="red", lty="dotted", lwd=3);
arrows(th, fs-sdev, th, fs+sdev, length=0.05, angle=90, code=3)
fs.correct = sum(fs-sdev < th & fs + sdev > th)/length(th)
sqrt(mean(sdev^2))
fs.correct

fs <- fscores(fit.Egrm, full.scores=T)[,"F1"]
sdev <- fscores(fit.Egrm, full.scores=T, full.scores.SE=T)[,"SE_F1"]
th <- Theta[,1]

c=cor(fs, th);
plot(fs ~ th, pch=19, xlim=c(-3,3), ylim=c(-3,3), main="egrm", xlab=paste("th :", round(c,2)));
abline(a=0, b=1, col="red", lty="dotted", lwd=3);
arrows(th, fs-sdev, th, fs+sdev, length=0.05, angle=90, code=3)
fs.correct = sum(fs-sdev < th & fs + sdev > th)/length(th)
sqrt(mean(sdev^2))
fs.correct


fs <- fscores(fit.grmSt, full.scores=T)[,"F1"]
sdev <- fscores(fit.grmSt, full.scores=T, full.scores.SE=T)[,"SE_F1"]
th <- Theta[,1]

c=cor(fs, th);
plot(fs ~ th, pch=19, xlim=c(-3,3), ylim=c(-3,3), main="grmSt", xlab=paste("th :", round(c,2)));
abline(a=0, b=1, col="red", lty="dotted", lwd=3);
arrows(th, fs-sdev, th, fs+sdev, length=0.05, angle=90, code=3)
fs.correct = sum(fs-sdev < th & fs + sdev > th)/length(th)
sqrt(mean(sdev^2))
fs.correct

fs <- fscores(fit.EgrmSt, full.scores=T)[,"F1"]
sdev <- fscores(fit.EgrmSt, full.scores=T, full.scores.SE=T)[,"SE_F1"]
th <- Theta[,1]

c=cor(fs, th);
plot(fs ~ th, pch=19, xlim=c(-3,3), ylim=c(-3,3), main="egrmSt", xlab=paste("th :", round(c,2)));
abline(a=0, b=1, col="red", lty="dotted", lwd=3);
arrows(th, fs-sdev, th, fs+sdev, length=0.05, angle=90, code=3)
fs.correct = sum(fs-sdev < th & fs + sdev > th)/length(th)
sqrt(mean(sdev^2))
fs.correct

# for normal distribtuion, -1SD/+1SD covers 68%




