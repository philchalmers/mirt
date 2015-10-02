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
n.subj = 500; n.item = 20; ncat=4; # we need this for following scripts
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

model <- paste("F1 = 1-15
               F2 = 1-15
               F3 = 2-15
               ", sep="")

model <- paste("F1 = 1-20
               F2 = 2-20
               ", sep="")

model <- paste("F1 = 1-20
               F2 = 1-20
               F3 = 2-20
               ", sep="")

model <- paste("F1 = 1-20
               F2 = 1-20
               ", sep="")

n.item=ncol(R2)

mirt(R2, model=mirt.model(model), itemtype=c('egrm10', rep('egrmSt2', n.item-1)), pars='values')

system.time(mirt.fit1e7 <- mirt(R2, model=mirt.model(model), itemtype=c('egrm10', rep('egrmSt2', n.item-1)),
                           TOL=1e-7, technical=list(NCYCLES=2000)))

system.time(mirt.fit <- mirt(R2, model=mirt.model(model), itemtype=c('egrm10', rep('egrmSt2', n.item-1)),
                          technical=list(NCYCLES=2000)))

system.time(mirt.fitegrm <- mirt(R, model=mirt.model(model), itemtype=c('egrm10', rep('egrm10', n.item-1)),
                                 technical=list(NCYCLES=2000))) # Check Model, Check input

system.time(mirt.fitgrmSt <- mirt(R2, model=mirt.model(model), itemtype=c('graded', rep('grmSt2', n.item-1)),
                                 technical=list(NCYCLES=2000))) # Check Model, Check input


system.time(mirt.fitMHRM <- mirt(R2, model=mirt.model(model), itemtype=c('egrm10', rep('egrmSt2', n.item-1)),
                              method="MHRM"))
# serious degeneration from EM!!





mirt.fit <- mirt(R2, model=mirt.model(model), itemtype=c('egrm10', rep('egrmSt2', n.item-1)),
                 TOL=1e-7,
                 optimizer = "nlminb")

# IT must be about optim!!
# or should have set lbound, ubound?!

######## CURRENT POSITION









system.time(mod2e5 <- mirt(R2,
                           TOL=1e-5,
                           model=mirt.model(model), itemtype=c('egrm', rep('egrmSt2', n.item-1)), technical=list(NCYCLES=1000)))
system.time(mod2e6 <- mirt(R2,
                           TOL=1e-6,
                           model=mirt.model(model), itemtype=c('graded', rep('grmSt2', n.item-1)), technical=list(NCYCLES=2000)))
system.time(mod2e5MHRM <- mirt(R2,
                               TOL=1e-5,
                               model=mirt.model(model), itemtype=c('graded', rep('grmSt2', n.item-1)),
                               method="MHRM", technical=list(NCYCLES=1000)))

#mod2e5 <- mod2

# transforming normal to jeffrey's prior (beta(1/2,1/2)))
hist((tanh(2*rnorm(10000))+1)/2, freq=F)
hist(norm2p(rnorm(10000)), freq=F)
curve(dbeta(x,1/2,1/2), xlim=c(0,1), add=T)

norm2p = function(x) {
    (tanh(2*x)+1)/2
}


#================= TESTING ESTIMATES
fit = mirt.fit
fit = mirt.fitMHRM
fit = mirt.fit1e7



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
fitgrm <- mirt(R, 1)
cor(fscores(fitgrm, full.scores=T)[,1], Theta[,1]);

plot(fscores(fit, full.scores=T)[,1], Theta[,1])
fs1 <- fscores(fit, full.scores=T, full.scores.SE=T)
fs2 <- fscores(fit, full.scores=T, full.scores.SE=T, rotate="none")
all.equal(fs1, fs2) # TRUE

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


fs <- fscores(fitgrm, full.scores=T)[,"F1"]
sdev <- fscores(fitgrm, full.scores=T, full.scores.SE=T)[,"SE_F1"]
th <- Theta[,1]

cor(fs, th);
plot(fs ~ th, pch=19, xlim=c(-3,3), ylim=c(-3,3)); abline(a=0, b=1, col="red", lty="dotted")
arrows(th, fs-sdev, th, fs+sdev, length=0.05, angle=90, code=3)


cor(fscores(fit, full.scores=T)[,2], Theta[,2]) # xi's
plot(fscores(fit, full.scores=T)[,2], Theta[,2])
cor(fscores(fit, full.scores=T)[,3], Theta[,3]) # st's
plot(fscores(fit, full.scores=T)[,3], Theta[,3])
cor(norm2p(fscores(fit, full.scores=T)[,3]), norm2p(Theta[,3])) # st's
plot(norm2p(fscores(fit, full.scores=T)[,3]), norm2p(Theta[,3]))


#plot(abs(fscores(fit2, full.scores=T)[,2]), abs(Theta[,2]))

#for n.subj=300, n.item=10, th1, cor=0.81 -> 0.74 (using grm for grmSt2 Data)
#for n.subj=300, n.item=30, th1,
#for n.subj=300, n.item=10, th1, cor=0.89 -> 0.87 (using grmSt for grm Data)
#for n.subj=1000, n.item=20, th1, cor=0.955 -> 0.954 (using grmSt for grm Data)
#                                     0.947 -> 0.940
fitgrm <- mirt(R,1)
cor(fscores(fitgrm, full.scores=T)[,1], Theta[,1])
plot(fscores(fit, full.scores=T)[,1], Theta[,1])
plot(fscores(fitgrm, full.scores=T)[,1], Theta[,1]); points(fscores(fit, full.scores=T)[,1], Theta[,1], col="red")

t1 <- Theta[,1]
t2 <- fscores(fitgrm, full.scores=T)[,1]
t3 <- fscores(fit, full.scores=T)[,1]

str(fit@Data$fulldata)
fit@Data$Freq
fit@Data$tabdata[15,]

pp1 <- which(apply(R2,1,function(x) all(x==c(2,6,6,6,6,6,6,6,6,6))))
fscores(fit, full.scores=T)[pp1,]
Theta[pp1,]
apply(Theta[pp1,], 2, mean)

#pp2 <- which(fscores(fit, full.scores=T)[,2] < -2.106 & fscores(fit, full.scores=T)[,2] > -2.107)
pp2 <- which(fscores(fit, full.scores=T)[,2] < 1.11 & fscores(fit, full.scores=T)[,2] > 1.10)
Theta[pp2,]
apply(Theta[pp2,], 2, mean)

str(R[2,])
str(c(2,6,6,6,6,6,6,6,6,6))
