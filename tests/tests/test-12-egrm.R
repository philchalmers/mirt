# Testing newly built package mirt with egrm item
# For simulation and estimation

# createR : BUG checked through R, R2
# create response with n.subj, n.item, ncat, and 1-dim latent var and xi, both from N(0,1)
createR=function(n.subj, n.item, ncat, Theta) {

    #

    a = runif(n.item, 1, 2)
    d = matrix(runif(n.item*(ncat-1), -3, 3), ncol=ncat-1)
    d = t(apply(d, 1, sort))

    #P2.egrm for creating response data R
    P2.egrm <- function(par, Theta, ncat) {
        th1 = Theta[,1]; xi1 = Theta[,2];
        a = par[1]
        d = par[2:ncat]
        d.mean=mean(d);
        D.star = matrix(exp(Theta[,2]), nrow=nrow(Theta), ncol=ncat-1) *
            matrix((d - d.mean) + d.mean, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
        TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
        A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
        P = 1/(1+exp(-1*(A*(TH1-D.star))))

        return(P)
    }

    # P.egrm for creating Items
    P.egrm <- function(par, Theta, ncat=4) {
        th1 = Theta[,1]; xi1 = Theta[,2];
        a = par[1]
        d = par[2:ncat]
        d.mean=mean(d);
        D.star = matrix(exp(Theta[,2]), nrow=nrow(Theta), ncol=ncat-1) *
            matrix((d - d.mean) + d.mean, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
        TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
        A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
        P = 1/(1+exp(-1*(A*(TH1-D.star))))
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
        P <- cbind(1, P2.egrm(c(a[i.item], d[i.item,]), Theta, ncat))
        R[, i.item]=apply(U-P <= 0 , 1, function(x) max(which(x)))
        R2[, i.item]=apply(U-P <=0, 1, sum)  # in this case, all categories should be number 1,2,3,4,...
    }

    return(R)
}

# Data(R) Generation
n.subj = 100; n.item = 10; ncat=4; # we need this for following scripts
Theta = matrix(rnorm(n.subj*2), ncol=2)
R <- createR(n.subj = n.subj, n.item = n.item, ncat=ncat, Theta=Theta)

#save(R, file= "test_egrm.RData")
load("test_egrm.RData")

#USING customItems EsTIMATE egrm model
name <- "c.egrm"
par <- c(a=1, d1=-1, d2=0, d3=1)
est <- c(T, T, T, T)
# P.egrm for creating Items
P.egrm <- function(par, Theta, ncat) {
    th1 = Theta[,1]; xi1 = Theta[,2];
    a = par[1]
    d = par[2:ncat]
    d.mean=mean(d);
    D.star = matrix(exp(Theta[,2]), nrow=nrow(Theta), ncol=ncat-1) *
        matrix((d - d.mean) + d.mean, nrow=nrow(Theta), ncol=ncat-1, byrow=T)
    TH1 = matrix(th1, nrow=nrow(Theta), ncol=ncat-1)
    A = matrix(a, nrow=nrow(Theta), ncol=ncat-1)
    P = 1/(1+exp(-1*(A*(TH1-D.star))))
    P.star=cbind(1, P)-cbind(P, 0)

    # Is this correct or justifiable?
    P.star <- ifelse(P.star < 1e-20, 1e-20, P.star)
    P.star <- ifelse(P.star > (1 - 1e-20), (1 - 1e-20), P.star)
    #

    return(P.star)
}
item.egrm <- createItem(name, par=par, est=est, P=P.egrm)

#mirt(R, 2, rep('egrm', n.item), customItems=list(egrm=item.egrm), pars='values')
system.time(mod1 <- mirt(R, 2, rep('c.egrm', n.item), customItems=list(c.egrm=item.egrm)))

#Are the estimates accurate?
cor(fscores(mod1, full.scores=T)[,1], Theta[,1])
cor(fscores(mod1, full.scores=T)[,2], Theta[,2])

#Checking built-in itemtype "egrm10" : which 1-dimensional factor model, with slope for factor 1 and xi
system.time(mod2 <- mirt(R, 2, itemtype=rep("egrm10",n.item)))

#Are the estimates accurate?
cor(fscores(mod2, full.scores=T)[,1], Theta[,1])
cor(fscores(mod2, full.scores=T)[,2], Theta[,2])



# Estimating with PRIORS
s.mirt.model1 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1), (1-10, d3, norm, 0, 1)"

s.mirt.model2 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a1, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1), (1-10, d3, norm, 0, 1)"

s.mirt.model3 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a1, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1), (2-10, d3, norm, 0, 1)"

s.mirt.model4 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a1, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1)"

system.time(mod1.prior <- mirt(R, mirt.model(s.mirt.model1), rep('c.egrm', n.item), customItems=list(c.egrm=item.egrm)))
# Works fine, but item1 has parameter "d3"
coef(mod1.prior)

#Are the estimates accurate?
cor(fscores(mod1.prior, full.scores=T)[,1], Theta[,1])
cor(fscores(mod1.prior, full.scores=T)[,2], Theta[,2])
#And estimates are poor

s.mirt.model11 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1), (1-10, d3, norm, 0, 1)
"

system.time(mod21.prior <- mirt(R, mirt.model(s.mirt.model11), rep('c.egrm', n.item), customItems=list(c.egrm=item.egrm)))
cor(fscores(mod21.prior, full.scores=T)[,1], Theta[,1])
cor(fscores(mod21.prior, full.scores=T)[,2], Theta[,2])
# Estimates distorted?

# Is it because of d3 parameter of item1?
s.mirt.model12 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1), (2-10, d3, norm, 0, 1)"

system.time(mod22.prior <- mirt(R, mirt.model(s.mirt.model12), rep('c.egrm', n.item), customItems=list(c.egrm=item.egrm)))
cor(fscores(mod22.prior, full.scores=T)[,1], Theta[,1])
cor(fscores(mod22.prior, full.scores=T)[,2], Theta[,2])

# Or maybe Starint point?
s.mirt.model13 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1), (1-10, d3, norm, 0, 1)
START = (1-10, a, 1.0), (1-10, d1, -1.0), (1-10, d2, 0), (1-10, d3, 1.0)"

mirt(R, mirt.model(s.mirt.model13), rep('c.egrm', n.item), customItems=list(c.egrm=item.egrm), pars='values')
temp <- mirt(R, mirt.model(s.mirt.model13), rep('c.egrm', n.item), customItems=list(c.egrm=item.egrm), technical=list(NCYCLES=1))
coef(temp)
# If you look at the result of coef(temp), START doesn't seem to work.

# THE LAST RESORT
parprior=list()
parnum=1:40
itemtype=rep("norm", 40); par1=rep(0, 40); par2=rep(1, 40);
itemtype[seq(1,37,4)]="lnorm";
for (i.parnum in parnum) {
    parprior=c(parprior, list(c(i.parnum, itemtype[i.parnum],par1[i.parnum],par2[i.parnum])))
}

system.time(mod23.prior <- mirt(R, 2, rep('c.egrm', n.item), customItems=list(c.egrm=item.egrm), parprior=parprior))
mirt(R, 2, rep('c.egrm', n.item), customItems=list(c.egrm=item.egrm), pars="values")

#Are the estimates accurate?
cor(fscores(mod23.prior, full.scores=T)[,1], Theta[,1])
cor(fscores(mod23.prior, full.scores=T)[,2], Theta[,2])

#Not so good, .16, .10
# maybe parameter d3 of item1 should be deleted?




system.time(mod2.prior <- mirt(R, 2, rep('egrm1', n.item)))
# ERROR : it seemes that every factor has to have a slope parameter

system.time(mod2.prior <- mirt(R, mirt.model(s.mirt.model2), rep('egrm10', n.item)))
# ERROR : yes, item 1 has only 3 categories and has no d3 parameter

system.time(mod3.prior <- mirt(R, mirt.model(s.mirt.model3), rep('egrm10', n.item)))
cor(fscores(mod3.prior, full.scores=T)[,1], Theta[,1])
cor(fscores(mod3.prior, full.scores=T)[,2], Theta[,2])


# but how about this
s.mirt.model2 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a1, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1), (2-10, d3, norm, 0, 1)"
system.time(mod2.prior <- mirt(R, mirt.model(s.mirt.model2), rep('egrm10', n.item)))


# adding starting point
s.mirt.model22 <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a1, lnorm, 0, 1), (1-10, d1, norm, 0, 1), (1-10, d2, norm, 0, 1), (2-10, d3, norm, 0, 1)
START = (1-10, a1, 1.0), (1-10, d1, -1.0), (1-10, d2, 0), (2-10, d3, 1.0)"
system.time(mod22.prior <- mirt(R, mirt.model(s.mirt.model22), rep('egrm10', n.item)))
# This sometimes works, but other times, ERROR, Error: Parameter 'd3' does not exist for item 1

system.time(mod2.prior <- mirt(R, mirt.model(s.mirt.model4), rep('egrm10', n.item)))
# Model has no d3 parameter for all items, no ERROR, ???????






