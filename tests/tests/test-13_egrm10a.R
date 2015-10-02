# testing itemtype="egrm10a"


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

    return(list(R,a,d))
}

# Data(R) Generation
n.subj = 100; n.item = 10; ncat=4; # we need this for following scripts
Theta = matrix(rnorm(n.subj*2), ncol=2)
Theta[,2] = Theta[,2]*.3
list.RespData <- createR(n.subj = n.subj, n.item = n.item, ncat=ncat, Theta=Theta)

R <- list.RespData[[1]]
a <- list.RespData[[2]]
d <- list.RespData[[3]]

#save(list.RespData, n.subj, n.item, ncat, file= "test_egrm.RData")
#load("test_egrm.RData")

apply(R, 2, function(x) length(unique(x)))

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
CONSTRAIN = (1-10, a.xi)
PRIOR = (1-10, a.xi, lnorm, 0, 0.3)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a.xi, lnorm, 0, 0.3)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a1, norm, 0, 1)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a1, lnorm, 0, 1)
" # I dont know why setting a1 PRIOS (lnorm, 0, 1) and a1* doesnt work and
# setting a1 PRIORS (norm, 0, 1) and exp(a1*...)

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a.xi, norm, 0.3, 0.1)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a.xi, norm, 0.3, 0.1)
CONSTRAIN = (1-10, a.xi)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a.xi, lnorm, -3, 0.1)
CONSTRAIN = (1-10, a.xi)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a.xi, norm, 0.184, 0.5)
CONSTRAIN = (1-10, a.xi)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
CONSTRAIN = (1-10, a.xi)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a.xi, norm, 0.184, 1)
CONSTRAIN = (1-10, a.xi)
"

s.mirt.model.a <- "F1 = 1-10
F2 = 1-10
PRIOR = (1-10, a.xi, norm, 0.299, 1)
CONSTRAIN = (1-10, a.xi)
"

#mirt(R, model=mirt.model(s.mirt.model.a), itemtype=rep('egrm10a', n.item), pars='values')
system.time(mod3a <- mirt(R, model=mirt.model(s.mirt.model.a), itemtype=rep('egrm10a', n.item)))

# sd of a.xi estimates
coef(mod)[[1]]["par","a.xi"]
sd(coef(mod)[[1]]["par","a.xi"]*fscores(mod, full.scores=T)[,2])
sd(fscores(mod, full.scores=T)[,2])
sd(Theta[,2])


# How accurate are the estimates?

mod = mod3a
a.est=c()
for (i in 1:n.item)
    a.est[i] = coef(mod)[[i]]["par","a1"]
plot(a.est,a); abline(a=0, b=1)
cor(a.est,a)

d.est = list()
for (i in 1:n.item)
    d.est[[i]] = coef(mod)[[i]]["par", c(-1,-2), drop=T]

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

#Are the estimates accurate?
cor(fscores(mod, full.scores=T)[,1], Theta[,1])
cor(fscores(mod, full.scores=T)[,2], Theta[,2])

plot(fscores(mod, full.scores=T)[,1], Theta[,1])
plot(fscores(mod, full.scores=T)[,2], Theta[,2])

cor(exp(coef(mod)[[1]]["par","a.xi"]*fscores(mod, full.scores=T)[,2]),
    Theta[,2])


# sd of a.xi estimates
sd(coef(mod)[[1]]["par","a.xi"]*fscores(mod, full.scores=T)[,2])
sd(fscores(mod, full.scores=T)[,2])
sd(Theta[,2])


cor(coef(mod)[[1]]["par","a.xi"]*fscores(mod, full.scores=T)[,2],
    Theta[,2])



cor(fscores(mod, full.scores=T)[,1], fscores(mod1, full.scores=T)[,1])
cor(fscores(mod, full.scores=T)[,2], fscores(mod1, full.scores=T)[,2])
