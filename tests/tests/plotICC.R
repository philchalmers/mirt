# THE MODEL PLOT

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
Theta[,2] = Theta[,2]*.1
list.RespData <- createR(n.subj = n.subj, n.item = n.item, ncat=ncat, Theta=Theta)

R <- list.RespData[[1]]
a <- list.RespData[[2]]
d <- list.RespData[[3]]

# for item 1

P.egrm <- function(par, Theta, ncat) {
    th1 = Theta[,1, drop=F]; xi1 = Theta[,2, drop=F];
    # drop=F is needed because Theta can be 1-row matrix
    a = par[1]
    d = par[2:ncat]
    d.mean=mean(d);
    D.star = matrix(exp(xi1), nrow=nrow(Theta), ncol=ncat-1) *
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


i.item=5; xi=0
P=matrix(NA, nrow=1, ncol=ncat)
for (i.th in seq(-4,4,0.1)) {
    P=rbind(P, P.egrm(c(a[i.item],d[i.item,]), matrix(c(i.th,xi), ncol=2), ncat=4))
}
P=P[-1,]

plot(seq(-4,4,0.1), P[,1], type="l", ylim=c(0,1))
lines(seq(-4,4,0.1), P[,2])
lines(seq(-4,4,0.1), P[,3])
lines(seq(-4,4,0.1), P[,4])

# Input var into P.egrm
par=c(a[i.item],d[i.item,])
Theta =  matrix(0, ncol=2)
ncat=4





par(mfcol=c(1,2))
#i.item=5; xi=0
P=matrix(NA, nrow=1, ncol=ncat);
for (i.th in seq(-4,4,0.1)) {
    P=rbind(P, P.egrmst(a= a[i.item], d=d[i.item,], g=g, u=u, th=i.th, xi=xi))
}
P=P[-1,]

plot(seq(-4,4,0.1), P[,1], type="l", ylim=c(0,1))
lines(seq(-4,4,0.1), P[,2])
lines(seq(-4,4,0.1), P[,3])
lines(seq(-4,4,0.1), P[,4])



P.egrmst <- function(a, d, g, u, th, xi) {
    d.mean=mean(d)
    d.new = exp(xi)*(d-d.mean)+d.mean
    prob=c(1,g+(u-g)/(1+exp(-(a*(th-d.new)))), 0)
    -diff(prob)}

plotit=function(xmin=-4, xmax=4) {
    P=matrix(NA, nrow=1, ncol=ncat);
    for (i.th in seq(xmin,xmax,0.1)) {
        P=rbind(P, P.egrmst(a= a[i.item], d=d[i.item,], g=g, u=u, th=i.th, xi=xi))
    }
    P=P[-1,]

    plot(seq(xmin,xmax,0.1), P[,1], type="l", ylim=c(0,1), col="red")
    lines(seq(xmin,xmax,0.1), P[,2], col="blue")
    lines(seq(xmin,xmax,0.1), P[,3], col="green")
    lines(seq(xmin,xmax,0.1), P[,4], col="black")
}

plotit2=function(xmin=-4, xmax=4) {
    P=matrix(NA, nrow=1, ncol=ncat);
    for (i.th in seq(xmin,xmax,0.1)) {
        q.straight = 1-sum(p.straight)
        P=rbind(P, p.straight + q.straight * P.egrmst(a= a[i.item], d=d[i.item,], g=0, u=1, th=i.th, xi=xi))
    }
    P=P[-1,]

    plot(seq(xmin,xmax,0.1), P[,1], type="l", ylim=c(0,1), col="red")
    lines(seq(xmin,xmax,0.1), P[,2], col="blue")
    lines(seq(xmin,xmax,0.1), P[,3], col="green")
    lines(seq(xmin,xmax,0.1), P[,4], col="black")
}

# calculation of probability : Same as P.egrm
#input par=c(a,d), th, xi, cov(item number of previouse response)

g=c(0.2,0,0); u=c(1,0.8,0.8); plotit(-8,8);
p.straight=c(0,0.2,0,0); plotit2(-8,8);





p.straight=c(0.2,0,0,0); plotit2(-8,8);




g=c(0.2,0,0); u=c(1,0.8,0.8); apply(P, 1, sum)
g=c(0.2,0.2,0); u=c(1,1,0.8); apply(P, 1, sum)
g=c(0.2,0.2,0.2); u=c(1, 1, 1); apply(P, 1,sum)


g=c(0,0,0); u=c(1,1,1); plotit(-8,8);
g=c(0.2,0,0); u=c(1,0.8,0.8); plotit(-8,8);

g=c(0,0,0); u=c(1,1,1); plotit(-8,8);
g=c(0,0,0); u=c(0.8,0.8,0.8); plotit(-8,8);

g=c(0,0,0); u=c(1,1,1); plotit(-8,8);
g=c(0.2,0.2,0); u=c(1,1,0.8); plotit(-8,8);

g=c(0,0,0); u=c(1,1,1); plotit(-8,8);
g=c(0.2,0,0); u=c(1,0.8,0.8); plotit(-8,8);

g=c(0.2,0.2,0); u=c(1,1,0.8); plotit(-8,8);

par(mfcol=c(1,3))
g=c(0,0,0); u=c(1,1,1); plotit(-8,8); apply(P, 1, sum)
g=c(0,0,0); u=c(0.7,0.7,0.7); plotit(-8,8); apply(P, 1, sum)
g=c(0.3,0.2,0.1); u=c(1,0.9, 0.8); plotit(-8,8); apply(P, 1, sum)

par(mfcol=c(1,3))
g=c(0,0,0); u=c(1,1,1); plotit(-8,8); apply(P, 1, sum)
g=c(0.3,0,0); u=c(1,0.7,0.7); plotit(-8,8); apply(P, 1, sum)
g=c(0.2,0.2,0.1); u=c(0.9,0.9, 0.8); plotit(-8,8); apply(P, 1, sum)


par(mfcol=c(1,3))
g=c(0,0,0); u=c(1,1,1); plotit(-8,8); apply(P, 1, sum)
g=c(1,0,0); u=c(1,0,0); plotit(-8,8); apply(P, 1, sum)
g=c(0.2,0.2,0.1); u=c(0.9,0.9, 0.8); plotit(-8,8); apply(P, 1, sum)

g=c(0.3, 0, 0); u=c(1, 0.7, 0.7);
g=c(0.3, 0, 0); u=c(1, 0.7, 0.7);

par(mfcol=c(1,2))
g=c(0,0,0); u=c(1,1,1); plotit(-8,8); apply(P, 1, sum)
g=c(0.3, 0, 0); u=c(1, 0.7, 0.7); plotit(-8,8); apply(P, 1, sum)

g=c(0.2, 0.1, 0.1); u=c(0.9, 0.8, 0.8); plotit(-8,8); apply(P, 1, sum)

g.all = 0.3 /1
g.c2 = -0.1

g=c(g.all+g.c2, g.all*2/3+g.c2, g.all/3);
u=c(1-g.all/3, 1-g.all*2/3, 1-g.all-g.c2)
plotit(-8,8); apply(P, 1, sum);

g.all = 0
g.c2 = 0.1

g=c(g.all+g.c2, g.all*2/3+g.c2, g.all/3);
u=c(1-g.all/3, 1-g.all*2/3, 1-g.all-g.c2)
plotit(-8,8); apply(P, 1, sum);



