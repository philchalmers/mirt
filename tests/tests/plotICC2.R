
P.egrmst <- function(a, d, g, u, th, xi) {
    d.mean=mean(d)
    d.new = exp(xi)*(d-d.mean)+d.mean
    prob=c(1,g+(u-g)/(1+exp(-(a*(th-d.new)))), 0)
    -diff(prob)}

plotit=function(xmin=-4, xmax=4) {
    P=matrix(NA, nrow=1, ncol=ncat);
    for (i.th in seq(xmin,xmax,0.1)) {
        P=rbind(P, P.egrmst(a= a, d=d, g=g, u=u, th=i.th, xi=xi))
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
        P=rbind(P, p.straight + q.straight * P.egrmst(a= a, d=d, g=0, u=1, th=i.th, xi=xi))
    }
    P=P[-1,]

    plot(seq(xmin,xmax,0.1), P[,1], type="l", ylim=c(0,1), col="red")
    lines(seq(xmin,xmax,0.1), P[,2], col="blue")
    lines(seq(xmin,xmax,0.1), P[,3], col="green")
    lines(seq(xmin,xmax,0.1), P[,4], col="black")
}

# ICC when probability of straight response of last response (in this case 2) is 0.2
# It can be parameterized as lower and upper asymptote of Pr(x>=k)
# or mixture model of Pr(straight) = 0.2 and Pr(GRM) = 0.8
par(mfcol=c(1,2))
xi=0; a=1; d=c(-1, 0, 1)
g=c(0.2,0,0); u=c(1,0.8,0.8); plotit(-4,4);
p.straight=c(0,0.2,0,0); plotit2(-4,4);



p.st <- rnorm(10000); a.st=1
p.repeat <- 2*abs(1/(1+exp(-a.st*p.st))-1/2)
#p.repeat
hist(p.repeat, xlim=c(0,1), f=F)



p.st=6




a.st <- 1/3; ncat=4; Theta=matrix(c(1,3,-2), ncol=1)
p.st <- matrix(Theta[,1], ncol=ncat, nrow=nrow(Theta))
p.repeat <- 2*abs(1/(1+exp(-a.st*p.st))-1/2)
q.repeat <- 1- p.repeat

P <- cbind(P.grm + p.repeat, P.grm * q.repeat)

