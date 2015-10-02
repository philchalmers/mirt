# TEST : itemtype = "grsm2", itemtype="grsm3"

R <-
    structure(c(1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 3, 1, 2, 1, 2, 1, 1,
                2, 1, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 2, 2, 1,
                1, 1, 2, 2, 3, 3, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 2, 1, 2, 4,
                1, 2, 2, 2, 1, 1, 4, 1, 1, 2, 1, 1, 3, 2, 1, 4, 1, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 2, 3,
                2, 4, 1, 3, 2, 2, 1, 1, 1, 3, 3, 2, 3, 2, 3, 3, 2, 1, 2, 3, 3,
                3, 1, 1, 1, 3, 3, 2, 2, 3, 1, 3, 3, 1, 2, 2, 3, 1, 3, 2, 2, 3,
                3, 1, 2, 3, 1, 3, 3, 2, 3, 2, 3, 3, 2, 1, 2, 4, 3, 2, 3, 2, 2,
                1, 4, 1, 3, 3, 1, 3, 3, 2, 3, 4, 1, 3, 3, 3, 2, 2, 3, 3, 3, 2,
                2, 3, 1, 2, 1, 1, 1, 2, 2, 3, 3, 3, 1, 3, 2, 3, 2, 1, 4, 2, 2,
                2, 1, 1, 1, 4, 1, 4, 4, 2, 3, 3, 2, 3, 2, 4, 3, 4, 4, 1, 1, 4,
                4, 3, 2, 3, 1, 3, 1, 1, 2, 2, 2, 1, 3, 2, 2, 2, 4, 1, 4, 1, 1,
                4, 3, 1, 3, 3, 4, 3, 2, 1, 2, 4, 3, 4, 4, 2, 1, 1, 4, 1, 3, 3,
                1, 1, 3, 2, 4, 4, 1, 3, 3, 3, 2, 4, 3, 4, 1, 2, 2, 2, 1, 1, 1,
                1, 1, 2, 2, 1, 2, 3, 1, 1, 2, 3, 2, 1, 3, 1, 2, 1, 1, 1, 1, 4,
                1, 4, 1, 2, 3, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 3, 1, 2, 2, 1,
                3, 1, 1, 2, 3, 1, 1, 2, 2, 3, 3, 4, 1, 2, 1, 1, 1, 1, 1, 4, 1,
                3, 3, 2, 1, 2, 1, 3, 4, 3, 1, 3, 1, 4, 1, 3, 1, 1, 1, 3, 2, 1,
                2, 1, 2, 2, 1, 2, 3, 1, 3, 1, 2, 2, 1, 1, 3, 1, 1, 1, 2, 2, 2,
                1, 2, 1, 3, 2, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 4, 4, 1, 4, 2, 2,
                1, 1, 1, 2, 1, 1, 4, 1, 1, 1, 1, 3, 1, 2, 3, 1, 1, 1, 1, 2, 3,
                3, 1, 2, 1, 3, 2, 4, 1, 2, 1, 1, 3, 1, 1, 4, 1, 3, 3, 2, 1, 2,
                1, 3, 4, 1, 1, 1, 1, 2, 1, 1, 3, 1, 1, 3, 2, 1, 1, 1, 4, 4, 2,
                2, 4, 1, 3, 3, 2, 2, 2, 1, 3, 1, 1, 1, 2, 2, 1, 1, 2, 1, 3, 2
    ), .Dim = c(100L, 5L))


ncat = function(R) {
    stopifnot(is.matrix(R))
    max(apply(R, 2, function(x) length(unique(x[!is.na(x)]))))
}

fit.grsm <- mirt(R, 1, itemtype="grsm")

n.item <- ncol(R); ncat <- ncat(R)
model <- paste("F1 = 1-",n.item,
               "\nCONSTRAIN = (1-",n.item,",d1)",
               paste(",(1-",n.item,",d",2:(ncat-1),")",sep="", collapse=""),
               sep="")
pars <- mod2values(fit.grsm);
str(pars[,"class"])
pars[,"class"]=factor(c(rep("grsm2.i1", ncat+1), rep("grsm2", (nrow(pars)-ncat-3)), rep("GroupPars",2)))
#pars[1:(ncat+1),"class"]="grsm2.i1"
#pars[(ncat+2):(nrow(pars)-2),"class"]="grsm2"
fit.grsm2 <- mirt(R, mirt.model(model), pars=pars, itemtype=c("grsm2.i1",rep("grsm2", n.item-1)),
                  TOL=1e-5)

fit.grsm2_ <- mirt(R, mirt.model(model), itemtype=c("grsm2.i1",rep("grsm2", n.item-1)))

pars[,"class"]=factor(c(rep("grsm3.i1", ncat+1), rep("grsm3", (nrow(pars)-ncat-3)), rep("GroupPars",2)))
fit.grsm3 <- mirt(R, mirt.model(model), pars=pars, itemtype=c("grsm3.i1",rep("grsm3", n.item-1)))

coef2mat = function(cf) {
    stopifnot(is.list(cf))
    do.call(rbind, cf[1:(length(cf)-1)])
}


coef2mat(coef(fit.grsm))
coef2mat(coef(fit.grsm2))
coef2mat(coef(fit.grsm2_))
coef2mat(coef(fit.grsm3))

fit.grsm@logLik
fit.grsm2@logLik
fit.grsm2_@logLik
fit.grsm3@logLik



# Do START and FIXED work?

model <- paste("F1 = 1-",n.item,
               "\nCONSTRAIN = (1-",n.item,",d1)",
               paste(",(1-",n.item,",d",2:(ncat-1),")",sep="", collapse=""),
               "\nSTART = (1, c, 0.0)",
               "\nFIXED = (1, c)",
               sep="")

fit.grsm3_ <- mirt(R, mirt.model(model), itemtype=rep("grsm3", n.item))

coef2mat(coef(fit.grsm3))
coef2mat(coef(fit.grsm3_))




















#===== IGNORE THE BELOW ==================


fit.grsm2 <- mirt(R, mirt.model(model), itemtype=c("grsm2.i1",rep("grsm2", n.item-1)),
                  TOL=1e-9)

fit.grsm3 <- mirt(R, mirt.model(model), itemtype=c("grsm3.i1",rep("grsm3", n.item-1)))




mod2values(fit.grsm)
mirt(R, 1, itemtype="grsm", pars="values")
mirt(R, 1, itemtype="grsm")
coef(fit.grsm)



mirt(R, mirt.model(model), itemtype="grsm2", pars="values")
fit.grsm2.constrain <- mirt(R, mirt.model(model), itemtype="grsm2")
#vals <- mirt(R, mirt.model(model), itemtype="graded.b", pars="values")
#mirt(R, mirt.model(model), itemtype="grsm2", pars="values")

#val <- mod2values(fit.grsm2)
#val[val$item=="Item.1" & val$name=="c","value"]=0
#val[val$item=="Item.1" & val$name=="c","est"]=F
#fit.grsm2 <- mirt(R, mirt.model(model), vals=val, itemtype="grsm2")
#fit.grsm2 <- mirt(R, mirt.model(model), itemtype="grsm2")
#fit.grsm2 <- mirt(R, mirt.model(model), itemtype=c("grsm2.i1",rep("grsm2",n.item-1)))

ncat = function(R) {
    stopifnot(is.matrix(R))
    max(apply(R, 2, function(x) length(unique(x[!is.na(x)]))))
}

n.item <- ncol(R); ncat <- ncat(R)
model <- paste("F1 = 1-",n.item,
               "\nCONSTRAIN = (1-",n.item,",d1)",
               paste(",(1-",n.item,",d",2:(ncat-1),")",sep="", collapse=""),
               "\nSTART = (1, c, 0.0)",
               "\nFIXED = (1, c)",
               sep="")
fit.grsm2 <- mirt(R, mirt.model(model), itemtype="grsm2", TOL=1e-5)

fit.grsm2 <- mirt(R, mirt.model(model), pars = mod2values(fit.grsm2), itemtype="grsm2", TOL=1e-7)
fit.grsm2 <- mirt(R, mirt.model(model), vals=val, itemtype=c("grsm2.i1",rep("grsm2", n.item-1)))
fit.grsm2 <- mirt(R, mirt.model(model), itemtype=c("grsm2.i1",rep("grsm2", n.item-1)))

val <- mirt(R, mirt.model(model), pars="values", itemtype=c("grsm2.i1",rep("grsm2", n.item-1)))
val[val$item!="Item.1" & val$name=="c","value"]=0
val[val$item=="Item.1" & val$name=="c","est"]=F

#vals <- mirt(R, mirt.model(model), itemtype="graded.b", pars="values")
#mirt(R, mirt.model(model), itemtype="grsm2", pars="values")

#val <- mod2values(fit.grsm2)
#



fit.grsm2a <- mirt(R, mirt.model(model), vals=val, itemtype=rep("grsm2", n.item))


coef(fit.grsm2)
#fit.grsm2 <- mirt(R, mirt.model(model), itemtype="grsm2")


coef(fit.grsm2)
coef(fit.grsm2a)

fit.grsm;
fit.grsm2;
mod2values(fit.grsm2)

grm(R)




a <- matrix(rlnorm(20,.2,.3))

# for the graded model, ensure that there is enough space between the intercepts,
# otherwise closer categories will not be selected often (minimum distance of 0.3 here)
diffs <- t(apply(matrix(runif(20*4, .3, 1), 20), 1, cumsum))
diffs <- -(diffs - rowMeans(diffs))
d <- diffs + rnorm(20)

R <- simdata(a, d, 500, itemtype = 'graded')
