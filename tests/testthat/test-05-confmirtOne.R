expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- mirt(fulldata, 1, verbose = FALSE, SE.type='MHRM', SE=TRUE)
    expect_class(onefact, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(onefact)))
    expect_equal(cfs, c(0.9879254,0.5978895,1.377961,1.85606,1.580046,2.132075,0,NA,NA,1,NA,NA,1.080885,0.772476,1.389293,0.8079786,0.6345122,0.981445,0,NA,NA,1,NA,NA,1.705801,0.8574762,2.554125,1.804219,1.308859,2.299578,0,NA,NA,1,NA,NA,0.7651853,0.4655423,1.064828,0.4859966,0.3374257,0.6345675,0,NA,NA,1,NA,NA,0.735798,0.4468093,1.024787,1.854513,1.62773,2.081295,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)
    names <- wald(onefact)
    L <- matrix(0, 1, length(names))
    L[1, c(1,3,5,7,9)] <- 1
    L2 <- matrix(0, 2, length(names))
    L2[1, 1] <- L2[2, 3] <- 1
    L2[1, 7] <- L2[2, 9] <- -1
    W1 <- wald(onefact, L)
    W2 <- wald(onefact, L2)
    expect_true(mirt:::closeEnough(W1$W - 178.6441, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(W2$W - 4.228611, -1e-2, 1e-2))
    onefactSEM <- mirt(fulldata, 1, verbose = FALSE, method='SEM')
    cfs <- as.numeric(do.call(c, coef(onefactSEM, verbose = FALSE)))
    expect_equal(cfs, c(0.97,1.843,0,1,1.08,0.805,0,1,1.741,1.819,0,1,
                        0.766,0.484,0,1,0.726,1.848,0,1,0,1), tolerance = 1e-2)


    fitonefact <- M2(onefact)
    expect_class(fitonefact, 'data.frame')
    expect_equal(fitonefact$M2, 11.93769, tolerance = 1e-2)
    twofact <- mirt(fulldata, 2, verbose = FALSE, draws = 10, method = 'MHRM')
    cfs <- as.numeric(do.call(c, coef(twofact, verbose = FALSE)))
    expect_equal(cfs, c(-1.389256,0.3038708,2.108244,0,1,-1.055821,-1.13757,0.9310239,0,1,-1.376235,-0.6974982,1.71311,0,1,-0.8597142,0.05825403,0.499093,0,1,-0.8210638,0,1.893461,0,1,0,0,1,0,1),
                 tolerance = 1e-2)
    expect_class(twofact, 'SingleGroupClass')
    modm7 <- mirt(fulldata, 1, '4PL', verbose=FALSE, parprior = list(c(3,7,11,15,19,'norm', -1.7, .25),
                                                                     c(4,8,12,16,20,'norm', 1.7, .25)), method = 'MHRM', draws = 10)
    expect_equal(extract.mirt(modm7, 'df'), 11)
    expect_class(modm7, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(modm7)))
    expect_equal(cfs, c(2.503,4.279,0.153,0.883,3.114,1.908,0.183,0.852,8.703,8.461,0.157,0.89,1.331,0.874,0.15,0.843,2.568,5.007,0.154,0.878,0,1), tolerance = 1e-2)
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_class(onefactmissing, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(onefactmissing, verbose = FALSE)))
    expect_equal(cfs, c(0.9366883,1.832658,0,1,1.070041,0.804908,0,1,1.770225,1.837219,0,1,0.7649542,0.4847653,0,1,0.7494625,1.859265,0,1,0,1),
                 tolerance = 1e-2)

    fs1 <- fscores(onefact, verbose = FALSE, mean=c(1), cov=matrix(2), full.scores=FALSE)
    expect_class(fs1, 'matrix')
    expect_true(mirt:::closeEnough(fs1[1:3,'F1'] - c(-2.1821, -1.6989, -1.6807), -1e-2, 1e-2))
    fs2 <- fscores(twofact, verbose = FALSE, full.scores=FALSE)
    expect_class(fs2, 'matrix')
    fs3 <- fscores(onefactmissing, verbose = FALSE, full.scores=FALSE)
    expect_class(fs3, 'matrix')
})

