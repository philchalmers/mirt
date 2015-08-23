context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- mirt(fulldata, 1, verbose = FALSE, SE.type='MHRM', SE=TRUE)
    expect_is(onefact, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(onefact, digits=4)))
    expect_equal(cfs, c(0.9879,0.6579,1.3179,1.8561,1.6074,2.1048,0,NA,NA,1,NA,NA,1.0809,0.7741,1.3877,0.808,0.6389,0.9771,0,NA,NA,1,NA,NA,1.7058,1.1881,2.2235,1.8042,1.4628,2.1457,0,NA,NA,1,NA,NA,0.7652,0.51,1.0203,0.486,0.3402,0.6318,0,NA,NA,1,NA,NA,0.7358,0.4549,1.0167,1.8545,1.6328,2.0762,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)
    names <- wald(onefact)
    L <- matrix(0, 1, ncol(names))
    L[1, c(1,3,5,7,9)] <- 1
    L2 <- matrix(0, 2, ncol(names))
    L2[1, 1] <- L2[2, 3] <- 1
    L2[1, 7] <- L2[2, 9] <- -1
    W1 <- wald(onefact, L)
    W2 <- wald(onefact, L2)
    expect_true(mirt:::closeEnough(W1$W - 187.1901, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(W2$W - 3.886675, -1e-2, 1e-2))

    fitonefact <- M2(onefact)
    expect_is(fitonefact, 'data.frame')
    expect_equal(fitonefact$M2, 11.93769, tolerance = 1e-2)
    twofact <- mirt(fulldata, 2, verbose = FALSE, draws = 10, method = 'MHRM')
    cfs <- as.numeric(do.call(c, coef(twofact, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(-1.4378,0.2262,2.1215,0,1,-1.0463,-1.439,1.0058,0,1,-1.291,-0.7683,1.6927,0,1,-0.8431,-0.0161,0.4919,0,1,-0.7941,0,1.8772,0,1,0,0,1,0,1),
                 tolerance = 1e-2)
    expect_is(twofact, 'SingleGroupClass')
    expect_message(modm7 <- mirt(fulldata, 1, '4PL', verbose=FALSE, parprior = list(c(3,7,11,15,19,'norm', -1.7, 1),
                                                                     c(4,8,12,16,20,'norm', 1.7, 1)), method = 'MHRM', draws = 10),
                            "MHRM terminated after 2000 iterations.")
    expect_equal(modm7@df, 11)
    expect_is(modm7, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(modm7)))
    expect_equal(cfs, c(5.103,7.255,0.148,0.895,8.396,1.902,0.323,0.891,10.745,6.593,0.356,0.926,5.189,5.404,0.148,0.696,2.352,3.802,0.338,0.897,0,1), tolerance = 1e-2)
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(onefactmissing, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(onefactmissing, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.969,1.8464,0,1,1.0833,0.806,0,1,1.6824,1.7834,0,1,0.7595,0.4829,0,1,0.7681,1.8658,0,1,0,1),
                 tolerance = 1e-2)

    fs1 <- fscores(onefact, verbose = FALSE, mean=c(1), cov=matrix(2))
    expect_is(fs1, 'matrix')
    expect_true(mirt:::closeEnough(fs1[1:3,'F1'] - c(-2.1821, -1.6989, -1.6807), -1e-2, 1e-2))
    fs2 <- fscores(twofact, verbose = FALSE)
    expect_is(fs2, 'matrix')
    fs3 <- fscores(onefactmissing, verbose = FALSE)
    expect_is(fs3, 'matrix')
})

