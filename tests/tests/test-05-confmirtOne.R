context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM', SE=TRUE)
    expect_is(onefact, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(onefact, digits=4)))
    expect_equal(cfs, c(0.9872,0.6067,1.3677,1.8552,1.5816,2.1289,0,NA,NA,1,NA,NA,1.0775,0.7526,1.4024,0.8071,0.619,0.9953,0,NA,NA,1,NA,NA,1.7159,0.8791,2.5526,1.8095,1.2797,2.3392,0,NA,NA,1,NA,NA,0.7506,0.4513,1.0499,0.4842,0.3336,0.6348,0,NA,NA,1,NA,NA,0.7582,0.4528,1.0636,1.8643,1.6372,2.0913,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)
    names <- wald(onefact)
    L <- matrix(0, 1, ncol(names))
    L[1, c(1,3,5,7,9)] <- 1
    L2 <- matrix(0, 2, ncol(names))
    L2[1, 1] <- L2[2, 3] <- 1
    L2[1, 7] <- L2[2, 9] <- -1
    W1 <- wald(onefact, L)
    W2 <- wald(onefact, L2)
    expect_true(mirt:::closeEnough(W1$W - 198.0846, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(W2$W - 3.502625, -1e-2, 1e-2))    
    
    fitonefact <- M2(onefact)
    expect_is(fitonefact, 'data.frame')
    expect_equal(fitonefact$M2, 11.92959, tolerance = 1e-2)
    twofact <- mirt(fulldata, 2, verbose = FALSE, draws = 10, method = 'MHRM')
    cfs <- as.numeric(do.call(c, coef(twofact, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.0023,-1.3362,2.0565,0,1,1.8471,0.1089,1.0121,0,1,1.0655,-0.6646,1.7122,0,1,0.185,-0.7199,0.4962,0,1,-0.0826,-0.9982,1.9631,0,1,0,0,1,-0.5326,1),
                 tolerance = 1e-2)
    expect_is(twofact, 'ExploratoryClass')
    expect_message(modm7 <- mirt(fulldata, 1, '4PL', verbose=FALSE, parprior = list(c(3,7,11,15,19,'norm', -1.7, 1),
                                                                     c(4,8,12,16,20,'norm', 1.7, 1)), method = 'MHRM', draws = 10),
                            "MHRM terminated after 2000 iterations.")
    expect_equal(modm7@df, 11)
    expect_is(modm7, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modm7)))
    expect_equal(cfs, c(3.2,4.795,0.125,0.902,8.599,1.951,0.331,0.886,3.978,2.637,0.313,0.942,2.893,3.095,0.118,0.715,1.92,3.628,0.147,0.905,0,1), tolerance = 1e-2)
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(onefactmissing, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(onefactmissing, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.9863,1.8592,0,1,1.0755,0.8079,0,1,1.6688,1.7831,0,1,0.7533,0.4845,0,1,0.7705,1.8698,0,1,0,1),
                 tolerance = 1e-2)

    fs1 <- fscores(onefact, verbose = FALSE, mean=c(1), cov=matrix(2))
    expect_is(fs1, 'matrix')
    expect_true(mirt:::closeEnough(fs1[1:3,'F1'] - c(-2.185196, -1.689718, -1.694422), -1e-2, 1e-2))
    fs2 <- fscores(twofact, verbose = FALSE)
    expect_is(fs2, 'matrix')
    fs3 <- fscores(onefactmissing, verbose = FALSE)
    expect_is(fs3, 'matrix')
})

