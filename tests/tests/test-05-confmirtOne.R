context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- mirt(fulldata, 1, verbose = FALSE, SE.type='MHRM', SE=TRUE)
    expect_is(onefact, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(onefact, digits=4)))
    expect_equal(cfs, c(0.9879,0.6413,1.3346,1.8561,1.6051,2.107,0,NA,NA,1,NA,NA,1.0809,0.7582,1.4036,0.808,0.635,0.9809,0,NA,NA,1,NA,NA,1.7058,1.1911,2.2205,1.8042,1.4441,2.1644,0,NA,NA,1,NA,NA,0.7652,0.4682,1.0622,0.486,0.3363,0.6356,0,NA,NA,1,NA,NA,0.7358,0.429,1.0426,1.8545,1.6205,2.0885,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)
    names <- wald(onefact)
    L <- matrix(0, 1, ncol(names))
    L[1, c(1,3,5,7,9)] <- 1
    L2 <- matrix(0, 2, ncol(names))
    L2[1, 1] <- L2[2, 3] <- 1
    L2[1, 7] <- L2[2, 9] <- -1
    W1 <- wald(onefact, L)
    W2 <- wald(onefact, L2)
    expect_true(mirt:::closeEnough(W1$W - 200.4343, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(W2$W - 3.237525, -1e-2, 1e-2))

    fitonefact <- M2(onefact)
    expect_is(fitonefact, 'data.frame')
    expect_equal(fitonefact$M2, 11.93769, tolerance = 1e-2)
    twofact <- mirt(fulldata, 2, verbose = FALSE, draws = 10, method = 'MHRM')
    cfs <- as.numeric(do.call(c, coef(twofact, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(-1.4813,0.2851,2.1789,0,1,-1.0415,-1.4008,0.9989,0,1,-1.2904,-0.7484,1.6994,0,1,-0.8552,-0.0268,0.505,0,1,-0.7969,0,1.8894,0,1,0,0,1,0,1),
                 tolerance = 1e-2)
    expect_is(twofact, 'SingleGroupClass')
    expect_message(modm7 <- mirt(fulldata, 1, '4PL', verbose=FALSE, parprior = list(c(3,7,11,15,19,'norm', -1.7, 1),
                                                                     c(4,8,12,16,20,'norm', 1.7, 1)), method = 'MHRM', draws = 10),
                            "MHRM terminated after 2000 iterations.")
    expect_equal(modm7@df, 11)
    expect_is(modm7, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(modm7)))
    expect_equal(cfs, c(3.682,5.396,0.135,0.899,8.241,1.506,0.339,0.896,6.646,3.843,0.359,0.937,3.94,4.184,0.121,0.704,1.934,3.644,0.153,0.904,0,1), tolerance = 1e-2)
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(onefactmissing, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(onefactmissing, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.9621,1.8471,0,1,1.0781,0.8082,0,1,1.7074,1.8052,0,1,0.7369,0.4818,0,1,0.753,1.862,0,1,0,1),
                 tolerance = 1e-2)

    fs1 <- fscores(onefact, verbose = FALSE, mean=c(1), cov=matrix(2))
    expect_is(fs1, 'matrix')
    expect_true(mirt:::closeEnough(fs1[1:3,'F1'] - c(-2.1821, -1.6989, -1.6807), -1e-2, 1e-2))
    fs2 <- fscores(twofact, verbose = FALSE)
    expect_is(fs2, 'matrix')
    fs3 <- fscores(onefactmissing, verbose = FALSE)
    expect_is(fs3, 'matrix')
})

