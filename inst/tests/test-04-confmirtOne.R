context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(onefact, 'ConfirmatoryClass')    
    cfs <- as.numeric(do.call(c, coef(onefact, digits=4)))
    expect_equal(cfs, c(0.9872, 0.5992, 1.3752, 1.8552, 1.5837, 2.1267, 0, NA, NA, 1, NA, NA, 1.0775, 0.6985, 1.4566, 0.8071, 0.6307, 0.9836, 0, NA, NA, 1, NA, NA, 1.7159, 0.9926, 2.4391, 1.8095, 1.3747, 2.2442, 0, NA, NA, 1, NA, NA, 0.7506, 0.4946, 1.0066, 0.4842, 0.322, 0.6464, 0, NA, NA, 1, NA, NA, 0.7582, 0.4457, 1.0706, 1.8643, 1.6067, 2.1218, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    names <- wald(onefact)
    L <- matrix(0, 1, length(names))
    L[1, c(1,3,5,7,9)] <- 1
    L2 <- matrix(0, 2, length(names))
    L2[1, 1] <- L2[2, 3] <- 1
    L2[1, 7] <- L2[2, 9] <- -1
    W1 <- wald(onefact, L)
    W2 <- wald(onefact, L2)    
    expect_is(W1, 'wald')
    expect_is(W2, 'wald')
    expect_true(mirt:::closeEnough(W1$W - 209.2077, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(W2$W - 1.515, -1e-2, 1e-2))
    fitonefact <- fitIndices(onefact)
    expect_is(fitonefact, 'data.frame')
    suppressWarnings(twofact <- mirt(fulldata, 2, verbose = FALSE, draws = 10, method = 'MHRM'))
    cfs <- as.numeric(do.call(c, coef(twofact, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.6554, 0.3271, 0.9837, -0.5402, -0.8868, -0.1937, 1.8729, 1.6057, 2.1401, 0, NA, NA, 1, NA, NA, 1.3934, 1.0633, 1.7234, 0.1951, -0.2438, 0.634, 0.8725, 0.71, 1.0351, 0, NA, NA, 1, NA, NA, 1.6312, 1.5436, 1.7188, -0.1309, -0.4142, 0.1525, 1.8022, 1.6739, 1.9305, 0, NA, NA, 1, NA, NA, 0.5695, 0.3362, 0.8027, -0.2449, -0.494, 0.0043, 0.4801, 0.3342, 0.626, 0, NA, NA, 1, NA, NA, -0.0045, -0.4506, 0.4417, -1.7778, -1.9886, -1.5669, 2.5053, 2.0813, 2.9292, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, -0.4354, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    expect_is(twofact, 'ExploratoryClass')
    modm7 <- suppressMessages(mirt(fulldata, 1, '4PL', verbose=FALSE, parprior = list(c(3,7,11,15,19,'norm', -1.7, 1), 
                                                                     c(4,8,12,16,20,'norm', 1.7, 1)), method = 'MHRM', draws = 10))
    expect_equal(modm7@df, 11)
    expect_is(modm7, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modm7)))
    expect_equal(cfs, c(3.2, 1.032, 5.368, 4.795, 2.938, 6.651, 0.125, 0.027, 0.422, 0.902, 0.868, 0.928, 8.599, 2.472, 14.726, 1.951, -0.24, 4.142, 0.331, 0.215, 0.472, 0.886, 0.642, 0.971, 3.978, -9.692, 17.648, 2.637, -2.853, 8.127, 0.313, 0.045, 0.815, 0.942, 0.854, 0.978, 2.893, 1.686, 4.1, 3.095, 1.287, 4.903, 0.118, 0.024, 0.424, 0.715, 0.676, 0.751, 1.92, 0.071, 3.769, 3.628, 0.705, 6.55, 0.147, 0.023, 0.559, 0.905, 0.832, 0.948, 0, NA, NA, 1, NA, NA), tollerance = 1e-2)
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')    
    expect_is(onefactmissing, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(onefactmissing, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.9863, 0.6398, 1.3328, 1.8592, 1.6213, 2.097, 0, NA, NA, 1, NA, NA, 1.0755, 0.7597, 1.3913, 0.8079, 0.6293, 0.9866, 0, NA, NA, 1, NA, NA, 1.6688, 1.1643, 2.1732, 1.7831, 1.4452, 2.121, 0, NA, NA, 1, NA, NA, 0.7533, 0.5069, 0.9997, 0.4845, 0.3279, 0.641, 0, NA, NA, 1, NA, NA, 0.7705, 0.42, 1.1211, 1.8698, 1.6063, 2.1333, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
        
    fs1 <- fscores(onefact, verbose = FALSE, mean=c(1), cov=matrix(2))
    expect_is(fs1, 'matrix')
    expect_true(mirt:::closeEnough(fs1[1:3,'F1'] - c(-2.148334, -1.681619, -1.687547), -1e-2, 1e-2)) 
    fs2 <- fscores(twofact, verbose = FALSE)
    expect_is(fs2, 'matrix') 
    fs3 <- fscores(onefactmissing, verbose = FALSE)
    expect_is(fs3, 'matrix')
})
 
