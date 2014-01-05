context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(onefact, 'ConfirmatoryClass')    
    cfs <- as.numeric(do.call(c, coef(onefact, digits=4)))
    expect_equal(cfs, c(0.9872, 0.6398, 1.3346, 1.8552, 1.5939, 2.1165, 0, NA, NA, 1, NA, NA, 1.0775, 0.7408, 1.4143, 0.8071, 0.6324, 0.9819, 0, NA, NA, 1, NA, NA, 1.7159, 1.0693, 2.3624, 1.8095, 1.4294, 2.1896, 0, NA, NA, 1, NA, NA, 0.7506, 0.5026, 0.9985, 0.4842, 0.3232, 0.6452, 0, NA, NA, 1, NA, NA, 0.7582, 0.4499, 1.0665, 1.8643, 1.609, 2.1195, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
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
    expect_true(mirt:::closeEnough(W1$W - 212.6877, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(W2$W - 2.080472, -1e-2, 1e-2))
    fitonefact <- fitIndices(onefact)
    expect_is(fitonefact, 'data.frame')
    suppressWarnings(twofact <- mirt(fulldata, 2, verbose = FALSE, draws = 10, method = 'MHRM'))
    cfs <- as.numeric(do.call(c, coef(twofact, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.6554, 0.3298, 0.9809, -0.5402, -0.8506, -0.2299, 1.8729, 1.6067, 2.1391, 0, NA, NA, 1, NA, NA, 1.3934, 1.1936, 1.5931, 0.1951, -0.1692, 0.5595, 0.8725, 0.7105, 1.0346, 0, NA, NA, 1, NA, NA, 1.6312, 1.4964, 1.7659, -0.1309, -0.4369, 0.1752, 1.8022, 1.7488, 1.8556, 0, NA, NA, 1, NA, NA, 0.5695, 0.3373, 0.8016, -0.2449, -0.4873, -0.0024, 0.4801, 0.3345, 0.6258, 0, NA, NA, 1, NA, NA, -0.0045, -0.4225, 0.4136, -1.7778, -1.9745, -1.581, 2.5053, 2.111, 2.8995, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, -0.4354, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    expect_is(twofact, 'ExploratoryClass')
    modm7 <- suppressMessages(mirt(fulldata, 1, '4PL', verbose=FALSE, parprior = list(c(3,7,11,15,19,'norm', -1.7, 1), 
                                                                     c(4,8,12,16,20,'norm', 1.7, 1)), method = 'MHRM', draws = 10))
    expect_equal(modm7@df, 11)
    expect_is(modm7, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modm7)))
    expect_equal(cfs, c(3.2, 1.055, 5.346, 4.795, 2.944, 6.645, 0.125, 0.028, 0.419, 0.902, 0.869, 0.928, 8.599, 4.618, 12.58, 1.951, 0.132, 3.77, 0.331, 0.215, 0.471, 0.886, 0.695, 0.963, 3.978, -7.554, 15.51, 2.637, -1.829, 7.103, 0.313, 0.056, 0.779, 0.942, 0.862, 0.976, 2.893, 1.685, 4.101, 3.095, 1.28, 4.91, 0.118, 0.024, 0.419, 0.715, 0.676, 0.75, 1.92, 0.072, 3.768, 3.628, 0.705, 6.55, 0.147, 0.023, 0.558, 0.905, 0.833, 0.948, 0, NA, NA, 1, NA, NA), tollerance = 1e-2)
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')    
    expect_is(onefactmissing, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(onefactmissing, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.9863, 0.6452, 1.3274, 1.8592, 1.6216, 2.0967, 0, NA, NA, 1, NA, NA, 1.0755, 0.7674, 1.3836, 0.8079, 0.6293, 0.9866, 0, NA, NA, 1, NA, NA, 1.6688, 1.1644, 2.1731, 1.7831, 1.4462, 2.12, 0, NA, NA, 1, NA, NA, 0.7533, 0.5088, 0.9978, 0.4845, 0.3287, 0.6402, 0, NA, NA, 1, NA, NA, 0.7705, 0.4204, 1.1206, 1.8698, 1.6072, 2.1324, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
        
    fs1 <- fscores(onefact, verbose = FALSE, mean=c(1), cov=matrix(2))
    expect_is(fs1, 'matrix')
    expect_true(mirt:::closeEnough(fs1[1:3,'F1'] - c(-2.148334, -1.681619, -1.687547), -1e-2, 1e-2)) 
    fs2 <- fscores(twofact, verbose = FALSE)
    expect_is(fs2, 'matrix') 
    fs3 <- fscores(onefactmissing, verbose = FALSE)
    expect_is(fs3, 'matrix')
})
 
