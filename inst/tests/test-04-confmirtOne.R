context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(onefact, 'ConfirmatoryClass')    
    cfs <- as.numeric(do.call(c, coef(onefact, digits=4)))
    expect_equal(cfs, c(0.9892, 0.5217, 1.4568, 1.8562, 1.5541, 2.1583, 0, NA, NA, 1, NA, NA, 1.0764, 0.4087, 1.7442, 0.8074, 0.5763, 1.0384, 0, NA, NA, 1, NA, NA, 1.7244, 0.8822, 2.5667, 1.815, 1.3459, 2.2841, 0, NA, NA, 1, NA, NA, 0.7553, 0.3079, 1.2026, 0.4849, 0.3339, 0.636, 0, NA, NA, 1, NA, NA, 0.7652, 0.4086, 1.1218, 1.8673, 1.6243, 2.1103, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
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
    expect_true(mirt:::closeEnough(W1$W - 96.686, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(W2$W - 1.728, -1e-2, 1e-2))
    fitonefact <- fitIndices(onefact)
    expect_is(fitonefact, 'data.frame')
    suppressWarnings(twofact <- mirt(fulldata, 2, verbose = FALSE, draws = 10, method = 'MHRM'))
    cfs <- as.numeric(do.call(c, coef(twofact, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.6568, 0.3086, 1.005, -0.5837, -1.3276, 0.1603, 1.8894, 1.6292, 2.1495, 0, NA, NA, 1, NA, NA, 1.3991, 0.9559, 1.8424, 0.1994, -0.1898, 0.5887, 0.8715, 0.6716, 1.0715, 0, NA, NA, 1, NA, NA, 1.5885, 0.7646, 2.4123, -0.1413, -0.7733, 0.4907, 1.7782, 1.369, 2.1875, 0, NA, NA, 1, NA, NA, 0.5764, 0.3311, 0.8217, -0.2592, -0.8007, 0.2823, 0.4802, 0.3313, 0.6291, 0, NA, NA, 1, NA, NA, -0.0103, -0.5703, 0.5497, -1.6973, -1.9067, -1.4879, 2.4348, 2.0092, 2.8604, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, -0.4373, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    expect_is(twofact, 'ExploratoryClass')
    modm7 <- suppressMessages(mirt(fulldata, 1, '4PL', verbose=FALSE, parprior = list(c(3,7,11,15,19,'norm', -1.7, 1), 
                                                                     c(4,8,12,16,20,'norm', 1.7, 1)), method = 'MHRM', draws = 10))
    expect_equal(modm7@df, 11)
    expect_is(modm7, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modm7)))
    expect_equal(cfs, c(5.322, 0.009, 10.636, 7.45, 0.371, 14.529, 0.136, 0.028, 0.462, 0.898, 0.861, 0.926, 11.756, -6.661, 30.172, 4.154, -3.63, 11.938, 0.33, 0.223, 0.457, 0.844, 0.763, 0.901, 5.234, -1.746, 12.213, 2.941, -1.111, 6.992, 0.362, 0.208, 0.552, 0.942, 0.8, 0.985, 1.851, 0.324, 3.377, 1.777, -0.672, 4.225, 0.114, 0.024, 0.399, 0.763, 0.565, 0.889, 1.777, 0.119, 3.435, 3.395, 0.827, 5.963, 0.151, 0.024, 0.56, 0.909, 0.835, 0.951, 0, NA, NA, 1, NA, NA), tollerance = 1e-2)
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- mirt(fulldata, 1, verbose = FALSE, draws = 10, method = 'MHRM')    
    expect_is(onefactmissing, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(onefactmissing, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.9741, 0.5695, 1.3787, 1.8525, 1.5673, 2.1377, 0, NA, NA, 1, NA, NA, 1.0568, 0.6981, 1.4155, 0.8032, 0.6244, 0.982, 0, NA, NA, 1, NA, NA, 1.7092, 0.5155, 2.9028, 1.8053, 1.1386, 2.4719, 0, NA, NA, 1, NA, NA, 0.7559, 0.4482, 1.0635, 0.4845, 0.3343, 0.6347, 0, NA, NA, 1, NA, NA, 0.7594, 0.4153, 1.1036, 1.8645, 1.6169, 2.112, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
        
    fs1 <- fscores(onefact, verbose = FALSE, mean=c(1), cov=matrix(2))
    expect_is(fs1, 'matrix')
    expect_true(mirt:::closeEnough(fs1[1:3,'F1'] - c(-2.148334, -1.681619, -1.687547), -1e-2, 1e-2)) 
    fs2 <- fscores(twofact, verbose = FALSE)
    expect_is(fs2, 'matrix') 
    fs3 <- fscores(onefactmissing, verbose = FALSE)
    expect_is(fs3, 'matrix')
})
 
