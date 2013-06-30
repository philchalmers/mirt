context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- confmirt(fulldata, 1, verbose = FALSE, draws = 10)
    expect_is(onefact, 'ConfirmatoryClass')    
    cfs <- as.numeric(do.call(c, coef(onefact, digits=4)))
    expect_equal(cfs, c(0.9892, 0.2385, 1.8562, 0.1541, 0, NA, 1, NA, 1.0764, 0.3407, 0.8074,
                        0.1179, 0, NA, 1, NA, 1.7244, 0.4297, 1.815, 0.2393, 0, NA, 1, NA, 0.7553, 
                        0.2282, 0.4849, 0.0771, 0, NA, 1, NA, 0.7652, 0.1819, 1.8673, 0.124, 0, NA, 
                        1, NA, 0, NA, 1, NA),
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
    expect_is(fitonefact, 'list')
    suppressWarnings(twofact <- confmirt(fulldata, 2, verbose = FALSE, draws = 10))
    cfs <- as.numeric(do.call(c, coef(twofact, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.4859, 0.1904, -0.7853, 0.1658, 1.9218, 0.1433, 0, NA, 1, NA, 1.4393, 
                        0.1143, 0.1872, 0.0945, 0.8808, 0.0921, 0, NA, 1, NA, 1.6235, 0.1142, 
                        -0.202, 0.2568, 1.8209, 0.1615, 0, NA, 1, NA, 0.4775, 0.1019, -0.381, 
                        0.1464, 0.4834, 0.072, 0, NA, 1, NA, -0.0591, 0.1688, -1.4161, 0.0637, 
                        2.2265, 0.122, 0, NA, 1, NA, 0, NA, 0, NA, 1, NA, -0.4901, NA, 1, NA),
                 tollerance = 1e-2)
    expect_is(twofact, 'ExploratoryClass')
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- confmirt(fulldata, 1, verbose = FALSE, draws = 10)    
    expect_is(onefactmissing, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(onefactmissing, digits=4, verbose = FALSE)))
    expect_equal(cfs, c(0.9741, 0.2064, 1.8525, 0.1455, 0, NA, 1, NA, 1.0568, 0.183, 0.8032, 
                        0.0912, 0, NA, 1, NA, 1.7092, 0.609, 1.8053, 0.3401, 0, NA, 1, NA, 0.7559, 
                        0.157, 0.4845, 0.0766, 0, NA, 1, NA, 0.7594, 0.1756, 1.8645, 0.1263, 0, 
                        NA, 1, NA, 0, NA, 1, NA),
                 tollerance = 1e-2)
        
    fs1 <- fscores(onefact, verbose = FALSE, mean=c(1), cov=matrix(2))
    expect_is(fs1, 'matrix')
    expect_true(mirt:::closeEnough(fs1[1:3,'F1'] - c(-2.148334, -1.681619, -1.687547), -1e-2, 1e-2)) 
    fs2 <- fscores(twofact, verbose = FALSE)
    expect_is(fs2, 'matrix') 
    fs3 <- fscores(onefactmissing, verbose = FALSE)
    expect_is(fs3, 'matrix')
})
 
