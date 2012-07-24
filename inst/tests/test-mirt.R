context('mirt')

test_that('dich', {
    data <- expand.table(LSAT7)        
    modm1 <- mirt(data, 1)
    modm2 <- mirt(data, 2)
    expect_is(modm1, 'mirtClass')          
    expect_is(modm1, 'mirtClass')
    
    fm1 <- fscores(modm1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    fm2 <- fscores(modm2, verbose = FALSE)
    expect_is(fm2, 'matrix')
})

test_that('poly', {
    modp1 <- mirt(Science, 1)
    modp2 <- mirt(Science, 2)
    expect_is(modp1, 'mirtClass')          
    expect_is(modp2, 'mirtClass')
    
    fm1 <- fscores(modp1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    fm2 <- fscores(modp2, rotate = 'oblimin', verbose = FALSE)
    expect_is(fm2, 'matrix')
})

