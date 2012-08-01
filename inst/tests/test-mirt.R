context('mirt')

test_that('dich', {
    data <- expand.table(LSAT7)        
    modm1 <- mirt(data, 1, SE=TRUE)
    modm2 <- mirt(data, 2)
    modm3 <- mirt(data, 1, rasch = TRUE, SE=TRUE)
    expect_is(modm1, 'mirtClass')          
    expect_is(modm2, 'mirtClass')
    expect_is(modm3, 'mirtClass')
    
    fm1 <- fscores(modm1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    fm2 <- fscores(modm2, verbose = FALSE)
    expect_is(fm2, 'matrix')
    fm3 <- fscores(modm3, verbose = FALSE)
    expect_is(fm3, 'matrix')
})

test_that('poly', {
    modp1 <- mirt(Science, 1, SE=TRUE)
    modp2 <- mirt(Science, 2)
    modp3 <- mirt(Science, 1, rasch = TRUE, SE=TRUE)
    expect_is(modp1, 'mirtClass')          
    expect_is(modp2, 'mirtClass')
    expect_is(modp3, 'mirtClass')
    
    fm1 <- fscores(modp1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    fm2 <- fscores(modp2, rotate = 'oblimin', verbose = FALSE)
    expect_is(fm2, 'matrix')
    fm3 <- fscores(modp3, verbose = FALSE)
    expect_is(fm3, 'matrix')
})

