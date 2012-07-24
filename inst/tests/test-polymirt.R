context('polymirt')

test_that('dich data', {
    data <- expand.table(LSAT7)        
    modm1 <- polymirt(data, 1, verbose = FALSE)
    modm2 <- polymirt(data, 2, verbose = FALSE)
    expect_is(modm1, 'polymirtClass')          
    expect_is(modm1, 'polymirtClass')
    fs <- fscores(modm2, rotate = 'promax')
    expect_is(fs, 'matrix')
})

test_that('poly data', {
    modp1 <- polymirt(Science, 1, verbose = FALSE)
    modp2 <- polymirt(Science, 2, verbose = FALSE)
    expect_is(modp1, 'polymirtClass')          
    expect_is(modp2, 'polymirtClass')
    fs <- fscores(modp2, rotate = 'oblimin')
    expect_is(fs, 'matrix')
})