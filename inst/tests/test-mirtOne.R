context('mirtOne')

test_that('dich', {
    data <- expand.table(LSAT7)        
    modm1 <- mirt(data, 1)
    expect_is(modm1, 'ConfirmatoryClass')          
    modm2 <- mirt(data, 2, SE = TRUE)
    expect_is(modm2, 'ExploratoryClass')
    modm3 <- mirt(data, 1, itemtype = 'Rasch', SE = FALSE)
    expect_is(modm3, 'ConfirmatoryClass')
    modm4 <- mirt(data, 1, itemtype = '1PL')    
    expect_is(modm4, 'ConfirmatoryClass')
    svalues <- mirt(data, 1, pars = 'values')
    svalues[22, 5] <- 2
    modm5 <- mirt(data, 1, pars = svalues)    
    expect_is(modm5, 'ConfirmatoryClass')
    data[1,1] <- data[2,2] <- NA
    modm6 <- mirt(data, 1)
    expect_is(modm6, 'ConfirmatoryClass')
    
    fm1 <- fscores(modm1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    fm2 <- fscores(modm2, method = 'MAP', verbose = FALSE)
    expect_is(fm2, 'matrix')
    fm3 <- fscores(modm3, method = 'ML', full.scores = TRUE, verbose = FALSE)
    expect_is(fm3, 'matrix')
    fm4 <- fscores(modm6, method = 'ML', full.scores = TRUE, verbose = FALSE)
    expect_is(fm4, 'matrix')
    fm5 <- fscores(modm6, method = 'ML', full.scores = FALSE, verbose = FALSE)
    expect_is(fm5, 'matrix')
    
    res1 <- residuals(modm1, verbose = FALSE)
    res2 <- residuals(modm2, verbose = FALSE)
    expect_is(res1, 'matrix')
    expect_is(res2, 'matrix')
    IP1 <- itemplot(modm1, 1)
    IP2 <- itemplot(modm2, 1)
    expect_is(IP1, 'NULL')
    expect_is(IP2, 'trellis')
    TP1 <- plot(modm1)
    TP2 <- plot(modm2)
    expect_is(TP1, 'trellis')    
    expect_is(TP2, 'trellis')
    ifit <- itemfit(modm1, X2 = TRUE)
    expect_is(ifit, 'data.frame')
})



