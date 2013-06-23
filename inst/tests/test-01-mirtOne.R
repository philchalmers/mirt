context('mirtOne')

test_that('dich', {
    data <- expand.table(LSAT7)        
    mod1 <- mirt(data, 1, verbose=FALSE)
    expect_is(mod1, 'ConfirmatoryClass')          
    expect_equal(mod1@df, 21)
    cfs <- as.numeric(do.call(c, coef(mod1)))
    expect_equal(cfs, c(0.988, 1.856, 0.000, 1.000, 1.081, 0.808, 0.000, 1.000, 1.706, 1.804, 0.000,
                        1.000, 0.765, 0.486, 0.000, 1.000, 0.736, 1.855, 0.000, 1.000, 0.000, 1.000),
                 tollerance = 1e-3)
    sv <- mod2values(mod1)
    sv$est <- FALSE
    moddummy <- mirt(data, 1, pars= sv, verbose=FALSE)
    expect_is(moddummy, 'ConfirmatoryClass')
    sv2 <- mod2values(moddummy)
    expect_equal(sv$value, sv2$value)
    modm1 <- mirt(data, 1, SE = TRUE, SE.type = 'SEM', verbose=FALSE)
    cfs <- as.numeric(do.call(c, coef(modm1)))
    expect_equal(cfs, c(0.988, 0.174, 1.856, 0.131, 0.000, NA, 1.000, NA, 1.081, 0.169, 
                        0.808, 0.087, 0.000, NA, 1.000, NA, 1.706, 0.297, 1.804, 0.199, 
                        0.000, NA, 1.000, NA, 0.765, 0.136, 0.486, 0.077, 0.000, NA, 
                        1.000, NA, 0.736, 0.149, 1.855, 0.112, 0.000, NA, 1.000, NA, 0.000, 
                        NA, 1.000, NA), tollerance = 1e-3)    
    expect_is(modm1, 'ConfirmatoryClass')          
    modm2 <- mirt(data, 1, SE = TRUE, SE.type = 'BL', verbose=FALSE)
    cfs <- as.numeric(do.call(c, coef(modm2)))
    expect_equal(cfs, c(0.988, 0.177, 1.856, 0.131, 0.000, NA, 1.000, NA, 1.081, 0.169, 0.808, 
                        0.091, 0.000, NA, 1.000, NA, 1.706, 0.320, 1.804, 0.204, 0.000,
                        NA, 1.000, NA, 0.765, 0.134, 0.486, 0.075, 0.000, NA, 1.000, NA, 0.736,
                        0.151, 1.855, 0.114, 0.000, NA, 1.000, NA, 0.000, NA, 1.000, NA), 
                 tollerance = 1e-3)
    expect_is(modm2, 'ConfirmatoryClass')          
    suppressWarnings(modm2 <- mirt(data, 1, SE = TRUE, SE.type = 'MHRM', verbose=FALSE))
    expect_is(modm2, 'ConfirmatoryClass')
    modm3 <- mirt(data, 1, itemtype = 'Rasch', verbose=FALSE)
    expect_is(modm3, 'ConfirmatoryClass')
    expect_equal(modm3@df, 25)
    modm3 <- mirt(data, 1, itemtype = 'Rasch', SE = TRUE, technical=list(TOL=1e-6), verbose=FALSE)
    expect_is(modm3, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modm3)))
    expect_equal(cfs, c(1.000, NA, 1.868, 0.098, 0.000, NA, 1.000, NA, 1.000, NA, 0.791, 0.080,
                        0.000, NA, 1.000, NA, 1.000, NA, 1.461, 0.094, 0.000, NA, 1.000, NA, 1.000,
                        NA, 0.522, 0.076, 0.000, NA, 1.000, NA, 1.000, NA, 1.993, 0.100, 0.000,
                        NA, 1.000, NA, 0.000, NA, 1.023, 0.344), 
                 tollerance = 1e-3)
    modm4 <- mirt(data, 1, itemtype = '1PL', verbose=FALSE)    
    expect_is(modm4, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modm4)))
    expect_equal(cfs, c(1.014, 1.868, 0.000, 1.000, 1.014, 0.791, 0.000, 1.000, 1.014, 1.461,
                        0.000, 1.000, 1.014, 0.522, 0.000, 1.000, 1.014, 1.993, 0.000, 1.000, 
                        0.000, 0.995), tollerance = 1e-3)
    svalues <- mirt(data, 1, pars = 'values', verbose=FALSE)
    svalues[22, 5] <- 2
    modm5 <- mirt(data, 1, pars = svalues, verbose=FALSE)    
    expect_is(modm5, 'ConfirmatoryClass')
    data[1,1] <- data[2,2] <- NA
    modm6 <- mirt(data, 1, verbose=FALSE)
    expect_equal(modm6@df, 23)
    expect_is(modm6, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modm6)))
    expect_equal(cfs, c(0.969, 1.851, 0.000, 1.000, 1.074, 0.808, 0.000, 1.000, 1.717, 1.811,
                        0.000, 1.000, 0.763, 0.486, 0.000, 1.000, 0.731, 1.852, 0.000, 1.000,
                        0.000, 1.000), tollerance = 1e-3)
    
    fm1 <- fscores(modm1, verbose = FALSE)
    expect_is(fm1, 'matrix')      
    expect_true(mirt:::closeEnough(fm1[1:6,'F1'] - c(-1.8665957, -1.5266920, -1.5134024,
                                                     -1.1852276, -1.0946830, -0.7666992), -1e-4, 1e-4))
    fm2 <- fscores(modm2, method = 'MAP', verbose = FALSE)
    expect_is(fm2, 'matrix')
    expect_true(mirt:::closeEnough(fm2[1:6,'F1'] - c(-1.8165552, -1.4946906, -1.4822982, 
                                                     -1.1789899, -1.0958928, -0.7951026), -1e-4, 1e-4))
    fm3 <- fscores(modm3, method = 'ML', full.scores = TRUE, verbose = FALSE)
    expect_is(fm3, 'matrix')
    expect_true(fm3[1, 'F1'] == -Inf && fm3[1000, 'F1'] == Inf)
    expect_true(mirt:::closeEnough(as.numeric(fm3[c(13,34,40),'F1'])
                                   - c(-2.812972, -1.769511, -2.812972), -1e-4, 1e-4))
    fm3 <- fscores(modm3, method = 'ML', full.scores = TRUE, verbose = FALSE, scores.only=TRUE)
    expect_is(fm3, 'matrix')
    fm4 <- fscores(modm6, method = 'ML', full.scores = TRUE, verbose = FALSE)
    expect_is(fm4, 'matrix')
    fm5 <- fscores(modm6, method = 'ML', full.scores = FALSE, verbose = FALSE)
    expect_is(fm5, 'matrix')
    fm6 <- fscores(modm1, method = 'EAPsum', full.scores = FALSE, verbose = FALSE)
    expect_is(fm6, 'data.frame')
    expect_true(mirt:::closeEnough(as.numeric(as.matrix(fm6)) -  
                 c(0.0000000,  1.0000000,  2.0000000,  3.0000000,  4.0000000,  5.0000000,
                   -1.8665957, -1.4314464, -0.9487476, -0.4131919,  0.1516851, 0.7269940,
                   0.6872827,  0.6831615,  0.6941894,  0.7210850,  0.7587511,  0.8004654), -1e-4, 1e-4))
    
    res1 <- residuals(modm1, verbose = FALSE)
    expect_equal(as.numeric(res1), c(NA, -0.452, -0.853,  2.574,  2.392,  0.021,     NA,  1.062,
                                     -0.267, -1.383,  0.029,  0.033,     NA, -0.153, -0.003,  0.051,
                                     0.016,  0.012, NA,  0.000,  0.049,  0.037,  0.002,  0.000, NA), 
                 tollerance = 1e-3)
    res2 <- residuals(modm2, verbose = FALSE)
    expect_is(res1, 'matrix')
    expect_is(res2, 'matrix')
    IP1 <- itemplot(modm1, 1)
    IP2 <- itemplot(modm2, 1)
    expect_is(IP1, 'trellis')
    expect_is(IP2, 'trellis')
    TP1 <- plot(modm1)
    TP2 <- plot(modm2)
    expect_is(TP1, 'trellis')    
    expect_is(TP2, 'trellis')
    ifit <- itemfit(modm1, X2 = TRUE)    
    expect_is(ifit, 'data.frame')
    expect_true(mirt:::closeEnough(as.numeric(ifit$Zh) - c(1.431838, 6.354917, 5.310844, 5.804449, 
                                                           0.696139), -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(as.numeric(ifit$X2) - c(15.691499, 39.888656, 23.843572, 
                                                           69.410816,  9.509994), -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(as.numeric(ifit$S_X2) - c(4.749440, 14.451071,  1.270381,
                                                             5.237400,  0.941125), -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(as.numeric(ifit$df) - c(5,5,5,6,5), -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(as.numeric(ifit$df.S_X2) - c(2,2,2,2,2), -1e-4, 1e-4))
    
    fitm1 <- fitIndices(modm1)
    expect_is(fitm1, 'list')
    expect_true(mirt:::closeEnough(fitm1$M2 - 11.45125, -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(fitm1$df.M2 - 5, -1e-4, 1e-4))
    fitm2 <- fitIndices(modm3)
    expect_is(fitm2, 'list')
    expect_true(mirt:::closeEnough(fitm2$M2 - 22.57281, -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(fitm2$df.M2 - 9, -1e-4, 1e-4))
    
    data <- expand.table(LSAT7)
    model <- confmirt.model('F1 = 1-3
        F2 = 3-5', quiet = TRUE)
    modm1 <- mirt(data, model, verbose=FALSE)
    expect_equal(modm1@df, 20)
    modm2 <- suppressMessages(mirt(data, model, itemtype=c('2PL','2PL', 'PC2PL','2PL', '2PL'), verbose=FALSE))
    expect_equal(modm2@df, 19)
    modm3 <- mirt(data, model, SE = TRUE, verbose=FALSE)
    expect_is(modm3, 'ConfirmatoryClass')
    
    fm1 <- fscores(modm1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    fm2 <- fscores(modm2, method = 'MAP', verbose = FALSE)
    expect_is(fm2, 'matrix')    
})

