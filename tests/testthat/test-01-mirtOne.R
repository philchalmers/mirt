expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('dich', {
    data <- expand.table(LSAT7)
    mod1 <- mirt(data, 1, verbose=FALSE)
    expect_class(mod1, 'SingleGroupClass')
    expect_equal(extract.mirt(mod1, 'df'), 21)
    cfs <- as.numeric(do.call(c, coef(mod1)))
    expect_equal(cfs, c(0.988, 1.8561, 0, 1, 1.081, 0.808, 0, 1, 1.706, 1.8043, 0, 1, 0.7651, 0.486, 0, 1, 0.7358, 1.8545, 0, 1, 0, 1),
                 tolerance = 1e-2)
    sv <- mod2values(mod1)
    sv$est <- FALSE
    moddummy <- mirt(data, 1, pars= sv, verbose=FALSE)
    expect_class(moddummy, 'SingleGroupClass')
    sv2 <- mod2values(moddummy)
    expect_equal(sv$value, sv2$value)
    modm1 <- mirt(data, 1, SE = TRUE, SE.type = 'SEM', verbose=FALSE)
    expect_true(extract.mirt(modm1, 'SEMconv'))
    cfs <- as.numeric(do.call(c, coef(modm1)))
    expect_equal(extract.mirt(modm1, 'condnum'), 30.12751, tolerance = 1e-4)
    expect_equal(cfs, c(0.9876, 0.6367, 1.3384, 1.8559, 1.5978, 2.1139, 0, NA, NA, 1, NA, NA, 1.0808, 0.7604, 1.4013, 0.808, 0.6335, 0.9825, 0, NA, NA, 1, NA, NA, 1.7075, 1.0868, 2.3281, 1.8052, 1.4028, 2.2076, 0, NA, NA, 1, NA, NA, 0.765, 0.5065, 1.0235, 0.486, 0.3114, 0.6606, 0, NA, NA, 1, NA, NA, 0.7357, 0.4246, 1.0467, 1.8545, 1.6332, 2.0757, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)
    expect_class(modm1, 'SingleGroupClass')
    modm2 <- mirt(data, 1, SE = TRUE, SE.type = 'Richardson', verbose=FALSE)
    cfs <- as.numeric(do.call(c, coef(modm2)))
    expect_equal(extract.mirt(modm2, 'condnum'), 30.24068, tolerance = 1e-3)
    expect_equal(cfs, c(0.988, 0.6406, 1.3354, 1.8561, 1.5984, 2.1138, 0, NA, NA, 1, NA, NA, 1.081, 0.7501, 1.4119, 0.808, 0.6291, 0.9869, 0, NA, NA, 1, NA, NA, 1.706, 1.0779, 2.334, 1.8043, 1.4036, 2.205, 0, NA, NA, 1, NA, NA, 0.7651, 0.5022, 1.028, 0.486, 0.3392, 0.6328, 0, NA, NA, 1, NA, NA, 0.7358, 0.4395, 1.032, 1.8545, 1.6302, 2.0787, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)
    expect_class(modm2, 'SingleGroupClass')
    modm2b <- mirt(data, 1, SE = TRUE, SE.type = 'Fisher', verbose=FALSE)
    expect_equal(extract.mirt(modm2b, 'condnum'), 29.0323, tolerance = 1e-3)
    modm3 <- mirt(data, 1, itemtype = 'Rasch', verbose=FALSE, SE=TRUE)
    expect_class(modm3, 'SingleGroupClass')
    expect_equal(extract.mirt(modm3, 'df'), 25)
    expect_equal(extract.mirt(modm3, 'condnum'), 4.488772, tolerance = 1e-4)
    LG <- lagrange(modm3, parnum = list(1, 5))
    expect_equal(LG$X2, c(0.37024340, 0.05429718), tolerance = 1e-4)
    LG2 <- lagrange(modm3, parnum = list(c(1, 5)))
    expect_equal(LG2$X2, 0.4816444, tolerance = 1e-4)
    dat <- expand.table(LSAT6)
    modm3 <- mirt(dat, 1, itemtype = 'Rasch', SE = TRUE, SE.type = 'SEM', verbose=FALSE)
    expect_class(modm3, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(modm3)))
    expect_equal(cfs, c(1,NA,NA,2.73,2.478,2.983,0,NA,NA,1,NA,NA,1,NA,NA,0.999,0.845,1.152,0,NA,NA,1,NA,NA,1,NA,NA,0.24,0.1,0.38,0,NA,NA,1,NA,NA,1,NA,NA,1.306,1.143,1.47,0,NA,NA,1,NA,NA,1,NA,NA,2.099,1.896,2.303,0,NA,NA,1,NA,NA,0,NA,NA,0.57,0.371,0.77),
                 tolerance = 1e-2)
    modm3b <- mirt(dat, 'F = 1-5
                   CONSTRAIN = (1-5, a1)', verbose=FALSE, SE=F)
    fitm2 <- M2(modm3b)
    expect_true(mirt:::closeEnough(fitm2$M2 - 5.292566, -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(fitm2$df.M2 - 9, -1e-4, 1e-4))
    model <- mirt.model('F = 1-5
                        CONSTRAIN = (1-5, a1)', quiet=TRUE)
    modm4 <- mirt(data, model, verbose = FALSE, SE=T, SE.type = 'crossprod')
    expect_equal(extract.mirt(modm4, 'condnum'), 5.171716, tolerance = 1e-4)
    cfs <- as.numeric(do.call(c, coef(modm4)))
    expect_equal(cfs, c(1.011, 0.885, 1.138, 1.868, 1.67, 2.067, 0, NA, NA, 1, NA, NA, 1.011, 0.885, 1.138, 0.791, 0.631, 0.951, 0, NA, NA, 1, NA, NA, 1.011, 0.885, 1.138, 1.461, 1.277, 1.644, 0, NA, NA, 1, NA, NA, 1.011, 0.885, 1.138, 0.521, 0.367, 0.676, 0, NA, NA, 1, NA, NA, 1.011, 0.885, 1.138, 1.993, 1.793, 2.193, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)
    svalues <- mirt(data, 1, pars = 'values', verbose=FALSE)
    svalues[22, 'value'] <- 2
    modm5 <- mirt(data, 1, pars = svalues, verbose=FALSE)
    expect_class(modm5, 'SingleGroupClass')
    expect_warning(modm7 <- mirt(data, 1, '4PL', verbose=FALSE, parprior = list(c(3,7,11,15,19,'norm', -1.7, .1),
                                                                 c(4,8,12,16,20,'norm', 1.7, .1))),
                   "EM cycles terminated after 500 iterations.")
    expect_equal(extract.mirt(modm7, 'df'), 11)
    expect_class(modm7, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(modm7)))
    expect_equal(cfs, c(5.12,8.557,0.154,0.859,5.741,3.595,0.16,0.843,10.215,10.937,0.155,0.861,1.284,0.854,0.153,0.845,4.629,8.905,0.154,0.859,0,1), tolerance = 1e-2)
    data[1,1] <- data[2,2] <- NA
    modm6 <- mirt(data, 1, verbose=FALSE)
    expect_equal(modm6@Fit$df, 21)
    expect_class(modm6, 'SingleGroupClass')
    cfs <- as.numeric(do.call(c, coef(modm6)))
    expect_equal(cfs, c(0.969, 1.851, 0, 1, 1.074, 0.808, 0, 1, 1.717, 1.811, 0, 1, 0.763, 0.486, 0, 1, 0.731, 1.852, 0, 1, 0, 1), tolerance = 1e-2)
    modideal <- mirt(data, 1, verbose=FALSE, itemtype='ideal')
    cfs <- as.numeric(do.call(c, coef(modideal)))
    expect_equal(cfs, c(0.2833761,-0.5685226,0.4186987,-0.8906244,0.5679343,-0.5649427,0.2920432,-0.999962,0.2066424,-0.5591496,0,1), tolerance = 1e-2)
    modspline <- mirt(data, 1, verbose=FALSE, itemtype=c(rep('2PL', 4), 'spline'))
    cfs <- as.numeric(do.call(c, coef(modspline)))
    expect_equal(cfs, c(0.963,1.848,0,1,1.08,0.809,0,1,1.719,1.812,0,1,0.758,0.485,0,1,-0.91,0.359,0.8,11.81,0,1), tolerance = 1e-2)
    expect_equal(logLik(modspline), -2657.558, tolerance = 1e-4)

    #QMCEM
    mod <- mirt(dat, 1, method = 'QMCEM', verbose=FALSE, optimizer='NR')
    expect_equal(extract.mirt(mod, 'logLik'), -2466.653, tolerance=1e-4)
    fs <- fscores(mod, QMC=TRUE, verbose=FALSE, full.scores=FALSE)
    expect_equal(fs[1:3,'F1'], c(-1.887714, -1.473829, -1.453797), tolerance=1e-4)
    m2 <- M2(mod, QMC=TRUE)
    expect_equal(m2$M2, 4.737141, tolerance=1e-5)
    ifit <- itemfit(mod, QMC=TRUE)
    expect_equal(ifit$p.S_X2[1], .7984, tolerance = 1e-2)
    rfit <- residuals(mod, QMC=TRUE, verbose=FALSE)
    expect_equal(as.numeric(rfit[,1]), c(NA, 0.0468866, 0.3906001, 0.2476980, 0.5195561), tolerance = 1e-2)

    fm1 <- fscores(modm1, verbose = FALSE, full.scores=FALSE)
    expect_class(fm1, 'matrix')
    expect_true(mirt:::closeEnough(fm1[1:6,'F1'] - c(-1.8665957, -1.5266920, -1.5134024,
                                                     -1.1852276, -1.0946830, -0.7666992), -1e-2, 1e-2))
    fm2 <- fscores(modm2, method = 'MAP', verbose = FALSE, full.scores=FALSE)
    expect_class(fm2, 'matrix')
    expect_true(mirt:::closeEnough(fm2[1:6,'F1'] - c(-1.8165552, -1.4946906, -1.4822982,
                                                     -1.1789899, -1.0958928, -0.7951026), -1e-2, 1e-2))
    fm3 <- fscores(modm4, method = 'ML', full.scores = TRUE, verbose = FALSE)
    expect_class(fm3, 'matrix')
    expect_true(fm3[1, 'F'] == -Inf && fm3[1000, 'F'] == Inf)
    expect_true(mirt:::closeEnough(as.numeric(fm3[c(13,34,40),'F'])
                                   - c(-2.783489, -1.750890, -2.783489), -1e-2, 1e-2))
    fm3 <- fscores(modm3, method = 'ML', full.scores = TRUE, verbose = FALSE, scores.only=TRUE)
    expect_class(fm3, 'matrix')
    fm4 <- fscores(modm6, method = 'ML', full.scores = TRUE, verbose = FALSE)
    expect_class(fm4, 'matrix')
    fm5 <- fscores(modm6, method = 'ML', full.scores = FALSE, verbose = FALSE)
    expect_class(fm5, 'matrix')
    fm6 <- fscores(modm1, method = 'EAPsum', full.scores = FALSE, verbose = FALSE)
    expect_class(fm6, 'data.frame')
    expect_true(mirt:::closeEnough(as.numeric(as.matrix(fm6)) - c(0,1,2,3,4,5,-1.86979,-1.431861,-0.9488463,-0.4131963,0.1517289,0.7271877,0.6927032,0.6838697,0.6942298,0.7210951,0.758772,0.8009335,12,40,114,205,321,308,10.08994,44.65882,109.773,207.7391,319.1854,308.5536,0.6013157,0.6971433,0.4034422,0.1900446,0.101567,0.0315184), -1e-2, 1e-2))
    expect_equal(as.numeric(attr(fm6, 'fit')['rxx_F1']), 0.4319948, tolerance = 1e-4)

    res1 <- residuals(modm1, verbose = FALSE)
    expect_equal(as.numeric(res1), c(NA,0.451213,0.8562096,2.577395,2.392183,-0.02124177,NA,1.053826,0.2662122,1.383089,-0.02926106,0.03246269,NA,0.1542321,0.002940504,0.05076805,-0.01631601,-0.01241902,NA,9.962506e-06,0.04890994,-0.0371899,-0.00171479,9.981236e-05,NA),
                 tolerance = 1e-2)
    res2 <- residuals(modm2, verbose = FALSE)
    expect_class(res1, 'matrix')
    expect_class(res2, 'matrix')
    IP1 <- itemplot(modm1, 1)
    IP2 <- itemplot(modm2, 1)
    expect_class(IP1, 'trellis')
    expect_class(IP2, 'trellis')
    TP1 <- plot(modm1)
    TP2 <- plot(modm2)
    expect_class(TP1, 'trellis')
    expect_class(TP2, 'trellis')
    ifit <- itemfit(modm1, c('S_X2', 'X2', 'Zh'))
    expect_class(ifit, 'data.frame')
    expect_true(mirt:::closeEnough(as.numeric(ifit$Zh) - c(1.431838, 6.354917, 5.310844, 5.804449,
                                                           0.696139), -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(as.numeric(ifit$X2) - c(91.71819, 390.07985, 145.39978, 329.48529, 129.49679), -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(as.numeric(ifit$S_X2) - c(4.749440, 14.451071,  1.270381,
                                                             5.237400,  0.941125), -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(as.numeric(ifit$df) - c(8,8,8,8,8), -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(as.numeric(ifit$df.S_X2) - c(2,2,2,2,2), -1e-4, 1e-4))

    fitm1 <- M2(modm1)
    expect_class(fitm1, 'data.frame')
    expect_true(mirt:::closeEnough(fitm1$M2 - 11.93841, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(fitm1$df.M2 - 5, -1e-4, 1e-4))
    fitm2 <- M2(modm3)
    expect_class(fitm2, 'data.frame')
    expect_true(mirt:::closeEnough(fitm2$M2 - 5.291576, -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(fitm2$df.M2 - 9, -1e-4, 1e-4))

    data <- expand.table(LSAT7)
    model <- 'F1 = 1-3
              F2 = 3-5'
    modm1 <- mirt(data, model, verbose=FALSE)
    expect_equal(extract.mirt(modm1, 'df'), 20)
    modm2 <- mirt(data, model, itemtype=c('2PL','2PL', 'PC2PL','2PL', '2PL'), TOL=1e-3, verbose=FALSE)
    cfs <- as.numeric(do.call(c, coef(modm2)))
    expect_equal(cfs, c(0.6514, 0, 1.7031, 0, 1, 1.4872, 0, 0.9174, 0, 1, 2.5151, 3.4949, 3.1337, 5.2577, 0, 1, 0, 0.7058, 0.4789, 0, 1, 0, 0.8524, 1.9092, 0, 1, 0, 0, 1, 0, 1),
                 tolerance = 1e-2)
    expect_equal(extract.mirt(modm2, 'df'), 19)
    modm3 <- mirt(data, model, SE = TRUE, verbose=FALSE)
    expect_class(modm3, 'SingleGroupClass')

    fm1 <- fscores(modm1, verbose = FALSE, full.scores=FALSE)
    expect_class(fm1, 'matrix')
    fm2 <- fscores(modm3, method = 'MAP', verbose = FALSE)
    expect_class(fm2, 'matrix')

    data[1,1] <- NA
    modm1 <- mirt(data, 1, verbose=FALSE)
    out1 <- M2(modm1, na.rm=TRUE)
    out2 <- itemfit(modm1, na.rm=TRUE)
    out3 <- fscores(modm1, na.rm=TRUE, method = 'EAPsum', full.scores=FALSE, verbose = FALSE)
    expect_equal(out1$M2, 11.76977, tolerance=1e-4)
    expect_equal(out2$S_X2[1], 4.8448539, tolerance=1e-4)
    expect_equal(out3$expected[1], 9.931098, tolerance=1e-4)

    # missing data
    set.seed(1234)
    pick <- sample(1:5000, 500)
    dat <- as.matrix(expand.table(LSAT7))
    dat[pick] <- NA
    mod <- mirt(dat, itemtype='Rasch', verbose=FALSE)
    syntax <- "F = 1-5
               START = (1, d, 0.0), (1, a1, 1.0)
               FIXED = (1, d), (1, a1)
               FREE = (GROUP, MEAN_1)"
    mod2 <-  mirt(dat, syntax,
                  itemtype='Rasch', verbose=FALSE)
    syntax <- 'F = 1-5
               CONSTRAIN = (1-5, a1)'
    mod3 <- mirt(dat, syntax, verbose=FALSE)
    expect_equal(2*logLik(mod) - logLik(mod2) - logLik(mod3),
                 0, 1e-2)
})

