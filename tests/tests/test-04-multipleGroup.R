context('multipleGroup')

test_that('one factor', {
    set.seed(12345)
    a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
    d <- matrix(rnorm(15,0,.7),ncol=1)
    itemtype <- rep('dich', nrow(a))
    N <- 1000
    dataset1 <- simdata(a, d, N, itemtype)
    dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
    dat <- rbind(dataset1, dataset2)
    group <- c(rep('D1', N), rep('D2', N))
    MGmodel1 <- 'F1 = 1-15'
    models <- mirt.model(MGmodel1, quiet = TRUE)

    mod_Rasch <- multipleGroup(dat, models, itemtype = 'Rasch', SE=TRUE, SE.type = 'crossprod',
                               group = group, verbose = FALSE, method = 'EM')
    cfs <- as.numeric(na.omit(do.call(rbind, coef(mod_Rasch, printSE=TRUE, as.data.frame=TRUE))))
    expect_equal(cfs, c(0.51134,-0.65255,-0.19084,0.86483,0.12324,0.76332,0.91934,-0.33708,-1.09074,-1.1665,1.24226,-0.2151,0.40659,0.44631,-0.06508,0.97838,0.66571,-0.49219,-0.13083,0.8908,0.3791,0.9721,1.03114,-0.33142,-0.90574,-1.22668,1.55031,-0.12033,0.37378,0.52922,-0.06789,1.55749,0.07894,0.07929,0.07721,0.08119,0.07786,0.07987,0.08224,0.07738,0.08299,0.08428,0.08508,0.0777,0.07804,0.07845,0.07726,0.06713,0.08612,0.08522,0.08467,0.08723,0.08485,0.08794,0.08911,0.08454,0.08803,0.09063,0.09475,0.08545,0.08438,0.08586,0.08411,0.10285),
                 tolerance = 1e-3)
    expect_equal(logLik(mod_Rasch), -17944.17, tolerance = 1e-4)
    EAP <- fscores(mod_Rasch, full.scores=TRUE)
    expect_equal(cor(EAP, rowSums(dat))[1], .99, tolerance = 1e-2)
    pf <- personfit(mod_Rasch, Theta=EAP)
    pffit <- c(as.numeric(as.matrix(head(pf))), as.numeric(as.matrix(tail(pf))))
    expect_equal(pffit, c(0.9509209,1.142614,0.598627,0.4390729,0.8847659,0.6006927,-0.1580301,0.4766065,-1.366012,-1.688975,-0.1716656,-1.35694,0.9657947,0.9039366,0.6798214,0.5300523,0.8204743,0.6796139,-0.1220695,-0.2348978,-1.36503,-1.825007,-0.5448904,-1.366104,0.1942602,0.1375231,1.282292,1.543588,0.5053905,1.280045,1.206,0.9890109,0.7871311,0.8204198,1.415266,1.078969,0.6701176,0.1288671,-0.5428167,-0.7281717,1.408353,0.3263531,0.8778794,1.013104,0.9167278,0.8353493,1.222656,0.7957533,-0.4393401,0.1491026,-0.268665,-0.8304078,1.080016,-0.3983634,0.1985518,-0.0408503,0.4460764,0.817872,-1.193991,0.3035265),
                 tolerance = 1e-3)
    # mod_QMCEM <- multipleGroup(dat, models, group=group, method = 'QMCEM', verbose=FALSE,
    #                            optimizer='NR')
    # expect_equal(extract.mirt(mod_QMCEM, 'logLik'), -17849.64, tolerance=1e-2)
    mod_configural <- multipleGroup(dat, models, SE=TRUE, SE.type = 'crossprod', optimizer='NR',
                                    group = group, verbose = FALSE, method = 'EM')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural)[[1L]]))
    cfs <- as.numeric(na.omit(cfs[cfs != 0 & cfs != 1]))
    expect_equal(cfs, c(1.0706,0.8462,1.295,0.524,0.3627,0.6852,1.217,0.9668,1.4672,-0.6995,-0.8745,-0.5244,0.9488,0.7445,1.1532,-0.1877,-0.3385,-0.037,0.8975,0.6948,1.1002,0.8423,0.681,1.0036,1.097,0.8781,1.316,0.1274,-0.0294,0.2841,0.5648,0.3916,0.7381,0.6824,0.5383,0.8266,1.2742,1.0146,1.5337,1.0016,0.8129,1.1904,0.9249,0.7242,1.1255,-0.3298,-0.4801,-0.1794,0.8903,0.6815,1.0992,-1.0591,-1.2284,-0.8897,0.723,0.5272,0.9188,-1.0828,-1.2463,-0.9192,0.8303,0.6252,1.0354,1.1879,1.0155,1.3604,1.4757,1.1986,1.7528,-0.252,-0.429,-0.0749,1.2905,1.0339,1.5471,0.4451,0.2742,0.616,1.0348,0.8239,1.2457,0.4526,0.2955,0.6096,0.8637,0.665,1.0624,-0.0617,-0.2093,0.0859),
                 tolerance = 1e-2)
    expect_equal(extract.mirt(mod_configural, 'df'), 32707)
    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), verbose = FALSE,
                                method = 'EM', optimizer = 'NR')
    expect_is(mod_metric, 'MultipleGroupClass')
    expect_equal(extract.mirt(mod_metric, 'df'), 32722)
    mod_scalar2 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_var','free_means'))
    cfs <- as.numeric(do.call(c, coef(mod_scalar2)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.1236,0.5623,1.1966,-0.6734,1.0347,-0.2125,0.9671,0.8262,1.1147,0.2057,0.4526,0.7108,1.1532,0.978,0.9278,-0.3708,0.9542,-1.0335,0.7012,-1.1231,0.871,1.296,1.5209,-0.2726,1.1443,0.3605,1.0558,0.4491,0.8774,-0.105),
                 tolerance = 1e-2)
    expect_is(mod_scalar2, 'MultipleGroupClass')
    expect_equal(extract.mirt(mod_scalar2, 'df'), 32735)
    newmodel <- mirt.model('F = 1-15
                            CONSTRAINB = (1-15, a1), (1,2,3-15,d)')
    mod_scalar1 <- multipleGroup(dat, newmodel, group = group, verbose = FALSE, invariance='free_var')
    expect_is(mod_scalar1, 'MultipleGroupClass')
    mod_EH <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                            dentype="empiricalhist", optimizer = 'NR')
    expect_is(mod_EH, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_EH)[[1L]]))
    expect_equal(cfs, c(1.0035,0.5413,0,1,1.186,-0.6822,0,1,0.8913,-0.1695,0,1,0.8556,0.8614,0,1,1.0557,0.1476,0,1,0.5402,0.694,0,1,1.2011,1.0217,0,1,0.8957,-0.3138,0,1,0.8306,-1.0388,0,1,0.6773,-1.072,0,1,0.7762,1.2018,0,1,1.4008,-0.2289,0,1,1.2235,0.4661,0,1,0.9764,0.4707,0,1,0.8229,-0.0455,0,1,0,1),
                 tolerance = 1e-2)

    dat[1,1] <- dat[2,2] <- NA
    mod_missing <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_var'))
    expect_is(mod_missing, 'MultipleGroupClass')
    expect_equal(extract.mirt(mod_missing, 'df'), 32736)

    fs1 <- fscores(mod_metric, verbose = FALSE, full.scores=FALSE)
    expect_true(mirt:::closeEnough(fs1[[1]][1:6, 'F1'] - c(-2.0826, -1.6822, -1.3988, -1.5287, -1.7450, -1.3712), -1e-2, 1e-2))
    fs2 <- fscores(mod_metric, full.scores = TRUE, full.scores.SE=TRUE, method = 'ML')
    expect_equal(as.numeric(head(fs2)), c(0.4035981,1.405389,0.7882771,1.259476,1.253976,1.027287,0.4782066,0.6323413,0.5177877,0.5994493,0.5982813,0.5545867),
                tolerance = 1e-2)
    fs3 <- fscores(mod_missing, verbose = FALSE, full.scores=FALSE)
    fs4 <- fscores(mod_missing, full.scores = TRUE)
    fs5 <- fscores(mod_metric, full.scores = TRUE, scores.only=TRUE)
    expect_is(fs1, 'list')
    expect_is(fs2, 'matrix')
    expect_is(fs3, 'list')
    expect_is(fs4, 'matrix')

    fit1 <- M2(mod_metric)
    expect_is(fit1, 'data.frame')
    expect_true(mirt:::closeEnough(fit1$M2 - c(213.1749), -1e-2, 1e-2))
    expect_equal(fit1$SRMSR.D1, 0.04152426, tolerance = 1e-4)
    expect_equal(fit1$TLI, .99758, tolerance = 1e-4)
    expect_true(mirt:::closeEnough(fit1$df - 195, -1e-4, 1e-4))
    fit2 <- itemfit(mod_metric, c('S_X2', 'Zh'))
    expect_is(fit2, 'list')
    expect_equal(as.numeric(fit2[[1]][1L,]), c(1.000000, 2.6646153, 8.1727058, 11.000000, 0.6977546),
                 tolerance = 1e-4)
    fit3 <- M2(mod_scalar2)
    expect_true(mirt:::closeEnough(fit3$M2 - c(198.6178), -1e-4, 1e-4))
    expect_equal(fit3$SRMSR.D1, 0.026854, tolerance = 1e-4)
    expect_equal(fit3$TLI, 1.001169, tolerance = 1e-4)
    expect_true(mirt:::closeEnough(fit3$df - 208, -1e-4, 1e-4))

    g1 <- extract.group(mod_metric, 1)
    expect_equal(as.numeric(coef(g1)[[1]]), c(1.272, 0.543, 0.000, 1.000), tolerance = 1e-2)
    fit3 <- M2(mod_metric)
    expect_is(fit1, 'data.frame')
    expect_true(mirt:::closeEnough(fit1$M2 - c(213.1749), -1e-2, 1e-2))
    expect_equal(fit1$SRMSR.D1, 0.04152426, tolerance = 1e-4)
    expect_equal(fit1$TLI, .99758, tolerance = 1e-4)
    expect_true(mirt:::closeEnough(fit1$df - 195, -1e-4, 1e-4))

    # missing by design
    dat[group == 'D1',1:2] <- NA
    dat[group == 'D2',14:15] <- NA
    mod <- multipleGroup(dat, 1, group, invariance = c('slopes', 'interecepts', 'free_means',
                                                       'free_var'), verbose=FALSE)
    cfs <- coef(mod, simplify=TRUE)
    expect_equal(as.vector(cfs$D1$items[1:3,1:2]), c(1.2071515, 1.1911074, 1.0377839, 1.7398616, 1.0476004, -0.1931686),
                 tolerance=1e-4)
    expect_equal(as.vector(fscores(mod)[1:3,]), c(0.7479307, 0.9566082, 0.4892834), tolerance=1e-4)
    expect_is(plot(mod, type = 'trace'), 'trellis')
    ifit <- itemfit(mod, 'X2')
    expect_equal(as.vector(ifit$D1$p.X2[1:4]), c(NaN,NaN,0.0001705316, 0.0079944214), tolerance=1e-4)

    #missing data
    set.seed(1234)
    Theta1 <- rnorm(1000, -1)
    Theta2 <- rnorm(1000, 1)
    Theta <- matrix(rbind(Theta1, Theta2))
    d <- rnorm(10,4)
    d <- cbind(d, d-1, d-2, d-3, d-4, d-5, d-6)
    a <- matrix(rlnorm(10, meanlog=.1))
    group <- factor(c(rep('g1',1000), rep('g2',1000)))

    dat <- simdata(a,d,2000, itemtype = rep('graded', 10), Theta=Theta)
    x <- multipleGroup(dat, 1, group=group, method='EM', verbose = FALSE)
    expect_is(x, 'MultipleGroupClass')
    out <- empirical_ES(x)
    expect_equal(as.numeric(out[1,]), c(-0.01984901,0.1397911,-0.02010119,0.1390092,-0.03393749,-2.687888,-0.4657692,3.53414,3.553989), tolerance=1e-4)
    out2 <- empirical_ES(x, DIF = FALSE)
    expect_equal(as.numeric(out2$Value), c(-0.8379638,1.252605,0.9030441,-0.05802242,-0.8532262,1.260376,0.9277066,-1.727569,-1.356197), tolerance=1e-4)

    dat[1,1] <- dat[2,2] <- NA
    x2 <- multipleGroup(dat, 1, group=group, method='EM', verbose = FALSE)
    expect_is(x2, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(x2)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_true(mirt:::closeEnough(cfs - c(0.54559,2.87963,2.09834,1.04898,0.08764,-1.01597,-1.88808,-2.80774,0.5721,3.66458,2.86098,2.07309,0.97108,-0.04591,-1.13943,-2.04077,2.19508,3.985,2.93076,1.92402,0.93888,-3e-04,-0.92533,-2.15966,2.75051,5.52755,4.54548,3.22867,2.2928,1.26823,0.37394,-0.49598,0.40521,2.24053,1.23686,0.3006,-0.65966,-1.68032,-2.71625,-3.65993,5.18399,3.21456,2.30067,1.19408,0.1821,-0.8039,-1.83942,-2.96744,2.50677,2.45572,1.40808,0.48888,-0.67763,-1.67403,-2.61442,-3.75363,1.90942,4.6516,3.56213,2.57738,1.45019,0.65221,-0.43677,-1.41043,2.0672,3.3566,2.50567,1.64636,0.6884,-0.28202,-1.35771,-2.31969,2.54578,2.31928,1.26861,0.29768,-0.87537,-1.86427,-2.76917,-3.71823), -1e-2, 1e-2))

    # three factor
    set.seed(12345)
    a <- matrix(c(abs(rnorm(5,1,.3)), rep(0,15),abs(rnorm(5,1,.3)),
                  rep(0,15),abs(rnorm(5,1,.3))), 15, 3)
    d <- matrix(rnorm(15,0,.7),ncol=1)
    mu <- c(-.4, -.7, .1)
    sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)
    itemtype <- rep('dich', nrow(a))
    N <- 1000
    dataset1 <- simdata(a, d, N, itemtype)
    dataset2 <- simdata(a, d, N, itemtype, mu = mu, sigma = sigma)
    dat <- rbind(dataset1, dataset2)
    group <- c(rep('D1', N), rep('D2', N))
    MGmodelg1 <- '
    F1 = 1-5
    F2 = 6-10
    F3 = 11-15'

    MGmodelg2 <- '
    F1 = 1-5
    F2 = 6-10
    F3 = 11-15
    COV = F1*F2, F1*F3, F2*F3'

    #group models
    model1 <- mirt.model(MGmodelg1, quiet = TRUE)
    model2 <- mirt.model(MGmodelg1, quiet = TRUE)
    models <- model1

    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), method = 'MHRM',
                                                 verbose = FALSE, draws = 10)
    expect_is(mod_metric, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_metric)[[1]]))[1:20]
    expect_equal(cfs, c(1.0639,0,0,0.5712,0,1,1.25,0,0,-0.5052,0,1,1.0689,0,0,-0.2563,0,1,0.9372,0),
                 tolerance = 1e-2)
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM', SE=TRUE,
                                    optimizer = 'NR')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural)[[1]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.1751,0.8707,1.4796,NA,NA,NA,NA,0.5905,0.4186,0.7624,NA,NA,NA,NA,1.4322,1.0538,1.8107,NA,NA,NA,NA,-0.5385,-0.7251,-0.3519,NA,NA,NA,NA,1.0549,0.7774,1.3325,NA,NA,NA,NA,-0.2574,-0.4131,-0.1017,NA,NA,NA,NA,0.9511,0.6867,1.2155,NA,NA,NA,NA,0.8776,0.708,1.0473,NA,NA,NA,NA,1.1473,0.8512,1.4433,NA,NA,NA,NA,0.2582,0.0973,0.4191,NA,NA,NA,NA,NA,NA,0.609,0.3694,0.8487,NA,NA,0.5172,0.3746,0.6597,NA,NA,NA,NA,NA,NA,1.59,0.9502,2.2299,NA,NA,1.1289,0.8318,1.426,NA,NA,NA,NA,NA,NA,1.0583,0.7093,1.4073,NA,NA,-0.5339,-0.6997,-0.3682,NA,NA,NA,NA,NA,NA,0.8615,0.5566,1.1664,NA,NA,-1.2163,-1.4071,-1.0255,NA,NA,NA,NA,NA,NA,0.5273,0.2866,0.768,NA,NA,-1.0236,-1.1802,-0.8669,NA,NA,NA,NA,NA,NA,NA,NA,1.0005,0.7058,1.2951,1.3027,1.1011,1.5042,NA,NA,NA,NA,NA,NA,NA,NA,1.4598,1.0451,1.8746,-0.5184,-0.7073,-0.3295,NA,NA,NA,NA,NA,NA,NA,NA,1.0278,0.7475,1.3082,0.474,0.3135,0.6345,NA,NA,NA,NA,NA,NA,NA,NA,1.4256,1.0218,1.8295,0.4234,0.2404,0.6063,NA,NA,NA,NA,NA,NA,NA,NA,0.6054,0.3957,0.8151,-0.152,-0.2888,-0.0152,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                 tolerance = 1e-2)

    fs1 <- fscores(mod_metric, verbose = FALSE, full.scores=FALSE)
    expect_is(fs1, 'list')
    expect_true(mirt:::closeEnough(fs1[[1L]][1:6, 'F3'] - c(-0.2929, -0.4214, -0.1471,  0.1641,  0.4519, -0.4214), -1e-3, 1e-3))
})
