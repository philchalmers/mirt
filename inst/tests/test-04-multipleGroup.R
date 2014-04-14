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

    mod_configural <- multipleGroup(dat, models, SE=TRUE, SE.type = 'crossprod', 
                                    group = group, verbose = FALSE, method = 'EM')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural, digits=4)[[1L]]))
    cfs <- as.numeric(na.omit(cfs[cfs != 0 & cfs != 1]))
    expect_equal(cfs, c(1.0693, 0.8484, 1.2901, 0.5541, 0.392, 0.7162, 1.278, 1.0225, 1.5335, -0.6918, -0.8705, -0.513, 0.8833, 0.6898, 1.0768, -0.1375, -0.2844, 0.0094, 1.1112, 0.8848, 1.3377, 0.8295, 0.6576, 1.0014, 1.2481, 1.0102, 1.4861, 0.3265, 0.1607, 0.4922, 0.476, 0.3118, 0.6402, 0.4796, 0.3432, 0.6161, 1.1617, 0.9322, 1.3912, 1.0847, 0.9013, 1.2682, 0.8586, 0.6638, 1.0535, -0.3852, -0.5338, -0.2367, 0.89, 0.6815, 1.0985, -1.048, -1.2174, -0.8787, 0.8085, 0.6115, 1.0056, -1.0908, -1.2579, -0.9237, 0.9013, 0.6958, 1.1069, 1.1642, 0.9902, 1.3382, 1.5832, 1.2885, 1.8779, -0.135, -0.3173, 0.0473, 1.4098, 1.1458, 1.6738, 0.6542, 0.4719, 0.8365, 1.0401, 0.826, 1.2542, 0.4073, 0.2506, 0.564, 0.8804, 0.6853, 1.0754, -0.0812, -0.2278, 0.0653),
                 tolerance = 1e-2)
    expect_equal(mod_configural@df, 32707)
#     dtf <- DTF(mod_configural, digits=7)
#     cfs <- as.numeric(c(dtf$signed, dtf$unsigned))
#     expect_equal(cfs, c(-0.0562109, -0.3747396, 0.1187070, 2.2969256, 4.4389564), tolerance=1e-3)
    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), verbose = FALSE,
                                method = 'EM')
    expect_is(mod_metric, 'MultipleGroupClass')
    expect_equal(mod_metric@df, 32722)
    mod_scalar2 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov','free_means'))
    cfs <- as.numeric(do.call(c, coef(mod_scalar2, digits=4)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.1424, 0.5623, 1.3257, -0.6508, 0.9936, -0.2008, 1.0489, 0.8867, 1.1449, 0.3383, 0.4314, 0.4965, 1.2256, 1.158, 0.916, -0.4197, 0.8163, -1.0164, 0.8011, -1.0888, 0.9486, 1.2348, 1.5887, -0.1893, 1.1991, 0.5387, 1.1291, 0.4329, 0.8934, -0.117),
                 tolerance = 1e-2)
    expect_is(mod_scalar2, 'MultipleGroupClass')
    expect_equal(mod_scalar2@df, 32735)
    newmodel <- mirt.model('F = 1-15
                            CONSTRAINB = (1-15, a1), (1,2,3-15,d)')
    mod_scalar1 <- multipleGroup(dat, newmodel, group = group, verbose = FALSE, invariance='free_var')
    expect_is(mod_scalar1, 'MultipleGroupClass')
    mod_EH <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM', empiricalhist=TRUE)
    expect_is(mod_EH, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_EH, digits=4)[[1L]]))
    expect_equal(cfs, c(0.9157, 0.5762, 0, 1, 1.1492, -0.6687, 0, 1, 0.7621, -0.1143, 0, 1, 0.962, 0.8544, 0, 1, 1.084, 0.3526, 0, 1, 0.409, 0.4941, 0, 1, 1.0155, 1.1149, 0, 1, 0.757, -0.3631, 0, 1, 0.7929, -1.027, 0, 1, 0.7202, -1.0713, 0, 1, 0.7833, 1.1886, 0, 1, 1.3881, -0.1058, 0, 1, 1.2098, 0.6779, 0, 1, 0.9116, 0.4329, 0, 1, 0.7701, -0.0585, 0, 1, 0, 1),
                 tolerance = 1e-2)

    dat[1,1] <- dat[2,2] <- NA
    mod_missing <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov'))
    expect_is(mod_missing, 'MultipleGroupClass')
    expect_equal(mod_missing@df, 32736)

    fs1 <- fscores(mod_metric, verbose = FALSE)
    expect_true(mirt:::closeEnough(fs1[[1]][1:6, 'F1'] - c(-2.084760, -1.683841, -1.412181,
                                                           -1.656879, -1.324689, -1.092169), -1e-2, 1e-2))
    fs2 <- fscores(mod_metric, full.scores = TRUE)
    fs3 <- fscores(mod_missing, verbose = FALSE)
    fs4 <- fscores(mod_missing, full.scores = TRUE)
    fs5 <- fscores(mod_metric, full.scores = TRUE, scores.only=TRUE)
    expect_is(fs1, 'list')
    expect_is(fs2, 'matrix')
    expect_is(fs3, 'list')
    expect_is(fs4, 'matrix')

    fit1 <- M2(mod_metric)
    expect_is(fit1, 'data.frame')
    expect_true(mirt:::closeEnough(fit1[1:2] - c(85.28706, 67.16565), -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(fit1$df.M2 - 195, -1e-4, 1e-4))
    fit2 <- itemfit(mod_metric)
    expect_is(fit2, 'list')

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

    dat[1,1] <- dat[2,2] <- NA
    x2 <- multipleGroup(dat, 1, group=group, method='EM', verbose = FALSE)
    expect_is(x2, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(x2, digits = 5)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_true(mirt:::closeEnough(cfs - c(0.67629, 2.87845, 1.89717, 0.98919, 0.11943, -0.95038, -2.02389, -3.14849, 0.61578, 3.97519, 2.86863, 1.70278, 0.77932, -0.11022, -1.06186, -2.12238, 2.11957, 4.01324, 2.92793, 1.87974, 0.91636, -0.14342, -1.10866, -2.06283, 2.75402, 5.36174, 4.40853, 3.28294, 2.38318, 1.2709, 0.29023, -0.82956, 0.46569, 2.4063, 1.42569, 0.45862, -0.58853, -1.63218, -2.68877, -3.68465, 4.89255, 3.07378, 2.03808, 1.12532, 0.08488, -0.88429, -1.835, -3.06312, 2.49489, 2.3397, 1.28589, 0.42082, -0.50011, -1.49587, -2.52183, -3.63643, 1.96523, 4.43821, 3.47756, 2.44846, 1.45835, 0.52316, -0.37892, -1.43195, 2.03157, 3.44315, 2.69415, 1.59947, 0.66758, -0.41574, -1.42144, -2.33709, 2.35308, 2.27145, 1.25765, 0.38389, -0.64434, -1.80528, -2.80008, -3.87182), -1e-2, 1e-2))

})

test_that('three factor', {
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

    suppressWarnings(mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), method = 'MHRM',
                                                 verbose = FALSE, draws = 10))
    expect_is(mod_metric, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_metric, digits=4)[[1]]))[1:20]
    expect_equal(cfs, c(1.3033, 1.0679, 1.5387, 0, NA, NA, 0, NA, NA, 0.6597, 0.4832, 0.8361, 0, NA, NA, 1, NA, NA, 1.2852, 1.0436),
                 tolerance = 1e-2)
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural, digits=4)[[1]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.2934, 0.655, 1.2387, -0.57, 0.9276, -0.1996, 0.8177, 0.7967, 1.0713, 0.2166, 0.4807, 0.61, 1.1778, 0.9948, 0.9453, -0.4464, 1.0761, -1.18, 0.8664, -1.1451, 0.8854, 1.3121, 1.4964, -0.3002, 1.0534, 0.4406, 1.0614, 0.4569, 0.8831, -0.1871),
                 tolerance = 1e-2)

    fs1 <- fscores(mod_metric, verbose = FALSE)
    expect_is(fs1, 'list')
    expect_true(mirt:::closeEnough(fs1[[1L]][1:6, 'F3'] - c(-1.4347448, -0.9987337, -0.8860504, -0.5417926, -0.8361526, -0.4013042), -1e-4, 1e-4))
})
