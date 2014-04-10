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
    dtf <- DTF(mod_configural, digits=7)
    cfs <- as.numeric(c(dtf$signed, dtf$unsigned))
    expect_equal(cfs, c(-0.0562109, -0.3747396, 0.1187070, 2.2969256, 4.4389564), tolerance=1e-3)
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
    expect_equal(cfs, c(0.9365, 0.5812, 0, 1, 1.1684, -0.6627, 0, 1, 0.7785, -0.1108, 0, 1, 0.9798, 0.8588, 0, 1, 1.1043, 0.3577, 0, 1, 0.4196, 0.4956, 0, 1, 1.0319, 1.1185, 0, 1, 0.7728, -0.3598, 0, 1, 0.8059, -1.0227, 0, 1, 0.7308, -1.0671, 0, 1, 0.7971, 1.1918, 0, 1, 1.4157, -0.0992, 0, 1, 1.2363, 0.6845, 0, 1, 0.9261, 0.4367, 0, 1, 0.7865, -0.0549, 0, 1, 0, 1),
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
    expect_true(mirt:::closeEnough(cfs - c(0.67621, 2.87682, 1.89559, 0.98763, 0.11786, -0.95196, -2.02549, -3.15013, 0.61568, 3.97356, 2.86709, 1.70132, 0.77785, -0.11162, -1.06331, -2.12382, 2.11904, 4.00787, 2.92268, 1.87459, 0.91128, -0.14847, -1.11366, -2.06776, 2.75261, 5.35379, 4.40072, 3.2754, 2.37591, 1.26395, 0.28351, -0.83603, 0.46562, 2.40519, 1.4246, 0.45753, -0.58961, -1.63328, -2.68991, -3.68585, 4.89091, 3.06113, 2.02537, 1.1127, 0.0724, -0.89656, -1.84694, -3.07453, 2.49477, 2.33401, 1.28, 0.4148, -0.50622, -1.50201, -2.52797, -3.64251, 1.96462, 4.43303, 3.47249, 2.44355, 1.45353, 0.51842, -0.38359, -1.43651, 2.03112, 3.43814, 2.68918, 1.5946, 0.66274, -0.42055, -1.42621, -2.34182, 2.35303, 2.26612, 1.25215, 0.37826, -0.6501, -1.8111, -2.80592, -3.87765), -1e-2, 1e-2))

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
