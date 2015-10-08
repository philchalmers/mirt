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
    cfs <- as.numeric(na.omit(do.call(rbind, coef(mod_Rasch, digits=5, printSE=TRUE, as.data.frame=TRUE))))
    expect_equal(cfs, c(0.545, -0.63665, -0.14376, 0.80569, 0.30492, 0.55009, 1.03783, -0.40416, -1.08672, -1.15675, 1.20202, -0.10968, 0.58072, 0.40407, -0.08536, 0.60372, -0.4505, -0.19841, 1.02121, 0.41623, 0.66369, 1.19678, -0.4028, -1.0233, -1.11901, 1.38872, -0.11527, 0.50138, 0.50138, -0.09971, 0.07974, 0.07977, 0.07809, 0.08083, 0.07878, 0.07937, 0.08405, 0.0783, 0.0834, 0.08447, 0.08523, 0.07865, 0.07975, 0.07838, 0.07795, 0.08454, 0.08442, 0.08268, 0.0882, 0.08401, 0.08539, 0.09027, 0.08389, 0.08684, 0.0885, 0.09122, 0.08414, 0.08454, 0.08447, 0.08317),
                 tolerance = 1e-3)
    EAP <- fscores(mod_Rasch, full.scores=TRUE)
    expect_equal(cor(EAP, rowSums(dat))[1], .99, tolerance = 1e-2)
    pf <- personfit(mod_Rasch, Theta=EAP)
    pffit <- c(as.numeric(as.matrix(head(pf))), as.numeric(as.matrix(tail(pf))))
    expect_equal(pffit, c(0.6388, 0.8582, 0.83267, 0.99725, 1.05999, 1.03414, -1.51153, -0.24137, -0.19464, 0.09861, 0.3459, 0.24062, 0.699, 0.89307, 0.84912, 1.08481, 0.98964, 1.05882, -1.48336, -0.26539, -0.31212, 0.40918, -0.00498, 0.40736, 1.39643, 0.34987, 0.38458, -0.21882, -0.05278, -0.30241, 1.15408, 1.54886, 0.69589, 1.64176, 1.054, 0.47791, 0.5427, 1.01894, -1.43216, 1.13883, 0.30973, -0.96363, 0.91257, 0.98546, 0.73961, 1.10462, 0.96747, 0.66029, -0.29309, 0.09788, -1.5242, 0.38551, -0.12126, -0.82876, 0.12927, -0.21665, 1.41519, -0.5449, 0.05218, 0.89503),
                 tolerance = 1e-3)
    mod_QMCEM <- multipleGroup(dat, models, group=group, method = 'QMCEM', verbose=FALSE,
                               optimizer='NR')
    expect_equal(mod_QMCEM@logLik, -17859.03, tolerance=1e-2)
    mod_configural <- multipleGroup(dat, models, SE=TRUE, SE.type = 'crossprod', optimizer='NR',
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
                                method = 'EM', optimizer = 'NR')
    expect_is(mod_metric, 'MultipleGroupClass')
    expect_equal(mod_metric@df, 32722)
    mod_scalar2 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_var','free_means'))
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
    mod_EH <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                            empiricalhist=TRUE, optimizer = 'NR')
    expect_is(mod_EH, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_EH, digits=4)[[1L]]))
    expect_equal(cfs, c(0.9339, 0.5807, 0, 1, 1.1652, -0.6633, 0, 1, 0.7763, -0.1112, 0, 1, 0.9771, 0.8583, 0, 1, 1.1012, 0.3572, 0, 1, 0.4185, 0.4954, 0, 1, 1.0291, 1.118, 0, 1, 0.7706, -0.3602, 0, 1, 0.8036, -1.0231, 0, 1, 0.7288, -1.0674, 0, 1, 0.7949, 1.1914, 0, 1, 1.4117, -0.0999, 0, 1, 1.2329, 0.6838, 0, 1, 0.9236, 0.4362, 0, 1, 0.7843, -0.0553, 0, 1, 0, 1),
                 tolerance = 1e-2)

    dat[1,1] <- dat[2,2] <- NA
    mod_missing <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_var'))
    expect_is(mod_missing, 'MultipleGroupClass')
    expect_equal(mod_missing@df, 32736)

    fs1 <- fscores(mod_metric, verbose = FALSE, full.scores=FALSE)
    expect_true(mirt:::closeEnough(fs1[[1]][1:6, 'F1'] - c(-2.084760, -1.683841, -1.412181,
                                                           -1.324478, -1.091952, -1.741399), -1e-2, 1e-2))
    fs2 <- fscores(mod_metric, full.scores = TRUE, full.scores.SE=TRUE, method = 'ML')
    expect_equal(as.numeric(head(fs2)), c(0.5531893,  1.1960187,  1.8287234,  1.1133180, -0.5164821, -0.2322618,  0.4968711,  0.6002451,  0.7688432,
                                          0.5831489,  0.4716423,  0.4613523),
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
    expect_true(mirt:::closeEnough(fit1[1:2] - c(85.28706, 67.16565), -1e-2, 1e-2))
    expect_equal(fit1$D1.SRMSR, 0.03606703, tolerance = 1e-4)
    expect_equal(fit1$TLI, 1.005683, tolerance = 1e-4)
    expect_true(mirt:::closeEnough(fit1$df - 195, -1e-4, 1e-4))
    fit2 <- itemfit(mod_metric, digits = 20)
    expect_is(fit2, 'list')
    expect_equal(as.numeric(fit2[[1]][1L,]), c(1.000000, 2.488648, 8.294027, 11.000000, 0.686700),
                 tolerance = 1e-4)

    g1 <- extract.group(mod_metric, 1)
    expect_equal(as.numeric(coef(g1)[[1]]), c(1.252, 0.575, 0.000, 1.000), tolerance = 1e-2)

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

    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), method = 'MHRM',
                                                 verbose = FALSE, draws = 10)
    expect_is(mod_metric, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_metric, digits=4)[[1]]))[1:20]
    expect_equal(cfs, c(1.1336,0,0,0.7018,0,1,1.3382,0,0,-0.566,0,1,0.8604,0,0,-0.2489,0,1,0.7673,0),
                 tolerance = 1e-2)
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM', SE=TRUE,
                                    optimizer = 'NR')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural, digits=4)[[1]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.4228, 1.0231, 1.8226, NA, NA, NA, NA, 0.7648, 0.5647, 0.9649, NA, NA, NA, NA, 1.151, 0.8359, 1.4661, NA, NA, NA, NA, -0.5379, -0.7068, -0.369, NA, NA, NA, NA, 0.8433, 0.5984, 1.0881, NA, NA, NA, NA, -0.25, -0.3965, -0.1034, NA, NA, NA, NA, 0.6891, 0.4615, 0.9167, NA, NA, NA, NA, 0.7509, 0.599, 0.9027, NA, NA, NA, NA, 1.3615, 0.9885, 1.7346, NA, NA, NA, NA, 0.2706, 0.0969, 0.4442, NA, NA, NA, NA, NA, NA, 0.454, 0.2234, 0.6846, NA, NA, 0.4949, 0.3577, 0.6322, NA, NA, NA, NA, NA, NA, 1.1874, 0.7276, 1.6472, NA, NA, 0.9519, 0.735, 1.1689, NA, NA, NA, NA, NA, NA, 0.8536, 0.5468, 1.1604, NA, NA, -0.5404, -0.696, -0.3848, NA, NA, NA, NA, NA, NA, 0.9625, 0.6074, 1.3175, NA, NA, -1.1942, -1.3984, -0.9899, NA, NA, NA, NA, NA, NA, 0.8721, 0.5555, 1.1888, NA, NA, -1.0898, -1.2742, -0.9053, NA, NA, NA, NA, NA, NA, NA, NA, 0.7713, 0.5174, 1.0252, 1.3637, 1.1783, 1.549, NA, NA, NA, NA, NA, NA, NA, NA, 1.4962, 1.058, 1.9343, -0.4238, -0.6112, -0.2364, NA, NA, NA, NA, NA, NA, NA, NA, 1.2049, 0.8765, 1.5334, 0.4431, 0.2751, 0.6111, NA, NA, NA, NA, NA, NA, NA, NA, 1.0857, 0.7886, 1.3827, 0.3515, 0.1927, 0.5104, NA, NA, NA, NA, NA, NA, NA, NA, 0.7394, 0.5085, 0.9704, -0.0684, -0.2096, 0.0729, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                 tolerance = 1e-2)

    fs1 <- fscores(mod_metric, verbose = FALSE, full.scores=FALSE)
    expect_is(fs1, 'list')
    expect_true(mirt:::closeEnough(fs1[[1L]][1:6, 'F3'] - c(-0.9750,  0.0475, -0.5315, -0.3341, 0.5062, -0.9750), -1e-3, 1e-3))
})
