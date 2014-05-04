context('mixedmirt')

test_that('mixed dich', {
    set.seed(1234)
    N <- 750
    a <- matrix(rlnorm(10,.2,.5),10,1)
    d <- matrix(rnorm(10), 10)
    Theta <- matrix(sort(rnorm(N)))
    pseudoIQ <- scale(Theta * 5 + 100  + rnorm(N, 0 , 5))
    group <- factor(rep(c('G1','G2','G3'), each = N/3))
    data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
    covdata <- data.frame(group, pseudoIQ)
    mixedmirt1 <- 'Theta = 1-10'
    model <- mirt.model(mixedmirt1, quiet = TRUE)

    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                       verbose = FALSE, draws = 10)
    expect_is(mod1, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -1.6963, -1.8782, -1.5145, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -2.1084, -2.296, -1.9208, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -1.7032, -1.8851, -1.5212, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -1.0479, -1.2259, -0.8698, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -0.3139, -0.4959, -0.1318, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -1.2438, -1.4224, -1.0653, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -1.6827, -1.8644, -1.501, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -2.1234, -2.3112, -1.9355, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, -2.0273, -2.2135, -1.8411, 0, NA, NA, 1, NA, NA, 1.1103, 0.9748, 1.2458, 2.2497, 2.1624, 2.337, 1, NA, NA, 1.6126, 1.3412, 1.884, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.1028, -0.0175, 0.223),
                 tolerance = 1e-2)
    names <- wald(mod1)
    L <- matrix(c(1, numeric(ncol(names) - 1L)), 1L)
    wld <- wald(mod1, L, C=as.numeric(L))
    expect_equal(wld$W[1], 2.547203, tolerance = 1e-4)

    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 10)
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, 0.0224, -0.2938, 0.3387, -1.7701, -1.9581, -1.5821, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, 0.2019, -0.0601, 0.464, -2.1949, -2.3966, -1.9933, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, 0.1026, -0.1209, 0.326, -1.7792, -1.9687, -1.5898, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, -2.4256, -3.4842, -1.367, -1.0825, -1.3595, -0.8055, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, 0.059, -0.202, 0.32, -0.3736, -0.5539, -0.1933, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, 0.0838, -0.1834, 0.3511, -1.3142, -1.4965, -1.132, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, 0.0108, -0.2815, 0.3031, -1.756, -1.9441, -1.5679, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, 0.312, 0.0342, 0.5897, -2.2213, -2.43, -2.0125, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, -0.0909, -0.3846, 0.2027, -2.1001, -2.2935, -1.9067, 0, NA, NA, 1, NA, NA, 1.1922, 1.055, 1.3294, 2.3777, 2.2458, 2.5097, -0.0866, -0.4855, 0.3124, 1.5568, 1.2878, 1.8259, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 10, verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1, NA, NA, -0.5542, -0.7049, -0.4035, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9713, -1.1252, -0.8174, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.5611, -0.7119, -0.4103, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.0964, -0.0489, 0.2418, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.8346, 0.6967, 0.9725, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.0999, -0.247, 0.0473, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.5405, -0.6911, -0.3898, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9865, -1.1405, -0.8325, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.8889, -1.0421, -0.7356, 0, NA, NA, 1, NA, NA, 1, NA, NA, 2.8345, 2.6455, 3.0235, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.3561, 0.2906, 0.4216, 0.6512, 0.3052, 0.9972),
                 tolerance = 1e-2)
})

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    mod <- mixedmirt(Science, covdat, model=model,
                                       fixed = ~ 0 + group, verbose = FALSE, draws = 10)
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(na.omit(do.call(c, coef(mod, digits=4))))
    expect_equal(cfs, c(-0.0455, -0.3298, 0.2387, 1, 0, 1, 2, 3, 0, 3.0579, 2.0316, 4.0841, 5.6485, 4.599, 6.6979, 4.2837, 3.1842, 5.3832, -0.0455, -0.3298, 0.2387, 1, 0, 1, 2, 3, 0, 1.8831, 1.4205, 2.3456, 2.7993, 2.263, 3.3355, 0.9734, 0.2594, 1.6874, -0.0455, -0.3298, 0.2387, 1, 0, 1, 2, 3, 0, 2.628, 1.9831, 3.273, 4.0487, 3.347, 4.7504, 2.9387, 2.1441, 3.7333, -0.0455, -0.3298, 0.2387, 1, 0, 1, 2, 3, 0, 2.4334, 1.8948, 2.9719, 3.3361, 2.7277, 3.9444, 2.0071, 1.2692, 2.7449, 0, 0.9545, 0.5554, 1.3537),
                 tolerance = 1e-2)

    mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE)
    expect_is(mod2, 'MixedClass')
    expect_equal(mod@df - mod2@df, 3)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod2, digits=4))))
    expect_equal(cfs, c(-0.1617, -0.3799, 0.0564, 0.8287, 0.4712, 1.1862, 0, 1, 2, 3, 0, 2.8314, 1.7069, 3.9558, 5.3694, 4.1408, 6.5979, 4.1331, 2.942, 5.3241, -0.1617, -0.3799, 0.0564, 0.8375, 0.5645, 1.1105, 0, 1, 2, 3, 0, 1.7767, 1.2475, 2.306, 2.7154, 2.0912, 3.3395, 1.0623, 0.3681, 1.7564, -0.1617, -0.3799, 0.0564, 2.5641, 0.9768, 4.1515, 0, 1, 2, 3, 0, 5.2528, 2.2768, 8.2289, 7.6976, 3.4654, 11.9299, 5.6962, 2.1506, 9.2417, -0.1617, -0.3799, 0.0564, 0.6946, 0.4353, 0.954, 0, 1, 2, 3, 0, 2.1363, 1.5668, 2.7059, 3.0002, 2.3542, 3.6462, 1.9153, 1.2202, 2.6103, 0, 1),
                 tolerance = 1e-2)

    mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE)
    expect_is(mod3, 'MixedClass')
    expect_equal(mod3@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod3, digits=4))))
    expect_equal(cfs, c(-0.198, -0.5356, 0.1396, 0.9961, 0.619, 1.3732, 4.9121, 3.9516, 5.8726, 2.7023, 2.2438, 3.1608, -1.3542, -1.7789, -0.9294, -0.198, -0.5356, 0.1396, 1.2157, 0.8826, 1.5487, 3.0075, 2.4655, 3.5495, 0.9901, 0.6019, 1.3784, -2.1652, -2.6131, -1.7174, -0.198, -0.5356, 0.1396, 2.52, 1.431, 3.609, 5.6312, 3.6983, 7.5641, 2.4423, 1.3505, 3.5342, -1.9982, -2.6425, -1.3539, -0.198, -0.5356, 0.1396, 1.0548, 0.718, 1.3917, 3.4081, 2.8535, 3.9626, 1.0723, 0.7202, 1.4243, -1.575, -2.0069, -1.1431, 0, 1),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(rmod1, digits=4))))
    expect_equal(cfs, c(0.994, 0.6734, 1.3147, 4.8198, 3.8831, 5.7565, 2.6108, 2.1982, 3.0235, -1.4401, -1.7495, -1.1308, 1.2413, 0.9069, 1.5757, 2.9377, 2.4641, 3.4113, 0.9084, 0.6162, 1.2006, -2.2677, -2.6606, -1.8748, 2.5726, 1.8406, 3.3045, 5.6288, 4.3821, 6.8755, 2.3998, 1.7592, 3.0403, -2.1022, -2.6991, -1.5053, 1.051, 0.714, 1.3879, 3.315, 2.7856, 3.8444, 0.9837, 0.7065, 1.2609, -1.6589, -1.9888, -1.3291, 0, 1, 1e-04, 0, 3e-04),
                 tolerance = 1e-2)

    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)

})
