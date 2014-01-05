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
    mod1 <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                       verbose = FALSE, draws = 10))
    expect_is(mod1, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -1.7032, -1.8947, -1.5117, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -2.1171, -2.3188, -1.9155, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -1.7101, -1.9018, -1.5184, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -1.0519, -1.2312, -0.8725, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -0.3147, -0.4864, -0.1429, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -1.2487, -1.4313, -1.0662, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -1.6895, -1.8807, -1.4983, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -2.1322, -2.3342, -1.9301, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, -2.0357, -2.2352, -1.8362, 0, NA, NA, 1, NA, NA, 1.1095, -0.2692, 2.4883, 2.2537, 0.8488, 3.6585, 1, NA, NA, 1.619, 1.3794, 1.8587, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.1254, -0.0429, 0.2936),
                 tolerance = 1e-2)

    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 10))
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.1835, -0.2476, 0.6146, -1.7122, -1.9023, -1.5221, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.4809, 0.0503, 0.9115, -2.1619, -2.3827, -1.9412, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.2635, -0.0436, 0.5706, -1.7236, -1.9155, -1.5318, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, -1.3626, -1.7702, -0.9551, -1.0484, -1.2759, -0.821, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.1472, -0.3208, 0.6152, -0.3258, -0.5081, -0.1436, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.1272, -0.2683, 0.5226, -1.2574, -1.439, -1.0759, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.2008, -0.2872, 0.6888, -1.6994, -1.8911, -1.5076, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.4034, -0.1106, 0.9175, -2.1656, -2.3837, -1.9476, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.2395, -0.1088, 0.5878, -2.0484, -2.2503, -1.8465, 0, NA, NA, 1, NA, NA, 1.1395, -0.2491, 2.5282, 2.2828, 0.8994, 3.6662, 0.0079, -0.4864, 0.5022, 1.5855, 1.313, 1.858, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- suppressMessages(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 10, verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1, NA, NA, -0.5481, -0.6984, -0.3978, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9631, -1.1198, -0.8063, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.555, -0.7054, -0.4046, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.0991, -0.0465, 0.2446, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.833, 0.6863, 0.9797, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.0962, -0.2426, 0.0503, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.5344, -0.6846, -0.3843, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9782, -1.1352, -0.8211, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.881, -1.0363, -0.7258, 0, NA, NA, 1, NA, NA, 1, NA, NA, 2.8197, 2.5906, 3.0489, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.3823, 0.2857, 0.4789, 0.5595, 0.3401, 0.7788),
                 tolerance = 1e-2)
})

test_that('item and group predictors', {
    data <- key2binary(SAT12,
                       key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
    model <- mirt.model('Theta = 1-32', quiet = TRUE)

    itemdesign <- data.frame(itemorder = factor(c(rep('easier', 16), rep('harder', 16))))
    fs <- scale(rowSums(data))
    covdata <- data.frame(gender=ifelse(fs > 1, 'M', 'F'))

    sv <- mixedmirt(data, covdata, model = model, fixed = ~ 0 + itemorder + gender,
                    itemdesign = itemdesign, pars = 'values')
    expect_is(sv, 'data.frame')
    suppressWarnings(LLTM <- mixedmirt(data, covdata, model = model, fixed = ~ 0 + itemorder + gender,
                      itemdesign = itemdesign, verbose = FALSE, draws = 10))
    expect_is(LLTM, 'MixedClass')
    cfs <- na.omit(as.numeric(do.call(c, coef(LLTM, digits=4))))[1:20]
    sv2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ 0 + itemorder * gender,
                     itemdesign = itemdesign, pars='values'))
    expect_is(sv2, 'data.frame')
    LLTM2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ 0 + itemorder * gender,
                       itemdesign = itemdesign, verbose = FALSE, draws = 10))
    expect_is(LLTM2, 'MixedClass')
    expect_equal(LLTM@df - LLTM2@df, 1)
    expect_equal(LLTM@logLik - LLTM2@logLik, 1.622451, tolerance = 1e-2)
    cfs <- na.omit(as.numeric(do.call(c, coef(LLTM2, digits=4))))[1:20]
    expect_equal(cfs, c(-0.0638, -1.7593, 1.6318, 0.2909, -1.7351, 2.317, 1.5022, -3.3885, 6.3929, 0.0027, -6.4064, 6.4118, 1, 0, 0, 1, -0.0638, -1.7593, 1.6318, 0.2909),
                 tolerance = 1e-2)
})

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    suppressWarnings(mod <- mixedmirt(Science, covdat, model=model,
                                       fixed = ~ 0 + group, verbose = FALSE, draws = 10))
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod, digits=4)))
    expect_equal(cfs, c(-0.0294, -1.0104, 0.9516, 1, NA, NA, 0, NA, NA, 3.0644, 2.0561, 4.0726, 5.6564, 4.6293, 6.6836, 4.2933, 3.1925, 5.3941, -0.0294, -1.0104, 0.9516, 1, NA, NA, 0, NA, NA, 1.8851, 1.4309, 2.3393, 2.8022, 2.2675, 3.337, 0.978, 0.29, 1.666, -0.0294, -1.0104, 0.9516, 1, NA, NA, 0, NA, NA, 2.6317, 2.0004, 3.263, 4.0533, 3.3625, 4.7441, 2.9451, 2.1507, 3.7395, -0.0294, -1.0104, 0.9516, 1, NA, NA, 0, NA, NA, 2.4361, 1.9078, 2.9643, 3.3397, 2.7354, 3.9439, 2.0126, 1.2857, 2.7395, 0, NA, NA, 0.9567, 0.7133, 1.2),
                 tolerance = 1e-2)

    suppressWarnings(mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE))
    expect_is(mod2, 'MixedClass')
    expect_equal(mod@df - mod2@df, 3)
    cfs <- as.numeric(do.call(c, coef(mod2, digits=4)))
    expect_equal(cfs, c(-0.1596, -1.0795, 0.7602, 0.8094, 0.5062, 1.1125, 0, NA, NA, 2.7928, 1.7047, 3.8809, 5.3159, 4.146, 6.4858, 4.0871, 2.9364, 5.2377, -0.1596, -1.0795, 0.7602, 0.8402, 0.5489, 1.1314, 0, NA, NA, 1.7784, 1.2527, 2.304, 2.7163, 2.1091, 3.3235, 1.0601, 0.409, 1.7112, -0.1596, -1.0795, 0.7602, 2.554, 1.0817, 4.0262, 0, NA, NA, 5.237, 2.4815, 7.9926, 7.6783, 3.8662, 11.4905, 5.6862, 2.5845, 8.788, -0.1596, -1.0795, 0.7602, 0.7038, 0.4562, 0.9515, 0, NA, NA, 2.1451, 1.5707, 2.7195, 3.0102, 2.3556, 3.6648, 1.9183, 1.2302, 2.6064, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    suppressWarnings(mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE))
    expect_is(mod3, 'MixedClass')
    expect_equal(mod3@df, 238)
    cfs <- as.numeric(do.call(c, coef(mod3, digits=4)))
    expect_equal(cfs, c(-0.1998, -1.3184, 0.9188, 1.0029, 0.8966, 1.1092, 4.9119, 3.9579, 5.8658, 2.7, 2.2554, 3.1446, -1.3591, -1.476, -1.2422, -0.1998, -1.3184, 0.9188, 1.2226, 0.8919, 1.5533, 3.0038, 2.6098, 3.3979, 0.9871, 0.8047, 1.1696, -2.1713, -2.4996, -1.843, -0.1998, -1.3184, 0.9188, 2.512, 1.1171, 3.9069, 5.5921, 3.2147, 7.9696, 2.4257, 1.1041, 3.7473, -1.9932, -2.4603, -1.526, -0.1998, -1.3184, 0.9188, 1.0722, 0.8211, 1.3234, 3.4151, 2.8754, 3.9549, 1.0718, 0.8235, 1.3201, -1.5858, -1.6473, -1.5244, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- suppressMessages(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 238)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(0.9922, 0.6573, 1.3272, 4.8202, 3.8824, 5.758, 2.6102, 2.1966, 3.0237, -1.4416, -1.7527, -1.1305, 1.2373, 0.9047, 1.57, 2.9358, 2.4465, 3.4251, 0.9069, 0.6079, 1.2059, -2.2689, -2.6277, -1.91, 2.5693, 1.7478, 3.3908, 5.6333, 4.2819, 6.9846, 2.3994, 1.7171, 3.0816, -2.1079, -2.7205, -1.4953, 1.0455, 0.7283, 1.3626, 3.3123, 2.7875, 3.8371, 0.9819, 0.709, 1.2547, -1.6596, -1.9733, -1.3459, 0, NA, NA, 1, NA, NA, 1e-04, -3e-04, 5e-04),
                 tolerance = 1e-2)

    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)

})
