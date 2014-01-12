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
    expect_equal(cfs, c(1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -1.6998, -1.8444, -1.5553, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -2.1137, -2.2648, -1.9626, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -1.7067, -1.8513, -1.5621, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -1.0486, -1.1933, -0.904, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -0.3117, -0.4686, -0.1548, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -1.2454, -1.3889, -1.102, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -1.6861, -1.8306, -1.5417, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -2.1287, -2.2801, -1.9772, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, -2.0322, -2.1816, -1.8828, 0, NA, NA, 1, NA, NA, 1.109, 0.9782, 1.2398, 2.2479, 2.1118, 2.3841, 1, NA, NA, 1.6213, 1.3591, 1.8835, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.124, 0.0656, 0.1823),
                 tolerance = 1e-2)

    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 10))
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, 0.1629, -0.1832, 0.5091, -1.7175, -1.8604, -1.5746, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, 0.4342, 0.2651, 0.6032, -2.1611, -2.3055, -2.0166, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, 0.2439, -0.0133, 0.5012, -1.7287, -1.8711, -1.5862, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, -1.4531, -2.0296, -0.8767, -1.056, -1.255, -0.857, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, 0.1389, -0.2221, 0.4999, -0.3306, -0.4895, -0.1717, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, 0.1208, -0.1438, 0.3855, -1.2631, -1.4038, -1.1224, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, 0.1702, -0.2803, 0.6208, -1.7041, -1.8465, -1.5617, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, 0.396, 0.0357, 0.7563, -2.1711, -2.3297, -2.0126, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, 0.1965, -0.1683, 0.5613, -2.0513, -2.2012, -1.9014, 0, NA, NA, 1, NA, NA, 1.1456, 1.1305, 1.1608, 2.2923, 2.1881, 2.3966, -0.0135, -0.541, 0.5141, 1.5809, 1.3111, 1.8507, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- suppressMessages(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 10, verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1, NA, NA, -0.5628, -0.7142, -0.4114, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9804, -1.1376, -0.8233, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.5697, -0.7212, -0.4182, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.089, -0.0559, 0.2339, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.8288, 0.6887, 0.9689, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.1077, -0.2544, 0.039, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.549, -0.7003, -0.3978, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9957, -1.153, -0.8383, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.8979, -1.0537, -0.7421, 0, NA, NA, 1, NA, NA, 1, NA, NA, 2.8328, 2.6225, 3.0431, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.3687, 0.253, 0.4844, 0.6415, 0.4113, 0.8718),
                 tolerance = 1e-2)
})

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    suppressWarnings(mod <- mixedmirt(Science, covdat, model=model,
                                       fixed = ~ 0 + group, verbose = FALSE, draws = 10))
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod, digits=4)))
    expect_equal(cfs, c(-0.0554, -0.3786, 0.2677, 1, NA, NA, 0, NA, NA, 3.0659, 2.0604, 4.0715, 5.6639, 4.7144, 6.6134, 4.3081, 3.483, 5.1332, -0.0554, -0.3786, 0.2677, 1, NA, NA, 0, NA, NA, 1.8901, 1.489, 2.2913, 2.8141, 2.6023, 3.0258, 0.9977, 0.6632, 1.3322, -0.0554, -0.3786, 0.2677, 1, NA, NA, 0, NA, NA, 2.6354, 2.0281, 3.2427, 4.0635, 3.5425, 4.5845, 2.9626, 2.7111, 3.2141, -0.0554, -0.3786, 0.2677, 1, NA, NA, 0, NA, NA, 2.4406, 1.9513, 2.9298, 3.3509, 2.9824, 3.7194, 2.0312, 1.7978, 2.2646, 0, NA, NA, 0.9521, 0.5819, 1.3223),
                 tolerance = 1e-2)

    suppressWarnings(mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE))
    expect_is(mod2, 'MixedClass')
    expect_equal(mod@df - mod2@df, 3)
    cfs <- as.numeric(do.call(c, coef(mod2, digits=4)))
    expect_equal(cfs, c(-0.1601, -0.4883, 0.1682, 0.8223, 0.5184, 1.1262, 0, NA, NA, 2.8303, 1.7347, 3.9259, 5.3709, 4.2186, 6.5232, 4.1417, 3.1204, 5.163, -0.1601, -0.4883, 0.1682, 0.8402, 0.5711, 1.1094, 0, NA, NA, 1.7885, 1.3101, 2.2669, 2.7338, 2.3147, 3.1529, 1.0821, 0.8821, 1.2821, -0.1601, -0.4883, 0.1682, 2.5353, 1.409, 3.6616, 0, NA, NA, 5.2316, 3.3755, 7.0876, 7.6801, 5.2949, 10.0652, 5.7051, 4.2908, 7.1194, -0.1601, -0.4883, 0.1682, 0.6943, 0.4433, 0.9453, 0, NA, NA, 2.1432, 1.5978, 2.6887, 3.0122, 2.4892, 3.5353, 1.9308, 1.6538, 2.2079, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    suppressWarnings(mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE))
    expect_is(mod3, 'MixedClass')
    expect_equal(mod3@df, 238)
    cfs <- as.numeric(do.call(c, coef(mod3, digits=4)))
    expect_equal(cfs, c(-0.1998, -0.4289, 0.0294, 1.0029, 0.6459, 1.3599, 4.9119, 3.9323, 5.8914, 2.7, 2.238, 3.162, -1.3591, -1.56, -1.1581, -0.1998, -0.4289, 0.0294, 1.2226, 0.8878, 1.5573, 3.0038, 2.5984, 3.4092, 0.9871, 0.7903, 1.184, -2.1713, -2.5108, -1.8318, -0.1998, -0.4289, 0.0294, 2.512, 1.1802, 3.8438, 5.5921, 3.9784, 7.2059, 2.4257, 1.8266, 3.0248, -1.9932, -2.9059, -1.0804, -0.1998, -0.4289, 0.0294, 1.0722, 0.747, 1.3975, 3.4151, 2.8775, 3.9527, 1.0718, 0.8395, 1.3041, -1.5858, -1.8141, -1.3576, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- suppressMessages(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 238)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(0.9922, 0.6575, 1.327, 4.8202, 3.8825, 5.7578, 2.6102, 2.1968, 3.0235, -1.4416, -1.7526, -1.1306, 1.2373, 0.9048, 1.5699, 2.9358, 2.4465, 3.4251, 0.9069, 0.608, 1.2059, -2.2689, -2.6271, -1.9106, 2.5693, 1.7479, 3.3908, 5.6333, 4.2822, 6.9844, 2.3994, 1.7172, 3.0815, -2.1079, -2.7195, -1.4962, 1.0455, 0.7283, 1.3626, 3.3123, 2.7875, 3.8371, 0.9819, 0.7091, 1.2546, -1.6596, -1.973, -1.3463, 0, NA, NA, 1, NA, NA, 1e-04, -3e-04, 5e-04),
                 tolerance = 1e-2)

    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)

})
