context('mixedmirt')

test_that('mixed dich', {
    set.seed(1234)
    N <- 750
    a <- matrix(rlnorm(10,.2,.5),10,1)
    d <- matrix(rnorm(10), 10)
    Theta <- matrix(sort(rnorm(N)))
    pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
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
    expect_equal(cfs, c(1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -1.6944, -1.8861, -1.5026, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -2.1067, -2.3075, -1.9059, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -1.7012, -1.8931, -1.5093, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -1.0454, -1.2295, -0.8613, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -0.3108, -0.4964, -0.1251, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -1.2415, -1.4272, -1.0559, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -1.6807, -1.8722, -1.4892, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -2.1216, -2.3228, -1.9204, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, -2.0255, -2.2243, -1.8268, 0, NA, NA, 1, NA, NA, 1.112, 0.9515, 1.2726, 2.2478, 2.1, 2.3956, 1, NA, NA, 1.6168, 1.3397, 1.8939, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.1057, -0.0086, 0.2201),
                 tollerance = 1e-2)
    
    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, 
                                        itemtype = '2PL', verbose = FALSE, draws = 10))
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 293) 
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, 0.0419, -0.2085, 0.2923, -1.7627, -1.955, -1.5703, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, 0.2383, 0.0166, 0.46, -2.1902, -2.3958, -1.9845, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, 0.1176, -0.0972, 0.3324, -1.772, -1.9657, -1.5783, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, -2.2564, -2.3686, -2.1442, -1.0794, -1.3754, -0.7835, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, 0.0714, -0.2436, 0.3865, -0.3678, -0.5506, -0.1849, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, 0.0955, -0.2731, 0.4642, -1.3073, -1.493, -1.1217, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, 0.0462, -0.2234, 0.3159, -1.7489, -1.9411, -1.5567, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, 0.3284, -0.0347, 0.6914, -2.2154, -2.4293, -2.0015, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, -0.0514, -0.2849, 0.1821, -2.0918, -2.2908, -1.8928, 0, NA, NA, 1, NA, NA, 1.1851, 1.0417, 1.3284, 2.3645, 2.217, 2.512, -0.0592, -0.5061, 0.3878, 1.5601, 1.2901, 1.83, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- suppressMessages(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group, 
                                        draws = 10, verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 303) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1, NA, NA, -0.5481, -0.6967, -0.3995, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9631, -1.1034, -0.8227, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.555, -0.7035, -0.4064, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.0991, -0.0383, 0.2364, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.833, 0.7715, 0.8945, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.0962, -0.2405, 0.0481, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.5344, -0.6831, -0.3857, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9782, -1.118, -0.8384, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.881, -1.0239, -0.7382, 0, NA, NA, 1, NA, NA, 1, NA, NA, 2.8197, 2.5741, 3.0654, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.3823, 0.3408, 0.4237, 0.5595, 0.3291, 0.7898),
                 tollerance = 1e-2)
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
    expect_equal(cfs, c(-0.062, -0.1115, -0.0125, 0.2291, 0.1795, 0.2787, 1.5033, 1.3819, 1.6247, 1, 0, 0, 1, -0.062, -0.1115, -0.0125, 0.2291, 0.1795, 0.2787, 1.5033),
                 tollerance = 1e-2)
    
    sv2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                     itemdesign = itemdesign, pars='values'))
    expect_is(sv2, 'data.frame')       
    LLTM2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                       itemdesign = itemdesign, verbose = FALSE, draws = 10))
    expect_is(LLTM2, 'MixedClass')
    expect_equal(LLTM@df - LLTM2@df, 1)
    cfs <- na.omit(as.numeric(do.call(c, coef(LLTM2, digits=4))))[1:20]
    expect_equal(cfs, c(-0.0611, -0.1075, -0.0147, 0.2902, 0.2277, 0.3528, 1.5006, 1.3504, 1.6508, 0.0032, -0.1969, 0.2033, 1, 0, 0, 1, -0.0611, -0.1075, -0.0147, 0.2902),
                 tollerance = 1e-2)
}) 

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    suppressWarnings(mod <- mixedmirt(Science, covdat, model=model,
                                       fixed = ~ 0 + group, verbose = FALSE, draws = 10))
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod, digits=4)))
    expect_equal(cfs, c(-0.0436, -0.3059, 0.2187, 1, NA, NA, 0, NA, NA, 3.0499, 2.0438, 4.056, 5.6391, 4.6066, 6.6716, 4.273, 3.1576, 5.3884, -0.0436, -0.3059, 0.2187, 1, NA, NA, 0, NA, NA, 1.8812, 1.4198, 2.3425, 2.7964, 2.2369, 3.3559, 0.9689, 0.2482, 1.6895, -0.0436, -0.3059, 0.2187, 1, NA, NA, 0, NA, NA, 2.6243, 1.9903, 3.2582, 4.0441, 3.3362, 4.7519, 2.9327, 2.1097, 3.7558, -0.0436, -0.3059, 0.2187, 1, NA, NA, 0, NA, NA, 2.4307, 1.8974, 2.964, 3.3324, 2.707, 3.9578, 2.002, 1.2438, 2.7603, 0, NA, NA, 0.955, 0.7229, 1.1871),
                 tollerance = 1e-2)
    
    suppressWarnings(mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE))    
    expect_is(mod2, 'MixedClass')   
    expect_equal(mod@df - mod2@df, 3)   
    cfs <- as.numeric(do.call(c, coef(mod2, digits=4)))
    expect_equal(cfs, c(-0.1624, -0.3909, 0.066, 0.8239, 0.4218, 1.226, 0, NA, NA, 2.8248, 1.6895, 3.9602, 5.3661, 4.1292, 6.6031, 4.14, 2.9473, 5.3327, -0.1624, -0.3909, 0.066, 0.843, 0.4582, 1.2278, 0, NA, NA, 1.7891, 1.2558, 2.3225, 2.7367, 2.1323, 3.341, 1.0878, 0.3586, 1.8169, -0.1624, -0.3909, 0.066, 2.5656, -0.3243, 5.4554, 0, NA, NA, 5.269, -0.1341, 10.672, 7.7369, -0.0726, 15.5464, 5.7573, -0.6232, 12.1378, -0.1624, -0.3909, 0.066, 0.6961, 0.3259, 1.0663, 0, NA, NA, 2.1428, 1.5483, 2.7372, 3.0133, 2.3483, 3.6783, 1.9344, 1.2194, 2.6494, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE))    
    expect_is(mod3, 'MixedClass')   
    expect_equal(mod3@df, 72) 
    cfs <- as.numeric(do.call(c, coef(mod3, digits=4)))
    expect_equal(cfs, c(-0.198, -0.5388, 0.1428, 0.9961, 0.6013, 1.3909, 4.9121, 3.9513, 5.8728, 2.7023, 2.2436, 3.161, -1.3542, -1.8017, -0.9066, -0.198, -0.5388, 0.1428, 1.2157, 0.8659, 1.5655, 3.0075, 2.4336, 3.5814, 0.9901, 0.5793, 1.401, -2.1652, -2.6144, -1.7161, -0.198, -0.5388, 0.1428, 2.52, 1.378, 3.6619, 5.6312, 3.5788, 7.6836, 2.4423, 1.2716, 3.613, -1.9982, -2.6426, -1.3539, -0.198, -0.5388, 0.1428, 1.0548, 0.6953, 1.4144, 3.4081, 2.8535, 3.9627, 1.0723, 0.7146, 1.4299, -1.575, -2.0349, -1.115, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- suppressMessages(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 72) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(0.9922, 0.6441, 1.3404, 4.8202, 3.8699, 5.7705, 2.6102, 2.1813, 3.0391, -1.4416, -1.7574, -1.1258, 1.2373, 0.8466, 1.6281, 2.9358, 2.4093, 3.4623, 0.9069, 0.5963, 1.2176, -2.2689, -2.6851, -1.8526, 2.5693, 1.3378, 3.8009, 5.6333, 3.7526, 7.514, 2.3994, 1.5434, 3.2553, -2.1079, -2.9177, -1.298, 1.0455, 0.6989, 1.392, 3.3123, 2.7619, 3.8627, 0.9819, 0.6997, 1.264, -1.6596, -2.0044, -1.3149, 0, NA, NA, 1, NA, NA, 1e-04, -3e-04, 5e-04),
                 tollerance = 1e-2)
    
    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)
    
}) 
