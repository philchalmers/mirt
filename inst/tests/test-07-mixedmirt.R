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
    expect_equal(cfs, c(1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -1.7032, -1.8969, -1.5095, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -2.1171, -2.3193, -1.915, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -1.7101, -1.9039, -1.5163, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -1.0519, -1.2393, -0.8645, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -0.3147, -0.5061, -0.1232, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -1.2487, -1.4371, -1.0603, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -1.6895, -1.883, -1.4961, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -2.1322, -2.3346, -1.9297, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, -2.0357, -2.2358, -1.8355, 0, NA, NA, 1, NA, NA, 1.1095, -0.4749, 2.694, 2.2537, 0.6226, 3.8847, 1, NA, NA, 1.619, 1.332, 1.9061, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.1254, 0.0025, 0.2482),
                 tollerance = 1e-2)
    
    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, 
                                        itemtype = '2PL', verbose = FALSE, draws = 10))
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 1001) 
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.1835, -0.219, 0.5859, -1.7122, -1.9059, -1.5185, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.4809, 0.0784, 0.8834, -2.1619, -2.3814, -1.9425, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.2635, -0.0149, 0.542, -1.7236, -1.918, -1.5293, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, -1.3626, -1.4658, -1.2595, -1.0484, -1.2751, -0.8218, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.1472, -0.3634, 0.6578, -0.3258, -0.5123, -0.1393, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.1272, -0.2891, 0.5434, -1.2574, -1.4412, -1.0737, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.2008, -0.2746, 0.6762, -1.6994, -1.892, -1.5068, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.4034, -0.1509, 0.9578, -2.1656, -2.3807, -1.9505, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.2395, 0.0038, 0.4752, -2.0484, -2.2496, -1.8472, 0, NA, NA, 1, NA, NA, 1.1395, -0.2813, 2.5603, 2.2828, 0.7932, 3.7724, 0.0079, -0.5987, 0.6145, 1.5855, 1.3083, 1.8626, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- suppressMessages(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group, 
                                        draws = 10, verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 1011) 
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
    expect_equal(cfs, c(-0.0688, -1.7315, 1.5939, 0.2223, -1.4719, 1.9165, 1.5018, -2.5686, 5.5723, 1, 0, 0, 1, -0.0688, -1.7315, 1.5939, 0.2223, -1.4719, 1.9165, 1.5018),
                 tollerance = 1e-2)
    
    sv2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                     itemdesign = itemdesign, pars='values'))
    expect_is(sv2, 'data.frame')       
    LLTM2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                       itemdesign = itemdesign, verbose = FALSE, draws = 10))
    expect_is(LLTM2, 'MixedClass')
    expect_equal(LLTM@df - LLTM2@df, 1)
    cfs <- na.omit(as.numeric(do.call(c, coef(LLTM2, digits=4))))[1:20]
    expect_equal(cfs, c(-0.0638, -1.7594, 1.6319, 0.2909, -1.7351, 2.317, 1.5022, -3.3888, 6.3932, 0.0027, -6.4064, 6.4118, 1, 0, 0, 1, -0.0638, -1.7594, 1.6319, 0.2909),
                 tollerance = 1e-2)
}) 

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    suppressWarnings(mod <- mixedmirt(Science, covdat, model=model,
                                       fixed = ~ 0 + group, verbose = FALSE, draws = 10))
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod, digits=4)))
    expect_equal(cfs, c(-0.0294, -1.0108, 0.952, 1, NA, NA, 0, NA, NA, 3.0644, 2.0561, 4.0726, 5.6564, 4.6292, 6.6837, 4.2933, 3.1917, 5.395, -0.0294, -1.0108, 0.952, 1, NA, NA, 0, NA, NA, 1.8851, 1.4306, 2.3395, 2.8022, 2.2657, 3.3388, 0.978, 0.2857, 1.6703, -0.0294, -1.0108, 0.952, 1, NA, NA, 0, NA, NA, 2.6317, 2.0004, 3.263, 4.0533, 3.3618, 4.7448, 2.9451, 2.1483, 3.742, -0.0294, -1.0108, 0.952, 1, NA, NA, 0, NA, NA, 2.4361, 1.9077, 2.9644, 3.3397, 2.7342, 3.9452, 2.0126, 1.2822, 2.743, 0, NA, NA, 0.9567, 0.7124, 1.2009),
                 tollerance = 1e-2)
    
    suppressWarnings(mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE))    
    expect_is(mod2, 'MixedClass')   
    expect_equal(mod@df - mod2@df, 3)   
    cfs <- as.numeric(do.call(c, coef(mod2, digits=4)))
    expect_equal(cfs, c(-0.1596, -1.0801, 0.7608, 0.8094, 0.4907, 1.1281, 0, NA, NA, 2.7928, 1.6957, 3.89, 5.3159, 4.1314, 6.5004, 4.0871, 2.9299, 5.2442, -0.1596, -1.0801, 0.7608, 0.8402, 0.5488, 1.1315, 0, NA, NA, 1.7784, 1.2513, 2.3055, 2.7163, 2.106, 3.3266, 1.0601, 0.4048, 1.7154, -0.1596, -1.0801, 0.7608, 2.554, 0.9379, 4.17, 0, NA, NA, 5.237, 2.2169, 8.2571, 7.6783, 3.49, 11.8667, 5.6862, 2.2896, 9.0829, -0.1596, -1.0801, 0.7608, 0.7038, 0.4343, 0.9734, 0, NA, NA, 2.1451, 1.5618, 2.7284, 3.0102, 2.3452, 3.6752, 1.9183, 1.2293, 2.6073, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE))    
    expect_is(mod3, 'MixedClass')   
    expect_equal(mod3@df, 238) 
    cfs <- as.numeric(do.call(c, coef(mod3, digits=4)))
    expect_equal(cfs, c(-0.1998, -1.3379, 0.9384, 1.0029, 0.8477, 1.1581, 4.9119, 3.9579, 5.8658, 2.7, 2.2547, 3.1453, -1.3591, -1.4661, -1.2521, -0.1998, -1.3379, 0.9384, 1.2226, 0.8916, 1.5535, 3.0038, 2.5953, 3.4124, 0.9871, 0.7722, 1.202, -2.1713, -2.5307, -1.8119, -0.1998, -1.3379, 0.9384, 2.512, 1.3263, 3.6977, 5.5921, 3.5636, 7.6206, 2.4257, 1.2906, 3.5608, -1.9932, -2.4935, -1.4928, -0.1998, -1.3379, 0.9384, 1.0722, 0.82, 1.3245, 3.4151, 2.8688, 3.9615, 1.0718, 0.801, 1.3426, -1.5858, -1.7324, -1.4393, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- suppressMessages(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 238) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(0.9922, 0.6441, 1.3404, 4.8202, 3.8699, 5.7705, 2.6102, 2.1813, 3.0391, -1.4416, -1.7574, -1.1258, 1.2373, 0.8466, 1.6281, 2.9358, 2.4093, 3.4623, 0.9069, 0.5963, 1.2176, -2.2689, -2.6851, -1.8526, 2.5693, 1.3378, 3.8009, 5.6333, 3.7526, 7.514, 2.3994, 1.5434, 3.2553, -2.1079, -2.9177, -1.298, 1.0455, 0.6989, 1.392, 3.3123, 2.7619, 3.8627, 0.9819, 0.6997, 1.264, -1.6596, -2.0044, -1.3149, 0, NA, NA, 1, NA, NA, 1e-04, -3e-04, 5e-04),
                 tollerance = 1e-2)
    
    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)
    
}) 
