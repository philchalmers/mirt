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
    expect_equal(cfs, c(1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -1.6946, -1.8863, -1.503, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -2.1071, -2.3076, -1.9065, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -1.7015, -1.8933, -1.5097, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -1.0456, -1.2293, -0.8618, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -0.3109, -0.4952, -0.1265, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -1.2417, -1.4271, -1.0564, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -1.681, -1.8724, -1.4895, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -2.122, -2.3229, -1.9211, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, -2.0259, -2.2244, -1.8274, 0, NA, NA, 1, NA, NA, 1.1116, 0.9623, 1.2609, 2.2465, 2.0999, 2.3931, 1, NA, NA, 1.6169, 1.3448, 1.8889, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.1067, 0.0791, 0.1343),
                 tollerance = 1e-2)
    
    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, 
                                        itemtype = '2PL', verbose = FALSE, draws = 10))
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 293) 
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, 0.0182, -0.2921, 0.3285, -1.7725, -1.9654, -1.5797, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, 0.1948, -0.0668, 0.4564, -2.1973, -2.4035, -1.9911, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, 0.098, -0.1065, 0.3024, -1.7817, -1.9755, -1.588, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, -2.4816, -2.6193, -2.344, -1.0787, -1.3685, -0.7889, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, 0.055, -0.2359, 0.3459, -0.3755, -0.5584, -0.1925, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, 0.0871, -0.2021, 0.3762, -1.3166, -1.5028, -1.1304, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, 0.004, -0.2585, 0.2665, -1.7584, -1.9509, -1.5659, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, 0.3122, 0.0081, 0.6163, -2.2245, -2.4329, -2.0161, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, -0.1009, -0.2916, 0.0897, -2.1025, -2.3019, -1.9031, 0, NA, NA, 1, NA, NA, 1.1944, 1.0501, 1.3387, 2.3817, 2.2332, 2.5303, -0.086, -0.4667, 0.2948, 1.5562, 1.2865, 1.8258, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- suppressMessages(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group, 
                                        draws = 10, verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 303) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1, NA, NA, -0.5536, -0.6897, -0.4176, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9713, -1.016, -0.9265, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.5605, -0.696, -0.4251, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.0982, -0.0309, 0.2273, 0, NA, NA, 1, NA, NA, 1, NA, NA, 0.838, 0.6664, 1.0096, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.0985, -0.2433, 0.0463, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.5399, -0.6771, -0.4026, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.9865, -1.0207, -0.9523, 0, NA, NA, 1, NA, NA, 1, NA, NA, -0.8887, -0.9679, -0.8096, 0, NA, NA, 1, NA, NA, 1, NA, NA, 2.8409, 2.2759, 3.4059, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0.3729, 0.3174, 0.4285, 0.6345, 0.2372),
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
    expect_equal(cfs, c(-0.0622, -0.1116, -0.0128, 0.2288, 0.1793, 0.2783, 1.5037, 1.3831, 1.6243, 1, 0, 0, 1, -0.0622, -0.1116, -0.0128, 0.2288, 0.1793, 0.2783, 1.5037),
                 tollerance = 1e-2)
    
    sv2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                     itemdesign = itemdesign, pars='values'))
    expect_is(sv2, 'data.frame')       
    LLTM2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                       itemdesign = itemdesign, verbose = FALSE, draws = 10))
    expect_is(LLTM2, 'MixedClass')
    expect_equal(LLTM@df - LLTM2@df, 1)
    cfs <- na.omit(as.numeric(do.call(c, coef(LLTM2, digits=4))))[1:20]
    expect_equal(cfs, c(-0.0623, -0.1122, -0.0125, 0.2903, 0.2277, 0.3529, 1.5023, 1.3559, 1.6486, 0.0031, -0.197, 0.2031, 1, 0, 0, 1, -0.0623, -0.1122, -0.0125, 0.2903),
                 tollerance = 1e-2)
}) 

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    suppressWarnings(mod <- mixedmirt(Science, covdat, model=model,
                                       fixed = ~ 0 + group, verbose = FALSE, draws = 10))
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod, digits=4)))
    expect_equal(cfs, c(-0.0444, -0.2772, 0.1884, 1, NA, NA, 0, NA, NA, 3.0618, 2.0323, 4.0912, 5.6507, 4.5894, 6.7121, 4.285, 3.179, 5.3909, -0.0444, -0.2772, 0.1884, 1, NA, NA, 0, NA, NA, 1.8823, 1.4152, 2.3493, 2.7962, 2.2531, 3.3392, 0.9698, 0.296, 1.6436, -0.0444, -0.2772, 0.1884, 1, NA, NA, 0, NA, NA, 2.6289, 1.9809, 3.2769, 4.0474, 3.3341, 4.7606, 2.9364, 2.1496, 3.7233, -0.0444, -0.2772, 0.1884, 1, NA, NA, 0, NA, NA, 2.4333, 1.8907, 2.9758, 3.3336, 2.7157, 3.9516, 2.004, 1.2906, 2.7174, 0, NA, NA, 0.9549, 0.5823, 1.3276),
                 tollerance = 1e-2)
    
    suppressWarnings(mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE))    
    expect_is(mod2, 'MixedClass')   
    expect_equal(mod@df - mod2@df, 3)   
    cfs <- as.numeric(do.call(c, coef(mod2, digits=4)))
    expect_equal(cfs, c(-0.1598, -0.3777, 0.0581, 0.8235, 0.2957, 1.3513, 0, NA, NA, 2.8221, 1.5391, 4.105, 5.3571, 3.8775, 6.8367, 4.1249, 2.7923, 5.4575, -0.1598, -0.3777, 0.0581, 0.8438, 0.5818, 1.1058, 0, NA, NA, 1.7845, 1.2541, 2.3149, 2.7255, 2.0982, 3.3528, 1.0697, 0.3939, 1.7455, -0.1598, -0.3777, 0.0581, 2.5755, 0.2193, 4.9317, 0, NA, NA, 5.2787, 1.1207, 9.4366, 7.7333, 1.7871, 13.6794, 5.7351, 1.1257, 10.3445, -0.1598, -0.3777, 0.0581, 0.6902, 0.3089, 1.0714, 0, NA, NA, 2.1313, 1.5121, 2.7505, 2.9944, 2.304, 3.6848, 1.9136, 1.2496, 2.5777, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE))    
    expect_is(mod3, 'MixedClass')   
    expect_equal(mod3@df, 72) 
    cfs <- as.numeric(do.call(c, coef(mod3, digits=4)))
    expect_equal(cfs, c(-0.1994, -0.5082, 0.1094, 0.9936, 0.7801, 1.2071, 4.9094, 3.9783, 5.8404, 2.6999, 2.2705, 3.1293, -1.3551, -1.543, -1.1673, -0.1994, -0.5082, 0.1094, 1.2207, 0.929, 1.5123, 3.0102, 2.6179, 3.4025, 0.9901, 0.7245, 1.2558, -2.1718, -2.5892, -1.7545, -0.1994, -0.5082, 0.1094, 2.5302, 0.9486, 4.1117, 5.6422, 3.095, 8.1894, 2.4468, 1.1376, 3.756, -2.011, -2.4075, -1.6144, -0.1994, -0.5082, 0.1094, 1.0506, 0.8756, 1.2257, 3.4052, 2.875, 3.9353, 1.0706, 0.7551, 1.386, -1.5746, -1.752, -1.3971, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- suppressMessages(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 72) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(0.9999, 0.6502, 1.3496, 4.825, 3.8768, 5.7732, 2.6141, 2.1908, 3.0374, -1.4432, -1.7709, -1.1156, 1.2486, 0.8766, 1.6206, 2.9426, 2.4399, 3.4453, 0.9099, 0.602, 1.2178, -2.2728, -2.6634, -1.8821, 2.525, 1.4283, 3.6218, 5.5601, 3.8662, 7.254, 2.3663, 1.5014, 3.2313, -2.0768, -2.7211, -1.4325, 1.0532, 0.7158, 1.3906, 3.3167, 2.7905, 3.843, 0.9839, 0.7125, 1.2552, -1.6614, -2.017, -1.3058, 0, NA, NA, 1, NA, NA, 1e-04, 1e-04),
                 tollerance = 1e-2)
    
    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)
    
}) 
