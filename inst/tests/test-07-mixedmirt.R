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
    expect_equal(cfs, c(1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -1.6946, 0.0978, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -2.1071, 0.1023, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -1.7015, 0.0979, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -1.0456, 0.0937, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -0.3109, 0.0941, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -1.2417, 0.0946, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -1.681, 0.0977, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -2.122, 0.1025, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, -2.0259, 0.1013, 0, NA, 1, NA, 1.1116, 0.0762, 2.2465, 0.0748, 1, NA, 1.6169, 0.1388, 0, NA, 1, NA, 0, NA, 0.1067, 0.0141),
                 tollerance = 1e-2)
    
    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, 
                                        itemtype = '2PL', verbose = FALSE, draws = 10))
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 293) 
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1944, 0.0736, 2.3817, 0.0758, 0.0182, 0.1583, -1.7725, 0.0984, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, 0.1948, 0.1335, -2.1973, 0.1052, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, 0.098, 0.1043, -1.7817, 0.0989, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, -2.4816, 0.0702, -1.0787, 0.1479, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, 0.055, 0.1484, -0.3755, 0.0933, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, 0.0871, 0.1475, -1.3166, 0.095, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, 0.004, 0.1339, -1.7584, 0.0982, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, 0.3122, 0.1552, -2.2245, 0.1063, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, -0.1009, 0.0973, -2.1025, 0.1018, 0, NA, 1, NA, 1.1944, 0.0736, 2.3817, 0.0758, -0.086, 0.1943, 1.5562, 0.1376, 0, NA, 1, NA, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- suppressMessages(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group, 
                                        draws = 10, verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 303) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1, NA, -0.5536, 0.0694, 0, NA, 1, NA, 1, NA, -0.9713, 0.0228, 0, NA, 1, NA, 1, NA, -0.5605, 0.0691, 0, NA, 1, NA, 1, NA, 0.0982, 0.0659, 0, NA, 1, NA, 1, NA, 0.838, 0.0875, 0, NA, 1, NA, 1, NA, -0.0985, 0.0739, 0, NA, 1, NA, 1, NA, -0.5399, 0.07, 0, NA, 1, NA, 1, NA, -0.9865, 0.0174, 0, NA, 1, NA, 1, NA, -0.8887, 0.0404, 0, NA, 1, NA, 1, NA, 2.8409, 0.2883, 0, NA, 1, NA, 0, NA, 0.3729, 0.0283, 0.6345, 0.2372),
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
    expect_equal(cfs, c(-0.0622, 0.0252, 0.2288, 0.0252, 1.5037, 0.0615, 1, 0, 0, 1, -0.0622, 0.0252, 0.2288, 0.0252, 1.5037, 0.0615, 1, 0, 0, 1),
                 tollerance = 1e-2)
    
    sv2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                     itemdesign = itemdesign, pars='values'))
    expect_is(sv2, 'data.frame')       
    LLTM2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                       itemdesign = itemdesign, verbose = FALSE, draws = 10))
    expect_is(LLTM2, 'MixedClass')
    expect_equal(LLTM@df - LLTM2@df, 1)
    cfs <- na.omit(as.numeric(do.call(c, coef(LLTM2, digits=4))))[1:20]
    expect_equal(cfs, c(-0.0623, 0.0254, 0.2903, 0.0319, 1.5023, 0.0747, 0.0031, 0.1021, 1, 0, 0, 1, -0.0623, 0.0254, 0.2903, 0.0319, 1.5023, 0.0747, 0.0031, 0.1021),
                 tollerance = 1e-2)
}) 

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    suppressWarnings(mod <- mixedmirt(Science, covdat, model=model,
                                       fixed = ~ 0 + group, verbose = FALSE, draws = 10))
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod, digits=4)))
    expect_equal(cfs, c(-0.0444, 0.1188, 1, NA, 0, NA, 3.0618, 0.5253, 5.6507, 0.5415, 4.285, 0.5643, -0.0444, 0.1188, 1, NA, 0, NA, 1.8823, 0.2383, 2.7962, 0.2771, 0.9698, 0.3438, -0.0444, 0.1188, 1, NA, 0, NA, 2.6289, 0.3306, 4.0474, 0.3639, 2.9364, 0.4015, -0.0444, 0.1188, 1, NA, 0, NA, 2.4333, 0.2768, 3.3336, 0.3153, 2.004, 0.364, 0, NA, 0.9549, 0.1901),
                 tollerance = 1e-2)
    
    suppressWarnings(mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE))    
    expect_is(mod2, 'MixedClass')   
    expect_equal(mod@df - mod2@df, 3)   
    cfs <- as.numeric(do.call(c, coef(mod2, digits=4)))
    expect_equal(cfs, c(-0.1598, 0.1112, 0.8235, 0.2693, 0, NA, 2.8221, 0.6546, 5.3571, 0.7549, 4.1249, 0.6799, -0.1598, 0.1112, 0.8438, 0.1337, 0, NA, 1.7845, 0.2706, 2.7255, 0.3201, 1.0697, 0.3448, -0.1598, 0.1112, 2.5755, 1.2022, 0, NA, 5.2787, 2.1214, 7.7333, 3.0338, 5.7351, 2.3518, -0.1598, 0.1112, 0.6902, 0.1945, 0, NA, 2.1313, 0.3159, 2.9944, 0.3523, 1.9136, 0.3388, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE))    
    expect_is(mod3, 'MixedClass')   
    expect_equal(mod3@df, 72) 
    cfs <- as.numeric(do.call(c, coef(mod3, digits=4)))
    expect_equal(cfs, c(-0.1994, 0.1575, 0.9936, 0.1089, 4.9094, 0.475, 2.6999, 0.2191, -1.3551, 0.0958, -0.1994, 0.1575, 1.2207, 0.1488, 3.0102, 0.2002, 0.9901, 0.1355, -2.1718, 0.2129, -0.1994, 0.1575, 2.5302, 0.8069, 5.6422, 1.2996, 2.4468, 0.668, -2.011, 0.2023, -0.1994, 0.1575, 1.0506, 0.0893, 3.4052, 0.2705, 1.0706, 0.1609, -1.5746, 0.0905, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- suppressMessages(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 72) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(0.9999, 0.1784, 4.825, 0.4838, 2.6141, 0.216, -1.4432, 0.1672, 1.2486, 0.1898, 2.9426, 0.2565, 0.9099, 0.1571, -2.2728, 0.1993, 2.525, 0.5596, 5.5601, 0.8642, 2.3663, 0.4413, -2.0768, 0.3287, 1.0532, 0.1721, 3.3167, 0.2685, 0.9839, 0.1384, -1.6614, 0.1814, 0, NA, 1, NA, 1e-04, 1e-04),
                 tollerance = 1e-2)
    
    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)
    
}) 
