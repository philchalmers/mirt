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
    model <- confmirt.model(mixedmirt1, quiet = TRUE)  
    
    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, 
                                       verbose = FALSE, draws = 10))
    expect_is(mod1, 'MixedClass')  
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -1.6974, 0.098, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -2.1101, 0.1026, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -1.7042, 0.0981, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -1.0478, 0.0937, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -0.3125, 0.0938, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -1.2441, 0.0947, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -1.6837, 0.0979, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -2.1251, 0.1028, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, -2.0289, 0.1016, 0, NA, 1, NA, 
                        1.1181, 0.0742, 2.2515, 0.0792, 1, NA, 1.6166, 0.1383, 0, NA, 1, NA, 0, 
                        NA, 0.1108, 0.0129),
                 tollerance = 1e-2)
    
    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- suppressWarnings(mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, 
                                        itemtype = '2PL', verbose = FALSE, draws = 10))
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 293) 
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.2178, 0.0749, 2.4259, 0.0768, 0.0153, 0.1275, -1.799, 0.0991, 0, NA, 1, NA,
                        1.2178, 0.0749, 2.4259, 0.0768, 0.1055, 0.134, -2.2187, 0.1045, 0, NA, 1, NA,
                        1.2178, 0.0749, 2.4259, 0.0768, 0.0269, 0.1054, -1.8061, 0.0992, 0, NA, 1, NA, 
                        1.2178, 0.0749, 2.4259, 0.0768, -3.5963, 0.421, -1.0678, 0.2463, 0, NA, 1, NA,
                        1.2178, 0.0749, 2.4259, 0.0768, 0.0266, 0.1191, -0.394, 0.0936, 0, NA, 1, NA,
                        1.2178, 0.0749, 2.4259, 0.0768, 0.0562, 0.1098, -1.3398, 0.0955, 0, NA, 1, NA,
                        1.2178, 0.0749, 2.4259, 0.0768, -0.0605, 0.1402, -1.7845, 0.0988, 0, NA, 1, NA,
                        1.2178, 0.0749, 2.4259, 0.0768, 0.2487, 0.1416, -2.2434, 0.1057, 0, NA, 1, NA,
                        1.2178, 0.0749, 2.4259, 0.0768, -0.1546, 0.1194, -2.1332, 0.1029, 0, NA, 1, NA,
                        1.2178, 0.0749, 2.4259, 0.0768, -0.1044, 0.1937, 1.5493, 0.138, 0, NA, 1, NA,
                        0, NA, 1, NA),
                 tollerance = 1e-2)
    
    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- suppressMessages(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group, 
                                        draws = 10, verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 303) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1, NA, -0.5526, 0.0745, 0, NA, 1, NA, 1, NA, -0.9644, 0.0808, 0, NA, 1, NA,
                        1, NA, -0.5594, 0.0746, 0, NA, 1, NA, 1, NA, 0.0897, 0.0552, 0, NA, 1, NA, 
                        1, NA, 0.8181, 0.0283, 0, NA, 1, NA, 1, NA, -0.1041, 0.0627, 0, NA, 1, NA, 
                        1, NA, -0.539, 0.0742, 0, NA, 1, NA, 1, NA, -0.9794, 0.0809, 0, NA, 1, NA, 
                        1, NA, -0.883, 0.0798, 0, NA, 1, NA, 1, NA, 2.7863, 0.0376, 0, NA, 1, NA, 0, 
                        NA, 0.2861, 0.0158, 0.6167, 0.323),
                 tollerance = 1e-2)
})   

test_that('item and group predictors', {    
    data <- key2binary(SAT12,
                       key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
    model <- confmirt.model('Theta = 1-32', quiet = TRUE)
    
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
    expect_equal(cfs, c(-0.0593,  0.0225,  0.2315,  0.0229,  1.4996,  0.0592,  1.0000,  0.0000, 
                        0.0000,  1.0000, -0.0593,  0.0225,  0.2315,  0.0229,  1.4996,
                        0.0592,  1.0000,  0.0000,  0.0000,  1.0000),
                 tollerance = 1e-2)
    
    sv2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                     itemdesign = itemdesign, pars='values'))
    expect_is(sv2, 'data.frame')       
    LLTM2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender,
                       itemdesign = itemdesign, verbose = FALSE, draws = 10))
    expect_is(LLTM2, 'MixedClass')
    expect_equal(LLTM@df - LLTM2@df, 1)
    cfs <- na.omit(as.numeric(do.call(c, coef(LLTM2, digits=4))))[1:20]
    expect_equal(cfs, c(-0.0606,  0.0176,  0.2867,  0.0214,  1.5025,  0.0553,  0.0101,  0.0722,
                        1.0000,  0.0000,  0.0000,  1.0000, -0.0606,  0.0176,  0.2867,
                         0.0214,  1.5025,  0.0553,  0.0101,  0.0722),
                 tollerance = 1e-2)
}) 

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- confmirt.model('F1 = 1-4', quiet = TRUE)
    suppressWarnings(mod <- mixedmirt(Science, covdat, model=model,
                                       fixed = ~ 0 + group, verbose = FALSE, draws = 10))
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod, digits=4)))
    expect_equal(cfs, c(-0.0344, 0.2983, 1, NA, 0, NA, 3.062, 0.5685, 5.6596, 0.619, 4.2985, 0.6017,
                        -0.0344, 0.2983, 1, NA, 0, NA, 1.8906, 0.2755, 2.8129, 0.3183, 0.9892, 0.3568,
                        -0.0344, 0.2983, 1, NA, 0, NA, 2.6349, 0.3742, 4.0624, 0.4298, 2.956, 0.4266,
                        -0.0344, 0.2983, 1, NA, 0, NA, 2.4406, 0.3161, 3.3497, 0.3652, 2.0238, 0.3779,
                        0, NA, 0.9621, 0.3357),
                 tollerance = 1e-2)
    
    suppressWarnings(mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE))    
    expect_is(mod2, 'MixedClass')   
    expect_equal(mod@df - mod2@df, 3)   
    cfs <- as.numeric(do.call(c, coef(mod2, digits=4)))
    expect_equal(cfs, c(-0.1647, 0.1309, 0.8226, 0.1634, 0, NA, 2.8261, 0.566, 5.3655, 0.6229, 
                        4.1373, 0.6336, -0.1647, 0.1309, 0.8462, 0.167, 0, NA, 1.7922, 0.2704, 
                        2.7389, 0.3201, 1.0843, 0.3952, -0.1647, 0.1309, 2.5442, 1.1668, 0, NA, 
                        5.2303, 2.1811, 7.677, 3.1151, 5.7004, 2.5926, -0.1647, 0.1309, 0.6893, 
                        0.1406, 0, NA, 2.1351, 0.2946, 3.0024, 0.3438, 1.9254, 0.3951, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE))    
    expect_is(mod3, 'MixedClass')   
    expect_equal(mod3@df, 72) 
    cfs <- as.numeric(do.call(c, coef(mod3, digits=4)))
    expect_equal(cfs, c(-0.1946, 0.1112, 1.0097, 0.1616, 4.9276, 0.5265, 2.7051, 0.1617, -1.3612,
                        0.1532, -0.1946, 0.1112, 1.2311, 0.1682, 3.0152, 0.2178, 0.9913, 0.1198, 
                        -2.1766, 0.2151, -0.1946, 0.1112, 2.5012, 0.3056, 5.5945, 0.6328, 2.422, 
                        0.067, -1.9886, 0.2961, -0.1946, 0.1112, 1.0689, 0.1736, 3.416, 0.2762, 
                        1.0757, 0.1185, -1.5828, 0.1721, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- suppressMessages(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 72) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1.0235, 0.1758, 4.8438, 0.4848, 2.6275, 0.2173, -1.4587, 0.1566, 1.2404, 
                        0.1975, 2.9387, 0.251, 0.9072, 0.1487, -2.2741, 0.2114, 2.3052, 0.3405, 
                        5.2567, 0.5408, 2.2194, 0.2729, -1.9701, 0.2844, 1.072, 0.1553, 3.3285, 
                        0.2697, 0.9862, 0.141, -1.6771, 0.1598, 0, NA, 1, NA, 1e-04, 1e-04),
                 tollerance = 1e-2)
    
    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)
    
}) 
