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
    expect_equal(cfs, c(1.2194, 0.0771, 2.4284, 0.0819, 0.0223, 0.132, -1.8006, 0.1008, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, 0.0998, 0.1398, -2.2197, 0.1058, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, 0.0239, 0.1075, -1.8075, 0.1009, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, -3.7769, 1.4659, -1.062, 0.2292, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, 0.0333, 0.1118, -0.3949, 0.0946, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, 0.0534, 0.1591, -1.3411, 0.0971, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, -0.0495, 0.1242, -1.7863, 0.1004, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, 0.2393, 0.1574, -2.2432, 0.1066, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, -0.1436, 0.137, -2.1353, 0.1042, 0, NA, 1, NA, 1.2194, 0.0771, 2.4284, 0.0819, -0.125, 0.2006, 1.5501, 0.1381, 0, NA, 1, NA, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- suppressMessages(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group, 
                                        draws = 10, verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 303) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1, NA, -0.5571, 0.0747, 0, NA, 1, NA, 1, NA, -0.9704, 0.0706, 0, NA, 1, NA, 1, NA, -0.5639, 0.0747, 0, NA, 1, NA, 1, NA, 0.0885, 0.0728, 0, NA, 1, NA, 1, NA, 0.8219, 0.0507, 0, NA, 1, NA, 1, NA, -0.1064, 0.0747, 0, NA, 1, NA, 1, NA, -0.5435, 0.0748, 0, NA, 1, NA, 1, NA, -0.9855, 0.0704, 0, NA, 1, NA, 1, NA, -0.8888, 0.0717, 0, NA, 1, NA, 1, NA, 2.8103, 0.0961, 0, NA, 1, NA, 0, NA, 0.3497, 0.0126, 0.5906, 0.2777),
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
    expect_equal(cfs, c(-0.1563, 0.1079, 0.822, 0.2338, 0, NA, 2.8092, 0.6237, 5.3378, 0.7025, 4.0997, 0.6456, -0.1563, 0.1079, 0.8452, 0.1357, 0, NA, 1.7802, 0.268, 2.7163, 0.314, 1.0521, 0.341, -0.1563, 0.1079, 2.5808, 1.0682, 0, NA, 5.2668, 1.8736, 7.7165, 2.6715, 5.6953, 2.0968, -0.1563, 0.1079, 0.6857, 0.1831, 0, NA, 2.1206, 0.316, 2.9776, 0.3548, 1.8936, 0.3436, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE))    
    expect_is(mod3, 'MixedClass')   
    expect_equal(mod3@df, 72) 
    cfs <- as.numeric(do.call(c, coef(mod3, digits=4)))
    expect_equal(cfs, c(-0.187, 0.1597, 1.0007, 0.1825, 4.9121, 0.4907, 2.7, 0.2279, -1.3589, 0.1982, -0.187, 0.1597, 1.2157, 0.1985, 3.0005, 0.2589, 0.985, 0.176, -2.1672, 0.2286, -0.187, 0.1597, 2.5263, 0.7053, 5.6266, 1.0828, 2.4402, 0.5751, -1.9992, 0.3701, -0.187, 0.1597, 1.07, 0.1718, 3.4151, 0.2832, 1.0721, 0.1672, -1.5855, 0.1997, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- suppressMessages(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE))
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 72) 
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1.0028, 0.1687, 4.827, 0.4797, 2.6162, 0.2109, -1.4461, 0.1523, 1.2369, 0.1639, 2.9335, 0.241, 0.905, 0.1453, -2.2685, 0.1903, 2.5053, 0.4909, 5.539, 0.742, 2.3496, 0.3533, -2.0704, 0.349, 1.0646, 0.1701, 3.325, 0.269, 0.9861, 0.1363, -1.669, 0.1635, 0, NA, 1, NA, 1e-04, 1e-04),
                 tollerance = 1e-2)
    
    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)
    
}) 
