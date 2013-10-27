context('multipleGroupTwo')

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
    
    suppressWarnings(mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), method = 'MHRM',
                                verbose = FALSE, draws = 10))
    expect_is(mod_metric, 'MultipleGroupClass') 
    cfs <- as.numeric(do.call(c, coef(mod_metric, digits=4)[[1]]))[1:20]
    expect_equal(cfs, c(1.3033, 0.8185, 1.788, 0, NA, NA, 0, NA, NA, 0.6597, 0.4677, 0.8516, 0, NA, NA, 1, NA, NA, 1.2852, 0.7614),
                 tollerance = 1e-2)
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural, digits=4)[[1]]))
    cfs <- cfs[cfs != 0 & cfs != 1]    
    expect_equal(cfs, c(1.2934, 0.655, 1.2387, -0.57, 0.9276, -0.1996, 0.8177, 0.7967, 1.0713, 0.2166, 0.4807, 0.61, 1.1778, 0.9948, 0.9453, -0.4464, 1.0761, -1.18, 0.8664, -1.1451, 0.8854, 1.3121, 1.4964, -0.3002, 1.0534, 0.4406, 1.0614, 0.4569, 0.8831, -0.1871),
                 tollerance = 1e-2)
    suppressWarnings(mod_scalar1 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'MHRM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov', draws = 10)))
    expect_is(mod_scalar1, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_scalar1, digits=4)[[1]]))[1:20]
    expect_equal(cfs, c(1.1907, 0.6966, 1.6848, 0, NA, NA, 0, NA, NA, 0.3851, 0.0887, 0.6815, 0, NA, NA, 1, NA, NA, 1.1926, 0.6811),
                 tollerance = 1e-2)
    
    fs1 <- fscores(mod_metric, verbose = FALSE)
    fs2 <- fscores(mod_scalar1, full.scores = TRUE)    
    expect_is(fs1, 'list')
    expect_is(fs2, 'data.frame')    
    expect_true(mirt:::closeEnough(fs2[1:6, 'F1'] - c(0.8076989,  0.8076989,  0.2915864, -0.1365438, -0.2653481,  0.8076989), -1e-4, 1e-4))     
})


