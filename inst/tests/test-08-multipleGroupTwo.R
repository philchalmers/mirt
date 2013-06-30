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
    model1 <- confmirt.model(MGmodelg1, quiet = TRUE)    
    model2 <- confmirt.model(MGmodelg1, quiet = TRUE)    
    models <- list(D1=model1, D2=model2)    
    
    suppressWarnings(mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), method = 'MHRM',
                                verbose = FALSE, draws = 10))
    expect_is(mod_metric, 'MultipleGroupClass') 
    cfs <- as.numeric(do.call(c, coef(mod_metric, digits=4)[[1]]))[1:20]
    expect_equal(cfs, c(1.3003, 0.1238, 0, NA, 0, NA, 0.6606, 0.0731, 0, NA, 1, NA, 
                        1.2537, 0.0371, 0, NA, 0, NA, -0.5673, 0.0748),
                 tollerance = 1e-2)
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural)[[1]]))
    cfs <- cfs[cfs != 0 & cfs != 1]    
    expect_equal(cfs, c(1.297,  0.655,  1.242, -0.570,  0.931, -0.200,  0.821,  0.797,  1.075,  0.217,  0.483,  
                        0.610,  1.181,  0.995,  0.948, -0.446,  1.080, -1.180,  0.870, -1.145,  0.890,  1.312, 
                        1.500, -0.300,  1.057,  0.441,  1.065,  0.457,  0.886, -0.187),
                 tollerance = 1e-2)
    suppressWarnings(mod_scalar1 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'MHRM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov', draws = 10)))
    expect_is(mod_scalar1, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_scalar1, digits=4)[[1]]))[1:20]
    expect_equal(cfs, c(1.214, 0.0852, 0, NA, 0, NA, 0.3883, 0.0513, 0, NA, 1, NA, 1.1956, 0.0775, 
                        0, NA, 0, NA, -0.8127, 0.0521),
                 tollerance = 1e-2)
    
    fs1 <- fscores(mod_metric, verbose = FALSE)
    fs2 <- fscores(mod_scalar1, full.scores = TRUE)    
    expect_is(fs1, 'list')
    expect_is(fs2, 'data.frame')    
    expect_true(mirt:::closeEnough(fs2[1:6, 'F1'] - c(0.8035599,  0.8035599,  0.2855899, 
                                                      -0.1416255, -0.2674919,  0.8035599), -1e-4, 1e-4))     
})


