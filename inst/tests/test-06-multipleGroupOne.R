context('multipleGroupOne')

test_that('one factor', {    
    set.seed(12345)
    a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
    d <- matrix(rnorm(15,0,.7),ncol=1)    
    itemtype <- rep('dich', nrow(a))
    N <- 1000    
    dataset1 <- simdata(a, d, N, itemtype)
    dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
    dat <- rbind(dataset1, dataset2)
    group <- c(rep('D1', N), rep('D2', N))    
    MGmodel1 <- 'F1 = 1-15'    
    models <- confmirt.model(MGmodel1, quiet = TRUE)
    
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]    
    expect_equal(cfs, c(1.069,  0.554,  1.278, -0.692,  0.883, -0.138,  1.111,  0.829,  1.248,  
                        0.326,  0.476,  0.480,  1.162,  1.085,  0.859, -0.385,  0.890,
                       -1.048,  0.809, -1.091,  0.901,  1.164,  1.583, -0.135,  1.410,  0.654,  
                        1.040,  0.407,  0.880, -0.081), tollerance = 1e-3)
    expect_equal(mod_configural@df, 1621)
    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), verbose = FALSE, 
                                method = 'EM')
    expect_is(mod_metric, 'MultipleGroupClass')
    expect_equal(mod_metric@df, 1636)
    mod_scalar2 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov','free_means'))
    cfs <- as.numeric(do.call(c, coef(mod_scalar2)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]    
    expect_equal(cfs, c(1.142,  0.563,  1.325, -0.651,  0.993, -0.201,  1.049,  0.887,  1.145, 
                        0.338,  0.431,  0.497,  1.226,  1.158,  0.916, -0.420,  0.816,
                       -1.016,  0.801, -1.089,  0.948,  1.235,  1.588, -0.189,  1.199,  0.539, 
                        1.129,  0.433,  0.893, -0.117), tollerance = 1e-3)
    expect_is(mod_scalar2, 'MultipleGroupClass')
    expect_equal(mod_scalar2@df, 1649)
    mod_scalar1 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'MHRM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov'), draws = 10)    
    expect_is(mod_scalar1, 'MultipleGroupClass')    
    dat[1,1] <- dat[2,2] <- NA
    mod_missing <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov'))    
    expect_is(mod_missing, 'MultipleGroupClass')
    expect_equal(mod_missing@df, 1651)
    
    fs1 <- fscores(mod_metric, verbose = FALSE)
    expect_true(mirt:::closeEnough(fs1[[1]][1:6, 'F1'] - c(-2.084760, -1.683841, -1.412181,
                                                           -1.656879, -1.324689, -1.092169), -1e-4, 1e-4))    
    fs2 <- fscores(mod_metric, full.scores = TRUE)
    fs3 <- fscores(mod_missing, verbose = FALSE)
    fs4 <- fscores(mod_missing, full.scores = TRUE)
    fs5 <- fscores(mod_metric, full.scores = TRUE, scores.only=TRUE)
    expect_is(fs1, 'list')
    expect_is(fs2, 'data.frame')
    expect_is(fs3, 'list')
    expect_is(fs4, 'data.frame') 
    
    fit1 <- fitIndices(mod_metric)
    expect_is(fit1, 'list')
    expect_true(mirt:::closeEnough(fit1$M2 - c(85.70548, 120.79436), -1e-4, 1e-4))
    expect_true(mirt:::closeEnough(fit1$df - 350, -1e-4, 1e-4))    
    fit2 <- itemfit(mod_metric)
    expect_is(fit2, 'list')

    #missing data
    set.seed(1234)
    Theta1 <- rnorm(1000, -1)
    Theta2 <- rnorm(1000, 1) 
    Theta <- matrix(rbind(Theta1, Theta2))
    d <- rnorm(10,4)
    d <- cbind(d, d-1, d-2, d-3, d-4, d-5, d-6)
    a <- matrix(rlnorm(10, meanlog=.1))
    group <- factor(c(rep('g1',1000), rep('g2',1000)))
    
    dat <- simdata(a,d,2000, itemtype = rep('graded', 10), Theta=Theta)
    x <- multipleGroup(dat, 1, group=group, method='EM', verbose = FALSE)
    expect_is(x, 'MultipleGroupClass')
    
    dat[1,1] <- dat[2,2] <- NA
    x2 <- multipleGroup(dat, 1, group=group, method='EM', verbose = FALSE)
    expect_is(x2, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(x2)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]    
    expect_equal(cfs, c(
        0.676,  2.876,  1.895,  0.987,  0.117, -0.953, -2.026, -3.151,  
        0.616,  3.973,  2.866,  1.701,  0.777, -0.112, -1.064, -2.124,  2.119,
        4.006,  2.920,  1.872,  0.909, -0.151, -1.116, -2.070,  2.752,
        5.350,  4.397,  3.272,  2.373,  1.261,  0.281, -0.839,  0.466,  2.405,
        1.424,  0.457, -0.590, -1.634, -2.690, -3.686,  4.890,  3.055, 
        2.019,  1.107,  0.067, -0.902, -1.852, -3.079,  2.495,  2.332,  1.278,
        0.412, -0.509, -1.505, -2.531, -3.645,  1.965,  4.431,  3.470,
        2.441,  1.451,  0.516, -0.386, -1.438,  2.031,  3.436,  2.687,  1.592,
        0.661, -0.423, -1.428, -2.344,  2.354,  2.264,  1.250,  0.376, 
        -0.653, -1.814, -2.809, -3.880), tollerance = 1e-3)
})