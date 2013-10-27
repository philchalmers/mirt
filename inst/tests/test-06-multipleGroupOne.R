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
    models <- mirt.model(MGmodel1, quiet = TRUE)
    
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural, digits=4)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]    
    expect_equal(cfs, c(1.0693, 0.5541, 1.278, -0.6918, 0.8833, -0.1375, 1.1112, 0.8295, 1.2481, 0.3265, 0.476, 0.4796, 1.1617, 1.0847, 0.8586, -0.3852, 0.89, -1.048, 0.8085, -1.0908, 0.9013, 1.1642, 1.5832, -0.135, 1.4098, 0.6542, 1.0401, 0.4073, 0.8804, -0.0812),
                 tollerance = 1e-2)
    expect_equal(mod_configural@df, 1621)
    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), verbose = FALSE, 
                                method = 'EM')
    expect_is(mod_metric, 'MultipleGroupClass')
    expect_equal(mod_metric@df, 1636)
    mod_scalar2 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov','free_means'))
    cfs <- as.numeric(do.call(c, coef(mod_scalar2, digits=4)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]    
    expect_equal(cfs, c(1.1424, 0.5623, 1.3257, -0.6508, 0.9936, -0.2008, 1.0489, 0.8867, 1.1449, 0.3383, 0.4314, 0.4965, 1.2256, 1.158, 0.916, -0.4197, 0.8163, -1.0164, 0.8011, -1.0888, 0.9486, 1.2348, 1.5887, -0.1893, 1.1991, 0.5387, 1.1291, 0.4329, 0.8934, -0.117),
                 tollerance = 1e-2)
    expect_is(mod_scalar2, 'MultipleGroupClass')
    expect_equal(mod_scalar2@df, 1649)
    newmodel <- mirt.model('F = 1-15
                            CONSTRAINB = (1-15, a1), (1,2,3-15,d)')
    mod_scalar1 <- multipleGroup(dat, newmodel, group = group, verbose = FALSE, invariance='free_var')   
    expect_is(mod_scalar1, 'MultipleGroupClass')  
    
    
    dat[1,1] <- dat[2,2] <- NA
    mod_missing <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov'))    
    expect_is(mod_missing, 'MultipleGroupClass')
    expect_equal(mod_missing@df, 1651)
    
    fs1 <- fscores(mod_metric, verbose = FALSE)
    expect_true(mirt:::closeEnough(fs1[[1]][1:6, 'F1'] - c(-2.084760, -1.683841, -1.412181,
                                                           -1.656879, -1.324689, -1.092169), -1e-2, 1e-2))    
    fs2 <- fscores(mod_metric, full.scores = TRUE)
    fs3 <- fscores(mod_missing, verbose = FALSE)
    fs4 <- fscores(mod_missing, full.scores = TRUE)
    fs5 <- fscores(mod_metric, full.scores = TRUE, scores.only=TRUE)
    expect_is(fs1, 'list')
    expect_is(fs2, 'data.frame')
    expect_is(fs3, 'list')
    expect_is(fs4, 'data.frame') 
    
    fit1 <- fitIndices(mod_metric)
    expect_is(fit1, 'data.frame')
    expect_true(mirt:::closeEnough(fit1[1:2] - c(1126.54, 2162.74), -1e-2, 1e-2))
    expect_true(mirt:::closeEnough(fit1$df.M2 - 350, -1e-4, 1e-4))    
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
    cfs <- as.numeric(do.call(c, coef(x2, digits = 5)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]    
    expect_true(mirt:::closeEnough(cfs - c(0.67643, 2.87461, 1.89346, 0.98554, 0.11572, -0.95404, -2.02763, -3.15235, 0.61587, 3.97145, 2.86507, 1.69939, 0.77596, -0.11359, -1.06517, -2.12574, 2.11956, 4.00101, 2.91579, 1.86772, 0.90442, -0.15526, -1.12035, -2.07435, 2.75164, 5.34259, 4.3898, 3.26498, 2.36583, 1.25436, 0.2745, -0.84438, 0.46576, 2.40369, 1.42319, 0.45604, -0.59102, -1.63473, -2.69143, -3.68742, 4.89341, 3.04587, 2.00935, 1.09628, 0.05593, -0.91307, -1.86358, -3.09127, 2.49587, 2.32619, 1.27206, 0.40676, -0.51432, -1.51021, -2.53623, -3.65075, 1.9649, 4.42635, 3.4659, 2.43702, 1.44714, 0.51205, -0.38967, -1.44249, 2.0317, 3.43162, 2.68263, 1.58805, 0.65621, -0.42705, -1.43263, -2.34819, 2.35415, 2.25881, 1.24474, 0.37074, -0.6577, -1.81886, -2.8138, -3.88558), -1e-2, 1e-2))
    
})
