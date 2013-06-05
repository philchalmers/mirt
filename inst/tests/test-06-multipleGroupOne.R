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
    expect_equal(mod_configural@df, 1166)
    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), verbose = FALSE, 
                                method = 'EM')
    expect_is(mod_metric, 'MultipleGroupClass')
    expect_equal(mod_metric@df, 1181)
    mod_scalar2 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov','free_means'))
    expect_is(mod_scalar2, 'MultipleGroupClass')
    expect_equal(mod_scalar2@df, 1194)
    mod_scalar1 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'MHRM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov'), draws = 10)    
    expect_is(mod_scalar1, 'MultipleGroupClass')    
    dat[1,1] <- dat[2,2] <- NA
    mod_missing <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_varcov'))    
    expect_is(mod_missing, 'MultipleGroupClass')
    expect_equal(mod_missing@df, 1195)
    
    fs1 <- fscores(mod_metric, verbose = FALSE)
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
    expect_is(x, 'MultipleGroupClass')
})