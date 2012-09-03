context('multipleGroup')

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
    models <- confmirt.model('confmods/MGmodel1', quiet = TRUE)
    
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE)
    expect_is(mod_configural, 'MultipleGroupClass')
    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), verbose = FALSE)
    expect_is(mod_metric, 'MultipleGroupClass')
    mod_scalar2 <- multipleGroup(dat, models, group = group, verbose = FALSE, 
                                 invariance=c('slopes', 'intercepts', 'free_varcov','free_means'))
    expect_is(mod_scalar2, 'MultipleGroupClass')
    mod_scalar1 <- multipleGroup(dat, models, group = group, verbose = FALSE,
                                 invariance=c('slopes', 'intercepts', 'free_varcov'))    
    expect_is(mod_scalar1, 'MultipleGroupClass')
    
    
})

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
    
    #group models
    model1 <- confmirt.model('confmods/MGmodelg1', quiet = TRUE)    
    model2 <- confmirt.model('confmods/MGmodelg2', quiet = TRUE)    
    models <- list(D1=model1, D2=model2)
    
    mod_configural <- multipleGroup(dat, models, group = group, verbose = FALSE)
    expect_is(mod_configural, 'MultipleGroupClass')
    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), verbose = FALSE)
    expect_is(mod_metric, 'MultipleGroupClass')    
    mod_scalar1 <- multipleGroup(dat, models, group = group, verbose = FALSE,
                                 invariance=c('slopes', 'intercepts', 'free_varcov'))    
    expect_is(mod_scalar1, 'MultipleGroupClass')
    
})


