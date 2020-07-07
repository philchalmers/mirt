context('mixedmirt')

test_that('mixed dich', {
    if(FALSE){
        rm(list=ls())
        set.seed(1234)
        N <- 750
        a <- matrix(rlnorm(10,.2,.5),10,1)
        d <- matrix(rnorm(10), 10)
        Theta <- matrix(sort(rnorm(N)))
        pseudoIQ <- scale(Theta * 5 + 100  + rnorm(N, 0 , 5))
        group <- factor(rep(c('G1','G2','G3'), each = N/3))
        data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
        covdata <- data.frame(group, pseudoIQ)
        save(data, covdata, file = 'tests/tests/testdata/mixed1.rds')
    }
    load('testdata/mixed1.rds')

    mixedmirt1 <- 'Theta = 1-10'
    model <- mirt.model(mixedmirt1, quiet = TRUE)

    #simple latent regression
    mod0 <- mirt(data, 1, 'Rasch', covdata=covdata, formula = ~ group + pseudoIQ, verbose=FALSE)
    expect_equal(mod0@Fit$logLik, -4058.968, tolerance = 1e-2)
    cfs <- coef(mod0)
    expect_equal(as.numeric(cfs$lr.betas), c(0.0000000, 0.8916977, 1.9757340, 0.2168971), tolerance=1e-4)
    set.seed(1234)
    plaus <- fscores(mod0, plausible.draws = 2)
    expect_equal(plaus[[1]][1:4], c(-0.6875802, -0.3252686, -0.3528408, -0.4759154),
                 tolerance = 1e-4)
    require(boot, quietly=TRUE, warn.conflicts=FALSE)
    set.seed(1)
    bs <- boot.mirt(mod0, R = 3)
    expect_is(bs, 'boot')
    fs <- fscores(mod0, full.scores.SE=TRUE, full.scores=FALSE)
    expect_equal(as.numeric(head(fs)), c(-0.3326607,-0.6466616,-0.4794936,-0.3070341,-0.4769525,-0.5488866,0.2940349,0.2963674,0.2951423,0.2938395,0.2951234,0.2956562),
                 tolerance=1e-4)

    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                       verbose = FALSE, draws = 1)
    expect_is(mod1, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod1)))
    expect_equal(cfs, c(1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-1.747882,-1.944582,-1.551183,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-2.066577,-2.270076,-1.863078,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-1.705846,-1.901775,-1.509917,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-1.043203,-1.230739,-0.8556669,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-0.2660702,-0.4541086,-0.07803179,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-1.377461,-1.568343,-1.186579,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-1.684929,-1.880486,-1.489373,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-2.127009,-2.331998,-1.922021,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,-2.007019,-2.209114,-1.804923,0,NA,NA,1,NA,NA,1.066421,0.9144826,1.218359,2.336545,2.178417,2.494673,1,NA,NA,1.546882,1.281149,1.812614,0,NA,NA,1,NA,NA,0,NA,NA,0.1554333,0.1170689,0.1937978),
                 tolerance = 1e-2)
    names <- wald(mod1)
    L <- matrix(c(1, numeric(length(names) - 1L)), 1L)
    wld <- wald(mod1, L, C=as.numeric(L))
    expect_equal(wld$W[1], .7341243, tolerance = 1e-4)
    set.seed(1)
    bs <- boot.mirt(mod1, R=2)
    expect_is(bs, 'boot')

    mod1a <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, SE=FALSE,
                      verbose = FALSE, draws = 1, internal_constraints = FALSE)
    cfs <- as.numeric(do.call(c, coef(mod1a)))
    expect_equal(cfs, c(0.7949969,1.593086,1,-1.361286,0,1,1.391169,2.930032,1,-2.458063,0,1,2.423806,4.488466,1,-3.154002,0,1,0.2993587,0.7774742,1,-0.287956,0,1,1.479693,3.058504,1,-0.547135,0,1,1.194509,3.038936,1,-1.655208,0,1,0.9844124,1.969708,1,-1.515223,0,1,1.422937,2.655203,1,-2.406514,0,1,0.8457056,2.057593,1,-1.808905,0,1,0.8300673,2.352666,1,1.61897,0,1,0,0.1639394),
                 tolerance = 1e-2)

    mod_items <- mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items,
                       verbose = FALSE, draws = 1)
    cfs <- c(coef(mod_items)[['GroupPars']], coef(mod_items)[['items']])
    expect_equal(cfs[1:3], c(0.000, .9937, .8372), tolerance = 1e-2)

    mod_items.group <- try(mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items:group,
                           verbose = FALSE, draws = 1), TRUE)
    cfs <- c(coef(mod_items.group)[['GroupPars']], coef(mod_items.group)[['items:group']])
    expect_equal(cfs[1:3], c(0.000, 0.1598, 2.1794), tolerance = 1e-3)

    #model using 2PL items instead of only Rasch, and with missing data
    covdataold <- covdata
    data[1,1] <- covdata[1,2] <- NA
    expect_warning(mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 1))
    expect_is(mod1b, 'MixedClass')
    expect_equal(extract.mirt(mod1b, 'df'), 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b)))
    expect_equal(cfs, c(1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.6801064,0.2472736,1.112939,-1.764637,-1.981088,-1.548186,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.3225088,-0.01222681,0.6572444,-2.04966,-2.252809,-1.846511,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.1730519,-0.08809251,0.4341962,-1.68213,-1.873847,-1.490413,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,-0.6751006,NaN,NaN,-1.027004,-1.217833,-0.8361756,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.442664,-0.003793412,0.8891215,-0.2375153,-0.4341578,-0.04087278,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.1315495,-0.1747718,0.4378707,-1.357544,-1.544107,-1.170982,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.2533267,-0.1310559,0.6377093,-1.664991,-1.856788,-1.473194,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.3229516,-0.1145078,0.760411,-2.110101,-2.317428,-1.902774,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.7662229,0.01146026,1.520985,-2.051964,-2.303527,-1.800401,0,NA,NA,1,NA,NA,1.055344,0.910874,1.199815,2.313518,2.1583,2.468737,0.4567099,-0.08615558,0.9995754,1.610767,1.276957,1.944577,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)

    covdataold$group <- factor(rep(paste0('G',1:50), each = nrow(data)/50))
    rmod1 <- try(mixedmirt(data, covdataold, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 1, verbose = FALSE), TRUE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(extract.mirt(rmod1, 'df'), 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1)))

    #polytomous
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    mod <- mixedmirt(Science, covdat, model=model, SE=FALSE,
                     fixed = ~ 0 + group, verbose = FALSE, draws = 1)
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(na.omit(do.call(c, coef(mod))))
    expect_equal(cfs, c(-0.08424103,1,0,1,2,3,0,3.097676,5.718138,4.364193,-0.08424103,1,0,1,2,3,0,1.918539,2.858765,1.037761,-0.08424103,1,0,1,2,3,0,2.664839,4.113727,3.013447,-0.08424103,1,0,1,2,3,0,2.469392,3.397754,2.075902,0,0.9854896),
                 tolerance = 1e-2)

    mod2 <- mixedmirt(Science, covdat, model=model, draws = 1,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE)
    expect_is(mod2, 'MixedClass')
    expect_equal(extract.mirt(mod, 'df') - extract.mirt(mod2, 'df'), 3)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod2))))
    expect_equal(cfs, c(-0.1764055,-0.4078213,0.05501025,0.8245979,0.5139055,1.13529,0,1,2,3,0,2.838464,1.720788,3.95614,5.387634,4.161899,6.613368,4.164414,2.956008,5.37282,-0.1764055,-0.4078213,0.05501025,0.844655,0.5531942,1.136116,0,1,2,3,0,1.800049,1.265608,2.33449,2.753831,2.133251,3.37441,1.105474,0.4461877,1.764761,-0.1764055,-0.4078213,0.05501025,2.48884,1.112024,3.865655,0,1,2,3,0,5.154023,2.685469,7.622577,7.582648,4.134116,11.03118,5.645098,2.915675,8.374521,-0.1764055,-0.4078213,0.05501025,0.7096472,0.447945,0.9713493,0,1,2,3,0,2.167818,1.586001,2.749634,3.048915,2.382779,3.715051,1.96501,1.268023,2.661996,0,1),
                 tolerance = 1e-2)

    mod3 <- mixedmirt(Science, covdat, model=model, draws = 1,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE)
    expect_is(mod3, 'MixedClass')
    expect_equal(extract.mirt(mod3, 'df'), 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod3))))
    expect_equal(cfs, c(-0.2481932,-0.5779923,0.08160588,0.9996557,0.63628,1.363031,4.951115,3.989263,5.912968,2.735017,2.288555,3.18148,-1.326709,-1.681992,-0.971426,-0.2481932,-0.5779923,0.08160588,1.21205,0.8732852,1.550815,3.033993,2.545489,3.522498,1.019672,0.6971102,1.342234,-2.132917,-2.553508,-1.712325,-0.2481932,-0.5779923,0.08160588,2.580794,1.250539,3.91105,5.76086,3.732607,7.789113,2.523845,1.49681,3.55088,-1.994636,-2.807197,-1.182076,-0.2481932,-0.5779923,0.08160588,1.075305,0.7278781,1.422733,3.456143,2.90494,4.007346,1.10973,0.7960739,1.423387,-1.553964,-1.924532,-1.183396,0,1),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- try(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE, SE=FALSE), TRUE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(extract.mirt(rmod1, 'df'), 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(rmod1))))
    expect_equal(cfs, c(1.032039,4.865985,2.63557,-1.47773,1.231335,2.936047,0.8979343,-2.2878,2.234105,5.191735,2.178876,-1.943804,1.059907,3.328831,0.9763014,-1.689545,0,1,0.01040838),
                 tolerance = 1e-4)
    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)


    ## latent regression
    if(FALSE){
        rm(list=ls())
        set.seed(1234)
        n <- 250
        Theta <- matrix(c(rnorm(n, -1, sd = sqrt(1/3)),
                          rnorm(n,0, sd = sqrt(1/3)),
                          rnorm(n, 1, sd = sqrt(1/3))))
        dat <- simdata(matrix(rlnorm(10)), matrix(rnorm(10)), N=n*3, Theta=Theta, itemtype = 'dich')
        covdata <- data.frame(group=rep(c('g1', 'g2', 'g3'), each=n))
        save(dat, covdata, Theta, file = 'tests/tests/testdata/mixed2.rds')
    }
    load('testdata/mixed2.rds')

    mod <- mixedmirt(dat, covdata = covdata, 1, itemtype = '2PL', fixed = ~ 0 + items,
                     lr.fixed = ~ group, verbose=FALSE, SE=TRUE)
    cfs <- coef(mod, printSE=TRUE)
    expect_equal(as.numeric(cfs$lr.betas)[-c(1:2)], c(1.2903165, 0.1097087, 2.9213031, 0.1191159),
                 tolerance=1e-4)
    expect_equal(extract.mirt(mod, 'logLik'), -4314.043, tolerance = 1e-2)
    expect_equal(extract.mirt(mod, 'df'), 1001)

    mod2 <- mixedmirt(dat, covdata = covdata, 1, itemtype = 'Rasch', fixed = ~ 0 + items,
                     lr.fixed = ~ group, verbose=FALSE, draws=1, SE=FALSE)
    so <- summary(mod2, verbose=FALSE)
    expect_equal(as.numeric(c(so$random$Theta, so$lr.out[,1])),
                 c(0.2094428,0,0.7084939,1.712383), tolerance=1e-3)
    expect_equal(extract.mirt(mod2, 'logLik'), -4647.105, tolerance = 1e-4)
    # set.seed(1)
    # bs <- boot.mirt(mod2, R = 3)
    # expect_is(bs, 'boot')

    #uncorrelated random slope
    covdata$theta <- Theta
    mod <- try(mixedmirt(dat, covdata = covdata, 1, fixed = ~ 0 + items, SE=FALSE,
                     random = ~ -1 + theta|group, verbose=FALSE, draws=1), TRUE)
    so <- summary(mod, verbose=FALSE)
    out <- try(as.numeric(diag(so$random$group)), TRUE)
    if(!is(out, 'try-error'))
        expect_equal(as.numeric(diag(so$random$group)), c(0.2859471, 0.6999401), tolerance = 1e-4)

})
