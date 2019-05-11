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
    MGmodel1 <- 'F1 = 1-15'
    models <- mirt.model(MGmodel1, quiet = TRUE)

    mod_Rasch <- multipleGroup(dat, models, itemtype = 'Rasch', SE=TRUE, SE.type = 'crossprod',
                               group = group, verbose = FALSE, method = 'EM')
    cfs <- as.numeric(na.omit(do.call(rbind, coef(mod_Rasch, printSE=TRUE, as.data.frame=TRUE))))
    expect_equal(cfs, c(0.5116698,-0.6522073,-0.1905068,0.8651635,0.1235676,0.7636483,0.919669,-0.3367473,-1.0904,-1.166159,1.242589,-0.2147728,0.4069149,0.4466409,-0.06474946,0.9783792,0.5880477,-0.4385498,-0.2643517,1.060335,0.2625935,0.7648064,0.9712697,-0.3590464,-0.9936327,-1.248645,1.270672,-0.09205005,0.4425358,0.4746364,-0.03480865,1.496147,0.07893867,0.0792895,0.07720629,0.08119347,0.07785767,0.0798737,0.08224574,0.07738347,0.08298824,0.08428115,0.08508337,0.07769511,0.07803545,0.07844748,0.07726093,0.0671328,0.08577857,0.08527442,0.08357395,0.08807036,0.08421204,0.0865046,0.08779314,0.08364731,0.08749941,0.0896606,0.09115384,0.08363036,0.08487376,0.08502025,0.0834237,0.1013833),
                 tolerance = 1e-4)
    expect_equal(logLik(mod_Rasch), -18007.78, tolerance = 1e-4)
    EAP <- fscores(mod_Rasch, full.scores=TRUE)
    expect_equal(cor(EAP, rowSums(dat))[1], .99, tolerance = 1e-2)
    pf <- personfit(mod_Rasch, Theta=EAP)
    pffit <- c(as.numeric(as.matrix(head(pf))), as.numeric(as.matrix(tail(pf))))
    expect_equal(pffit, c(0.950921,1.142593,0.5986253,0.4390725,0.8847643,0.6006911,-0.1580319,0.4765587,-1.366033,-1.688994,-0.1716739,-1.356961,0.9657947,0.9039326,0.6798182,0.5300499,0.820473,0.6796108,-0.12207,-0.2349138,-1.365052,-1.825028,-0.5448986,-1.366125,0.1942611,0.1375458,1.282308,1.543602,0.5053975,1.280061,0.9206344,1.182217,1.033595,1.287318,0.4910474,0.365361,-0.2859851,0.6103204,0.2228756,0.8827224,-1.29217,-0.998842,0.9182596,1.20314,1.062435,1.285586,0.6090764,0.5207685,-0.37187,0.8019668,0.4074629,1.176899,-1.291869,-1.025544,0.4027076,-0.6934658,-0.2963764,-1.118536,1.201121,0.9925202),
                 tolerance = 1e-3)
    # mod_QMCEM <- multipleGroup(dat, models, group=group, method = 'QMCEM', verbose=FALSE,
    #                            optimizer='NR')
    # expect_equal(extract.mirt(mod_QMCEM, 'logLik'), -17849.64, tolerance=1e-2)
    mod_configural <- multipleGroup(dat, models, SE=TRUE, SE.type = 'crossprod', optimizer='NR',
                                    group = group, verbose = FALSE, method = 'EM')
    expect_is(mod_configural, 'MultipleGroupClass')
    expect_equal(logLik(mod_configural), -17903.76, tolerance = 1e-4)
    cfs <- as.numeric(do.call(c, coef(mod_configural)[[1L]]))
    cfs <- as.numeric(na.omit(cfs[cfs != 0 & cfs != 1]))
    expect_equal(cfs, c(1.0706,0.8462,1.295,0.524,0.3627,0.6852,1.217,0.9668,1.4672,-0.6995,-0.8745,-0.5244,0.9488,0.7445,1.1532,-0.1877,-0.3385,-0.037,0.8975,0.6948,1.1002,0.8423,0.681,1.0036,1.097,0.8781,1.316,0.1274,-0.0294,0.2841,0.5648,0.3916,0.7381,0.6824,0.5383,0.8266,1.2742,1.0146,1.5337,1.0016,0.8129,1.1904,0.9249,0.7242,1.1255,-0.3298,-0.4801,-0.1794,0.8903,0.6815,1.0992,-1.0591,-1.2284,-0.8897,0.723,0.5272,0.9188,-1.0828,-1.2463,-0.9192,0.8303,0.6252,1.0354,1.1879,1.0155,1.3604,1.4757,1.1986,1.7528,-0.252,-0.429,-0.0749,1.2905,1.0339,1.5471,0.4451,0.2742,0.616,1.0348,0.8239,1.2457,0.4526,0.2955,0.6096,0.8637,0.665,1.0624,-0.0617,-0.2093,0.0859),
                 tolerance = 1e-2)
    expect_equal(extract.mirt(mod_configural, 'df'), 32707)
    mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), verbose = FALSE,
                                method = 'EM', optimizer = 'NR')
    expect_is(mod_metric, 'MultipleGroupClass')
    expect_equal(extract.mirt(mod_metric, 'df'), 32722)
    mod_scalar2 <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_var','free_means'))
    expect_equal(logLik(mod_scalar2), -17915.33, tolerance = 1e-4)
    cfs <- as.numeric(do.call(c, coef(mod_scalar2)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.104323,0.5385656,1.20104,-0.6295554,1.061246,-0.2649867,0.8824239,0.9004537,1.081875,0.1640334,0.4445889,0.636301,1.180488,0.9758879,0.977318,-0.3768236,0.8786338,-1.034887,0.6847666,-1.117687,0.948221,1.213217,1.389262,-0.2233858,1.261964,0.4288545,1.153485,0.4531765,0.7745967,-0.06987377),
                 tolerance = 1e-2)
    expect_is(mod_scalar2, 'MultipleGroupClass')
    expect_equal(extract.mirt(mod_scalar2, 'df'), 32735)
    newmodel <- mirt.model('F = 1-15
                            CONSTRAINB = (1-15, a1), (1,2,3-15,d)')
    mod_scalar1 <- multipleGroup(dat, newmodel, group = group, verbose = FALSE, invariance='free_var')
    expect_is(mod_scalar1, 'MultipleGroupClass')
    mod_EH <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                            dentype="empiricalhist", optimizer = 'NR')
    expect_is(mod_EH, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_EH)[[1L]]))
    expect_equal(cfs, c(1.019061,0.5052184,0,1,1.209095,-0.7255769,0,1,0.9008162,-0.2009865,0,1,0.8723441,0.8312948,0,1,1.076235,0.1099079,0,1,0.5510721,0.6746446,0,1,1.219269,0.9782172,0,1,0.9121723,-0.3459642,0,1,0.8412548,-1.067533,0,1,0.6879252,-1.096485,0,1,0.7896133,1.173946,0,1,1.419879,-0.2785068,0,1,1.247844,0.4227584,0,1,0.9913394,0.4356204,0,1,0.8361253,-0.074975,0,1,0,1),
                 tolerance = 1e-2)
    set.seed(1)
    mod_mixture <- suppressWarnings(multipleGroup(dat, 1, itemtype = 'Rasch', GenRandomPars = TRUE,
                                 verbose = FALSE, dentype = 'mixture-2', SE=TRUE))
    expect_equal(extract.mirt(mod_mixture, 'condnum'), 112.2736, tolerance=1e-4)
    so <- summary(mod_mixture, verbose=FALSE)
    expect_equal(so[[1]]$class_proportion, .512121, tolerance=1e-4)

    dat[1,1] <- dat[2,2] <- NA
    mod_missing <- multipleGroup(dat, models, group = group, verbose = FALSE, method = 'EM',
                                 invariance=c('slopes', 'intercepts', 'free_var'))
    expect_is(mod_missing, 'MultipleGroupClass')
    expect_equal(extract.mirt(mod_missing, 'df'), 32736)
    out1 <- M2(mod_missing, na.rm=TRUE)
    out2 <- itemfit(mod_missing, na.rm=TRUE)
    out3 <- fscores(mod_missing, na.rm=TRUE, method = 'EAPsum', full.scores=FALSE, verbose = FALSE)
    expect_equal(out1$M2, 166.4414, tolerance=1e-4)
    expect_equal(out2$D1$S_X2[1], 7.633155, tolerance=1e-4)
    expect_equal(out3$D1$expected[1], 5.258545, tolerance=1e-2)

    fs1 <- fscores(mod_metric, verbose = FALSE, full.scores=FALSE)
    expect_true(mirt:::closeEnough(fs1[[1]][1:6, 'F1'] - c(-2.083008, -1.653961, -1.405526, -1.573013, -1.723711, -1.324450), -1e-2, 1e-2))
    fs2 <- fscores(mod_metric, full.scores = TRUE, full.scores.SE=TRUE, method = 'ML')
    expect_equal(as.numeric(head(fs2)), c(0.3827467,1.358339,0.83162,1.289526,1.286425,1.067579,0.4827119,0.6292271,0.531026,0.6137522,0.613074,0.5693166),
                tolerance = 1e-2)
    fs3 <- fscores(mod_missing, verbose = FALSE, full.scores=FALSE)
    fs4 <- fscores(mod_missing, full.scores = TRUE)
    fs5 <- fscores(mod_metric, full.scores = TRUE, scores.only=TRUE)
    expect_is(fs1, 'list')
    expect_is(fs2, 'matrix')
    expect_is(fs3, 'list')
    expect_is(fs4, 'matrix')

    fit1 <- M2(mod_metric)
    expect_is(fit1, 'data.frame')
    expect_true(mirt:::closeEnough(fit1$M2 - 184.5311, -1e-2, 1e-2))
    expect_equal(fit1$SRMSR.D1, 0.04055562, tolerance = 1e-4)
    expect_equal(fit1$TLI, 1.001405, tolerance = 1e-4)
    expect_true(mirt:::closeEnough(fit1$df - 195, -1e-4, 1e-4))
    fit2 <- itemfit(mod_metric, c('S_X2', 'Zh'))
    expect_is(fit2, 'list')
    expect_equal(as.numeric(fit2[[1]][1L,]), c(1.000000, 2.733099, 7.851266, 11.000000, 0, 0.726562),
                 tolerance = 1e-4)
    fit3 <- M2(mod_scalar2)
    expect_true(mirt:::closeEnough(fit3$M2 - 165.1392, -1e-4, 1e-4))
    expect_equal(fit3$SRMSR.D1, 0.02754769, tolerance = 1e-4)
    expect_equal(fit3$TLI, 1.005392, tolerance = 1e-4)
    expect_true(mirt:::closeEnough(fit3$df - 208, -1e-4, 1e-4))

    # missing by design
    dat[group == 'D1',1:2] <- NA
    dat[group == 'D2',14:15] <- NA
    mod <- multipleGroup(dat, 1, group, invariance = c('slopes', 'intercepts', 'free_means',
                                                       'free_var'), verbose=FALSE)
    cfs <- coef(mod, simplify=TRUE)
    expect_equal(as.vector(cfs$D1$items[1:3,1:2]), c(1.139968,1.192491,1.093887,0.5647073,-0.5384472,-0.2615739),
                 tolerance=1e-4)
    expect_equal(as.vector(fscores(mod)[1:3,]), c(0.7325525, 0.9279133, 0.5071795), tolerance=1e-4)
    expect_is(plot(mod, type = 'trace'), 'trellis')
    ifit <- itemfit(mod, 'X2')
    expect_equal(as.vector(ifit$D1$p.X2[1:4]), c(NaN,NaN,0.0002541653, 0.0009022802), tolerance=1e-4)

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
    out <- empirical_ES(x)
    expect_equal(as.numeric(out[1,]), c(-0.01984901,0.1397911,-0.02010119,0.1390092,-0.03393749,-2.687888,-0.4657692,3.53414,3.553989), tolerance=1e-4)
    out2 <- empirical_ES(x, DIF = FALSE)
    expect_equal(as.numeric(out2$Value), c(-0.8379638,1.252605,0.9030441,-0.05802242,-0.8532262,1.260376,0.9277066,-1.727569,-1.356197), tolerance=1e-4)

    dat[1,1] <- dat[2,2] <- NA
    x2 <- multipleGroup(dat, 1, group=group, method='EM', verbose = FALSE)
    expect_is(x2, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(x2)[[1L]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_true(mirt:::closeEnough(cfs - c(0.54559,2.87963,2.09834,1.04898,0.08764,-1.01597,-1.88808,-2.80774,0.5721,3.66458,2.86098,2.07309,0.97108,-0.04591,-1.13943,-2.04077,2.19508,3.985,2.93076,1.92402,0.93888,-3e-04,-0.92533,-2.15966,2.75051,5.52755,4.54548,3.22867,2.2928,1.26823,0.37394,-0.49598,0.40521,2.24053,1.23686,0.3006,-0.65966,-1.68032,-2.71625,-3.65993,5.18399,3.21456,2.30067,1.19408,0.1821,-0.8039,-1.83942,-2.96744,2.50677,2.45572,1.40808,0.48888,-0.67763,-1.67403,-2.61442,-3.75363,1.90942,4.6516,3.56213,2.57738,1.45019,0.65221,-0.43677,-1.41043,2.0672,3.3566,2.50567,1.64636,0.6884,-0.28202,-1.35771,-2.31969,2.54578,2.31928,1.26861,0.29768,-0.87537,-1.86427,-2.76917,-3.71823), -1e-2, 1e-2))

    # three factor
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

    mod_metric <- multipleGroup(dat, MGmodelg2, group = group, invariance=c('slopes'), method = 'MHRM',
                                                 verbose = FALSE, draws = 10)
    expect_is(mod_metric, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_metric)[[1]]))[1:20]
    expect_equal(cfs, c(1.223844,0,0,0.5986703,0,1,1.478027,0,0,-0.542961,0,1,1.090491,0,0,-0.2585476,0,1,0.8359449,0),
                 tolerance = 1e-2)
    mod_configural <- multipleGroup(dat, MGmodelg1, group = group, verbose = FALSE, method = 'EM', SE=TRUE,
                                    optimizer = 'NR')
    expect_is(mod_configural, 'MultipleGroupClass')
    cfs <- as.numeric(do.call(c, coef(mod_configural)[[1]]))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(as.numeric(na.exclude(cfs)), c(1.175122,0.8723134,1.47793,0.5905031,0.4203281,0.7606781,1.432219,1.051584,1.812853,-0.5385373,-0.7236171,-0.3534576,1.054938,0.7785073,1.331368,-0.2574002,-0.4120382,-0.1027621,0.9510795,0.6902863,1.211873,0.8776461,0.7088456,1.046447,1.147247,0.8538526,1.440642,0.2582053,0.09907628,0.4173343,0.6090205,0.3766148,0.8414262,0.5171881,0.3763434,0.6580327,1.590855,0.9590698,2.22264,1.129202,0.8384696,1.419934,1.058027,0.7076615,1.408393,-0.5338947,-0.6989874,-0.368802,0.8613521,0.5565073,1.166197,-1.216261,-1.405812,-1.026709,0.5272607,0.2895619,0.7649595,-1.023556,-1.178142,-0.8689692,1.00044,0.707042,1.293839,1.302659,1.102789,1.502529,1.459818,1.046385,1.873251,-0.5183716,-0.7054787,-0.3312645,1.027813,0.7488713,1.306755,0.4739741,0.3155796,0.6323687,1.425603,1.026227,1.824979,0.4233488,0.2429515,0.6037461,0.6053667,0.397785,0.8129484,-0.1520328,-0.2869917,-0.01707401),
                 tolerance = 1e-2)

    fs1 <- fscores(mod_metric, verbose = FALSE, full.scores=FALSE)
    expect_is(fs1, 'list')
    expect_true(mirt:::closeEnough(fs1[[1L]][1:6, 'F3'] - c(-0.3666644,-0.4531451,-0.194526,0.07536214,0.3359495,-0.4426196), -1e-3, 1e-3))
})
