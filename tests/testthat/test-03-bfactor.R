expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('dich data', {
    key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
    data <- key2binary(SAT12, key)
    specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
    mod1 <- bfactor(data, specific, verbose=FALSE, SE = TRUE, SE.type = 'crossprod')
    expect_class(mod1, 'SingleGroupClass')
    expect_equal(extract.mirt(mod1, 'df'), 4294967199)
    cfs <- do.call(c, coef(mod1))
    cfs <- as.numeric(na.omit(cfs[cfs != 0 & cfs != 1]))
    expect_equal(cfs, c(0.7861, 0.5159, 1.0563, 0.4426, 0.0102, 0.875, -1.0729, -1.3116, -0.8343, 1.4909, 1.066, 1.9157, 0.8135, 0.3039, 1.3231, 0.471, 0.1911, 0.7509, 1.1476, 0.8265, 1.4687, -0.1526, -0.5584, 0.2532, -1.1727, -1.4361, -0.9093, 0.5277, 0.296, 0.7595, 0.5813, 0.1573, 1.0053, -0.5576, -0.7686, -0.3465, 0.968, 0.6684, 1.2677, 0.5136, 0.0972, 0.9299, 0.6283, 0.3893, 0.8672, 1.1389, 0.8043, 1.4734, 0.5835, 0.0236, 1.1435, -2.1287, -2.5185, -1.7389, 1.0383, 0.6538, 1.4228, 0.9505, 0.3252, 1.5758, 1.5624, 1.1836, 1.9412, 0.6731, 0.3993, 0.9468, 0.524, 0.0133, 1.0347, -1.57, -1.8547, -1.2853, 0.4473, 0.051, 0.8437, 1.0716, 0.2209, 1.9223, 2.4985, 1.8584, 3.1386, 1.0389, 0.7071, 1.3706, 0.8009, 0.2857, 1.316, -0.4025, -0.6423, -0.1626, 1.5827, 0.1568, 3.0087, 0.8549, -0.9418, 2.6517, 5.4019, 3.4336, 7.3702, 0.1221, -0.0922, 0.3364, 0.2729, -0.0654, 0.6113, -0.3504, -0.532, -0.1689, 1.0958, 0.7507, 1.4409, 0.5884, 0.1265, 1.0502, 0.8847, 0.629, 1.1404, 1.0908, 0.6922, 1.4894, 1.034, 0.4259, 1.642, 1.3605, 1.0081, 1.713, 1.3277, 0.9036, 1.7518, 0.5674, 0.0418, 1.0931, 2.0162, 1.6303, 2.4022, 0.7375, 0.4793, 0.9957, 0.3912, 0.0054, 0.777, -0.3909, -0.5967, -0.185, 1.5223, 0.7735, 2.271, 0.3034, -0.545, 1.1518, 4.1836, 3.3304, 5.0369, 1.7575, 1.3192, 2.1959, 0.1929, -0.3447, 0.7304, -0.8748, -1.1651, -0.5845, 0.8645, 0.5785, 1.1506, 0.0331, -0.3091, 0.3754, 0.2388, 0.0304, 0.4472, 1.5332, 0.9808, 2.0856, 0.4104, -0.1434, 0.9642, 2.6505, 2.1186, 3.1823, 0.5263, 0.0806, 0.9721, 0.6773, -0.1142, 1.4687, 2.6609, 2.1603, 3.1615, 1.6737, 0.8655, 2.4818, -0.0344, -0.7411, 0.6722, 3.6147, 2.7752, 4.4542, 0.6087, 0.3653, 0.8522, 0.4945, 0.0858, 0.9031, -0.8815, -1.1074, -0.6556, 1.2319, 0.8525, 1.6114, 0.2423, -0.1535, 0.6382, 1.2882, 1.0023, 1.5741, 0.7354, 0.476, 0.9947, 0.6426, 0.2139, 1.0712, -0.5995, -0.8249, -0.3742, 1.4898, 1.0991, 1.8806, 0.4859, 0.0468, 0.9251, -0.1725, -0.4279, 0.0828, 1.9058, 1.3075, 2.504, 0.415, -0.169, 0.9991, 2.8064, 2.2378, 3.375, 1.0572, 0.761, 1.3534, 0.149, -0.2358, 0.5338, 0.1725, -0.0414, 0.3863, 1.2276, 0.0476, 2.4075, 2.0402, -1.1072, 5.1875, -1.1821, -2.307, -0.0571, 0.4342, 0.1946, 0.6738, -0.1733, -0.5266, 0.1801, -0.2521, -0.4383, -0.0658, 2.5816, 1.5551, 3.6081, -0.235, -0.9097, 0.4398, 2.9944, 2.0731, 3.9156, 0.1326, -0.1185, 0.3838, 0.0294, -0.3714, 0.4302, -1.6521, -1.8929, -1.4113),
                 tolerance = 1e-2)
    fs <- suppressWarnings(fscores(mod1, verbose = FALSE, full.scores=FALSE, method='MAP'))
    expect_class(fs, 'matrix')
    expect_equal(fs[1:6,'G'], c(-1.0103793, -0.7807428, -1.2088350, -1.4099888, -2.0074429, -1.2226505), tolerance = 1e-2)
    cof <- coef(mod1, verbose = FALSE)
    expect_class(cof, 'list')
    sum <- summary(mod1, verbose = FALSE)
    expect_class(sum, 'list')
    fs <- suppressWarnings(fscores(mod1, method = 'EAPsum', verbose = FALSE, full.scores=FALSE))
    expect_equal(fs[1:3,'G'], c(-3.153030, -2.973976, -2.786622), tolerance = 1e-4)
    expect_equal(fs[1:3,'S1'], c(-1.319310, -1.097647, -0.893891), tolerance = 1e-4)
    fit <- suppressWarnings(M2(mod1))
    expect_equal(fit$M2, 553.6781, tolerance = 1e-2)
    expect_equal(fit$df, 432, tolerance = 1e-2)
    pfit1 <- suppressWarnings(personfit(mod1))
    expect_class(pfit1, 'data.frame')
    ifit <- suppressWarnings(itemfit(mod1))
    expect_class(ifit, 'data.frame')
    mod2 <- bfactor(data, specific, pars=mod2values(mod1), TOL=NaN,
                    verbose=FALSE, SE = TRUE, SE.type = 'Oakes')
    expect_equal(extract.mirt(mod2, 'condnum'), 619.3619, tolerance = 1e-4)

    #nestlogit
    scoredSAT12 <- data
    scoredSAT12[,1:5] <- as.matrix(SAT12[,1:5])
    scoredSAT12[scoredSAT12 == 8] <- NA
    nestmod <- mirt(scoredSAT12, 1, c(rep('2PLNRM',5),rep('2PL', 27)), key=key, verbose=FALSE)
    cfs <- as.numeric(do.call(c, coef(nestmod)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(0.8075,-1.0421,-0.5669,-0.5737,-3.0614,0.2125,0.0698,-5.3709,1.5103,0.444,-0.9065,-0.3884,-0.2221,-2.9412,-1.3382,-0.6365,1.0693,-1.1192,0.1294,-0.1866,-0.332,0.1263,0.4042,-0.624,0.576,-0.5181,0.0332,-0.0055,0.1145,-0.0994,0.0384,-0.2303,0.9778,0.6075,0.4063,0.0481,-0.0747,0.6153,0.0278,-0.7043,1.1276,-2.0376,1.0278,1.3932,0.7003,-1.5112,0.5251,2.1404,1.003,-0.3598,1.7028,5.2068,0.1643,-0.3455,1.0956,0.8494,1.0374,1.1747,1.3068,1.9325,0.7252,-0.3815,1.4946,4.1057,1.7048,-0.8517,0.8311,0.2367,1.5509,2.6202,0.6068,2.5177,1.5574,3.4923,0.6458,-0.8513,1.2023,1.2685,0.7725,-0.5667,1.5348,-0.1701,1.8923,2.7513,1.0725,0.1742,0.832,-0.7494,0.3767,-0.2478,2.4171,2.8476,0.1256,-1.6512),
                 tolerance = 1e-2)
    expect_equal(extract.mirt(nestmod, 'logLik'), -11626.06, tolerance = 1e-4)

    #simulate data
    if(FALSE){
        rm(list=ls())

        set.seed(1234)
        a <- matrix(c(
            1,0.5,NA,
            1,0.5,NA,
            1,0.5,NA,
            1,0.5,NA,
            1,0.5,NA,
            1,0.5,NA,
            1,0.5,NA,
            1,NA,0.5,
            1,NA,0.5,
            1,NA,0.5,
            1,NA,0.5,
            1,NA,0.5,
            1,NA,0.5,
            1,NA,0.5),ncol=3,byrow=TRUE)

        d <- matrix(c(
            -1.0,NA,NA,
            -1.5,NA,NA,
            1.5,NA,NA,
            0.0,NA,NA,
            0.0,-1.0,1.5,
            0.0,2.0,-0.5,
            3.0,2.0,-0.5,
            3.0,2.0,-0.5,
            2.5,1.0,-1,
            2.0,0.0,NA,
            -1.0,NA,NA,
            -1.5,NA,NA,
            1.5,NA,NA,
            0.0,NA,NA),ncol=3,byrow=TRUE)

        nominal <- matrix(NA, nrow(d), ncol(d))
        nominal[5, ] <- c(0,1.2,2)
        sigma <- diag(3)
        set.seed(1234)
        items <- itemtype <- c(rep('dich', 4), 'nominal', 'gpcm', rep('graded',4),rep('dich', 4))
        dataset <- simdata(a,d,3000,itemtype, sigma=sigma, nominal=nominal)

        save(dataset, items, file = 'tests/tests/testdata/bfactor1.rds')
    }
    load('testdata/bfactor1.rds')
    specific <- c(rep(1,7),rep(2,7))
    items[items == 'dich'] <- '2PL'
    simmod <- bfactor(dataset, specific, itemtype = items, verbose=FALSE, TOL=3e-3)
    expect_class(simmod, 'SingleGroupClass')
    expect_equal(extract.mirt(simmod, 'df'), 442315)
    cfs <- as.numeric(do.call(c, coef(simmod)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(0.9336,0.6975,-0.9818,0.9105,0.4485,-1.4454,1.0832,0.5553,1.495,0.9833,0.4433,0.0846,1.062,0.5301,1.0833,2,-0.9409,1.4898,1.0732,0.4348,2,1.9718,-0.5391,0.9738,0.4775,2.8693,1.8961,-0.5221,1.0478,0.4884,3.0804,2.0606,-0.5033,0.9921,0.5407,2.5673,1.1235,-0.9224,0.955,0.5072,1.9868,-0.0063,1.115,0.4878,-0.9547,1.0441,0.655,-1.4294,1.008,0.4946,1.4922,1.0458,0.5287,0.0506),
                 tolerance = 1e-2)
    specific[1] <- NA
    simmod2 <- bfactor(dataset, specific, itemtype = items, verbose=FALSE, TOL=3e-3)
    expect_class(simmod2, 'SingleGroupClass')
    expect_equal(extract.mirt(simmod2, 'df'), 442316)
    cfs <- as.numeric(do.call(c, coef(simmod2)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.0905,-0.9601,1.0223,-0.0369,-1.4489,1.2195,-0.0022,1.4959,1.092,-0.014,0.0848,1.2755,0.4096,1.0859,2,-0.8057,1.6183,1.2266,-0.2,2,2.0288,-0.5539,1.0768,0.333,2.9025,1.9204,-0.5287,0.9394,0.6695,3.0785,2.0591,-0.5033,0.8864,0.7024,2.5682,1.1235,-0.9231,0.8554,0.6607,1.9864,-0.0067,0.9996,0.6903,-0.9545,0.9342,0.8017,-1.4291,0.9083,0.6546,1.4906,0.932,0.7122,0.0501),
                 tolerance = 1e-2)
    fs <- fscores(simmod, verbose = FALSE, full.scores=FALSE)
    expect_true(mirt:::closeEnough(fs[1:6,'G'] - c(-2.6461, -2.0541, -2.3838, -2.4061, -1.7492, -1.8682), -1e-2, 1e-2))
    expect_class(fs, 'matrix')

    res <- residuals(simmod, verbose = FALSE, QMC=TRUE)
    expect_class(res, 'matrix')
    expect_equal(res[2,1], 0.08305725, tolerance = 1e-2)
    sum <- summary(simmod, verbose = FALSE)
    expect_class(sum, 'list')

    #two-tier
    if(FALSE){
        rm(list=ls())
        set.seed(1234)
        a <- matrix(c(
            0,1,0.5,NA,NA,
            0,1,0.5,NA,NA,
            0,1,0.5,NA,NA,
            0,1,0.5,NA,NA,
            0,1,0.5,NA,NA,
            0,1,NA,0.5,NA,
            0,1,NA,0.5,NA,
            0,1,NA,0.5,NA,
            1,0,NA,0.5,NA,
            1,0,NA,0.5,NA,
            1,0,NA,0.5,NA,
            1,0,NA,NA,0.5,
            1,0,NA,NA,0.5,
            1,0,NA,NA,0.5,
            1,0,NA,NA,0.5,
            1,0,NA,NA,0.5),ncol=5,byrow=TRUE)

        d <- matrix(rnorm(16))
        items <- rep('dich', 16)

        sigma <- diag(5)
        sigma[1,2] <- sigma[2,1] <- .7
        dataset <- simdata(a,d,2000,itemtype=items,sigma=sigma)

        save(dataset, items, file = 'tests/tests/testdata/bfactor2.rds')
    }
    load('testdata/bfactor2.rds')
    specific <- c(rep(1,5),rep(2,6),rep(3,5))
    model <- mirt.model('
        G1 = 1-8
        G2 = 9-16
        COV = G1*G2')

    simmod <- bfactor(dataset, specific, model, quadpts = 11, TOL = 1e-2, verbose=FALSE)
    expect_class(simmod, 'SingleGroupClass')
    expect_equal(extract.mirt(simmod, 'df'), 65486)
    cfs <- as.numeric(do.call(c, coef(simmod)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.014,0.3654,-1.1222,1.0028,0.3867,0.3248,0.9736,0.4333,1.0799,0.9276,0.6708,-2.2706,1.002,0.7353,0.5082,0.9243,0.5424,0.4103,1.1482,0.5269,-0.659,0.9209,0.5934,-0.6195,0.8717,0.4396,-0.5551,1.0235,0.6562,-0.9677,1.1029,0.4787,-0.5439,1.0948,0.6413,-0.987,1.0947,0.4871,-0.706,1.0222,0.4093,-0.0148,0.9406,0.4599,0.8405,0.9952,0.6163,-0.074,0.6361),
                 tolerance = 1e-2)
})
