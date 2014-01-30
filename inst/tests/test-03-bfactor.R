context('bfactor')

test_that('dich data', {
    key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
    data <- key2binary(SAT12, key)
    specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
    mod1 <- bfactor(data, specific, verbose=FALSE)
    expect_is(mod1, 'ConfirmatoryClass')
    expect_equal(mod1@df, 4294967199)
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(0.7851, 0.4358, -1.0717, 1.4896, 0.8156, 0.4711, 1.1475, -0.1538, -1.1727, 0.5263, 0.5811, -0.5573, 0.9675, 0.5152, 0.6284, 1.1361, 0.5778, -2.1255, 1.0436, 0.9261, 1.5563, 0.6731, 0.5188, -1.5687, 0.4541, 1.108, 2.5245, 1.0477, 0.8196, -0.4047, 1.5918, 0.8201, 5.3882, 0.1206, 0.2755, -0.3506, 1.1001, 0.5777, 0.8844, 1.0876, 1.0359, 1.3602, 1.3332, 0.534, 2.0112, 0.7362, 0.387, -0.3906, 1.5293, 0.2701, 4.1835, 1.758, 0.1904, -0.8748, 0.8626, 0.0379, 0.2387, 1.5305, 0.4172, 2.6498, 0.5335, 0.6558, 2.6531, 1.6817, -0.0683, 3.6239, 0.6068, 0.4983, -0.8817, 1.2375, 0.2173, 1.2886, 0.7326, 0.6474, -0.5997, 1.4874, 0.4916, -0.1725, 1.9096, 0.4021, 2.807, 1.0559, 0.1523, 0.1724, 1.2577, 2.1328, -1.2139, 0.4331, -0.1688, -0.252, 2.599, -0.273, 3.0115, 0.1328, 0.028, -1.6521),
                 tolerance = 1e-2)
    fs <- fscores(mod1, verbose = FALSE)
    expect_is(fs, 'matrix')
    expect_true(mirt:::closeEnough(fs[1:6,'F1'] - c(-0.72059353, -0.07439544, -1.91235316,
                                                    -1.99421796, -2.00030284, -1.92556420), -1e-2, 1e-2))
    cof <- coef(mod1, verbose = FALSE)
    expect_is(cof, 'list')
    sum <- summary(mod1, verbose = FALSE)
    expect_is(sum, 'list')
    pfit1 <- personfit(mod1)
    expect_is(pfit1, 'data.frame')
    ifit <- itemfit(mod1)
    expect_is(ifit, 'data.frame')
    
    #nestlogit
    scoredSAT12 <- data
    scoredSAT12[,1:5] <- as.matrix(SAT12[,1:5])
    nestmod <- mirt(scoredSAT12, 1, c(rep('2PLNRM',5),rep('2PL', 27)), key=key, verbose=FALSE)
    cfs <- as.numeric(do.call(c, coef(nestmod, digits=4)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(0.7928, -1.0428, -0.5715, -0.5787, -3.1137, 1.4113, 0.2128, 0.07, -5.4259, -5.7802, 1.5217, 0.4431, -0.919, -0.4077, -0.2366, -1.6949, -2.9472, -1.3499, -0.6444, -6.4068, 1.063, -1.1381, 0.1255, -0.1845, -0.3254, -0.2015, 0.1254, 0.4052, -0.6207, -2.5713, 0.5877, -0.5306, 0.0333, -0.0066, 0.1137, -1.1324, -0.0994, 0.0382, -0.2306, -3.6532, 0.983, 0.605, 0.4158, 0.0505, -0.0696, -1.3233, 0.6194, 0.0295, -0.7007, -5.3263, 1.1315, -2.04, 1.0201, 1.3905, 0.7031, -1.5122, 0.5175, 2.1378, 0.9976, -0.3594, 1.7057, 5.2077, 0.1694, -0.3456, 1.0861, 0.847, 1.0424, 1.177, 1.327, 1.9449, 0.7268, -0.3816, 1.4883, 4.0989, 1.6847, -0.8479, 0.8226, 0.2359, 1.5806, 2.6433, 0.6012, 2.5153, 1.5424, 3.4789, 0.649, -0.8519, 1.2108, 1.2726, 0.7751, -0.567, 1.5423, -0.17, 1.8603, 2.7282, 1.076, 0.1748, 0.8315, -0.7493, 0.3816, -0.2479, 2.3798, 2.8218, 0.1269, -1.6513),
                 tolerance = 1e-2)
    expect_equal(nestmod@logLik, -11715.17, tolerance = 1e-4)

    #simulate data
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

    specific <- c(rep(1,7),rep(2,7))
    items[items == 'dich'] <- '2PL'
    simmod <- suppressMessages(bfactor(dataset, specific, itemtype = items, verbose=FALSE))
    expect_is(simmod, 'ConfirmatoryClass')
    expect_equal(simmod@df, 442315)
    cfs <- as.numeric(do.call(c, coef(simmod, digits=4)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.1041, 0.08, -1.0015, 1.1639, -0.5016, -1.5207, 1.1132, 0.6341, 1.6349, 0.9834, 0.0885, 0.0694, 1.1464, 0.2551, 1.1043, 2, -0.9361, 1.5983, 1.1543, -0.1421, 2, 2.0817, -0.4029, 1.1246, 0.0932, 3.0943, 2.0269, -0.4557, 0.9981, 0.6139, 3.0811, 2.0627, -0.4541, 0.8432, 0.6361, 2.4515, 0.9724, -0.9835, 1.018, 0.5648, 2.0248, -0.0309, 0.8408, 0.7581, -1.0135, 0.9375, 0.5238, -1.3633, 0.8808, 0.4932, 1.5456, 0.951, 0.7649, 0.0282),
                 tolerance = 1e-2)
    specific[1] <- NA
    simmod2 <- suppressMessages(bfactor(dataset, specific, itemtype = items, verbose=FALSE))
    expect_is(simmod2, 'ConfirmatoryClass')
    expect_equal(simmod2@df, 442316)
    cfs <- as.numeric(do.call(c, coef(simmod2, digits=4)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.1073, -1.0015, 1.153, -0.5001, -1.5158, 1.1239, 0.6421, 1.6415, 0.9848, 0.0714, 0.0694, 1.15, 0.2443, 1.1028, 2, -0.9349, 1.5996, 1.1502, -0.1402, 2, 2.0788, -0.4023, 1.1264, 0.0702, 3.0944, 2.027, -0.4557, 0.997, 0.6156, 3.0811, 2.0627, -0.454, 0.8424, 0.6371, 2.4515, 0.9724, -0.9835, 1.0168, 0.5666, 2.0248, -0.0308, 0.8398, 0.7591, -1.0135, 0.938, 0.5235, -1.3634, 0.8795, 0.4952, 1.5456, 0.9501, 0.7663, 0.0282),
                 tolerance = 1e-2)
    fs <- fscores(simmod, verbose = FALSE)
    expect_true(mirt:::closeEnough(fs[1:6,'F1'] - c(-2.713717, -2.440282, -2.177029,
                                                    -2.265682, -2.249449, -2.416284), -1e-2, 1e-2))
    expect_is(fs, 'matrix')

    res <- residuals(simmod, verbose = FALSE)
    expect_is(res, 'matrix')
    sum <- summary(simmod, verbose = FALSE)
    expect_is(sum, 'list')
})
