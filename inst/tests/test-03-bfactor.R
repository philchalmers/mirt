context('bfactor')

test_that('dich data', {
    data <- key2binary(SAT12,
                       key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
    specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
    mod1 <- bfactor(data, specific, verbose=FALSE)    
    expect_is(mod1, 'ConfirmatoryClass')
    expect_equal(mod1@df, 501)
    cfs <- as.numeric(do.call(c, coef(mod1)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(0.785,  0.438, -1.072,  1.490,  0.816,  0.471,  1.147, -0.153, -1.173,  
                        0.526,  0.581, -0.557,  0.968,  0.515,  0.628,  1.136,  0.580,
                       -2.126,  1.044,  0.925,  1.556,  0.673,  0.520, -1.569,  0.455,  1.110,
                        2.526,  1.048,  0.821, -0.405,  1.592,  0.818,  5.388,  0.121,
                        0.276, -0.351,  1.100,  0.577,  0.884,  1.087,  1.036,  1.360,  1.334,
                        0.532,  2.011,  0.736,  0.389, -0.391,  1.530,  0.268,  4.184,
                        1.758,  0.190, -0.875,  0.862,  0.038,  0.239,  1.530,  0.418,  2.650,
                        0.534,  0.655,  2.653,  1.682, -0.070,  3.624,  0.607,  0.499,
                       -0.882,  1.238,  0.216,  1.289,  0.732,  0.648, -0.600,  1.487,  0.492,
                        -0.172,  1.910,  0.401,  2.807,  1.056,  0.153,  0.172,  1.248,
                        2.108, -1.205,  0.433, -0.169, -0.252,  2.600, -0.275,  3.013,  
                        0.133,  0.029, -1.652), tollerance = 1e-2)
    fs <- fscores(mod1, verbose = FALSE)
    expect_is(fs, 'matrix')        
    expect_true(mirt:::closeEnough(fs[1:6,'F1'] - c(-0.72059353, -0.07439544, -1.91235316,
                                                    -1.99421796, -2.00030284, -1.92556420), -1e-4, 1e-4))
    cof <- coef(mod1, verbose = FALSE)
    expect_is(cof, 'list')
    sum <- summary(mod1, verbose = FALSE)
    expect_is(sum, 'list')    
    pfit1 <- personfit(mod1)
    expect_is(pfit1, 'data.frame')    
    ifit <- itemfit(mod1)
    expect_is(ifit, 'data.frame')

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
    expect_equal(simmod@df, 2442)
    cfs <- as.numeric(do.call(c, coef(simmod)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.104,  0.080, -1.001,  1.163, -0.500, -1.520,  1.113,  0.636,  1.635, 
                        0.983,  0.088,  0.069,  1.146,  0.255,  1.104,  2.000, -0.936,
                        1.598,  1.154, -0.142,  2.082, -0.403,  1.125,  0.094,  3.094,  2.027, 
                        -0.456,  0.998,  0.614,  3.081,  2.063, -0.454,  0.843,  0.636,
                        2.451,  0.972, -0.983,  1.018,  0.565,  2.025, -0.031,  0.841,  
                        0.758, -1.014,  0.938,  0.524, -1.363,  0.881,  0.493,  1.546,  0.951,
                        0.765,  0.028), tollerance = 1e-2)   
    specific[1] <- NA
    simmod2 <- suppressMessages(bfactor(dataset, specific, itemtype = items, verbose=FALSE))
    expect_is(simmod2, 'ConfirmatoryClass')              
    expect_equal(simmod2@df, 2443)
    cfs <- as.numeric(do.call(c, coef(simmod2)))
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.107, -1.002,  1.153, -0.502, -1.516,  1.124,  0.641,  1.641,  0.985,  
                        0.071,  0.069,  1.150,  0.244,  1.103,  2.000, -0.935,  1.600,
                        1.150, -0.141,  2.079, -0.402,  1.126,  0.070,  3.094,  2.027, -0.456, 
                        0.997,  0.616,  3.081,  2.063, -0.454,  0.842,  0.637,  2.451,
                        0.972, -0.983,  1.017,  0.567,  2.025, -0.031,  0.840,  0.759, -1.013, 
                        0.938,  0.524, -1.363,  0.879,  0.495,  1.546,  0.950,  0.766,  0.028),
                 tollerance = 1e-2)   
    fs <- fscores(simmod, verbose = FALSE)
    expect_true(mirt:::closeEnough(fs[1:6,'F1'] - c(-2.713717, -2.440282, -2.177029,
                                                    -2.265682, -2.249449, -2.416284), -1e-4, 1e-4))
    expect_is(fs, 'matrix')
    
    res <- residuals(simmod, verbose = FALSE)
    expect_is(res, 'matrix')
    fit <- fitted(simmod)
    expect_is(fit, 'matrix')  
    sum <- summary(simmod, verbose = FALSE)
    expect_is(sum, 'list')
})
