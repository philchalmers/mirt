context('mirt.model')

test_that('syntax', {
    data <- expand.table(LSAT7)
    group <- rep(c('male', 'female'), length.out=1000)
    group2 <- c('male', rep(c('male', 'female', 'other'), each = 333))
    model0 <- 'F = 1-5'
    model1 <- mirt.model('F = 1-5')
    model2 <- mirt.model('F = 1-5
                   CONSTRAIN = (1,2,3-5,a1)')
    model3 <- mirt.model('F = 1-5
                   CONSTRAIN = (2,3-5,a1)
                   PRIOR = (2,3-5, a1, lnorm, .2, .2), (1-2, d, norm, 0, 2)')
    model4 <- mirt.model('F = 1-5
                   CONSTRAIN = (1-2, d)
                   CONSTRAINB = (2-4,5,a1), (1, a1)')
    model5 <- mirt.model('F = 1-5
                   CONSTRAIN [male] = (1-5, d)
                   CONSTRAINB = (1-4,5,a1)
                   PRIOR [female] = (1-5, d, norm, 0, 2)')
    model6 <- 'F1 = 1-2
               F2 = 3-5
               CONSTRAIN = (3-5, a2), (1-2, a1)
               COV = F1*F2'
    model7 <- mirt.model('F1 = 1-2
                         F2 = 3-5
                         START = (2, a2, 1.5), (4,a1,-1)')
    model8 <- mirt.model('F1 = 1-2
                         F2 = 3-5
                         CONSTRAIN = (1, 3, a1, a2), (5, 2, a2, a1), (1-3, d)')
    model9 <- mirt.model('F1 = 1-5
                         LBOUND = (1-3, g, 0.2), (4,5, g, 0.2)')
    model10 <- mirt.model('F1 = 1-5
                          START = (1,3-4, a1, 1)
                          FIXED = (1-3, a1)')
    model11 <- mirt.model('F1 = 1-5
                          PRIOR = (1-5, g, expbeta, 10, 40)')
    model12 <- mirt.model('F1 = 1-5
                           F2 = 2-5
                          CONSTRAIN = (1, 2, a1, a2)')
    model13 <- mirt.model('F = 1-5
                   CONSTRAINB = (1-5,a1)
                   CONSTRAINB [male, female]= (1-5,d)')
    model14 <- mirt.model("F = 1-5
                          FIXED = (1-5, a1)
                          START = (1-5, a1, 1.0)
                          FREE = (GROUP, COV_11)")

    mod0 <- mirt(data, model0, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod0)$value, c(0.9879254,1.85606,0,1,1.080885,0.8079786,0,1,1.705801,1.804219,0,1,0.7651853,0.4859966,0,1,0.735798,1.854513,0,1,0,1),
                 tolerance = 1e-2)
    mod1 <- mirt(data, model1, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod1)$value, c(0.9879254,1.85606,0,1,1.080885,0.8079786,0,1,1.705801,1.804219,0,1,0.7651853,0.4859966,0,1,0.735798,1.854513,0,1,0,1),
                 tolerance = 1e-2)
    mod2 <- mirt(data, model2, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod2)$value, c(1.01114,1.86813,0,1,1.01114,0.7909392,0,1,1.01114,1.460878,0,1,1.01114,0.5214695,0,1,1.01114,1.992827,0,1,0,1),
                 tolerance = 1e-2)
    mod3 <- mirt(data, model3, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod3)$value, c(1.108067,1.912278,0,1,1.046681,0.7965584,0,1,1.046681,1.472833,0,1,1.046681,0.5258648,0,1,1.046681,2.008483,0,1,0,1),
                 tolerance = 1e-2)
    mod4 <- multipleGroup(data, model4, group=group, verbose = FALSE)
    expect_equal(mod2values(mod4)$value, c(0.6407961,1.630026,0,1,2.96864,1.630026,0,1,1.189039,1.572577,0,1,0.558073,0.4682583,0,1,0.5097135,1.804266,0,1,0,1,0.6407961,1.621704,0,1,2.96864,1.621704,0,1,1.189039,1.578345,0,1,0.558073,0.4883855,0,1,0.5097135,1.759616,0,1,0,1),
                 tolerance = 1e-2)
    mod5 <- multipleGroup(data, model5, group=group, verbose = FALSE)
    expect_equal(mod2values(mod5)$value, c(0.6711941,1.6994,0,1,1.366601,0.8824884,0,1,1.074449,1.473972,0,1,1.487291,0.5832004,0,1,0.4453427,1.760829,0,1,0,1,0.6711941,1.388356,0,1,1.366601,1.388356,0,1,1.074449,1.388356,0,1,1.487291,1.388356,0,1,0.4453427,1.388356,0,1,0,1),
                 tolerance = 1e-2)
    mod6 <- mirt(data, model6, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod6)$value, c(1.074887,0,1.902959,0,1,1.074887,0,0.8065663,0,1,0,1.00348,1.458001,0,1,0,1.00348,0.520351,0,1,0,1.00348,1.989104,0,1,0,0,1,0.939999,1),
                 tolerance = 1e-2)
    mod7 <- mirt(data, model7, verbose=FALSE, calcNull=FALSE)
    expect_equal(as.numeric(coef(mod7, simplify=TRUE)$items), c(-1.153815,-0.2728293,0,-1,0,0,1.5,1.740508,0.5660805,0.5774055,1.945502,0.9276134,1.818513,0.5439836,1.789067,0,0,0,0,0,1,1,1,1,1),
                 tolerance = 1e-2)
    mod8 <- mirt(data, model8, verbose=FALSE, calcNull=FALSE)
    expect_equal(as.numeric(coef(mod8, simplify=TRUE)$items), c(0.5501291,3.146379,0,0,0,0,0,0.5501291,0.4096282,3.146379,1.46666,1.46666,1.46666,0.4535099,3.685736,0,0,0,0,0,1,1,1,1,1),
                 tolerance = 1e-2)
    mod9 <- mirt(data, model9, '3PL', verbose=FALSE, calcNull=FALSE)
    expect_equal(as.numeric(coef(mod9, simplify=TRUE)$items), c(1.09262,1.819549,2.095646,0.8938963,0.8182848,1.587373,0.1118206,1.542427,0.0396478,1.595971,0.2,0.2901414,0.2,0.2,0.2,1,1,1,1,1),
                 tolerance = 1e-2)
    mod10 <- mirt(data, model10, '3PL', pars = 'values')
    expect_equal(mod10$value[mod10$name == 'a1'], c(1, 0.851, 1, 1, .851), tolerance = 1e-4)
    expect_equal(mod10$est[mod10$name == 'a1'], c(FALSE, FALSE, FALSE, TRUE, TRUE))
    mod11 <- mirt(data, model11, '3PL', verbose=FALSE)
    expect_equal(as.vector(unname(coef(mod11)[[1]])),
                 c(1.0767651, 1.6027628, 0.1871268, 1.0000000), tolerance = 1e-4)
    mod12 <- mirt(data, model12, verbose=FALSE)
    expect_equal(as.vector(unname(c(coef(mod12)[[1]], coef(mod12)[[2]]))),
                 c(1.397911,0,2.093339,0,1,0.807397,1.397911,0.9542606,0,1), tolerance = 1e-4)
    mod13 <- multipleGroup(data, model13, group=group2, verbose = FALSE)
    expect_equal(mod2values(mod13)$value, c(-0.01562703,1.057122,0,1,1.243694,-0.063944,0,1,0.3680028,0.6760485,0,1,-1.259412,-0.4828693,0,1,0.05418581,1.400352,0,1,0,1,-0.01562703,1.057122,0,1,1.243694,-0.063944,0,1,0.3680028,0.6760485,0,1,-1.259412,-0.4828693,0,1,0.05418581,1.400352,0,1,0,1,-0.01562703,13.08788,0,1,1.243694,12.94597,0,1,0.3680028,13.22963,0,1,-1.259412,12.96195,0,1,0.05418581,2.512473,0,1,0,1),
                 tolerance = 1e-2)
    mod14 <- mirt(data, model14, verbose = FALSE)
    expect_equal(mod2values(mod14)$value, c(1,1.868016,0,1,1,0.7908857,0,1,1,1.46078,0,1,1,0.5214175,0,1,1,1.99271,0,1,0,1.021912),
                 tolerance = 1e-2)

    data(data.read, package = 'sirt')
    dat <- data.read

    # syntax with variable names
    mirtsyn2 <- "
            F1 = A1,B2,B3,C4
            F2 = A1-A4,C2,C4
            MEAN = F1
            COV = F1*F1, F1*F2
            CONSTRAIN=(A2-A4,a2),(A3,C2,d)
            PRIOR = (C3,A2-A4,a2,lnorm, .2, .2),(B3,d,norm,0,.0001)"
    # create a mirt model
    mirtmodel <- mirt.model(mirtsyn2, itemnames=dat)
    # or equivelently:
    mirtmodel2 <- mirt.model(mirtsyn2, itemnames=colnames(dat))

    expect_true(all(mirtmodel$x == mirtmodel2$x))
    got <- matrix(c(c('F1', 'F2', "MEAN", 'COV', 'CONSTRAIN', 'PRIOR'),
                    c("1,6,7,12", "1-4,10,12","F1","F1*F1,F1*F2","(2-4,a2),(3,10,d)",
                      "(11,2-4,a2,lnorm,.2,.2),(7,d,norm,0,.0001)")), nrow = 6)
    expect_true(all(mirtmodel$x == got))

    mod <- mirt(dat, mirtsyn2, TOL = NaN)
    sv <- mod2values(mod)
    expect_true(all(sv$est == c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE)))
    expect_true(all(as.character(sv$prior.type) == c("none","none","none","none","none","none","lnorm","none","none","none","none","lnorm","none","none","none","none","lnorm","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","norm","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","lnorm","none","none","none","none","none","none","none","none","none","none","none","none","none")))

})

