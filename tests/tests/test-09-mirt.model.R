context('mirt.model')

test_that('syntax', {
    data <- expand.table(LSAT7)
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

    set.seed(1234)
    group <- sample(c('male', 'female'), 1000, TRUE)
    group2 <- sample(c('male', 'female', 'other'), 1000, TRUE)

    mod0 <- mirt(data, model0, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod0)$value, c(0.987973787231699, 1.85608912732841, 0, 1, 1.08103954211169, 0.808007534786952, 0, 1, 1.70595475896956, 1.80426768080187, 0, 1, 0.765076394253259, 0.486005938565521, 0, 1, 0.735771996169788, 1.85448564531374, 0, 1, 0, 1),
                 tolerance = 1e-2)
    mod1 <- mirt(data, model1, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod1)$value, c(0.987973787231699, 1.85608912732841, 0, 1, 1.08103954211169, 0.808007534786952, 0, 1, 1.70595475896956, 1.80426768080187, 0, 1, 0.765076394253259, 0.486005938565521, 0, 1, 0.735771996169788, 1.85448564531374, 0, 1, 0, 1),
                 tolerance = 1e-2)
    mod2 <- mirt(data, model2, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod2)$value, c(1.01052705474606, 1.86793304554809, 0, 1, 1.01052705474606, 0.790899504740889, 0, 1, 1.01052705474606, 1.46073200902294, 0, 1, 1.01052705474606, 0.521457445159032, 0, 1, 1.01052705474606, 1.99261823763434, 0, 1, 0, 1),
                 tolerance = 1e-2)
    mod3 <- mirt(data, model3, verbose=FALSE, calcNull=FALSE)
    expect_equal(mod2values(mod3)$value, c(1.10784174676958, 1.91239633478765, 0, 1, 1.04684043398884, 0.797880409701628, 0, 1, 1.04684043398884, 1.47310584064986, 0, 1, 1.04684043398884, 0.526104887330641, 0, 1, 1.04684043398884, 2.00877254688793, 0, 1, 0, 1),
                 tolerance = 1e-2)
    mod4 <- multipleGroup(data, model4, group=group, verbose = FALSE)
    expect_equal(mod2values(mod4)$value, c(0.6329101,1.460023,0,1,3.05139,1.460023,0,1,1.172012,1.489809,0,1,0.5481464,0.3904254,0,1,0.4969134,1.650933,0,1,0,1,0.6329101,1.809123,0,1,3.05139,1.809123,0,1,1.172012,1.635588,0,1,0.5481464,0.5571859,0,1,0.4969134,1.903352,0,1,0,1),
                 tolerance = 1e-2)
    mod5 <- multipleGroup(data, model5, group=group, verbose = FALSE)
    expect_equal(mod2values(mod5)$value, c(0.6146542,1.474662,0,1,1.433285,0.8417281,0,1,1.091704,1.425829,0,1,1.546467,0.5011313,0,1,0.4210692,1.614541,0,1,0,1,0.6146542,1.508259,0,1,1.433285,1.508259,0,1,1.091704,1.508259,0,1,1.546467,1.508259,0,1,0.4210692,1.508259,0,1,0,1),
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
    expect_equal(mod2values(mod13)$value, c(0.9911819,1.936458,0,1,1.08577,0.7493856,0,1,1.705515,1.770363,0,1,0.7650503,0.5171894,0,1,0.7388148,1.960346,0,1,0,1,0.9911819,1.936458,0,1,1.08577,0.7493856,0,1,1.705515,1.770363,0,1,0.7650503,0.5171894,0,1,0.7388148,1.960346,0,1,0,1,0.9911819,1.706259,0,1,1.08577,0.9337099,0,1,1.705515,1.874088,0,1,0.7650503,0.4229254,0,1,0.7388148,1.660745,0,1,0,1),
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

