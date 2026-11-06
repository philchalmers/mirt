expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('High Correlations', {
    set.seed(1234)
    a <- matrix(c(rlnorm(10, .2, .2), numeric(30), rlnorm(10, .2, .2),  numeric(30),rlnorm(10, .2, .2)), ncol=3)
    d <- rnorm(30)
    sigma <- matrix(.95, 3, 3)
    diag(sigma) <- 1
    dat <- simdata(a, d, 5000, itemtype = '2PL', sigma=sigma)
    model <- "F1 = 1-10
                  F2 = 11-20
                  F3 = 21-30
                  COV = F1*F2*F3"

    mod.MHRM <- mirt(dat, model, method='MHRM', technical=list(warn=FALSE))
    expect_equal(logLik(mod.MHRM), -82020.94, tolerance=1e-4)

    mod.EM <- mirt(dat, model, method='EM', TOL=.001, quadpts=26, technical=list(warn=FALSE))
    expect_equal(logLik(mod.EM), -82017.93, tolerance=1e-4)

    mod.EM_start <- mirt(dat, model, method='EM', TOL=.001, quadpts=31,
                         technical=list(warn=FALSE, MAXQUAD=30000),
                         pars = mod2values(mod.MHRM))
    expect_equal(logLik(mod.EM_start), -82016.68, tolerance=1e-4)

    mod.QMCEM <- mirt(dat, model, method='QMCEM', TOL=.001, technical=list(warn=FALSE))
    expect_equal(logLik(mod.QMCEM), -82007.65, tolerance=1e-4)

    set.seed(4321)
    mod.MCEM <- mirt(dat, model, method='MCEM', technical=list(warn=FALSE, NCYCLES=100))
    expect_equal(logLik(mod.MCEM), -82022.53, tolerance=1e-4)




})