context('DCIRT')

test_that('DCIRT', {

    set.seed(1234)
    N <- 1000
    P <- 25

    # First, sample item parameters:
    a <- matrix(rlnorm(25, .2, .3))
    b <- matrix(rnorm(25, 0, 1.2))
    d <- -a*b # IRT -> FA (mirt)

    # Then, sample latent traits and simulate data:
    bimodal.woodslin <- c(rnorm(N*.6, mean = -.70, sd = .50), rnorm(N*.4, mean = 1.05, sd = .54))
    dat_bm <- simdata(a, d, itemtype = 'dich', Theta = as.matrix(bimodal.woodslin))

    mod <- mirt(dat_bm, 1, dentype = 'empiricalhist_Woods', verbose=FALSE,
                technical = list(zeroExtreme = TRUE))
    expect_equal(extract.mirt(mod, 'logLik'), -12709.16, tolerance=1e-4)
    fs1 <- fscores(mod)
    fs2 <- fscores(mod, use_dentype_estimate = TRUE)
    expect_equal(fs1[1:3], c(-1.5345966, -0.8099821, -1.2411356), tolerance=1e-4)
    expect_equal(fs2[1:3], c(-1.3336196, -0.9165531, -1.0881438), tolerance=1e-4)
    fs1 <- fscores(mod, method = 'EAPsum')
    fs2 <- fscores(mod, method = 'EAPsum', use_dentype_estimate = TRUE)
    expect_equal(fs1[1:3], c(-1.4499389, -0.8284078, -1.4499389), tolerance=1e-4)
    expect_equal(fs2[1:3], c(-1.236668, -0.919067, -1.236668), tolerance=1e-4)
    resid <- residuals(mod, verbose=FALSE)
    resid2 <- residuals(mod, use_dentype_estimate=TRUE, verbose=FALSE)
    expect_equal(unname(resid[2:4,1]), c(3.454542, 5.869551, 0.119278), tolerance=1e-4)
    expect_equal(unname(resid2[2:4,1]), c(2.1035571, 1.9256263, 0.1226077), tolerance=1e-4)

    pp <- plot(mod, type = 'empiricalhist')
    expect_is(pp, 'trellis')

    mod2 <- mirt(dat_bm, 1, dentype = 'EHW', verbose=FALSE)
    expect_equal(extract.mirt(mod2, 'logLik'), -12758.689, tolerance=1e-4)

    res_bm <- mirt(dat_bm, model = 1, dentype='Davidian-6', verbose=FALSE)
    expect_equal(extract.mirt(res_bm, 'logLik'), -12769.95, tolerance=1e-4)
    cfs <- coef(res_bm)$GroupPars
    expect_equal(as.vector(cfs), c(0,1,1.331314,0.1005421,-0.397783,0.3760957,-0.7119315,-0.8868246), tolerance=1e-4)

    cfs2 <- coef(res_bm, simplify=TRUE)
    expect_equal(c(0, 1, cfs2$Davidian_phis), as.vector(cfs))

    fs1 <- fscores(res_bm)
    fs2 <- fscores(res_bm, use_dentype_estimate = TRUE)
    expect_equal(fs1[1:3], c(-1.5086792, -0.8190973, -1.2681388), tolerance=1e-4)
    expect_equal(fs2[1:3], c(-1.3837530, -0.8473187, -1.2122867), tolerance=1e-4)
    fs1 <- fscores(res_bm, method = 'EAPsum')
    fs2 <- fscores(res_bm, method = 'EAPsum', use_dentype_estimate = TRUE)
    expect_equal(fs1[1:2], c(-1.4400090, -0.8319593), tolerance=1e-4)
    expect_equal(fs2[1:2], c(-1.3283695, -0.8591202), tolerance=1e-4)

    out <- itemfit(res_bm)
    out2 <- itemfit(res_bm, use_dentype_estimate=TRUE)
    expect_equal(out$S_X2[1:3], c(25.04989, 14.21401, 14.81014), tolerance=1e-4)
    expect_equal(out2$S_X2[1:3], c(25.44769, 14.06283, 13.88460), tolerance=1e-4)

    pp <- plot(res_bm, type = 'Davidian')
    expect_is(pp, 'trellis')
})

test_that('DCIRT Option Errors and Warnings', {

    set.seed(1234)
    N <- 1000
    P <- 25

    # First, sample item parameters:
    a <- matrix(rlnorm(25, .2, .3))
    b <- matrix(rnorm(25, 0, 1.2))
    d <- -a*b # IRT -> FA (mirt)

    # Then, sample latent traits and simulate data:
    bimodal.woodslin <- c(rnorm(N*.6, mean = -.70, sd = .50), rnorm(N*.4, mean = 1.05, sd = .54))
    dat_bm <- simdata(a, d, itemtype = 'dich', Theta = as.matrix(bimodal.woodslin))

    # Err from in dcurver:: ?
    expect_error(mirt(dat_bm, model = 1, dentype='Davidian-12', verbose=FALSE))

    # Estimation methods should fail:
    expect_error(mirt(dat_bm, model = 1, dentype='Davidian-6', method = "QMCEM", verbose=FALSE))
    expect_error(mirt(dat_bm, model = 1, dentype='Davidian-6', method = "MCEM", verbose=FALSE))
    expect_error(mirt(dat_bm, model = 1, dentype='Davidian-6', method = "MHRM", verbose=FALSE))
    expect_error(mirt(dat_bm, model = 1, dentype='Davidian-6', method = "SEM", verbose=FALSE))
    expect_error(mirt(dat_bm, model = 1, dentype='Davidian-6', method = "BL", verbose=FALSE))

})

# test_that('DCIRT Convergence', {
#
#     # An inconvergent data is needed.
#     set.seed(1234)
#     N <- 1000
#     P <- 25
#
#     # First, sample item parameters:
#     a <- matrix(rlnorm(25, .2, .3))
#     b <- matrix(rnorm(25, 0, 1.2))
#     d <- -a*b # IRT -> FA (mirt)
#
#     # Then, sample latent traits and simulate data:
#     bimodal.woodslin <- c(rnorm(N*.6, mean = -.70, sd = .50), rnorm(N*.4, mean = 1.05, sd = .54))
#     dat_bm <- simdata(a, d, itemtype = 'dich', Theta = as.matrix(bimodal.woodslin))
#
#     # Check inconvergent DCs
#
#     # Use alternative/random starting points
#
# })

test_that('DCIRT-MG', {
    set.seed(1234)
    N <- 1000
    P <- 25

    # First, sample item parameters:
    a <- matrix(rlnorm(25, .2, .3))
    b <- matrix(rnorm(25, 0, 1.2))
    d <- -a*b # IRT -> FA (mirt)

    # Then, sample latent traits and simulate data:
    bimodal.woodslin <- c(rnorm(N*.6, mean = -.70, sd = .50), rnorm(N*.4, mean = 1.05, sd = .54))
    dat_bm <- simdata(a, d, itemtype = 'dich', Theta = as.matrix(bimodal.woodslin))
    bimodal.woodslin2 <- c(rnorm(N*.6, mean = -.70, sd = .50), rnorm(N*.4, mean = 1.05, sd = .54))*.75 + .5
    dat_bm2 <- simdata(a, d, itemtype = 'dich', Theta = as.matrix(bimodal.woodslin2))
    dat <- rbind(dat_bm, dat_bm2)
    colnames(dat) <- paste0("Item", 1:ncol(dat))
    group <- rep(c('G1', 'G2'), each = nrow(dat_bm))

    mod_configural <- multipleGroup(dat, 1, group = group, dentype='Davidian-6',
                                    verbose = FALSE)
    # plot(mod_configural)
    expect_equal(extract.mirt(mod_configural, 'logLik'), -23882.37, tolerance=1e-4)
    cfs <- coef(mod_configural)$G2$GroupPars
    expect_equal(as.vector(cfs), c(0,1,1.414418,0.08941072,-0.2462128,0.5657355,-0.8215269,-0.9705975), tolerance=1e-1)

    cfs2 <- coef(mod_configural, simplify=TRUE)
    expect_equal(c(0, 1, cfs2$G2$Davidian_phis), as.vector(cfs))

    fs1 <- fscores(mod_configural)
    fs2 <- fscores(mod_configural, use_dentype_estimate = TRUE)
    expect_equal(fs1[1:3], c(-1.510665,-0.8180199,-1.268329), tolerance=1e-1)
    expect_equal(fs2[1:3], c(-1.385374,-0.8463518,-1.212446), tolerance=1e-1)
    fs1 <- fscores(mod_configural, method = 'EAPsum')
    fs2 <- fscores(mod_configural, method = 'EAPsum', use_dentype_estimate = TRUE)
    expect_equal(fs1[1:2], c(-1.441303,-0.8309972), tolerance=1e-1)
    expect_equal(fs2[1:2], c(-1.3294050, -0.8582679), tolerance=1e-1)

    out <- itemfit(mod_configural)
    out2 <- itemfit(mod_configural, use_dentype_estimate=TRUE)
    expect_equal(out$G2$S_X2[1:3], c(10.11545,15.17313,8.989033), tolerance=1e-1)
    expect_equal(out2$G2$S_X2[1:3], c(10.2253,15.41447,9.010604), tolerance=1e-1)

    pp <- plot(mod_configural, type = 'Davidian')
    expect_is(pp, 'trellis')

    # not equated
    mod_scalar0 <- multipleGroup(dat, 1, group = group, verbose=FALSE, dentype='Davidian-6',
                                invariance=colnames(dat)[1:5])
    # coef(mod_scalar0, simplify=TRUE)
    expect_equal(extract.mirt(mod_scalar0, 'logLik'), -23968.56, tolerance=1e-4)
    expect_equal(extract.mirt(mod_scalar0, 'df'), 67108760)
    # plot(mod_scalar0)
    pp <- plot(mod_scalar0, type = 'Davidian')
    expect_is(pp, 'trellis')

    # equated
    expect_error(multipleGroup(dat, 1, group = group, verbose=FALSE, dentype='Davidian-6',
                                 invariance=c(colnames(dat)[1:10], 'free_var','free_means')))
    # plot(mod_scalar)
    # expect_equal(extract.mirt(mod_scalar, 'df'), 33554337)
    # expect_equal(extract.mirt(mod_scalar, 'logLik'), -48247.41, tolerance=1e-4)
    # expect_equal(M2(mod_scalar)$M2, 600.9787, tolerance=1e-4)
    # cfs <- as.vector(unname(coef(mod_scalar)$G2$GroupPars))
    # expect_equal(cfs, c(0.4745419,0.6083134,1.200818,2.058516,2.109649,3.667984,-0.6940906,0.7198784), tolerance=1e-4)
    #
    # pp <- plot(mod_scalar, type = 'Davidian')
    # expect_is(pp, 'trellis')

    # equated EHW
    mod_scalarEHW <- multipleGroup(dat, 1, group = group, verbose=FALSE, dentype='EHW',
                                invariance=c(colnames(dat)[1:10], 'free_var','free_means'))
    # plot(mod_scalarEHW)
    expect_equal(extract.mirt(mod_scalarEHW, 'df'), 67108544)
    expect_equal(extract.mirt(mod_scalarEHW, 'logLik'), -23877.63, tolerance=1e-4)
    cfs <- as.vector(unname(coef(mod_scalarEHW)$G2$GroupPars))
    expect_equal(cfs, c(0.5093032, 0.6000460), tolerance=1e-4)

    pp <- plot(mod_scalarEHW, type = 'empiricalhist')
    expect_is(pp, 'trellis')


})


