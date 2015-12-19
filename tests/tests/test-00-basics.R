context('basics')

test_that('basics', {
    set.seed(1)
    group <- sample(c('G1', 'G2'), nrow(Science), TRUE)
    gmod <- multipleGroup(Science, 1, group=group, TOL = NaN)
    mod <- extract.group(gmod, 1)
    item1 <- extract.item(gmod, 2, group=1)
    item2 <- extract.item(mod, 2)
    Theta <- matrix(-3:3)
    E <- expected.item(item1, Theta)
    P <- probtrace(item1, Theta)
    expect_equal(E, c(0.7007354,1.050795,1.411307,1.74344,2.036256,2.312032,2.56886), tolerance=1e-6)
    expect_equal(P[1:8], c(0.4615411,0.2679321,0.1351536,0.06255336,0.02770248,0.01201941,
                           0.005167734,0.3854359), tolerance=1e-6)
    theta_se <- fscores(mod, full.scores.SE = TRUE)
    expect_equal(unname(empirical_rxx(theta_se)), 0.509281, tolerance = 1e-4)
    tscore <- expected.test(mod, Theta)
    expect_equal(tscore, c(7.935196,9.270222,10.53705,11.69166,12.7803,13.83405,14.73532), tolerance=1e-6)
    IG <- itemGAM(Science[group == 'G1',1], theta_se[,1, drop=FALSE])
    expect_is(IG, 'itemGAM')

    info <- iteminfo(item1, Theta)
    expect_equal(as.vector(info), c(0.2031261,0.2118485,0.2069816,0.187747,0.1793368,0.1902587,0.1761818),
                 tolerance=1e-6)
    info2 <- iteminfo(item1, Theta, total.info = FALSE)
    expect_equal(as.vector(info2[1:2, ]), c(0.09691134,0.1039888,0.02656885,0.0002719432,
                 0.07306797,0.09274217,0.006577983,0.01484559), tolerance = 1e-4)
    rxx <- marginal_rxx(mod)
    expect_equal(rxx, 0.4209159, tolerance=1e-6)
    MD <- unname(MDISC(mod)[1])
    expect_equal(MD, 0.724201, tolerance=1e-6)
    MDF <- unname(MDIFF(mod)[1L,])
    expect_equal(MDF, c(-6.061399, -3.566437, 2.031588), tolerance=1e-6)
})

