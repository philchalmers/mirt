expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('LLTM', {
    set.seed(42)
    a <- matrix(rep(1,30))
    d <- rep(c(1,0, -1),each = 10)  # first easy, then medium, last difficult
    dat <- simdata(a, d, 1000, itemtype = '2PL')
    itemdesign <- data.frame(difficulty =
        factor(c(rep('easy', 10), rep('medium', 10), rep('hard', 10))))
    rownames(itemdesign) <- colnames(dat)

    lltm <- mirt(dat, itemtype = 'Rasch', SE=TRUE, verbose=FALSE,
        item.formula = ~ 0 + difficulty, itemdesign=itemdesign)
    expect_equal(as.numeric(coef(lltm, simplify=TRUE)$items[1:2,1]), c(0.9587247, 0.9587247), tolerance=1e-2)
    expect_equal(as.numeric(coef(lltm, printSE=TRUE)[[1]][2, 1:3]), c(0.03886556, 0.03898911, 0.03769205), tolerance=1e-2)
    expect_equal(extract.mirt(lltm, 'condnum'), 6.633144, tolerance=1e-2)

    # additional information for LLTM
    oo <- plot(lltm)
    expect_class(oo, 'trellis')
    ifit <- itemfit(lltm)
    expect_equal(ifit$S_X2[1:3], c(20.06072, 20.90161, 23.48163), tolerance=1e-2)
    eap <- fscores(lltm)
    expect_equal(eap[1:3], c(1.0141828, -0.2427042, -0.2427042), tolerance=1e-2)

    # using unconditional modeling for first four items
    itemdesign.sub <- itemdesign[5:nrow(itemdesign), , drop=FALSE]
    lltm.4 <- mirt(dat, itemtype = 'Rasch', SE=TRUE, verbose=FALSE,
        item.formula = ~ 0 + difficulty, itemdesign=itemdesign.sub)
    cfs <- coef(lltm.4, simplify=TRUE)$items
    expect_equal(as.vector(cfs[c(1,5), 1:3]), c(0, .9063842, numeric(4)), tolerance=1e-2)
    expect_equal(anova(lltm, lltm.4)$p[2], 0.04288353, tolerance=1e-2)
    m2 <- M2(lltm.4)
    expect_equal(m2$TLI, 1.001862, tolerance=1e-2)

    group <- factor(rep(c('G1', 'G2'), each=500))
    lltm.G <- multipleGroup(dat, group=group, itemtype = 'Rasch', SE=TRUE, verbose=FALSE,
                            item.formula = ~ 0 + difficulty, itemdesign=itemdesign)
    cfsG1 <- coef(lltm.G, simplify=TRUE)$G1$items
    cfsG2 <- coef(lltm.G, simplify=TRUE)$G2$items
    expect_equal(as.vector(cfsG1[1, 1:3]), c(1.009664, 0.000000, 0.000000), tolerance=1e-2)
    expect_equal(as.vector(cfsG2[1, 1:3]), c(0.9074375, 0.000000, 0.000000), tolerance=1e-2)
    fs <- fscores(lltm.G)
    expect_equal(fs[1:3], c(1.0170938, -0.2457873, -0.2457873), tolerance=1e-2)


})

test_that('MLTM', {
    set.seed(42)

    as <- matrix(rep(1,60), ncol=2)
    as[11:18,1] <- as[1:9,2] <- 0
    d1 <- rep(c(3,1),each = 6)  # first easy, then medium, last difficult for first trait
    d2 <- rep(c(0,1,2),times = 4)    # difficult to easy
    d <- rnorm(18)
    ds <- rbind(cbind(d1=NA, d2=d), cbind(d1, d2))
    dat <- simdata(as, ds, 2500,
                   itemtype = c(rep('dich', 18), rep('partcomp', 12)))

    # unconditional model
    syntax <- "theta1 = 1-9, 19-30
           theta2 = 10-30
           COV = theta1*theta2"
    itemtype <- c(rep('Rasch', 18), rep('PC1PL', 12))
    mod <- mirt(dat, syntax, itemtype=itemtype, verbose=FALSE)
    expect_equal(as.numeric(coef(mod, simplify=TRUE)$items[19:21, c('d1', 'd2')]),
                 c(2.86919500, 3.71533102, 3.23797689, 0.01323629, 0.83160425, 1.90030392), tolerance=1e-2)
    expect_equal(logLik(mod), -43860.17, tolerance=1e-2)

    # MLTM design only for PC1PL items
    itemdesign <- data.frame(t1_difficulty= factor(d1, labels=c('medium', 'easy')),
                             t2_difficulty=factor(d2, labels=c('hard', 'medium', 'easy')))
    rownames(itemdesign) <- colnames(dat)[19:30]

    # fit MLTM design, leaving first 18 items as 'Rasch' type
    mltm <- mirt(dat, syntax, itemtype=itemtype, itemdesign=itemdesign,
                 item.formula = list(theta1 ~ 0 + t1_difficulty,
                                     theta2 ~ 0 + t2_difficulty), SE=TRUE, verbose=FALSE)
    expect_equal(extract.mirt(mltm, 'condnum'), 36.33511, tolerance=1e-2)
    expect_equal(anova(mltm, mod)$p[2], 0.592, tolerance=1e-2)
    cfs <- coef(mltm, simplify=TRUE)$items
    expect_equal(sort(unique(as.vector(cfs[,1:5]))),
                 c(-0.07838095, 0.00000000,  0.92424686,  1.03069049,  1.85666146,  3.18998660), tolerance=1e-2)
    fs <- fscores(mltm)
    expect_equal(as.vector(fs[1:3,]), c(-2.0019607,1.138449,0.149316,
                                        -0.4814751,0.3697978,1.7928627), tolerance=1e-2)
    m2 <- M2(mltm)
    expect_equal(m2$CFI, 0.9758302, tolerance=1e-2)


})

test_that('MLTM', {
    set.seed(42)

    as <- matrix(rep(1,60), ncol=2)
    as[11:18,1] <- as[1:9,2] <- 0
    d1 <- rep(c(3,1),each = 6)  # first easy, then medium, last difficult for first trait
    d2 <- rep(c(0,1,2),times = 4)    # difficult to easy
    d <- rep(-1:1, each=6)
    ds <- rbind(cbind(d1=NA, d2=d), cbind(d1, d2))
    dat <- simdata(as, ds, 2500,
                   itemtype = c(rep('dich', 18), rep('partcomp', 12)))

    # MLTM design only for PC1PL items
    itemdesign <- data.frame(t1_difficulty= factor(c(rep(NA, 18), d1), labels=c('medium', 'easy')),
                             t2_difficulty=factor(c(rep(NA, 18), d2), labels=c('hard', 'medium', 'easy')),
                             difficulty=factor(c(d, rep(NA, 12)), labels=c('hard', 'medium', 'easy')))

    # mixed <- mirt(dat, syntax, itemtype=itemtype, itemdesign=itemdesign,
    #              item.formula = list(theta1 ~ 0 + t1_difficulty,
    #                                  theta2 ~ 0 + t2_difficulty,
    #                                  ~ difficulty), SE=TRUE, verbose=FALSE)


})