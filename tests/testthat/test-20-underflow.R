test_that('E-step likelihood underflow (long tests)', {

    # Independent marginal log-likelihood per response pattern via log-sum-exp
    indep_LL <- function(mod){
        Theta <- mod@Model$Theta
        pr <- mod@Internals$Prior[[1L]]
        if(is.matrix(pr)) pr <- pr[1L, ]
        ll <- mod@Data$fulldata[[1L]] %*% t(log(probtrace(mod, Theta)))
        apply(ll, 1L, function(x){
            x <- x + log(pr)
            m <- max(x)
            m + log(sum(exp(x - m)))
        })
    }

    # 400 steep items; person 1 misses every easy item and passes every hard one,
    # so their likelihood is < 1e-308 at every quadrature node (log-lik ~ -885)
    set.seed(1234)
    J <- 400
    d <- seq(-4, 4, length.out=J)
    dat <- simdata(matrix(rep(4, J)), matrix(d), 100, itemtype='dich')
    dat[1L, ] <- as.integer(d < 0)
    sv <- mirt(dat, 1, pars='values')
    sv$value[sv$name == 'a1'] <- 4
    sv$value[sv$name == 'd'] <- d
    sv$est <- FALSE

    mod <- mirt(dat, 1, pars=sv, verbose=FALSE)
    LL <- indep_LL(mod)
    expect_true(LL[1L] < log(.Machine$double.xmin))
    expect_true(is.finite(logLik(mod)))
    expect_equal(logLik(mod), sum(LL), tolerance=1e-8)
    expect_true(is.finite(extract.mirt(mod, 'G2')))
    expect_true(is.finite(extract.mirt(mod, 'AIC')))

    # full-data E-step (latent regression path)
    covdata <- data.frame(X=rep(0:1, 50))
    sv2 <- mirt(dat, 1, pars='values', covdata=covdata, formula=~X)
    pick <- sv2$item %in% colnames(dat)
    sv2$value[pick & sv2$name == 'a1'] <- 4
    sv2$value[pick & sv2$name == 'd'] <- d
    sv2$est <- FALSE
    mod2 <- mirt(dat, 1, pars=sv2, covdata=covdata, formula=~X, verbose=FALSE)
    expect_true(is.finite(logLik(mod2)))
    expect_equal(logLik(mod2), sum(LL), tolerance=1e-8)
})
