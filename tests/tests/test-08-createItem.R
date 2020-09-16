context('createItem')

test_that('old2PL', {
    name <- 'old2PL'
    par <- c(a = .5, b = -2)
    est <- c(TRUE, TRUE)
    P.old2PL <- function(par,Theta,ncat){
        a <- par[1]
        b <- par[2]
        P1 <- 1 / (1 + exp(-1*a*(Theta - b)))
        cbind(1-P1, P1)
    }
    lbound <- c(-Inf, -Inf)
    ubound <- c(Inf, Inf)

    x <- createItem(name, par=par, est=est, lbound=lbound, ubound=ubound, P=P.old2PL)

    dat <- expand.table(LSAT7)
    sv <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), pars = 'values', verbose=FALSE)
    expect_is(sv, 'data.frame')
    mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose=FALSE)
    expect_is(mod, 'SingleGroupClass')
    expect_is(coef(mod), 'list')
    expect_equal(logLik(mod), -2658.805, tolerance = 1e-4)
    mod2 <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod2, 'SingleGroupClass')
    expect_is(coef(mod2), 'list')
    mod3 <- mirt(dat, 'F = 1-5
                       PRIOR = (5, b, norm, 0, 1)', c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose=FALSE)
    expect_equal(logLik(mod3), -2659.1, tolerance = 1e-4)


    #' #nonlinear
    name <- 'nonlin'
    par <- c(a1 = .5, a2 = .1, d = 0)
    est <- c(TRUE, TRUE, TRUE)
    P.nonlin <- function(par,Theta,ncat){
      a1 <- par[1]
      a2 <- par[2]
      d <- par[3]
      P1 <- 1 / (1 + exp(-1*(a1*Theta + a2*Theta^2 + d)))
      cbind(1-P1, P1)
    }

    x2 <- createItem(name, par=par, est=est, P=P.nonlin)
    mod <- mirt(dat, 1, c(rep('2PL',4), 'nonlin'), customItems=list(nonlin=x2), verbose=FALSE)
    expect_is(mod, 'SingleGroupClass')
    expect_is(coef(mod), 'list')

    fs <- fscores(mod)
    expect_equal(unname(fs[1,]), c(-1.828444), tolerance = 1e-4)
})


