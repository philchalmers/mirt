library(testthat)
library(mirt)
library(numDeriv)
options(error = utils::recover, browserNLdisabled=TRUE)

i.count <- 8
data <- simdata(as.matrix(rlnorm(i.count, sdlog=.5)),
                as.matrix(rnorm(i.count)),
                100, 'dich')
spoint <- c(1,0,.2,.9)

deriv <- genD(function(param) {
  tech <- list(debugderiv=list(param=1:4, value=param))
  fit <- mirt(data, 1, rep('4PL',i.count), D=1, calcNull=FALSE, technical=tech)
  fit@MstepLogLik
}, spoint, method.args=list(eps=0.01, d=0.01, r=2))

tech <- list(debugderiv=list(param=1:4, value=spoint))
fit <- mirt(data, 1, rep('4PL',i.count), D=1, calcNull=FALSE, technical=tech)

if (0) {
  for (ix in 1:i.count) {
    ii <- extract.item(fit, ix)
    print(ii@par)
  }
}

if (1) {
  ii <- extract.item(fit, 1)

  np <- length(spoint)
  expect_equal(-ii@gradient, deriv$D[1:np], tolerance=1e-4)
  
  hess <- matrix(NA, nrow=np, ncol=np)
  dx <- np+1
  for (hr in 1:np) {
    hess[hr,1:hr] <- deriv$D[dx:(dx+hr-1)]
    dx <- dx + hr
  }
  ii@hessian[is.na(hess)] <- NA
  expect_equal(hess, -ii@hessian, tolerance=1e-4)
}
