\name{DTF}
\alias{DTF}
\title{Differential test functioning}
\usage{
DTF(mod, MI = NULL, CI = 0.95, npts = 200, digits = 4,
  theta_lim = c(-6, 6))
}
\arguments{
  \item{mod}{a multipleGroup object which estimated only 2
  groups}

  \item{MI}{a number indicating how many draws to take to
  form a suitable multiple imputation for the expected test
  scores (100 or more). Requires an estimated parameter
  information matrix}

  \item{CI}{range of condfidince interval when using MI}

  \item{npts}{number of points to use in the integration,
  ranging from -6 to 6}

  \item{digits}{number of digits to round result to}
}
\description{
Function performs various omnibus differential test
functioning proceduces on an object estimated with
\code{multipleGroup()}. If the latent means/covariances are
suspected to differ then the input object should contain a
set of 'anchor' items to ensure that only differential test
features are being detected rather than group differences.
Returns signed (average area above and below) and unsigned
(total area) statistics, with descriptives such as the
percent average bias between group total scores.
}
\examples{
\dontrun{
set.seed(1234)
n <- 30
N <- 500

# only first 5 items as anchors
model <- mirt.model('F = 1-30
                    CONSTRAINB = (1-5, a1), (1-5, d)')

a <- matrix(1, n)
d <- matrix(rnorm(n), n)
group <- c(rep('G1', N), rep('G2', N))

## -------------
# groups completely equal
dat1 <- simdata(a, d, N, itemtype = 'dich')
dat2 <- simdata(a, d, N, itemtype = 'dich')
dat <- rbind(dat1, dat2)
mod <- multipleGroup(dat, model, group=group, SE=TRUE, SE.type='crossprod',
                     invariance=c('free_means', 'free_var'))
plot(mod, type = 'score')

DTF(mod)
DTF(mod, MI = 100) #95\% C.I. for DTI containins 0

## -------------
## random slopes and intercepts for 15 items, and latent mean difference
##    (no systematic DTF should exist, but DIF will be present)
dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
dat2 <- simdata(a + c(numeric(15), sign(rnorm(n-15))*runif(n-15, .25, .5)),
                d + c(numeric(15), sign(rnorm(n-15))*runif(n-15, .5, 1)), N, itemtype = 'dich')
dat <- rbind(dat1, dat2)
mod1 <- multipleGroup(dat, 1, group=group)
plot(mod1, type = 'score') #does not account for group differences! Need anchors

mod2 <- multipleGroup(dat, model, group=group, SE=TRUE, SE.type = 'crossprod',
                      invariance=c('free_means', 'free_var'))
plot(mod2, type = 'score')

mirtCluster()
#significant DIF in multiple items....
DIF(mod2, which.par=c('a1', 'd'), items2test=16:30)
DTF(mod2) #...but not substantial DTF due to randomness of DIF
DTF(mod2, MI=100)

## -------------
## systematic differing slopes and intercepts (clear DTF)
dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)), d + c(numeric(15), rnorm(n-15, 1, .5)),
                N, itemtype = 'dich')
dat <- rbind(dat1, dat2)
mod3 <- multipleGroup(dat, model, group=group, SE=TRUE, SE.type='crossprod',
                      invariance=c('free_means', 'free_var'))
plot(mod3, type = 'score') #visable DTF happening

DIF(mod3, c('a1', 'd'), items2test=16:30)
DTF(mod3) #huge unsigned bias. Signed bias indicates group 2 scores generally lower
DTF(mod3, MI=100)

}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
\code{\link{multipleGroup}}, \code{\link{DIF}}
}
\keyword{DTF}

