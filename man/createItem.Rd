\name{createItem}
\alias{createItem}
\title{Create a user defined item with correct generic functions}
\usage{
  createItem(name, par, est, P, gr = NULL, hss = NULL,
    lbound = NULL, ubound = NULL)
}
\arguments{
  \item{name}{a character indicating the item class name to
  be defined}

  \item{makeitem}{logical; after initializing the methods,
  create an instantiation?}

  \item{par}{a named vector of the starting values for the
  parameters}

  \item{est}{a logical vector indicating which parameters
  should be freely estimated by default}

  \item{P}{the probability trace function for all
  categories (first column is category 1, second category
  two, etc). Must return a \code{matrix} object}

  \item{gr}{gradient function (vector of first derivatives)
  used in estimation. If not specified a numeric
  approximation will be used}

  \item{hss}{hessian function (matrix of second
  derivatives) used in estimation. If not specified a
  numeric approximation will be used (required for the
  MH-RM algorithm only)}

  \item{lbound}{optional vector indicating the lower bounds
  of the parameters. If not specified then the bounds will
  be set to -Inf}

  \item{ubound}{optional vector indicating the lower bounds
  of the parameters. If not specified then the bounds will
  be set to Inf}
}
\description{
  Initializes the proper S4 class and methods necessary for
  mirt functions to use in estimation. To use the defined
  objects pass to the \code{mirt(..., customItems =
  list())} command, and ensure that the classes are
  properly labelled and unique in the list.
}
\details{
  Addionally, the \code{summary()} function will not return
  proper standardized loadings since the function is not
  sure how to handle them (no slopes could be defined at
  all!). Instead loadings of .001 are filled in as
  placeholders.
}
\examples{
\dontrun{

name <- 'old2PL'
par <- c(a = .5, b = -2)
est <- c(TRUE, TRUE)
P.old2PL <- function(par,Theta){
     a <- par[1]
     b <- par[2]
     P1 <- 1 / (1 + exp(-1.702*a*(Theta - b)))
     cbind(1-P1, P1)
}
lbound <- c(-Inf, -Inf)
ubound <- c(Inf, Inf)

x <- createItem(name, par=par, est=est, lbound=lbound, ubound=ubound, P=P.old2PL)

#So, let's estimate it!
dat <- expand.table(LSAT7)
sv <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), pars = 'values')
tail(sv) #looks good
mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose = TRUE)
coef(mod)
mod2 <- confmirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose = TRUE)
coef(mod2)
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{createItem}

