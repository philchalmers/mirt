\docType{data}
\name{LSAT7}
\alias{LSAT7}
\title{Description of LSAT7 data}
\description{
Data from Bock & Lieberman (1970); contains 5 dichotomously
scored items obtained from the Law School Admissions Test,
section 7.
}
\examples{
\dontrun{
dat <- expand.table(LSAT7)
head(dat)
(mod <- mirt(dat, 1))
coef(mod)
}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
Bock, R. D., & Lieberman, M. (1970). Fitting a response
model for \emph{n} dichotomously scored items.
\emph{Psychometrika, 35}(2), 179-197.
}
\keyword{data}

