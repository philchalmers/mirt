\name{specify.model}
\alias{specify.model}
\title{Specify confirmatory loadings}
\description{
The \code{specify.model} function scans user input to specify the confirmatory factor loadings. This function is nearly identical to \code{\link[sem]{specify.model}}, except that residuals terms are omitted.
}

\usage{
specify.model()
}

\value{
  Returns a model specification object to be used in \code{\link{confmirt}}.
}

\examples{
\dontrun{
model <- specify.model()
    F1 -> Item_1,   lam11, NA
    F1 -> Item_2,   lam12, NA
    F1 -> Item_3,   lam13, NA
    F1 -> Item_4,   lam14, NA
    F2 -> Item_4,   lam24, NA
    F2 -> Item_5,   lam25, NA
    F2 -> Item_6,   lam26, NA
    F2 -> Item_7,   lam27, NA
    F2 -> Item_8,   lam28, NA
    F1 <-> F2,      phi12, NA
    F1 <-> F1,      NA,     1
    F2 <-> F2,      NA,     1

    }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}


