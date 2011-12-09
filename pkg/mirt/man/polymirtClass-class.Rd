\name{polymirtClass-class}
\docType{class}
\alias{polymirtClass-class}
\alias{coef,polymirtClass-method}
\alias{plot,polymirtClass,missing-method}
\alias{print,polymirtClass-method}
\alias{residuals,polymirtClass-method}
\alias{show,polymirtClass-method}
\alias{summary,polymirtClass-method}
\alias{anova,polymirtClass-method}
\alias{fitted,polymirtClass-method}

\title{Class "polymirtClass"}
\description{
Defines the object returned from \code{\link{polymirt}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("polymirtClass", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{pars}:}{Object of class \code{"matrix"}, parameters  }
    \item{\code{guess}:}{Object of class \code{"numeric"}, guessing parameters  }
    \item{\code{SEpars}:}{Object of class \code{"matrix"}, standard errors for pars }
    \item{\code{cycles}:}{Object of class \code{"numeric"}, number of MHRM cycles }
    \item{\code{Theta}:}{Object of class \code{"matrix"}, final MH draw }
    \item{\code{fulldata}:}{Object of class \code{"matrix"}, complete binary data }
    \item{\code{data}:}{Object of class \code{"matrix"}, complete data with 99 as missing }
    \item{\code{K}:}{Object of class \code{"numeric"}, number of item categories  }
    \item{\code{F}:}{Object of class \code{"matrix"}, factor loadings }
    \item{\code{h2}:}{Object of class \code{"numeric"}, commonalities }
    \item{\code{logLik}:}{Object of class \code{"numeric"}, log-likelihood }
    \item{\code{SElogLik}:}{Object of class \code{"numeric"}, std error for log-likelihood }
    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
    \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian Information Criteria }
    \item{\code{df}:}{Object of class \code{"integer"}, degrees of freedom }
    \item{\code{itemloc}:}{Object of class \code{"numeric"}, item locations in fulldata }
    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
    \item{\code{G2}:}{Object of class \code{"numeric"}, G2 model fit statistic }
    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for G2 }
    \item{\code{tabdata}:}{Object of class \code{"matrix"}, tabulated data with expected values }
    \item{\code{estGuess}:}{Object of class \code{"logical"}, vector indicating which guessing parameters were estimated }
    \item{\code{Call}:}{Object of class \code{"call"}, call }
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "polymirtClass")}}
    \item{plot}{\code{signature(x = "polymirtClass", y = "missing")} }
    \item{print}{\code{signature(x = "polymirtClass")} }
    \item{residuals}{\code{signature(object = "polymirtClass")} }
    \item{show}{\code{signature(object = "polymirtClass")} }
    \item{summary}{\code{signature(object = "polymirtClass")}}
    \item{logLik}{\code{signature(object = "polymirtClass")} }
    \item{anova}{\code{signature(object = "polymirtClass")} }
	\item{fitted}{\code{signature(object = "polymirtClass")} }
	 }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{classes}
