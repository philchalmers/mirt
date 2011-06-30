\name{confmirtClass-class}
\docType{class}
\alias{confmirtClass-class}
\alias{coef,confmirtClass-method}
\alias{print,confmirtClass-method}
\alias{residuals,confmirtClass-method}
\alias{show,confmirtClass-method}
\alias{summary,confmirtClass-method}
\alias{anova,confmirtClass-method}

\title{Class "confmirtClass"}
\description{
Defines the object returned from \code{\link{confmirt}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("confmirtClass", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{pars}:}{Object of class \code{"matrix"}, parameters }
    \item{\code{guess}:}{Object of class \code{"numeric"}, guessing parameters }
    \item{\code{SEpars}:}{Object of class \code{"matrix"},  standard errors for pars}
    \item{\code{SEg}:}{Object of class \code{"numeric"}, standard errors for guess }
    \item{\code{gpars}:}{Object of class \code{"list"}, group parameters }
    \item{\code{SEgpars}:}{Object of class \code{"list"}, standard errors for gpars }
    \item{\code{estpars}:}{Object of class \code{"list"}, indicators for gpars }
    \item{\code{cycles}:}{Object of class \code{"numeric"},  number of MHRM cycles}
    \item{\code{Theta}:}{Object of class \code{"matrix"}, final MH draw }
    \item{\code{fulldata}:}{Object of class \code{"matrix"},  binary data}
    \item{\code{data}:}{Object of class \code{"matrix"}, data with 99 as missing }
    \item{\code{K}:}{Object of class \code{"numeric"},  number of item categories}
    \item{\code{itemloc}:}{Object of class \code{"numeric"}, item location in fulldata }
    \item{\code{logLik}:}{Object of class \code{"numeric"}, log-likelihood }
    \item{\code{SElogLik}:}{Object of class \code{"numeric"}, std error for log-likelihood }
    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
    \item{\code{df}:}{Object of class \code{"integer"}, degrees of freedom }
    \item{\code{h2}:}{Object of class \code{"numeric"}, communalities }
    \item{\code{F}:}{Object of class \code{"matrix"}, factor loadings }
    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
    \item{\code{nconstvalues}:}{Object of class \code{"integer"}, fixed parameters count }
    \item{\code{G2}:}{Object of class \code{"numeric"}, G2 model fit statistic }
    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for G2 }
    \item{\code{Call}:}{Object of class \code{"call"},  call }
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "confmirtClass")} }
    \item{print}{\code{signature(x = "confmirtClass")} }
    \item{residuals}{\code{signature(object = "confmirtClass")}}
    \item{show}{\code{signature(object = "confmirtClass")} }
    \item{summary}{\code{signature(object = "confmirtClass")} }
    \item{logLik}{\code{signature(object = "confmirtClass")} }
    \item{anova}{\code{signature(object = "confmirtClass")} }
	 }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{classes}
