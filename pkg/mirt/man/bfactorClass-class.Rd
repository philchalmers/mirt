\name{bfactorClass-class}
\docType{class}
\alias{bfactorClass-class}
\alias{coef,bfactorClass-method}
\alias{fitted,bfactorClass-method}
\alias{print,bfactorClass-method}
\alias{residuals,bfactorClass-method}
\alias{show,bfactorClass-method}
\alias{summary,bfactorClass-method}

\title{Class "bfactorClass"}
\description{
Defines the object returned from \code{\link{bfactor}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bfactorClass", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{EMiter}:}{Object of class \code{"numeric"}, number of EM iterations  }
    \item{\code{pars}:}{Object of class \code{"matrix"}, estimated parameters }
    \item{\code{guess}:}{Object of class \code{"numeric"}, guessing parameters }
    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
    \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian information criteria }
    \item{\code{X2}:}{Object of class \code{"numeric"}, chi-squared value }
    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom }
    \item{\code{logLik}:}{Object of class \code{"numeric"}, observed log-likelihood }
    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for chi-squared value}
    \item{\code{F}:}{Object of class \code{"matrix"},  factor loadings}
    \item{\code{h2}:}{Object of class \code{"numeric"}, commonalities }
    \item{\code{itemnames}:}{Object of class \code{"character"},  item names}
    \item{\code{tabdata}:}{Object of class \code{"matrix"},  tabulated data}
    \item{\code{N}:}{Object of class \code{"numeric"},  sample size}
    \item{\code{Pl}:}{Object of class \code{"numeric"},  likelihood-values for tabdata}
    \item{\code{Theta}:}{Object of class \code{"matrix"}, theta grid }
    \item{\code{fulldata}:}{Object of class \code{"matrix"},  complete binary data}
    \item{\code{logicalfact}:}{Object of class \code{"matrix"}, logical factor indicators }
    \item{\code{facility}:}{Object of class \code{"numeric"},  item easiness}
    \item{\code{specific}:}{Object of class \code{"numeric"}, specific item indicators }
    \item{\code{cormat}:}{Object of class \code{"matrix"}, quasi-correlation matrix }
    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
    \item{\code{par.prior}:}{Object of class \code{"matrix"},  prior parameter distribution specifications}
    \item{\code{quadpts}:}{Object of class \code{"numeric"}, number for quadrature points }
    \item{\code{Call}:}{Object of class \code{"call"},  call }
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "bfactorClass")} }
    \item{fitted}{\code{signature(object = "bfactorClass")} }
    \item{print}{\code{signature(x = "bfactorClass")} }
    \item{residuals}{\code{signature(object = "bfactorClass")}}
    \item{show}{\code{signature(object = "bfactorClass")} }
    \item{summary}{\code{signature(object = "bfactorClass")} }
	 }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{classes}
