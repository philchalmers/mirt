\name{mirtClass-class}
\docType{class}
\alias{mirtClass-class}
\alias{anova,mirtClass-method}
\alias{coef,mirtClass-method}
\alias{fitted,mirtClass-method}
\alias{plot,mirtClass,missing-method}
\alias{print,mirtClass-method}
\alias{residuals,mirtClass-method}
\alias{show,mirtClass-method}
\alias{summary,mirtClass-method}

\title{Class "mirtClass"}
\description{
Defines the object returned from \code{\link{mirt}}. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mirtClass", ...).}.
}
\section{Slots}{
  \describe{
    \item{\code{EMiter}:}{Object of class \code{"numeric"}, number of iterations  }
    \item{\code{pars}:}{Object of class \code{"matrix"}, parameters  }
    \item{\code{guess}:}{Object of class \code{"numeric", guessing parameters}  }
    \item{\code{X2}:}{Object of class \code{"numeric"}, Chi-squared  }
    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom for X2 }
    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for X2  }
    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
    \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian Information Criteria }
    \item{\code{log.lik}:}{Object of class \code{"numeric"}, log-likelihood }
    \item{\code{F}:}{Object of class \code{"matrix"}, factor loadings }
    \item{\code{h2}:}{Object of class \code{"numeric"}, communalities }
    \item{\code{tabdata}:}{Object of class \code{"matrix"}, tabulated data }
    \item{\code{Theta}:}{Object of class \code{"matrix"}, ability grid }
    \item{\code{Pl}:}{Object of class \code{"numeric"}, p's for tabulated response}
    \item{\code{fulldata}:}{Object of class \code{"matrix"}, complete data }
    \item{\code{cormat}:}{Object of class \code{"matrix"}, quasi-correlation matrix }
    \item{\code{facility}:}{Object of class \code{"numeric"}, item easiness }
    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
    \item{\code{quadpts}:}{Object of class \code{"numeric"}, number of quadrature points }
    \item{\code{Call}:}{Object of class \code{"call"}, call }
  }
}
\section{Methods}{
  \describe{
    \item{anova}{\code{signature(object = "mirtClass")}}
    \item{coef}{\code{signature(object = "mirtClass")}}
    \item{fitted}{\code{signature(object = "mirtClass")}}
    \item{plot}{\code{signature(x = "mirtClass", y = "missing")}}
    \item{print}{\code{signature(x = "mirtClass")} }
    \item{residuals}{\code{signature(object = "mirtClass")}}
    \item{show}{\code{signature(object = "mirtClass")} }
    \item{summary}{\code{signature(object = "mirtClass")}}
	 }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{classes}
