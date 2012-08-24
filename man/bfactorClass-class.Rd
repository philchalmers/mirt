\name{bfactorClass-class}
\docType{class}
\alias{bfactorClass-class}
\alias{coef,bfactorClass-method}
\alias{anova,bfactorClass-method}
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
    \item{\code{pars}:}{Object of class \code{"list"}, estimated parameter objects list }        
    \item{\code{K}:}{Object of class \code{"numeric", number of item categories}  }    
    \item{\code{itemloc}:}{Object of class \code{"numeric", index for tabdata}  }     
    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
    \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian information criteria }
    \item{\code{G2}:}{Object of class \code{"numeric"}, G squared stat }
    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom }
    \item{\code{RMSEA}:}{Object of class \code{"numeric"}, root mean-square error of approximation}
    \item{\code{TLI}:}{Object of class \code{"numeric"}, Tucker-Lewis index}
    \item{\code{logLik}:}{Object of class \code{"numeric"}, observed log-likelihood }
    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for chi-squared value}
    \item{\code{F}:}{Object of class \code{"matrix"},  factor loadings}
    \item{\code{h2}:}{Object of class \code{"numeric"}, commonalities }
    \item{\code{itemnames}:}{Object of class \code{"character"},  item names}
    \item{\code{tabdata}:}{Object of class \code{"matrix"},  tabulated data}
	\item{\code{tabdatalong}:}{Object of class \code{"matrix"}, dichotomous version of tabulated data }		
    \item{\code{N}:}{Object of class \code{"numeric"},  sample size}
    \item{\code{Pl}:}{Object of class \code{"numeric"},  likelihood-values for tabdata}
    \item{\code{Theta}:}{Object of class \code{"matrix"}, theta grid }
    \item{\code{data}:}{Object of class \code{"matrix"},  complete binary data}
    \item{\code{logicalfact}:}{Object of class \code{"matrix"}, logical factor indicators }    
    \item{\code{specific}:}{Object of class \code{"numeric"}, specific item indicators }
    \item{\code{cormat}:}{Object of class \code{"matrix"}, quasi-correlation matrix }
    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }    
    \item{\code{quadpts}:}{Object of class \code{"numeric"}, number for quadrature points }
    \item{\code{vcov}:}{Object of class \code{"matrix"}, covariance matrix for estimated parameters }
    \item{\code{null.mod}:}{Object of class \code{"S4"}, null model}
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
