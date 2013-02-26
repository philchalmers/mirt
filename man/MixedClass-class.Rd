\name{MixedClass-class}
\docType{class}
\alias{MixedClass-class}
\alias{show,MixedClass-method}
\alias{print,MixedClass-method}

\title{Class "MixedClass"}
\description{
Defines the object returned from \code{\link{mixedmirt}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MixedClass", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{iter}:}{Object of class \code{"numeric"}, number of iterations  }
    \item{\code{pars}:}{Object of class \code{"list"}, estimated parameter objects list }    
    \item{\code{model}:}{Object of class \code{"list"}, list containing original model }    
    \item{\code{K}:}{Object of class \code{"numeric", number of item categories}  }
    \item{\code{itemloc}:}{Object of class \code{"numeric", index for tabdata}  }        
    \item{\code{G2}:}{Object of class \code{"numeric"}, G squared stat }
    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom for X2 }
    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for X2  }
    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
    \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian information criteria }
    \item{\code{RMSEA}:}{Object of class \code{"numeric"}, root mean-square error of approximation}
    \item{\code{TLI}:}{Object of class \code{"numeric"}, Tucker-Lewis index}
    \item{\code{logLik}:}{Object of class \code{"numeric"}, observed log-likelihood }
    \item{\code{SElogLik}:}{Object of class \code{"numeric"}, Monte Carlo standard error for log-likelihood }
    \item{\code{F}:}{Object of class \code{"matrix"}, unrotated factor loadings }
    \item{\code{h2}:}{Object of class \code{"numeric"}, commonalities }
    \item{\code{tabdata}:}{Object of class \code{"matrix"}, tabulated data }
    \item{\code{tabdatalong}:}{Object of class \code{"matrix"}, dichotomous version of tabulated data }    
    \item{\code{Theta}:}{Object of class \code{"matrix"}, ability grid }
    \item{\code{Pl}:}{Object of class \code{"numeric"}, normed likelihoods for tabulated response}
    \item{\code{prodlist}:}{Object of class \code{"list"}, list containing product combination of factors }
    \item{\code{data}:}{Object of class \code{"matrix"}, complete data }          
    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
    \item{\code{quadpts}:}{Object of class \code{"numeric"}, number of quadrature points }
    \item{\code{esttype}:}{Object of class \code{"character"}, indicates whether estimation was 'EM' or 'MHRM'}
    \item{\code{mixedlist}:}{Object of class \code{"list"}, contains info about the mixed effect design}
    %\item{\code{vcov}:}{Object of class \code{"matrix"}, covariance matrix for estimated parameters } 
    \item{\code{null.mod}:}{Object of class \code{"ExploratoryClass"}, null model}        
    \item{\code{Call}:}{Object of class \code{"call"}, call }
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "MixedClass")}}
    \item{print}{\code{signature(x = "MixedClass")} }
    \item{residuals}{\code{signature(object = "MixedClass")}}
    \item{show}{\code{signature(object = "MixedClass")} }
    \item{summary}{\code{signature(object = "MixedClass")} }
    \item{logLik}{\code{signature(object = "MixedClass")} }
    \item{anova}{\code{signature(object = "MixedClass")} }
    \item{fitted}{\code{signature(object = "MixedClass")} }
	 }
}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\keyword{classes}
