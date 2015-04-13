setClass("AllModelClass",
         representation(pars='list',
                        Data='list',
                        shortpars='numeric',
                        model='list',
                        G2='numeric',
                        df='numeric',
                        p='numeric',
                        AIC='numeric',
                        AICc='numeric',
                        K='numeric',
                        F='matrix',
                        h2='numeric',
                        Theta='matrix',
                        converge='numeric',
                        itemloc = 'numeric',
                        Prior='list',
                        BIC='numeric',
                        SABIC='numeric',
                        RMSEA='numeric',
                        null.mod = 'S4',
                        TLI = 'numeric',
                        logLik='numeric',
                        SElogLik='numeric',
                        Call='call',
                        esttype='character',
                        nest='integer',
                        iter='numeric',
                        quadpts='numeric',
                        theta_lim='numeric',
                        nfact='numeric',
                        prodlist='list',
                        constrain='list',
                        parprior='list',
                        information='matrix',
                        infomethod='character',
                        factorNames='character',
                        method='character',
                        accelerate='character',
                        Moptim='character',
                        itemtype='character',
                        time='numeric',
                        condnum='numeric',
                        secondordertest='logical',
                        empiricalhist='logical',
                        CFI='numeric',
                        CUSTOM.IND='integer',
                        SLOW.IND='integer',
                        collectLL='numeric',
                        TOL='numeric',
                        lrPars='S4',
                        exploratory='logical',
                        'VIRTUAL'),
             validity = function(object) return(TRUE)
)

#------------------------------------------------------------------------------
#' Class "SingleGroupClass"
#'
#' Defines the object returned from \code{\link{mirt}} when model is exploratory.
#'
#' @section Slots:
#'
#' \describe{
#'     \item{\code{Data}:}{Object of class \code{"list"}, contains various data matricies and properties }
#'     \item{\code{iter}:}{Object of class \code{"numeric"}, number of iterations  }
#'     \item{\code{pars}:}{Object of class \code{"list"}, estimated parameter objects list }
#'     \item{\code{shortpars}:}{Object of class \code{"numeric"}, unique estimated parameters}
#'     \item{\code{model}:}{Object of class \code{"list"}, list containing original model }
#'     \item{\code{K}:}{Object of class \code{"numeric", number of item categories}  }
#'     \item{\code{itemloc}:}{Object of class \code{"numeric", index for tabdata}  }
#'     \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
#'     \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian information criteria }
#'     \item{\code{G2}:}{Object of class \code{"numeric"}, G squared stat }
#'     \item{\code{p}:}{Object of class \code{"numeric"}, p-value for G2  }
#'     \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom }
#'     \item{\code{RMSEA}:}{Object of class \code{"numeric"}, root mean-square error of approximation for G2}
#'     \item{\code{TLI}:}{Object of class \code{"numeric"}, Tucker-Lewis index for G2}
#'     \item{\code{CFI}:}{Object of class \code{"numeric"}, CFI for G2}
#'     \item{\code{logLik}:}{Object of class \code{"numeric"}, observed log-likelihood }
#'     \item{\code{SElogLik}:}{Object of class \code{"numeric"}, Monte Carlo standard error for log-likelihood }
#'     \item{\code{F}:}{Object of class \code{"matrix"}, unrotated factor loadings }
#'     \item{\code{h2}:}{Object of class \code{"numeric"}, commonalities }
#'     \item{\code{Theta}:}{Object of class \code{"matrix"}, ability grid }
#'     \item{\code{Pl}:}{Object of class \code{"numeric"}, normed likelihoods for tabulated response}
#'     \item{\code{prodlist}:}{Object of class \code{"list"}, list containing product combination of factors }
#'     \item{\code{rotate}:}{Object of class \code{"character"}, type of rotation to be used in \code{summary}}
#'     \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
#'     \item{\code{quadpts}:}{Object of class \code{"numeric"}, number of quadrature points }
#'     \item{\code{esttype}:}{Object of class \code{"character"}, indicates whether estimation was 'EM' or 'MHRM'}
#'     \item{\code{null.mod}:}{Object of class \code{"SingleGroupClass"}, null model}
#'     \item{\code{Target}:}{Object of class \code{"numeric"}, dummy rotation matrix}
#'     \item{\code{condnum}:}{Object of class \code{"numeric"}, condition number of information matrix}
#'     \item{\code{secondordertest}:}{Object of class \code{"logical"}, indicate whether information matrix passes
#'       second-order test}
#'     \item{\code{bfactor}:}{Object of class \code{"list"}, an empty list}
#'     \item{\code{infomethod}:}{Object of class \code{"character"}, indiciates which information estimation method was used}
#'     \item{\code{TOL}:}{Object of class \code{"numeric"}, tollerence stopping criteria}
#'     \item{\code{CUSTOM.IND}:}{Object of class \code{"integer"}, an internal index}
#'     \item{\code{SLOW.IND}:}{Object of class \code{"integer"}, an internal index}
#'     \item{\code{random}:}{Object of class \code{"list"}, typicall null, except for internal mixed model usage}
#'     \item{\code{exploratory}:}{Object of class \code{"logical"}, indicates whether model contains
#'       slopes that should be rotated}
#'     \item{\code{Call}:}{Object of class \code{"call"}, call }
#' }
#' @section Methods:
#'
#' \describe{
#'     \item{anova}{\code{signature(object = "SingleGroupClass")}}
#'     \item{coef}{\code{signature(object = "SingleGroupClass")}}
#'     \item{plot}{\code{signature(x = "SingleGroupClass", y = "missing")}}
#'     \item{print}{\code{signature(x = "SingleGroupClass")} }
#'     \item{residuals}{\code{signature(object = "SingleGroupClass")}}
#'     \item{show}{\code{signature(object = "SingleGroupClass")} }
#'     \item{summary}{\code{signature(object = "SingleGroupClass")}}
#' }
#'
#' @name SingleGroupClass-class
#' @rdname SingleGroupClass-class
#' @exportClass SingleGroupClass
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords classes
# @examples
setClass(
    Class = 'SingleGroupClass', contains = 'AllModelClass',
    representation = representation(Pl='numeric',
                                    Target='numeric',
                                    rotate='character',
                                    bfactor='list',
                                    random='list'),
    validity = function(object) return(TRUE)
)

#------------------------------------------------------------------------------
#' Class "MultipleGroupClass"
#'
#' Defines the object returned from \code{\link{multipleGroup}}.
#'
#' @section Slots:
#'
#'  \describe{
#'    \item{\code{iter}:}{Object of class \code{"numeric"}, number of iterations  }
#'    \item{\code{pars}:}{Object of class \code{"list"}, estimated parameter objects list }
#'    \item{\code{shortpars}:}{Object of class \code{"numeric"}, unique estimated parameters}
#'    \item{\code{model}:}{Object of class \code{"list"}, list containing original model }
#'    \item{\code{K}:}{Object of class \code{"numeric", number of item categories}  }
#'    \item{\code{itemloc}:}{Object of class \code{"numeric", index for tabdata}  }
#'    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom}
#'    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
#'    \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian information criteria }
#'    \item{\code{G2}:}{Object of class \code{"numeric"}, G squared stat }
#'    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for G2  }
#'    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom }
#'    \item{\code{RMSEA}:}{Object of class \code{"numeric"}, root mean-square error of approximation for G2}
#'    \item{\code{TLI}:}{Object of class \code{"numeric"}, Tucker-Lewis index for G2}
#'    \item{\code{CFI}:}{Object of class \code{"numeric"}, CFI for G2}
#'    \item{\code{logLik}:}{Object of class \code{"numeric"}, observed log-likelihood }
#'    \item{\code{SElogLik}:}{Object of class \code{"numeric"}, Monte Carlo standard error for log-likelihood }
#'    \item{\code{Prior}:}{Object of class \code{"numeric"}, prior distribution used during estimation. Empty unless
#'        \code{empiricalhist = TRUE}}
#'    \item{\code{F}:}{Object of class \code{"matrix"}, unrotated factor loadings }
#'    \item{\code{h2}:}{Object of class \code{"numeric"}, commonalities }
#'    \item{\code{Theta}:}{Object of class \code{"matrix"}, ability grid }
#'    \item{\code{Pl}:}{Object of class \code{"numeric"}, normed likelihoods for tabulated response}
#'    \item{\code{prodlist}:}{Object of class \code{"list"}, list containing product combination of factors }
#'    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
#'    \item{\code{quadpts}:}{Object of class \code{"numeric"}, number of quadrature points }
#'    \item{\code{esttype}:}{Object of class \code{"character"}, indicates whether estimation was 'EM' or 'MHRM'}
#'    \item{\code{constrain}:}{Object of class \code{"list"}, list of constraints}
#'    \item{\code{invariance}:}{Object of class \code{"character"}, invariance input}
#'    \item{\code{null.mod}:}{Object of class \code{"SingleGroupClass"}, null model}
#'    \item{\code{condnum}:}{Object of class \code{"numeric"}, condition number of information matrix}
#'     \item{\code{bfactor}:}{Object of class \code{"list"}, contains information from bfactor() estimation}
#'    \item{\code{secondordertest}:}{Object of class \code{"logical"}, indicate whether information matrix passes
#'       second-order test}
#'     \item{\code{infomethod}:}{Object of class \code{"character"}, indiciates which information estimation method was used}
#'     \item{\code{TOL}:}{Object of class \code{"numeric"}, tollerence stopping criteria}
#'     \item{\code{CUSTOM.IND}:}{Object of class \code{"integer"}, an internal index}
#'     \item{\code{SLOW.IND}:}{Object of class \code{"integer"}, an internal index}
#'    \item{\code{Call}:}{Object of class \code{"call"}, call }
#'  }
#' @section Methods:
#'
#' \describe{
#'    \item{coef}{\code{signature(object = "MultipleGroupClass")}}
#'    \item{print}{\code{signature(x = "MultipleGroupClass")} }
#'    \item{show}{\code{signature(object = "MultipleGroupClass")} }
#'    \item{anova}{\code{signature(object = "MultipleGroupClass")} }
#' }
#'
#' @name MultipleGroupClass-class
#' @rdname MultipleGroupClass-class
#' @exportClass MultipleGroupClass
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords classes
# @examples
# @keywords classes
setClass(
    Class = 'MultipleGroupClass', contains = 'AllModelClass',
    representation = representation(Pl='list',
                                    invariance='character',
                                    bfactor='list'),
    validity = function(object) return(TRUE)
)

#------------------------------------------------------------------------------
#' Class "MixedClass"
#'
#' Defines the object returned from \code{\link{mixedmirt}}.
#'
#' @section Slots:
#'
#'  \describe{
#'    \item{\code{Data}:}{Object of class \code{"list"}, contains various data matricies and properties }
#'    \item{\code{iter}:}{Object of class \code{"numeric"}, number of iterations  }
#'    \item{\code{pars}:}{Object of class \code{"list"}, estimated parameter objects list }
#'    \item{\code{shortpars}:}{Object of class \code{"numeric"}, unique estimated parameters}
#'    \item{\code{model}:}{Object of class \code{"list"}, list containing original model }
#'    \item{\code{K}:}{Object of class \code{"numeric", number of item categories}  }
#'    \item{\code{itemloc}:}{Object of class \code{"numeric", index for tabdata}  }
#'    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom}
#'    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
#'    \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian information criteria }
#'    \item{\code{G2}:}{Object of class \code{"numeric"}, G squared stat }
#'    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for G2  }
#'    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom }
#'    \item{\code{RMSEA}:}{Object of class \code{"numeric"}, root mean-square error of approximation for G2}
#'    \item{\code{TLI}:}{Object of class \code{"numeric"}, Tucker-Lewis index for G2}
#'    \item{\code{CFI}:}{Object of class \code{"numeric"}, CFI for G2}
#'    \item{\code{logLik}:}{Object of class \code{"numeric"}, observed log-likelihood }
#'    \item{\code{SElogLik}:}{Object of class \code{"numeric"}, Monte Carlo standard error for log-likelihood }
#'    \item{\code{F}:}{Object of class \code{"matrix"}, unrotated factor loadings }
#'    \item{\code{h2}:}{Object of class \code{"numeric"}, commonalities }
#'    \item{\code{Theta}:}{Object of class \code{"matrix"}, ability grid }
#'    \item{\code{Pl}:}{Object of class \code{"numeric"}, normed likelihoods for tabulated response}
#'    \item{\code{prodlist}:}{Object of class \code{"list"}, list containing product combination of factors }
#'    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
#'    \item{\code{quadpts}:}{Object of class \code{"numeric"}, number of quadrature points }
#'    \item{\code{esttype}:}{Object of class \code{"character"}, indicates whether estimation was 'EM' or 'MHRM'}
#'    \item{\code{cand.t.var}:}{Object of class \code{"numeric"}, parameter used to control the MH sampler for Theta}
#'    \item{\code{random}:}{Object of class \code{"list"}, typicall null, except for internal mixed model usage}
#'    \item{\code{null.mod}:}{Object of class \code{"SingleGroupClass"}, null model}
#'    \item{\code{condnum}:}{Object of class \code{"numeric"}, condition number of information matrix}
#'    \item{\code{secondordertest}:}{Object of class \code{"logical"}, indicate whether information matrix passes
#'       second-order test}
#'     \item{\code{infomethod}:}{Object of class \code{"character"}, indiciates which information estimation method was used}
#'     \item{\code{TOL}:}{Object of class \code{"numeric"}, tollerence stopping criteria}
#'     \item{\code{CUSTOM.IND}:}{Object of class \code{"integer"}, an internal index}
#'     \item{\code{SLOW.IND}:}{Object of class \code{"integer"}, an internal index}
#'     \item{\code{formulas}:}{Object of class \code{"list"}, list of formula}
#'     \item{\code{covdata}:}{Object of class \code{"data.frame"}, covariate data}
#'     \item{\code{itemdesign}:}{Object of class \code{"data.frame"}, item design data}
#'    \item{\code{Call}:}{Object of class \code{"call"}, call }
#'  }
#' @section Methods:
#'
#'  \describe{
#'    \item{coef}{\code{signature(object = "MixedClass")}}
#'    \item{print}{\code{signature(x = "MixedClass")} }
#'    \item{residuals}{\code{signature(object = "MixedClass")}}
#'    \item{show}{\code{signature(object = "MixedClass")} }
#'    \item{summary}{\code{signature(object = "MixedClass")} }
#'    \item{logLik}{\code{signature(object = "MixedClass")} }
#'    \item{anova}{\code{signature(object = "MixedClass")} }
#'   }
#'
#' @name MixedClass-class
#' @rdname MixedClass-class
#' @exportClass MixedClass
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords classes
# @examples
# @keywords classes
setClass(
    Class = 'MixedClass', contains = 'AllModelClass',
    representation = representation(Pl='numeric',
                                    random='list',
                                    cand.t.var='numeric',
                                    formulas='list',
                                    covdata='data.frame',
                                    itemdesign='data.frame'),
    validity = function(object) return(TRUE)
)

#------------------------------------------------------------------------------
#' Class "DiscreteClass"
#'
#' Defines the object returned from \code{\link{mdirt}}.
#'
#' @section Slots:
#'
#'  \describe{
#'    \item{\code{iter}:}{Object of class \code{"numeric"}, number of iterations  }
#'    \item{\code{pars}:}{Object of class \code{"list"}, estimated parameter objects list }
#'    \item{\code{shortpars}:}{Object of class \code{"numeric"}, unique estimated parameters}
#'    \item{\code{model}:}{Object of class \code{"list"}, list containing original model }
#'    \item{\code{K}:}{Object of class \code{"numeric", number of item categories}  }
#'    \item{\code{itemloc}:}{Object of class \code{"numeric", index for tabdata}  }
#'    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom}
#'    \item{\code{AIC}:}{Object of class \code{"numeric"}, Akaike's information criteria }
#'    \item{\code{BIC}:}{Object of class \code{"numeric"}, Bayesian information criteria }
#'    \item{\code{G2}:}{Object of class \code{"numeric"}, G squared stat }
#'    \item{\code{p}:}{Object of class \code{"numeric"}, p-value for G2  }
#'    \item{\code{df}:}{Object of class \code{"numeric"}, degrees of freedom }
#'    \item{\code{RMSEA}:}{Object of class \code{"numeric"}, root mean-square error of approximation for G2}
#'    \item{\code{TLI}:}{Object of class \code{"numeric"}, Tucker-Lewis index for G2}
#'    \item{\code{CFI}:}{Object of class \code{"numeric"}, CFI for G2}
#'    \item{\code{logLik}:}{Object of class \code{"numeric"}, observed log-likelihood }
#'    \item{\code{SElogLik}:}{Object of class \code{"numeric"}, Monte Carlo standard error for log-likelihood }
#'    \item{\code{Prior}:}{Object of class \code{"numeric"}, prior distribution used during estimation. Empty unless
#'        \code{empiricalhist = TRUE}}
#'    \item{\code{F}:}{Object of class \code{"matrix"}, unrotated factor loadings }
#'    \item{\code{h2}:}{Object of class \code{"numeric"}, commonalities }
#'    \item{\code{Theta}:}{Object of class \code{"matrix"}, ability grid }
#'    \item{\code{Pl}:}{Object of class \code{"numeric"}, normed likelihoods for tabulated response}
#'    \item{\code{prodlist}:}{Object of class \code{"list"}, list containing product combination of factors }
#'    \item{\code{converge}:}{Object of class \code{"numeric"}, convergence diagnostic }
#'    \item{\code{quadpts}:}{Object of class \code{"numeric"}, number of quadrature points }
#'    \item{\code{esttype}:}{Object of class \code{"character"}, indicates whether estimation was 'EM' or 'MHRM'}
#'    \item{\code{constrain}:}{Object of class \code{"list"}, list of constraints}
#'    \item{\code{null.mod}:}{Object of class \code{"SingleGroupClass"}, null model}
#'    \item{\code{condnum}:}{Object of class \code{"numeric"}, condition number of information matrix}
#'     \item{\code{bfactor}:}{Object of class \code{"list"}, contains information from bfactor() estimation}
#'    \item{\code{secondordertest}:}{Object of class \code{"logical"}, indicate whether information matrix passes
#'       second-order test}
#'     \item{\code{infomethod}:}{Object of class \code{"character"}, indiciates which information estimation method was used}
#'     \item{\code{TOL}:}{Object of class \code{"numeric"}, tollerence stopping criteria}
#'     \item{\code{CUSTOM.IND}:}{Object of class \code{"integer"}, an internal index}
#'     \item{\code{SLOW.IND}:}{Object of class \code{"integer"}, an internal index}
#'    \item{\code{Call}:}{Object of class \code{"call"}, call }
#'  }
#' @section Methods:
#'
#' \describe{
#'    \item{print}{\code{signature(x = "DiscreteClass")} }
#'    \item{show}{\code{signature(object = "DiscreteClass")} }
#'    \item{anova}{\code{signature(object = "DiscreteClass")} }
#'    \item{coef}{\code{signature(x = "DiscreteClass")} }
#'    \item{summary}{\code{signature(object = "DiscreteClass")} }
#'    \item{residuals}{\code{signature(object = "DiscreteClass")} }
#' }
#'
#' @name DiscreteClass-class
#' @rdname DiscreteClass-class
#' @exportClass DiscreteClass
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords classes
# @examples
# @keywords classes
setClass(
    Class = 'DiscreteClass', contains = 'AllModelClass',
    representation = representation(Pl='list',
                                    bfactor='list'),
    validity = function(object) return(TRUE)
)
