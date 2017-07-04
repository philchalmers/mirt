setClass("AllModelClass",
         representation(Call='call',
                        Data='list',
                        Options='list',
                        Fit='list',
                        Model='list',
                        ParObjects='list',
                        OptimInfo='list',
                        Internals='list',
                        vcov='matrix',
                        time='numeric',
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
#'     \item{\code{Call}:}{function call }
#'     \item{\code{Data}:}{list of data, sometimes in different forms }
#'     \item{\code{Options}:}{list of estimation options}
#'     \item{\code{Fit}:}{a list of fit information }
#'     \item{\code{Model}:}{a list of model-based information }
#'     \item{\code{ParObjects}:}{a list of the S4 objects used during estimation}
#'     \item{\code{OptimInfo}:}{a list of arguments from the optimization process}
#'     \item{\code{Internals}:}{a list of internal arguments for secondary computations (inspecting this
#'       object is generally not required)}
#'     \item{\code{vcov}:}{a matrix represented the asymtotic covariance matrix of the parameter estimates}
#'     \item{\code{time}:}{a data.frame indicating the breakdown of computation times in seconds}
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
#' @references
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
# @examples
setClass(
    Class = 'SingleGroupClass', contains = 'AllModelClass',
    representation = representation(),
    validity = function(object) return(TRUE)
)

#------------------------------------------------------------------------------
#' Class "MultipleGroupClass"
#'
#' Defines the object returned from \code{\link{multipleGroup}}.
#'
#' @section Slots:
#'
#' \describe{
#'     \item{\code{Call}:}{function call }
#'     \item{\code{Data}:}{list of data, sometimes in different forms }
#'     \item{\code{Options}:}{list of estimation options}
#'     \item{\code{Fit}:}{a list of fit information }
#'     \item{\code{Model}:}{a list of model-based information }
#'     \item{\code{ParObjects}:}{a list of the S4 objects used during estimation}
#'     \item{\code{OptimInfo}:}{a list of arguments from the optimization process}
#'     \item{\code{Internals}:}{a list of internal arguments for secondary computations (inspecting this
#'       object is generally not required)}
#'     \item{\code{vcov}:}{a matrix represented the asymtotic covariance matrix of the parameter estimates}
#'     \item{\code{time}:}{a data.frame indicating the breakdown of computation times in seconds}
#' }
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
#' @references
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
# @examples
# @keywords classes
setClass(
    Class = 'MultipleGroupClass', contains = 'AllModelClass',
    representation = representation(),
    validity = function(object) return(TRUE)
)

#------------------------------------------------------------------------------
#' Class "MixedClass"
#'
#' Defines the object returned from \code{\link{mixedmirt}}.
#'
#' @section Slots:
#'
#' \describe{
#'     \item{\code{Call}:}{function call }
#'     \item{\code{Data}:}{list of data, sometimes in different forms }
#'     \item{\code{Options}:}{list of estimation options}
#'     \item{\code{Fit}:}{a list of fit information }
#'     \item{\code{Model}:}{a list of model-based information }
#'     \item{\code{ParObjects}:}{a list of the S4 objects used during estimation}
#'     \item{\code{OptimInfo}:}{a list of arguments from the optimization process}
#'     \item{\code{Internals}:}{a list of internal arguments for secondary computations (inspecting this
#'       object is generally not required)}
#'     \item{\code{vcov}:}{a matrix represented the asymtotic covariance matrix of the parameter estimates}
#'     \item{\code{time}:}{a data.frame indicating the breakdown of computation times in seconds}
#' }
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
#' @references
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
# @examples
# @keywords classes
setClass(
    Class = 'MixedClass', contains = 'AllModelClass',
    representation = representation(),
    validity = function(object) return(TRUE)
)

#------------------------------------------------------------------------------
#' Class "DiscreteClass"
#'
#' Defines the object returned from \code{\link{mdirt}}.
#'
#' @section Slots:
#'
#' \describe{
#'     \item{\code{Call}:}{function call }
#'     \item{\code{Data}:}{list of data, sometimes in different forms }
#'     \item{\code{Options}:}{list of estimation options}
#'     \item{\code{Fit}:}{a list of fit information }
#'     \item{\code{Model}:}{a list of model-based information }
#'     \item{\code{ParObjects}:}{a list of the S4 objects used during estimation}
#'     \item{\code{OptimInfo}:}{a list of arguments from the optimization process}
#'     \item{\code{Internals}:}{a list of internal arguments for secondary computations (inspecting this
#'       object is generally not required)}
#'     \item{\code{vcov}:}{a matrix represented the asymtotic covariance matrix of the parameter estimates}
#'     \item{\code{time}:}{a data.frame indicating the breakdown of computation times in seconds}
#' }
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
#' @references
#'
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
# @examples
# @keywords classes
setClass(
    Class = 'DiscreteClass', contains = 'AllModelClass',
    representation = representation(),
    validity = function(object) return(TRUE)
)
