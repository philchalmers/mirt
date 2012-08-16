# Class "mirtClass"
# 
# Defines the object returned from \code{\link{mirt}}.
# 
# 
# @name mirtClass-class
# @aliases mirtClass-class anova,mirtClass-method coef,mirtClass-method
# fitted,mirtClass-method plot,mirtClass,missing-method print,mirtClass-method
# residuals,mirtClass-method show,mirtClass-method summary,mirtClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("mirtClass", ...).}.
# @method Emiter number of EM iterations
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass mirtClass
# @keywords classes
setClass(
    Class = 'mirtClass',
    representation = representation(EMiter='numeric', 
                                    pars='list', 
                                    guess='numeric', 
                                    upper='numeric',
                                    K='numeric', 
                                    parsSE='list', 
                                    X2='numeric', 
                                    df='numeric', 
                                    p='numeric', 
                                    AIC='numeric', 
                                    F='matrix', 
                                    h2='numeric', 
                                    tabdata='matrix', 
                                    tabdatalong='matrix', 
                                    Theta='matrix', 
                                    Pl='numeric',
                                    data='matrix', 
                                    cormat='matrix', 
                                    facility='numeric', 
                                    converge='numeric', 
                                    itemloc = 'numeric',
                                    quadpts='numeric', 
                                    BIC='numeric', 
                                    vcov='matrix', 
                                    RMSEA='numeric', 
                                    rotate='character', 
                                    null.mod = 'S4', 
                                    TLI = 'numeric', 
                                    Target='numeric', 
                                    logLik='numeric',                                    
                                    Call='call'),	
    validity = function(object) return(TRUE)
)	

#------------------------------------------------------------------------------
# Class "bfactorClass"
# 
# Defines the object returned from \code{\link{bfactor}}.
# 
# 
# @name bfactorClass-class
# @aliases bfactorClass-class coef,bfactorClass-method
# fitted,bfactorClass-method print,bfactorClass-method
# residuals,bfactorClass-method show,bfactorClass-method
# summary,bfactorClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("bfactorClass", ...)}.
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass bfactorClass
# @keywords classes
setClass(
    Class = 'bfactorClass',
    representation = representation(EMiter = 'numeric', 
                                    pars = 'list', 
                                    upper='numeric',
                                    guess = 'numeric', 
                                    parsSE='list', 
                                    AIC = 'numeric', 
                                    X2 = 'numeric', 
                                    df = 'numeric', 
                                    logLik = 'numeric', 
                                    p = 'numeric', 
                                    F = 'matrix', 
                                    h2 = 'numeric', 
                                    itemnames = 'character', 
                                    tabdata = 'matrix', 
                                    N = 'numeric', 
                                    K='numeric',
                                    Pl = 'numeric', 
                                    Theta = 'matrix', 
                                    data = 'matrix', 
                                    itemloc = 'numeric',
                                    logicalfact = 'matrix', 
                                    facility = 'numeric', 
                                    specific = 'numeric', 
                                    tabdatalong='matrix',
                                    BIC = 'numeric', 
                                    cormat = 'matrix', 
                                    converge = 'numeric', 
                                    RMSEA = 'numeric',
                                    par.prior = 'matrix', 
                                    quadpts = 'numeric', 
                                    vcov = 'matrix', 
                                    null.mod = 'S4', 
                                    TLI = 'numeric',                                    
                                    Call = 'call'),	
    validity = function(object) return(TRUE)
)	

#------------------------------------------------------------------------------
# Class "confmirtClass"
# 
# Defines the object returned from \code{\link{confmirt}}.
# 
# 
# @name confmirtClass-class
# @aliases confmirtClass-class coef,confmirtClass-method
# print,confmirtClass-method residuals,confmirtClass-method
# show,confmirtClass-method summary,confmirtClass-method
# anova,confmirtClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("confmirtClass", ...)}.
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass confmirtClass
# @keywords classes
setClass(
    Class = 'confmirtClass',
    representation = representation(pars = 'list', 
                                    parsprint = 'matrix', 
                                    guess = 'numeric', 
                                    SEpars = 'matrix', 
                                    SEup='numeric',
                                    SEg='numeric', 
                                    gpars = 'list', 
                                    SEgpars = 'list', 
                                    estpars = 'list',
                                    cycles = 'numeric', 
                                    Theta = 'matrix', 
                                    fulldata = 'matrix', 
                                    data = 'matrix', 
                                    K = 'numeric', 
                                    itemloc = 'numeric', 
                                    h2 = 'numeric',
                                    F = 'matrix', 
                                    converge = 'numeric', 
                                    logLik = 'numeric',
                                    SElogLik = 'numeric', 
                                    df = 'integer', 
                                    AIC = 'numeric', 
                                    nconstvalues = 'integer', 
                                    G2 = 'numeric', 
                                    p = 'numeric',
                                    tabdata = 'matrix', 
                                    BIC = 'numeric', 
                                    estComp = 'logical', 
                                    prodlist = 'list', 
                                    upper = 'numeric', 
                                    RMSEA = 'numeric', 
                                    null.mod = 'S4', 
                                    TLI = 'numeric', 
                                    constrain='list',
                                    nfact='integer',
                                    Call = 'call'),	
    validity = function(object) return(TRUE)
)	
