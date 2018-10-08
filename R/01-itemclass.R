#Generic Item class

setClass("AllItemsClass",
         representation(par='numeric',
                        SEpar='numeric',
                        parnames='character',
                        est='logical',
                        dps='function',
                        dps2='function',
                        constr='logical',
                        itemclass='integer',
                        parnum='numeric',
                        nfact='integer',
                        nfixedeffects='numeric',
                        fixed.design='matrix',
                        dat='matrix',
                        ncat='integer',
                        gradient='numeric',
                        hessian='matrix',
                        itemtrace='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric',
                        'VIRTUAL')
)

#--------------------------------------------------------------------------

#Generics

setGeneric('ProbTrace', function(x, Theta, ...) standardGeneric("ProbTrace"))

setGeneric('ExtractLambdas', function(x) standardGeneric("ExtractLambdas"))

setGeneric('ExtractZetas', function(x) standardGeneric("ExtractZetas"))

setGeneric('Deriv', function(x, Theta, ...) standardGeneric("Deriv"))

setGeneric('DerivTheta', function(x, Theta) standardGeneric("DerivTheta"))

setGeneric('dP', function(x, Theta) standardGeneric("dP"))

setGeneric('calcLogLik', function(object, ...) standardGeneric("calcLogLik"))

setGeneric('set_null_model', function(x) standardGeneric('set_null_model'))

setGeneric("itemplot.internal",  function(object, ...) standardGeneric("itemplot.internal"))

setGeneric("fscores.internal", function(object, ...) standardGeneric("fscores.internal"))

setGeneric("summary", function(object, ...) standardGeneric("summary"))

setGeneric("vcov", function(object, ...) standardGeneric("vcov"))

setGeneric("logLik", function(object, ...) standardGeneric("logLik"))

setGeneric('DrawValues', function(x, Theta, ...) standardGeneric("DrawValues"))

setGeneric('GenRandomPars', function(x) standardGeneric("GenRandomPars"))

setGeneric('CheckIntercepts', function(x) standardGeneric("CheckIntercepts"))
