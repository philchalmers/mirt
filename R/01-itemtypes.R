#Classes

setClass("GroupPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        parnum='numeric',
                        nfact='numeric',
                        gradient='numeric',
                        hessian='matrix',
                        itemtrace='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        n.prior.mu='numeric',
                        n.prior.sd='numeric',
                        b.prior.alpha='numeric',
                        b.prior.beta='numeric')
)

setClass("AllItemsClass",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        constr='logical',
                        parnum='numeric',
                        nfact='numeric',
                        nfixedeffects='numeric', #number of fixed effect predictors
                        dat='matrix',
                        ncat='numeric',
                        rs='matrix',
                        gradient='numeric',
                        hessian='matrix',
                        itemtrace='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        n.prior.mu='numeric',
                        n.prior.sd='numeric',
                        b.prior.alpha='numeric',
                        b.prior.beta='numeric',
                        D='numeric', #scaling correction
                        'VIRTUAL')
)

setClass("dich", contains = 'AllItemsClass')

setClass("graded", contains = 'AllItemsClass')

setClass("rating", contains = 'AllItemsClass')

setClass("gpcm", contains = 'AllItemsClass')

setClass("rsm", contains = 'AllItemsClass')

setClass("nominal", contains = 'AllItemsClass')

setClass("partcomp", contains = 'AllItemsClass')

setClass("mcm", contains = 'AllItemsClass')

setClass("nestlogit", contains = 'AllItemsClass',
         representation = representation(correctcat='integer'))

setClass('custom', contains = 'AllItemsClass',
         representation = representation(name='character',
                                         P='function',
                                         gr='function',
                                         usegr='logical',
                                         hss='function',
                                         usehss='logical',
                                         userdata='matrix',
                                         useuserdata='logical'))

#--------------------------------------------------------------------------

#Generics

setGeneric('ProbTrace', function(x, Theta, ...) standardGeneric("ProbTrace"))

setGeneric('LogLik', function(x, Theta, ...) standardGeneric("LogLik"))

setGeneric('ExtractLambdas', function(x) standardGeneric("ExtractLambdas"))

setGeneric('ExtractZetas', function(x) standardGeneric("ExtractZetas"))

setGeneric('Deriv', function(x, Theta, ...) standardGeneric("Deriv"))

setGeneric('DerivTheta', function(x, Theta) standardGeneric("DerivTheta"))

setGeneric('calcLogLik', function(object, ...) standardGeneric("calcLogLik"))

setGeneric("itemplot.internal",  function(object, ...) standardGeneric("itemplot.internal"))

setGeneric("fscores.internal", function(object, ...) standardGeneric("fscores.internal"))

setGeneric("summary", function(object, ...) standardGeneric("summary"))
