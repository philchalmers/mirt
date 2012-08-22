#Classes

setClass("GroupPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        parnum='numeric',
                        nfact='numeric',
                        gradient='numeric',
                        hessian='matrix',
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
                        dat='matrix',
                        rs='matrix',
                        bfactor='logical',
                        gradient='numeric',
                        hessian='matrix',
                        method='character',
                        lbound='numeric',
                        ubound='numeric',
                        n.prior.mu='numeric',
                        n.prior.sd='numeric',
                        b.prior.alpha='numeric',
                        b.prior.beta='numeric',
                        'VIRTUAL')
)

setClass("NullModel", contains = 'AllItemsClass')

setClass("dich", contains = 'AllItemsClass')

setClass("graded", contains = 'AllItemsClass',
         representation(ncat='numeric'))

setClass("gpcm", contains = 'AllItemsClass',
         representation(ncat='numeric'))

setClass("nominal", contains = 'AllItemsClass',
         representation(ncat='numeric'))

setClass("partcomp", contains = 'AllItemsClass')

setClass("mcm", contains = 'AllItemsClass',
         representation(ncat='numeric'))

#--------------------------------------------------------------------------

#Generics

setGeneric('ProbTrace', function(x, Theta, ...) standardGeneric("ProbTrace"))

setGeneric('LogLik', function(x, Theta) standardGeneric("LogLik"))

setGeneric('ExtractLambdas', function(x) standardGeneric("ExtractLambdas"))

setGeneric('ExtractZetas', function(x) standardGeneric("ExtractZetas"))

setGeneric('Deriv', function(x, Theta) standardGeneric("Deriv"))

setGeneric('calcLogLik', function(object, ...) standardGeneric("calcLogLik"))         

setGeneric('ItemInfo', function(x, A, Theta) standardGeneric("ItemInfo"))
