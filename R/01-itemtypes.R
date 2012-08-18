#Classes

setClass("GroupPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        parnum='numeric',
                        nfact='numeric',
                        gradient='numeric',
                        hessian='matrix')
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
                        lbound='numeric',
                        ubound='numeric',
                        'VIRTUAL')
)

setClass("NullModel", contains = 'AllItemsClass')

setClass("dich", contains = 'AllItemsClass',
         representation(a.prior='numeric', 
                        d.prior='numeric',
                        g.prior='numeric',
                        u.prior='numeric')
)

setClass("graded", contains = 'AllItemsClass',
         representation(ncat='numeric',
                        a.prior='numeric', 
                        d.prior='numeric')
)

setClass("gpcm", contains = 'AllItemsClass',
         representation(ncat='numeric',
                        a.prior='numeric', 
                        d.prior='numeric')
)

setClass("nominal", contains = 'AllItemsClass',
         representation(ncat='numeric',
                        a.prior='numeric', 
                        d.prior='numeric')
)

setClass("partcomp", contains = 'AllItemsClass',
         representation(a.prior='numeric', 
                        d.prior='numeric')
)

#--------------------------------------------------------------------------

#Generics

setGeneric('ProbTrace', function(x, Theta, ...) standardGeneric("ProbTrace"))

setGeneric('LogLik', function(x, Theta) standardGeneric("LogLik"))

setGeneric('ExtractLambdas', function(x) standardGeneric("ExtractLambdas"))

setGeneric('ExtractZetas', function(x) standardGeneric("ExtractZetas"))

setGeneric('Deriv', function(x, Theta) standardGeneric("Deriv"))

setGeneric('calcLogLik', function(object, ...) standardGeneric("calcLogLik"))         

setGeneric('ItemInfo', function(x, A, Theta) standardGeneric("ItemInfo"))
