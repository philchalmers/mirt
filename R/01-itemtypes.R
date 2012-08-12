#Classes

setClass("AllItemsClass",
         representation(par='numeric',
                        est='logical', 
                        constr='logical',
                        parnum='numeric',
                        nfact='numeric',
                        dat='matrix',
                        rs='matrix',
                        bfactor='logical',
                        'VIRTUAL')
)

setClass("NullModel", contains = 'AllItemsClass')

setClass("dich", contains = 'AllItemsClass',
         representation(a.prior='numeric', 
                        d.prior='numeric',
                        g.prior='numeric',
                        u.prior='numeric')
)

setClass("grad", contains = 'AllItemsClass',
         representation(ncat='numeric',
                        a.prior='numeric', 
                        d.prior='numeric')
)

setClass("gpcm", contains = 'AllItemsClass',
         representation(ncat='numeric',
                        a.prior='numeric', 
                        d.prior='numeric')
)

setClass("nom", contains = 'AllItemsClass',
         representation(ncat='numeric',
                        a.prior='numeric', 
                        d.prior='numeric')
)

#--------------------------------------------------------------------------

#Generics

setGeneric('ProbTrace', function(x, Theta) standardGeneric("ProbTrace"))

setGeneric('LogLik', function(x, Theta) standardGeneric("LogLik"))

setGeneric('ExtractLambdas', function(x) standardGeneric("ExtractLambdas"))

         